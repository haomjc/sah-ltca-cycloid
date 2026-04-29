function results = runTCA_analysis(params, load_params, analysis_flags)
% RUNTCA_ANALYSIS 运行齿面接触分析 (TCA) 和承载接触分析 (LTCA)
% 增加多因素对比分析支持 (通过 analysis_flags)
%
% 输入:
%   params: 包含摆线针轮几何参数的结构体 (来自 cycloidPin)
%   load_params: 包含工况参数的结构体
%   analysis_flags: 控制摩擦、刚度、位移的开关

if nargin < 3
    analysis_flags.use_kg = true;
    analysis_flags.use_pin_disp = true;
end

%% 1. 参数提取与适配
% 几何参数
Drz = params.Drz;
r_zo = params.r_zo;
Rp_profile = params.Rz;
if isfield(params, 'Rrz')
    Rp_pin_center = params.Rrz;
else
    Rp_pin_center = Drz/2 - r_zo;
end
Rrp = params.Rc;
Rrp_pin = params.actual_pin_diameter / 2;
Nc = params.Za;
Np = Nc + 1;
A = params.A;

gear_ratio = abs(Np - Nc) / Nc;        % 理论传动比大小

% 负载与摩擦参数
Tin = load_params.Tin * 1000;          % Nm -> Nmm
E = load_params.E;
PR = load_params.PR;
B = load_params.B;
if isfield(load_params, 'L_pin')
    L_pin = load_params.L_pin;
else
    L_pin = B;
end
if isfield(load_params, 'mu_fric')
    mu_fric = load_params.mu_fric;
else
    mu_fric = 0.05;
end

% 误差参数
dRp = 0.00;
dRrp = 0.00;
dA = 0.00;

Ap = 2*pi/Np;

% 角度参数
start_angle = load_params.angle_start;
end_angle = load_params.angle_end;
step_angle = params.deltatheta;
range = end_angle - start_angle;
num_steps = floor(range/step_angle) + 1;

mod_params = buildModParams(params);

%% 1.5 初始化并行池与进度条机制
% 检查是否已有并池，如果没有则在后台启动
pool = gcp('nocreate');
if isempty(pool)
    hWaitInit = waitbar(0, '首次运行正在启动多核并行引擎 (Parpool)，请稍候...');
    parpool;
    close(hWaitInit);
end

hWait = waitbar(0, '准备并行计算...');
dq = parallel.pool.DataQueue;
total_tasks = num_steps * 3; % DCA(1) + 刚度(1) + LTCA(1)
tasks_completed = 0;

% 定义回调函数，用于在并行过程中安全更新进度条
afterEach(dq, @updateWaitbar);

    function updateWaitbar(~)
        tasks_completed = tasks_completed + 1;
        if isvalid(hWait)
            waitbar(tasks_completed / total_tasks, hWait, ...
                sprintf('多核并行计算中... 总体进度: %d / %d', tasks_completed, total_tasks));
        end
    end

%% 2. 预分配变量 (为 parfor 准备)
IN = zeros(1, num_steps);
OUT_theo = zeros(1, num_steps);

for j = 1:num_steps
    input_angle_deg = start_angle + (j-1) * step_angle;
    IN(j) = deg2rad(input_angle_deg + 90);
    OUT_theo(j) = IN(j) * gear_ratio;
end

TO = zeros(num_steps, Np);
TC = zeros(num_steps, Np);
TP = zeros(num_steps, Np);
BL = zeros(num_steps, Np);

%% 3. 齿面接触分析 (DCA) - 并行执行
fprintf('\n===== DCA 并行计算开始 =====\n');

% fsolve 选项预设
opt_dca = optimoptions('fsolve', 'MaxFunEvals', 30, 'MaxIter', 10, 'Display', 'off');

parfor j = 1:num_steps
    theo_out_j = OUT_theo(j);
    IN_j = IN(j);
    
    row_TO = zeros(1, Np);
    row_TC = zeros(1, Np);
    row_TP = zeros(1, Np);

    for i = 1:Np
        pin_angle = (i-1) * Ap;
        iTP = atan2(A*Np*sin(IN_j)-(Rp_pin_center)*sin(pin_angle), ...
            A*Np*cos(IN_j)-(Rp_pin_center)*cos(pin_angle));

        iTC = fun_findtc([IN_j, theo_out_j, iTP, A, Rp_profile, Rrp, Nc, i-1, dRp, dRrp, Rrp_pin, Rp_pin_center]);

        initial = [theo_out_j, iTC, iTP];
        y_args = [IN_j, A, Rp_profile, dRp, Rrp, Nc, i-1, dA, Rrp_pin, Rp_pin_center];
        
        [S, ~, exitflag] = fsolve(@(x) fun_output_mod(x, y_args, mod_params), initial, opt_dca);

        if exitflag > 0
            row_TO(i) = S(1);
            row_TC(i) = S(2);
            row_TP(i) = S(3);
        else
            row_TO(i) = NaN;
            row_TC(i) = NaN;
            row_TP(i) = NaN;
        end
    end
    
    TO(j,:) = row_TO;
    TC(j,:) = row_TC;
    TP(j,:) = row_TP;
    
    send(dq, j); % 发送进度信号
end

%% 4. 计算背隙/过盈与几何误差 (串行极快，无需并行)
OUT_geo = zeros(1, num_steps);
Err_geo = zeros(1, num_steps);
Clearance = max(r_zo - Rrp_pin, 0);

for j = 1:num_steps
    for i = 1:Np
        if ~isnan(TO(j,i))
            BL(j,i) = TO(j,i) - OUT_theo(j);
        else
            BL(j,i) = NaN;
        end
    end

    valid_TO_indices = find(~isnan(TO(j,:)));
    if ~isempty(valid_TO_indices)
        [min_val, min_idx] = min(TO(j, valid_TO_indices));
        actual_pin_idx = valid_TO_indices(min_idx);

        OUT_geo(j) = min_val;

        tp_val = TP(j, actual_pin_idx);
        m_vec = tan(tp_val);
        X0 = (A+dA)*cos(IN(j));
        Y0 = (A+dA)*sin(IN(j));

        pin_angle_rad = (actual_pin_idx - 1) * Ap;
        XS_geo = Rp_pin_center * cos(pin_angle_rad);
        YS_geo = Rp_pin_center * sin(pin_angle_rad);

        num = abs(m_vec*X0 - Y0 + YS_geo - m_vec*XS_geo);
        den = sqrt(m_vec^2 + 1);
        lp_val = num / den;

        if lp_val < 1e-3, lp_val = Rp_pin_center; end

        clearance_angle_rad = Clearance / lp_val;
        Err_geo(j) = (radtodeg(OUT_theo(j) - OUT_geo(j)) + radtodeg(clearance_angle_rad)) * 3600;
    else
        OUT_geo(j) = OUT_theo(j);
        Err_geo(j) = 0;
    end
end

%% 4.5 预计算齿体刚度 (能量法) - 并行执行
fprintf('\n===== 齿体能量法刚度并行计算 =====\n');
Kg_matrix = zeros(num_steps, Np);
Ap_vec = (0:Np-1) * Ap;

parfor j = 1:num_steps
    IN_j = IN(j);
    row_Kg = zeros(1, Np);
    
    for i = 1:Np
        if ~isnan(TP(j,i)) && ~isnan(TO(j,i)) && ~isnan(TC(j,i))
            theta_ci = mod(TC(j,i) * Nc, 2*pi);
            if theta_ci > pi, theta_ci = 2*pi - theta_ci; end
            phi0 = theta_ci / Nc;
            row_Kg(i) = fun_tooth_stiffness(phi0, Rp_profile, A, Nc, Np, Rrp, B, E, PR, mod_params);
        end
    end
    Kg_matrix(j,:) = row_Kg;
    send(dq, j); % 发送进度信号
end

%% 5. 承载接触分析 (LTCA) - 并行执行
fprintf('\n===== LTCA 并行计算开始 =====\n');

FP = zeros(num_steps, Np);
Hb = zeros(num_steps, Np);
k = zeros(num_steps, Np);
Hz = zeros(num_steps, Np);
dphi = zeros(num_steps, 1);
OUT_loaded = zeros(1, num_steps);

opt_ltca = optimset('Display','off','TolX',1e-12);

parfor j = 1:num_steps
    IN_j = IN(j);
    OUT_j = OUT_theo(j);
    
    force_params = [IN_j, OUT_j, A, Rp_profile, Rrp, Nc, Tin, E, B, PR, ...
        dRp, dRrp, dA, r_zo, L_pin, Rrp_pin, mu_fric, Rp_pin_center];

    TP_j = TP(j,:);
    BL_j = BL(j,:);
    Kg_j = Kg_matrix(j,:);
    
    F_obj = @(x) fun_force_analysis(x, TP_j, BL_j, force_params, 0, Kg_j, analysis_flags);

    dOut_scan = linspace(-1e-4, 1e-4, 500);
    F_scan = arrayfun(F_obj, dOut_scan);
    sign_idx = find(F_scan(1:end-1) .* F_scan(2:end) <= 0, 1);

    if isempty(sign_idx)
        dOut_scan = linspace(-0.01, 0.01, 2000);
        F_scan = arrayfun(F_obj, dOut_scan);
        sign_idx = find(F_scan(1:end-1) .* F_scan(2:end) <= 0, 1);
    end

    if ~isempty(sign_idx)
        bracket = [dOut_scan(sign_idx), dOut_scan(sign_idx+1)];
        dOut_sol = fzero(F_obj, bracket, opt_ltca);
    else
        warning('runTCA:LTCA failed to find root at step %d', j);
        dOut_sol = 0;
    end

    FP(j,:) = fun_force_analysis(dOut_sol, TP_j, BL_j, force_params, 1, Kg_j, analysis_flags);
    Hb(j,:) = fun_force_analysis(dOut_sol, TP_j, BL_j, force_params, 2, Kg_j, analysis_flags);
    k(j,:) = fun_force_analysis(dOut_sol, TP_j, BL_j, force_params, 3, Kg_j, analysis_flags);
    Hz(j,:) = fun_force_analysis(dOut_sol, TP_j, BL_j, force_params, 4, Kg_j, analysis_flags);

    OUT_loaded(j) = OUT_j - dOut_sol;
    dphi(j) = radtodeg(dOut_sol) * 3600;
    
    send(dq, j); % 发送进度信号
end

if isvalid(hWait)
    close(hWait);
end

%% 6. 整理输出结果
Err_loaded = radtodeg(OUT_theo - OUT_loaded) * 3600;

dphi_safe = abs(dphi);
dphi_safe(dphi_safe == 0) = NaN;
Ktm = Tin ./ (dphi_safe/60) / 1000;

results.angles_deg = rad2deg(IN);
results.Err_geo = Err_geo;
results.Err_loaded = Err_loaded;
results.dphi = dphi;
results.FP = FP;
results.Hz = Hz;
results.Hb = Hb;
results.k = k;
results.Ktm = Ktm;
results.BL = BL;
results.TO = TO;
results.TP = TP;
results.IN = IN;
results.OUT = OUT_loaded;

fprintf('===== TCA 分析完成 =====\n');

end

%% ================= 辅助函数 ================= %%

function mod_params = buildModParams(params)
mod_params.type = 0;
mod_params.theta = [];
mod_params.gamma = [];

if isfield(params, 'shape_type')
    switch params.shape_type
        case 1
            mod_params.type = 1;
            mod_params.theta = [params.theta1_3, params.theta2_3, params.theta3_3, params.theta4_3];
            mod_params.gamma = [params.gamma1_3, params.gamma2_3];
        case 2
            mod_params.type = 2;
            mod_params.theta = [params.theta1_4, params.theta2_4, params.theta3_4, params.theta4_4, params.theta5_4, params.theta6_4];
            mod_params.gamma = [params.gamma1_4, params.gamma2_4];
        case 3
            mod_params.type = 3;
            mod_params.theta = [params.theta1_5, params.theta2_5, params.theta3_5, params.theta4_5, params.theta5_5, params.theta6_5, params.theta7_5];
            mod_params.gamma = [params.gamma1_5, params.gamma2_5, params.gamma3_5];
    end
end
end
