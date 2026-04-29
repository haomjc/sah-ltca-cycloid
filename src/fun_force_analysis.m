function F = fun_force_analysis(x,y1,y2,z,o,varargin)
dOut = x;

if nargin >= 6
    Kg_vec = varargin{1};
else
    Kg_vec = [];
end

if nargin >= 7
    flags = varargin{2};
else
    flags.use_kg = true;
    flags.use_pin_disp = true;
end

TP = y1;
BL = y2;
In = z(1);
Out = z(2);
A = z(3);
Rp_profile = z(4);
Rrp = z(5); % Generating Radius for Cycloid Profile
Nc = z(6);

Tin = z(7);
E = z(8);
B = z(9);
PR = z(10);
dRp = z(11);
dRrp = z(12);
dA= z(13);
r_zo = z(14); % 针齿槽半径
L_pin = z(15); % 针齿有效长度
Rrp_pin = z(16); % 实际针齿半径 (Operating Radius)
mu_fric = z(17); % 摩擦系数
if length(z) >= 18
    Rp = z(18); % Operating pin/hole center radius
else
    Rp = Rp_profile;
end

% ===== SAH Precondition Contracts =====
% These assertions encode the formal specification from sah-ltca.schema.json
% (AlgorithmPassport section). They do not alter computation logic.
assert(length(z) >= 17, 'SAH_PRE: parameter vector z must have >= 17 elements');
assert(all([A, Rp_profile, Rrp, Nc, Tin, B, r_zo, L_pin, Rrp_pin] > 0), ...
    'SAH_PRE: geometry and load params must be strictly positive');
assert(E > 0, 'SAH_PRE: elastic modulus E must be positive');
assert(PR >= 0 && PR <= 0.5, 'SAH_PRE: Poisson ratio PR must be in [0, 0.5]');
assert(mu_fric >= 0, 'SAH_PRE: friction coefficient must be non-negative');
assert(~isempty(TP) && ~isempty(BL), 'SAH_PRE: TP and BL arrays must not be empty');
assert(length(TP) == length(BL), 'SAH_PRE: TP and BL must have same length');
% ===== End SAH Preconditions =====

Obj = o; % 0 = 平衡方程式 1=针齿受力 2=接触半宽 3=赫兹接触应力
Np = Nc+1;
Ap = 2*pi/Np;
k1=A*Np/(Rp_profile+dRp);
Ee = 1/(2*(1-PR^2)/E);
T = @(x)(abs(x)+x)/2; %判断接触函数
C1 = pi/4*Ee*B;
EQ1 = Tin;

% 向量化计算代逻辑
idx = 1:Np;
Ap_vec = (idx - 1) * Ap;  % 从 0 开始编号

% 有效性掩码 (TP ~= NaN)
valid_mask = ~isnan(TP);

% 1. 预计算不随迭代变化的量
XS = Rp * cos(Ap_vec);
YS = Rp * sin(Ap_vec);
m_vec = tan(TP);
X0 = (A+dA)*cos(In);
Y0 = (A+dA)*sin(In);

% 力臂初始值 (基于孔中心, 后续在迭代中用实际针销位置更新)
lp = zeros(1, Np);
denominator = sqrt(m_vec.^2 + 1);

% 2. 浮动针销几何与力学模型
%
% 物理过程:
% Step 1: 计算从节点到孔中心的射线方向
% Step 2: 沿射线方向将针销偏移到刚好接触孔壁 (几何计算)
% Step 3: 在该位置计算摆线轮-针销接触力
% Step 4: 根据针销受力平衡,计算沿射线方向的额外位移 (针销压入孔壁)
% Step 5: 在新位置重新计算摆线轮-针销接触力

% 针孔参数
Clearance = r_zo - Rrp_pin; % 单边游隙 (孔半径 - 针销半径)
if Clearance < 0, Clearance = 0; end % 物理保护

% 节点位置 (Pitch Point / Node)
% 节点与偏心同向，距离为 Np*A (与 cycloidPin 一致: pitch_point = Zb*A)
node_x = Np * A * cos(In);  % 节点X (同向，倍率 Np = Zb)
node_y = Np * A * sin(In);  % 节点Y

% 针齿槽孔中心位置 (固定在针齿壳上)
hole_center_x = Rp * cos(Ap_vec);
hole_center_y = Rp * sin(Ap_vec);

% 从节点到孔中心的方向向量 (射线方向)
ray_dir_x = hole_center_x - node_x;
ray_dir_y = hole_center_y - node_y;
ray_mag = sqrt(ray_dir_x.^2 + ray_dir_y.^2);
ray_dir_x = ray_dir_x ./ ray_mag;  % 单位向量
ray_dir_y = ray_dir_y ./ ray_mag;

% Step 1: 几何偏移 - 针销沿射线方向偏移到接触孔壁
% 针销中心从孔中心偏移 Clearance 到接触孔壁
pin_offset = Clearance;  % 针销中心到孔中心的距离 = 游隙

% 针销中心的新位置 (已偏移): 沿 ray_dir (向外) 贴靠孔壁
pin_center_x = hole_center_x + ray_dir_x * pin_offset;
pin_center_y = hole_center_y + ray_dir_y * pin_offset;

% 初始化变量
fp = zeros(1, Np);
delta_ph = zeros(1, Np);
force_projection = zeros(1, Np);

% 迭代计算力和变形
% 动态收敛检测
tol = 1e-6;  % 收敛容差
max_iter = 50; % 最大迭代次数

for iter = 1:max_iter

    % Step 2: 计算当前针销位置下的摆线轮-针销接触力
    %
    % ===== 标准 LTCA 标量公式 + 2D叉积判断工作侧 =====
    % 有效变形 = (dOut - BL) * lp
    % 工作侧/背侧由 2D 叉积 (r × F) 的符号自然决定:
    %   力矩 > 0 → 工作侧 (接近针销, 提供正力矩)
    %   力矩 < 0 → 背侧 (远离针销, 不应承载)

    % 2.1 计算力臂 lp (基于当前针销位置, 始终为正)
    numerator = abs(m_vec.*X0 - Y0 + pin_center_y - m_vec.*pin_center_x);
    lp = zeros(1, Np);
    lp(valid_mask) = numerator(valid_mask) ./ denominator(valid_mask);

    % 接触点 K 坐标 (近似为针销中心沿法线向外)
    K_x = pin_center_x - Rrp_pin .* cos(TP);
    K_y = pin_center_y - Rrp_pin .* sin(TP);
    
    % 向量 PK (瞬心指向接触点)
    PK_x = K_x - node_x;
    PK_y = K_y - node_y;
    
    % 摆线轮在 K 点的滑动速度方向判定 (基于瞬心旋转)
    V_proj = PK_y .* sin(TP) + (-PK_x) .* (-cos(TP));
    S_slip = sign(V_proj);
    S_slip(abs(V_proj) < 1e-12) = 0;
    S_slip(~valid_mask) = 0;

    % 静摩擦方向: 按极限静摩擦近似，方向与相对滑动趋势相反。
    % t_vec 为法向 F_x/F_y 逆时针旋转 90 度后的切向单位方向。
    t_x = sin(TP);
    t_y = -cos(TP);
    fric_unit_x = -S_slip .* t_x;
    fric_unit_y = -S_slip .* t_y;
    fric_unit_x(~valid_mask) = 0;
    fric_unit_y(~valid_mask) = 0;

    % 摩擦力对摆线轮中心的附加力臂，稍后与法向力臂一起用于扭矩平衡。
    rK_x = K_x - X0;
    rK_y = K_y - Y0;
    friction_moment_arm = mu_fric .* (rK_x .* fric_unit_y - rK_y .* fric_unit_x);
    friction_moment_arm(~valid_mask) = 0;

    % 2.2 用 2D 叉积判断力矩方向 (r × F)
    % r = 针销中心 - 摆线轮中心 (力臂向量)
    % F = 针销对摆线轮的接触力方向 (从针销指向摆线轮, 即 TP 反方向)
    % 叉积 > 0: 针销力矩与 Tin 同向 → 工作侧 (接近侧)
    % 叉积 < 0: 针销力矩与 Tin 反向 → 背侧 (远离侧)
    r_x = pin_center_x - X0;
    r_y = pin_center_y - Y0;
    F_x = -cos(TP);   % 针销对摆线轮的力 (从针销指向摆线轮)
    F_y = -sin(TP);
    torque_dir = r_x .* F_y - r_y .* F_x;  % 2D 叉积

    % 沿用原模型约定: [-cos(TP), -sin(TP)] 表示法向接触力有效方向。
    % 刚体游隙只沿法向接近量扣减；孔壁受力投影则由法向力+静摩擦合力决定。
    normal_projection = F_x .* ray_dir_x + F_y .* ray_dir_y;
    normal_projection = max(normal_projection, 0);
    normal_projection(~valid_mask) = 0;

    friction_projection = fric_unit_x .* ray_dir_x + fric_unit_y .* ray_dir_y;
    resultant_projection = normal_projection + mu_fric .* friction_projection;
    resultant_projection = max(resultant_projection, 0);
    resultant_projection(~valid_mask) = 0;

    % 工作侧: 针销力矩方向与输入扭矩同号的针齿
    % Tin > 0 时, 工作侧的 torque_dir 应 > 0
    work_mask = (torque_dir > 0) & valid_mask;

    normal_moment_arm = torque_dir;
    normal_moment_arm(~work_mask) = 0;
    friction_moment_arm(~work_mask) = 0;
    lp_eff = normal_moment_arm + friction_moment_arm;

    % 2.3 标量变形计算
    % 总接近量 = 转角造成的接近量 - 针销在孔内贴壁前的刚体游隙损失。
    % 孔壁弹性压缩 delta_ph 作为串联柔度 C_ph 进入下方力迭代，避免重复扣减。
    clearance_loss = Clearance .* normal_projection;
    effective_displacement = (dOut - BL) .* lp - clearance_loss;

    % 只有正变形才产生接触力 (压力接触)
    effective_displacement = max(effective_displacement, 0);

    % 只有工作侧针齿才能承载 (背侧力矩为负, 排除)
    effective_displacement(~work_mask) = 0;

    effective_deformation = effective_displacement;

    % 调试输出已抑制 (scan+bisection 会调用数百次, 仅在需要时启用)

    % 计算力: 使用非线性刚度迭代 (Hertz + 齿体串联)
    for pin_idx = 1:Np
        if work_mask(pin_idx) && effective_deformation(pin_idx) > 0
            if ~isempty(Kg_vec) && Kg_vec(pin_idx) > 0
                C_struct = 1 / Kg_vec(pin_idx);
            else
                C_struct = 0;
            end
            
            f_temp = max(fp(pin_idx), 1); % 初始猜测
            
            for f_iter = 1:5
                cos_t = cos(Ap_vec(pin_idx) - In);
                s_t = 1 + k1^2 - 2*k1*cos_t;
                denom_t = k1*(1+Np)*cos_t - (1+Np*k1^2);
                p_cur = Rp_profile * s_t^1.5 / denom_t + Rrp;
                p0_cur = abs(p_cur * Rrp_pin / (p_cur + Rrp_pin));
                
                Hb_cur = (4 * f_temp * p0_cur / (pi * B * Ee))^0.5;
                if Hb_cur > 1e-9
                    log_v = log(4*abs(p_cur)/Hb_cur) + log(4*Rrp_pin/Hb_cur) - 1;
                    C_hertz = (2*(1-PR^2)/E / (pi * B)) * max(log_v, 0.1);
                else
                    C_hertz = 1e-12;
                end
                
                normal_proj_i = normal_projection(pin_idx);
                resultant_proj_i = resultant_projection(pin_idx);
                if iter > 1 && force_projection(pin_idx) > 1e-9 && normal_proj_i > 0 && resultant_proj_i > 0
                    C_ph = delta_ph(pin_idx) / force_projection(pin_idx) * normal_proj_i * resultant_proj_i;
                else
                    C_ph = 0;
                end
                
                C_total = C_struct + C_hertz + C_ph;
                if C_total > 1e-15
                    f_temp = effective_deformation(pin_idx) / C_total;
                end
            end
            fp(pin_idx) = f_temp;
        else
            fp(pin_idx) = 0;
        end
    end

    % Step 3: 针销受力平衡 - 计算针销压入孔壁的额外位移
    %
    % 摆线轮-针销接触力 fp 沿接触方向 (约等于TP方向)
    % 该力的分量沿射线方向会将针销压入孔壁
    % 孔壁提供反力, 刚度由赫兹接触公式决定

    % 力在射线方向的投影 (仅工作侧针齿有力)
    force_projection = fp .* resultant_projection;
    force_projection(~work_mask) = 0;

    % Step 4: 计算针孔变形 (赫兹接触)
    E_star = E / (2 * (1 - PR^2));
    R_eff_inv = max(1/Rrp_pin - 1/r_zo, 1e-12);  % 有效曲率半径倒数

    % 避免除零
    safe_force = max(force_projection, 1e-9);

    % 赫兹半宽
    b_width = sqrt(4 .* safe_force ./ (pi * L_pin * E_star * R_eff_inv));

    % 赫兹变形 (修正为标准圆柱赫兹接近量公式)
    valid_idx = (b_width > 1e-9) & work_mask;
    delta_ph_new = zeros(1, Np);

    if any(valid_idx) && flags.use_pin_disp
        factor = 2 .* safe_force(valid_idx) ./ (pi * L_pin * E_star);
        % 标准公式: ln(4*R1/b) + ln(4*R2/b) - 1
        term_log = log(4 * Rrp_pin ./ b_width(valid_idx)) + log(4 * r_zo ./ b_width(valid_idx));
        delta_ph_new(valid_idx) = factor .* max(term_log - 1, 0);
    end

    % 更新针销中心位置 (孔壁弹性变形使针销进一步外移)
    pin_center_x = hole_center_x + ray_dir_x .* (Clearance + delta_ph_new);
    pin_center_y = hole_center_y + ray_dir_y .* (Clearance + delta_ph_new);

    % 更新针孔变形
    diff_ph = max(abs(delta_ph_new - delta_ph));
    delta_ph = delta_ph_new;
    if diff_ph < tol
        break;
    end

    if iter >= max_iter
        warning('SAH_CONV: pin-hole iteration did not converge (max_iter=%d, residual=%.2e)', max_iter, diff_ph);
    end


end

% 3. 平衡方程剩余力矩
% 使用带方向的力臂: 背侧针齿已被 work_mask 排除 (fp=0)
% 保留 lp 的正值 (work_mask 保证只有正力矩方向的针齿有 fp)
EQ1 = Tin - sum(fp .* lp_eff);

% 4. 赫兹应力相关 (摆线-针销侧)
cos_term = cos(Ap_vec - In);
s = 1 + k1^2 - 2*k1*cos_term;
denom = k1*(1+Np)*cos_term - (1+Np*k1^2);

% 摆线齿廓曲率半径 (等距曲线)
% 基础摆线曲率 + 等距量Rrp (生成半径)
% Rrp 决定齿廓几何, 与实际针销尺寸无关
p = Rp_profile * s.^1.5 ./ denom + Rrp;

% 当量曲率半径 (赫兹接触: 摆线齿廓 vs 实际针销)
p0 = abs(p .* Rrp_pin ./ (p + Rrp_pin));

% 赫兹接触半宽 Hb
Hb = (4 * fp .* p0 / (pi * B * Ee)).^0.5;

% 赫兹接触应力 Hz
Hz = zeros(1, Np);
active_idx = find(fp > 0 & Hb > 0);
if ~isempty(active_idx)
    Hz(active_idx) = 2 * fp(active_idx) ./ (pi * Hb(active_idx) * B);
end

% 接触刚度 k (赫兹接近量)
% 注意: 第二个 log 项应使用实际针销半径 Rrp_pin
k = zeros(1, Np);
if ~isempty(active_idx)
    term = 2*(1-PR^2)/E;
    log_val = log(4*abs(p(active_idx))./Hb(active_idx)) + log(4*Rrp_pin./Hb(active_idx)) - 1;
    k(active_idx) = abs(pi * B ./ (term * log_val));
end


% ===== SAH Postcondition Contracts =====
% Verify physical invariants after convergence.
assert(all(fp >= -1e-6), 'SAH_POST: pin forces fp must be non-negative');
assert(all(fp(~work_mask) == 0), 'SAH_POST: backside (non-working) pins must carry zero force');
assert(isfinite(EQ1), 'SAH_POST: torque balance residual EQ1 must be finite');
assert(all(Hb >= 0), 'SAH_POST: Hertz half-width Hb must be non-negative');
assert(all(Hz >= 0), 'SAH_POST: Hertz contact stress Hz must be non-negative');
% ===== End SAH Postconditions =====

if Obj == 0
    F = (EQ1);
elseif Obj ==1
    F = fp;
elseif Obj ==2
    F =Hb;
elseif Obj ==3
    F =k;
elseif Obj ==4
    F =Hz;
elseif Obj == 5
    F = sum(fp .* lp_eff);
end

end
