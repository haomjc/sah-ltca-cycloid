% run_parametric_study.m
% 论文多因素对比计算脚本
clear; clc; close all;

%% 1. 基础参数设置 (根据你现有的默认参数进行设定)
params.A   = 2.2;           % 偏心距
params.Za  = 59;            % 摆线齿数
params.Zb  = 60;            % 针齿数
params.Rz  = 168.0104;      % 针销PCR
params.Rc  = 4;             % 针齿半径
params.Drz = 343.968;       % 顶根距
params.r_zo = 4.004;        % 针齿槽半径
params.deltatheta = 0.3;    % 计算角度增量
params.actual_pin_diameter = 7.98;  % 实际使用针齿直径
params.Rrz = (params.Drz - 2*params.r_zo) / 2;

load_params.Tin = 1100;     % 输入扭矩 Nm
load_params.E = 206000;     % 弹性模量
load_params.PR = 0.3;
load_params.B = 14;         % 齿宽
load_params.angle_start = 0;
load_params.angle_end = 0; % 仅计算一个角度 (0度) 以快速验证

%% 2. 定义 5 种研究工况
% Case 1: 理想刚体基准 (无摩擦, 无刚度, 无位移)
cases(1).name = 'Ideal Rigid';
cases(1).mu_fric = 0.0;
cases(1).flags.use_kg = false;
cases(1).flags.use_pin_disp = false;

% Case 2: 仅摩擦力
cases(2).name = 'Only Friction';
cases(2).mu_fric = 0.05;
cases(2).flags.use_kg = false;
cases(2).flags.use_pin_disp = false;

% Case 3: 仅摆线轮刚度
cases(3).name = 'Only Gear Stiffness';
cases(3).mu_fric = 0.0;
cases(3).flags.use_kg = true;
cases(3).flags.use_pin_disp = false;

% Case 4: 仅针销孔壁变形(位移)
cases(4).name = 'Only Pin Disp';
cases(4).mu_fric = 0.0;
cases(4).flags.use_kg = false;
cases(4).flags.use_pin_disp = true;

% Case 5: 全耦合 (摩擦+刚度+位移)
cases(5).name = 'Full Coupled';
cases(5).mu_fric = 0.05;
cases(5).flags.use_kg = true;
cases(5).flags.use_pin_disp = true;

%% 3. 批量运行并保存结果
all_results = cell(1, length(cases));

for i = 1:length(cases)
    fprintf('\n======================================================\n');
    fprintf('正在运行 Case %d/5: %s\n', i, cases(i).name);
    fprintf('======================================================\n');
    
    % 设置摩擦力
    current_load_params = load_params;
    current_load_params.mu_fric = cases(i).mu_fric;
    
    % 运行改造后的分析程序
    res = runTCA_analysis(params, current_load_params, cases(i).flags);
    
    % 存储结果
    all_results{i} = res;
end

% 保存到工作区数据文件
save('parametric_study_results.mat', 'cases', 'all_results');
disp('✅ 所有组合计算完成，数据已保存至 parametric_study_results.mat');

% 自动调用数据导出脚本 (用于 LaTeX)
export_analysis_data;
