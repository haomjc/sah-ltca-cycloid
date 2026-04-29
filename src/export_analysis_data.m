% export_analysis_data.m
% 导出用于 LaTeX Pgfplots 绘图的 .dat 文本数据文件
clear; clc;
load('parametric_study_results.mat'); % 加载上面脚本跑出来的结果

output_dir = fullfile(pwd, 'manuscript', 'data');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% 指定要观察的特征针齿编号 (例如 1号齿)
target_pin = 1;

%% 1. 导出 传动误差 (Loaded Transmission Error) 对比数据
filename_te = fullfile(output_dir, 'compare_TE.dat');
fid_te = fopen(filename_te, 'w');
fprintf(fid_te, 'Angle\tCase1\tCase2\tCase3\tCase4\tCase5\n'); % 表头

angles = all_results{1}.angles_deg;
for i = 1:length(angles)
    fprintf(fid_te, '%.3f\t', angles(i));
    for c = 1:5
        fprintf(fid_te, '%.4f', all_results{c}.Err_loaded(i));
        if c < 5, fprintf(fid_te, '\t'); else, fprintf(fid_te, '\n'); end
    end
end
fclose(fid_te);
disp(['✅ 传动误差数据已导出至: ', filename_te]);

%% 2. 导出 最大单齿接触力 (Force on specific pin) 对比数据
filename_force = fullfile(output_dir, 'compare_Force_Pin1.dat');
fid_f = fopen(filename_force, 'w');
fprintf(fid_f, 'Angle\tCase1\tCase2\tCase3\tCase4\tCase5\n');

for i = 1:length(angles)
    fprintf(fid_f, '%.3f\t', angles(i));
    for c = 1:5
        fprintf(fid_f, '%.4f', all_results{c}.FP(i, target_pin));
        if c < 5, fprintf(fid_f, '\t'); else, fprintf(fid_f, '\n'); end
    end
end
fclose(fid_f);
disp(['✅ 针齿受力对比数据已导出至: ', filename_force]);

%% 3. 导出 最大赫兹接触应力对比数据 (取每一时刻所有齿中的最大值)
filename_hz = fullfile(output_dir, 'compare_MaxHz.dat');
fid_hz = fopen(filename_hz, 'w');
fprintf(fid_hz, 'Angle\tCase1\tCase2\tCase3\tCase4\tCase5\n');

for i = 1:length(angles)
    fprintf(fid_hz, '%.3f\t', angles(i));
    for c = 1:5
        max_hz = max(all_results{c}.Hz(i, :));
        fprintf(fid_hz, '%.4f', max_hz);
        if c < 5, fprintf(fid_hz, '\t'); else, fprintf(fid_hz, '\n'); end
    end
end
fclose(fid_hz);
disp(['✅ 最大接触应力对比数据已导出至: ', filename_hz]);

%% 4. 导出 扭转刚度对比数据 (取平均值输出到命令行，因为通常是单值)
disp('--- 扭转刚度 (平均值) 对比 ---');
for c = 1:5
    valid_ktm = all_results{c}.Ktm(~isnan(all_results{c}.Ktm) & ~isinf(all_results{c}.Ktm));
    avg_ktm = mean(valid_ktm);
    disp([cases(c).name, ': ', num2str(avg_ktm), ' Nm/arcmin']);
end
