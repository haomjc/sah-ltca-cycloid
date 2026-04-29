function Kg = fun_tooth_stiffness(phi0, Rz, A, Nc, Np, Rc, B, E, PR, mod_params)
% FUN_TOOTH_STIFFNESS 基于能量法计算摆线轮齿体结构刚度
% 采用修型后的实际齿廓坐标，并用数值差分计算齿廓导数。
%
% 输入:
%   phi0      - 接触点齿内参数 (rad)
%   Rz        - 针齿分布圆半径 / 齿廓生成半径 (mm)
%   A         - 偏心距 (mm)
%   Nc        - 摆线齿数 zg
%   Np        - 针齿数 zb (= Nc + 1)
%   Rc        - 名义针齿半径 / 等距半径 (mm)
%   B         - 齿宽 (mm)
%   E         - 弹性模量 (MPa = N/mm^2)
%   PR        - 泊松比
%   mod_params - 修型参数结构体，可省略
%
% 输出:
%   Kg        - 摆线轮齿体结构刚度 (N/mm)

if nargin < 10 || isempty(mod_params)
    mod_params = struct('type', 0, 'theta', [], 'gamma', []);
end

Zb = Np;
Zg = Nc;
K1 = A * Zb / Rz;
beta = pi / Zg;
G = E / (2 * (1 + PR));

% 避开齿根和齿顶附近的奇异区。phi0 过小时等效为极短梁，刚度极大。
phi_min = pi / (Zg^2);
phi_max = (Zg - 1) * pi / (Zg^2);
phi0 = min(max(phi0, phi_min), phi_max);

[x2_0, y2_0, dx2_0, dy2_0] = profile_local(phi0);
if ~isfinite(x2_0) || ~isfinite(y2_0) || ~isfinite(dx2_0) || ~isfinite(dy2_0)
    Kg = 1e12;
    return;
end

alpha0 = pi/2 - atan2(dy2_0, dx2_0);
cos_a0 = cos(alpha0);
sin_a0 = sin(alpha0);

if phi0 <= phi_min + 1e-9
    Kg = 1e12;
    return;
end

try
    total_inv = integral(@integrand, phi_min, phi0, 'RelTol', 1e-6, 'ArrayValued', true);
    if total_inv > 0 && isfinite(total_inv)
        Kg = 1 / total_inv;
    else
        Kg = 1e12;
    end
catch
    Kg = 1e12;
end

    function val = integrand(phi)
        [x2, y2, ~, dy2_dphi] = profile_local(phi);

        x2_abs = max(abs(x2), 1e-9);
        dy_abs = abs(dy2_dphi);

        % 论文式(11): 压力角取接触点 alpha0，截面半高以 |x2| 代入。
        moment = (y2_0 - y2) .* cos_a0 + x2_0 .* sin_a0;
        inv_Kb = 3 .* moment.^2 .* dy_abs ./ (2 * E * B .* x2_abs.^3);
        inv_Ks = 1.2 * cos_a0^2 .* dy_abs ./ (2 * G * B .* x2_abs);
        inv_Kc = sin_a0^2 .* dy_abs ./ (2 * E * B .* x2_abs);

        val = inv_Kb + inv_Ks + inv_Kc;
        val(~isfinite(val)) = 0;
    end

    function [x2, y2, dx2_dphi, dy2_dphi] = profile_local(phi)
        [x1, y1] = profile_xy(phi);

        h = max(1e-6, 1e-4 * max(abs(phi0), phi_min));
        phi_m = max(phi - h, phi_min * 0.1);
        phi_p = min(phi + h, phi_max * 1.05);
        if phi_p <= phi_m
            phi_p = phi + h;
            phi_m = phi - h;
        end

        [x1_p, y1_p] = profile_xy(phi_p);
        [x1_m, y1_m] = profile_xy(phi_m);
        dx1_dphi = (x1_p - x1_m) ./ (phi_p - phi_m);
        dy1_dphi = (y1_p - y1_m) ./ (phi_p - phi_m);

        x2 = x1 .* cos(beta) - y1 .* sin(beta);
        y2 = x1 .* sin(beta) + y1 .* cos(beta) - (Rz - Rc - A) .* cos(beta);

        dx2_dphi = dx1_dphi .* cos(beta) - dy1_dphi .* sin(beta);
        dy2_dphi = dx1_dphi .* sin(beta) + dy1_dphi .* cos(beta);
    end

    function [x1, y1] = profile_xy(phi)
        denom = sqrt(1 + K1^2 - 2 * K1 .* cos(Zg .* phi));
        denom = max(denom, 1e-12);

        x0 = Rz .* (sin(phi) - (K1 / Zb) .* sin(Zb .* phi));
        y0 = Rz .* (cos(phi) - (K1 / Zb) .* cos(Zb .* phi));

        normal_x = (K1 .* sin(Zb .* phi) - sin(phi)) ./ denom;
        normal_y = (-K1 .* cos(Zb .* phi) + cos(phi)) ./ denom;

        Rc_eff = Rc + modification_amount(phi);
        x1 = x0 + Rc_eff .* normal_x;
        y1 = y0 - Rc_eff .* normal_y;
    end

    function dRc = modification_amount(phi)
        dRc = zeros(size(phi));
        if ~isfield(mod_params, 'type') || mod_params.type == 0
            return;
        end

        theta_deg = rad2deg(phi .* Zg);
        theta_norm = mod(theta_deg, 360);
        theta_norm(theta_norm > 180) = theta_norm(theta_norm > 180) - 360;
        abs_theta = abs(theta_norm);

        t = mod_params.theta;
        g = mod_params.gamma;

        switch mod_params.type
            case 1
                idx = abs_theta >= t(1) & abs_theta < t(2);
                dRc(idx) = g(1) .* (t(2) - abs_theta(idx)) ./ (t(2) - t(1));
                idx = abs_theta >= t(3) & abs_theta <= t(4);
                dRc(idx) = g(2) .* (abs_theta(idx) - t(3)) ./ (t(4) - t(3));

            case 2
                idx = abs_theta >= t(1) & abs_theta < t(2);
                dRc(idx) = g(1);
                idx = abs_theta >= t(2) & abs_theta < t(3);
                dRc(idx) = g(1) .* (t(3) - abs_theta(idx)) ./ (t(3) - t(2));
                idx = abs_theta >= t(4) & abs_theta < t(5);
                dRc(idx) = g(2) .* (abs_theta(idx) - t(4)) ./ (t(5) - t(4));
                idx = abs_theta >= t(5) & abs_theta <= t(6);
                dRc(idx) = g(2);

            case 3
                idx = abs_theta >= t(1) & abs_theta < t(2);
                dRc(idx) = g(1);
                idx = abs_theta >= t(2) & abs_theta < t(3);
                dRc(idx) = g(1) .* (t(3) - abs_theta(idx)) ./ (t(3) - t(2));
                idx = abs_theta >= t(4) & abs_theta < t(5);
                dRc(idx) = g(2) .* (abs_theta(idx) - t(4)) ./ (t(5) - t(4));
                idx = abs_theta >= t(5) & abs_theta < t(6);
                dRc(idx) = g(2) + (g(3) - g(2)) .* (abs_theta(idx) - t(5)) ./ (t(6) - t(5));
                idx = abs_theta >= t(6) & abs_theta <= t(7);
                dRc(idx) = g(3);
        end
    end
end
