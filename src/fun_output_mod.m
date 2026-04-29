function F = fun_output_mod(x, y, mod_params)
% 齿面接触分析函数 - 支持动态修形
% x: [Out, TC, TP] - 求解变量
% y: [In, A, Rp_profile, dRp, Rrp, Nc, N, dA, Rrp_pin, Rp_pin_center]
% mod_params: 修形参数结构体
%   .type: 修形类型 (1=3段, 2=4段, 3=5段)
%   .theta: 角度节点数组
%   .gamma: 修形量数组
%
% 修形量 dRrp 由 TC 动态计算

syms TP_sym Alpha_sym;

Out = x(1);
TC = x(2);

In = y(1);
A = y(2);
Rp = y(3);
dRp = y(4);
Rrp = y(5);
Nc = y(6);
N = y(7);
dA = y(8);
if length(y) >= 9
    Rrp_pin = y(9);
else
    Rrp_pin = Rrp; % Fallback
end
if length(y) >= 10
    Rp_pin_center = y(10);
else
    Rp_pin_center = Rp; % Fallback
end

Np = Nc + 1;
Ap = 2*pi/Np;

% ===== 计算动态修形量 dRrp =====
% TC -> theta_deg 映射
% TC (rad) 范围 0..2π 对应整个摆线轮
% 一个齿距 = 2π/Nc, 对应修形角度 0..360°
theta_deg = rad2deg(TC * Nc);

% 映射到 -180..180 区间
theta_norm = mod(theta_deg, 360);
if theta_norm > 180
    theta_norm = theta_norm - 360;
end
abs_theta = abs(theta_norm);

% 根据修形类型计算 gamma
dRrp = 0;
if isfield(mod_params, 'type')
    switch mod_params.type
        case 1 % 3段修形
            t = mod_params.theta; % [t1, t2, t3, t4]
            g = mod_params.gamma; % [g1, g2]
            if abs_theta >= t(1) && abs_theta < t(2)
                dRrp = g(1) * (t(2) - abs_theta) / (t(2) - t(1));
            elseif abs_theta >= t(2) && abs_theta < t(3)
                dRrp = 0;
            elseif abs_theta >= t(3) && abs_theta <= t(4)
                dRrp = g(2) * (abs_theta - t(3)) / (t(4) - t(3));
            end

        case 2 % 4段修形
            t = mod_params.theta; % [t1, t2, t3, t4, t5, t6]
            g = mod_params.gamma; % [g1, g2]
            if abs_theta >= t(1) && abs_theta < t(2)
                dRrp = g(1);
            elseif abs_theta >= t(2) && abs_theta < t(3)
                dRrp = g(1) * (t(3) - abs_theta) / (t(3) - t(2));
            elseif abs_theta >= t(3) && abs_theta < t(4)
                dRrp = 0;
            elseif abs_theta >= t(4) && abs_theta < t(5)
                dRrp = g(2) * (abs_theta - t(4)) / (t(5) - t(4));
            elseif abs_theta >= t(5) && abs_theta <= t(6)
                dRrp = g(2);
            end

        case 3 % 5段修形
            t = mod_params.theta; % [t1, t2, t3, t4, t5, t6, t7]
            g = mod_params.gamma; % [g1, g2, g3]
            if abs_theta >= t(1) && abs_theta < t(2)
                dRrp = g(1);
            elseif abs_theta >= t(2) && abs_theta < t(3)
                dRrp = g(1) * (t(3) - abs_theta) / (t(3) - t(2));
            elseif abs_theta >= t(3) && abs_theta < t(4)
                dRrp = 0;
            elseif abs_theta >= t(4) && abs_theta < t(5)
                dRrp = g(2) * (abs_theta - t(4)) / (t(5) - t(4));
            elseif abs_theta >= t(5) && abs_theta < t(6)
                dRrp = g(2) + (g(3) - g(2)) * (abs_theta - t(5)) / (t(6) - t(5));
            elseif abs_theta >= t(6) && abs_theta <= t(7)
                dRrp = g(3);
            end
    end
end

% ===== 以下与原始 fun_output 完全一致 =====
% rC (Cycloid Profile) uses Generating Radius (Rrp + dRrp)
rC = [(Rp+dRp)*cos(TC)+(Rrp+dRrp)*cos(Alpha_sym-TC)-A*cos(Np*TC);
    (Rp+dRp)*sin(TC)-(Rrp+dRrp)*sin(Alpha_sym-TC)-A*sin(Np*TC);
    0;
    1];

T = [cos(Out)   sin(Out)    0   (A+dA)*cos(In);
    -sin(Out)   cos(Out)    0   (A+dA)*sin(In);
    0          0           1   0;
    0          0           0   1];

K = [0 0 1];
% frR (Pin Surface) uses operating center radius and actual pin radius
frR = [Rp_pin_center*cos(N*Ap)+Rrp_pin*cos(TP_sym), Rp_pin_center*sin(N*Ap)+Rrp_pin*sin(TP_sym), 0];
frC = T*rC;

dfrR = diff(frR);
drC = -diff(rC');
drC = drC(1:3);

nC = cross(drC,K)/norm(cross(drC,K));
fnR = cross(dfrR,K)/norm(cross(dfrR,K));
fnC = T(1:3,1:3)*nC';

TP_val = x(3);
Alpha_val = atan2(-sin(Nc*TC), (cos(Nc*TC)-(Rp+dRp)/(A*Np)));

frC = double(subs(frC, [TP_sym, Alpha_sym], [TP_val, Alpha_val]));
fnC = double(subs(fnC, [TP_sym, Alpha_sym], [TP_val, Alpha_val]));
frR = double(subs(frR, [TP_sym, Alpha_sym], [TP_val, Alpha_val]));
fnR = double(subs(fnR, [TP_sym, Alpha_sym], [TP_val, Alpha_val]));

F = [frC(1)-frR(1); frC(2)-frR(2); atan2(fnC(2),fnC(1))-atan2(fnR(2),fnR(1))];
end
