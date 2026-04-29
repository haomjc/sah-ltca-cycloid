function F = fun_findtc(y)
% 寻找TC参数的函数
In = y(1);
Out = y(2);
TP = y(3);
A = y(4);
Rp = y(5);
Rrp = y(6); % Generating Radius
Nc = y(7);
N = y(8);
dRp = y(9);
dRrp = y(10);
if length(y) >= 11
    Rrp_pin = y(11); % Operating Radius
else
    Rrp_pin = Rrp;
end
if length(y) >= 12
    Rp_pin_center = y(12); % Operating pin center radius
else
    Rp_pin_center = Rp;
end

Np = Nc+1;
Ap = 2*pi/Np;

% F = (Cycloid_Profile_Projection) - (Pin_Surface_Projection)
% Cycloid uses Generating Params: (Rp+dRp), (Rrp+dRrp)
% Pin uses Operating Params: Rp_pin_center, Rrp_pin
% (Note: dRp applies to profile generation, not operating pin position)

F = @(x)(( (Rp+dRp)*cos(x) + (Rrp+dRrp)*cos(atan2(-sin(Nc*x),(cos(Nc*x)-(Rp+dRp)/(A*Np)))-x) - A*cos(Np*x) )*cos(Out)+...
    ( (Rp+dRp)*sin(x) - (Rrp+dRrp)*sin(atan2(-sin(Nc*x),(cos(Nc*x)-(Rp+dRp)/(A*Np)))-x) - A*sin(Np*x) )*sin(Out)+...
    A*cos(In) - (Rp_pin_center*cos(N*Ap) + Rrp_pin*cos(TP)));

% Initial guess calculation - simplified
iTC = atan2(Rp_pin_center*sin(N*Ap)+Rrp_pin*sin(TP)-A*sin(In), ...
    Rp_pin_center*cos(N*Ap)+Rrp_pin*cos(TP)-A*cos(In))+Out;

option = optimoptions('fsolve','MaxFunEvals',100,'TolFun',1e-10,'Display','off');
S = fsolve(F,(iTC),option);
F = double(S);
end
