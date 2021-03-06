close all;
clear;

%変換式
wrm2rpm = 60/(2*pi);
rpm2wrm = (2*pi)/60;

%定格
fpwm = 2e+04;

% Motor parameter
% Resistance (Ohm)
Ra = 0.1197;
% Inductance of d-axis (H)
Ld = 1.09e-4;
% Inductance of q-axis (H)
Lq = 1.19e-4;
% EMF constant (Wb)
Phia = 4.5e-3;
% Pole pairs number (-)
Pn = 4;
% Inertia of moment
Jm = 0.0001;
% Damping coefficient
Dm = 1;
%負荷トルク
taul = 0;

%電流制御
wc = 2*pi*fpwm/3/5;
Kdp = Ld*wc;
Kdi = Ra*wc;
Kqp = Lq*wc;
Kqi = Ra*wc;
% Kdp = 1;
% Kdi = 3;
% Kqp = 1;
% Kqi = 3;

%速度制御
wsc = wc/10;
Ksp = Jm*wsc;
Ksi = Jm*wsc*(wsc/5);
Ksa = 3/Ksp;
