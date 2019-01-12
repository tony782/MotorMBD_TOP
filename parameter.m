close all;
clear;

% Motor parameter
% Resistance (Ohm)
Ra = 0.1197;
% Inductance of d-axis (H)
Ld = 0.97e-3;
% Inductance of q-axis (H)
Lq = 2.03e-3;
% EMF constant (Wb)
Phia = 0.0045;
% Pole pairs number (-)
Pn = 4;
% Inertia of moment
Jm = 1;
% Damping coefficient
Dm = 1;
%���׃g���N
taul = 1;

% �V�~���[���[�V�������s
sim( 'MotorMBD_TOP' );