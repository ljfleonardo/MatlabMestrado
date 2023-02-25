%Preditor de Kalman
%Item 3.1 (Lima, B. M. Lima, D. M. Normey-Rico, J. E. A robust predictor for dead-time systems based on the Kalman filter. IFAC, 2018)
clc, clear all

%---------------Definição do sistema--------------
s = tf('s');

Ts = 0.2; %Dado no trabalho
% Ts = 0.1;
Y = (3*s+1)/((s+1)*((2*s+1)^2));
Y.inputdelay = 5;
% Y = 1/(3*s+1)*exp(-2*s);
Y_d = c2d(Y,Ts);
Y_dss = ss(Y_d);

temp = Y_d;
temp.inputdelay = 0;

%--------------- Todas as variâncias = 1 ---------

% Better disturbance rejection: if you wish to improve
% the disturbance rejection, make Qs = I, Qd > 1
% an R < 1. An interesting initial choice is R = 1
% and/or Qd = 10. 

Q  = 2;

R  = 1; %pode ir de 0.01 a 100^4, depende de cada caso, influencia no sinal de controle

%---------------Controlador--------------------
z = tf('z',Ts);

% C = 4.3423*(((s+1)*(s+0.5)^2)/(s*(s^2+0.49)));
% Cz = c2d(C,Ts,'tustin');
Cz = (12.479*(z-0.7531)^2)/(z*(z-1));
% t = feedback(Cz*Fz*Y_d,1);       %MF do sistema 

iteracoes = 40;

% [fkalman,Fr2,KalObs] = calcKalmanPredictor(Y_dss,Q,[]);%Filtro de erro de predição

[Fr,L,P] = kalman_leonardo(Y_dss,Q,R,iteracoes);

% Fr=Fr2;

Sr    = temp-Y_d*Fr;

t = feedback(Cz*temp,1);
t.inputdelay = 20;
y_q = temp*(1-t*Fr); %TF perturbação para saída
y_n = (1-t*Fr);      %TF ruído para saída


figure
step(t,y_q,y_n)
legend('saída','perturbação/saída','ruído/saída')
title({['R is ',num2str(R)];['Q is ', num2str(Q)]})