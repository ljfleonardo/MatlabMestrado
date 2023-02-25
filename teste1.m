clc, clear all
%---------------Variáveis toolbox CasADi--------------------
import casadi.*
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x3 = SX.sym('x3');
x4 = SX.sym('x4');
u1 = SX.sym('u1');
u2 = SX.sym('u2');
d = SX.sym('d');

s = tf('s');

global Ts;
Ts = 0.01;

%---------------Sistema exemplo CSTR--------------------
%x1 = Ca; x2 = Cb; u = F/V; d = Caf

%Representação no contínuo
% dx1 = -x2*sin(x1)-(1/2)*x1+(2/3)*u;
% dx2 = (1/2)*x1*sin(x1)-x2-3*x1-(1/2)*x2;
delay1 = 2;
delay2 = 1.5;
delay3 = 0.5;
delay4 = 1;

dx1 = -x1+(1+(sin(x1))^2)*x2+delay1^2;
dx2 = x1*x2+x3+x4+((1+(sin(x1))^2)+0.5*(cos(x4))^2)*u1+delay2;
dx3 = -x3+x4+delay3;
dx4 = (x1+x3)*x4-x1*u1+((2+(sin(u1))^2)-sin(x3*x4-x1))*u2+delay1*delay4;


%Representação não-linear discreta (aproximação Euler)
% dx1_discreto = x1 + Ts*(-x2*sin(x1)-(1/2)*x1+(2/3)*u);
% dx2_discreto = x2 + Ts*((1/2)*x1*sin(x1)-x2-3*x1-(1/2)*x2);

dx1_discreto = x1 + Ts*(dx1);
dx2_discreto = x2 + Ts*(dx2);
dx3_discreto = x2 + Ts*(dx3);
dx4_discreto = x2 + Ts*(dx4);

fun_ax_ext = [dx1_discreto;dx2_discreto;dx3_discreto;dx4_discreto;d]; %Sistema aumentado discreto
fun_x = Function('fun_x',{x1,x2,x3,x4,d,u1,u2},{fun_ax_ext});

fun_yx_ext = x1;
fun_y = Function('fun_y',{x1,x2,x3,x4,d,u1,u2},{fun_yx_ext});

na = size(fun_ax_ext,1);                      %Dimensão vetor de estados aumentado
m  = size(fun_yx_ext,1);                      %Dimensão vetor de saídas
