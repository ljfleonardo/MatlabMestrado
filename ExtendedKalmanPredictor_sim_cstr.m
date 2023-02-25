%Preditor de Kalman - Extended Kalman
clc, clear all
%---------------Variáveis toolbox CasADi--------------------
import casadi.*
x1 = SX.sym('x1');
x2 = SX.sym('x2');
u = SX.sym('u');
d = SX.sym('d');

s = tf('s');

global k1;
global k2;
global k3;
global Ts;

k1 = 6.01;
k2 = 0.8433;
k3 = 0.1123;
% Ts = 0.01;
Ts = 0.02;

%---------------Sistema exemplo CSTR--------------------
%x1 = Ca; x2 = Cb; u = F/V; d = Caf

%Representação no contínuo
dx1 = -k1*x1 - k3*x1^2 -x1*u;
dx2 = k1*x1 - k2*x2 - x2*u;

%Representação não-linear discreta
dx1_discreto = x1 + Ts*(-k1*x1-k3*x1^2+(d-x1)*u);
dx2_discreto = x2 + Ts*(k1*x1-k2*x2-x2*u);

fun_ax_ext = [dx1_discreto;dx2_discreto;d]; %Sistema aumentado discreto
fun_x = Function('fun_x',{x1,x2,d,u},{fun_ax_ext});

fun_yx_ext = x2;
fun_y = Function('fun_y',{x1,x2,d,u},{fun_yx_ext});
na = size(fun_ax_ext,1);                      %Dimensão vetor de estados aumentado
m  = size(fun_yx_ext,1);                      %Dimensão vetor de saídas
%% !!!!!!!!!!!!!! ATRASO !!!!!!!!!!!!!!!!!!
delay_real = 3;
delay_modelado = 4;

%% ---------------Inicialização das variáveis--------------------
iteracoes   = 650;

saidas      = zeros(2,iteracoes);            %inicia vetor de saídas
entradas    = zeros(1,iteracoes);            %inicia vetor de entradas
ref_cb      = zeros(1,iteracoes);
saidas(1,:) = 0.7192;
saidas(2,:) = 2.3450;
caf         = zeros(1,iteracoes);            %inicia vetor de perturbação
ruido = zeros(1,iteracoes);

entradas(1:end)   = 1;                       %inicia vetor de entradas com 1
ref_cb(1:end)     = 2.3450;
ref_cb(200:end)   = 1.6;
caf(1:end)        = 5.1;                     %inicia vetor de perturbação com 5.1
caf(400:end)      = 5.3;                     %degrau na perturbação de 0.2 a partir de 400 it

ca_estimado = 0.7192;                        %inicia ca estimado próximo do ponto de operação
cb_estimado = 2.3450;                        %inicia cb estimado próximo do ponto de operação
x_estimado  = [ca_estimado;cb_estimado];     %cria vetor de variáveis estimadas

d_estimado = caf(1);                         %inicia variável da perturbação estimada
y_estimado = saidas(2,1);                    %inicia variável da saída estimada
P_estimado    = 0;                           %inicia variável de erro de covariância
P_estimado_at = eye(3);

d_estimado_vect = zeros(1,iteracoes);        %inicia vetor de pert. estimada (usado apenas para plotar)
d_estimado_vect(:) = d_estimado;
y_estimado_vect = zeros(1,iteracoes);        %inicia vetor de saída estimada (usado apenas para plotar)
y_estimado_vect(:) = y_estimado;
x_estimado_vect     = zeros(na,iteracoes);
x_estimado_vect(1,:) = ca_estimado;
x_estimado_vect(2,:) = cb_estimado;
x_estimado_vect(3,:) = d_estimado;

x_a_pred = zeros(na-1,iteracoes);
x_pred_vect = zeros(na-1,iteracoes);
x_pred_vect(1,:)      = ca_estimado;
x_pred_vect(2,:)      = cb_estimado;

x_a_estim = [x_estimado;d_estimado];

A_jacobian = jacobian(fun_ax_ext,[x1 x2 d]); %cálculo do jacobiano para matriz A
B_jacobian = jacobian(fun_ax_ext,u);         %cálculo do jacobiano para matriz B
C_jacobian = jacobian(fun_yx_ext,[x1 x2 d]); %cálculo do jacobiano para matriz C

Q  = 1;                                      %inicia variável da ponderação dos estados
R  = 1e-2;                                    %inicia variável da ponderação da saída
% R = 1e-3;

%Utilizando CasADi (otimização para diferenciação), cria função para realizar jacobiano
fun_a = Function('fun_a',{x1,x2,d,u},{A_jacobian});
fun_b = Function('fun_b',{x1,x2,d,u},{B_jacobian});
fun_c = Function('fun_a',{x1,x2,d,u},{C_jacobian});

%% ---------------Controle---------------
for controle=1
Aa = full(fun_a(saidas(1,1),saidas(2,1),d_estimado,entradas(1,1))); %Transforma para double
Ba = full(fun_b(saidas(1,1),saidas(2,1),d_estimado,entradas(1,1)));
Ca = [0 1];
Ba_x = [Ba(1);Ba(2)];
Aa_x = [Aa(1,1) Aa(1,2); Aa(2,1) Aa(2,2)];

A_exp = [Aa_x zeros(size(Aa_x,1),1);
    -Ca 1];
B_exp = [Ba_x;0];
C_exp = [Ca 0];
D_exp = 0;
sistema = ss(A_exp,B_exp,C_exp,D_exp,Ts);                            %Sistema expandido

Co = ctrb(sistema.A,sistema.B);                                      %Matriz de Controlabilidade
rank(Co);
Q_lqr = diag([0,10,1]);                                              %Matrizes de ponderação
R_lqr = 10*eye(size(sistema.B,2));
[K,P,eig] = dlqr(sistema.A,sistema.B,Q_lqr,R_lqr);                   %Ganho do controlador via dlqr

%inicializa integrador no ponto de operação
integrador_anterior = (entradas(1)+K(1)*ca_estimado+K(2)*cb_estimado)/(-K(3));
integrador_atual = 0;
end

%% ---------------Simulação--------------------
tic
for k = 2+delay_real+delay_modelado:iteracoes
%     if(k>450)
%     ruido(k) = 0.05*randn;
%     end
    %------simulação do processo------
    saidas(:,k) = modeloCSTR(saidas(1,k-1),saidas(2,k-1),caf(k-1),(entradas(k-1-delay_real)+ruido(k)));
    
    Aa = full(fun_a(x_a_estim(1),x_a_estim(2),...
        x_a_estim(3),entradas(k-1-delay_modelado)));                %substitui os valores reais em A_jacobian
    Ba = full(fun_b(x_a_estim(1),x_a_estim(2),...
        x_a_estim(3),entradas(k-1-delay_modelado)));                %substitui os valores reais em B_jacobian
    Ca = full(fun_c(x_a_estim(1),x_a_estim(2),...
        x_a_estim(3),entradas(k-1-delay_modelado)));                %substitui os valores reais em C_jacobian
    
    %------estimação------

    [x_a_estim,P_estimado] = ekf(fun_x,fun_y,saidas(2,k),...
    entradas(k-1-delay_modelado),x_a_estim,P_estimado_at,Aa,Ca,Q,R);

%     x_a_estim = x_estimado;
    %------predição------
    %primeira iteração
    x_a_pred = modeloCSTR(x_a_estim(1),x_a_estim(2),x_a_estim(3),entradas(k-delay_modelado-1));
    if(k<=(k-delay_modelado+j-1))
        %segunda iteração em diante:
        for j=2:delay_modelado
            x_a_pred = modeloCSTR(x_a_pred(1),x_a_pred(2),x_a_estim(3),entradas(k-delay_modelado+j-1));
        end
    end
    
    %------controle------
    erro = ref_cb(k)-x_a_pred(2);
    integrador_atual = integrador_anterior+erro;
    
    integrador_anterior = integrador_atual;
    
    nova_entrada = -K*[x_a_pred(1);x_a_pred(2);integrador_atual];
   
    
    % saturação
    if(nova_entrada<0)
        nova_entrada = 0;
    elseif(nova_entrada>10)
        nova_entrada = 10;
    end
    
    entradas(k) = nova_entrada;
    
    %------atualização das variáveis------
    x_estimado_vect(:,k) = x_a_estim(1:end);
    x_pred_vect(:,k) = x_a_pred(1:end);
    P_estimado_at = P_estimado;
%     y_estimado_vect(k) = y_estimado;
end
Elapsed_time = toc;
disp('Tempo gasto na simulação:');
disp(Elapsed_time);
%% ---------------plot---------------
plot_cstr; %função para plotar os gráficos