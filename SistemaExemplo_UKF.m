%Preditor de Kalman - Unscented Kalman
clc, clear all
%---------------Variáveis toolbox CasADi--------------------
import casadi.*
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x3 = SX.sym('x3');
u = SX.sym('u');
d = SX.sym('d');

global Ts

Ts = 0.02;

dx1_discreto = x1 + Ts*(-x1+x2+d);
dx2_discreto = x2 + Ts*(x1-x2-x1*x3+u);
dx3_discreto = x3 + Ts*(x1+x1*x2-2*x3);


fun_ax_ext = [dx1_discreto;dx2_discreto;dx3_discreto;d]; %Sistema aumentado discreto
fun_x = Function('fun_x',{x1,x2,x3,d,u},{fun_ax_ext});

fun_yx_ext = x1;
fun_y = Function('fun_y',{x1,x2,x3,d,u},{fun_yx_ext});

% A_jacobian = jacobian(fun_ax_ext,[x1 x2 x3]); %cálculo do jacobiano para matriz A
% B_jacobian = jacobian(fun_ax_ext,u);         %cálculo do jacobiano para matriz B
% C_jacobian = jacobian(fun_yx_ext,[x1 x2 x3]); %cálculo do jacobiano para matriz C

na = size(fun_ax_ext,1);                      %Dimensão vetor de estados aumentado
m  = size(fun_yx_ext,1);                      %Dimensão vetor de saídas
%% ---------------Inicialização das variáveis---------------
iteracoes   = 850;

saidas      = zeros(3,iteracoes);            %inicia vetor de saídas
entradas    = zeros(1,iteracoes);            %inicia vetor de entradas
ref_cb      = zeros(1,iteracoes);
pert        = zeros(1,iteracoes);

saidas(1,:) = 0.998;
saidas(2,:) = 1.002;
saidas(3,:) = 1;

entradas (:) = 1;
entradas(100:end) = 2;
pert(400:end) = 1;

x1_estimado = saidas(1,1);                   %inicia ca estimado próximo do ponto de operação
x2_estimado = saidas(2,1);                   %inicia cb estimado próximo do ponto de operação
x3_estimado = saidas(3,1);
x_estimado  = [x1_estimado;x2_estimado;x3_estimado];     %cria vetor de variáveis estimadas

d_estimado    = pert(1);                         %inicia variável da perturbação estimada
y_estimado    = saidas(3,1);                    %inicia variável da saída estimada
P_estimado    = 0;                           %inicia variável de erro de covariância
P_estimado_at = 0.1*eye(na);

d_estimado_vect       = zeros(1,iteracoes);        %inicia vetor de pert. estimada (usado apenas para plotar)
d_estimado_vect(:)    = d_estimado;
y_estimado_vect       = zeros(1,iteracoes);        %inicia vetor de saída estimada (usado apenas para plotar)
y_estimado_vect(:)    = y_estimado;
x_pred_vect           = zeros(3,iteracoes);
x_estimado_vect       = zeros(na,iteracoes);
x_estimado_vect(1,:)  = x1_estimado;
x_estimado_vect(2,:)  = x2_estimado;
x_estimado_vect(3,:)  = x3_estimado;
x_pred_vect(1,:)      = x1_estimado;
x_pred_vect(2,:)      = x2_estimado;
x_pred_vect(3,:)      = x3_estimado;

x_a_estim = [x_estimado;d_estimado];
x_a_pred  = zeros(3,iteracoes);

% fun_a = Function('fun_a',{x1,x2,x3,u},{A_jacobian});
% fun_b = Function('fun_b',{x1,x2,x3,u},{B_jacobian});
% fun_c = Function('fun_a',{x1,x2,x3,u},{C_jacobian});

% Aa = full(fun_a(saidas(1,1),saidas(2,1),saidas(3,1),entradas(1,1))); %Transforma para double
% Ba = full(fun_b(saidas(1,1),saidas(2,1),saidas(3,1),entradas(1,1)));
% Ca = [1 0 0];


% Q = 0.1^2*eye(na-1);                               %inicia variável da ponderação dos estados
Q = diag(0.1*[ones(1,na-1),100]);
% Q = diag([1,1,1,100]);
% Q = 10*eye(na);
R = diag(ones(1,m));
% R = 0.1^2*eye(m);                                  %inicia variável da ponderação da saída
% R = 1e2*eye(m);
P_estimado_at = 1.1*Q;

%% !!!!!!!!!!!!!! ATRASO !!!!!!!!!!!!!!!!!!
delay_real     = 0;
delay_modelado = 0;
%% ---------------Simulação---------------

tic
for k = 2+delay_real:iteracoes
%     ruido(k) = 0.05*randn;
    %------Simulação do processo------
    saidas(:,k) = full(modelo_exemplo(saidas(1,k-1),saidas(2,k-1),saidas(3,k-1),pert(k-1),(entradas(k-1-delay_real))));
    
    %------Estimação------
    [x_estimado,P_estimado] = ukf(fun_x,fun_y,saidas(1,k),...
        entradas(k-1-delay_modelado),x_a_estim,P_estimado_at,Q,R,na,m);
    
    x_a_estim = x_estimado;
%     ------Predição------
    %primeira iteração
    x_a_pred = modelo_exemplo(x_a_estim(1),x_a_estim(2),x_a_estim(3),x_a_estim(4),entradas(k-delay_modelado-1));
    if(k<=(k-delay_modelado+j-1))
        %segunda iteração em diante:
        for j=2:delay_modelado
            x_a_pred = modelo_exemplo(x_a_pred(1),x_a_pred(2),x_a_pred(3),x_a_estim(4),entradas(k-delay_modelado+j-1));
        end
    end

    
    %------Controle------
%     erro = ref_cb(k)-x_a_pred(2);
%     integrador_atual = integrador_anterior+erro;
%     
%     integrador_anterior = integrador_atual;
%     
%     nova_entrada = -K*[x_a_pred(1);x_a_pred(2);integrador_atual];
    
%     %------Saturação------
%     if(nova_entrada<0)
%         nova_entrada = 0;
%     elseif(nova_entrada>10)
%         nova_entrada = 10;
%     end
%     entradas(k) = nova_entrada;
%     
    %------Atualização das variáveis------
%     x1_estimado_vect(k) = x_a_estim(1);
%     x2_estimado_vect(k) = x_a_estim(2);
    x_estimado_vect(:,k) = x_a_estim(1:end);
    x_pred_vect(:,k) = x_a_pred;
    y_estimado_vect(k) = y_estimado;
%     d_estimado_vect(k) = d_estimado;
    P_estimado_at = P_estimado;
end

Elapsed_time = toc;
disp('Tempo gasto na simulação:');
disp(Elapsed_time);
%% plot
plot_exemplo
