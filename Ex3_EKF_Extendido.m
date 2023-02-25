clc, clear all, close all

%Exemplo ECE5550-Notes06 Universidade do colorado slide 8
Ew = 1;
Ev = 2;
% syms xk wk vk
%
% fk = sqrt(5+xk)+wk;
% hk = xk^3+vk;
%
% A_hat = diff(fk,xk);
% B_hat = diff(fk,wk);
% C_hat = diff(hk,xk);
% D_hat = diff(hk,vk);

it = 40;%N�mero de itera��es desejadas

x_true_inicial = 2+randn(1); %Estado real
x_hat_inicial = 2;           %Estado estimado
Ex_inicial = 1;              %Covariancia

x_store = zeros(it+1,length(x_true_inicial)); 
x_store(1,:) = x_true_inicial;
x_hat_store = zeros(it,length(x_hat_inicial));
Ex_store = zeros(it,length(x_hat_inicial)^2);

for (k = 1:it)
    A_hat = (1/2)*(sqrt(5+x_hat_inicial));
    B_hat = 1;
    
    %Passo 1.1 - Estima��o de estados = xk = f_(k-1) (xhat_(k-1), u_(k-1),w\_(k-1))
    x_hat_atual = sqrt(5+x_hat_inicial); %Sistema com wk=0
    
    %Passo 1.2 - Estima��o do erro de covari�ncia = Ex
    Ex_atual = A_hat*Ex_inicial*A_hat' + B_hat*Ew*B_hat';
    
    %Defini��o das entradas desconhecidas e aleat�rias
    w = chol(Ew)'*randn(1);   %Ru�do aleat�rio
    v = chol(Ev)'*randn(1);   %Perturba��o de sa�da alat�ria
    
    %C�lculo da sa�da e estados verdadeiros
    z_true_atual = x_true_inicial^3 + v;
    x_true_atual = sqrt(5+x_true_inicial) + w;
    
    C_hat = 3*x_hat_atual^2+v;
    D_hat = 1;
    
    %Passo 1.3 - Estima��o de sa�da zk = hk(xhat_k,u_k,v\_k)
    z_hat_atual = x_hat_atual^3 + v;
    
    %Passo 1.4 - Ganho do filtro L = Ex*C'*inv(C*Ex*C' + D*Ev*D')
    L = Ex_atual*C_hat'*inv(C_hat*Ex_atual*C_hat' + D_hat*Ev*D_hat');
    
    %Passo 1.5 - Estimar estado futuro (atualiza��o) x_hat_+ = x_hat_- + L*(z_k - z_hat_-)
    x_hat_prox = x_hat_atual + L*(z_true_atual-z_hat_atual);
    x_hat_prox = max(-5,x_hat_prox);
    
    %Passo 1.6 - Estimar erro de covari�ncia (atualiza��o) Ex = (I-L*C)*Ex
    Ex_prox = (eye(length(L))-L*C_hat)*Ex_atual;
    
    
    x_store(k+1,:) = x_true_atual;
    x_hat_store(k,:) = x_hat_prox;
    Ex_store(k,:) = Ex_prox(:);
    
    %Atualiza��o das vari�veis
    x_true_inicial = x_true_atual;
    x_hat_inicial = x_hat_prox;
    Ex_inicial = Ex_prox;
    
end

figure(1); 
plot(0:it-1,x_store(1:it),'k-',0:it-1,x_hat_store,'b--', ...
    0:it-1,x_hat_store+3*sqrt(Ex_store),'m-.',...
    0:it-1,x_hat_store-3*sqrt(Ex_store),'m-.'); 
grid;
legend('True','Estimate','Bounds'); 
xlabel('Iteration'); 
ylabel('State');
title('Extended Kalman filter in action');


figure(2); 
plot(0:it-1,x_store(1:it)-x_hat_store,'b-',0:it-1, ...
    3*sqrt(Ex_store),'m--',0:it-1,-3*sqrt(Ex_store),'m--');
grid; 
legend('Error','Bounds',0);
title('EKF Error with bounds');
xlabel('Iteration'); 
ylabel('Estimation Error');