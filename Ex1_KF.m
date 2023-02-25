clc, clear all, close all

%Exemplo ECE5550-Notes04 Universidade do colorado slide 28

A = 1;
B = 1;
C = 1;
D = 0;

Ew = 1;   %Covariancia do ruído do processo
Ev = 0.1; %Covariancia da perturbação de saída = ruido de medição

%Sem incertezas iniciais
it = 40;%Número de iterações desejadas

x_true = 0; %Cond. inicial
Ex = 0;     %Variância da cond. inicial

x_hat = 0;


%Assumindo as condições inicias:
u = 0;    %Entrada inicial

xstore = zeros(length(x_true),it+1); xstore(:,1) = x_hat;
xhatstore = zeros(length(x_hat),it);
Exstore = zeros(length(x_hat)^2,it);

for(k = 1:it)
    
    %Estimação
    x_hat = A*x_hat+B*u;   %Estado estimado
    Ex    = A*Ex*A'+Ew;    %Erro da covariancia
    
    u = 0.5*randn(1);                     %Entrada aleatória
    w = chol(Ew)'*randn(length(x_hat));   %Ruído aleatório
    v = chol(Ev)'*randn(length(C*x_hat)); %Perturbação de saída alatória
    
    %Chol é a decomposição de Cholesky que pega a matriz triangular superior
    
    z_true = C*x_true + D*u + v; %z baseado em x e u
    x_true = A*x_true + B*u + w; %próximo x baseado na entrada atual
    
    %Atualização
    z_hat  = C*x_hat + D*u;             %Estimação da saída
    L      = Ex*C'*inv(C*Ex*C'+Ev);     %Ganho do filtro de Kalman
    x_hat  = x_hat + L*(z_true-z_hat);  %Atualização da medição da estimativa do estado
    Ex     = (eye(length(A))-L*C)*Ex;   %Atualização da medição do erro da covariância
    
    xstore(:,k+1) = x_true; xhatstore(:,k) = x_hat;
    Exstore(:,k) = Ex(:);
    
end

figure(1);
plot(0:it-1,xstore(1:it)','k-',...               %Real
    0:it-1,xhatstore','b--', ...                 %Estimado
    0:it-1,xhatstore'+3*sqrt(Exstore'),'m-.',... %Limite superior
    0:it-1,xhatstore'-3*sqrt(Exstore'),'m-.');   %Limite inferior
grid;

legend('true','estimate','bounds');
title('Kalman filter in action'); 
xlabel('Iteration'); 
ylabel('State');

figure(2); 
clf;
plot(0:it-1,xstore(1:it)'-xhatstore','b-',... %Erro
    0:it-1,3*sqrt(Exstore'),'m--',...         %Limite superior
    0:it-1,-3*sqrt(Exstore'),'m--');          %Limite inferior
grid;
legend('Error','Bounds',0);
title('Error w/ bounds'); 
xlabel('Iteration'); 
ylabel('Estimation Error');

