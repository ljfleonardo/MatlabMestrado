clc, clear all, close all

%Exemplo ECE5550-Notes05 Universidade do colorado slide 15
SR_Ew = chol(1,'lower'); %Raíz quadrada da covariância da perturbação de saída
SR_Ev = chol(1,'lower'); %Raíz quadrada da covariância do ruído do sensor
A = 1;
B = 1;
C = 1;
D = 0;

it = 40;%Número de iterações desejadas

x_true = 0; %Cond. inicial
x_hat = 0;

Ex = 0.1;   %Covariância do filtro de kalman
SR_Ex = chol(Ex,'lower');
u = 0;    %Entrada inicial

xstore = zeros(it+1,length(x_true)); xstore(1,:) = x_true;
xhatstore = zeros(it,length(x_hat));
SR_Exstore = zeros(it,length(x_hat));

for(k = 1:it)
    
    %Estimação
    x_hat = A*x_hat+B*u;              %Estado estimado
    SR_Ex = qr([A*SR_Ex, SR_Ew]');    %Erro da covariancia
    %Decomposição QR -> Q ortogonal e R triangular superior (relacionado com cholesky)
    SR_Ex = tril(SR_Ex(1:length(x_hat),1:length(x_hat)));%Triangular inferior
    
    u = 0.5*randn(1);                  %Entrada aleatória
    w = SR_Ew*randn(length(x_true));   %Ruído aleatório
    v = SR_Ev*randn(length(C*x_true)); %Perturbação de saída alatória
        
    z_true = C*x_true + D*u + v; %z baseado em x e u
    x_true = A*x_true + B*u + w; %próximo x baseado na entrada atual
    
    %Atualização
    z_hat  = C*x_hat + D*u;             %Estimação da saída
    
    SR_Ez = qr([C*SR_Ex, SR_Ev]');
    SR_Ez = tril(SR_Ez(1:length(z_hat),1:length(z_hat)));
    
    L      = (SR_Ex*SR_Ex')*C'*inv(SR_Ez')*inv(SR_Ez);     %Ganho do filtro de Kalman
    x_hat  = x_hat + L*(z_true-z_hat);  %Atualização da medição da estimativa do estado
    
    %Como Sx = Sx-L*Sz*Sz'*L', não é possível usar QR por causa do sinal -
    Sx_ = SR_Ex';
    cov_update_vectors = L*SR_Ez;
    
    for j=1:length(z_hat),
        Sx_ = cholupdate(Sx_,cov_update_vectors(:,j),'-');
    end
    SR_Ex = Sx_';
    xstore(k+1,:) = x_true; 
    xhatstore(k,:) = x_hat;
    SR_Exstore(k,:) = diag(SR_Ex*SR_Ex');
    
end

figure(1);
plot(0:it-1,xstore(1:it)','k-',...                  %Real
    0:it-1,xhatstore','b--', ...                    %Estimado
    0:it-1,xhatstore'+3*sqrt(SR_Exstore'),'m-.',... %Limite superior
    0:it-1,xhatstore'-3*sqrt(SR_Exstore'),'m-.');   %Limite inferior
grid;

legend('True','Estimate','Bounds');
title('Kalman filter in action');
xlabel('Iteration');
ylabel('State');

figure(2);
clf;
plot(0:it-1,xstore(1:it)'-xhatstore','b-',...    %Erro
    0:it-1,3*sqrt(SR_Exstore'),'m--',...         %Limite superior
    0:it-1,-3*sqrt(SR_Exstore'),'m--');          %Limite inferior
grid;
legend('Error','Bounds',0);
title('Error w/ bounds');
xlabel('Iteration');
ylabel('Estimation Error');

