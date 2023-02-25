clc, clear all, close all

%Dimensões
nx = 1;  %Sistema
nxa = 3; %Sistema aumentado, contendo w e v
nz = 1;  %Saída

%Central Difference Kalman Filter - uma das ab15ordagens do sigma point (implementação mais simples)
% h = sqrt(3);
% gama=h;
% Wmx(1) = (h^2-nxa)/(h^2);
% Wmx(2) = 1/(2*h*h);
% Wcx=Wmx;
% Wmxz = [Wmx(1) repmat(Wmx(2),[1 2*nxa])]';

%Unscented Kalman Filter - outra abordagem para sigma point
%[1] - Wan, E. A. an der Merwe, R. The unscented Kalman filter for nonlinear estimation

%alpha determines the spread of the sigma points around and is usually
%set to a small positive value (e.g., 1e-3) [1]

alpha = 1e-3; 

%lambda = alpha^2(L+k)-L onde L é a dimensão do sistema [1]
%k is a secondary scaling parameter which is usually set to 0 [1]
k=0;
lambda = (alpha^2)*(nx+k)-nx;% alfa^2(L+k)-L onde L é a dimensão do sistema

%beta = 2 para ruidos gaussianos [1]
beta = 2;
gamma = sqrt(nxa+lambda);

Wmx(1) = lambda/(nxa+lambda);
Wmx(2) = 1/(2*(nxa+lambda));
Wcx(1) = (lambda/(nxa+lambda))+(1-(alpha)^2+beta);
Wcx(2) = Wmx(2);
Wmxz   = [Wmx(1) repmat(Wmx(2),[1 2*nxa])]';

Ew = 1;
Ev = 2;

it = 40;
x_true_inicial = 2+randn(1); %Estado real
x_hat_inicial = 2;           %Estado estimado
Ex_inicial = 1;              %Covariancia
u_inicial = 0;               %Entrada atual desconhecida = 0;

x_store = zeros(it+1,length(x_true_inicial));
x_store(1,:) = x_true_inicial;
x_hat_store = zeros(it,length(x_hat_inicial));
Ex_store = zeros(it,length(x_hat_inicial)^2);

for (k = 1:it)
    
    %Passo 1.1 - Estimação de estados = xk = f_(k-1) (xhat_(k-1), u_(k-1),w\_(k-1))
    %a) calcular a matriz aumentada
    x_hat_a = [x_hat_inicial;0;0];   %Ruído do processo e do sensor
    
    %b) Fator de Cholesky
    Pxa = blkdiag(Ex_inicial,Ew,Ev); %Cria um bloco diagonal das covariancias
    sPxa = chol(Pxa,'lower');        %Matriz triangular inferior
    
    %c) calcular sigma points
    X = x_hat_a(:,ones([1 2*nxa+1])) + gamma*[zeros([nxa 1]), sPxa, -sPxa];
    
    %d) Calcular as equações de cada elemento
    Xx = sqrt(5+X(1,:))+X(2,:);
    x_hat_atual = Xx*Wmxz;
    
    %Passo 1.2 Estimação do erro de covariância
    Xs = ((Xx(:,2:end)) - x_hat_atual(:,ones([1 2*nxa])))*sqrt(Wcx(2));
    Xs1 = Xx(:,1) - x_hat_atual;
    Ex_atual = Xs*Xs' + Wcx(1)*Xs1*Xs1';
    
    %Definição das entradas desconhecidas e aleatórias
    w = chol(Ew)'*randn(1);   %Ruído aleatório
    v = chol(Ev)'*randn(1);   %Perturbação de saída alatória
    
    %Cálculo da saída e estados verdadeiros
    z_true_atual = x_true_inicial^3 + v;
    x_true_atual = sqrt(5+x_true_inicial) + w;
    
    %Passo 1.3 Estimação da saída
    Z = Xx.^3+X(3,:);%A saída do sistema é definida por xk^3+vk
    z_hat_atual = Z*Wmxz;
    
    %Passo 1.4 Ganho do filtro
    Zs = (Z(:,2:end) - z_hat_atual*ones([1 2*nxa])) * sqrt(Wcx(2));
    Zs1 = Z(:,1) - z_hat_atual;
    Exz = Xs*Zs' + Wcx(1)*Xs1*Zs1';
    Ez  = Zs*Zs' + Wcx(1)*Zs1*Zs1';
    L = Exz*inv(Ez);
    
    %Passo 1.5 Estimação do estado futuro
    x_hat_prox = x_hat_atual + L*(z_true_atual-z_hat_atual);
    
    %Passo 1.6 Estimação do erro de covariâncnia futuro
    Ex_prox = Ex_atual - L*Ez*L;
    
    x_store(k+1,:) = x_true_atual;
    x_hat_store(k,:) = x_hat_prox;
    Ex_store(k,:) = Ex_prox(:);
    
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
title('Sigma-point Kalman filter in action');


figure(2); 
plot(0:it-1,x_store(1:it)-x_hat_store,'b-',0:it-1, ...
    3*sqrt(Ex_store),'m--',0:it-1,-3*sqrt(Ex_store),'m--');
grid; 
legend('Error','Bounds',0);
title('SPKF Error with bounds');
xlabel('Iteration'); 
ylabel('Estimation Error');