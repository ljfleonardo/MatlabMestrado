
Qy = diag([1 1]);
Qu = diag(1);
K_1 = inv(G_PNMPC'*Qy*G_PNMPC+Qu)*G_PNMPC'*Qy;
K = K_1(1,:);

tic
for k = 2+delay_total:iteracoes

    %Compila vetor de saídas e perturbações
    xd = [saidas(:,k-1)' perturbacoes(:,k-1)']; 
    
    %Compila as entradas com atraso real
    entrada_atrasada_real = entradas(1,k-1-delay_real(1));
    
    %Compila as entradas com atraso modelado
    entrada_atrasada_mod  = entradas(1,k-1-delay_modelado(1));
    
    for t=1:N1
        for kk=1:N % N periodos de amostragem que a planta se estabiliza, N1 é o horizonte
            %---- Resposta Livre ----
            %f_1 = f(y(k+n-1|k),u(k+n-1))
            f_1 = modeloCSTR(saidas(:,t+kk-1),entradas(:,t+kk-1));
        end
        %f_2 = y(k)
        f_2 = modeloCSTR(saidas(:,t),entradas(:,t));
        %f_3 = f(y(k-1),u(k-1))
        f_3 = modeloCSTR(saidas(:,t-1),entradas(:,t-1));
        %y(k+n|k) = f(y(k+n-1|k),u(k+n-1))+y(k)-f(y(k-1),u(k-1))
        free = f_1(2,t)+f_2(2,t)-f_3(2,t);%saídas é apenas x2
                        %f(kk) = yp(k) + vect_g*u_livre'; %Saída livre= saída do processo + g*u_passado
    end
       
    %---- Função Custo ----
    delta_u = K*(ref_cb(k)-free); %deltaU = K(ref - respostaLivre)
    
    if k==1               %Não tenho o anterior, logo mantem o próprio delta_u
        entradas(:,k) = delta_u;
    else
        entradas(:,k) = entradas(:,k-1) + delta_u;
    end
    
    x_a_pred = [x_a_pred(1:2);x_a_estim(3)]'; %Compila estados e perturbações preditas em um vetor
    
    %---- Controle ----
    if controle == 1
        %Integrador para altura do tanque - h
        erro_1                = ref_cb(k) - x_a_pred(1);
        integrador_atual_1    = integrador_anterior_1 + erro_1;
        integrador_anterior_1 = integrador_atual_1;
        
        nova_entrada = -K*[x_a_pred(1);x_a_pred(2);x_a_pred(3); integrador_atual_1];
        entradas(:,k) = nova_entrada;
    end
    
    %---- Atualização das variáveis ----
    x_estimado_vect(:,k) = x_a_estim(1:end);
    x_pred_vect(:,k)     = x_a_pred';
%     P_estimado_at        = P_estimado;
    
end
Elapsed_time = toc;
disp('Tempo gasto na simulação:');
disp(Elapsed_time);

%---- Plot de gráficos ----
plot_cstr_FSP