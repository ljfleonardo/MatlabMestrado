%% plot
figure

subplot(3,2,1);
plot(saidas(1,:),'linewidth',2)
hold on
plot(x_pred_vect(1,:),'linewidth',2)
legend('Saída Real', 'Saída Predita')
title('X_1');
xlim([0 iteracoes]);
grid

subplot(3,2,2);
plot(saidas(2,:),'linewidth',2)
hold on
plot(x_pred_vect(2,:),'linewidth',2)
legend('Saída Real', 'Saída Predita')
title('X_2');
xlim([0 iteracoes]);
grid

subplot(3,2,3);
plot(saidas(3,:),'linewidth',2)
hold on
plot(x_pred_vect(3,:),'linewidth',2)
legend('Saída Real', 'Saída Predita')
title('X_3');
xlim([0 iteracoes]);
grid

subplot(3,2,4);
plot(pert,'linewidth',2)
hold on
plot(x_estimado_vect(3,:),'linewidth',2)
legend('Perturbação Real', 'Perturbação Predita')
title('Perturbação');
xlim([0 iteracoes]);
grid

% figure
subplot(3,2,[5,6]);
plot(entradas,'k','linewidth',2);
title('Sinal de controle');
xlim([0 iteracoes]);
grid