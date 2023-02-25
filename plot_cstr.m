%% plot
h=figure;
subplot(3,1,1);
plot(saidas(2,:),'r','linewidth',2);
hold on
plot(x_pred_vect(2,:),'color','g','linestyle','-','linewidth',2);
% plot(x2_estimado_vect,'color','g','linestyle','-','linewidth',2);
% plot(x_pred_vect,'color','b','linestyle','-','linewidth',2);
plot(ref_cb,'color','k','linewidth',2);
ylabel('C_b [mol/l]','FontSize',12);
xlabel('Iterações [-]','FontSize',12);
if((delay_real-delay_modelado)~=0)
    title({'Sistema com erro de modelagem'; ['Delay  Real = ',num2str(delay_real), ' Delay Modelado = ',num2str(delay_modelado)]})
end
% title(['Ajustes de sintonia: Q = ',num2str(Q), ', R = ',num2str(R)],'FontSize',14);
%ylim([min(y_estimado_vect)-0.5 max(y_estimado_vect)+0.5]);
% legend({'Saída Real','Referência'},'Location','best','FontSize',12);
legend({'Saída Real','Saída Predita','Referência'},'Location','best','FontSize',12);
% legend({'Saída Real','Estimado','Predito','Referencia'},'Location','best','FontSize',12);
xlim([0 iteracoes]);
grid

subplot(3,1,2);
plot(caf','r','linewidth',2);
hold on
plot(x_estimado_vect(3,:),'g','linewidth',2);
ylabel('C_{af} [mol/l]', 'FontSize',12);
xlabel('Iterações [-]','FontSize',12);
% title('Perturbação','FontSize',14);
% ylim([min(d_estimado_vect)-0.5 max(d_estimado_vect)+0.5]);
legend({'Perturbação Real','Perturbação Estimada'},'Location','southeast','FontSize',12);
xlim([0 iteracoes]);
grid

subplot(3,1,3);
hold on
plot(entradas','k','linewidth',2);
ylabel('F [l/min]','FontSize',12);
xlabel('Iterações [-]','FontSize',12);
legend({'Sinal de Controle'},'Location','southeast','FontSize',12);
% title('Sinal de Controle','FontSize',14);
% ylim([min(entradas)-0.5 max(entradas)+0.5]);
grid
xlim([0 iteracoes]);
%}