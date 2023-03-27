tamLetra = 10;
tamTitulo = 12;
espes = 3;

h = figure();
h.WindowState = 'maximized';

subplot(2,1,1);
plot(saidas(1,1:iteracoes),'r','linewidth',espes);
hold on
% plot(ref_h(1:iteracoes),'b--','linewidth',espes);
plot(x_pred_vect(1,1:iteracoes),'k--','linewidth',espes);
legend({'Real','Predição'},'FontSize',tamLetra);
% legend({'Real','Referência','Predição'},'FontSize',tamLetra);
xlim([0 iteracoes])
% ylim([0.92 1.2])
xlabel('Iterações (-)','FontSize',tamLetra);
ylabel('Ca (mol/l)','FontSize',tamLetra);
title('$C_a$ - Concentra\c{c}\~{a}o do produto A','interpreter','latex','FontSize',tamTitulo)
grid

subplot(2,1,2);
plot(saidas(2,1:iteracoes),'r','linewidth',espes);
hold on
% plot(ref_h(1:iteracoes),'b--','linewidth',espes);
plot(x_pred_vect(2,1:iteracoes),'k--','linewidth',espes);
legend({'Real','Predição'},'FontSize',tamLetra);
% legend({'Real','Referência','Predição'},'FontSize',tamLetra);
xlim([0 iteracoes])
% ylim([0.92 1.2])
xlabel('Iterações (-)','FontSize',tamLetra);
ylabel('C_b (mol/l)','FontSize',tamLetra);
title('$C_b$ - Concentra\c{c}\~{a}o do produto B','interpreter','latex','FontSize',tamTitulo)
grid