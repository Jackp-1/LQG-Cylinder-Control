function BODE_plot(transfer,U_inf)

colors = {[0.062000000000000   0.258000000000000   0.501000000000000]};
w=logspace(-5, 1,30000);
[mag, phase, wout] = bode(transfer,w);
mag=squeeze(mag);
magdb = 20*log10(mag);
phase=squeeze(phase);
St = wout(:,1)/(U_inf*2*pi);
p=semilogx(St,magdb,'LineWidth',1);
ylabel('Magnitude (dB)','Interpreter','latex','FontSize',16,'Rotation',90)
xlim([10^-4 10])
xlabel('$\mathrm{St}$','Interpreter','latex','FontSize',16)
grid on
hold on
colororder("earth")
legend('11','12','13','14','15')
set(gcf,'Position',[0 0 900 300])