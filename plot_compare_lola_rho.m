function []=plot_compare_lola_rho(T,EV1,EV2,EV3,H1,H2,H3,EH1,EH2,EH3)

figure
hold on
plot(T,H1,'r');
plot(T,H2,'b');
plot(T,H3,'k');
make_label('Time (\mus)','Horizontal magnetic field (A/m)');
legend('rho_g=0.5','rho_g=0.8','rho_g=1');
hold off

figure
hlod on
plot(T,EV1,'r');
plot(T,EV2,'b');
plot(T,EV3,'k');
make_label('Time (\mus)','Vertical electric field (V/m');
legend('rho_g=0.5','rho_g=0.8','rho_g=1');
hold off

figure
hold on
plot(T,EH1,'r');
plot(T,EH2,'b');
plot(T,EH3,'k');
make_label('Time (\mus)','Horizontal electric field (V/m)');
legend('rho_g=0.5','rho_g=0.8','rho_g=1');
hold off

end

