function [] = plot_currents_fields(T, I0_T, I1_T, I2_T, EV_TOTAL, HPHI_TOTAL, EH_TOTAL, T_carlos, I0_T_carlos, I1_T_carlos, I2_T_carlos, EV_TOTAL_carlos, HPHI_TOTAL_carlos, EH_TOTAL_carlos)

cl = 'r';

figure

hold on

plot(T,I0_T,cl,'LineWidth',2);
plot(T_carlos,I0_T_carlos,'b');

plot(T,I1_T,cl,'LineWidth',2);
plot(T_carlos,I1_T_carlos,'b');

plot(T,I2_T,cl,'LineWidth',2);
plot(T_carlos,I2_T_carlos,'b');

make_label('Time (\mus)','Current (kA)');

hold off

figure
hold on
plot(T, EV_TOTAL,cl,'LineWidth',2);
plot(T_carlos, EV_TOTAL_carlos,'b');

make_label('Time (\mus)','Vertical electric field (V/m)');

figure
hold on
plot(T, HPHI_TOTAL,cl,'LineWidth',2);
plot(T_carlos, HPHI_TOTAL_carlos,'b');

make_label('Time (\mus)','Horizontal magnetic field (A/m)');

figure
hold on
plot(T, EH_TOTAL,cl,'LineWidth',2);
plot(T_carlos, EH_TOTAL_carlos,'b');
make_label('Time (\mus)','Horizontal electric field (V/m)');


