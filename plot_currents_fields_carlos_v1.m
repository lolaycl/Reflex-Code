function [] = plot_currents_fields_carlos_v1(T, I0_T, I1_T, I2_T, EV_TOTAL, HPHI_TOTAL, EH_TOTAL)

cl = 'r';

figure

hold on

plot(T,I0_T,'k');

plot(T,I1_T,'r');

plot(T,I2_T,'b');

make_label('Time (\mus)','Current (kA)');

hold off

figure

plot(T, EV_TOTAL,cl);

make_label('Time (\mus)','Vertical electric field (V/m)');

figure

plot(T, HPHI_TOTAL,cl);

make_label('Time (\mus)','Horizontal magnetic field (A/m)');

figure

plot(T, EH_TOTAL,cl);

make_label('Time (\mus)','Horizontal electric field (V/m)');


