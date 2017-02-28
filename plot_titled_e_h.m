tilte_anlge_d = 1;

tilte_angle_r = tilte_anlge_d*pi/180;

EH_TOTAL_M = EH_TOTAL*cos(tilte_angle_r)-EV_TOTAL*sin(tilte_angle_r);

figure

plot(T,EH_TOTAL,'k','LineWidth',4)

hold on

plot(T,EH_TOTAL_M,'b','LineWidth',4)

make_label('Time (\mus)','Horizontal electric field (A/m)');