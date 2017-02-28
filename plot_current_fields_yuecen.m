function [] = plot_current_fields_yuecen(T,I0,I0_T,I1_T,I2_T,EV_TOTAL,HPHI_TOTAL,EH_TOTAL)

cl='r';

% figure
% plot(T,I0,cl);
% make_label('Time (\mus)','Undisturbed current (kA)');

% figure
% hold on
% plot(T,I0_T,cl);
% plot(T,I1_T,'b');
% plot(T,I2_T,'k');
% make_label('Time (\mus)','Current (kA)');
% hold off

figure

plot(T, EV_TOTAL,cl);

make_label('Time (\mus)','Vertical electric field (V/m)');

figure

plot(T, HPHI_TOTAL,cl);

make_label('Time (\mus)','Horizontal magnetic field (A/m)');

% figure
% 
% plot(T, EH_TOTAL,cl);
% 
% make_label('Time (\mus)','Horizontal electric field (V/m)');

end