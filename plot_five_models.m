function []= plot_five_models(T_TL,I_TL,I2_TL,EV_TOTAL_TL,Hh_TOTAL_TL,Eh_TOTAL_TL,T_MTLL,I_MTLL,I2_MTLL,EV_TOTAL_MTLL,Hh_TOTAL_MTLL,Eh_TOTAL_MTLL,T_MTLE,I_MTLE,I2_MTLE,EV_TOTAL_MTLE,Hh_TOTAL_MTLE,Eh_TOTAL_MTLE,T_BG,I_BG,I2_BG,EV_TOTAL_BG,Hh_TOTAL_BG,Eh_TOTAL_BG,T_TCS,I_TCS,I2_TCS,EV_TOTAL_TCS,Hh_TOTAL_TCS,Eh_TOTAL_TCS)

figure
hold on
plot(T_TL,I_TL,'r');
plot(T_MTLL,I_MTLL,'b');
plot(T_MTLE,I_MTLE,'k');
plot(T_BG,I_BG,'g');
plot(T_TCS,I_TCS,'c');
make_label('Time (\mus)','Current (kA)');
legend('TL','MTLL','MTLE','BG','TCS');
hold off

figure 
hold on
plot(T_TL,I2_TL,'r');
plot(T_MTLL,I2_MTLL,'b');
plot(T_MTLE,I2_MTLE,'k');
plot(T_BG,I2_BG,'g');
plot(T_TCS,I2_TCS,'c');
make_label('Time (\mus)','Current (kA)');
legend('TL','MTLL','MTLE','BG','TCS');
hold off

figure
hold on
plot(T_TL,EV_TOTAL_TL,'r');
plot(T_MTLL,EV_TOTAL_MTLL,'b');
plot(T_MTLE,EV_TOTAL_MTLE,'k');
plot(T_BG,EV_TOTAL_BG,'g');
plot(T_TCS,EV_TOTAL_TCS,'c');
make_label('Time (\mus)','Vertical electric field (V/m');
legend('TL','MTLL','MTLE','BG','TCS');
hold off

figure
hold on
plot(T_TL,Hh_TOTAL_TL,'r');
plot(T_MTLL,Hh_TOTAL_MTLL,'b');
plot(T_MTLE,Hh_TOTAL_MTLE,'k');
plot(T_BG,Hh_TOTAL_BG,'g');
plot(T_TCS,Hh_TOTAL_TCS,'c');
make_label('Time (\mus)','Horizontal magnetic field (A/m)');
legend('TL','MTLL','MTLE','BG','TCS');
hold off

figure
hold on
plot(T_TL,Eh_TOTAL_TL,'r');
plot(T_MTLL,Eh_TOTAL_MTLL,'b');
plot(T_MTLE,Eh_TOTAL_MTLE,'k');
plot(T_BG,Eh_TOTAL_BG,'g');
plot(T_TCS,Eh_TOTAL_TCS,'c');
make_label('Time (\mus)','Horizontal electric field (V/m)');
legend('TL','MTLL','MTLE','BG','TCS');
hold off

end


