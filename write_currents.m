function [str_current] = write_currents(I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B)

str1 = mat2str([I01 ALPHA BETA]);

str2 = mat2str([I02 GAMMA DELTA]);

str3 = mat2str([IW1 TAU11 TAU21 NW1]);

str4 = mat2str([IW2 TAU12 TAU22 NW2]);

str5 = mat2str([I0B T0B]);

str_current = strvcat(str1(2:end-1),str2(2:end-1),str3(2:end-1),str4(2:end-1),str5(2:end-1));