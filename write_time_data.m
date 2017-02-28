function [str_time_data] = write_time_data(DELTATE1,TINT,DTS,TMAX)

str_txt1 = mat2str([DELTATE1 TINT DTS TMAX]);

str_time_data = str_txt1(2:end-1);