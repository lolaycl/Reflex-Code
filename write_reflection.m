function [str_reflection] = write_reflection(HR,RR,RG,ALPHA_TOWER)

str_txt1 = mat2str([HR RR RG ALPHA_TOWER]);

str_reflection = str_txt1(2:end-1);
