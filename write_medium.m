function [str_medium] = write_medium(ZS,SIG,EPSR)

str_txt1 = mat2str([ZS SIG EPSR]);

str_medium = str_txt1(2:end-1);