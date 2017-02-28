function [str_fs] = write_fs(F1,F2,F3,F4)

str_txt1 = mat2str([F1]);

str_txt2 = mat2str([F2]);

str_txt3 = mat2str([F3]);

str_txt4 = mat2str([F4]);

str_fs = strvcat(str_txt1,str_txt2,str_txt3,str_txt4);