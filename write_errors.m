function [str_errors] = write_errors(ERR,ERR2)

str_txt1 = mat2str([ERR ERR2]);

str_errors = str_txt1(2:end-1);