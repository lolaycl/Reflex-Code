function [str_scaling_factors] = write_scaling_factors(QEV,QH,QEH)

str_txt1 = mat2str([QEV QH QEH]);

str_scaling_factors = str_txt1(2:end-1);