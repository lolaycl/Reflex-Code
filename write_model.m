function [str_model] = write_model(MPROBLEM,MODEL)

str_txt1 = mat2str([MODEL MPROBLEM]);

str_model = str_txt1(2:end-1);