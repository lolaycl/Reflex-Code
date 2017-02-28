function [str_locations] = write_current_locations(HAB,HAM)

str_txt1 = mat2str([HAB HAM]);

str_locations = str_txt1(2:end-1);