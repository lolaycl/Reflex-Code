function [str_keys] = write_keys(KEY1,KEY2)

str_txt1 = mat2str([KEY1 KEY2]);

str_keys = str_txt1(2:end-1);