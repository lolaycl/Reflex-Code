function [str_channel] = write_speed_height_decay(V,W,H,ZLAM)

str1 = mat2str([V W H ZLAM]);

str_channel = str1(2:end-1);