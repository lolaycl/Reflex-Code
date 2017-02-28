function [str_distances] = write_distances(D1,D2,NSIMP)

str_txt1 = mat2str([D1 D2 NSIMP]);

str_distances = str_txt1(2:end-1);