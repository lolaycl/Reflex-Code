function [] = write_infile(MPROBLEM,MODEL,I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B,V,W,H,ZLAM,ALPHA_TOWER,HR,RR,RG,HAB,HAM,ERR,...
    ERR2,ZS,SIG,EPSR,D1,D2,NSIMP,DELTATE1,TINT,DTS,TMAX,QEV,QH,QEH,KEY1,KEY2,F1,F2,F3,F4)

[str_model] = write_model(MPROBLEM,MODEL);

[str_current] = write_currents(I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B);

[str_channel] = write_speed_height_decay(V,W,H,ZLAM);

[str_reflection] = write_reflection(HR,RR,RG,ALPHA_TOWER);

[str_locations] = write_current_locations(HAB,HAM);

[str_errors] = write_errors(ERR,ERR2);

[str_medium] = write_medium(ZS,SIG,EPSR);

[str_distances] = write_distances(D1,D2,NSIMP);

[str_time_data] = write_time_data(DELTATE1,TINT,DTS,TMAX);

[str_scaling_factors] = write_scaling_factors(QEV,QH,QEH);

[str_keys] = write_keys(KEY1,KEY2);

[str_fs] = write_fs(F1,F2,F3,F4);

Infiletxt = strvcat(str_model,str_current,str_channel,str_reflection,str_locations,str_errors,str_medium,str_distances,...
    str_time_data,str_scaling_factors,str_keys,str_fs);

Infiletxt_size=size(Infiletxt);

fid2=fopen('infile.dat','wt+');

for k=1:Infiletxt_size(1,1)
    
    po=size(str2num(Infiletxt(k,:)));
    
    if (po == [0 0] & k~=1)
        
        fprintf(fid2,'\n');
        
        fprintf(fid2, '% s ',Infiletxt(k,:));
        
    else
        
        fprintf(fid2, '% s ',Infiletxt(k,:));
        
        fprintf(fid2,'\n');
        
    end
    
end

fclose(fid2);

clear Infiletxt;