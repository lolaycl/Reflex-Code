function [ezhr] = obtain_ez_eh_ratio(T,EV,EH)

td = 1;

Nt = max(size(T));

EH = abs(EH);

EV = abs(EV);

jj = 1;

while ((jj > 0) && (jj < Nt))
    
    if (EH(jj,1) > EH(jj+1,1))
        
        break;
        
    end
    
    jj = jj + 1;
    
end

ih_max = jj;

eh_max = EH(ih_max,1);

t_max = T(ih_max,1);

t_max_b = t_max + td;

t_min_b = t_max - td;

if (t_min_b < 0)
    
    t_min_b = 0;
    
end

[y,i_min_b] = min(abs(T-t_min_b));

[y,i_max_b] = min(abs(T-t_max_b));

ev_max = max(EV(i_min_b:i_max_b,1));

ezhr = ev_max/eh_max;
