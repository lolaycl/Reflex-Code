function [tz] = extract_zero_cross_time(T,Ezt)

Nd = size(T,1);

pt = 2;

while ((pt > 0) && (pt<(Nd-1)))
    
    if ((Ezt(pt-1,1) > 0) && (Ezt(pt+1,1) < 0))
        
        break;
        
    end
    
    pt = pt + 1;
    
end

tz = T(pt,1);

