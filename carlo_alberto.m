function [pw] = carlo_alberto(I02,T3,T4,I01,T1,T2,n)

% Subsequent return stroke

t1 = [0:0.001e-6:5e-6];
t2 = [5e-6:0.1e-6:500e-6];
t = [t1 t2];

N = size(t,2);
It = zeros(1,N);

eta = exp(-T1/T2*(n*T2/T1)^(1/n));

for jj = 1:N
    tt = t(1,jj);
    It(1,jj) = (I01/eta*((tt/T1)^n))/(1+(tt/T1)^n)*exp(-tt/T2)+I02*(exp(-tt/T3)-exp(-tt/T4));
end

pw_criterion = 0.05*max(It);

pt = 1;
start_point = N;
end_point = N;
while ((pt > 0) && (pt < N))
    if ((It(1,pt) > pw_criterion) && (pt < start_point))
        start_point = pt;
    end
    if(((pt > start_point) && (It(1,pt) < pw_criterion)))
        end_point = pt;
        break
    end
    pt = pt + 1;
end

pw = (t(1,end_point) - t(1,start_point))/1e-6;

% plot(t,It)
% hold on