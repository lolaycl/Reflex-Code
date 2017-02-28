% function [tot] = sum_expo(t,z,vf,rho_g)
% c = 3e8;
% k = (c-vf)/(c+vf);
% 
% tot  = zeros();
% tot1 = zeros();
% tot2 = zeros();
% 
% n_max = 1e5;
% for n= 1:n_max
%     update1 = (-1)^(n+1)*rho_g^n*channel_base_current(k^(n-1)*(t-2*z/c)).*heaviside(t-z/vf);
%     tot1= tot1+update1;
%     
%    update2 = (-1)^n*rho_g^n*channel_base_current(k^n*t).*heaviside(t-z/vf);
%    tot2 = tot2+update2;
%    
%    tot = tot1+tot2;
%    if norm(update1)<(1e-5)&&norm(update2)<(1e-6)
%     n
%        break;
%    end
% end


function [tot] = sum_expo(t,z,vf,rho_g)
c = 3e8;
k = (c-vf)/(c+vf);

tot  = zeros();
tot1 = zeros();
tot2 = zeros();

n_max = 1e5;
for n= 1:30
    update1 = (-1)^(n+1)*rho_g^n*channel_base_current(k^(n-1)*(t-z/c)).*heaviside(t-z/vf);
    tot1= tot1+update1;
    
   update2 = (-1)^n*rho_g^n*channel_base_current(k^n*(t+z/c)).*heaviside(t-z/vf);
   tot2 = tot2+update2;
   
   tot = tot1+tot2;
%    if norm(update1)<(1e-5)&&norm(update2)<(1e-6)
%     n
%        break;
%    end
end
