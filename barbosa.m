function [It] = barbosa(t)

Nt = size(t,2);
It = zeros(1,Nt);
I01 = 8.5e3;
T11 = 0.5e-6;
T21 = 100000e-6;
n1 = 2;
I02 = 3.2e3;
T12 = 1e-6;
T22 = 60000e-6;
n2 = 2;

eta1 = exp(-T11/T21*(n1*T21/T11)^(1/n1));
eta2 = exp(-T12/T22*(n2*T22/T12)^(1/n2));
It = (I01/eta1*((t/T11).^n1))./(1+(t/T11).^n1).*exp(-t/T21)+(I02/eta2*((t/T12).^n2))./(1+(t/T12).^n2).*exp(-t/T22);