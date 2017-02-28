function [] = plot_ez_eh_energy_ratio(T,EV_TOTAL,EH_TOTAL,NL)

Nt = size(T,1);

tilte_anlge_d_vector = linspace(0,1,NL);

tilte_angle_r_vector = tilte_anlge_d_vector*pi/180;

energy_non_contaminated = 0;

energy_contaminated = zeros(1,NL);

energy_error = zeros(1,NL);

peak_error = zeros(1,NL);

for jj = 1:Nt-1
    
    energy_non_contaminated = energy_non_contaminated + 0.5*(EH_TOTAL(jj+1,1)^2+EH_TOTAL(jj,1)^2)*(T(jj+1,1)-T(jj,1));
    
end

for ii = 1:NL
    
    tilte_angle_r = tilte_angle_r_vector(1,ii);
    
    EH_TOTAL_M = EH_TOTAL - (EH_TOTAL*cos(tilte_angle_r)-EV_TOTAL*sin(tilte_angle_r));
    
    for jj = 1:Nt-1
        
        energy_contaminated(1,ii) = energy_contaminated(1,ii) + 0.5*(EH_TOTAL_M(jj+1,1)^2+EH_TOTAL_M(jj,1)^2)*(T(jj+1,1)-T(jj,1));
        
    end
    
end

energy_contaminated = sqrt(energy_contaminated);

energy_non_contaminated = sqrt(energy_non_contaminated);

energy_error = (energy_contaminated-energy_non_contaminated)/energy_non_contaminated*100;

figure;

plot(tilte_anlge_d_vector,abs(energy_error),'k');

make_label('Title angle (degree)','Energy error');
