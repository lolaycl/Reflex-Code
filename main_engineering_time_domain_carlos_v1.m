clc;

clear all;

c = 3e8;

current_type = 'farhad_subsequent_return_stroke';
% 'farhad_first_return_stroke'
% 'farhad_subsequent_return_stroke'
% 'farhad_subsequent_return_stroke_decay'
% 'du_first_return_stroke'
% 'du_subsequent_return_stroke'
% 'carlo_alberto_subsequent_return_stroke'
% 'fredholm_subsequent_return_stroke'
% 'marcos_horizontal_electric_field'
% 'barbosa_linear_ramp_current_original'
% 'abbas'

[I01 TAU31 TAU41 ALPHA BETA I02 TAU32 TAU42 GAMMA DELTA IW1 TAU11 TAU21 NW1 IW2 TAU12 TAU22 NW2 I0B T0B] = current_wave(current_type);

V = 1e8; % V is the return stroke speed
% V = 1.5e8;
% V =1.9e8;

W = c; % w equal to speed of light (never used actually)


H = 10e3; % H is the channel height

ZLAM = 2000; % ZLAM is lambda for the MTLE model
% ZLAM = 5000;

ALPHA_TOWER = 0; % ALPHA_TOWER is upward current attenuation rate inside the tower. It is zero for no loss tower.

HR = 0; % we don't consider strike object

RR = 0; % RR=(Zt-Zch)/(Zt+Zg)  Zt=ZcH RR is the current reflection coefficient at the tower top
 
% RG = 0.5; % RG is the current reflection coefficient at the tower base
RG = 0.8;
% RG = 1;

HAB = 10; % HAB,HAM are heights at which the current is observed

HAM = 500; %500; % HAB,HAM are heights at which the current is observed %this is the value that changes everything

ERR = 1.0e-7; % ERR is parameter for the precision of the numerical integration

ERR2 = 1.0e-7; % ERR2 is parameter for the precision of the numerical integration

ZS = 0.0; % ZS is height of the observation point over ground

SIG = 0.002e100; % SIG is conductivity of the soil

EPSR = 1.0; % EPSR is dielectric constant of the soil

% D1 = 50; % D1 is distance of the observation point (projection on X axis)
D1 = 5e3;
% D1 = 100e3;

D2 = 0.0; % D2 is distance of the observation point (projection on Y axis)

NSIMP = 10; % NSIMP is parameter for precision of the integration routine of Simpson

DELTATE1 = 5e-8; % DELTATE1 is the time step of the interval [0:TINT]

TINT = 2e-6; % TINT is the first interval limit

DTS = 1.0; % DTS*DELTATE1 is the time step for the period [TINT:TMAX]

TMAX = H/V; % TMAx is the second time interval limit

QEV = 1.e0;

QH = 1.e0;

QEH = 1.e0;

KEY1 = 0;

KEY2 = 0;

F1 = 1.0e0; % E-FIELD SCALE [RESULTS IN 1.e-3=kV/m 1=V/m]

F2 = 1.0e-3; % CURRENT SCALE [RESULTS IN 1.e-3=kA 1=A]

F3 = 1.0e6; % TIME SCALE [RESULTS IN 1.e6=microsec 1=sec]

F4 = 1.0e0; % H-FIELD SCALE [RESULTS IN 1.e-3=kA/m 1=A/m]

model_type = 'MTLE'; %model_type = 'TL';

if strcmp(model_type,'TL')
    
    MODEL = 1;
    
elseif strcmp(model_type,'MTLL')
    
    MODEL = 2;
    
elseif strcmp(model_type,'MTLE')
    
    MODEL = 3;
    
elseif strcmp(model_type,'BG')
    
    MODEL = 4;
    
elseif strcmp(model_type,'TCS')
    
    MODEL = 5;
    
end

problem_type = 'single_point'; % 'single_point' - 'multi_point'

if strcmp(problem_type,'single_point')
    
    MPROBLEM = 1;
    
elseif strcmp(problem_type,'multi_point')
    
    multi_type = 'distance_ez_er_ratio'; % 'decay' - 'height' - 'speed' - 'pulse_width' - 'distance' - 'distance_ez_er_ratio'
    
    MPROBLEM = 2;
    
    NL = 21;
    
end

if strcmp(problem_type,'single_point')
    
    write_infile(MPROBLEM,MODEL,I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B,V,W,H,ZLAM,ALPHA_TOWER,HR,RR,RG,HAB,HAM,ERR,...
        ERR2,ZS,SIG,EPSR,D1,D2,NSIMP,DELTATE1,TINT,DTS,TMAX,QEV,QH,QEH,KEY1,KEY2,F1,F2,F3,F4);
    
tic  
     [T_lola I0 I0_T_lola I1_T_lola I2_T_lola EV_TOTAL_lola HPHI_TOTAL_lola EH_TOTAL_lola]=Lola(MPROBLEM,MODEL,I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B,V,W,H,ZLAM,ALPHA_TOWER,HR,RR,RG,HAB,HAM,ERR,ERR2,ZS,SIG,EPSR,D1,D2,NSIMP,DELTATE1,TINT,DTS,TMAX,QEV,QH,QEH,KEY1,KEY2,F1,F2,F3,F4);
toc
plot_current_fields_yuecen(T_lola,I0,I0_T_lola,I1_T_lola,I2_T_lola,EV_TOTAL_lola,HPHI_TOTAL_lola,EH_TOTAL_lola);

% tic  
%      [T_yliu I0_yliu I0_T_yliu I1_T_yliu I2_T_yliu EV_TOTAL_yliu HPHI_TOTAL_yliu EH_TOTAL_yliu]=run_engineering_model_matlab_current_vectorized_yuecen(MPROBLEM,MODEL,I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B,V,W,H,ZLAM,ALPHA_TOWER,HR,RR,RG,HAB,HAM,ERR,ERR2,ZS,SIG,EPSR,D1,D2,NSIMP,DELTATE1,TINT,DTS,TMAX,QEV,QH,QEH,KEY1,KEY2,F1,F2,F3,F4);
% toc
%  plot_current_fields_yuecen(T_yliu,I0_yliu,I0_T_yliu,I1_T_yliu,I2_T_yliu,EV_TOTAL_yliu,HPHI_TOTAL_yliu,EH_TOTAL_yliu);

% tic  %_current_vectorized
%      [T_carlos I0_carlos I0_T_carlos I1_T_carlos I2_T_carlos EV_TOTAL_carlos HPHI_TOTAL_carlos EH_TOTAL_carlos]=run_engineering_model_matlab_current_vectorized(MPROBLEM,MODEL,I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B,V,W,H,ZLAM,ALPHA_TOWER,HR,RR,RG,HAB,HAM,ERR,ERR2,ZS,SIG,EPSR,D1,D2,NSIMP,DELTATE1,TINT,DTS,TMAX,QEV,QH,QEH,KEY1,KEY2,F1,F2,F3,F4);
% toc
% plot_current_fields_yuecen(T_carlos,I0_carlos,I0_T_carlos,I1_T_carlos,I2_T_carlos,EV_TOTAL_carlos,HPHI_TOTAL_carlos,EH_TOTAL_carlos);

% tic
%     run_engineering_model()
%     [T I0_T I1_T I2_T EV_TOTAL HPHI_TOTAL EH_TOTAL] = read_out_file();
% toc  

%  plot_currents_fields_carlos_v1(T_yliu,I0_T_yliu,I1_T_yliu,I2_T_yliu,EV_TOTAL_yliu,HPHI_TOTAL_yliu,EH_TOTAL_yliu);
%  plot_currents_fields_carlos_v1(T_carlos,I0_T_carlos,I1_T_carlos,I2_T_carlos,EV_TOTAL_carlos,HPHI_TOTAL_carlos,EH_TOTAL_carlos);
%  plot_currents_fields_carlos_v1(T,I0_T,I1_T,I2_T,EV_TOTAL,HPHI_TOTAL,EH_TOTAL);

%  plot_current_fields_yuecen(T_yliu,I0_yliu,I0_T_yliu,I1_T_yliu,I2_T_yliu,EV_TOTAL_yliu,HPHI_TOTAL_yliu,EH_TOTAL_yliu);

%  plot_currents_fields(T_yliu,I0_T_yliu,I1_T_yliu,I2_T_yliu,EV_TOTAL_yliu,HPHI_TOTAL_yliu,EH_TOTAL_yliu, T_carlos,I0_T_carlos,I1_T_carlos,I2_T_carlos,EV_TOTAL_carlos,HPHI_TOTAL_carlos,EH_TOTAL_carlos);
%  plot_currents_fields(T, I0_T, I1_T, I2_T, EV_TOTAL, HPHI_TOTAL, EH_TOTAL, I0_T_carlos, I1_T_carlos, I2_T_carlos, EV_TOTAL_carlos, HPHI_TOTAL_carlos, EH_TOTAL_carlos);
%  plot_currents_fields(T_lola,I0_T_lola,I1_T_lola,I2_T_lola,EV_TOTAL_lola,HPHI_TOTAL_lola,EH_TOTAL_lola, T_carlos,I0_T_carlos,I1_T_carlos,I2_T_carlos,EV_TOTAL_carlos,HPHI_TOTAL_carlos,EH_TOTAL_carlos);
% plot_currents_fields(T_lola,I0_T_lola,I1_T_lola,I2_T_lola,EV_TOTAL_lola,HPHI_TOTAL_lola,EH_TOTAL_lola, T_yliu,I0_T_yliu,I1_T_yliu,I2_T_yliu,EV_TOTAL_yliu,HPHI_TOTAL_yliu,EH_TOTAL_yliu);
    
% % % % %     t = T';
% % % % %     
% % % % %     ez = EV_TOTAL';
% % % % %     
% % % % %     er = EH_TOTAL';
% % % % %     
% % % % %     hphi = HPHI_TOTAL';
% % % % %     
% % % % %     save ('t','t');
% % % % %     
% % % % %     save ('ez_100','ez');
% % % % %     
% % % % %     save ('er_100','er');
% % % % %     
% % % % %     save ('hphi_100','hphi');
% % % % % 
% % % % %     save ('hphi_100_z_00','hphi');

% elseif strcmp(problem_type,'multi_point')
% 
%     if strcmp(multi_type,'decay')
% 
%         dlam = linspace(1500,4500,NL);
%         
%         tz = zeros(1,NL);
% 
%         for jj = 1:NL
% 
%             jj
% 
%             ZLAM = dlam(1,jj);
% 
%             write_infile(MPROBLEM,MODEL,I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B,V,W,H,ZLAM,ALPHA_TOWER,HR,RR,RG,HAB,HAM,ERR,...
%                 ERR2,ZS,SIG,EPSR,D1,D2,NSIMP,DELTATE1,TINT,DTS,TMAX,QEV,QH,QEH,KEY1,KEY2,F1,F2,F3,F4);
% 
%             run_engineering_model()
% 
%             [T I0_T I1_T I2_T EV_TOTAL HPHI_TOTAL EH_TOTAL] = read_out_file();
% 
%             tz(1,jj) = extract_zero_cross_time(T,EV_TOTAL);
%             
%         end
% 
%         figure;
%         
%         plot(dlam,tz)
%         
%         grid off
%         
%         make_label('\lambda','Time of zero crossing (\mus)');
% 
%     elseif strcmp(multi_type,'height')
% 
%         dheight = linspace(5e3,15e3,NL);
% 
%         for jj = 1:NL
% 
%             jj
% 
%             H = dheight(1,jj);
%             
%             TMAX = H/V;
% 
%             write_infile(MPROBLEM,MODEL,I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B,V,W,H,ZLAM,ALPHA_TOWER,HR,RR,RG,HAB,HAM,ERR,...
%                 ERR2,ZS,SIG,EPSR,D1,D2,NSIMP,DELTATE1,TINT,DTS,TMAX,QEV,QH,QEH,KEY1,KEY2,F1,F2,F3,F4);
% 
%             run_engineering_model()
% 
%             [T I0_T I1_T I2_T EV_TOTAL HPHI_TOTAL EH_TOTAL] = read_out_file();
% 
%             tz(1,jj) = extract_zero_cross_time(T,EV_TOTAL);
%             
%         end
% 
%         figure;
%         
%         plot(dheight/1e3,tz)
%         
%         grid off
%         
%         make_label('H (km)','Time of zero crossing (\mus)');
% 
%     elseif strcmp(multi_type,'speed')
% 
%         dspeed = linspace(c/3,c/2,NL);
%         
%         tz = zeros(1,NL);
% 
%         for jj = 1:NL
% 
%             jj
% 
%             V = dspeed(1,jj);
%             
%             TMAX = H/V;
% 
% 
%             write_infile(MPROBLEM,MODEL,I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B,V,W,H,ZLAM,ALPHA_TOWER,HR,RR,RG,HAB,HAM,ERR,...
%                 ERR2,ZS,SIG,EPSR,D1,D2,NSIMP,DELTATE1,TINT,DTS,TMAX,QEV,QH,QEH,KEY1,KEY2,F1,F2,F3,F4);
% 
%             run_engineering_model()
% 
%             [T I0_T I1_T I2_T EV_TOTAL HPHI_TOTAL EH_TOTAL] = read_out_file();
% 
%             tz(1,jj) = extract_zero_cross_time(T,EV_TOTAL);
%             
%         end
% 
%         figure;
%         
%         plot(dspeed/1e8,tz)
%         
%         grid off
%         
%         make_label('v (\times10^8 m/s)','Time of zero crossing (\mus)');
% 
%     elseif strcmp(multi_type,'pulse_width')
%         
%         dTAU31 = linspace(50e-6,150e-6,NL);
% 
%         tz = zeros(1,NL);
%         
%         pw = zeros(1,NL);
% 
%         for jj = 1:NL
% 
%             jj
% 
%             TAU31 = dTAU31(1,jj);
%             
%             ALPHA = 1/TAU31;
% 
%             [pw(1,jj)] = carlo_alberto(I01,TAU31,TAU41,IW1,TAU11,TAU21,NW1);
% 
%             write_infile(MPROBLEM,MODEL,I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B,V,W,H,ZLAM,ALPHA_TOWER,HR,RR,RG,HAB,HAM,ERR,...
%                 ERR2,ZS,SIG,EPSR,D1,D2,NSIMP,DELTATE1,TINT,DTS,TMAX,QEV,QH,QEH,KEY1,KEY2,F1,F2,F3,F4);
% 
%             run_engineering_model()
% 
%             [T I0_T I1_T I2_T EV_TOTAL HPHI_TOTAL EH_TOTAL] = read_out_file();
% 
%             tz(1,jj) = extract_zero_cross_time(T,EV_TOTAL);
%             
%         end
% 
%         figure;
%         
%         plot(pw,tz)
%         
%         make_label('Pulse width (\mus)','Time of zero crossing (\mus)');
% 
%     elseif strcmp(multi_type,'distance')
% 
%         dD1 = linspace(100e3,200e3,NL);
% 
%         tz = zeros(1,NL);
% 
%         for jj = 1:NL
% 
%             jj
% 
%             D1 = dD1(1,jj);
% 
%             write_infile(MPROBLEM,MODEL,I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B,V,W,H,ZLAM,ALPHA_TOWER,HR,RR,RG,HAB,HAM,ERR,...
%                 ERR2,ZS,SIG,EPSR,D1,D2,NSIMP,DELTATE1,TINT,DTS,TMAX,QEV,QH,QEH,KEY1,KEY2,F1,F2,F3,F4);
% 
%             run_engineering_model()
% 
%             [T I0_T I1_T I2_T EV_TOTAL HPHI_TOTAL EH_TOTAL] = read_out_file();
% 
%             tz(1,jj) = extract_zero_cross_time(T,EV_TOTAL);
%             
%         end
% 
%         figure;
%         
%         plot(dD1/1e3,tz)
%         
%         make_label('Distance from the channel base','Time of zero crossing (\mus)');
%         
%     elseif strcmp(multi_type,'distance_ez_er_ratio')
%         
%         dD1 = linspace(10,1500,NL);
% 
%         ezhr = zeros(1,NL);
% 
%         for jj = 1:NL
% 
%             jj
% 
%             D1 = dD1(1,jj);
% 
%             write_infile(MPROBLEM,MODEL,I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B,V,W,H,ZLAM,ALPHA_TOWER,HR,RR,RG,HAB,HAM,ERR,...
%                 ERR2,ZS,SIG,EPSR,D1,D2,NSIMP,DELTATE1,TINT,DTS,TMAX,QEV,QH,QEH,KEY1,KEY2,F1,F2,F3,F4);
% 
%             run_engineering_model()
% 
%             [T I0_T I1_T I2_T EV_TOTAL HPHI_TOTAL EH_TOTAL] = read_out_file();
% 
%             [ezhr(1,jj)] = obtain_ez_eh_ratio(T,EV_TOTAL,EH_TOTAL);
%             
%         end
% 
%         figure;
%         
%         plot(dD1,ezhr)
%         
%         make_label('r (m)','e_z/e_r');
%         
%     end
%     
end

