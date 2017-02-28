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

% V = 2.999999e8; % V is the return stroke speed
V =1.9e8;

W = c; % w equal to speed of light (never used actually)


H = 10e3; % H is the channel height

ZLAM = 2000; % ZLAM is lambda for the MTLE model

ALPHA_TOWER = 0; % ALPHA_TOWER is upward current attenuation rate inside the tower. It is zero for no loss tower.

HR = 0; % we don't consider strike object
% HR = 553; %553; % HR is the height of the tower

RR = 0; % RR=(Zt-Zch)/(Zt+Zg)  Zt=ZcH
% RR = -0.5; % RR is the current reflection coefficient at the tower top

% RG = 0.48; % RG is the current reflection coefficient at the tower base
% RG = 0.8;
% RG = 1;

HAB = 10; % HAB,HAM are heights at which the current is observed

HAM = 500; %500; % HAB,HAM are heights at which the current is observed %this is the value that changes everything

ERR = 1.0e-7; % ERR is parameter for the precision of the numerical integration

ERR2 = 1.0e-7; % ERR2 is parameter for the precision of the numerical integration

ZS = 0.0; % ZS is height of the observation point over ground

SIG = 0.002e100; % SIG is conductivity of the soil

EPSR = 1.0; % EPSR is dielectric constant of the soil

% D1 = 2e3; % D1 is distance of the observation point (projection on X axis)
% D1 = 50;
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

problem_type = 'single_point'; % 'single_point' - 'multi_point'

if strcmp(problem_type,'single_point')
    
    MPROBLEM = 1;
    
elseif strcmp(problem_type,'multi_point')
    
    multi_type = 'distance_ez_er_ratio'; % 'decay' - 'height' - 'speed' - 'pulse_width' - 'distance' - 'distance_ez_er_ratio'
    
    MPROBLEM = 2;
    
    NL = 21;
    
end

% if strcmp(model_type,'TL')
%     
%     MODEL = 1;
%     
% elseif strcmp(model_type,'MTLL')
%     
%     MODEL = 2;
%     
% elseif strcmp(model_type,'MTLE')
%     
%     MODEL = 3;
%     
% elseif strcmp(model_type,'BG')
%     
%     MODEL = 4;
%     
% elseif strcmp(model_type,'TCS')
%     
%     MODEL = 5;
%     
% end

% MODEL=1; %model_type = 'TL';
% if strcmp(problem_type,'single_point')
%     
%     write_infile(MPROBLEM,MODEL,I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B,V,W,H,ZLAM,ALPHA_TOWER,HR,RR,RG,HAB,HAM,ERR,...
%         ERR2,ZS,SIG,EPSR,D1,D2,NSIMP,DELTATE1,TINT,DTS,TMAX,QEV,QH,QEH,KEY1,KEY2,F1,F2,F3,F4);
%     
% tic  
%      [T_yliu_TL I0_yliu_TL I0_T_yliu_TL I1_T_yliu_TL I2_T_yliu_TL EV_TOTAL_yliu_TL HPHI_TOTAL_yliu_TL EH_TOTAL_yliu_TL]=run_engineering_model_matlab_current_vectorized_yuecen(MPROBLEM,MODEL,I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B,V,W,H,ZLAM,ALPHA_TOWER,HR,RR,RG,HAB,HAM,ERR,ERR2,ZS,SIG,EPSR,D1,D2,NSIMP,DELTATE1,TINT,DTS,TMAX,QEV,QH,QEH,KEY1,KEY2,F1,F2,F3,F4);
% toc
% 
% end

% MODEL=2; %model_type = 'MTLL';
% if strcmp(problem_type,'single_point')
%     
%     write_infile(MPROBLEM,MODEL,I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B,V,W,H,ZLAM,ALPHA_TOWER,HR,RR,RG,HAB,HAM,ERR,...
%         ERR2,ZS,SIG,EPSR,D1,D2,NSIMP,DELTATE1,TINT,DTS,TMAX,QEV,QH,QEH,KEY1,KEY2,F1,F2,F3,F4);
%     
% tic  
%      [T_yliu_MTLL I0_yliu_MTLL I0_T_yliu_MTLL I1_T_yliu_MTLL I2_T_yliu_MTLL EV_TOTAL_yliu_MTLL HPHI_TOTAL_yliu_MTLL EH_TOTAL_yliu_MTLL]=run_engineering_model_matlab_current_vectorized_yuecen(MPROBLEM,MODEL,I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B,V,W,H,ZLAM,ALPHA_TOWER,HR,RR,RG,HAB,HAM,ERR,ERR2,ZS,SIG,EPSR,D1,D2,NSIMP,DELTATE1,TINT,DTS,TMAX,QEV,QH,QEH,KEY1,KEY2,F1,F2,F3,F4);
% toc
% 
% end
%  plot_current_fields_yuecen(T_yliu_MTLL, I0_yliu_MTLL ,I0_T_yliu_MTLL, I1_T_yliu_MTLL, I2_T_yliu_MTLL, EV_TOTAL_yliu_MTLL, HPHI_TOTAL_yliu_MTLL ,EH_TOTAL_yliu_MTLL);
 
MODEL=3; %model_type = 'MTLE'; 
if strcmp(problem_type,'single_point')
    
    write_infile(MPROBLEM,MODEL,I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B,V,W,H,ZLAM,ALPHA_TOWER,HR,RR,RG,HAB,HAM,ERR,...
        ERR2,ZS,SIG,EPSR,D1,D2,NSIMP,DELTATE1,TINT,DTS,TMAX,QEV,QH,QEH,KEY1,KEY2,F1,F2,F3,F4);
    
tic  
     [T_lola I0 I0_T_lola I1_T_lola I2_T_lola EV_TOTAL_lola HPHI_TOTAL_lola EH_TOTAL_lola]=Lola(MPROBLEM,MODEL,I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B,V,W,H,ZLAM,ALPHA_TOWER,HR,RR,RG,HAB,HAM,ERR,ERR2,ZS,SIG,EPSR,D1,D2,NSIMP,DELTATE1,TINT,DTS,TMAX,QEV,QH,QEH,KEY1,KEY2,F1,F2,F3,F4);
toc

end
% plot_current_fields_yuecen(T_yliu_MTLE, I0_yliu_MTLE ,I0_T_yliu_MTLE, I1_T_yliu_MTLE, I2_T_yliu_MTLE, EV_TOTAL_yliu_MTLE, HPHI_TOTAL_yliu_MTLE ,EH_TOTAL_yliu_MTLE);

% 
% MODEL=4; %model_type = 'BG';
% if strcmp(problem_type,'single_point')
%     
%     write_infile(MPROBLEM,MODEL,I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B,V,W,H,ZLAM,ALPHA_TOWER,HR,RR,RG,HAB,HAM,ERR,...
%         ERR2,ZS,SIG,EPSR,D1,D2,NSIMP,DELTATE1,TINT,DTS,TMAX,QEV,QH,QEH,KEY1,KEY2,F1,F2,F3,F4);
%     
% tic  
%      [T_yliu_BG I0_yliu_BG I0_T_yliu_BG I1_T_yliu_BG I2_T_yliu_BG EV_TOTAL_yliu_BG HPHI_TOTAL_yliu_BG EH_TOTAL_yliu_BG]=run_engineering_model_matlab_current_vectorized_yuecen(MPROBLEM,MODEL,I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B,V,W,H,ZLAM,ALPHA_TOWER,HR,RR,RG,HAB,HAM,ERR,ERR2,ZS,SIG,EPSR,D1,D2,NSIMP,DELTATE1,TINT,DTS,TMAX,QEV,QH,QEH,KEY1,KEY2,F1,F2,F3,F4);
% toc
% 
% end
% %  plot_current_fields_yuecen(T_yliu_BG, I0_yliu_BG ,I0_T_yliu_BG, I1_T_yliu_BG, I2_T_yliu_BG, EV_TOTAL_yliu_BG, HPHI_TOTAL_yliu_BG ,EH_TOTAL_yliu_BG);
 
% MODEL=5; %model_type = 'TCS';
% if strcmp(problem_type,'single_point')
%     
%     write_infile(MPROBLEM,MODEL,I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B,V,W,H,ZLAM,ALPHA_TOWER,HR,RR,RG,HAB,HAM,ERR,...
%         ERR2,ZS,SIG,EPSR,D1,D2,NSIMP,DELTATE1,TINT,DTS,TMAX,QEV,QH,QEH,KEY1,KEY2,F1,F2,F3,F4);
%     
% tic  
%      [T_yliu_TCS I0_yliu_TCS I0_T_yliu_TCS I1_T_yliu_TCS I2_T_yliu_TCS EV_TOTAL_yliu_TCS HPHI_TOTAL_yliu_TCS EH_TOTAL_yliu_TCS]=run_engineering_model_matlab_current_vectorized_yuecen(MPROBLEM,MODEL,I01,ALPHA,BETA,I02,GAMMA,DELTA,IW1,TAU11,TAU21,NW1,IW2,TAU12,TAU22,NW2,I0B,T0B,V,W,H,ZLAM,ALPHA_TOWER,HR,RR,RG,HAB,HAM,ERR,ERR2,ZS,SIG,EPSR,D1,D2,NSIMP,DELTATE1,TINT,DTS,TMAX,QEV,QH,QEH,KEY1,KEY2,F1,F2,F3,F4);
% toc
% 
% end

% plot_five_models(T_yliu_TL,I0_T_yliu_TL,I2_T_yliu_TL, EV_TOTAL_yliu_TL, HPHI_TOTAL_yliu_TL, EH_TOTAL_yliu_TL, T_yliu_MTLL,I0_T_yliu_MTLL,I2_T_yliu_MTLL,EV_TOTAL_yliu_MTLL,HPHI_TOTAL_yliu_MTLL,EH_TOTAL_yliu_MTLL,T_yliu_MTLE,I0_T_yliu_MTLE, I2_T_yliu_MTLE,EV_TOTAL_yliu_MTLE, HPHI_TOTAL_yliu_MTLE, EH_TOTAL_yliu_MTLE, T_yliu_BG,I0_T_yliu_BG , I2_T_yliu_BG, EV_TOTAL_yliu_BG ,HPHI_TOTAL_yliu_BG, EH_TOTAL_yliu_BG ,T_yliu_TCS , I0_T_yliu_TCS , I2_T_yliu_TCS ,EV_TOTAL_yliu_TCS, HPHI_TOTAL_yliu_TCS, EH_TOTAL_yliu_TCS);



