function [T I0_T I1_T I2_T EV_TOTAL HPHI_TOTAL EH_TOTAL]=run_engineering_model_matlab_current...
    (MPROBLEM_in,MODEL_in,I01_in,ALPHA_in,BETA_in,I02_in,GAMMA_in,DELTA_in,IW1_in,TAU11_in,TAU21_in,...
     NW1_in,IW2_in,TAU12_in,TAU22_in,NW2_in,I0B_in,T0B_in,V_in,W_in,H_in,ZLAM_in,ALPHA_TOWER_in,HR_in,RR_in,RG_in,HAB,HAM,ERR_in,...
     ERR2_in,ZS_in,SIG,EPSR,D1,D2,NSIMP_in,DELTATE1,TINT,DTS,TMAX,QEV_in,QH_in,QEH_in,KEY1_in,KEY2_in,F1,F2,F3,F4)
% clear
global MU C PI EPS IFRONT % /COM1/
global V VC W H ZLAM HC  % /COM2/
global ZS R R1 FV        % /COM3/ on other hand: R R1 FV are calculated inside the program (down)
global IW1 TAU11 TAU21 IW2 TAU12 TAU22 NW2 NW1 ETA ETA2 %/COM4/ %ETA AND ETA2 ARE CALCULATED INSIDE
global ALPHA BETA I01 GAMMA DELTA I02 %/COM5/
global HR RR RG ALPHA_TOWER MODEL MPROBLEM NC %/COM6/, NC is calculated inside
global ERR ERR2 %/COM7/
global NSIMP    %/COM8/
global RZ %/COM9/
global TT %/COMT/ 
global I0B T0B %/COMBARB/ INPUT VALUES
global EV1 EV2 EV3 EV4 EV5 EV6 EV7 EV8 EV9 EV10 EV11 EV12 EV1TO EV2TO %/EV/ CALCULATED VALUES
global EH1 EH2 EH3 EH4 EH5 EH6 EH7 EH8 EH9 EH10 EH11 EH12 EH1TO EH2TO %/EH/ CALCULATED VALUES
global HP1 HP2 HP3 HP4 HP5 HP6 HP7 HP8 HP1TO HP2TO %/HP/
global QEV QH QEH %/SCALE/ INPUT VALUES
global MAXLIM1 KEY1 MAXLIM2 KEY2 %/INTEG/, KEY1 and KEY2 are INPUT VALUES

% load('test_wksp1')
% MPROBLEM_in=MPROBLEM; MODEL_in=MODEL; I01_in=I01; ALPHA_in=ALPHA; BETA_in=BETA; I02_in=I02; GAMMA_in=GAMMA; DELTA_in=DELTA;
% IW1_in=IW1; TAU11_in=TAU11; TAU21_in=TAU21; NW1_in=NW1; IW2_in=IW2; TAU12_in=TAU12; TAU22_in=TAU22; NW2_in=NW2; I0B_in=I0B;
% T0B_in=T0B; V_in=V; W_in=W; H_in=H; ZLAM_in=ZLAM; ALPHA_TOWER_in=ALPHA_TOWER; HR_in=HR; RR_in=RR; RG_in=RG; HAB=HAB; HAM=HAM;
% ERR_in=ERR; ERR2_in=ERR2; ZS_in=ZS; SIG=SIG; EPSR=EPSR; D1=D1; D2=D2; NSIMP_in=NSIMP; DELTATE1=DELTATE1; TINT=TINT; DTS=DTS;
% TMAX=TMAX; QEV_in=QEV; QH_in=QH; QEH_in=QEH; KEY1_in=KEY1; KEY2_in=KEY2; F1=F1; F2=F2; F3=F3; F4=F4;
% clc

%%
%     I01,I02,ALPHA,BETA,GAMMA,DELTA: parameters of bi-exp function 
%     IW,NW,TAU1,TAU2: parameters of Heidler function
%     V: return-stroke speed
%     W: equal to speed of light (never used actually)
%     H: total height of the lightning channel
%     ZLAM: current attenuation factor with height for MTLE model
%     HR: height of the tower
%     RR: current reflection coefficient at the tower top
%     RG: current reflection coefficient at the tower base
%	  HAB,HAM: heights at which the current is observed
%	  ERR,ERR2: parameters for the precision of the numerical integration 
%     D1: distance of the observation point (projection on X axis)
%	  D1: distance of the observation point (projection on Y axis)
%     ZS: height of the observation point over ground
%	  SIG: conductivity of the soil
%	  EPSR: dielectric constant of the soil
%	  DELTATE1,TINT,DTS,TMAX: parameters for time resolution
%	  NSIMP: parameter for precision of the integration routine of Simpson

%%
%   Variables definition
	C=3e8; EPS=8.8542e-12; PI=pi;  MU=4*pi*1e-7;
    IFRONT=0;  
    MAXLIM1=0;
    MAXLIM2=0;   
    MX1=0; MX2=0;
    V=V_in; W=W_in; H=H_in; ZLAM=ZLAM_in;
    ZS=ZS_in;
    IW1=IW1_in;TAU11=TAU11_in;TAU21=TAU21_in;IW2=IW2_in;TAU12=TAU12_in;TAU22=TAU22_in;NW2=NW2_in;NW1=NW1_in;
    ALPHA=ALPHA_in;BETA=BETA_in;I01=I01_in;GAMMA=GAMMA_in;DELTA=DELTA_in;I02=I02_in;
    MODEL=MODEL_in; MPROBLEM=MPROBLEM_in; ALPHA_TOWER=ALPHA_TOWER_in;
    NC=0;
    HR=HR_in; RR=RR_in; RG=RG_in;
    ERR=ERR_in; ERR2=ERR2_in;
    NSIMP=NSIMP_in;
    I0B=I0B_in; T0B=T0B_in;
%   NW=n for the heidler function
    ETA=exp(-(TAU11/TAU21)*((NW1*TAU21/TAU11)^(1/(NW1))));
    ETA2=exp(-(TAU12/TAU22)*((NW2*TAU22/TAU12)^(1/(NW2))));
    EV1=0; EV2=0;  EV3=0;  EV4=0;  EV5=0;  EV6=0;  EV7=0;  EV8=0;  EV9=0;  EV10=0;  EV11=0;  EV12=0;  EV1TO=0;  EV2TO=0;  %/EV/ CALCULATED VALUES
    EH1=0; EH2=0;  EH3=0;  EH4=0;  EH5=0;  EH6=0;  EH7=0;  EH8=0;  EH9=0;  EH10=0;  EH11=0;  EH12=0;  EH1TO=0;  EH2TO=0;  %/EH/ CALCULATED VALUES
    HP1=0;  HP2=0;  HP3=0;  HP4=0;  HP5=0;  HP6=0;  HP7=0;  HP8=0;  HP1TO=0;  HP2TO=0;  %/HP/    
    QEV=QEV_in; QH=QH_in; QEH=QEH_in;
    KEY1=KEY1_in; KEY2=KEY2_in;
% Put correct names on the characters      
      if (F2==1e-3) 
          MCUR='   KA';
      else
          MCUR='    A';
      end
      
      if (F4==1e-3) 
          MHP='KA/(M';
      else
          MHP=' A/(M';
      end
      if (F1==1e-3) 
          MEVH='KV/(M';
      else
          MEVH=' V/(M';
      end

      if (F3==1e6) 
          MDCUR=strcat(MCUR,'/MICROS');
          MDHP=strcat(MHP,'*MICROS)');
          MDEV=strcat(MEVH,'*MICROS)');
      else
          MDCUR=strcat(MCUR,'/S');
          MDHP=strcat(MHP,'*S)');
          MDEV=strcat(MEVH,'*S)');
      end

%     ==============================================
%     SELECTION OF CURRENT-WAVE SPEED
%     ==============================================
%	VC: current-wave speed
%	C: speed of light

% 'TL'     MODEL = 1;
% 'MTLL'   MODEL = 2;
% 'MTLE'   MODEL = 3;
% 'BG'     MODEL = 4;
% 'TCS'    MODEL = 5;
    
    if ((MODEL==1) || (MODEL==2) || (MODEL==3)) 
        VC=V;
    elseif (MODEL==4) 
        VC=1e15; %infinite?
    elseif 	 (MODEL==5) 
        VC=-C;
    else
        VC='error';
    end
%%  ZLAM:current attenuation factor with height for MTLE model
    if (ZLAM~=0) 
        ZLAM=1/ZLAM; 
    end
    R=sqrt(D1^2+D2^2);
    DELTATE2=DELTATE1*DTS;
    NBE1=round(TINT/DELTATE1);
    NBE2=round((TMAX-TINT)/DELTATE2);
    NBE=NBE1+NBE2;
    

%% Preallocation of variables (we add one element at the end to compensate)
    EVO=zeros(1,NBE1+NBE2);
	EEL=zeros(1,NBE1+NBE2);
	EIND=zeros(1,NBE1+NBE2);
	ERAD=zeros(1,NBE1+NBE2);
	EV_TO=zeros(1,NBE1+NBE2);
	EHO=zeros(1,NBE1+NBE2);
	EHEL=zeros(1,NBE1+NBE2);
	EHIND=zeros(1,NBE1+NBE2);
	EHRAD=zeros(1,NBE1+NBE2);
	EH_TO=zeros(1,NBE1+NBE2);
	EELtower=zeros(1,NBE1+NBE2);
	EVtower=zeros(1,NBE1+NBE2);
	HPtower=zeros(1,NBE1+NBE2);
	HPHI=zeros(1,NBE1+NBE2);
	HIND=zeros(1,NBE1+NBE2);
	HRAD=zeros(1,NBE1+NBE2);
	HP_TO=zeros(1,NBE1+NBE2);
	ECO=zeros(1,NBE1+NBE2);
	CURT=zeros(1,NBE1+NBE2);
	DCURT=zeros(1,NBE1+NBE2);
	CURC=zeros(1,NBE1+NBE2);
	DCURC=zeros(1,NBE1+NBE2);
	CURB=zeros(1,NBE1+NBE2);
	DCURB=zeros(1,NBE1+NBE2);
	CURR1=zeros(1,NBE1+NBE2);
	DCURR1=zeros(1,NBE1+NBE2);
	TIMEE=zeros(1,NBE1+NBE2);
    
%%    =============================================      
%     CURRENT CALCULATION
%     =============================================
%	TE: time step
    TE=0.0;
    
    for N=1:1:NBE

          if (MPROBLEM == 1) 
            disp(strcat('N=',num2str(N),' NC=',num2str(NC)))
          end

          if (N<=NBE1) 
            DELTATE=DELTATE1;
          else
            DELTATE=DELTATE2;
          end
            TIMEE(N)=TE;

          if (HR>0) 
            NC=TE/(2*HR/C);
          else
            NC=0;
          end
    %   Actual wavefront current
        RZ=0;
        HC=HR+V*TE; %HR: height of the tower + V*time_step
        [CURC(N),DCURC(N)]=CURRENT(TE,HC,CURC(N),DCURC(N), RZ, HR, RR, RG, ALPHA_TOWER, NC, VC,C, MODEL);
        ifRONT=CURC(N);

    %	Turn-On Term fields
        ZS=abs(ZS);
        T01=sqrt(R*R+(HR-ZS)*(HR-ZS))/C;
        TT=TE+T01;
        A=1/(V*V)-1/(C*C);
        B=2*ZS/(C*C)-2*HR/(V*V)-2*TT/V;
        CC=TT*TT+2*TT*HR/V+HR*HR/(V*V)-R*R/(C*C);
        RADIC=abs(B*B-4*A*CC);
        Z1max=(-B-sqrt(RADIC))/(2*A);
        R1=sqrt(R*R+(ZS-Z1max)^2);
        RZ=R1;
        
        [CURR1(N),DCURR1(N)]=CURRENT(TT,Z1max,CURR1(N),DCURR1(N), RZ, HR, RR, RG, ALPHA_TOWER, NC, VC,C, MODEL);
        ifRONT=CURR1(N);
        FV=1/(1/V-(ZS-Z1max)/C/R1);


        ZS=-abs(ZS);
        TT=TE+T01;
        A=1/(V*V)-1/(C*C);
        B=2*ZS/(C*C)-2*HR/(V*V)-2*TT/V;
        % it necessary to calculate the effective distance, given the speed
        % of light propagation, (speed of light is assumed i think), to the
        % observation point
        CC=TT*TT+2*TT*HR/V+HR*HR/(V*V)-R*R/(C*C);
        RADIC=abs(B*B-4*A*CC);
        Z1max=(-B-sqrt(RADIC))/(2*A);
        R1=sqrt(R*R+(ZS-Z1max)^2);
        RZ=R1;
        [CURR1(N),DCURR1(N)]=CURRENT(TT,Z1max,CURR1(N),DCURR1(N), RZ, HR, RR, RG, ALPHA_TOWER, NC, VC,C, MODEL);
        ifRONT=CURR1(N);
        FV=1/(1/V-(ZS-Z1max)/C/R1);
        EV2TO=-ifRONT*R*R*FV/(4*PI*EPS*C*C*R1^3);
        HP2TO=ifRONT*R*FV/(4*PI*C*R1^2);
        EH2TO=-ifRONT*R*(ZS-Z1max)*FV/(4*PI*EPS*C*C*R1^3);

        RZ=0;
        [CURT(N),DCURT(N)]=CURRENT(TE,HAM,CURT(N),DCURT(N), RZ, HR, RR, RG, ALPHA_TOWER, NC, VC,C, MODEL); %current at the top 
        [CURB(N),DCURB(N)]=CURRENT(TE,HAB,CURB(N),DCURB(N), RZ, HR, RR, RG, ALPHA_TOWER, NC, VC,C, MODEL); %current at base


        TE=TE+DELTATE;
    end
 
        TM=TIMEE;
        I0_T=IPULSE(TM)*F2;
        DI0_T=DIPULSE(TM)*F2/F3;

        EHF=0; %Erase just for test current calc (Carlos)
        DHPHI=0; %Erase just for test current calc (Carlos)
        DEVO=0; %Erase just for test current calc (Carlos)

        T=TM*F3; EV=-EVO*F1; EVEL=-EEL*F1; EVIN=-EIND*F1; EVRD=-ERAD*F1; EV_TO=-EV_TO*F1; EHO=EHO*F1;
        EHF=EHF*F1; EHEL=EHEL*F1; EHIN=EHIND*F1; EHRD=EHRAD*F1; EH_TO=EH_TO*F1; HPHI=-HPHI*F4; HPIND=-HIND*F4;HPRAD=-HRAD*F4;
        HP_TO=HP_TO*F4; DEV=-DEVO*F1/F3; DHPHI=DHPHI*F4/F3; IT_T=CURT*F2; IM_T=CURC*F2; IB_T=CURB*F2; IFRONT_T= CURR1*F2;
        DIT_T=DCURT*F2/F3; DIM_T=DCURC*F2/F3; DIB_T=DCURB*F2/F3;
        EVtower_T=-EVtower*F1;HPtower_T=-HPtower*F4; EELtower_T=-EELtower*F1;
% 
%         % Now do the same as the inside of read_out_file
        I1_T=IT_T;
        I2_T=IB_T;
        
        EV_TOTAL = EV + EV_TO; %Vertical Eletric field
        HPHI_TOTAL = HPHI + HP_TO; 
        EH_TOTAL = EHF; 
        
%current at height QUO at time TY and its derivative    
function [CUR_out,DCUR_out]=CURRENT(TY,QUO,CUR,DCUR, RZ, HR, RR, RG, ALPHA_TOWER, NC, VC,C, MODEL)

  CUR=0;DCUR=0;
  B1=0;B2=0;B3=0; 
  J=0;
      
     
  CUR=0;
  DCUR=0;

    if (QUO<=HR)  
%  ===========1	
        for J=0:1:NC
            B1=TY-RZ/C-(HR-QUO)/C-2.*J*HR/C;
%             B1=B1+1; %I had to add this to make compatible with 1-indexing
            if (B1>=0) 
               CUR=CUR+IPULSE(B1)*(RR^J)*(RG^J)*exp(-J*HR*ALPHA_TOWER);  %I had to subs -1 to B1 this to make compatible with 1-indexing
               DCUR=DCUR+DIPULSE(B1)*(RR^J)*(RG^J)*exp(-J*HR*ALPHA_TOWER); %I had to subs -1 to B1 this to make compatible with 1-indexing
               if or(not(isreal(CUR)),not(isreal(DCUR)))
                 disp('Non-real value on the current routine') 
               end
            end
            B2=B1-2.*QUO/C;
%             B2=B2+1; %I had to add this to make compatible with 1-indexing
            if (B2>=0) 
               CUR=CUR+IPULSE(B2)*(RR^J)*(RG^(J+1))*exp(-(J*HR+QUO)*ALPHA_TOWER); %I had to subs -1 to B2 this to make compatible with 1-indexing
               DCUR=DCUR+DIPULSE(B2)*(RR^J)*(RG^(J+1))*exp(-(J*HR+QUO)*ALPHA_TOWER); %I had to subs -1 to B2 this to make compatible with 1-indexing
               if or(not(isreal(CUR)),not(isreal(DCUR)))
                 disp('Non-real value on the current routine') 
               end
            end
        end

        CUR=CUR*(1-RR);
        DCUR=DCUR*(1-RR);
    else
% ===========Z GREATER THAN HR 
        B1=TY-RZ/C-(QUO-HR)/VC;
%         B1=B1+1; %I had to add this to make compatible with 1-indexing
	    if (B1>=0) 
	     if (MODEL == 5) 
% B1 = B1 + (QUO-HR)/VC + (QUO-HR)/C
         end 
		 CUR=CUR+IPULSE(B1)*P(QUO);
	     DCUR=DCUR+DIPULSE(B1)*P(QUO);
         if or(not(isreal(CUR)),not(isreal(DCUR)))
            disp('Non-real value on the current routine') 
         end
        end
% C===========2
        B1=TY-RZ/C-(QUO-HR)/C;
%         B1=B1+1; %I had to add this to make compatible with 1-indexing
	    if (B1>=0) 
             CUR=CUR+IPULSE(B1)*(-RR);
             DCUR=DCUR+DIPULSE(B1)*(-RR);
             if or(not(isreal(CUR)),not(isreal(DCUR)))
                disp('Non-real value on the current routine') 
             end
        end
% C===========3
        for J=0:1:NC
          B1=TY-RZ/C-(QUO+HR)/C-2.*J*HR/C;
%           B1=B1+1; %I had to add this to make compatible with 1-indexing
              if (B1>=0)
                 CUR=CUR+(1-RR)*(1+RR)*(RR^J)*(RG^(J+1))*IPULSE(B1)*exp(-(J+1)*HR*ALPHA_TOWER);
                 DCUR=DCUR+(1-RR)*(1+RR)*(RR^J)*(RG^(J+1))*DIPULSE(B1)*exp(-(J+1)*HR*ALPHA_TOWER);
                 if or(not(isreal(CUR)),not(isreal(DCUR)))
                    disp('Non-real value on the current routine') 
                end
              end
        end
	
    end
CUR_out=CUR; DCUR_out=DCUR;
% disp(J)
%     ==============================================================
%                       CURRENT FUNCTIONS
%     ==============================================================

function IPULSE_calc=IPULSE(TA)
%	============================== 2Heidler+2exp
      global IW1 ETA TAU11 TAU21 IW2 ETA2 TAU12 TAU22 NW2 NW1
      global ALPHA BETA I01 GAMMA DELTA I02
	  global I0B T0B
%tau1=2-5e-07 tau2=2.5e-06
      FK=(TA./TAU11).^NW1./(1+((TA./TAU11).^NW1));
      FKK=(TA./TAU12).^NW2./(1+((TA./TAU12).^NW2));
      IPULSE_calc=(IW1/ETA).*FK.*exp(-TA./TAU21);
      IPULSE_calc=IPULSE_calc+(IW2/ETA2).*FKK.*exp(-TA./TAU22);
      IPULSE_calc=IPULSE_calc+I01*(exp(-ALPHA.*TA)-exp(-BETA.*TA));
	  IPULSE_calc=IPULSE_calc+I02*(exp(-GAMMA.*TA)-exp(-DELTA.*TA));
	if (TA<=T0B) 
        IPULSE_calc=IPULSE_calc+I0B.*TA./T0B;
    else
        IPULSE_calc=IPULSE_calc+I0B;
    end
    if not(isreal(IPULSE_calc))
         disp('Pulse imaginary') 
    end

function DIPULSE_calc=DIPULSE(TA)
%===============================d(2Heidler+2exp)/dt

      global IW1 ETA TAU11 TAU21 IW2 ETA2 TAU12 TAU22 NW2 NW1
      global ALPHA BETA I01 GAMMA DELTA I02
	  global I0B T0B

      %index_zero=(TA<=0); %Carlos added for solving a problem
      FK2=NW1./(1+((TA./TAU11).^NW1))-TA./TAU21;
      FK1=((TA./TAU11).^(NW1-1))./(1+((TA./TAU11).^NW1));
      FKK2=NW2./(1+((TA./TAU12).^NW2))-TA./TAU22;
      FKK1=((TA./TAU12).^(NW2-1))./(1+((TA./TAU12).^NW2));
      DIPULSE_calc=(IW1./(ETA.*TAU11)).*FK1.*exp(-TA./TAU21).*FK2;
      DIPULSE_calc=DIPULSE_calc+(IW2/(ETA2*TAU12)).*FKK1.*exp(-TA./TAU22).*FKK2;
      DIPULSE_calc=DIPULSE_calc+I01.*(BETA*exp(-BETA.*TA)-ALPHA.*exp(-ALPHA.*TA));
      DIPULSE_calc=DIPULSE_calc+I02.*(DELTA*exp(-DELTA.*TA)-GAMMA.*exp(-GAMMA.*TA));
      if (TA<=T0B) 
            DIPULSE_calc=DIPULSE_calc+I0B/T0B;
      end
      if not(isreal(DIPULSE_calc))
         disp('Pulse derivative imaginary') 
      end
      if or(isinf(DIPULSE_calc),isnan(DIPULSE_calc))
        disp('Infinite Number or NaN');
      end
     DIPULSE_calc=real(DIPULSE_calc);


function P_calc=P(ZZ)
%	============================== CHOOSING THE ATTENUATION FACTOR P(z) =====
      global H ZLAM
      global HR MODEL 


	if ((MODEL==1)||(MODEL==4)||(MODEL==5)) 
        P_calc=1;
    end
    if (MODEL==2) 
        P_calc=(1-(ZZ-HR)./H);
    end
 	if (MODEL==3) 
        P_calc=exp(-(ZZ-HR).*ZLAM);
    end
	
