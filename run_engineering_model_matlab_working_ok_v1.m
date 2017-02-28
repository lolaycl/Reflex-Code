function [T I0_T I1_T I2_T EV_TOTAL HPHI_TOTAL EH_TOTAL]=run_engineering_model_matlab...
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
% IW1_in=IW1; TAU11_in=TAU11; TAU21_in=TAU21; NW1_in=TAU21; IW2_in=TAU21; TAU12_in=TAU21; TAU22_in=TAU22; NW2_in=NW2; I0B_in=I0B;
% T0B_in=T0B; V_in=V; W_in=W; H_in=H; ZLAM_in=ZLAM; ALPHA_TOWER_in=ALPHA_TOWER; HR_in=HR; RR_in=RR; RG_in=RG; HAB=HAB; HAM=HAM;
% ERR_in=ERR; ERR2_in=ERR2; ZS_in=ZS; SIG=SIG; EPSR=EPSR; D1=D1; D2=D2; NSIMP_in=NSIMP; DELTATE1=DELTATE1; TINT=TINT; DTS=DTS;
% TMAX=TMAX; QEV_in=QEV; QH_in=QH; QEH_in=QEH; KEY1_in=KEY1; KEY2_in=KEY2; F1=F1; F2=F2; F3=F3; F4=F4;
% clc
% function run_engineering_model_matlab
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
%       REAL*8 MU,C,PI,EPS,V,VC,W,H,ZLAM,R,ZS,R1,FV,IFRONT,Z1max,RZ             
%       REAL*8 IW,ETA,TAU1,TAU2,IW2,ETA2,TAU12,TAU22,I0B,T0B
%       REAL*8 ALPHA,BETA,I01,GAMMA,DELTA,I02      
%       REAL*8 HR,RR,RG,ALPHA_TOWER,ERR,ERR2                           
%       REAL*8 DELTATE,DELTATE1,DELTATE2
%       REAL*8 TE,DTS,TMAX,D1,D2,TINT,TT,T01,A,B,CC,RADIC
%       REAL*8 EV1,EV2,EV3,EV4,EV5,EV6,EV7,EV8,EV9,EV10,EV11,EV12,QEV
%       REAL*8 EH1,EH2,EH3,EH4,EH5,EH6,EH7,EH8,EH9,EH10,EH11,EH12,QEH
%       REAL*8 HP1,HP2,HP3,HP4,HP5,HP6,HP7,HP8,QH
% 	    REAL*8 EV1TO,EV2TO,EH1TO,EH2TO,HP1TO,HP2TO
%       REAL*8 SIG,EPSR,MK1,MK,FF
% 	    REAL*8 F1,F2,F3,F4,TM,HAB,HAM,HC
%       REAL*8 IPULSE,DIPULSE,P
%       INTEGER NBE,NBE1,NBE2,NW,NW2,NC,N,NSIMP,MAXLIM1,KEY1,MAXLIM2,KEY2
%       INTEGER MX1,MX2,K,S,MODEL,MPROBLEM
%       REAL*8 DEVO,DHPHI,EHF
%       REAL*8 ,ALLOCATABLE::EVO(:),EEL(:),EIND(:),ERAD(:)
%       REAL*8 ,ALLOCATABLE::EHO(:),EHEL(:),EHIND(:),EHRAD(:)    
%       REAL*8 ,ALLOCATABLE::HPHI(:),HIND(:),HRAD(:),EVtower(:),HPtower(:)
% 	    REAL*8 ,ALLOCATABLE::EV_TO(:),EH_TO(:),HP_TO(:),EELtower(:)
%       REAL*8 ,ALLOCATABLE::ECO(:)
% 	    REAL*8 ,ALLOCATABLE::CURT(:),DCURT(:),CURC(:),DCURC(:),CURB(:),DCURB(:),CURR1(:),DCURR1(:)
%                    
%       REAL*8 ,ALLOCATABLE::TIMEE(:)
%       
%       CHARACTER*1 CHO,RE
%       CHARACTER*13 MDCUR,MDEV,MDHP
%       CHARACTER*5 MCUR,MEVH,MHP 
%       CHARACTER*15 NOME
% %      L'instruction global permet le regroupement de zones mémoires pouvant 
% %      être partagées par différentes unités de programme (fonctions, procédures).
% %      When a subroutine or function is exited, local variables become undefined.
% %      The same thing may happen with the variables stored in named COMMON blocks. 
% %      Therefore, it is possible to SAVE an entire named COMMON block 
% %      (but not individual variables in the block) in a procedure with the command
% %      SAVE /named common block1/, /named common block2/, ..., /named common blockn/
% %      Variables in a blank COMMON block never become undefined since the it
% %      is declared in the main program. Similarly, 
% %      if a named COMMON block is declared in the main program, 
% %      then it is unecessary to use the SAVE command in other program units. 
% %      Note that you cannot use the SAVE command on a blank COMMON block.
    
% 	  global MU C PI EPS IFRONT % /COM1/
%       global V VC W H ZLAM HC  % /COM2/
%       global ZS R R1 FV        % /COM3/ on other hand: R R1 FV are calculated inside the program (down)
%       global IW1 TAU11 TAU21 IW2 TAU12 TAU22 NW2 NW1 ETA ETA2 %/COM4/ %ETA AND ETA2 ARE CALCULATED INSIDE
%       global ALPHA BETA I01 GAMMA DELTA I02 %/COM5/
%       global HR RR RG ALPHA_TOWER MODEL MPROBLEM NC %/COM6/, NC is calculated inside
%       global ERR ERR2 %/COM7/
%       global NSIMP    %/COM8/
% 	  global RZ %/COM9/
% 	  global TT %/COMT/ 
% 	  global I0B T0B %/COMBARB/ INPUT VALUES
%       global EV1 EV2 EV3 EV4 EV5 EV6 EV7 EV8 EV9 EV10 EV11 EV12 EV1TO EV2TO %/EV/ CALCULATED VALUES
%       global EH1 EH2 EH3 EH4 EH5 EH6 EH7 EH8 EH9 EH10 EH11 EH12 EH1TO EH2TO %/EH/ CALCULATED VALUES
%       global HP1 HP2 HP3 HP4 HP5 HP6 HP7 HP8 HP1TO HP2TO %/HP/
%       global QEV QH QEH %/SCALE/ INPUT VALUES
%       global MAXLIM1 KEY1 MAXLIM2 KEY2 %/INTEG/, KEY1 and KEY2 are INPUT VALUES

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
%     FIELDS CALCULATION 
%     =============================================
%	TE: time step
    TE=0.0;
    
    for N=1:1:NBE
          if (MAXLIM1>=MX1) 
              MX1=MAXLIM1;
          end
          if (MAXLIM2>=MX2) 
              MX2=MAXLIM2;
          end

          if (MPROBLEM == 1) 
            disp(strcat('N=',num2str(N),' NC=',num2str(NC)))
          end
    
          MAXLIM1=0;
          MAXLIM2=0;
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
        HC=HR+V*TE;
        [CURC(N),DCURC(N)]=CURRENT(TE,HC,CURC(N),DCURC(N));
        ifRONT=CURC(N);

% Now the fields calculation

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
        [CURR1(N),DCURR1(N)]=CURRENT(TT,Z1max,CURR1(N),DCURR1(N));
        ifRONT=CURR1(N);
        FV=1/(1/V-(ZS-Z1max)/C/R1);
        EV1TO=-ifRONT*R*R*FV/(4*PI*EPS*C*C*R1^3);
        HP1TO=ifRONT*R*FV/(4*PI*C*R1^2);
        EH1TO=ifRONT*R*(ZS-Z1max)*FV/(4*PI*EPS*C*C*R1^3);


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
        [CURR1(N),DCURR1(N)]=CURRENT(TT,Z1max,CURR1(N),DCURR1(N));
        ifRONT=CURR1(N);
        FV=1/(1/V-(ZS-Z1max)/C/R1);
        EV2TO=-ifRONT*R*R*FV/(4*PI*EPS*C*C*R1^3);
        HP2TO=ifRONT*R*FV/(4*PI*C*R1^2);
        EH2TO=-ifRONT*R*(ZS-Z1max)*FV/(4*PI*EPS*C*C*R1^3);


        HZONTAL(TE) %call to horizontal function
        EHIND(N)=EH1+EH4+EH7+EH10;
        EHRAD(N)=EH2+EH5+EH8+EH11;
        EH_TO(N)=EH1TO+EH2TO;
        EHEL(N)=EH3+EH6+EH9+EH12;
        EHO(N)=EHIND(N)+EHRAD(N)+EHEL(N);
    %     =============================================================
    %	OUTPUT VARIABLES OF ROUTINE HZONTAL
    %     EH1  INDUCTION  COMPONENT OF Er DUE TO THE TOWER FOR ZS>0 
    %     EH2  RADIAL	    COMPONENT OF Er DUE TO THE TOWER FOR ZS>0
    %     EH3  ELECTROST. COMPONENT OF Er DUE TO THE TOWER FOR ZS>0   

    %     EH4  INDUCTION  COMPONENT OF Er DUE TO THE CHANNEL FOR ZS>0 
    %     EH5  RADIAL	    COMPONENT OF Er DUE TO THE CHANNEL FOR ZS>0
    %     EH6  ELECTROST. COMPONENT OF Er DUE TO THE CHANNEL FOR ZS>0      

    %     EH7  INDUCTION  COMPONENT OF Er DUE TO THE TOWER FOR ZS<0 
    %     EH8  RADIAL	    COMPONENT OF Er DUE TO THE TOWER FOR ZS<0
    %     EH9  ELECTROST. COMPONENT OF Er DUE TO THE TOWER FOR ZS<0   

    %     EH10  INDUCTION  COMPONENT OF Er DUE TO THE CHANNEL FOR ZS<0 
    %     EH11  RADIAL	 COMPONENT OF Er DUE TO THE CHANNEL FOR ZS<0
    %     EH12  ELECTROST. COMPONENT OF Er DUE TO THE CHANNEL FOR ZS<0      
    %     =============================================================

        VERTICAL(TE) %call to vertical function  
        EIND(N)=EV1+EV4+EV7+EV10;
        ERAD(N)=EV2+EV5+EV8+EV11;
        EV_TO(N)=EV1TO+EV2TO; %check, there is a problem with the BG-TCS models
        EELtower(N)=EV3+EV9;
        EEL(N)=EV3+EV6+EV9+EV12;
        EVtower(N)=EV1+EV2+EV3+EV7+EV8+EV9;
        EVO(N)=EIND(N)+ERAD(N)+EEL(N);

        HIND(N)=(HP1+HP3+HP5+HP7)*(-1);
        HRAD(N)=(HP2+HP4+HP6+HP8)*(-1);
        HP_TO(N)=HP1TO+HP2TO;
        HPtower(N)=(HP1+HP2+HP5+HP6)*(-1);
        HPHI(N)=HIND(N)+HRAD(N);

    %	=============================================================
    %	OUTPUT VARIABLES OF ROUTINE VERTICAL
    %
    %     EV1  INDUCTION  COMPONENT OF Ez DUE TO THE TOWER FOR ZS>0 
    %     EV2  RADIAL	    COMPONENT OF Ez DUE TO THE TOWER FOR ZS>0
    %     EV3  ELECTROST. COMPONENT OF Ez DUE TO THE TOWER FOR ZS>0   

    %     EV4  INDUCTION  COMPONENT OF Ez DUE TO THE CHANNEL FOR ZS>0 
    %     EV5  RADIAL	    COMPONENT OF Ez DUE TO THE CHANNEL FOR ZS>0
    %     EV6  ELECTROST. COMPONENT OF Ez DUE TO THE CHANNEL FOR ZS>0      

    %     HP1  INDUCTION  COMPONENT OF H DUE TO THE TOWER FOR ZS>0  
    %     HP2  RADIAL	    COMPONENT OF H DUE TO THE TOWER FOR ZS>0  
    %     HP3  INDUCTION  COMPONENT OF H DUE TO THE CHANNEL FOR ZS>0  
    %     HP4  RADIAL	    COMPONENT OF H DUE TO THE CHANNEL FOR ZS>0  
    %------------------------------------------------------------------
    %     EV7  INDUCTION  COMPONENT OF Ez DUE TO THE TOWER FOR ZS<0 
    %     EV8  RADIAL	    COMPONENT OF Ez DUE TO THE TOWER FOR ZS<0
    %     EV9  ELECTROST. COMPONENT OF Ez DUE TO THE TOWER FOR ZS<0   

    %     EV10  INDUCTION  COMPONENT OF Ez DUE TO THE CHANNEL FOR ZS<0 
    %     EV11  RADIAL	 COMPONENT OF Ez DUE TO THE CHANNEL FOR ZS<0
    %     EV12  ELECTROST. COMPONENT OF Ez DUE TO THE CHANNEL FOR ZS<0      

    %     HP5  INDUCTION  COMPONENT OF H DUE TO THE TOWER FOR ZS<0  
    %     HP6  RADIAL	    COMPONENT OF H DUE TO THE TOWER FOR ZS<0  
    %     HP7  INDUCTION  COMPONENT OF H DUE TO THE CHANNEL FOR ZS<0  
    %     HP8  RADIAL	    COMPONENT OF H DUE TO THE CHANNEL FOR ZS<0  
    %==================================================================      

 % End of fields calculation
 
        RZ=0;
        [CURT(N),DCURT(N)]=CURRENT(TE,HAM,CURT(N),DCURT(N));
        [CURB(N),DCURB(N)]=CURRENT(TE,HAB,CURB(N),DCURB(N));


        TE=TE+DELTATE;
    end
 

%     ==================================================================
%     CORRECTIVE TERM FOR HORIZONTAL FIELD (APPROACH Rubinstein)
%     ==================================================================
	ECO = zeros(NBE,1);
	ECO(1)=0.0;
	ECO(2)=HPHI(2)/TIMEE(2)*FF(TIMEE(2),SIG,EPSR);
      for N=3:NBE
          for K=2:N-1
              MK=(HPHI(K+1)-HPHI(K))/(TIMEE(K+1)-TIMEE(K));
              MK1=(HPHI(K)-HPHI(K-1))/(TIMEE(K)-TIMEE(K-1));
              ECO(N)=ECO(N)+(MK-MK1)*FF(TIMEE(N)-TIMEE(K),SIG,EPSR);
          end
          ECO(N)=ECO(N)+HPHI(2)/TIMEE(2)*FF(TIMEE(N),SIG,EPSR);
      end

%     ===================================================================
%     WRITING IN THE OUTPUT FILE : RESULT.DAT
%     ===================================================================
%       OPEN (11, FILE ='result.dat')
%       WRITE(11,500)'T','(-)EV','EVEL','EVIN','EVRD','EV_TO','EHO',
%      1'EHF','EHEL','EHIN','EHRD','EH_TO','HPHI','HPIND','HPRAD',
%      2'HP_TO','(-)DEV','DHPHI','I0(T)','IT(T)','IM(T)',
%      3'IB(T)','IFRONT(T)','DI0(T)','DIT(T)','DIM(T)','DIB(T)',
%      4'EVtower(T)','HPtower(T)','EELtower(T)'
% 
% 500   FORMAT(30(A11,1X))

      for N=2:NBE %had to change index to make it compatible with fortran
          if (N==1)
              DEVO(N)=0.0;
              DHPHI(N)=0.0;
          else
              DEVO(N)=(EVO(N)-EVO(N-1))/(TIMEE(N)-TIMEE(N-1));
              DHPHI(N)=(HPHI(N)-HPHI(N-1))/(TIMEE(N)-TIMEE(N-1));
              EHF(N)=EHO(N)+ECO(N); 
          end
          
      end
% 
%         T=TM*F3; EV=-EVO(N)*F1; EVEL=-EEL(N)*F1; EVIN=-EIND(N)*F1; EVRD=-ERAD(N)*F1; EV_TO=1-EV_TO(N)*F1; EHO=EHO(N)*F1;
%         EHF=EHF*F1; EHEL=EHEL(N)*F1; EHIN=EHIND(N)*F1; EHRD=EHRAD(N)*F1; EH_TO=EH_TO(N)*F1; HPHI=-HPHI(N)*F4; HPIND=-HIND(N)*F4;HPRAD=-HRAD(N)*F4;
%         HP_TO=HP_TO(N)*F4; DEV=-DEVO*F1/F3; DHPHI=DHPHI*F4/F3; I0_T=IPULSE(TM)*F2; IT_T=CURT(N)*F2; IM_T=CURC(N)*F2; IB_T=CURB(N)*F2; IFRONT_T= CURR1(N)*F2;
%         DI0_T=DIPULSE(TM)*F2/F3; DIT_T=DCURT(N)*F2/F3; DIM_T=DCURC(N)*F2/F3; DIB_T=DCURB(N)*F2/F3;
%         EVtower_T=6-EVtower(N)*F1;HPtower_T=-HPtower(N)*F4; EELtower_T=-EELtower(N)*F1;
% 
         TM=TIMEE;
         I0_T=IPULSE(TM)*F2;
         DI0_T=DIPULSE(TM)*F2/F3;

        T=TM*F3; EV=-EVO*F1; EVEL=-EEL*F1; EVIN=-EIND*F1; EVRD=-ERAD*F1; EV_TO=-EV_TO*F1; EHO=EHO*F1;
        EHF=EHF*F1; EHEL=EHEL*F1; EHIN=EHIND*F1; EHRD=EHRAD*F1; EH_TO=EH_TO*F1; HPHI=-HPHI*F4; HPIND=-HIND*F4;HPRAD=-HRAD*F4;
        HP_TO=HP_TO*F4; DEV=-DEVO*F1/F3; DHPHI=DHPHI*F4/F3; IT_T=CURT*F2; IM_T=CURC*F2; IB_T=CURB*F2; IFRONT_T= CURR1*F2;
        DIT_T=DCURT*F2/F3; DIM_T=DCURC*F2/F3; DIB_T=DCURB*F2/F3;
        EVtower_T=-EVtower*F1;HPtower_T=-HPtower*F4; EELtower_T=-EELtower*F1;

        % Now do the same as the inside of read_out_file
        I1_T=IT_T;
        I2_T=IB_T;
        
        EV_TOTAL = EV + EV_TO; %Vertical Eletric field
        HPHI_TOTAL = HPHI + HP_TO; 
        EH_TOTAL = EHF; 
        
%current at height QUO at time TY and its derivative    
function [CUR_out,DCUR_out]=CURRENT(TY,QUO,CUR,DCUR)

  CUR=0;DCUR=0;
  B1=0;B2=0;B3=0; 
  J=0;
      
   global RZ
   global HR RR RG ALPHA_TOWER NC 
   global VC
   global C
   global MODEL
      
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
            end
            B2=B1-2.*QUO/C;
%             B2=B2+1; %I had to add this to make compatible with 1-indexing
            if (B2>=0) 
               CUR=CUR+IPULSE(B2)*(RR^J)*(RG^(J+1))*exp(-(J*HR+QUO)*ALPHA_TOWER); %I had to subs -1 to B2 this to make compatible with 1-indexing
               DCUR=DCUR+DIPULSE(B2)*(RR^J)*(RG^(J+1))*exp(-(J*HR+QUO)*ALPHA_TOWER); %I had to subs -1 to B2 this to make compatible with 1-indexing
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
        end
% C===========2
        B1=TY-RZ/C-(QUO-HR)/C;
%         B1=B1+1; %I had to add this to make compatible with 1-indexing
	    if (B1>=0) 
             CUR=CUR+IPULSE(B1)*(-RR);
             DCUR=DCUR+DIPULSE(B1)*(-RR);
        end
% C===========3
        for J=0:1:NC
          B1=TY-RZ/C-(QUO+HR)/C-2.*J*HR/C;
%           B1=B1+1; %I had to add this to make compatible with 1-indexing
              if (B1>=0)
                 CUR=CUR+(1-RR)*(1+RR)*(RR^J)*(RG^(J+1))*IPULSE(B1)*exp(-(J+1)*HR*ALPHA_TOWER);
                 DCUR=DCUR+(1-RR)*(1+RR)*(RR^J)*(RG^(J+1))*DIPULSE(B1)*exp(-(J+1)*HR*ALPHA_TOWER);
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
% 	  REAL*8 MU TAU22  ALPHA BETA I01 GAMMA DELTA I02 I0B T0B
%       REAL*8 FK FKK TA
%       INTEGER NW2,NW1
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

function DIPULSE_calc=DIPULSE(TA)
%===============================d(2Heidler+2exp)/dt
% 	REAL*8 MU,%,PI,EPS,V,VC,W,H,ZLAM,IW,ETA,TAU1,TAU2,IW2,ETA2,TAU12
%      1       ,TAU22,ALPHA,BETA,I01,GAMMA,DELTA,I02,I0B,T0B
%       REAL*8 FK1,FK2,FKK1,FKK2,TA
%       INTEGER NW2,NW1

      global IW1 ETA TAU11 TAU21 IW2 ETA2 TAU12 TAU22 NW2 NW1
      global ALPHA BETA I01 GAMMA DELTA I02
	  global I0B T0B


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

function INTPULSE_calc=INTPULSE(A,B)
%===============================integrale(2Heidler+2exp)dt 0-->t
	
%       REAL*8 IW,ETA,TAU1,TAU2,IW2,ETA2,TAU12,TAU22
%       REAL*8 ALPHA,BETA,I01,GAMMA,DELTA,I02,I0B,T0B
%       REAL*8 V,VC,W,H,ZLAM
%       REAL*8 SUM,SUM1,SUM2
%       REAL*8 PULSE,A,B,HS,E1,TA
%       INTEGER NN,K,K1,K2,NW2,NW1,NSIMP,NN1,NN2,K12,K22

      global IW1 ETA TAU11 TAU21 IW2 ETA2 TAU12 TAU22 NW2 NW1
      global ALPHA BETA I01 GAMMA DELTA I02
%       global V VC W H ZLAM
      global NSIMP
	  global I0B T0B
% an anonymous function is defined here
      PULSE = @(TA) (IW1/ETA).*(TA./TAU11).^NW1./(1+((TA./TAU11).^NW1)).*exp(-TA./TAU21)+(IW2/ETA2).*(TA./TAU12).^NW2./(1+((TA./TAU12).^NW2)).*exp(-TA./TAU22);
%original      PULSE(TA)=(IW1/ETA)*(TA/TAU11)^NW1/(1+((TA/TAU11)^NW1))*exp(-TA/TAU21)+(IW2/ETA2)*(TA/TAU12)^NW2/(1+((TA/TAU12)^NW2))*exp(-TA/TAU22);
      HS=(B-A)/double(NSIMP);
      NN=NSIMP/2-1;
      SUM=PULSE(A)+4*PULSE(A+(NSIMP-1)*HS);
      SUM1=0;
      SUM2=0;
      for K=1:1:NN
          K1=2*K;
          K2=2*K-1;
          SUM1=SUM1+2*PULSE(A+(double(K1))*HS);
          SUM2=SUM2+4*PULSE(A+(K2)*HS);
      end
      SUM=SUM+SUM1+SUM2+PULSE(B);	

      INTPULSE_calc=HS.*SUM/3;
      INTPULSE_calc=INTPULSE_calc+I01*((exp(-ALPHA*A)-exp(-ALPHA*B))./ALPHA-(exp(-BETA*A)-exp(-BETA*B))./BETA)+I02*((exp(-GAMMA*A)-exp(-GAMMA*B))./GAMMA-(exp(-DELTA*A)-exp(-DELTA*B))./DELTA);

	if any(and((A<T0B),(B<=T0B))) 
        INTPULSE_calc=INTPULSE_calc+0.5*I0B*(B.^2-A.^2)./T0B;
    elseif any(and((A<=T0B),(B>T0B))) 
        INTPULSE_calc=INTPULSE_calc+0.5*I0B*(T0B^2-A.^2)./T0B+I0B*(B-T0B);
    elseif any(and((A>=T0B),(B>=T0B))) 
        INTPULSE_calc=INTPULSE_calc+I0B*(B-A);
	end


function FF_calc=FF(T,SIGMA,EPSILR)
%     ========================================
%       FORMULA FOR CORRECTIVE TERM FOR HORIZONTAL FIELD (APPROACH Rubinstein)

    global MU EPS
      
    P=SIGMA/(2.*EPS*EPSILR);
    Z0=sqrt(MU/EPS);
    X=P*T;
    Y=X/3.75;
%      GOTO 444 old formula RUBINSTEIN
	if (abs(X)<3.75)

        S0=1+3.5156229*(Y^2)+3.0899424*(Y^4)+1.2067492*(Y^6)+0.2659732*(Y^8)+0.0360768*(Y^10)+0.0045813*(Y^12);
        S1=X*(0.5+0.87890594*Y*Y+0.51498869*(Y^4)+0.15084934*(Y^6)+0.02658733*(Y^8)+0.00301532*(Y^10)+0.00032411*(Y^12));
        FF_calc=Z0*T*exp(-X)*(S0+S1)/sqrt(EPSILR);
    else
        S0=(0.39894228+0.01328592/Y+0.00225319/(Y^2)-0.00157565/(Y^3)+0.00916281/(Y^4)-0.02057706/(Y^5)+0.02635537/(Y^6)-0.01647633/(Y^7)+0.00392377/(Y^8))/sqrt(X);
        S1=(0.39894228-0.03988024/Y-0.00362018/(Y*Y)+0.00163801/(Y*Y*Y)-.01031555/(Y^4)+0.02282967/(Y^5)-0.02895312/(Y^6)+0.01787654/(Y^7)-0.00420059/(Y^8))/sqrt(X);
        FF_calc=Z0*T*(S0+S1)/sqrt(EPSILR);
    end
% C444   FF=2.0*T*sqrt(MU/(PI*SIGMA*T))
	

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
	

%     ==========================================================
%     CALCULATION OF VERTICAL ELECTRIC FIELD AND MAGNETIC FIELD
%     ==========================================================
function VERTICAL(TE)

%   REAL*8 MU
%   REAL*8 TT,D01AHF
%   REAL*8 HR,RR,RG,ALPHA_TOWER,ERR,ERR2
%   REAL*8 EV1,EV2,EV3,EV4,EV5,EV6,EV7,EV8,EV9,EV10,EV11,EV12
%   REAL*8 QEV,QH,QEH
%   REAL*8 HP1,HP2,HP3,HP4,HP5,HP6,HP7,HP8
%   REAL*8 EV1TO,EV2TO,HP1TO,HP2TO
%   REAL*8 PHI,PHI1,PHI2,DELTAP
%   REAL*8 TE,T01,H0,H1
%   REAL*8 C1,C3,RELLER,ER0
%   REAL*8 EELP,EINDP,ERADP,BINDP,BRADP
%   INTEGER NC,IFAIL,NPTS,NLIMIT,MAXLIM1,KEY1,MAXLIM2,KEY2,MODEL
%   INTEGER MPROBLEM

  NTPS=0; 
  
  global  C PI EPS 
  global  V  H 
  global  R ZS 
  global  TT
  global  HR 
  global  ERR 
  global  EV1 EV2 EV3 EV4 EV5 EV6 EV7 EV8 EV9 EV10 EV11 EV12
  global  HP1 HP2 HP3 HP4 HP5 HP6 HP7 HP8
  global  QEV QH 
  global  MAXLIM1 KEY1 

%   EXTERNAL EELP,EINDP,ERADP,BINDP,BRADP
%   EXTERNAL D01AHF

  ZS=abs(ZS);
  T01=sqrt(R*R+(HR-ZS)*(HR-ZS))/C;

%===========T01: OBSERVATION TIME DELAY
  TT=TE+T01;
%===========TT: ACTUAL TIME

  if (TT==T01)
      EV1=0;
      EV2=0;
      EV3=0;
      EV4=0;
      EV5=0;
      EV6=0;
      EV7=0;
      EV8=0;
      EV9=0;
      EV10=0;
      EV11=0;
      EV12=0;
      HP1=0;
      HP2=0;
      HP3=0;
      HP4=0;
      HP5=0;
      HP6=0;
      HP7=0;
      HP8=0;
      return
  end
  C1=R/(4.*PI);
  C3=1./(4.*PI*EPS);

%===========DETERMINATION OF THE UPPER INTEGRATION LIMIT
  PHI=1./(V*V)-1./(C*C);
  PHI1=ZS/(C*C)-(TT/V)-HR/(V*V);
  PHI2=(TT+HR/V)^2.-(R/C)^2.-(ZS/C)^2;
  DELTAP=abs(sqrt(PHI1*PHI1-PHI*PHI2));
  H1=(-PHI1-DELTAP)/PHI;
  if (H1<HR) 
      H1=HR;
  end
  if (H1>H) 
      H1=H;
  end
%===========DETERMINATION OF THE LOWER INTEGRATION LIMIT
  PHI1=(ZS/(C*C)+(TT/C)-HR/(C*C))*2;
  PHI2=-(TT-HR/C)^2+(R/C)^2+(ZS/C)^2;
  H0=PHI2/PHI1;
  if (H0>HR) 
      H0=HR;
  end
  if (H0<0) 
      H0=0;
  end

  ER0=ERR;
  IFAIL=0;
  NLIMIT=KEY1;
% the function is called as: D01AHF(A,B,EPR,NPTS,RELERR,F,NL,IFAIL)
%THIS FUNCTION ROUTINE PERFORMS AUTOMATIC INTEGRATION OVER A
%     FINITE INTERVAL USING THE BASIC INTEGRATION ALGORITHMS D01AHY
%     AND D01AHX, TOGETHER WITH, IF NECESSARY, AN ADAPTIVE
%     SUBDIVISION PROCESS.
%
%     INPUT ARGUMENTS
%     ----- ----------
%     A,B     -  LOWER AND UPPER INTEGRATION LIMITS.
%     EPR     -  REQUIRED RELATIVE ACCURACY.
%     NL      -  APPROXIMATE LIMIT ON NUMBER OF INTEGRAND
%                EVALUATIONS. IF SET NEGATIVE OR ZERO THE
%                DEFAULT IS 10000.
%     F       -  THE USER NAMED AND PREPARED FUNCTION  F(X)
%                GIVES THE VALUE OF THE INTEGRAND AT X.
%     IFAIL      INTEGER VARIABLE
%             - 0  FOR HARD FAIL REPORT
%             - 1  FOR SOFT FAIL REPORT
%    OUTPUT ARGUMENTS
%     ------ ----------
%     NPTS    -  NUMBER OF INTEGRAND EVALUATIONS USED IN OBTAINING
%                THE RESULT.
%     RELERR  -  ROUGH ESTIMATE OF RELATIVE ACCURACY ACHIEVED.
%     IFAIL   -  VALUE INDICATES THE OUTCOME OF THE INTEGRATION -
%                IFAIL  = 0  CONVERGED
%                IFAIL  = 1  INTEGRAND EVALUATIONS EXCEEDED  NL.
%                            THE RESULT WAS OBTAINED BY CONTINUING
%                            BUT IGNORING ANY NEED TO SUBDIVIDE.
%                            RESULT LIKELY TO BE INACCURATE.
%                IFAIL  = 2  DURING THE SUBDIVISION PROCESS
%                            THE STACK BECAME FULL
%                            (PRESENTLY SET TO HOLD 20
%                            LEVELS OF INFORMATION.  MAY BE
%                            INCREASED BY  ALTERING  ISMAX
%                            AND THE DIMENSIONS OF STACK
%                            AND ISTACK). RESULT IS
%                            OBTAINED BY CONTINUING BUT
%                            IGNORING CONVERGENCE FAILURES
%                            ON INTERVALS  WHICH CANNOT BE
%                            ACCOMMODATED ON THE STACKS.
%                            RESULT LIKELY TO BE
%                            INACCURATE.
%                IFAIL  = 3  INVALID ACCURACY REQUEST.
%
[integrated_out,NPTS,RELER,IFAIL]=D01AHF(H0,HR,ER0,@EINDP,NLIMIT,IFAIL);
%ORIGINAL CODE LINE  EV1=(C3/C)*D01AHF(H0,HR,ER0,NPTS,RELLER,EINDP,NLIMIT,IFAIL)/QEV; 
  EV1=(C3/C)*integrated_out/QEV; 
  %EINDP is the Elecric field (induction component of Ez and its a function
  %(should be defined as an abstract one in Matlab)
  if (NPTS>=MAXLIM1) %APARENTEMENTE ESTA LINEA NUNCA SIRVE PARA NIMIER
      MAXLIM1=NPTS;  %APARENTEMENTE ESTA LINEA NUNCA SIRVE PARA NIMIER
  end                %APARENTEMENTE ESTA LINEA NUNCA SIRVE PARA NIMIER, porque la variable MAXLIM1 no se utiliza

  NLIMIT=KEY1;
  [integrated_out,NPTS,RELER,IFAIL]=D01AHF(H0,HR,ER0,@ERADP,NLIMIT,IFAIL);
%   EV2=-(C3*R*R/(C*C))*D01AHF(H0,HR,ER0,NPTS,RELLER,ERADP,NLIMIT,IFAIL)/QEV;
  EV2=-(C3*R*R/(C*C))*integrated_out/QEV;
  if (NPTS>=MAXLIM1) 
      MAXLIM1=NPTS;
  end

  NLIMIT=KEY1;
%   EV3=C3*D01AHF(H0,HR,ER0,NPTS,RELLER,EELP,NLIMIT,IFAIL)/QEV;
  [integrated_out,NPTS,RELER,IFAIL]=D01AHF(H0,HR,ER0,@EELP,NLIMIT,IFAIL);
  EV3=C3*integrated_out/QEV;
  if (NPTS>=MAXLIM1) 
      MAXLIM1=NPTS;
  end

  NLIMIT=KEY1;
%   EV4=(C3/C)*D01AHF(HR,H1,ER0,NPTS,RELLER,EINDP,NLIMIT,IFAIL)/QEV;
  [integrated_out,NPTS,RELER,IFAIL]=D01AHF(HR,H1,ER0,@EINDP,NLIMIT,IFAIL);
  EV4=(C3/C)*integrated_out/QEV;
  if (NPTS>=MAXLIM1) 
      MAXLIM1=NPTS;
  end

  NLIMIT=KEY1;
%   EV5=-(C3*R*R/(C*C))*D01AHF(HR,H1,ER0,NPTS,RELLER,ERADP,NLIMIT,IFAIL)/QEV;
  [integrated_out,NPTS,RELER,IFAIL]=D01AHF(HR,H1,ER0,@ERADP,NLIMIT,IFAIL);
  EV5=-(C3*R*R/(C*C))*integrated_out/QEV;
  if (NPTS>=MAXLIM1) 
      MAXLIM1=NPTS;
  end

  NLIMIT=KEY1;
%   EV6=C3*D01AHF(HR,H1,ER0,NPTS,RELLER,EELP,NLIMIT,IFAIL)/QEV;
  [integrated_out,NPTS,RELER,IFAIL]=D01AHF(HR,H1,ER0,@EELP,NLIMIT,IFAIL);
  EV6=C3*integrated_out/QEV;
  if (NPTS>=MAXLIM1) 
      MAXLIM1=NPTS;
  end

  NLIMIT=KEY1;
%   HP1=(C1)*D01AHF(H0,HR,ER0,NPTS,RELLER,BINDP,NLIMIT,IFAIL)/QH;
  [integrated_out,NPTS,RELER,IFAIL]=D01AHF(H0,HR,ER0,@BINDP,NLIMIT,IFAIL);
  HP1=(C1)*integrated_out/QH;
  if (NPTS>=MAXLIM1) 
      MAXLIM1=NPTS;
  end

  NLIMIT=KEY1;
%   HP2=(C1/C)*D01AHF(H0,HR,ER0,NPTS,RELLER,BRADP,NLIMIT,IFAIL)/QH;
  [integrated_out,NPTS,RELER,IFAIL]=D01AHF(H0,HR,ER0,@BRADP,NLIMIT,IFAIL);
  HP2=(C1/C)*integrated_out/QH;
  if (NPTS>=MAXLIM1) 
      MAXLIM1=NPTS;
  end

  NLIMIT=KEY1;
%   HP3=(C1)*D01AHF(HR,H1,ER0,NPTS,RELLER,BINDP,NLIMIT,IFAIL)/QH;
  [integrated_out,NPTS,RELER,IFAIL]=D01AHF(HR,H1,ER0,@BINDP,NLIMIT,IFAIL);
  HP3=(C1)*integrated_out/QH;
  if (NPTS>=MAXLIM1) 
      MAXLIM1=NPTS;
  end

  NLIMIT=KEY1;
%   HP4=(C1/C)*D01AHF(HR,H1,ER0,NPTS,RELLER,BRADP,NLIMIT,IFAIL)/QH;
  [integrated_out,NPTS,RELER,IFAIL]=D01AHF(HR,H1,ER0,@BRADP,NLIMIT,IFAIL);
  HP4=(C1/C)*integrated_out/QH;
  if (NPTS>=MAXLIM1) 
      MAXLIM1=NPTS;
  end


%==========DETERMINATION OF THE UPPER INTEGRATION LIMIT
  ZS=-abs(ZS);
  PHI=1./(V*V)-1./(C*C);
  PHI1=ZS/(C*C)-(TT/V)-HR/(V*V);
  PHI2=(TT+HR/V)^2.-(R/C)^2.-(ZS/C)^2;
  DELTAP=abs(sqrt(PHI1*PHI1-PHI*PHI2));
  H1=(-PHI1-DELTAP)/PHI;
  if (H1<HR) 
      H1=HR;
  end
  if (H1>H) 
      H1=H;
  end

%==========DETERMINATION OF THE LOWER INTEGRATION LIMIT
  PHI1=(ZS/(C*C)+(TT/C)-HR/(C*C))*2;
  PHI2=-(TT-HR/C)^2+(R/C)^2+(ZS/C)^2;
  H0=PHI2/PHI1;
  if (H0>HR) 
      H0=HR;
  end
  if (H0<0) 
      H0=0;
  end

  NLIMIT=KEY1;
%   EV7=(C3/C)*D01AHF(H0,HR,ER0,NPTS,RELLER,EINDP,NLIMIT,IFAIL)/QEV;
  [integrated_out,NPTS,RELER,IFAIL]=D01AHF(H0,HR,ER0,@EINDP,NLIMIT,IFAIL);
  EV7=(C3/C)*integrated_out/QEV;
  if (NPTS>=MAXLIM1) 
      MAXLIM1=NPTS;
  end

  NLIMIT=KEY1;
%   EV8=-(C3*R*R/(C*C))*D01AHF(H0,HR,ER0,NPTS,RELLER,ERADP,NLIMIT,IFAIL)/QEV;
  [integrated_out,NPTS,RELER,IFAIL]=D01AHF(H0,HR,ER0,@ERADP,NLIMIT,IFAIL);
  EV8=-(C3*R*R/(C*C))*integrated_out/QEV;
  if (NPTS>=MAXLIM1) 
      MAXLIM1=NPTS;
  end

  NLIMIT=KEY1;
%   EV9=C3*D01AHF(H0,HR,ER0,NPTS,RELLER,EELP,NLIMIT,IFAIL)/QEV;
  [integrated_out,NPTS,RELER,IFAIL]=D01AHF(H0,HR,ER0,@EELP,NLIMIT,IFAIL);
  EV9=C3*integrated_out/QEV;
  if (NPTS>=MAXLIM1) 
      MAXLIM1=NPTS;
  end

  NLIMIT=KEY1;
%   EV10=(C3/C)*D01AHF(HR,H1,ER0,NPTS,RELLER,EINDP,NLIMIT,IFAIL)/QEV;
  [integrated_out,NPTS,RELER,IFAIL]=D01AHF(HR,H1,ER0,@EINDP,NLIMIT,IFAIL);
  EV10=(C3/C)*integrated_out/QEV;
  if (NPTS>=MAXLIM1) 
      MAXLIM1=NPTS;
  end

  NLIMIT=KEY1;
%   EV11=-(C3*R*R/(C*C))*D01AHF(HR,H1,ER0,NPTS,RELLER,ERADP,NLIMIT,IFAIL)/QEV;
  [integrated_out,NPTS,RELER,IFAIL]=D01AHF(HR,H1,ER0,@ERADP,NLIMIT,IFAIL);
  EV11=-(C3*R*R/(C*C))*integrated_out/QEV;
  if (NPTS>=MAXLIM1) 
      MAXLIM1=NPTS;
  end

  NLIMIT=KEY1;
%   EV12=C3*D01AHF(HR,H1,ER0,NPTS,RELLER,EELP,NLIMIT,IFAIL)/QEV;
  [integrated_out,NPTS,RELER,IFAIL]=D01AHF(HR,H1,ER0,@EELP,NLIMIT,IFAIL);
  EV12=C3*integrated_out/QEV;
  if (NPTS>=MAXLIM1) 
      MAXLIM1=NPTS;
  end

  NLIMIT=KEY1;
%   HP5=(C1)*D01AHF(H0,HR,ER0,NPTS,RELLER,BINDP,NLIMIT,IFAIL)/QH;
  [integrated_out,NPTS,RELER,IFAIL]=D01AHF(H0,HR,ER0,@BINDP,NLIMIT,IFAIL);
  HP5=(C1)*integrated_out/QH;
  if (NPTS>=MAXLIM1) 
      MAXLIM1=NPTS;
  end

  NLIMIT=KEY1;
%   HP6=(C1/C)*D01AHF(H0,HR,ER0,NPTS,RELLER,BRADP,NLIMIT,IFAIL)/QH;
  [integrated_out,NPTS,RELER,IFAIL]=D01AHF(H0,HR,ER0,@BRADP,NLIMIT,IFAIL); 
  HP6=(C1/C)*integrated_out/QH;
  if (NPTS>=MAXLIM1) 
      MAXLIM1=NPTS;
  end

  NLIMIT=KEY1;
%   HP7=(C1)*D01AHF(HR,H1,ER0,NPTS,RELLER,BINDP,NLIMIT,IFAIL)/QH;
  [integrated_out,NPTS,RELER,IFAIL]=D01AHF(HR,H1,ER0,@BINDP,NLIMIT,IFAIL);
  HP7=(C1)*integrated_out/QH;
  if (NPTS>=MAXLIM1) 
      MAXLIM1=NPTS;
  end

  NLIMIT=KEY1;
%   HP8=(C1/C)*D01AHF(HR,H1,ER0,NPTS,RELLER,BRADP,NLIMIT,IFAIL)/QH;
  [integrated_out,NPTS,RELER,IFAIL]=D01AHF(HR,H1,ER0,@BRADP,NLIMIT,IFAIL);
  HP8=(C1/C)*integrated_out/QH;
  if (NPTS>=MAXLIM1) 
      MAXLIM1=NPTS;
  end
%====================================================================

ZS=abs(ZS);


%	============================= ELECTROSTATIC COMPONENT OF Ez =========
function y=EELP(x)
%	=============================
%   IMPLICIT NONE
% 	REAL*8 INTPULSE,P
%   REAL*8 MU,%,PI,EPS,V,VC,W,H,ZLAM,R,ZS,TT,HR,RR,RG,ALPHA_TOWER
%   REAL*8 ICI,ICII,x,RZ,QEV,QH,QEH
% 	REAL*8 A1,A2,A3,B1,B2,B3
% 	INTEGER NC,J,MODEL,MPROBLEM
%   EXTERNAL INTPULSE,P

      global  C
      global  V VC
      global  R ZS
      global  TT
      global  HR RR RG ALPHA_TOWER NC
      global  QEV

	  RZ=sqrt((R*R)+((x-ZS).^2));
      ICI=0;
      ICII=0;
      y=0; %add to eliminate the Matlab message
	if any(x<HR) 
%===========1	
        for J=0:1:NC		
           B1=TT-RZ/C-(HR-x)/C-2.*J*HR/C;
              if any(B1>=0)
               A1=0;
               ICII=INTPULSE(A1,B1);
               ICI=ICI+ICII.*(RR.^J)*(RG.^J).*exp(-J.*HR.*ALPHA_TOWER);
              end
              B2=B1-2.*x/C;
              if any(B2>=0)
               A2=0;
               ICII=INTPULSE(A2,B2);
               ICI=ICI+ICII.*(RR.^J)*(RG.^(J+1)).*exp(-(J.*HR+x).*ALPHA_TOWER);
              end
        end
	ICI=ICI.*(1-RR);
   else
%===========x GREATER THAN HR
	B1=TT-RZ/C-(x-HR)/V;
	    if any(B1>=0)
	      A3=(x-HR)/V-(x-HR)/VC;
	      B3=TT-RZ/C-(x-HR)/VC;
	      ICII=INTPULSE(A3,B3);
		  ICI=ICI+ICII.*P(x);
        end
%===========2
	B2=TT-RZ/C-(x-HR)/C;
	    if any(B2>=0)
	     A1=0;
		 ICII=INTPULSE(A1,B2);
		 ICI=ICI+ICII*(-RR);
        end
%===========3
	   for J=0:1:NC
          B3=TT-RZ/C-(x+HR)/C-2.*J*HR/C;
          if any(B3>=0)
	         A1=0;
		     ICII=INTPULSE(A1,B3);
		     ICI=ICI+ICII.*(1-RR).*(1+RR).*(RR^J)*(RG^(J+1)).*exp(-(J+1)*HR.*ALPHA_TOWER);
          end
        end
%================================	
	
    end
y=(2.*(ZS-x).*(ZS-x)-R*R).*ICI./(RZ.^5)*QEV;
    
    

%	============================= INDUCTION COMPONENT OF Ez ==========
function y=EINDP(x)
%	=============================
%       IMPLICIT NONE
% 	    REAL*8 IPULSE,P
%       REAL*8 MU,%,PI,EPS,V,VC,W,H,ZLAM,R,ZS,TT,HR,RR,RG,ALPHA_TOWER
%       REAL*8 ICI,ICII,x,RZ,QEV,QH,QEH
% 	    REAL*8 B1,B2,B3
%       INTEGER NC,J,MODEL,MPROBLEM
%       EXTERNAL IPULSE,P

      global  C
      global  V VC
      global  R ZS
      global  TT
      global  HR RR RG ALPHA_TOWER NC
      global  QEV

      RZ=sqrt((R*R)+((x-ZS).^2));
      ICI=0;
      ICII=0;
      y=0; %added to eliminate Matlab message
	if any(x<HR) %any act on a vector
%===========1	
	 for J=0:1:NC
       B1=TT-RZ/C-(HR-x)/C-2.*J*HR/C;
          if any(B1>=0) 
            ICII=IPULSE(B1);
            ICI=ICI+ICII*(RR^J)*(RG^J).*exp(-J*HR*ALPHA_TOWER);
          end
          B2=B1-2.*x/C;
          if any(B2>=0) 
            ICII=IPULSE(B2);
            ICI=ICI+ICII*(RR^J)*(RG^(J+1)).*exp(-(J*HR+x)*ALPHA_TOWER);
          end
     end
	ICI=ICI*(1-RR);
    else
%===========x GREATER THAN HR
	B1=TT-RZ/C-(x-HR)/V;
	B3=TT-RZ/C-(x-HR)/VC;
	    if any(B1>=0) 
	    ICII=IPULSE(B3);
		ICI=ICI+ICII.*P(x);
        end
%===========2
	    B1=TT-RZ/C-(x-HR)/C;
	    if any(B1>=0) 
		 ICII=IPULSE(B1);
		 ICI=ICI+ICII*(-RR);
        end
%===========3
	   for J=0:1:NC
          B1=TT-RZ/C-(x+HR)/C-2.*J*HR/C;
          if any(B1>=0) 
            ICII=IPULSE(B1);
            ICI=ICI+ICII*(1-RR)*(1+RR)*(RR^J)*(RG^(J+1))*exp(-(J+1)*HR*ALPHA_TOWER);
          end
       end
%================================	
	
    end
y=ICI.*(2*(ZS-x).*(ZS-x)-R*R)./(RZ.^4)*QEV;
%     

%	============================= RADIATION COMPONENT OF Ez ===============
function y=ERADP(x)
%	=============================
%       IMPLICIT NONE
% 	    REAL*8 DIPULSE,P
%       REAL*8 MU,%,PI,EPS,V,VC,W,H,ZLAM,R,ZS,TT,HR,RR,RG,ALPHA_TOWER
%       REAL*8 ICI,ICII,x,RZ,QEV,QH,QEH
% 	    REAL*8 B1,B2,B3
% 	    INTEGER NC,J,MODEL,MPROBLEM
%       EXTERNAL DIPULSE,P

      global  C
      global  V VC
      global  R ZS
      global  TT
      global  HR RR RG ALPHA_TOWER NC
      global  QEV

      RZ=sqrt((R*R)+((x-ZS).^2));
      ICI=0;
      ICII=0;
      y=0; %added to eliminate Matlab message
	if any(x<HR) 
%===========1	
	 for J=0:1:NC
          B1=TT-RZ/C-(HR-x)/C-2.*J*HR/C;
          if any(B1>=0)
          ICII=DIPULSE(B1);
          ICI=ICI+ICII.*(RR.^J).*(RG.^J).*exp(-J.*HR.*ALPHA_TOWER);
          end
          B2=B1-2.*x/C;
          if (B2>=0)
          ICII=DIPULSE(B2);
          ICI=ICI+ICII.*(RR.^J).*(RG.^(J+1)).*exp(-(J.*HR+x).*ALPHA_TOWER);
          end
       end
	ICI=ICI.*(1-RR);
    else
%===========x GREATER THAN HR
	B1=TT-RZ/C-(x-HR)/V;
	B3=TT-RZ/C-(x-HR)/VC;
	    if any(B1>=0)
	    ICII=DIPULSE(B3);
		ICI=ICI+ICII.*P(x);
		end
%===========2
	    B1=TT-RZ/C-(x-HR)/C;
	    if any(B1>=0) 
		 ICII=DIPULSE(B1);
		 ICI=ICI+ICII*(-RR);
		end
%===========3
	    for J=0:1:NC
            B1=TT-RZ/C-(x+HR)/C-2.*J*HR/C;
            if any(B1>=0)
              ICII=DIPULSE(B1);
              ICI=ICI+ICII*(1-RR)*(1+RR)*(RR^J)*(RG^(J+1))*exp(-(J+1)*HR*ALPHA_TOWER);
            end
        end
%================================	
	
	end

      y=ICI./(RZ.^3).*QEV;

%	============================= INDUCTION COMPONENT OF H ===========
function y=BINDP(x)
%	=============================
%   IMPLICIT NONE
% 	REAL*8 IPULSE,P
%   REAL*8 MU,%,PI,EPS,V,VC,W,H,ZLAM,R,ZS,TT,HR,RR,RG,ALPHA_TOWER
%   REAL*8 ICI,ICII,x,RZ,QEV,QH,QEH
% 	REAL*8 B1,B2,B3
% 	INTEGER NC,J,MODEL,MPROBLEM
%   EXTERNAL IPULSE,P

      global  C
      global  V VC
      global  R ZS
      global  TT
      global  HR RR RG ALPHA_TOWER NC
      global  QH
%	--------------------------------------------------
	RZ=sqrt((R*R)+((x-ZS).^2));
	ICI=0;
    ICII=0;
    y=0; %added to eliminate Matlab message
	if any(x<HR) 
%===========1	
	 for J=0:1:NC
          B1=TT-RZ/C-(HR-x)/C-2.*J*HR/C;
          if any(B1>=0) 
              ICII=IPULSE(B1);
              ICI=ICI+ICII.*(RR^J).*(RG.^J).*exp(-J*HR.*ALPHA_TOWER);
          end
          B2=B1-2.*x/C;
          if any(B2>=0) 
              ICII=IPULSE(B2);
              ICI=ICI+ICII.*(RR.^J).*(RG.^(J+1)).*exp(-(J*HR+x).*ALPHA_TOWER);
          end
     end
	ICI=ICI.*(1-RR);
    else
%===========x GREATER THAN HR
        B1=TT-RZ/C-(x-HR)/V;
        B3=TT-RZ/C-(x-HR)/VC;
	    if any(B1>=0) 
            ICII=IPULSE(B3);
            ICI=ICI+ICII.*P(x);
		end
%===========2
	B1=TT-RZ/C-(x-HR)/C;
	    if any(B1>=0) 
		 ICII=IPULSE(B1);
		 ICI=ICI+ICII*(-RR);
		end
%===========3
	   for J=0:1:NC
          B1=TT-RZ/C-(x+HR)/C-2.*J*HR/C;
          if any(B1>=0) 
             ICII=IPULSE(B1);
             ICI=ICI+ICII*(1-RR)*(1+RR)*(RR^J)*(RG^(J+1))*exp(-(J+1)*HR*ALPHA_TOWER);
          end
       end
%================================	
	
	end

    y=ICI./(RZ.^3).*QH;


%	============================= RADIATION COMPONENT OF H =================
function y=BRADP(x)
%	============================
%  IMPLICIT NONE
%  REAL*8 DIPULSE,P
%  REAL*8 MU,%,PI,EPS,V,VC,W,H,ZLAM,R,ZS,TT,HR,RR,RG,ALPHA_TOWER
%  REAL*8 ICI,ICII,x,RZ,QEV,QH,QEH
%  REAL*8 B1,B2,B3
%  INTEGER NC,J,MODEL,MPROBLEM
%  EXTERNAL DIPULSE,P

      global  C
      global  V VC
      global  R ZS
      global  TT
      global  HR RR RG ALPHA_TOWER NC
      global  QH
%     ------------------------------------------------
      RZ=sqrt((R*R)+((x-ZS).^2));
	  ICI=0;
      ICII=0;
      y=0; %added to eliminate Matlab message
	if any(x<HR) 
%===========1	
	 for J=0:1:NC
       B1=TT-RZ/C-(HR-x)/C-2.*J*HR/C;
          if any(B1>=0) 
          ICII=DIPULSE(B1);
          ICI=ICI+ICII*(RR.^J).*(RG.^J).*exp(-J.*HR.*ALPHA_TOWER);
          end
       B2=B1-2.*x/C;
          if any(B2>=0) 
          ICII=DIPULSE(B2);
          ICI=ICI+ICII*(RR.^J).*(RG.^(J+1)).*exp(-(J.*HR+x).*ALPHA_TOWER);
          end
     end
	ICI=ICI.*(1-RR);
    else
%===========x GREATER THAN HR
	B1=TT-RZ/C-(x-HR)/V;
	B3=TT-RZ/C-(x-HR)/VC;
	    if any(B1>=0) 
	    ICII=DIPULSE(B3);
		ICI=ICI+ICII.*P(x);
		end
%===========2
	  B1=TT-RZ/C-(x-HR)/C;
	  if any(B1>=0) 
		 ICII=DIPULSE(B1);
		 ICI=ICI+ICII*(-RR);
	  end
%===========3
	  for J=0:1:NC
          B1=TT-RZ/C-(x+HR)/C-2.*J*HR/C;
          if any(B1>=0) 
             ICII=DIPULSE(B1);
             ICI=ICI+ICII.*(1-RR).*(1+RR).*(RR.^J).*(RG.^(J+1)).*exp(-(J+1).*HR.*ALPHA_TOWER);
          end
      end

%================================	
	end

    y=ICI./((RZ.^2)).*QH;


%     ==========================================
%     CALCULATION OF RADIAL ELECTRIC FIELD Er
%     ==========================================
function HZONTAL(TE)


% 	IMPLICIT NONE
%   REAL*8 MU,%,PI,EPS,V,VC,W,H,ZLAM,R1,FV,IFRONT
%   REAL*8 R,ZS,TT,HR,RR,RG,ALPHA_TOWER,ERR,ERR2
%   REAL*8 EH1,EH2,EH3,EH4,EH5,EH6,EH7,EH8,EH9,EH10,EH11,EH12
%   REAL*8 QEV,QH,QEH,D01AHF
%   REAL*8 TE,T01, H0,H1,C1,C3,ER0,RELLER
%   REAL*8 PHI,DELTAP,PHI1,PHI2
%   REAL*8 EELRP,EINDRP,ERADRP,EH1TO,EH2TO
%   INTEGER NPTS,NLIMIT,IFAIL,NC,MAXLIM1,KEY1,MAXLIM2,KEY2,MODEL
% 	INTEGER MPROBLEM

      global  C PI EPS
      global  V H
      global  R ZS
      global  TT
      global  HR 
      global  ERR2
      global  EH1 EH2 EH3 EH4 EH5 EH6 EH7 EH8 EH9 EH10 EH11 EH12
      global  QEH
      global  MAXLIM2 KEY2

%       EXTERNAL EELRP,EINDRP,ERADRP
%       EXTERNAL D01AHF

      ZS=abs(ZS);
      T01=sqrt(R*R+(HR-ZS)*(HR-ZS))/C;

      TT=TE+T01;

      if (TT==T01) 
          EH1=0;
          EH2=0;
          EH3=0;
          EH4=0;
          EH5=0;
          EH6=0;
          EH7=0;
          EH8=0;
          EH9=0;
          EH10=0;
          EH11=0;
          EH12=0;
          return
      end
      if (ZS==(0.0))
          EH1=0;
          EH2=0;
          EH3=0;
          EH4=0;
          EH5=0;
          EH6=0;
          EH7=0;
          EH8=0;
          EH9=0;
          EH10=0;
          EH11=0;
          EH12=0;
          return
      end
	  C1=R/(4*PI);
      C3=1./(4.*PI*EPS);
      ER0=ERR2;


%===========DETERMINATION OF THE UPPER INTEGRATION LIMIT
      PHI=1./(V*V)-1./(C*C);
      PHI1=ZS/(C*C)-(TT/V)-HR/(V*V);
      PHI2=(TT+HR/V)^2.-(R/C)^2.-(ZS/C)^2;
      DELTAP=abs(sqrt(PHI1*PHI1-PHI*PHI2));
      H1=(-PHI1-DELTAP)/PHI;
      if (H1<HR)
          H1=HR;
      end
      if (H1>H) 
          H1=H;
      end
%============DETERMINATION OF THE LOWER INTEGRATION LIMIT
      PHI1=(ZS/(C*C)+(TT/C)-HR/(C*C))*2;
      PHI2=-(TT-HR/C).^2+(R/C)^2+(ZS/C).^2;
      H0=PHI2/PHI1;
      if (H0>HR) 
          H0=HR;
      end
      if (H0<0) 
          H0=0;
      end

      IFAIL=0;

      NLIMIT=KEY2;
      EH1=(C3/C)*D01AHF(H0,HR,ER0,NPTS,RELLER,EINDRP,NLIMIT,IFAIL)*3.*R/QEH;
      if (NPTS>=MAXLIM2) 
          MAXLIM2=NPTS;
      end

      NLIMIT=KEY2;
      EH4=(C3/C)*D01AHF(HR,H1,ER0,NPTS,RELLER,EINDRP,NLIMIT,IFAIL)*3.*R/QEH;
      if (NPTS>=MAXLIM2) 
          MAXLIM2=NPTS;
      end

      NLIMIT=KEY2;
      EH2=(C3/(C*C))*D01AHF(H0,HR,ER0,NPTS,RELLER,ERADRP,NLIMIT,IFAIL)*R/QEH;
      if (NPTS>=MAXLIM2) 
          MAXLIM2=NPTS;
      end

      NLIMIT=KEY2;
      EH5=(C3/(C*C))*D01AHF(HR,H1,ER0,NPTS,RELLER,ERADRP,NLIMIT,IFAIL)*R/QEH;
      if (NPTS>=MAXLIM2) 
          MAXLIM2=NPTS;
      end

      NLIMIT=KEY2;
      EH3=C3*D01AHF(H0,HR,ER0,NPTS,RELLER,EELRP,NLIMIT,IFAIL)*3.*R/QEH;
      if (NPTS>=MAXLIM2) 
          MAXLIM2=NPTS;
      end

      NLIMIT=KEY2;
      EH6=C3*D01AHF(HR,H1,ER0,NPTS,RELLER,EELRP,NLIMIT,IFAIL)*3.*R/QEH;
      if (NPTS>=MAXLIM2) 
          MAXLIM2=NPTS;
      end

%============DETERMINATION OF THE UPPER INTEGRATION LIMIT
      ZS=-abs(ZS);
      PHI=1./(V*V)-1./(C*C);
      PHI1=ZS/(C*C)-(TT/V)-HR/(V*V);
      PHI2=(TT+HR/V)^2-(R/C)^2-(ZS/C)^2;
      DELTAP=abs(sqrt(PHI1*PHI1-PHI*PHI2));
      H1=(-PHI1-DELTAP)/PHI;
      if (H1<HR) 
          H1=HR;
      end
      if (H1>H) 
          H1=H;
      end

%=============DETERMINATION OF THE LOWER INTEGRATION LIMIT
      PHI1=(ZS/(C*C)+(TT/C)-HR/(C*C))*2;
      PHI2=-(TT-HR/C)^2+(R/C)^2+(ZS/C).^2;
      H0=PHI2/PHI1;
      if (H0>HR) 
          H0=HR;
      end
      if (H0<0) 
          H0=0;
      end

      NLIMIT=KEY2;
      EH7=-(C3/C)*D01AHF(H0,HR,ER0,NPTS,RELLER,EINDRP,NLIMIT,IFAIL)*3.*R/QEH;
      if (NPTS>=MAXLIM2) 
          MAXLIM2=NPTS;
      end

      NLIMIT=KEY2;
      EH10=-(C3/C)*D01AHF(HR,H1,ER0,NPTS,RELLER,EINDRP,NLIMIT,IFAIL)*3.*R/QEH;
      if (NPTS>=MAXLIM2) 
          MAXLIM2=NPTS;
      end

      NLIMIT=KEY2;
      EH8=-(C3/(C*C))*D01AHF(H0,HR,ER0,NPTS,RELLER,ERADRP,NLIMIT,IFAIL)*R/QEH;
      if (NPTS>=MAXLIM2) 
          MAXLIM2=NPTS;
      end

      NLIMIT=KEY2;
      EH11=-(C3/(C*C))*D01AHF(HR,H1,ER0,NPTS,RELLER,ERADRP,NLIMIT,IFAIL)*R/QEH;
      if (NPTS>=MAXLIM2) 
          MAXLIM2=NPTS;
      end

      NLIMIT=KEY2;
      EH9=-C3*D01AHF(H0,HR,ER0,NPTS,RELLER,EELRP,NLIMIT,IFAIL)*3.*R/QEH;
      if (NPTS>=MAXLIM2) 
          MAXLIM2=NPTS;
      end

      NLIMIT=KEY2;
      EH12=-C3*D01AHF(HR,H1,ER0,NPTS,RELLER,EELRP,NLIMIT,IFAIL)*3.*R/QEH;
      if (NPTS>=MAXLIM2) 
          MAXLIM2=NPTS;
      end

      ZS=abs(ZS);




%	============================= ELECTROSTATIC COMPONENT OF Er ==============
function y=EELRP(x)
%	===========================
%       IMPLICIT NONE
% 	    REAL*8 INTPULSE,P
%       REAL*8 MU,%,PI,EPS,V,VC,W,H,ZLAM,R,ZS,TT,HR,RR,RG,ALPHA_TOWER
%       REAL*8 ICI,ICII,x,RZ,QEV,QH,QEH
% 	    REAL*8 A1,A2,A3,B1,B2,B3
% 	    INTEGER NC,J,MODEL,MPROBLEM
%       EXTERNAL INTPULSE,P

      global  C
      global  V VC
      global  R ZS
      global  TT
      global  HR RR RG ALPHA_TOWER NC
      global  QEH

      RZ=sqrt((R*R)+((x-ZS)^2));
      ICI=0;
      ICII=0;
	if any(x<HR) 
%===========1	
	 for J=0:1:NC
          B1=TT-RZ/C-(HR-x)/C-2.*J*HR/C;
          if any(B1>0) 
              A1=0;
              ICII=INTPULSE(A1,B1);
              ICI=ICI+ICII*(RR^J)*(RG^J)*exp(-J*HR*ALPHA_TOWER);
          end
          B2=B1-2.*x/C;
          if any(B2>0) 
             A2=0;
             ICII=INTPULSE(A2,B2);
             ICI=ICI+ICII*(RR^J)*(RG^(J+1))*exp(-(J*HR+x)*ALPHA_TOWER);
          end
    end
	ICI=ICI*(1-RR);
    else
%===========x GREATER THAN HR
	B1=TT-RZ/C-(x-HR)/V;
	      if any(B1>0) 
	       A3=(x-HR)/VC-(x-HR)/V;
	       B3=TT-RZ/C-(x-HR)/VC;
	       ICII=INTPULSE(A3,B3);
		   ICI=ICI+ICII.*P(x);
		  end
%===========2
	B1=TT-RZ/C-(x-HR)/C;
	    if any(B1>0) 
	     A1=0;
		 ICII=INTPULSE(A1,B1);
		 ICI=ICI+ICII*(-RR);
		end
%===========3
	   for J=0:1:NC
           B1=TT-RZ/C-(x+HR)/C-2.*J*HR/C;
           if any(B1>0) 
                 A1=0;
                 ICII=INTPULSE(A1,B1);
                 ICI=ICI+ICII*(1-RR)*(1+RR)*(RR^J)*(RG^(J+1))*exp(-(J+1)*HR*ALPHA_TOWER);
           end
       end
%================================	
	
	end

    y=(ZS-x).*ICI./(RZ.^5).*QEH;
      

%	============================= INDUCTION COMPONENT OF Er =========
function y=EINDRP(x)
%	==============================
%       IMPLICIT NONE
% 	    REAL*8 IPULSE,P
%       REAL*8 MU,%,PI,EPS,V,VC,W,H,ZLAM,R,ZS,TT,HR,RR,RG,ALPHA_TOWER
%       REAL*8 ICI,ICII,x,B1,B2,RZ,QEV,QH,QEH
% 	    INTEGER NC,J,MODEL,MPROBLEM
%       EXTERNAL IPULSE,P

      global  C
      global  VC
      global  R ZS
      global  TT
      global  HR RR RG ALPHA_TOWER NC
      global  QEH

      RZ=sqrt((R*R)+((x-ZS)^2));
      ICI=0;
      ICII=0;
	if any(x<HR) 
%===========1	
	 for J=0:1:NC
          B1=TT-RZ/C-(HR-x)/C-2.*J*HR/C;
          if any(B1>0) 
              ICII=IPULSE(B1);
              ICI=ICI+ICII*(RR^J)*(RG^J)*exp(-J*HR*ALPHA_TOWER);
          end
          B2=B1-2.*x/C;
          if any(B2>0) 
              ICII=IPULSE(B2);
              ICI=ICI+ICII*(RR^J)*(RG^(J+1))*exp(-(J*HR+x)*ALPHA_TOWER);
          end
     end
	ICI=ICI*(1-RR);
    else
%===========x GREATER THAN HR
	    B1=TT-RZ/C-(x-HR)/VC;
	    if any(B1>0) 
            ICII=IPULSE(B1);
            ICI=ICI+ICII.*P(x);
		end
%===========2
	    B1=TT-RZ/C-(x-HR)/C;
	    if any(B1>0) 
             ICII=IPULSE(B1);
             ICI=ICI+ICII*(-RR);
		end
%===========3
	  for J=0:1:NC
          B1=TT-RZ/C-(x+HR)/C-2.*J*HR/C;
          if any(B1>0) 
             ICII=IPULSE(B1);
             ICI=ICI+ICII*(1-RR)*(1+RR)*(RR^J)*(RG^(J+1))*exp(-(J+1)*HR*ALPHA_TOWER);
          end
      end
%================================	
	
	end

    y=(ZS-x).*ICI./(RZ.^4)*QEH;


%	============================= RADIATION COMPONENT OF Er ========
function y=ERADRP(x)
%	=============================
%       IMPLICIT NONE
% 	  REAL*8 DIPULSE,P
%       REAL*8 MU,%,PI,EPS,V,VC,W,H,ZLAM,R,ZS,TT,HR,RR,RG,ALPHA_TOWER
%       REAL*8 ICI,ICII,x,B1,B2,RZ,QEV,QH,QEH
% 	  INTEGER NC,J,MODEL,MPROBLEM
%       EXTERNAL DIPULSE,P

      global  C
      global  VC
      global  R ZS
      global  TT
      global  HR RR RG ALPHA_TOWER NC
      global  QEH

      RZ=sqrt((R*R)+((x-ZS)^2));
      ICI=0;
      ICII=0;
	if any(x<HR) 
%===========1	
	for J=0:1:NC
          B1=TT-RZ/C-(HR-x)/C-2.*J*HR/C;
          if any(B1>0) 
              ICII=DIPULSE(B1);
              ICI=ICI+ICII*(RR^J)*(RG^J)*exp(-J*HR*ALPHA_TOWER);
          end
          B2=B1-2.*x/C;
          if any(B2>0) 
              ICII=DIPULSE(B2);
              ICI=ICI+ICII*(RR^J)*(RG^(J+1))*exp(-(J*HR+x)*ALPHA_TOWER);
          end
    end
	ICI=ICI*(1-RR);
    else
%===========x GREATER THAN HR
	    B1=TT-RZ/C-(x-HR)/VC;
	    if any(B1>0) 
            ICII=DIPULSE(B1);
            ICI=ICI+ICII.*P(x);
		end
%===========2
	    B1=TT-RZ/C-(x-HR)/C;
	    if any(B1>0) 
             ICII=DIPULSE(B1);
             ICI=ICI+ICII*(-RR);
		end
%===========3
	  for J=0:1:NC
          B1=TT-RZ/C-(x+HR)/C-2.*J*HR/C;
          if any(B1>0) 
             ICII=DIPULSE(B1);
             ICI=ICI+ICII*(1-RR)*(1+RR)*(RR^J)*(RG^(J+1))*exp(-(J+1)*HR*ALPHA_TOWER);
          end
      end
%================================	
	
	end

      y=(ZS-x).*ICI./(RZ.^3)*QEH;
    