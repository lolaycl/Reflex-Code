      PROGRAM REFLEX_15

c	Added variables HPtower and EELtower


      IMPLICIT NONE

	
	REAL*8 MU,C,PI,EPS,V,VC,W,H,ZLAM,R,ZS,R1,FV,IFRONT,Z1max,RZ             
      REAL*8 IW,ETA,TAU1,TAU2,IW2,ETA2,TAU12,TAU22,I0B,T0B
      REAL*8 ALPHA,BETA,I01,GAMMA,DELTA,I02      
      REAL*8 HR,RR,RG,ALPHA_TOWER,ERR,ERR2                           
      REAL*8 DELTATE,DELTATE1,DELTATE2
      REAL*8 TE,DTS,TMAX,D1,D2,TINT,TT,T01,A,B,CC,RADIC
      REAL*8 EV1,EV2,EV3,EV4,EV5,EV6,EV7,EV8,EV9,EV10,EV11,EV12,QEV
      REAL*8 EH1,EH2,EH3,EH4,EH5,EH6,EH7,EH8,EH9,EH10,EH11,EH12,QEH
      REAL*8 HP1,HP2,HP3,HP4,HP5,HP6,HP7,HP8,QH
	REAL*8 EV1TO,EV2TO,EH1TO,EH2TO,HP1TO,HP2TO
      REAL*8 SIG,EPSR,MK1,MK,FF
	REAL*8 F1,F2,F3,F4,TM,HAB,HAM,HC
      REAL*8 IPULSE,DIPULSE,P
      INTEGER NBE,NBE1,NBE2,NW,NW2,NC,N,NSIMP,MAXLIM1,KEY1,MAXLIM2,KEY2
      INTEGER MX1,MX2,K,S,MODEL,MPROBLEM
      REAL*8 DEVO,DHPHI,EHF
      REAL*8 ,ALLOCATABLE::EVO(:),EEL(:),EIND(:),ERAD(:)
c     The alocatable variables reserve memory space for the vectors
      REAL*8 ,ALLOCATABLE::EHO(:),EHEL(:),EHIND(:),EHRAD(:)    
      REAL*8 ,ALLOCATABLE::HPHI(:),HIND(:),HRAD(:),EVtower(:),HPtower(:)
	REAL*8 ,ALLOCATABLE::EV_TO(:),EH_TO(:),HP_TO(:),EELtower(:)
      REAL*8 ,ALLOCATABLE::ECO(:)
	REAL*8 ,ALLOCATABLE::CURT(:),DCURT(:),CURC(:),DCURC(:),CURB(:)
	1                    ,DCURB(:),CURR1(:),DCURR1(:)
      REAL*8 ,ALLOCATABLE::TIMEE(:)
      
      CHARACTER*1 CHO,RE
      CHARACTER*13 MDCUR,MDEV,MDHP
      CHARACTER*5 MCUR,MEVH,MHP 
      CHARACTER*15 NOME
c     L'instruction COMMON permet le regroupement de zones mémoires pouvant 
c     être partagées par différentes unités de programme (fonctions, procédures).
	COMMON /COM1/ MU,C,PI,EPS,IFRONT
      COMMON /COM2/ V,VC,W,H,ZLAM,HC
      COMMON /COM3/ R,ZS,R1,FV
      COMMON /COM4/ IW,ETA,TAU1,TAU2,IW2,ETA2,TAU12,TAU22,NW2,NW
      COMMON /COM5/ ALPHA,BETA,I01,GAMMA,DELTA,I02
      COMMON /COM6/ HR,RR,RG,ALPHA_TOWER,NC,MODEL,MPROBLEM
      COMMON /COM7/ ERR,ERR2
      COMMON /COM8/ NSIMP
	COMMON /COM9/ RZ
	COMMON /COMT/ TT
	COMMON /COMBARB/ I0B,T0B
      COMMON /EV/ EV1,EV2,EV3,EV4,EV5,EV6,EV7,EV8,EV9,EV10,EV11,EV12,
	1			EV1TO,EV2TO
      COMMON /EH/ EH1,EH2,EH3,EH4,EH5,EH6,EH7,EH8,EH9,EH10,EH11,EH12,
	1			EH1TO,EH2TO
      COMMON /HP/ HP1,HP2,HP3,HP4,HP5,HP6,HP7,HP8,HP1TO,HP2TO
      COMMON /SCALE/ QEV,QH,QEH
      COMMON /INTEG/ MAXLIM1,KEY1,MAXLIM2,KEY2
      
	EXTERNAL IPULSE,FF,DIPULSE,P
      DATA C,PI,EPS /3.D8,3.141592654D0,8.8542D-12/
	     
      
      MAXLIM1=0
      
      MAXLIM2=0
C     =========================================      
C     READING INPUT FILE NAME
C     =========================================
      MU=4.*PI*1.D-7
      RE=CHAR(13)
C      PRINT*,'INPUT FILE NAME?    ( MAX 15)'
C      READ*,NOME
	NOME = 'INFILE.DAT'
	OPEN (UNIT=9,FILE=NOME,STATUS='OLD',IOSTAT=S)
	IF(S.GT.0)THEN
C     .GT.   Greater than
      PRINT*,'ERROR - PRESS RETURN!'
	READ(*,'(A1)')CHO   
	STOP
	ENDIF
C17    PRINT*,'MODEL?'
C      PRINT*,'1 - TL'
C      PRINT*,'2 - MTLL'
C      PRINT*,'3 - MTLE'
C	PRINT*,'4 - BG'
C	PRINT*,'5 - TCS'
      READ(9,*) MODEL,MPROBLEM
	IF (MODEL.GT.5) THEN
      PRINT*,'Model not valid'
      READ(*,'(A1)')CHO   
	STOP
	END IF


	IF (MODEL.LT.1) THEN
	PRINT*,'Model not valid'
      READ(*,'(A1)')CHO   
	STOP
      END IF


C     ============================     
C     READING INPUT DATA
C     ============================
C     I01,I02,ALPHA,BETA,GAMMA,DELTA: parameters of bi-EXP
C      function
C     IW,NW,TAU1,TAU2: parameters of Heidler function
C     V: return-stroke speed
C     W: equal to speed of light (never used actually)
C     H: total height of the lightning channel
C     ZLAM: current attenuation factor with height for MTLE model
C     HR: height of the tower
C     RR: current reflection coefficient at the tower top
C     RG: current reflection coefficient at the tower base
C	HAB,HAM: heights at which the current is observed
C	ERR,ERR2: parameters for the precision of the numerical integration 
C     D1: distance of the observation point (projection on X axis)
C	D1: distance of the observation point (projection on Y axis)
C     ZS: height of the observation point over ground
C	SIG: conductivity of the soil
C	EPSR: dielectric constant of the soil
C	DELTATE1,TINT,DTS,TMAX: parameters for time resolution
C	NSIMP: parameter for precision of the integration routine of Simpson
C    
    	READ(9,*)I01,ALPHA,BETA
      READ(9,*)I02,GAMMA,DELTA
C      READ(9,'(A1)')CHO
	READ(9,*)IW,TAU1,TAU2,NW
      READ(9,*)IW2,TAU12,TAU22,NW2
      READ(9,*)I0B,T0B
C      READ(9,'(A1)')CHO
	READ(9,*)V,W,H,ZLAM
C      READ(9,'(A1)')CHO
	READ(9,*)HR,RR,RG,ALPHA_TOWER
C	READ(9,'(A1)')CHO
	READ(9,*)HAB,HAM
C	READ(9,'(A1)')CHO
	READ(9,*)ERR,ERR2
C      READ(9,'(A1)')CHO
	ETA=EXP(-(TAU1/TAU2)*((NW*TAU2/TAU1)**(1.D0/FLOAT(NW))))
      ETA2=EXP(-(TAU12/TAU22)*((NW2*TAU22/TAU12)**
     1     (1.D0/FLOAT(NW2))))
      
      READ(9,*)ZS,SIG,EPSR
C      READ(9,'(A1)')CHO
	READ(9,*)D1,D2,NSIMP
      READ(9,*)DELTATE1,TINT,DTS,TMAX
C      READ(9,'(A1)')CHO
	READ(9,*)QEV,QH,QEH
      READ(9,*)KEY1,KEY2
      READ(9,*)F1
      READ(9,*)F2
      READ(9,*)F3
      READ(9,*)F4
      CLOSE(9)                   

      IF (F2.EQ.1D-3) THEN  
                      MCUR='   KA'
                      ELSE
                      MCUR='    A'
      ENDIF
      IF (F4.EQ.1D-3) THEN
                      MHP='KA/(M'
                      ELSE
                      MHP=' A/(M'
      ENDIF
      IF (F1.EQ.1D-3) THEN
                      MEVH='KV/(M'
                      ELSE
                      MEVH=' V/(M'
      END IF
      IF (F3.EQ.1.D6) THEN 
                      MDCUR=MCUR//'/MICROS'
                      MDHP=MHP//'*MICROS)'
                      MDEV=MEVH//'*MICROS)'
                      ELSE
                      MDCUR=MCUR//'/S'
                      MDHP=MHP//'*S)'
                      MDEV=MEVH//'*S)'
      ENDIF

C     ==============================================
C     SELECTION OF CURRENT-WAVE SPEED
C     ==============================================
C	VC: current-wave speed
C	C: speed of light

      IF ((MODEL.EQ.1).OR.(MODEL.EQ.2).OR.(MODEL.EQ.3)) VC=V
	IF (MODEL.EQ.4) VC=1.D15
	IF (MODEL.EQ.5) VC=-C 

C     ====================================      
C     INPUT DATA DISPLAY 
C     ====================================      
C	WRITE(*,30)'DISPLAY OF INPUT DATA'
30    FORMAT('*',/,26X,A26,/,26X,26('-'))
C      WRITE(*,39)'I01','ALPHA','BETA','I02','GAMMA','DELTA'
39    FORMAT(1X,6(5X,A5,3X))
C      WRITE(*,40)I01,ALPHA,BETA,I02,GAMMA,DELTA
40    FORMAT(1X,6(E10.3,3X))      
C      WRITE(*,34)'IW','TAU1','TAU2','NW','ETA'
34    FORMAT(1X,5(5X,A5,3X))      
C      WRITE(*,33)IW,TAU1,TAU2,NW,ETA  
33    FORMAT(1X,3(E10.3,3X),I10,3X,E10.3)      
C      WRITE(*,34)'IW2','TAU12','TAU22','NW2','ETA2'
C      WRITE(*,33)IW2,TAU12,TAU22,NW2,ETA2 
C      WRITE(*,39)'ZLAM','V_CAN','H_CAN','HR_OG','SIGMA','H_OSM'
35    FORMAT(1X,5(5X,A5,3X))
C      WRITE(*,40)ZLAM,V,H,HR,SIG,HAM
36    FORMAT(1X,4(E10.3,3X))
      IF (ZLAM.NE.0) ZLAM=1/ZLAM
C	WRITE(*,39)'DELT1','TINT','DTS','TMAX','EPSR','H_OSB'
C      WRITE(*,40)DELTATE1,TINT,DTS,TMAX,EPSR,HAB
29    FORMAT(1X,5(E10.3,3X))

C      WRITE(*,31)'CR_RR','CR_RG','VEL_W','ERR','ERR2'
31    FORMAT(1X,5(5X,A5,3X))
C      WRITE(*,38)RR,RG,W,ERR,ERR2
38    FORMAT(1X,5(E10.3,3X))
C      WRITE(*,39)'H_LIN','D1','D2','NSIMP','NLIM1','NLIM2'
C      WRITE(*,32)ZS,D1,D2,NSIMP,KEY1,KEY2
32    FORMAT(1X,3(E10.3,3X),3(I10,3X))

      R=DSQRT(D1**2+D2**2)
      DELTATE2=DELTATE1*DTS
      NBE1=INT(TINT/DELTATE1)
      NBE2=INT((TMAX-TINT)/DELTATE2)
      NBE=NBE1+NBE2

C      WRITE(*,34)'QEV','QH','QEH','NBE1','NBE2'
C      WRITE(*,37)QEV,QH,QEH,NBE1,NBE2
37    FORMAT(1X,3(E10.1,3X),2(I10,3X))
C      WRITE(*,'(4X,A9,1X,A5,8X,2(2X,A9,1X,A5,A1,7X))')'CUR IN',
C     1  MCUR,'EVH IN',MEVH,')','HPH IN',MHP,')'
C      WRITE(*,'(2X,3(2X,A9,1X,A13))')'DCUR IN',MDCUR,'DEV/DT IN',MDEV,
C     1   'DHP/DT IN',MDHP
C      PRINT*,'CORRECT?    (Y/N)'
C      READ(*,'(A1)')CHO
C      IF ((CHO.EQ.'N').OR.(CHO.EQ.'n')) STOP


C     ==========================================   
C     DINAMIC ALLOCATION OUTPUT ARRAYS
C     ==========================================    
	ALLOCATE(EVO(0:NBE1+NBE2))
	ALLOCATE(EEL(0:NBE1+NBE2))
	ALLOCATE(EIND(0:NBE1+NBE2))
	ALLOCATE(ERAD(0:NBE1+NBE2))
	ALLOCATE(EV_TO(0:NBE1+NBE2))

	
	ALLOCATE(EHO(0:NBE1+NBE2))
	ALLOCATE(EHEL(0:NBE1+NBE2))
	ALLOCATE(EHIND(0:NBE1+NBE2))
	ALLOCATE(EHRAD(0:NBE1+NBE2))
	ALLOCATE(EH_TO(0:NBE1+NBE2))


	ALLOCATE(EELtower(0:NBE1+NBE2))
	ALLOCATE(EVtower(0:NBE1+NBE2))
	ALLOCATE(HPtower(0:NBE1+NBE2))
	ALLOCATE(HPHI(0:NBE1+NBE2))
	ALLOCATE(HIND(0:NBE1+NBE2))
	ALLOCATE(HRAD(0:NBE1+NBE2))
	ALLOCATE(HP_TO(0:NBE1+NBE2))

	ALLOCATE(ECO(0:NBE1+NBE2))
	
	ALLOCATE(CURT(0:NBE1+NBE2))
	ALLOCATE(DCURT(0:NBE1+NBE2))
	ALLOCATE(CURC(0:NBE1+NBE2))
	ALLOCATE(DCURC(0:NBE1+NBE2))
	ALLOCATE(CURB(0:NBE1+NBE2))
	ALLOCATE(DCURB(0:NBE1+NBE2))
	ALLOCATE(CURR1(0:NBE1+NBE2))
	ALLOCATE(DCURR1(0:NBE1+NBE2))


	ALLOCATE(TIMEE(0:NBE1+NBE2))


C     =============================================      
C     FIELDS CALCULATION 
C     =============================================
C	TE: time step


      TE=0.0
C      WRITE(*,41)NBE,'MAXL1','MAXL2','MX1','MX2'
41    FORMAT(' NMAX=',I3,4(2X,A6),/)
     
      DO 5 N=0,NBE
      IF (MAXLIM1.GE.MX1) MX1=MAXLIM1
      IF (MAXLIM2.GE.MX2) MX2=MAXLIM2
      
      IF (MPROBLEM .EQ. 1) THEN
      WRITE(*,42)N,NC
	ENDIF
42    FORMAT('N=',I5, '            NC=',I3,/)  
	
      MAXLIM1=0
      MAXLIM2=0
      IF (N.LE.NBE1) THEN
      DELTATE=DELTATE1
      ELSE
      DELTATE=DELTATE2
      END IF
      TIMEE(N)=TE

      IF (HR.GT.0.) THEN
       NC=TE/(2*HR/C)
      ELSE
       NC=0
      END IF


C     actual wavefront current

	RZ=0
	HC=HR+V*TE

      CALL CURRENT(TE,HC,CURC(N),DCURC(N))
	IFRONT=CURC(N)

C	Turn-On Term fields

      ZS=DABS(ZS)
      T01=DSQRT(R*R+(HR-ZS)*(HR-ZS))/C
      TT=TE+T01
	A=1/(V*V)-1/(C*C)
	B=2*ZS/(C*C)-2*HR/(V*V)-2*TT/V
	CC=TT*TT+2*TT*HR/V+HR*HR/(V*V)-R*R/(C*C)
	RADIC=DABS(B*B-4*A*CC)
	Z1max=(-B-DSQRT(RADIC))/(2*A)
	R1=DSQRT(R*R+(ZS-Z1max)**2)
	RZ=R1
	CALL CURRENT(TT,Z1max,CURR1(N),DCURR1(N))
	IFRONT=CURR1(N)
	FV=1/(1/V-(ZS-Z1max)/C/R1)
	EV1TO=-IFRONT*R*R*FV/(4*PI*EPS*C*C*R1**3)
	HP1TO=IFRONT*R*FV/(4*PI*C*R1**2)
	EH1TO=IFRONT*R*(ZS-Z1max)*FV/(4*PI*EPS*C*C*R1**3)


      ZS=-DABS(ZS)
      TT=TE+T01
	A=1/(V*V)-1/(C*C)
	B=2*ZS/(C*C)-2*HR/(V*V)-2*TT/V
	CC=TT*TT+2*TT*HR/V+HR*HR/(V*V)-R*R/(C*C)
	RADIC=DABS(B*B-4*A*CC)
	Z1max=(-B-DSQRT(RADIC))/(2*A)
	R1=DSQRT(R*R+(ZS-Z1max)**2)
	RZ=R1
	CALL CURRENT(TT,Z1max,CURR1(N),DCURR1(N))
	IFRONT=CURR1(N)
	FV=1/(1/V-(ZS-Z1max)/C/R1)
	EV2TO=-IFRONT*R*R*FV/(4*PI*EPS*C*C*R1**3)
	HP2TO=IFRONT*R*FV/(4*PI*C*R1**2)
	EH2TO=-IFRONT*R*(ZS-Z1max)*FV/(4*PI*EPS*C*C*R1**3)


      CALL HZONTAL(TE)
      EHIND(N)=EH1+EH4+EH7+EH10
      EHRAD(N)=EH2+EH5+EH8+EH11
	EH_TO(N)=EH1TO+EH2TO
      EHEL(N)=EH3+EH6+EH9+EH12
      EHO(N)=EHIND(N)+EHRAD(N)+EHEL(N)
C     =============================================================
C	OUTPUT VARIABLES OF ROUTINE HZONTAL
C     EH1  INDUCTION  COMPONENT OF Er DUE TO THE TOWER FOR ZS>0 
C     EH2  RADIAL	    COMPONENT OF Er DUE TO THE TOWER FOR ZS>0
C     EH3  ELECTROST. COMPONENT OF Er DUE TO THE TOWER FOR ZS>0   
   
C     EH4  INDUCTION  COMPONENT OF Er DUE TO THE CHANNEL FOR ZS>0 
C     EH5  RADIAL	    COMPONENT OF Er DUE TO THE CHANNEL FOR ZS>0
C     EH6  ELECTROST. COMPONENT OF Er DUE TO THE CHANNEL FOR ZS>0      

C     EH7  INDUCTION  COMPONENT OF Er DUE TO THE TOWER FOR ZS<0 
C     EH8  RADIAL	    COMPONENT OF Er DUE TO THE TOWER FOR ZS<0
C     EH9  ELECTROST. COMPONENT OF Er DUE TO THE TOWER FOR ZS<0   

C     EH10  INDUCTION  COMPONENT OF Er DUE TO THE CHANNEL FOR ZS<0 
C     EH11  RADIAL	 COMPONENT OF Er DUE TO THE CHANNEL FOR ZS<0
C     EH12  ELECTROST. COMPONENT OF Er DUE TO THE CHANNEL FOR ZS<0      
C     =============================================================

      CALL VERTICAL(TE)
      EIND(N)=EV1+EV4+EV7+EV10
      ERAD(N)=EV2+EV5+EV8+EV11
	EV_TO(N)=EV1TO+EV2TO
	EELtower(N)=EV3+EV9
      EEL(N)=EV3+EV6+EV9+EV12
	EVtower(N)=EV1+EV2+EV3+EV7+EV8+EV9
      EVO(N)=EIND(N)+ERAD(N)+EEL(N)

      HIND(N)=(HP1+HP3+HP5+HP7)*(-1.D0) 
      HRAD(N)=(HP2+HP4+HP6+HP8)*(-1.D0)
	HP_TO(N)=HP1TO+HP2TO
	HPtower(N)=(HP1+HP2+HP5+HP6)*(-1.D0)
	HPHI(N)=HIND(N)+HRAD(N)

C	=============================================================
C	OUTPUT VARIABLES OF ROUTINE VERTICAL
C
C     EV1  INDUCTION  COMPONENT OF Ez DUE TO THE TOWER FOR ZS>0 
C     EV2  RADIAL	    COMPONENT OF Ez DUE TO THE TOWER FOR ZS>0
C     EV3  ELECTROST. COMPONENT OF Ez DUE TO THE TOWER FOR ZS>0   

C     EV4  INDUCTION  COMPONENT OF Ez DUE TO THE CHANNEL FOR ZS>0 
C     EV5  RADIAL	    COMPONENT OF Ez DUE TO THE CHANNEL FOR ZS>0
C     EV6  ELECTROST. COMPONENT OF Ez DUE TO THE CHANNEL FOR ZS>0      

C     HP1  INDUCTION  COMPONENT OF H DUE TO THE TOWER FOR ZS>0  
C     HP2  RADIAL	    COMPONENT OF H DUE TO THE TOWER FOR ZS>0  
C     HP3  INDUCTION  COMPONENT OF H DUE TO THE CHANNEL FOR ZS>0  
C     HP4  RADIAL	    COMPONENT OF H DUE TO THE CHANNEL FOR ZS>0  
C------------------------------------------------------------------
C     EV7  INDUCTION  COMPONENT OF Ez DUE TO THE TOWER FOR ZS<0 
C     EV8  RADIAL	    COMPONENT OF Ez DUE TO THE TOWER FOR ZS<0
C     EV9  ELECTROST. COMPONENT OF Ez DUE TO THE TOWER FOR ZS<0   

C     EV10  INDUCTION  COMPONENT OF Ez DUE TO THE CHANNEL FOR ZS<0 
C     EV11  RADIAL	 COMPONENT OF Ez DUE TO THE CHANNEL FOR ZS<0
C     EV12  ELECTROST. COMPONENT OF Ez DUE TO THE CHANNEL FOR ZS<0      

C     HP5  INDUCTION  COMPONENT OF H DUE TO THE TOWER FOR ZS<0  
C     HP6  RADIAL	    COMPONENT OF H DUE TO THE TOWER FOR ZS<0  
C     HP7  INDUCTION  COMPONENT OF H DUE TO THE CHANNEL FOR ZS<0  
C     HP8  RADIAL	    COMPONENT OF H DUE TO THE CHANNEL FOR ZS<0  
C==================================================================      
      
	RZ=0
	CALL CURRENT(TE,HAM,CURT(N),DCURT(N))
	CALL CURRENT(TE,HAB,CURB(N),DCURB(N))
	
		
	TE=TE+DELTATE


5     CONTINUE


C     ==================================================================     
C     CORRECTIVE TERM FOR HORIZONTAL FIELD (APPROACH Rubinstein)
C     ==================================================================
	ECO = 0.0
	ECO(0)=0.0
	ECO(1)=HPHI(1)/TIMEE(1)*FF(TIMEE(1),SIG,EPSR)
      DO 88 N=2,NBE
      DO 98 K=1,N-1
	MK=(HPHI(K+1)-HPHI(K))/(TIMEE(K+1)-TIMEE(K))
      MK1=(HPHI(K)-HPHI(K-1))/(TIMEE(K)-TIMEE(K-1))
      ECO(N)=ECO(N)+(MK-MK1)*FF(TIMEE(N)-TIMEE(K),SIG,EPSR)
98    CONTINUE
      ECO(N)=ECO(N)+HPHI(1)/TIMEE(1)*FF(TIMEE(N),SIG,EPSR)
88    CONTINUE
C     ===================================================================      
C     WRITING IN THE OUTPUT FILE : RESULT.DAT
C     ===================================================================
      OPEN (11, FILE ='result.dat')
      WRITE(11,500)'T','(-)EV','EVEL','EVIN','EVRD','EV_TO','EHO',
     1'EHF','EHEL','EHIN','EHRD','EH_TO','HPHI','HPIND','HPRAD',
     2'HP_TO','(-)DEV','DHPHI','I0(T)','IT(T)','IM(T)',
     3'IB(T)','IFRONT(T)','DI0(T)','DIT(T)','DIM(T)','DIB(T)',
     4'EVtower(T)','HPtower(T)','EELtower(T)'

500   FORMAT(30(A11,1X))
      
      DO 400 N=0,NBE             
      IF (N.EQ.0) THEN
      DEVO=0.0
      DHPHI=0.0
      ELSE
      DEVO=(EVO(N)-EVO(N-1))/(TIMEE(N)-TIMEE(N-1))
      DHPHI=(HPHI(N)-HPHI(N-1))/(TIMEE(N)-TIMEE(N-1))
      EHF=EHO(N)+ECO(N)
	END IF
      
	TM=TIMEE(N)

	WRITE(11,300)TM*F3,-EVO(N)*F1,-EEL(N)*F1,-EIND(N)*F1,-ERAD(N)*F1,
     1-EV_TO(N)*F1,EHO(N)*F1,EHF*F1,EHEL(N)*F1,EHIND(N)*F1,
     2EHRAD(N)*F1,EH_TO(N)*F1,-HPHI(N)*F4,-HIND(N)*F4,-HRAD(N)*F4,
     3HP_TO(N)*F4,-DEVO*F1/F3,DHPHI*F4/F3,
     4IPULSE(TM)*F2,CURT(N)*F2, CURC(N)*F2, CURB(N)*F2, CURR1(N)*F2,
     5DIPULSE(TM)*F2/F3,DCURT(N)*F2/F3,DCURC(N)*F2/F3,DCURB(N)*F2/F3,
     6-EVtower(N)*F1,-HPtower(N)*F4,-EELtower(N)*F1

300   FORMAT(30(1E11.4,1X))
400   CONTINUE      
      CLOSE(11)
      STOP
      END
C     ==========================================================  
C     CALCULATION OF VERTICAL ELECTRIC FIELD AND MAGNETIC FIELD
C     ==========================================================
      SUBROUTINE VERTICAL(TE)
	
	IMPLICIT NONE     
      REAL*8 MU,C,PI,EPS,V,VC,W,H,ZLAM,R,ZS,R1,FV,IFRONT
      REAL*8 TT,D01AHF
      REAL*8 HR,RR,RG,ALPHA_TOWER,ERR,ERR2
      REAL*8 EV1,EV2,EV3,EV4,EV5,EV6,EV7,EV8,EV9,EV10,EV11,EV12
      REAL*8 QEV,QH,QEH
      REAL*8 HP1,HP2,HP3,HP4,HP5,HP6,HP7,HP8
	REAL*8 EV1TO,EV2TO,HP1TO,HP2TO
      REAL*8 PHI,PHI1,PHI2,DELTAP
      REAL*8 TE,T01,H0,H1
      REAL*8 C1,C3,RELLER,ER0
      REAL*8 EELP,EINDP,ERADP,BINDP,BRADP
      INTEGER NC,IFAIL,NPTS,NLIMIT,MAXLIM1,KEY1,MAXLIM2,KEY2,MODEL
	INTEGER MPROBLEM

      COMMON /COM1/ MU,C,PI,EPS,IFRONT
      COMMON /COM2/ V,VC,W,H,ZLAM
      COMMON /COM3/ R,ZS,R1,FV
      COMMON /COMT/ TT
      COMMON /COM6/ HR,RR,RG,ALPHA_TOWER,NC,MODEL,MPROBLEM
      COMMON /COM7/ ERR,ERR2
      COMMON /EV/ EV1,EV2,EV3,EV4,EV5,EV6,EV7,EV8,EV9,EV10,EV11,EV12,
	1			EV1TO,EV2TO
      COMMON /HP/ HP1,HP2,HP3,HP4,HP5,HP6,HP7,HP8,HP1TO,HP2TO
      COMMON /SCALE/ QEV,QH,QEH
      COMMON /INTEG/ MAXLIM1,KEY1,MAXLIM2,KEY2
      
      EXTERNAL EELP,EINDP,ERADP,BINDP,BRADP
      EXTERNAL D01AHF
      
      ZS=DABS(ZS)
      T01=DSQRT(R*R+(HR-ZS)*(HR-ZS))/C

C===========T01: OBSERVATION TIME DELAY
      TT=TE+T01
C===========TT: ACTUAL TIME

	IF (TT.EQ.T01) GOTO667
      C1=R/(4.*PI)
      C3=1./(4.*PI*EPS)

C===========DETERMINATION OF THE UPPER INTEGRATION LIMIT
      PHI=1./(V*V)-1./(C*C)
      PHI1=ZS/(C*C)-(TT/V)-HR/(V*V)
      PHI2=(TT+HR/V)**2.-(R/C)**2.-(ZS/C)**2.
      DELTAP=DABS(DSQRT(PHI1*PHI1-PHI*PHI2))
      H1=(-PHI1-DELTAP)/PHI
      IF (H1.LT.HR) H1=HR
      IF (H1.GT.H) H1=H
C===========DETERMINATION OF THE LOWER INTEGRATION LIMIT
      PHI1=(ZS/(C*C)+(TT/C)-HR/(C*C))*2.
      PHI2=-(TT-HR/C)**2.+(R/C)**2.+(ZS/C)**2.
      H0=PHI2/PHI1
      IF (H0.GT.HR) H0=HR
      IF (H0.LT.0.) H0=0.

      ER0=ERR
      IFAIL=0
      NLIMIT=KEY1
      EV1=(C3/C)*D01AHF(H0,HR,ER0,NPTS,RELLER,EINDP,NLIMIT,IFAIL)/QEV
      IF (NPTS.GE.MAXLIM1) MAXLIM1=NPTS
            
      NLIMIT=KEY1
      EV2=-(C3*R*R/(C*C))*D01AHF(H0,HR,ER0,NPTS,RELLER,ERADP,NLIMIT
     1      ,IFAIL)/QEV
      IF (NPTS.GE.MAXLIM1) MAXLIM1=NPTS

      NLIMIT=KEY1
      EV3=C3*D01AHF(H0,HR,ER0,NPTS,RELLER,EELP,NLIMIT,IFAIL)/QEV
      IF (NPTS.GE.MAXLIM1) MAXLIM1=NPTS

      NLIMIT=KEY1
      EV4=(C3/C)*D01AHF(HR,H1,ER0,NPTS,RELLER,EINDP,NLIMIT,IFAIL)/QEV
      IF (NPTS.GE.MAXLIM1) MAXLIM1=NPTS

      NLIMIT=KEY1
      EV5=-(C3*R*R/(C*C))*D01AHF(HR,H1,ER0,NPTS,RELLER,ERADP,NLIMIT
     1       ,IFAIL)/QEV
      IF (NPTS.GE.MAXLIM1) MAXLIM1=NPTS

      NLIMIT=KEY1
      EV6=C3*D01AHF(HR,H1,ER0,NPTS,RELLER,EELP,NLIMIT,IFAIL)/QEV
      IF (NPTS.GE.MAXLIM1) MAXLIM1=NPTS

      NLIMIT=KEY1
      HP1=(C1)*D01AHF(H0,HR,ER0,NPTS,RELLER,BINDP,NLIMIT,IFAIL)/QH
      IF (NPTS.GE.MAXLIM1) MAXLIM1=NPTS

      NLIMIT=KEY1
      HP2=(C1/C)*D01AHF(H0,HR,ER0,NPTS,RELLER,BRADP,NLIMIT,IFAIL)/QH
      IF (NPTS.GE.MAXLIM1) MAXLIM1=NPTS

      NLIMIT=KEY1
      HP3=(C1)*D01AHF(HR,H1,ER0,NPTS,RELLER,BINDP,NLIMIT,IFAIL)/QH
      IF (NPTS.GE.MAXLIM1) MAXLIM1=NPTS

      NLIMIT=KEY1
      HP4=(C1/C)*D01AHF(HR,H1,ER0,NPTS,RELLER,BRADP,NLIMIT,IFAIL)/QH
      IF (NPTS.GE.MAXLIM1) MAXLIM1=NPTS


C==========DETERMINATION OF THE UPPER INTEGRATION LIMIT
      ZS=-DABS(ZS)
      PHI=1./(V*V)-1./(C*C)
      PHI1=ZS/(C*C)-(TT/V)-HR/(V*V)
      PHI2=(TT+HR/V)**2.-(R/C)**2.-(ZS/C)**2.
      DELTAP=DABS(DSQRT(PHI1*PHI1-PHI*PHI2))
      H1=(-PHI1-DELTAP)/PHI
      IF (H1.LT.HR) H1=HR      
      IF (H1.GT.H) H1=H

C==========DETERMINATION OF THE LOWER INTEGRATION LIMIT
      PHI1=(ZS/(C*C)+(TT/C)-HR/(C*C))*2.
      PHI2=-(TT-HR/C)**2.+(R/C)**2.+(ZS/C)**2.
      H0=PHI2/PHI1
      IF (H0.GT.HR) H0=HR      
      IF (H0.LT.0.) H0=0.



      NLIMIT=KEY1
      EV7=(C3/C)*D01AHF(H0,HR,ER0,NPTS,RELLER,EINDP,NLIMIT,IFAIL)/QEV
      IF (NPTS.GE.MAXLIM1) MAXLIM1=NPTS    

      NLIMIT=KEY1
      EV8=-(C3*R*R/(C*C))*D01AHF(H0,HR,ER0,NPTS,RELLER,ERADP
     1     ,NLIMIT,IFAIL)/QEV
      IF (NPTS.GE.MAXLIM1) MAXLIM1=NPTS

      NLIMIT=KEY1
      EV9=C3*D01AHF(H0,HR,ER0,NPTS,RELLER,EELP,NLIMIT,IFAIL)/QEV
      IF (NPTS.GE.MAXLIM1) MAXLIM1=NPTS

      NLIMIT=KEY1
      EV10=(C3/C)*D01AHF(HR,H1,ER0,NPTS,RELLER,EINDP,NLIMIT  
     1             ,IFAIL)/QEV
      IF (NPTS.GE.MAXLIM1) MAXLIM1=NPTS

      NLIMIT=KEY1
      EV11=-(C3*R*R/(C*C))*D01AHF(HR,H1,ER0,NPTS,RELLER,ERADP
     1              ,NLIMIT,IFAIL)/QEV
      IF (NPTS.GE.MAXLIM1) MAXLIM1=NPTS

      NLIMIT=KEY1
      EV12=C3*D01AHF(HR,H1,ER0,NPTS,RELLER,EELP,NLIMIT,IFAIL)/QEV
      IF (NPTS.GE.MAXLIM1) MAXLIM1=NPTS


      NLIMIT=KEY1
      HP5=(C1)*D01AHF(H0,HR,ER0,NPTS,RELLER,BINDP,NLIMIT,IFAIL)/QH
      IF (NPTS.GE.MAXLIM1) MAXLIM1=NPTS

      NLIMIT=KEY1
      HP6=(C1/C)*D01AHF(H0,HR,ER0,NPTS,RELLER,BRADP,NLIMIT,IFAIL)/QH
      IF (NPTS.GE.MAXLIM1) MAXLIM1=NPTS

      NLIMIT=KEY1
      HP7=(C1)*D01AHF(HR,H1,ER0,NPTS,RELLER,BINDP,NLIMIT,IFAIL)/QH
      IF (NPTS.GE.MAXLIM1) MAXLIM1=NPTS

      NLIMIT=KEY1
      HP8=(C1/C)*D01AHF(HR,H1,ER0,NPTS,RELLER,BRADP,NLIMIT,IFAIL)/QH
      IF (NPTS.GE.MAXLIM1) MAXLIM1=NPTS
C====================================================================     
      
	ZS=DABS(ZS)
      RETURN
667   EV1=0.   
      EV2=0.
      EV3=0.
      EV4=0.
      EV5=0.
      EV6=0.
      EV7=0.
      EV8=0.
      EV9=0.
      EV10=0.
      EV11=0.
      EV12=0.
      HP1=0.
      HP2=0.
      HP3=0.
      HP4=0.
      HP5=0.
      HP6=0.
      HP7=0.
      HP8=0.
      RETURN
      END

C	============================= ELECTROSTATIC COMPONENT OF Ez =========
      REAL*8 FUNCTION EELP(Z)
C	=============================
      IMPLICIT NONE
	REAL*8 INTPULSE,P
      REAL*8 MU,C,PI,EPS,V,VC,W,H,ZLAM,R,ZS,TT,HR,RR,RG,ALPHA_TOWER
      REAL*8 ICI,ICII,Z,RZ,QEV,QH,QEH
	REAL*8 A1,A2,A3,B1,B2,B3   
	INTEGER NC,J,MODEL,MPROBLEM
      EXTERNAL INTPULSE,P

      COMMON /COM1/ MU,C,PI,EPS
      COMMON /COM2/ V,VC,W,H,ZLAM
      COMMON /COM3/ R,ZS
      COMMON /COMT/ TT
      COMMON /COM6/ HR,RR,RG,ALPHA_TOWER,NC,MODEL,MPROBLEM
      COMMON /SCALE/ QEV,QH,QEH
      
	RZ=DSQRT((R*R)+((Z-ZS)**2.))
     	ICI=0.
      ICII=0.
	IF (Z.LT.HR) THEN 
C===========1	  
	DO 10 J=0,NC,1		 
       B1=TT-RZ/C-(HR-Z)/C-2.*J*HR/C
          IF (B1.GE.0.) THEN
		 A1=0.
           ICII=INTPULSE(A1,B1)
           ICI=ICI+ICII*(RR**J)*(RG**J)*EXP(-J*HR*ALPHA_TOWER)
          END IF
       B2=B1-2.*Z/C
          IF (B2.GE.0.) THEN
	     A2=0.
           ICII=INTPULSE(A2,B2)
           ICI=ICI+ICII*(RR**J)*(RG**(J+1))*EXP(-(J*HR+Z)*ALPHA_TOWER)
          END IF
10    CONTINUE
	ICI=ICI*(1-RR)
	ELSE
C===========Z GREATER THAN HR 
	B1=TT-RZ/C-(Z-HR)/V
	    IF (B1.GE.0.) THEN
	      A3=(Z-HR)/V-(Z-HR)/VC
	      B3=TT-RZ/C-(Z-HR)/VC
	      ICII=INTPULSE(A3,B3)
		  ICI=ICI+ICII*P(Z)
		END IF	    
C===========2
	B2=TT-RZ/C-(Z-HR)/C
	    IF (B2.GE.0.) THEN
	     A1=0.
		 ICII=INTPULSE(A1,B2)
		 ICI=ICI+ICII*(-RR)
		END IF  
C===========3
	   DO 20 J=0,NC,1
      B3=TT-RZ/C-(Z+HR)/C-2.*J*HR/C
          IF (B3.GE.0.) THEN
	     A1=0.
		 ICII=INTPULSE(A1,B3)     
		 ICI=ICI+ICII*(1-RR)*(1+RR)*(RR**J)*(RG**(J+1))*EXP(-(J+1)*HR*ALPHA_TOWER)
          END IF
20       CONTINUE
C================================	
	
	END IF
      EELP=(2.*(ZS-Z)*(ZS-Z)-R*R)*ICI/(RZ**5)*QEV
      RETURN
      END

C	============================= INDUCTION COMPONENT OF Ez ==========
      REAL*8 FUNCTION EINDP(Z)
C	=============================
      IMPLICIT NONE
	REAL*8 IPULSE,P
      REAL*8 MU,C,PI,EPS,V,VC,W,H,ZLAM,R,ZS,TT,HR,RR,RG,ALPHA_TOWER              
      REAL*8 ICI,ICII,Z,RZ,QEV,QH,QEH
	REAL*8 B1,B2,B3 
      INTEGER NC,J,MODEL,MPROBLEM
      EXTERNAL IPULSE,P
      
      COMMON /COM1/ MU,C,PI,EPS
      COMMON /COM2/ V,VC,W,H,ZLAM
      COMMON /COM3/ R,ZS
      COMMON /COMT/ TT
      COMMON /COM6/ HR,RR,RG,ALPHA_TOWER,NC,MODEL,MPROBLEM      
      COMMON /SCALE/ QEV,QH,QEH
      
      RZ=DSQRT((R*R)+((Z-ZS)**2.))
      ICI=0.
      ICII=0.
	IF (Z.LT.HR) THEN 
C===========1	
	 DO 10 J=0,NC,1
       B1=TT-RZ/C-(HR-Z)/C-2.*J*HR/C
          IF (B1.GE.0.) THEN
          ICII=IPULSE(B1)
          ICI=ICI+ICII*(RR**J)*(RG**J)*EXP(-J*HR*ALPHA_TOWER)
          END IF
       B2=B1-2.*Z/C
          IF (B2.GE.0.) THEN
          ICII=IPULSE(B2)
          ICI=ICI+ICII*(RR**J)*(RG**(J+1))*EXP(-(J*HR+Z)*ALPHA_TOWER)
          END IF
10    CONTINUE
	ICI=ICI*(1-RR)
	ELSE
C===========Z GREATER THAN HR
	B1=TT-RZ/C-(Z-HR)/V
	B3=TT-RZ/C-(Z-HR)/VC
	    IF (B1.GE.0.) THEN
	    ICII=IPULSE(B3)
		ICI=ICI+ICII*P(Z)
		END IF
C===========2
	B1=TT-RZ/C-(Z-HR)/C
	    IF (B1.GE.0.) THEN
		 ICII=IPULSE(B1)
		 ICI=ICI+ICII*(-RR)
		END IF  
C===========3
	   DO 20 J=0,NC,1
      B1=TT-RZ/C-(Z+HR)/C-2.*J*HR/C
          IF (B1.GE.0.) THEN
		 ICII=IPULSE(B1)     
		 ICI=ICI+ICII*(1-RR)*(1+RR)*(RR**J)*(RG**(J+1))*EXP(-(J+1)*HR*ALPHA_TOWER)
          END IF
20       CONTINUE
C================================	
	
	END IF

      EINDP=ICI*(2*(ZS-Z)*(ZS-Z)-R*R)/(RZ**4)*QEV
      RETURN
      END


C	============================= RADIATION COMPONENT OF Ez ===============
      REAL*8 FUNCTION ERADP(Z)
C	=============================
      IMPLICIT NONE
	REAL*8 DIPULSE,P
      REAL*8 MU,C,PI,EPS,V,VC,W,H,ZLAM,R,ZS,TT,HR,RR,RG,ALPHA_TOWER              
      REAL*8 ICI,ICII,Z,RZ,QEV,QH,QEH
	REAL*8 B1,B2,B3
	INTEGER NC,J,MODEL,MPROBLEM
      EXTERNAL DIPULSE,P
      
      COMMON /COM1/ MU,C,PI,EPS
      COMMON /COM2/ V,VC,W,H,ZLAM
      COMMON /COM3/ R,ZS
      COMMON /COMT/ TT
      COMMON /COM6/ HR,RR,RG,ALPHA_TOWER,NC,MODEL,MPROBLEM
      COMMON /SCALE/ QEV,QH,QEH
      
      RZ=DSQRT((R*R)+((Z-ZS)**2.))
      ICI=0.
      ICII=0.
	IF (Z.LT.HR) THEN 
C===========1	
	 DO 10 J=0,NC,1
       B1=TT-RZ/C-(HR-Z)/C-2.*J*HR/C
          IF (B1.GE.0.) THEN
          ICII=DIPULSE(B1)
          ICI=ICI+ICII*(RR**J)*(RG**J)*EXP(-J*HR*ALPHA_TOWER)
          END IF
       B2=B1-2.*Z/C
          IF (B2.GE.0.) THEN
          ICII=DIPULSE(B2)
          ICI=ICI+ICII*(RR**J)*(RG**(J+1))*EXP(-(J*HR+Z)*ALPHA_TOWER)
          END IF
10    CONTINUE
	ICI=ICI*(1-RR)
	ELSE
C===========Z GREATER THAN HR
	B1=TT-RZ/C-(Z-HR)/V
	B3=TT-RZ/C-(Z-HR)/VC
	    IF (B1.GE.0.) THEN
	    ICII=DIPULSE(B3)
		ICI=ICI+ICII*P(Z)
		END IF
C===========2
	B1=TT-RZ/C-(Z-HR)/C
	    IF (B1.GE.0.) THEN
		 ICII=DIPULSE(B1)
		 ICI=ICI+ICII*(-RR)
		END IF  
C===========3
	   DO 20 J=0,NC,1
      B1=TT-RZ/C-(Z+HR)/C-2.*J*HR/C
          IF (B1.GE.0.) THEN
		 ICII=DIPULSE(B1)     
		 ICI=ICI+ICII*(1-RR)*(1+RR)*(RR**J)*(RG**(J+1))*EXP(-(J+1)*HR*ALPHA_TOWER)
          END IF
20       CONTINUE
C================================	
	
	END IF

      ERADP=ICI/(RZ**3)*QEV
      RETURN
      END


C	============================= INDUCTION COMPONENT OF H ===========
      REAL*8 FUNCTION BINDP(Z)
C	=============================
      IMPLICIT NONE
	REAL*8 IPULSE,P
      REAL*8 MU,C,PI,EPS,V,VC,W,H,ZLAM,R,ZS,TT,HR,RR,RG,ALPHA_TOWER             
      REAL*8 ICI,ICII,Z,RZ,QEV,QH,QEH
	REAL*8 B1,B2,B3      
	INTEGER NC,J,MODEL,MPROBLEM
      EXTERNAL IPULSE,P
      
      COMMON /COM1/ MU,C,PI,EPS
      COMMON /COM2/ V,VC,W,H,ZLAM
      COMMON /COM3/ R,ZS
      COMMON /COMT/ TT
      COMMON /COM6/ HR,RR,RG,ALPHA_TOWER,NC,MODEL,MPROBLEM
      COMMON /SCALE/ QEV,QH,QEH
C	--------------------------------------------------      
	RZ=DSQRT((R*R)+((Z-ZS)**2.))  
	ICI=0.
      ICII=0.
	IF (Z.LT.HR) THEN 
C===========1	
	 DO 10 J=0,NC,1
       B1=TT-RZ/C-(HR-Z)/C-2.*J*HR/C
          IF (B1.GE.0.) THEN
          ICII=IPULSE(B1)
          ICI=ICI+ICII*(RR**J)*(RG**J)*EXP(-J*HR*ALPHA_TOWER)
          END IF
       B2=B1-2.*Z/C
          IF (B2.GE.0.) THEN
          ICII=IPULSE(B2)
          ICI=ICI+ICII*(RR**J)*(RG**(J+1))*EXP(-(J*HR+Z)*ALPHA_TOWER)
          END IF
10    CONTINUE
	ICI=ICI*(1-RR)
	ELSE
C===========Z GREATER THAN HR 
	B1=TT-RZ/C-(Z-HR)/V
	B3=TT-RZ/C-(Z-HR)/VC
	    IF (B1.GE.0.) THEN
	    ICII=IPULSE(B3)
		ICI=ICI+ICII*P(Z)
		END IF
C===========2
	B1=TT-RZ/C-(Z-HR)/C
	    IF (B1.GE.0.) THEN
		 ICII=IPULSE(B1)
		 ICI=ICI+ICII*(-RR)
		END IF  
C===========3
	   DO 20 J=0,NC,1
      B1=TT-RZ/C-(Z+HR)/C-2.*J*HR/C
          IF (B1.GE.0.) THEN
		 ICII=IPULSE(B1)     
		 ICI=ICI+ICII*(1-RR)*(1+RR)*(RR**J)*(RG**(J+1))*EXP(-(J+1)*HR*ALPHA_TOWER)
          END IF
20       CONTINUE
C================================	
	
	END IF

      BINDP=ICI/(RZ**3)*QH
      RETURN
      END


C	============================= RADIATION COMPONENT OF H =================
       REAL*8 FUNCTION BRADP(Z)
C	============================
      IMPLICIT NONE
	REAL*8 DIPULSE,P
      REAL*8 MU,C,PI,EPS,V,VC,W,H,ZLAM,R,ZS,TT,HR,RR,RG,ALPHA_TOWER           
      REAL*8 ICI,ICII,Z,RZ,QEV,QH,QEH
	REAL*8 B1,B2,B3      
	INTEGER NC,J,MODEL,MPROBLEM
      EXTERNAL DIPULSE,P
      
      COMMON /COM1/ MU,C,PI,EPS
      COMMON /COM2/ V,VC,W,H,ZLAM
      COMMON /COM3/ R,ZS
      COMMON /COMT/ TT
      COMMON /COM6/ HR,RR,RG,ALPHA_TOWER,NC,MODEL,MPROBLEM
      COMMON /SCALE/ QEV,QH,QEH
C     ------------------------------------------------      
      RZ=DSQRT((R*R)+((Z-ZS)**2.)) 
	ICI=0.
      ICII=0.
	IF (Z.LT.HR) THEN 
C===========1	
	 DO 10 J=0,NC,1
       B1=TT-RZ/C-(HR-Z)/C-2.*J*HR/C
          IF (B1.GE.0.) THEN
          ICII=DIPULSE(B1)
          ICI=ICI+ICII*(RR**J)*(RG**J)*EXP(-J*HR*ALPHA_TOWER)
          END IF
       B2=B1-2.*Z/C
          IF (B2.GE.0.) THEN
          ICII=DIPULSE(B2)
          ICI=ICI+ICII*(RR**J)*(RG**(J+1))*EXP(-(J*HR+Z)*ALPHA_TOWER)
          END IF
10    CONTINUE
	ICI=ICI*(1-RR)
	ELSE
C===========Z GREATER THAN HR
	B1=TT-RZ/C-(Z-HR)/V
	B3=TT-RZ/C-(Z-HR)/VC
	    IF (B1.GE.0.) THEN
	    ICII=DIPULSE(B3)
		ICI=ICI+ICII*P(Z)
		END IF
C===========2
	B1=TT-RZ/C-(Z-HR)/C
	    IF (B1.GE.0.) THEN
		 ICII=DIPULSE(B1)
		 ICI=ICI+ICII*(-RR)
		END IF  
C===========3
	   DO 20 J=0,NC,1
      B1=TT-RZ/C-(Z+HR)/C-2.*J*HR/C
          IF (B1.GE.0.) THEN
		 ICII=DIPULSE(B1)     
		 ICI=ICI+ICII*(1-RR)*(1+RR)*(RR**J)*(RG**(J+1))*EXP(-(J+1)*HR*ALPHA_TOWER)
          END IF
20       CONTINUE

C================================	
	
	END IF

      BRADP=ICI/((RZ**2))*QH
      RETURN
      END


C     ==========================================
C     CALCULATION OF RADIAL ELECTRIC FIELD Er
C     ==========================================
      SUBROUTINE HZONTAL(TE) 
    

	IMPLICIT NONE
      REAL*8 MU,C,PI,EPS,V,VC,W,H,ZLAM,R1,FV,IFRONT
      REAL*8 R,ZS,TT,HR,RR,RG,ALPHA_TOWER,ERR,ERR2
      REAL*8 EH1,EH2,EH3,EH4,EH5,EH6,EH7,EH8,EH9,EH10,EH11,EH12
      REAL*8 QEV,QH,QEH,D01AHF
      REAL*8 TE,T01, H0,H1,C1,C3,ER0,RELLER
      REAL*8 PHI,DELTAP,PHI1,PHI2
      REAL*8 EELRP,EINDRP,ERADRP,EH1TO,EH2TO
      INTEGER NPTS,NLIMIT,IFAIL,NC,MAXLIM1,KEY1,MAXLIM2,KEY2,MODEL
	INTEGER MPROBLEM

      COMMON /COM1/ MU,C,PI,EPS,IFRONT
      COMMON /COM2/ V,VC,W,H,ZLAM
      COMMON /COM3/ R,ZS,R1,FV
      COMMON /COMT/ TT
      COMMON /COM6/ HR,RR,RG,ALPHA_TOWER,NC,MODEL,MPROBLEM
      COMMON /COM7/ ERR,ERR2
      COMMON /EH/ EH1,EH2,EH3,EH4,EH5,EH6,EH7,EH8,EH9,EH10,EH11,EH12,
	1			EH1TO,EH2TO
      COMMON /SCALE/ QEV,QH,QEH
      COMMON /INTEG/ MAXLIM1,KEY1,MAXLIM2,KEY2
      
      EXTERNAL EELRP,EINDRP,ERADRP
      EXTERNAL D01AHF
      
      ZS=DABS(ZS)
      T01=DSQRT(R*R+(HR-ZS)*(HR-ZS))/C

      TT=TE+T01

      IF (TT.EQ.T01) GOTO 668
      IF (ZS.EQ.(0.0)) GOTO 668
	C1=R/(4*PI)
      C3=1./(4.*PI*EPS)
      ER0=ERR2
      
      
C===========DETERMINATION OF THE UPPER INTEGRATION LIMIT
      PHI=1./(V*V)-1./(C*C)
      PHI1=ZS/(C*C)-(TT/V)-HR/(V*V)
      PHI2=(TT+HR/V)**2.-(R/C)**2.-(ZS/C)**2.
      DELTAP=DABS(DSQRT(PHI1*PHI1-PHI*PHI2))
      H1=(-PHI1-DELTAP)/PHI
      IF (H1.LT.HR)H1=HR
      IF (H1.GT.H) H1=H
C============DETERMINATION OF THE LOWER INTEGRATION LIMIT
      PHI1=(ZS/(C*C)+(TT/C)-HR/(C*C))*2.
      PHI2=-(TT-HR/C)**2.+(R/C)**2.+(ZS/C)**2.
      H0=PHI2/PHI1
      IF (H0.GT.HR) H0=HR
      IF (H0.LT.0.) H0=0.
      
      IFAIL=0

      NLIMIT=KEY2
      EH1=(C3/C)*D01AHF(H0,HR,ER0,NPTS,RELLER,EINDRP
     1                   ,NLIMIT,IFAIL)*3.*R/QEH
      IF (NPTS.GE.MAXLIM2) MAXLIM2=NPTS
      
      NLIMIT=KEY2
      EH4=(C3/C)*D01AHF(HR,H1,ER0,NPTS,RELLER,EINDRP
     1                   ,NLIMIT,IFAIL)*3.*R/QEH
      IF (NPTS.GE.MAXLIM2) MAXLIM2=NPTS

      NLIMIT=KEY2
      EH2=(C3/(C*C))*D01AHF(H0,HR,ER0,NPTS,RELLER,ERADRP,NLIMIT
     1                     ,IFAIL)*R/QEH
      IF (NPTS.GE.MAXLIM2) MAXLIM2=NPTS

      NLIMIT=KEY2
      EH5=(C3/(C*C))*D01AHF(HR,H1,ER0,NPTS,RELLER,ERADRP,NLIMIT
     1                     ,IFAIL)*R/QEH
      IF (NPTS.GE.MAXLIM2) MAXLIM2=NPTS

      NLIMIT=KEY2
      EH3=C3*D01AHF(H0,HR,ER0,NPTS,RELLER,EELRP,NLIMIT
     1               ,IFAIL)*3.*R/QEH
      IF (NPTS.GE.MAXLIM2) MAXLIM2=NPTS

      NLIMIT=KEY2
      EH6=C3*D01AHF(HR,H1,ER0,NPTS,RELLER,EELRP,NLIMIT
     1               ,IFAIL)*3.*R/QEH
      IF (NPTS.GE.MAXLIM2) MAXLIM2=NPTS

C============DETERMINATION OF THE UPPER INTEGRATION LIMIT
      ZS=-DABS(ZS)
      PHI=1./(V*V)-1./(C*C)
      PHI1=ZS/(C*C)-(TT/V)-HR/(V*V)
      PHI2=(TT+HR/V)**2.-(R/C)**2.-(ZS/C)**2.
      DELTAP=DABS(DSQRT(PHI1*PHI1-PHI*PHI2))
      H1=(-PHI1-DELTAP)/PHI
      IF (H1.LT.HR) H1=HR
      IF (H1.GT.H) H1=H

C=============DETERMINATION OF THE LOWER INTEGRATION LIMIT
      PHI1=(ZS/(C*C)+(TT/C)-HR/(C*C))*2.
      PHI2=-(TT-HR/C)**2.+(R/C)**2.+(ZS/C)**2.
      H0=PHI2/PHI1
      IF (H0.GT.HR) H0=HR
      IF (H0.LT.0.) H0=0.
      


      NLIMIT=KEY2
      EH7=-(C3/C)*D01AHF(H0,HR,ER0,NPTS,RELLER,EINDRP
     1            ,NLIMIT,IFAIL)*3.*R/QEH
      IF (NPTS.GE.MAXLIM2) MAXLIM2=NPTS

      NLIMIT=KEY2
      EH10=-(C3/C)*D01AHF(HR,H1,ER0,NPTS,RELLER,EINDRP
     1            ,NLIMIT,IFAIL)*3.*R/QEH
      IF (NPTS.GE.MAXLIM2) MAXLIM2=NPTS

      NLIMIT=KEY2
      EH8=-(C3/(C*C))*D01AHF(H0,HR,ER0,NPTS,RELLER,ERADRP
     1            ,NLIMIT,IFAIL)*R/QEH
      IF (NPTS.GE.MAXLIM2) MAXLIM2=NPTS

      NLIMIT=KEY2
      EH11=-(C3/(C*C))*D01AHF(HR,H1,ER0,NPTS,RELLER,ERADRP
     1            ,NLIMIT,IFAIL)*R/QEH
      IF (NPTS.GE.MAXLIM2) MAXLIM2=NPTS
 
      NLIMIT=KEY2
      EH9=-C3*D01AHF(H0,HR,ER0,NPTS,RELLER,EELRP
     1          ,NLIMIT,IFAIL)*3.*R/QEH
      IF (NPTS.GE.MAXLIM2) MAXLIM2=NPTS

      NLIMIT=KEY2
      EH12=-C3*D01AHF(HR,H1,ER0,NPTS,RELLER,EELRP
     1          ,NLIMIT,IFAIL)*3.*R/QEH
      IF (NPTS.GE.MAXLIM2) MAXLIM2=NPTS

      ZS=DABS(ZS)
      
	RETURN
668   EH1=0.
      EH2=0.
      EH3=0.
      EH4=0.  
      EH5=0.
      EH6=0.
      EH7=0.
      EH8=0.
      EH9=0.
      EH10=0.
      EH11=0.
      EH12=0.
      RETURN
      END


C	============================= ELECTROSTATIC COMPONENT OF Er ==============
      REAL*8 FUNCTION EELRP(Z)
C	===========================
      IMPLICIT NONE
	REAL*8 INTPULSE,P
      REAL*8 MU,C,PI,EPS,V,VC,W,H,ZLAM,R,ZS,TT,HR,RR,RG,ALPHA_TOWER             
      REAL*8 ICI,ICII,Z,RZ,QEV,QH,QEH
	REAL*8 A1,A2,A3,B1,B2,B3
	INTEGER NC,J,MODEL,MPROBLEM
      EXTERNAL INTPULSE,P
      
      COMMON /COM1/ MU,C,PI,EPS
      COMMON /COM2/ V,VC,W,H,ZLAM
      COMMON /COM3/ R,ZS
      COMMON /COMT/ TT
      COMMON /COM6/ HR,RR,RG,ALPHA_TOWER,NC,MODEL,MPROBLEM
      COMMON /SCALE/ QEV,QH,QEH
      
      RZ=DSQRT((R*R)+((Z-ZS)**2.))
      ICI=0.
      ICII=0.
	IF (Z.LT.HR) THEN 
C===========1	
	 DO 10 J=0,NC,1
       B1=TT-RZ/C-(HR-Z)/C-2.*J*HR/C
          IF (B1.GT.0.) THEN
	    A1=0.
          ICII=INTPULSE(A1,B1)
          ICI=ICI+ICII*(RR**J)*(RG**J)*EXP(-J*HR*ALPHA_TOWER)
          END IF
       B2=B1-2.*Z/C
          IF (B2.GT.0.) THEN
	    A2=0.
          ICII=INTPULSE(A2,B2)
          ICI=ICI+ICII*(RR**J)*(RG**(J+1))*EXP(-(J*HR+Z)*ALPHA_TOWER)
          END IF
10    CONTINUE
	ICI=ICI*(1-RR)
	ELSE
C===========Z GREATER THAN HR
	B1=TT-RZ/C-(Z-HR)/V
	      IF (B1.GT.0.) THEN
	       A3=(Z-HR)/VC-(Z-HR)/V
	       B3=TT-RZ/C-(Z-HR)/VC
	       ICII=INTPULSE(A3,B3)
		   ICI=ICI+ICII*P(Z)
		  END IF
C===========2
	B1=TT-RZ/C-(Z-HR)/C
	    IF (B1.GT.0.) THEN
	     A1=0.
		 ICII=INTPULSE(A1,B1)
		 ICI=ICI+ICII*(-RR)
		END IF  
C===========3
	   DO 20 J=0,NC,1
      B1=TT-RZ/C-(Z+HR)/C-2.*J*HR/C
          IF (B1.GT.0.) THEN
	     A1=0.
		 ICII=INTPULSE(A1,B1)     
		 ICI=ICI+ICII*(1-RR)*(1+RR)*(RR**J)*(RG**(J+1))*EXP(-(J+1)*HR*ALPHA_TOWER)
          END IF
20       CONTINUE
C================================	
	
	END IF

      EELRP=(ZS-Z)*ICI/(RZ**5)*QEH
      RETURN
      END


C	============================= INDUCTION COMPONENT OF Er =========
      REAL*8 FUNCTION EINDRP(Z)
C	==============================
      IMPLICIT NONE
	REAL*8 IPULSE,P
      REAL*8 MU,C,PI,EPS,V,VC,W,H,ZLAM,R,ZS,TT,HR,RR,RG,ALPHA_TOWER             
      REAL*8 ICI,ICII,Z,B1,B2,RZ,QEV,QH,QEH
	INTEGER NC,J,MODEL,MPROBLEM
      EXTERNAL IPULSE,P
      
      COMMON /COM1/ MU,C,PI,EPS
      COMMON /COM2/ V,VC,W,H,ZLAM
      COMMON /COM3/ R,ZS
      COMMON /COMT/ TT
      COMMON /COM6/ HR,RR,RG,ALPHA_TOWER,NC,MODEL,MPROBLEM
      COMMON /SCALE/ QEV,QH,QEH
      
      RZ=DSQRT((R*R)+((Z-ZS)**2.))
      ICI=0.
      ICII=0.
	IF (Z.LT.HR) THEN 
C===========1	
	 DO 10 J=0,NC,1
       B1=TT-RZ/C-(HR-Z)/C-2.*J*HR/C
          IF (B1.GT.0.) THEN
          ICII=IPULSE(B1)
          ICI=ICI+ICII*(RR**J)*(RG**J)*EXP(-J*HR*ALPHA_TOWER)
          END IF
       B2=B1-2.*Z/C
          IF (B2.GT.0.) THEN
          ICII=IPULSE(B2)
          ICI=ICI+ICII*(RR**J)*(RG**(J+1))*EXP(-(J*HR+Z)*ALPHA_TOWER)
          END IF
10    CONTINUE
	ICI=ICI*(1-RR)
	ELSE
C===========Z GREATER THAN HR 
	B1=TT-RZ/C-(Z-HR)/VC
	    IF (B1.GT.0.) THEN
	    ICII=IPULSE(B1)
		ICI=ICI+ICII*P(Z)
		END IF
C===========2
	B1=TT-RZ/C-(Z-HR)/C
	    IF (B1.GT.0.) THEN
		 ICII=IPULSE(B1)
		 ICI=ICI+ICII*(-RR)
		END IF  
C===========3
	   DO 20 J=0,NC,1
      B1=TT-RZ/C-(Z+HR)/C-2.*J*HR/C
          IF (B1.GT.0.) THEN
		 ICII=IPULSE(B1)     
		 ICI=ICI+ICII*(1-RR)*(1+RR)*(RR**J)*(RG**(J+1))*EXP(-(J+1)*HR*ALPHA_TOWER)
          END IF
20       CONTINUE
C================================	
	
	END IF

      EINDRP=(ZS-Z)*ICI/(RZ**4)*QEH
      RETURN
      END


C	============================= RADIATION COMPONENT OF Er ========
      REAL*8 FUNCTION ERADRP(Z)
C	=============================
      IMPLICIT NONE
	REAL*8 DIPULSE,P
      REAL*8 MU,C,PI,EPS,V,VC,W,H,ZLAM,R,ZS,TT,HR,RR,RG,ALPHA_TOWER                
      REAL*8 ICI,ICII,Z,B1,B2,RZ,QEV,QH,QEH
	INTEGER NC,J,MODEL,MPROBLEM
      EXTERNAL DIPULSE,P
     
      COMMON /COM1/ MU,C,PI,EPS
      COMMON /COM2/ V,VC,W,H,ZLAM
      COMMON /COM3/ R,ZS
      COMMON /COMT/ TT
      COMMON /COM6/ HR,RR,RG,ALPHA_TOWER,NC,MODEL,MPROBLEM
      COMMON /SCALE/ QEV,QH,QEH
      
      RZ=DSQRT((R*R)+((Z-ZS)**2.))
      ICI=0.
      ICII=0.
	IF (Z.LT.HR) THEN 
C===========1	
	DO 10 J=0,NC,1
       B1=TT-RZ/C-(HR-Z)/C-2.*J*HR/C
          IF (B1.GT.0.) THEN
          ICII=DIPULSE(B1)
          ICI=ICI+ICII*(RR**J)*(RG**J)*EXP(-J*HR*ALPHA_TOWER)
          END IF
       B2=B1-2.*Z/C
          IF (B2.GT.0.) THEN
          ICII=DIPULSE(B2)
          ICI=ICI+ICII*(RR**J)*(RG**(J+1))*EXP(-(J*HR+Z)*ALPHA_TOWER)
          END IF
10    CONTINUE
	ICI=ICI*(1-RR)
	ELSE
C===========Z GREATER THAN HR 
	B1=TT-RZ/C-(Z-HR)/VC
	    IF (B1.GT.0.) THEN
	    ICII=DIPULSE(B1)
		ICI=ICI+ICII*P(Z)
		END IF
C===========2
	B1=TT-RZ/C-(Z-HR)/C
	    IF (B1.GT.0.) THEN
		 ICII=DIPULSE(B1)
		 ICI=ICI+ICII*(-RR)
		END IF  
C===========3
	   DO 20 J=0,NC,1
      B1=TT-RZ/C-(Z+HR)/C-2.*J*HR/C
          IF (B1.GT.0.) THEN
		 ICII=DIPULSE(B1)     
		 ICI=ICI+ICII*(1-RR)*(1+RR)*(RR**J)*(RG**(J+1))*EXP(-(J+1)*HR*ALPHA_TOWER)
          END IF
20       CONTINUE
C================================	
	
	END IF

      ERADRP=(ZS-Z)*ICI/(RZ**3)*QEH
      RETURN
      END


C     ==============================================================      
C                       CURRENT FUNCTIONS 
C     ==============================================================

      REAL*8 FUNCTION IPULSE(TA)
C	============================== 2Heidler+2exp
      IMPLICIT NONE
	REAL*8 MU,C,PI,EPS,V,VC,W,H,ZLAM,IW,ETA,TAU1,TAU2,IW2,ETA2,TAU12
     1       ,TAU22,ALPHA,BETA,I01,GAMMA,DELTA,I02,I0B,T0B
      REAL*8 FK,FKK,TA
      INTEGER NW2,NW
      
      COMMON /COM1/ MU,C,PI,EPS
      COMMON /COM2/ V,VC,W,H,ZLAM
      COMMON /COM4/ IW,ETA,TAU1,TAU2,IW2,ETA2,TAU12,TAU22,NW2,NW
      COMMON /COM5/ ALPHA,BETA,I01,GAMMA,DELTA,I02
	COMMON /COMBARB/ I0B,T0B

      
      FK=(TA/TAU1)**NW/(1.+((TA/TAU1)**NW))
      FKK=(TA/TAU12)**NW2/(1.+((TA/TAU12)**NW2)) 
      IPULSE=(IW/ETA)*FK*DEXP(-TA/TAU2)
      IPULSE=IPULSE+(IW2/ETA2)*FKK*DEXP(-TA/TAU22)
      IPULSE=IPULSE+I01*(DEXP(-ALPHA*TA)-DEXP(-BETA*TA))
	IPULSE=IPULSE+I02*(DEXP(-GAMMA*TA)-DEXP(-DELTA*TA))
	IF (TA.LE.T0B) THEN
	IPULSE=IPULSE+I0B*TA/T0B
	ELSE
	IPULSE=IPULSE+I0B
	ENDIF
	      
      RETURN
      END

      REAL*8 FUNCTION DIPULSE(TA)
C===============================d(2Heidler+2exp)/dt
      IMPLICIT NONE
	REAL*8 MU,C,PI,EPS,V,VC,W,H,ZLAM,IW,ETA,TAU1,TAU2,IW2,ETA2,TAU12
     1       ,TAU22,ALPHA,BETA,I01,GAMMA,DELTA,I02,I0B,T0B
      REAL*8 FK1,FK2,FKK1,FKK2,TA
      INTEGER NW2,NW

      COMMON /COM1/ MU,C,PI,EPS
      COMMON /COM2/ V,VC,W,H,ZLAM
      COMMON /COM4/ IW,ETA,TAU1,TAU2,IW2,ETA2,TAU12,TAU22,NW2,NW
      COMMON /COM5/ ALPHA,BETA,I01,GAMMA,DELTA,I02
	COMMON /COMBARB/ I0B,T0B

      
      FK2=NW/(1.+((TA/TAU1)**NW))-TA/TAU2
      FK1=((TA/TAU1)**(NW-1))/(1.+((TA/TAU1)**NW))
      FKK2=NW2/(1.+((TA/TAU12)**NW2))-TA/TAU22
      FKK1=((TA/TAU12)**(NW2-1))/(1.+((TA/TAU12)**NW2))
      DIPULSE=(IW/(ETA*TAU1))*FK1*DEXP(-TA/TAU2)*FK2
      DIPULSE=DIPULSE+(IW2/(ETA2*TAU12))*FKK1*DEXP(-TA/TAU22)*FKK2
      DIPULSE=DIPULSE+I01*(BETA*DEXP(-BETA*TA)-ALPHA*DEXP(-ALPHA*TA))
      DIPULSE=DIPULSE+I02*(DELTA*DEXP(-DELTA*TA)-GAMMA*DEXP(-GAMMA*TA))
	IF (TA.LE.T0B) THEN
	DIPULSE=DIPULSE+I0B/T0B
	ENDIF

      RETURN
      END


       
      REAL*8 FUNCTION INTPULSE(A,B)
C===============================integrale(2Heidler+2exp)dt 0-->t 
	IMPLICIT NONE
      REAL*8 IW,ETA,TAU1,TAU2,IW2,ETA2,TAU12,TAU22
      REAL*8 ALPHA,BETA,I01,GAMMA,DELTA,I02,I0B,T0B
      REAL*8 V,VC,W,H,ZLAM
      REAL*8 SUM,SUM1,SUM2
      REAL*8 PULSE,A,B,HS,E1,TA 
      INTEGER NN,K,K1,K2,NW2,NW,NSIMP,NN1,NN2,K12,K22

      COMMON /COM4/ IW,ETA,TAU1,TAU2,IW2,ETA2,TAU12,TAU22,NW2,NW
      COMMON /COM5/ ALPHA,BETA,I01,GAMMA,DELTA,I02
      COMMON /COM2/ V,VC,W,H,ZLAM
      COMMON /COM8/ NSIMP
	COMMON /COMBARB/ I0B,T0B
      
      PULSE(TA)=(IW/ETA)*(TA/TAU1)**NW/(1.+((TA/TAU1)**NW))
     1       *DEXP(-TA/TAU2)+(IW2/ETA2)
     1       *(TA/TAU12)**NW2/(1.+((TA/TAU12)**NW2))*DEXP(-TA/TAU22) 
	
      HS=(B-A)/DFLOAT(NSIMP)
      NN=NSIMP/2-1
      SUM=PULSE(A)+4.D0*PULSE(A+(NSIMP-1)*HS)
      SUM1=0.D0
      SUM2=0.D0
      DO 109 K=1,NN
      K1=2*K
      K2=2*K-1
      SUM1=SUM1+2.D0*PULSE(A+(DFLOAT(K1))*HS)
      SUM2=SUM2+4.D0*PULSE(A+(K2)*HS)
109   CONTINUE
      SUM=SUM+SUM1+SUM2+PULSE(B)	

      INTPULSE=HS*SUM/3.D0
      INTPULSE=INTPULSE+I01*((DEXP(-ALPHA*A)-DEXP(-ALPHA*B))/ALPHA-
	1                       (DEXP(-BETA*A)-DEXP(-BETA*B))/BETA)
     1			     +I02*((DEXP(-GAMMA*A)-DEXP(-GAMMA*B))/GAMMA-
	1		               (DEXP(-DELTA*A)-DEXP(-DELTA*B))/DELTA)

	IF ((A.LT.T0B).AND.(B.LE.T0B)) THEN
	INTPULSE=INTPULSE+0.5D0*I0B*(B**2-A**2)/T0B
	ELSEIF ((A.LE.T0B).AND.(B.GT.T0B)) THEN
	INTPULSE=INTPULSE+0.5D0*I0B*(T0B**2-A**2)/T0B+I0B*(B-T0B)
	ELSEIF((A.GE.T0B).AND.(B.GE.T0B)) THEN
	INTPULSE=INTPULSE+I0B*(B-A)
	ENDIF


      RETURN
      END


C======================= current at height QUO at time TY and its derivative

      SUBROUTINE CURRENT (TY,QUO,CUR,DCUR)
      
	IMPLICIT NONE
      REAL*8 HR,RR,RG,ALPHA_TOWER,V,VC,W,H,ZLAM,RZ
      REAL*8 CUR,DCUR,TY,QUO
	REAL*8 B1,B2,B3 
      REAL*8 DIPULSE,IPULSE,P
      REAL*8 MU,C,PI,EPS
	INTEGER J,NC,MODEL,MPROBLEM
      
	COMMON /COM9/ RZ
	COMMON /COM6/ HR,RR,RG,ALPHA_TOWER,NC,MODEL,MPROBLEM
      COMMON /COM2/ V,VC,W,H,ZLAM
	COMMON /COM1/ MU,C,PI,EPS
	EXTERNAL IPULSE,DIPULSE,P
      
	CUR=0.D0
	DCUR=0.D0

      IF (QUO.LE.HR) THEN 
C===========1	
	DO 10 J=0,NC,1
       B1=TY-RZ/C-(HR-QUO)/C-2.*J*HR/C
          IF (B1.GE.0.) THEN
           CUR=CUR+IPULSE(B1)*(RR**J)*(RG**J)*EXP(-J*HR*ALPHA_TOWER)
           DCUR=DCUR+DIPULSE(B1)*(RR**J)*(RG**J)*EXP(-J*HR*ALPHA_TOWER)
		END IF
       B2=B1-2.*QUO/C
          IF (B2.GE.0.) THEN
           CUR=CUR+IPULSE(B2)*(RR**J)*(RG**(J+1))*EXP(-(J*HR+QUO)*ALPHA_TOWER)
           DCUR=DCUR+DIPULSE(B2)*(RR**J)*(RG**(J+1))*EXP(-(J*HR+QUO)*ALPHA_TOWER)
		END IF
10    CONTINUE
	CUR=CUR*(1-RR)
	DCUR=DCUR*(1-RR)
	ELSE
C===========Z GREATER THAN HR 
	B1=TY-RZ/C-(QUO-HR)/VC
	    IF (B1.GE.0) THEN
	     IF (MODEL .EQ. 5) THEN
C	      B1 = B1 + (QUO-HR)/VC + (QUO-HR)/C
	     END IF
		 CUR=CUR+IPULSE(B1)*P(QUO) 
	     DCUR=DCUR+DIPULSE(B1)*P(QUO)
		END IF
C===========2
	B1=TY-RZ/C-(QUO-HR)/C
	    IF (B1.GE.0.) THEN
		 CUR=CUR+IPULSE(B1)*(-RR)
		 DCUR=DCUR+DIPULSE(B1)*(-RR)
		END IF  
C===========3
	DO 20 J=0,NC,1
      B1=TY-RZ/C-(QUO+HR)/C-2.*J*HR/C
          IF (B1.GE.0.) THEN
		 CUR=CUR+(1-RR)*(1+RR)*(RR**J)*(RG**(J+1))*IPULSE(B1)*EXP(-(J+1)*HR*ALPHA_TOWER)    
		 DCUR=DCUR+(1-RR)*(1+RR)*(RR**J)*(RG**(J+1))*DIPULSE(B1)*EXP(-(J+1)*HR*ALPHA_TOWER)
          END IF
20    CONTINUE
	
	END IF
999   RETURN
	END
	
	REAL*8 FUNCTION FF(T,SIGMA,EPSILR)
C     ========================================
C       FORMULA FOR CORRECTIVE TERM FOR HORIZONTAL FIELD (APPROACH Rubinstein)
	IMPLICIT NONE
      COMMON /COM1/ MU,C0,PI,EPS
      REAL*8  Y,X,S0,P,T,S1,SIGMA,EPSILR,EPS,
	1	    PI,MU,Z0,C0

        P=SIGMA/(2.*EPS*EPSILR)
        Z0=DSQRT(MU/EPS)
        X=P*T
        Y=X/3.75
C      GOTO 444 old formula RUBINSTEIN
	IF (DABS(X).LT.3.75)THEN

        S0=1+3.5156229*(Y**2)+3.0899424*(Y**4)+
	1	1.2067492*(Y**6)+
	1	0.2659732*(Y**8)+
	1	0.0360768*(Y**10)+
	1	0.0045813*(Y**12)
        S1=X*(0.5+0.87890594*Y*Y+0.51498869*
	1	(Y**4)+0.15084934*(Y**6)
	1	+0.02658733*(Y**8)
	1	+0.00301532*(Y**10)
	1	+0.00032411*(Y**12))
        FF=Z0*T*DEXP(-X)*(S0+S1)/DSQRT(EPSILR)
        ELSE
        S0=(0.39894228+
	1	0.01328592/Y
	1	+0.00225319/(Y**2)
	1	-0.00157565
	1	/(Y**3)+0.00916281/(Y**4)
	1	-0.02057706/(Y**5)
	1	+0.02635537/(Y**6)
	1	-0.01647633/(Y**7)
	1	+0.00392377/(Y**8))/DSQRT(X)


        S1=(0.39894228-0.03988024/Y-
	1	0.00362018/(Y*Y)+0.00163801
	1	/(Y*Y*Y)-.01031555/(Y**4)
	1	+0.02282967/(Y**5)-
	1	.02895312/(Y**6)+0.01787654/(Y**7)
	1	-0.00420059/(Y**8))/DSQRT(X)


        FF=Z0*T*(S0+S1)/DSQRT(EPSILR)
        END IF
C444   FF=2.0*T*DSQRT(MU/(PI*SIGMA*T))     
	  
	  RETURN
        END


	REAL*8 FUNCTION P(ZZ)
C	============================== CHOOSING THE ATTENUATION FACTOR P(z) =====
      IMPLICIT NONE
      
      REAL*8 MU,C,PI,EPS,V,VC,W,H,ZLAM,R,ZS,TT,HR,RR,RG,ALPHA_TOWER              
      REAL*8 Z,A1,A2,RZ,ZZ
      INTEGER NC,J,MODEL,MPROBLEM
      EXTERNAL IPULSE
      
      COMMON /COM1/ MU,C,PI,EPS
      COMMON /COM2/ V,VC,W,H,ZLAM
      COMMON /COM3/ R,ZS
      COMMON /COMT/ TT
      COMMON /COM6/ HR,RR,RG,ALPHA_TOWER,NC,MODEL,MPROBLEM
    

	IF ((MODEL.EQ.1).OR.(MODEL.EQ.4).OR.(MODEL.EQ.5)) P=1
      IF (MODEL.EQ.2) P=(1-(ZZ-HR)/H)
 	IF (MODEL.EQ.3) P=DEXP(-(ZZ-HR)*ZLAM)
	
	RETURN
      END


	   	  	
