
C ======================= current at height QUO at time TY and its derivative

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
