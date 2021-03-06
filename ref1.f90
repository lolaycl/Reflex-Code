
     SUBROUTINE ref1(X)

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

   

        RETURN
END SUBROUTINE ref1


