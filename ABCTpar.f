      SUBROUTINE ABCTpar(IIparal)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INTEGER IIparal

      CHARACTER*1  cpar
      DIMENSION JTEMP(40),JSHIFT(40)
     
      INCLUDE 'PRIVCOM21.FOR'      

      
      INCLUDE 'PRIVCOM.FOR'         
      INCLUDE 'PRIVCOM1.FOR'
      INCLUDE 'PRIVCOM6.FOR'
      INCLUDE 'PRIVCOM8.FOR'
      INCLUDE 'PRIVCOM9.FOR'
      INCLUDE 'PRIVCOM13.FOR'
      INCLUDE 'PRIVCOM14.FOR'  
C     Input commons (do not change with energy) 
      INCLUDE 'PRIVCOM10.FOR'
      INCLUDE 'PRIVCOM12.FOR'
      INCLUDE 'PRIVCOM15.FOR'
      INCLUDE 'PRIVCOM16D.FOR'
      INCLUDE 'PRIVCOM20.FOR'

      
      CHARACTER*8 PNAME
      INTEGER TID
!$    INTEGER OMP_GET_THREAD_NUM

      TID = 0
!$    TID = OMP_GET_THREAD_NUM()
!$    IF(MEPRI.NE.98)
!$   &    PRINT *, 'Thread',TID,' starting...','IIparal=',IIparal

      EN=EE(IIparal)

  777 MECHA=MCHAE(IIparal)
      
      PNAME='NEUTRONS'
      IF(MECHA.NE.0) PNAME=' PROTONS'       
  
C     CREATING LEVELS FOR (P,N) ANALOG STATES CALCULATIONS
     
      NURRR=NUR
            
      NURC=0
       HW=HWIS(1)
       AMB0=AMB0IS(1)
       AMG0=AMG0IS(1)
       GAM0=GAM0IS(1)
       BET0=BET0IS(1)
       BET4=BET4IS(1)
       BB42=BB42IS(1)
       GAMG=GAMGIS(1)
       DELG=DELGIS(1)
       BET3=BET3IS(1)
       ETO=ETOIS(1)
       AMUO=AMUOIS(1)
       HWO=HWOIS(1)
       BB32=BB32IS(1)
       GAMDE=GAMDIS(1)
       DPAR=DPARIS(1)
      GSHAPE=GSHAEIS(1)
      IF(MEPOT.GT.1) GO TO 638
      
      JSHIFT=0
      
      DO 601 I=1, NUR
      IF(MECHA.EQ.0.AND. NCAC(I).NE.NCAC(1)) GO TO 601
      NURC=NURC+1
      EL(NURC)=ELC(I)
      JO(NURC)=JOC(I)
      NPO(NURC)=NPOC(I)
      KO(NURC)=KOC(I)
      NCA(NURC)=NCAC(I)
      NUMB(NURC)=NUMBC(I)
      BETB(NURC)=BETBC(I)
      AIGS(NURC)=AGSIC(I)
      NTU(NURC)=NTUC(I)
      NNB(NURC)=NNBC(I)
      NNG(NURC)=NNGC(I)
      NNO(NURC)=NNOC(I) 
      ES(NURC)=EL(NURC) 
      JU(NURC)=JO(NURC)/2
      NPI(NURC)=NPO(NURC)
      IF (NNO(NURC).EQ.1) JSHIFT(NURC)=1
  601 CONTINUE
      NUR=NURC
      
      !JBASE=NINT(DBLE(JO(1))/4.0)*2
      
      IF(MOD(JO(1),2).GT.0) THEN
          JTEMP=JU
          !JU=NINT(DBLE(JO)/4.0)*2!-JBASE
          JU=NINT(DBLE(JO-JO(1))/4.0)*2+JSHIFT ! FOR GS, BETA, GAMMA, AND INV PARITY BANDS
          !!! ABNORMAL BAND SHOULD BE ASSIGNED SEPARATELY!!!
      END IF     
            
      EFFDEF=0.d0
      IF(MEDEF.GT.0.OR.MEAXI.EQ.1.OR.MEVOL.EQ.1) CALL OVLOPT
 
       DO IID=1,NUR
         DO JJD=IID,NUR
C             EFFDIS(IIS,IID,JJD,:)=EFFDEF(IID,JJD,:)
             EFFDEF(JJD,IID,:)=EFFDEF(IID,JJD,:)
         END DO
       END DO

      IF(MOD(JO(1),2).GT.0)  THEN
c          NUMBGS=NUMB(1)
c           DO IID=1,NUR
c             DO JJD=1,NUR
c                IF(NUMB(IID).NE.NUMBGS.OR.NUMB(JJD).NE.NUMBGS)
c     *                  EFFDEF(JJD,IID,:)=0.0
c             END DO
c           END DO          
          JU=JTEMP
      END IF         
  
      GO TO 639
  638 DO 602 I=1, NUR
      IF(MECHA.EQ.0.AND. NCAC(I).NE.NCAC(1)) GO TO 602
      NURC=NURC+1
      EL(NURC)=ELC(I)
      JO(NURC)=JOC(I)
      NTU(NURC)=NTUC(I)
      NNB(NURC)=NNBC(I)
      NNG(NURC)=NNGC(I)
      NNO(NURC)=NNOC(I)      
      NPO(NURC)=NPOC(I)
      NCA(NURC)=NCAC(I)
  602 CONTINUE
      NUR=NURC

  639 CONTINUE 
  
       IF(MEHAM.GT.2) CALL PREQU 
  
       ANEU=1.008664924
       IF(MECHA.EQ.1) ANEU=1.007825032
       AMI=939.56536D0
       IF(MECHA.EQ.1) AMI=938.272029D0
       
       IF(MECHA.EQ.1.AND.IDINT(ASP+0.1D0).EQ.1) AMI=1875.612859D0 
       IF(MECHA.EQ.1.AND.IDINT(ASP+0.1D0).EQ.1) ANEU=2.013553212712D0 
       
      REL=(EN+AMI)/AMI
      IF(MEREL.EQ.0) REL=1. 
      ENC=EN*AT/(AT+ANEU*REL)
      DO 5 I1=1,NUR
      IF(ENC-EL(I1))6,6,5
    5 CONTINUE
      NMAX=NUR
      GO TO 7
    6 NMAX=I1-1
    7 CONTINUE
C     KODMA=KOD
      IF(NMAX.LT.NUR) KODMA=0
   
      CALL RIPAT
      CALL ASFUT
      
      IF(MEPOT.GT.1) CALL KNCOE
        
      CALL QUANT
              
      IF(MTET.EQ.0) GO TO 15
      NUI=1
      NUF=NMAX
      
      KEYAP=1
      CALL ANPOW                       !NOT WORKING YET !!!!!
      IF(MEPRI.NE.98) THEN
      IF(MEPRI.LT.98) PRINT 300, PNAME
      WRITE(21,300) PNAME
  300 FORMAT(/23X,'ANALYZING POWERS FOR SCATTERED ',A8/)
      WRITE(327,'(F10.6,2I3)') EN,MTET,NMAX
      DO 314 M=1,MTET
      IF(MEPRI.LT.98) PRINT 11,TET(M),(DISC(K,M),K=1,NMAX)
      WRITE(21,11)TET(M),(DISC(K,M),K=1,NMAX)
      WRITE(327,111)TET(M),(DISC(K,M),K=1,NMAX)
  314 CONTINUE  
      ENDIF
    
      KEYAP=2
      CALL ANPOW                          !NOT WORKING YET !!!!!

      IF(MEPRI.NE.98) THEN
      IF(MEPRI.LT.98) PRINT 338, PNAME
      WRITE(21,338) PNAME
  338 FORMAT(/29X,'POLARIZATION FOR SCATTERED ',A8/)
      WRITE(328,'(F10.6,2I3)') EN,MTET,NMAX
      DO 315 M=1,MTET
      IF(MEPRI.LT.98) PRINT 11,TET(M),(DISC(K,M),K=1,NMAX)
      WRITE(21,11)TET(M),(DISC(K,M),K=1,NMAX)
      WRITE(328,111)TET(M),(DISC(K,M),K=1,NMAX)      
  315 CONTINUE
      ENDIF

      CALL DISCA
      
      IF(MEPRI.NE.98) THEN
      IF(MEPRI.LT.98) PRINT 110, PNAME
      WRITE(21,110) PNAME
  110 FORMAT(/23X,'ANGULAR DISTRIBUTIONS OF SCATTERED ',A8/
     *'    ANGLES[deg]      gs          1st       ...      '  )
      WRITE(27,'(F10.6,2I3)') EN,MTET,NMAX
      DO 14 M=1,MTET
      IF(MEPRI.LT.98) PRINT 11,TET(M),(DISC(K,M),K=1,NMAX)
      WRITE(21,11)TET(M),(DISC(K,M),K=1,NMAX)
      WRITE(27,111)TET(M),(DISC(K,M),K=1,NMAX)
  111 FORMAT(E13.6,(8E13.5),(12E13.5))
   14 CONTINUE
   11 FORMAT(2X,F10.3,3X,(8E11.3),(12E11.3))
      ENDIF

      DO 25 L=1,40
      IF(COEFR(L).EQ.0.AND.COEFI(L).EQ.0.) LCOUL=L-1
      IF(COEFR(L).EQ.0.AND.COEFI(L).EQ.0.) GO TO 26
   25 CONTINUE
   26 IF(MEPRI.NE.98) THEN
        IF(MEPRI.LT.98) PRINT 101, EN,LKK,LCOUL
        WRITE(26,101)EN,LKK,LCOUL
        WRITE(26,131)EN,CST,CSR,(CSN(K),K=1,NMAX)
  101 FORMAT (/3X,'EN=',F9.3,3X,'LKK=',I3,3X,'LCOUL=',I3)
        DO 9 N=1,NMAX
          IF(MEPRI.LT.98) PRINT 112,(COEF(N,L),L=1,LKK)
          WRITE(26,112)(COEF(N,L),L=1,LKK)
          IF(ETA.EQ.0.D0.OR.N.GT.1) GO TO 9
          IF(MEPRI.LT.98) PRINT 113, (COEFR(L),COEFI(L),L=1,LCOUL)
          WRITE(26,113) (COEFR(L),COEFI(L),L=1,LCOUL)
    9     CONTINUE
 112  FORMAT(/14X,'LEGENDRE COEFFICIENTS FOR  SCATTERED NUCLEONS'/
     *16X,'ANGULAR DISTRIBUTIONS - NUCLEAR AMPLITUDE'/(1X,6E15.7))
 113  FORMAT(/14X,'LEGENDRE COEFFICIENTS FOR  SCATTERED PROTONS'/
     *16X,'ANGULAR DISTRIBUTIONS - COULOMB AMPLITUDE'/(1X,6E15.7))
      ENDIF 
   15 IF(MECHA.NE.0) GO TO 102
C RCN
      IF(MEPRI.LT.98) PRINT 1210,EN,CST
      IF(EN.LT.2.d0) THEN
        IF(MEPRI.LT.98) PRINT 12,EN,CST,CSN(1),CSR,CST-CSN(1),
     *    SQRT((CST-CSR)/0.125663706D0)
        WRITE(21,12)EN,CST,CSN(1),CSR,CST-CSN(1),
     *    SQRT((CST-CSR)/0.125663706D0)
      ELSE
        IF(MEPRI.LT.98) PRINT 120,EN,CST,CSN(1),CSR,CST-CSN(1)
        WRITE(21,120)EN,CST,CSN(1),CSR,CST-CSN(1)
      ENDIF
C RCN
      IF(MEPRI.NE.98) WRITE(25,33)EN
      GO TO 103

  102 IF(MEPRI.NE.98) THEN
         IF(MEPRI.LT.98) PRINT 104,EN,CST,CSN(1),CSR,CST-CSN(1)
         IF(MEPRI.LT.98) PRINT 1214,EN,CST
         WRITE(21,104)EN,CST,CSN(1),CSR,CST-CSN(1)
      ENDIF

      IF(MEPRI.NE.98) WRITE(25,34)EN
  103 IF(MEPRI.NE.98) WRITE(23,131)EN,CST,CSN(1),CST-CSN(1),
     *                       CSR,(CSN(K),K=2,NURRR)
  131 FORMAT(1P50E14.5)
      IF(MEPRI.NE.98) THEN 
      DO 31 M=1,MTET
      IF(MEPRI.LT.98) PRINT 11,TET(M),(DISC(K,M),K=1,NMAX)
      WRITE(25,11)TET(M),(DISC(K,M),K=1,NMAX)
   31 CONTINUE
      ENDIF
   33 FORMAT(///1X,'NEUTRON ENERGY =',F10.6)
   34 FORMAT(///1X,'PROTON  ENERGY =',F10.6)
C     IF(MEPRI.NE.99) PRINT 130,(K,CSN(K),K=1,NMAX)
C     WRITE(21,130)(K,CSN(K),K=1,NMAX)
      IF(MEPRI.LT.98) PRINT 130,(K,EL(k),0.5*JO(k),
     *                cpar(NPO(k)),CSN(K),K=1,NMAX)
      WRITE(21,130)(K,EL(k),0.5*JO(k),cpar(NPO(k)),CSN(K),K=1,NMAX)
      IF(MEPRI.LT.98) PRINT 129,SF0,SF1,SF2
      WRITE(21,129)SF0,SF1,SF2
c
c [GN] 11/2015 add strenght function from SPRT+ESW
c
      WRITE(21,229)SFR0,SFR1,SFR2
      WRITE(21,329)RRPRIME0,RRPRIME1,RRPRIME2
c [GN] end      
C RCN
C
C     CROSS SECTION FILES
C

      IF(MEPRI.EQ.98) THEN  
        open(unit=93,file=TRIM(fname)//'.CS')
        open(unit=98,file=TRIM(fname)//'.ICS')

       IF(ETA.EQ.0.D0) THEN
C        WRITE(93,'(10H<CROSS-S.>,F10.2,F10.5,F10.2,2I5)') 
         WRITE(93,1009) ANEU,EN,AT,NINT(0.5*JO(1)),3
       ELSE
C        WRITE(93,'(10H<CROSS-S.>,F10.2,F10.5,F10.2,2I5)') 
         WRITE(93,1009) ANEU,EN,AT,NINT(0.5*JO(1)),1
       ENDIF
 1009  FORMAT ('<CROSS-S.>',F10.2,1P,D20.8,0P,F10.2,2I5)                 RESU-604
 1010  FORMAT ('<INE.C.S.>',F10.2,1P,D20.8,0P,F10.2,2I5)                 RESU-605
C      WRITE(98,'(10H<INE.C.S.>,F10.2,F10.5,F10.2,2I5)') 
       WRITE(98,1010) ANEU,EN,AT,NINT(0.5*JO(1)),NMAX-1

       IF(ETA.EQ.0.D0) THEN
C
C        TOTAL
C        WRITE(93,'(1X,E14.8)') CST*1000.
         WRITE(93,1012) CST*1000.d0
 1012    FORMAT (1P,D12.5)                                               RESU-607
C
C        INELASTIC TO LEVELS
         DO K=2,NMAX
C          WRITE(98,'(1X,E14.8)') CSN(K)*1000.
           WRITE(98,1012) CSN(K)*1000.d0
         ENDDO
C
C        REACTION + INELASTIC TO LEVELS
C        WRITE(93,'(1X,E14.8)') (CST - CSN(1))*1000.
         WRITE(93,1012) (CST - CSN(1))*1000.d0
C
C        ELASTIC
C        WRITE(93,'(1X,E14.8)') CSN(1)*1000.
         WRITE(93,1012) CSN(1)*1000.d0

       ELSE
C
C        INELASTIC TO LEVELS
         DO K=2,NMAX
C          WRITE(98,'(1X,E14.8)') CSN(K)*1000.
           WRITE(98,1012) CSN(K)*1000.d0
         ENDDO
C
C        REACTION + INELASTIC TO LEVELS
C        WRITE(93,'(1X,E14.8)') (CST - CSN(1))*1000.
         WRITE(93,1012) (CST - CSN(1))*1000.d0

       ENDIF
       close(93)
       close(98)

      ENDIF     
       
      IF(MEPRI.EQ.98) THEN  
C
        open(unit=96,file=TRIM(fname)//'.LEG')
        open(unit=97,file=TRIM(fname)//'.ANG')
C
c       IF(ETA.EQ.0.) WRITE(96,'(10H<LEGENDRE>,F10.2,F10.5,F10.2,2I5)') 
        IF(ETA.EQ.0.D0) WRITE(96,1000) ANEU,EN,AT,NINT(0.5*JO(1)),NMAX
 1000   FORMAT ('<LEGENDRE>',F10.2,1P,D20.8,0P,F10.2,2I5)               
c       WRITE(97,'(10H<ANG.DIS.>,F10.2,F10.5,F10.2,2I5)')
        WRITE(97,1008) ANEU,EN,AT,NINT(0.5*JO(1)),NMAX
 1008   FORMAT ('<ANG.DIS.>',F10.2,1P,D20.8,0P,F10.2,2I5)               
        DO K=1,NMAX
          IF(ETA.EQ.0.D0) 
     *     WRITE(96,'(2I5,'' COUPLED LEVEL, NUMBER OF VALUES'')') K, LKK
          IF(NPO(K).eq.+1)
     *      WRITE(97,'(I5,F5.1,A1,I4,I5)') K,0.5*JO(K),'+',1,MTET
C                                                          Cross section printed
C    *      WRITE(97,'(I5,F5.1,A1,I4,I5)') K,0.5*JO(K),'+',0,MTET  
          IF(NPO(K).eq.-1)
     *      WRITE(97,'(I5,F5.1,A1,I4,I5)') K,0.5*JO(K),'-',1,MTET
C                                                          Cross section printed
C    *      WRITE(97,'(I5,F5.1,A1,I4,I5)') K,0.5*JO(K),'-',0,MTET  
          IF(ETA.EQ.0.D0) THEN
            DO L=1,LKK
C             WRITE(96,'(2I5,1P,D20.10)') K,L-1,COEF(K,L)
C                                               OUTPUT IN MB 
              WRITE(96,'(2I5,1P,D20.10)') K,L-1,1000.D0*COEF(K,L) 
            ENDDO
          ENDIF
          DO M=1,MTET
            WRITE(97,1038) 0,TET(M),DISC(K,M)*1000.d0,'CROSS-SECTION   '
C           WRITE(97,1038)   TET(M),DISC(K,M)*1000.d0,'CROSS-SECTION   '
C1038       FORMAT (I3,1P,2D12.5,5X,4A4,A2)                             RESU-648
 1038       FORMAT (I3,1P,2D12.5,5X,A16)     
          ENDDO
        ENDDO
        close(96)
        close(97)
      ENDIF
 1210 FORMAT( 1X,'NEUTRON ENERGY =',F10.6,2X,'TOTAL CR-SECT.=',F10.6)
   12 FORMAT(/1X,'NEUTRON ENERGY =',F10.6/1X,'TOTAL  CR-SECT.=',F10.6/
     *1X,'ELASTIC  CR-SECT. =',F10.6/     
     *1X,'REACTION CR-SECT. =',F10.6/
     *1X,'REACTION CR-SECT. incl. coupled levels =',F10.6/
     *1X,'SCATTERING RADIUS =',F10.6)
  120 FORMAT(/1X,'NEUTRON ENERGY =',F10.6/1X,'TOTAL  CR-SECT.=',F10.6/
     *1X,'ELASTIC  CR-SECT. =',F10.6/     
     *1X,'REACTION CR-SECT. =',F10.6/
     *1X,'REACTION CR-SECT. incl. coupled levels =',F10.6)
C RCN
 1214 FORMAT( 1X,'PROTON  ENERGY =',F10.6,2X,'TOTAL CR-SECT.=',F10.6)
  104 FORMAT(/1X,'PROTON  ENERGY =',F10.6/1X,'TOTAL  CR-SECT.=',F10.6/
     *1X,'ELASTIC  CR-SECT. =',F10.6/ 
     *1X,'REACTION CR-SECT. =',F10.6/
     *1X,'REACTION CR-SECT. incl. coupled levels =',F10.6)
C 130 FORMAT(/3X,'Nlev',17X,'CR-SECT. OF LEVEL EXCITATION '
C    */(1X,I5,25X,F10.6))
  130 FORMAT(
     */2x,'Nlev',4X,'Elev',3x,'Jpi',9x,'CR-SECT(Nlev)'
     */(2X,I2,3X,F7.4,2x,F4.1,A1,10X,F10.6))
  129 FORMAT(/30X,'STRENGTH  FUNCTIONS'
     */1X,'SF0=',E15.7,8X,'SF1=',E15.7,8X,'SF2=',E15.7/)
c
c [GN] 11/2015 add strenght function from SPRT+ESW
c
  229 FORMAT(/25X,'STRENGTH  FUNCTIONS FROM ESW'
     */1X,'S0  =',E15.7,8X,'S1  =',E15.7,8X,'S2  =',E15.7)
  329 FORMAT(1X,'R0  =',E15.7,8X,'R1  =',E15.7,8X,'R2  =',E15.7)
c [GN] end

      NUR=NURRR
      TID = 0
!$    TID = OMP_GET_THREAD_NUM()
!$    IF(MEPRI.NE.98) 
!$   &  PRINT *, 'Thread',TID,' finished.',' IIparal=',IIparal

      RETURN
      END

