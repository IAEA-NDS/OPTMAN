      SUBROUTINE DATETpar(IE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INTEGER IE
      
      INTEGER NNTTii,MEISii,IIS,IIIS
      REAL*8 FUii,FUU
      CHARACTER*1 cpar
      DIMENSION JTEMP(40), JSHIFT(40)
C---------------------------------
C     These commons are also used and initialized in ABCT 
      INCLUDE 'PRIVCOM10.FOR'
      INCLUDE 'PRIVCOM12.FOR'
      INCLUDE 'PRIVCOM9.FOR'
      INCLUDE 'PRIVCOM13.FOR'
      INCLUDE 'PRIVCOM1.FOR'
      INCLUDE 'PRIVCOM.FOR'         
      INCLUDE 'PRIVCOM15.FOR'
      INCLUDE 'PRIVCOM6.FOR'
      INCLUDE 'PRIVCOM20.FOR'
C--------------------------------- 
C
      INCLUDE 'PRIVCOM8.FOR'

      INCLUDE 'PRIVCOM16D.FOR'  
      INCLUDE 'PRIVCOM17D.FOR'
C---------------------------------
C     These common is used FOR initialization CCOULii <-> CCOUL
      INCLUDE 'PRIVCOM18D.FOR'
  

      INCLUDE 'PRIVCOM22.FOR'
      CHARACTER*8 PNAME


      INTEGER TID
!$    INTEGER OMP_GET_THREAD_NUM

      TID = 0
!$    TID = OMP_GET_THREAD_NUM()


 
      
      FUii   = 0.d0
      NNTTii = 0
      MEISii = 1

      NNRA=0
      WCST=0.D0
      SUMWEI=0.D0      
      IIIS=IIS  
      CCOULi=CCCOUL  
C     write(*,*) 'Nuc.Index=',IIS,' AT=',NINT(ATIS(IIS))
 740  EFERMN=EFISN(IIIS)
      EFERMP=EFISP(IIIS)
      ZNUC=ZNUCIS(IIIS)      
      AT=ATIS(IIIS)
      ASQ=AT**(1.d0/3.d0)
      NRESN=NRES(IIIS)
      WEIGHT=WEI(IIIS)
      MESOL=MESOIS(IIIS)
      MEDEF=MEDEIS(IIIS)
      MEAXI=MEAXIS(IIIS)

C---- Only two quantitites that depend on IIS (not IIIS) 
      EN=EEIS(IIS,IE)
      MECHA=MCHAIS(IIS,IE)
      
      PNAME='NEUTRONS'
      IF(MECHA.NE.0) PNAME=' PROTONS'      
C----
      WRITE (21,118) NINT(AT),NINT(ZNUC),EN
  118 FORMAT(//1X,'ADJUSTING TO EXPERIMENTAL DATA FOR NUCLEUS
     * WITH MASS NUMBER=',I3,' AND CHARGE=',I3,' EN=',D12.6/)

C     PRINT 126, 'Thread ',TID,'start IE=',IE,' En=',EN,' AT=',INT(AT),
C    *           ' IIS=',IIS,' IIIS=',IIIS,' MECHA=',MECHA,
C    *           ' MEIS=',MEISii
C126  FORMAT(/1x,A7,I2,A10,I2,A4,F8.4,A4,I3,A5,I2,A6,I2,
C    *          A7,I2,A6,I1)

      IF(MEPOT.EQ.1.AND.MEDEF.EQ.0.AND.MEAXI.EQ.0.AND.MEVOL.EQ.0)
     *     GO TO 701
      HW=HWIS(IIIS)
      AMB0=AMB0IS(IIIS)
      AMG0=AMG0IS(IIIS)
      GAM0=GAM0IS(IIIS)
      BET0=BET0IS(IIIS)
      BET4=BET4IS(IIIS)
      BB42=BB42IS(IIIS)
      GAMG=GAMGIS(IIIS)
      DELG=DELGIS(IIIS)
      BET3=BET3IS(IIIS)
      ETO=ETOIS(IIIS)
      AMUO=AMUOIS(IIIS)
      HWO=HWOIS(IIIS)
      BB32=BB32IS(IIIS)
      GAMDE=GAMDIS(IIIS)
      DPAR=DPARIS(IIIS)
      GSHAPE=GSHAEIS(IIIS)
C     GO TO 702
      BET2SUM=0.d0
  701 DO I=2,NPD,2
         BET(I)=BETIS(IIIS,I)
         BET2SUM=BET2SUM+BET(I)**2
      END DO
      !BET2SUM=BET2SUM+BET3**2
      RCORR=1.d0
      IF(MERAD.EQ.1) RCORR=1.d0-BET2SUM*7.9577471546d-2 ! 1-bet2sum/(4*pi)
C      BETB(MELEV)=BETBIS(IIIS,MELEV)
C
  702 CONTINUE
C
C     dtmp = DBLE(NINT(ATIS(IIIS)-ATIS(1)))
      dtmp = ATIS(IIIS)-ATIS(1)
      VRLA=VRG+CAVR*dtmp
      RR=(RRG*RCORR+CARR*dtmp)*ASQ
      RC=RCG*RCORR*ASQ
      RD=(RDG*RCORR+CARD*dtmp)*ASQ
      !RD=RDG*RCORR*ASQ
      RW=RWG*RCORR*ASQ
      RS=RSG*RCORR*ASQ
      RZ=RZG*RCORR*ASQ 
      AR0=ARG+CAAR*dtmp
      AC0=ACG+CAAC*dtmp
      !AC0=ACG
      
      WDBW=WDBWG+CAWD*dtmp
      WDWID=WDWIDG+CAWDW*dtmp
       
C     PRINT 131, 'Thread ',TID,' VRG=',VRG,' RRG=',RRG,
C    * ' RDG=',RDG,' ARG=',ARG,' ACG=',ACG,'  IIS=',IIS,' IIIS=',IIIS 
C     WRITE(21,131) 'Thread ',TID,' VRG=',VRG,' RRG=',RRG,
C    * ' RDG=',RDG,' ARG=',ARG,' ACG=',ACG
C131  FORMAT(1x,A7,I2,5(A5,d12.6),2(A6,I2))
C     PRINT 132, 'Thread ',TID,' ATIS(1)=',INT(ATIS(1)),' AT-ATIS(1)=',
C    *  INT(AT-ATIS(1)),' VRLA=',VRLA,' RR=',RR,' RC=', RC,' RD=', RD,
C    *  ' RW=', RW,' RS=', RS,' RZ=', RZ,' AR0=',AR0,' AC0=',AC0
C132  FORMAT(1x,A7,I2,A9,I2,A12,I2,A6,d12.6,3(A4,d12.6)/
C    *   10x,3(A4,d12.6),2(A5,d12.6))
C
      NUR=NURIS(IIIS)
      NURRR=NUR
      NURC=0

       IF(MEPOT.GT.1) GO TO 638
      
      JSHIFT=0 
       
      DO 601 I=1, NUR
      IF(MECHA.EQ.0.AND. NCAIS(IIIS,I).NE.NCAIS(IIIS,1)) GO TO 601
      NURC=NURC+1
      EL(NURC)=ELIS(IIIS,I)
      JO(NURC)=JOIS(IIIS,I)
      NPO(NURC)=NPOIS(IIIS,I)
      KO(NURC)=KOIS(IIIS,I)
      NCA(NURC)=NCAIS(IIIS,I)
      NUMB(NURC)=NUMBIS(IIIS,I)
      BETB(NURC)=BETBIS(IIIS,I)
      NUMB(NURC)=NUMBIS(IIIS,I)
      BETB(NURC)=BETBIS(IIIS,I)
      AIGS(NURC)=AIGSIS(IIIS,I)
      NTU(NURC)=NTUIS(IIIS,I) 
      NNB(NURC)=NNBIS(IIIS,I)
      NNG(NURC)=NNGIS(IIIS,I)
      NNO(NURC)=NNOIS(IIIS,I)
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
       IF(MEDEF.GT.0.OR.MEAXI.EQ.1.OR.MEVOL.GE.1) CALL OVLOPT 
       DO IID=1,NUR
         DO JJD=IID,NUR
             EFFDEF(JJD,IID,:)=EFFDEF(IID,JJD,:)
             
         END DO
       END DO

       
       IF(MOD(JO(1),2).GT.0)  THEN
c          NUMBGS=NUMB(1)
C            DO IID=1,NUR
C              DO JJD=1,NUR
C                IF(NUMB(IID).NE.NUMBGS.OR.NUMB(JJD).NE.NUMBGS)
C      *                  EFFDEF(JJD,IID,:)=0.0
C              END DO
C            END DO          
           JU=JTEMP
       END IF
       
      DEFNUL=0.D0
      DEFNUL=SUM(EFFDEF*EFFDEF)     

      GO TO 639
  638 DO 602 I=1, NUR
      IF(MECHA.EQ.0.AND. NCAIS(IIIS,I).NE.NCAIS(IIIS,1)) GO TO 602
      NURC=NURC+1
      EL(NURC)=ELIS(IIIS,I)
      JO(NURC)=JOIS(IIIS,I)
      NTU(NURC)=NTUIS(IIIS,I)
      NNB(NURC)=NNBIS(IIIS,I)
      NNG(NURC)=NNGIS(IIIS,I)
      NNO(NURC)=NNOIS(IIIS,I)      
      NPO(NURC)=NPOIS(IIIS,I)
      NCA(NURC)=NCAIS(IIIS,I)
  602 CONTINUE
      NUR=NURC
      !IF(MEISii.EQ.1 .AND.MEDEF.NE.0) CALL PREQU
      IF(MEHAM.GT.2) CALL PREQU 
      MEISii=0

  639 CONTINUE 
      IF(NRESN.EQ.0) GO TO 604
      DO 603 I=1,NRESN
      ERN(I)=ERIS(IIIS,I)
      GNN(I)=GNIS(IIIS,I)
      GREN(I)=GREIS(IIIS,I)
      LON(I)=LOIS(IIIS,I)
      JMN(I)=JMIS(IIIS,I)
      JCON(I)=JCOIS(IIIS,I)
  603 NEL(I)=NELA(IIIS,I)
C     CREATING LEVELS FOR (P,N) ANALOG STATES CALCULATIONS
  604 CONTINUE
C      IF(EEIS(IIIS,IE).LT.EL(NURRR)*(AT+1.007825032D0)/AT+0.5D0)
C     *GO TO 779
C      IF(MEHAM.GE.1.AND.MCHAIS(IIIS,IE).EQ.1) NUR=NURRR
C      IF(MEHAM.GE.1.AND.MCHAIS(IIIS,IE).EQ.1) GO TO 777
C  779 DO 778 ILEV=1,NURRR
C      IF(NCAIS(IIIS,ILEV).EQ.NCAIS(IIIS,1)) NUR=ILEV
C 778 CONTINUE
C 777 CONTINUE
 
      ANEU=1.008664924D0
      IF(MECHA.EQ.1) ANEU=1.007825032D0
      AMI=939.56536D0
      IF(MECHA.EQ.1) AMI=938.272029D0
      
      IF(MECHA.EQ.1.AND.IDINT(ASP+0.1D0).EQ.1) AMI=1875.612859D0 
      IF(MECHA.EQ.1.AND.IDINT(ASP+0.1D0).EQ.1) ANEU=2.013553212712D0 
      
      REL=(EN+AMI)/AMI
      IF(MEREL.EQ.0) REL=1.d0
      ENC=EN*AT/(AT+ANEU*REL)
      DO 29 I1=1,NUR
      IF(ENC-EL(I1)) 30,30,29
   29 CONTINUE
      NMAX=NUR
      GO TO 31
   30 NMAX=I1-1
   31 CONTINUE
      IF(NMAX.LT.NUR) KODMA=0
C
C
C
C     IF(IIIS.EQ.IIS) THEN
C     WRITE(21,90) TID,IE,EN
C     PRINT    90, TID,IE,EN
C 90  FORMAT(/10X,'POTENTIAL PARAMETERS V(R), TID= ',I3,
C    *                  ' IE=',I3,' E=',D12.6/)
C     PRINT    91, VR0,VR1,VR2,RR,AR0,AR1,WD0,WD1,VR3,RD,AD0,
C     WRITE(21,91) VR0,VR1,VR2,RR,AR0,AR1,WD0,WD1,VR3,RD,AD0,
C    *AD1,WC0,WC1,RC,AC0,AC1,RW,AW0,AW1,VS,RS,AS0,AS1,ALF,ANEU,RZ,
C    *AZ,BNDC,WDA1,WCA1,CCOUL,CISO,WCISO,WS0,WS1,VRLA,ALAVR,
C    *WCBW,WCWID,WDBW,WDWID,ALAWD,EFERMN,EFERMP,ALASO,PDIS,
C    *WSBW,WSWID,RRBWC,RRWID,RZBWC,RZWID,EA,WDISO,WDSHI,WDWID2,
C    *ALFNEW,VRD,CAVR,CARR,CAAR,CARD,CAAC,AT,ASP,ENC,REL,KODMA,NMAX,
C    *MECHA,LAS,NPD,NUR
C     ENDIF
   91 FORMAT(/1x,'VR0=',D13.7,5X,'VR1=',D13.7,5X,
     *'VR2=',D13.7,2X,'RR=',D13.7,5X,'AR0=',D13.7,5X,'AR1=',D13.7
     */1X,'WD0=',D13.7,5X,'WD1=',D13.7,5X,'VR3=',D13.7,2X,
     *'RD=',D13.7,5X,'AD0=',D13.7,5X,'AD1=',D13.7
     */1X,'WC0=',D13.7,5X,'WC1=',D13.7,24X,
     *'RC=',D13.7,5X,'AC0=',D13.7,5X,'AC1=',D13.7
     */45X,'RW=',D13.7,3X,'AW0=',D13.7,4X,'AW1=',D13.7
     */1X,'VSO=',D13.7,5X,'RS=',D13.7,6X,'AS0=',D13.7,2X,'AS1=',D13.7
     */1X,'ALF=',D13.7,5X,'ANEU=',D13.7,23X,'RZ=',D13.7,5X,'AZ0=',D13.7,
     */1X,'BNDC=',D13.7,4X,'WDA1=',D13.7,4X,'WCA1=',D13.7,4X,
     *'CCOUL=',D13.7,5X,'CISO=',D13.7,4X,'WCISO=',D13.7
     */1X,'WS0=',D13.7,5X,'WS1=',D13.7,5X,'VRLA=',D13.7
     *,4X,'ALAVR=',D13.7,5X,'WCBW=',D13.7,4X,'WCWID=',D13.7,/1X,'WDBW='
     *,D13.7,4X,'WDWID=',D13.7,3X,'ALAWD=',D13.7
     *,3X,'EFERMN=',D13.7,4X,'EFERMP=',D13.7,2X,'ALASO=',D13.7,
     */1X,'PDIS=',D13.7,4X,'WSBW=',D13.7,4X,'WSWID=',D13.7,3X,
     *'RRBWC=',D13.7,5X,'RRWID=',D13.7,3X,'RZBWC=',D13.7,
     */1X,'RZWID=',D13.7,3X,'EA=',D13.7,6X,'WDISO=',D13.7,
     *3X,'WDSHI=',D13.7,5X,'WDWID2=',D13.7,2X,'ALFNEW=',D13.7,
     */1X,'VRD=',D13.7,5X,'CAVR=',D13.7,4X,'CARR=',D13.7,
     *4X,'CAAR=',D13.7,6X,'CARD=',D13.7,4X,'CAAC=',D13.7,
     */1X,'AT=',D13.7,6X,'ASP=',D13.7,5X,'ENC=',D13.7,5X,'REL=',D13.7,
     */1X,'KODMA=',I2,'  NMAX=',I3,'  MECHA=',I2,'  LAS=',I2,' NPD=',I2,
     *    '  NUR=',I3 /)
C
C
C
      CCOUL=CCOULii
      CALL RIPAT
      CALL ASFUT
C      PRINT*,CCOULii,CCOUL,3
 
      IF(MEPOT.GT.1) CALL KNCOE
      CALL QUANT
C      Write (21,99999)EFERMN,EFERMP,ZNUC,AT,ASQ,EN
C99999 format(/10e11.4)
C      Write (21,99999)VRLA,RR,RC,RD,RW,RS,RZ,AZ,AR0,AC0
C      Write (21,99999)VRG,CAVR,RRG,CARR,RCG,RDG,CARD,RWG,RSG,RZG,ARG,
C     *CAAR,ACG,CAAC
C      Write (21,99997)CDE,CCOUL,MECHA,MECUL,VISO,CISO,WISO,WCISO,WVISO,
C     *WDISO
C99997 format(/2e11.4,2I3,6e11.4)      
C      Write (21,99998)EL(18),NUR,JO(18),NPO(18),KO(18),NCA(18),
C     *NUMB(18),BETB(18),NUMB(18),BETB(18),AIGS(18)
C      Write (21,99998)EL(20),NUR,JO(20),NPO(20),KO(20),NCA(20),
C     *NUMB(20),BETB(20),NUMB(20),BETB(20),AIGS(20)
     
C99998 format(/e11.4,6i3,e11.4,i3,2e11.4) 
C      Write (21,99999) BET(2),BET(4),BET(6)
C      Write (21,99999)(BETB(KK),AIGS(KK),KK=1,NUR)
C
      IF(NNRA.EQ.1.) GO TO 788 ! IF RATIO goto 788
C
      IF(MECHA.NE.0) GO TO 102
      IF(EN.GT.2.D0) THEN
        IF(MEPRI.LT.98) PRINT 112,EN,CST,CSR,(CST-CSR)
        WRITE(21,112)EN,CST,CSR,(CST-CSR)
      ELSE
        IF(MEPRI.LT.98) PRINT 92,EN,CST,CSR,(CST-CSR),
     *                  SQRT((CST-CSR)/0.125663706D0)
        WRITE(21,92)EN,CST,CSR,(CST-CSR),SQRT((CST-CSR)/0.125663706D0)
      ENDIF
C     PRINT 1292,EN,CST,TID
      GO TO 103
  102 IF(MEPRI.LT.98) PRINT 104,EN,CST,CSR
C     PRINT 1214,EN,CST,TID
      WRITE(21,104)EN,CST,CSR
  103 IF(MEPRI.LT.98) PRINT 130,   
     *    (K,EL(k),0.5*JO(k),cpar(NPO(k)),CSN(K),K=1,NMAX)
      WRITE(21,130)(K,EL(k),0.5*JO(k),cpar(NPO(k)),CSN(K),K=1,NMAX)

C 103 IF(MEPRI.LT.98) PRINT 130,(K,CSN(K),K=1,NMAX)
C     WRITE(21,130)(K,CSN(K),K=1,NMAX)
      IF(EN.GT.0.75) GOTO 33
      IF(MEPRI.LT.98) PRINT 129,SF0,SF1,SF2
      WRITE(21,129)SF0,SF1,SF2
c
c [GN] 11/2015 add strenght function from sprt+ESW 
c
      WRITE(21,229)SFR0,SFR1,SFR2
      WRITE(21,329)RRPRIME0,RRPRIME1,RRPRIME2
c [GN] end     
 1292 FORMAT( 1X,'NEUTRON ENERGY =',F10.6,2X,'TOTAL CR-SECT.=',F10.6,
     * ' TID=',I2)
   92 FORMAT(/1X,'NEUTRON ENERGY =',F10.6/1X,'TOTAL  CR-SECT.=',F10.6/
     *1X,'REACTION CR-SECT. =',F10.6/
     *1X,'TOTAL DIRECT CR-SECT.(ELASTIC + DIR.LEV EXCIT.) =',F10.6/
     *1X,'SCATTERING RADIUS =',F10.6)
  112 FORMAT(/1X,'NEUTRON ENERGY =',F10.6/1X,'TOTAL  CR-SECT.=',F10.6/
     *1X,'REACTION CR-SECT. =',F10.6/
     *1X,'TOTAL DIRECT CR-SECT.(ELASTIC + DIR.LEV EXCIT.) =',F10.6)
 1214 FORMAT( 1X,'PROTON  ENERGY =',F10.6,2X,'TOTAL CR-SECT.=',F10.6,
     * ' TID=',I2)
  104 FORMAT(/1X,'PROTON  ENERGY =',F10.6/1X,'TOTAL  CR-SECT.=',F10.6/
     *1X,'REACTION CR-SECT. =',F10.6)
C 130 FORMAT(/3X,'NMAX',17X,'CR-SECT. OF LEVEL EXCITATION '
C    */(1X,I5,25X,F10.6))
  130 FORMAT(
     */2x,'Nlev',4X,'Elev',3x,'Jpi',9x,'CR-SECT(Nlev)'
     */(2X,I2,3X,D13.7,2x,F4.1,A1,10X,F10.6))
  129 FORMAT(/30X,'STRENGTH  FUNCTIONS'
     */1X,'SF0=',E15.7,8X,'SF1=',E15.7,8X,'SF2=',E15.7)
c
c [GN] 11/2015 add strenght function from SPRT+ESW
c
  229 FORMAT(/25X,'STRENGTH  FUNCTIONS FROM ESW'
     */1X,'S0  =',E15.7,8X,'S1  =',E15.7,8X,'S2  =',E15.7)
  329 FORMAT(1X,'R0  =',E15.7,8X,'R1  =',E15.7,8X,'R2  =',E15.7)
c [GN] end
      IF(NSF1(IIS,IE).NE.1) GO TO 32
      NNTTii=NNTTii+1
C     FU=FU+((SE1(IIS,IE)-SF0)/DS1(IIS,IE))**2
      FUU=((SE1(IIS,IE)-SF0)/DS1(IIS,IE))**2
      FUii=FUii+FUU
      WRITE (21,149)FUU,SF0,SE1(IIS,IE),DS1(IIS,IE)
  149 FORMAT(/1X,'FU FOR S0 STRENGTH FUNCTION=',E14.7, 
     *'    CALC :',E14.7,' EXP :',E14.7,' +/- ',E14.7/)      
   32 IF(NSF2(IIS,IE).NE.1) GO TO 33
      NNTTii=NNTTii+1
C     FU=FU+((SE2(IIS,IE)-SF1)/DS2(IIS,IE))**2
      FUU=((SE2(IIS,IE)-SF1)/DS2(IIS,IE))**2
      FUii=FUii+FUU
      WRITE (21,152)FUU,SF1,SE2(IIS,IE),DS2(IIS,IE)
  152 FORMAT(/1X,'FU FOR S1 STRENGTH FUNCTION=',E14.7, 
     *'    CALC :',E14.7,' EXP :',E14.7,' +/- ',E14.7/)   
   33 IF(NT(IIS,IE).NE.1) GO TO 2
      NNTTii=NNTTii+1
C     FU=FU+((STE(IIS,IE)-CST)/DST(IIS,IE))**2
      FUU=((STE(IIS,IE)-CST)/DST(IIS,IE))**2
      FUii=FUii+FUU
      WRITE (21,150)FUU,CST,STE(IIS,IE),DST(IIS,IE)
  150 FORMAT(/1X,'FU FOR TOTAL CS=',E14.7, 
     *'    CALC :',F10.4,' EXP :',F10.4,' +/- ',F9.4/)
    2 IF(NR(IIS,IE).NE.1) GO TO 3
      NNTTii=NNTTii+1
C     FU=FU+((SRE(IIS,IE)-CSR)/DSR(IIS,IE))**2
      FUU=((SRE(IIS,IE)-CSR)/DSR(IIS,IE))**2
      FUii=FUii+FUU
      WRITE (21,153)FUU,CSR,SRE(IIS,IE),DSR(IIS,IE)
  153 FORMAT(/1X,'FU FOR REACTION CS=',E14.7, 
     *'    CALC :',F10.4,' EXP :',F10.4,' +/- ',F9.4/)
    3 NG=NGN(IIS,IE)
      IF(NG.EQ.0) GO TO 6
      DO 4 KG=1,NG
      NNTTii=NNTTii+1
      NUI=NIN(IIS,IE,KG)
      NUF=NFN(IIS,IE,KG)
      SNG=0.
      DO 5 I=NUI,NUF
    5 SNG=SNG+CSN(I)
      FUU=((SNE(IIS,IE,KG)-SNG)/DSN(IIS,IE,KG))**2
      IF (NUI.EQ.NUF.AND.NUI.NE.1) GO TO 4
      WRITE (21,156)FUU,SNG,SNE(IIS,IE,KG),DSN(IIS,IE,KG)
  156 FORMAT(/1X,'FU FOR ELASTIC SCATTERING CS "LOW ENERGY" R =',E14.7, 
     *'    CALC :',F10.4,' EXP :',F10.4,' +/- ',F9.4/)
    4 FUii=FUii+FUU
    6 NG=NGD(IIS,IE)
      IF(NG.EQ.0) GO TO 1
      DO 7 KG=1,NG
      NNTTii=NNTTii+1
      NUI=NID(IIS,IE,KG)
      NUF=NFD(IIS,IE,KG)
      KEYAP=0
      IF(NUI.GT.40) KEYAP=1
      IF(NUI.GT.80) KEYAP=2 
      IF(KEYAP.EQ.1) NUI=NID(IIS,IE,KG)-40
      IF(KEYAP.EQ.1) NUF=NID(IIS,IE,KG)-40
      IF(KEYAP.EQ.2) NUI=NID(IIS,IE,KG)-80
      IF(KEYAP.EQ.2) NUF=NID(IIS,IE,KG)-80
      MTET=MTD(IIS,IE,KG)
      DO 8 M=1,MTET
    8 TET(M)=TED(IIS,IE,KG,M)
      IF(KEYAP.EQ.0)CALL DISCA
      IF(KEYAP.GT.0)CALL ANPOW
      DO 9 M=1,MTET
      DISG(M)=0.D0
      DO 9 I=NUI,NUF
    9 DISG(M)=DISG(M)+DISC(I,M)
      
      
      
      IF(KEYAP.EQ.0 .AND. MEPRI.LT.98) PRINT 100,
     *                        PNAME,TID,IIS,IE,EEIS(IIS,IE),KG,NUI,NUF
      IF(KEYAP.EQ.0) WRITE (21,100) 
     *                        PNAME,TID,IIS,IE,EEIS(IIS,IE),KG,NUI,NUF
  100 FORMAT(/23X,'ANGULAR DISTRIBUTIONS OF SCATTERED ',A8/
     *        19X,'THREAD ',I2,' IIS=',I2,' IE=',I2,' E=',F8.4,
     *            '  KG=',I2,' NUI=',I2,' NUF=',I2/)
      IF(KEYAP.EQ.0 .AND. MEPRI.LT.98) 
     *  PRINT 39,(M,TET(M),SNGD(IIS,IE,KG,M),DISG(M),M=1,MTET)
      IF(KEYAP.LE.0) WRITE(21,39)(M,TET(M),SNGD(IIS,IE,KG,M),DISG(M),
     *M=1,MTET)
   39 FORMAT(1X,2('MTET',2X,'ANGL(CENT)',1X,'EXP. C.-S.',1X,
     *'CALC. C.-S. ')/(1X,2(I3,3D12.5)))
      IF(KEYAP.EQ.1 .AND. MEPRI.LT.98) PRINT 300
      IF(KEYAP.EQ.1) WRITE (21,300)
  300 FORMAT(/23X,'ANALYZING POWERS FOR SCATTERED PARTICLES'/)
      IF(KEYAP.EQ.1 .AND. MEPRI.LT.98) 
     *  PRINT 339,(M,TET(M),SNGD(IIS,IE,KG,M),DISG(M),M=1,MTET)
      IF(KEYAP.EQ.1) WRITE(21,339)(M,TET(M),SNGD(IIS,IE,KG,M),DISG(M),
     *M=1,MTET)
  339 FORMAT(1X,2('MTET',2X,'ANGL(CENT)',1X,'EXP. AN.P.',1X,
     *'CALC. AN.P. ')/(1X,2(I3,3D12.5)))
     
      IF(KEYAP.EQ.2 .AND. MEPRI.LT.98) PRINT 328
      IF(KEYAP.EQ.2) WRITE (21,328)
  328 FORMAT(/20X,'POLARIZATION FOR SCATTERED PARTICLES'/)
      IF(KEYAP.EQ.2 .AND. MEPRI.LT.98) 
     *  PRINT 340,(M,TET(M),SNGD(IIS,IE,KG,M),DISG(M),M=1,MTET)
      IF(KEYAP.EQ.2) WRITE(21,340)(M,TET(M),SNGD(IIS,IE,KG,M),DISG(M),
     *M=1,MTET)
  340 FORMAT(1X,2('MTET',2X,'ANGL(CENT)',1X,'EXP. POLA.',1X,
     *'CALC. POLA. ')/(1X,2(I3,3D12.5)))

      dtmp=0.D0
      DO 10 M=1,MTET
   10 dtmp=dtmp+((SNGD(IIS,IE,KG,M)-DISG(M))/DSD(IIS,IE,KG,M))**2
      FUU=dtmp/MTET
      WRITE (21,151)FUU,MTET
  151 FORMAT(/1X,'FU FOR THIS ANGULAR DATA GROUP IS=',E14.7,
     *           ' # OF EXP.ANGLES=',I4/)
      FUii=FUii+FUU
    7 CONTINUE

    1 NUR=NURRR

      IF(NRAT(IIS,IE).EQ.0.AND.NNAT(IIS,IE).EQ.0) THEN
C===============================================================
C       TID = OMP_GET_THREAD_NUM()
        FU = FU + FUii
        NNTT = NNTT + NNTTii
        PRINT 125, 'Thread ',TID,' done. AT=',NINT(AT),' E=',EN, 
     *  ' IE=',IE,' IIS=',IIS,' ChisqR =',FUii/DBLE(NNTTii),
     *  ' ChisqT =',FU/DBLE(NNTT),' INT'
 125    FORMAT(1x,A7,I2,A10,I3,A3,D12.6,A4,I2,A5,I2,2(A9,D14.8),A4)
C===============================================================
        RETURN
      ENDIF

      IF(NRAT(IIS,IE).NE.0.AND.NNAT(IIS,IE).NE.0.AND.
     *NRAT(IIS,IE).NE.NNAT(IIS,IE)) WRITE(21,158)
  158 FORMAT(1X,'BOTH NRAT AND NNAT NON-EQUAL 0, THAN SHOULD BE EQUAL!')
      IF(NRAT(IIS,IE).NE.0.AND.NNAT(IIS,IE).NE.0.AND.
     *NRAT(IIS,IE).NE.NNAT(IIS,IE)) STOP
C
C     RATIO calculation
C      
  788 RATIOS=2.D0*(CSTR-CST)/(CSTR+CST) 
      IF(NNRA.EQ.1.AND.NRAT(IIS,IE).NE.0)WRITE (21,155)EN,CSTR, CST
C  155 FORMAT(/1X,'CS FOR RATIO AT THIS ENERGY ARE
C     *=',E14.7,2X,'AND',2X,E14.7/) 
     
  155 FORMAT(/1X,'CS FOR RATIO AT THE ENERGY=',E14.7,2X,'ARE
     *=',E14.7,2X,'AND',2X,E14.7/) 
     
      CSTR=CST
      WCST=WCST+WEIGHT*CSTR
      SUMWEI=SUMWEI+WEIGHT
            
      IF(NRAT(IIS,IE).NE.0) IIIS=NRAT(IIS,IE)
      IF(NNAT(IIS,IE).NE.0) IIIS=NNAT(IIS,IE)      
      NNRA=NNRA+1 
      MEISii=1
      IF(NNRA.LT.2) GO TO 740
       
      IF(NRAT(IIS,IE).EQ.0)GO TO 741
      NNTTii=NNTTii+1
      FUU=((RATIOS-RATIO(IIS,IE))/DRAT(IIS,IE))**2
      WRITE (21,154) FUU,RATIOS,RATIO(IIS,IE),
     *DRAT(IIS,IE)
      FUii=FUii+FUU
     
  741 IF(NNAT(IIS,IE).EQ.0)GO TO 742 
      NNTTi=NNTTii+1
      CALNAT= WCST/SUMWEI    
      FUU=((CALNAT-CSNAT(IIS,IE))/DCSNAT(IIS,IE))**2
      WRITE (21,157) FUU,CALNAT,CSNAT(IIS,IE),
     *DCSNAT(IIS,IE)
      FUii=FUii+FUU
            
  154 FORMAT(1X,'FU FOR RATIO IS=',E14.7,
     *'  CALC :',E14.7,' EXP :',E14.7,' +/- ',E14.7/) 
     
  157 FORMAT(1X,'FU FOR NATURAL TOTAL IS=',E14.7,
     *'  CALC :',E14.7,' EXP :',E14.7,' +/- ',E14.7/) 
     
 742  CONTINUE
C===============================================================
      FU = FU + FUii
      NNTT = NNTT + NNTTii
      PRINT 125, 'Thread ',TID,' done. AT=',NINT(AT),' E=',EN, 
     *  ' IE=',IE,' IIS=',IIS,' ChisqR =',FUii/DBLE(NNTTii),
     *  ' ChisqT =',FU/DBLE(NNTT),' RATIO'
C===============================================================
      RETURN
      END

