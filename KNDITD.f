C     *******************************************************
C     START of knditd
C     *******************************************************
C     ****************************************************************
      SUBROUTINE KNCOE
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM1.FOR'
      INCLUDE 'PRIVCOM9.FOR'
      INCLUDE 'PRIVCOM10.FOR' 
      INCLUDE 'PRIVCOM12.FOR'
      INCLUDE 'PRIVCOM13.FOR'    

      DIMENSION NLA(5)
      DATA NLA/0,1,3,6,10/

C     RCN June 2014
c     Following variables were being used in an uninitialized state
c     SDH, May 2014
      real*8 :: ak3 = 0.D0
      real*8 :: ak4 = 0.D0
      real*8 :: bk3 = 0.D0
      real*8 :: bk5 = 0.D0
      real*8 :: bk8 = 0.D0

      real*8, parameter :: SQ2 = 1.414213562373095D0
      real*8, parameter :: SQ32 = 2.D0*1.414213562373095D0
C     RCN June 2014

      C1=SQ2      
C     SQ32=2.828427D0
      C95=1.8D0
      CQ95=1.34164D0
      C75=7.D0/5.D0
      CSQ75=SQRT(C75)
      C6325=63.D0/25.D0
      C6325=SQRT(C6325)

      DO 11 I=1,15
      ALK11(I)=0.D0
      ALK21(I)=0.D0
      ALK22(I)=0.D0
      ALK31(I)=0.D0
      ALK32(I)=0.D0
      ALK41(I)=0.D0
      ALK42(I)=0.D0
      ALK43(I)=0.D0
      B0LK1(I)=0.D0
      B0LK2(I)=0.D0
      B1LK21(I)=0.D0
      B1LK22(I)=0.D0
      B2LK31(I)=0.D0
      B2LK32(I)=0.D0
      B2LK33(I)=0.D0
      B0O1(I)=0.D0
      B1O21(I)=0.D0
      B1O22(I)=0.D0
      BO2(I)=0.D0
      B0O2(I)=0.D0
      B2O31(I)=0.D0
      B2O32(I)=0.D0
      B2O33(I)=0.D0
   11 CONTINUE
      IF(MESHA.EQ.1)GO TO 5
      BET42=BET4*BET4
      IF(MESHA.GT.2) GO TO 82
      A40=1.D0
      A42=0.D0
      A44=0.D0
      GO TO 5
   82 C7=7.D0/12.D0
      C5=5.D0/12.D0
      CQ5=SQRT(C5)
      CQ7=SQRT(C7)
      IF(MESHA.EQ.4) CSDG=COS(DELG)
      IF(MESHA.EQ.4) GO TO 7
      GAMG=GAM0
      CSDG=CQ7*COS(3.D0*GAM0)
    7 SSDG=SQRT(1.D0-CSDG*CSDG)
      CSGG=COS(GAMG)
      SSGG=SIN(GAMG)
      A40=CQ7*CSDG+CQ5*SSDG*CSGG
      A42=-SSDG*SSGG/SQ2
      A44=(CQ5*CSDG-CQ7*SSDG*CSGG)/SQ2
    5 ALK11(2)=1.D0
      ALK11(3)=1.D0/SQ2
      IF(LAS.EQ.1) GO TO 21
      DO 15 JJ2=1,3
      JI2=4*(JJ2-1)
      J1=4
      J2=4
      J=JI2
      M1=0
      M2=0
      M=0
      CALL KLEGO
      AK1=AKG
      DO 16 KK2=1,JJ2
      LK=NLA(JJ2)+KK2
      GO TO (17,18,19),KK2
   17 M1=4
      M2=-4
      CALL KLEGO
      AK2=AKG
      ALK21(LK)=AK1*AK1
      ALK22(LK)=AK1*AK2
      GO TO 16
   18 M1=0
      M2=4
      M=4
      CALL KLEGO
      AK4=AKG
      ALK21(LK)=AK1*AK4*SQ2
      GO TO 16
   19 M1=4
      M=8
      CALL KLEGO
      AK3=AKG
      ALK21(LK)=AK1*AK3/2.D0
   16 CONTINUE
      IF(LAS.EQ.2) GO TO 15
      J3M=IABS(JJ2-2)+1
      J3P=JJ2+1
      DO 46 JJ3=J3M,J3P
      JI3=4*(JJ3-1)
      J2=JI2
      M1=0
      M2=0
      J=JI3
      M=0
      CALL KLEGO
      BK1=AKG
      ABK=AK1*BK1
      DO 45 KK3=1,JJ3
      LK=NLA(JJ3)+KK3
      GO TO (36,37,38,39), KK3
   36 M1=-4
      M2=4
      CALL KLEGO
      BK2=AKG
      ALK31(LK)=ALK31(LK)+ABK*AK1*BK1
      ALK32(LK)=ALK32(LK)+ABK*(BK1*AK2+AK4*BK2)
      GO TO 45
   37 M1=0
      M=4
      CALL KLEGO
      BK6=AKG
      M1=4
      M2=0
      CALL KLEGO
      BK7=AKG
      M1=-4
      M2=8
      CALL KLEGO
      BK9=AKG
      ALK31(LK)=ALK31(LK)+ABK*(SQ2*AK4*BK6+AK1*BK7/SQ2)
      ALK32(LK)=ALK32(LK)+ABK*(AK3*BK9/SQ32+AK2*BK7/SQ2)
      GO TO 45
   38 M1=0
      M=8
      CALL KLEGO
      BK5=AKG
      M1=4
      M2=4
      CALL KLEGO
      BK3=AKG
      ALK31(LK)=ALK31(LK)+ABK*(AK3*BK5/2.D0+AK4*BK3)
      GO TO 45
   39 M2=8
      M=12
      CALL KLEGO
      BK8=AKG
      ALK31(LK)=ALK31(LK)+ABK*AK3*BK8
   45 CONTINUE
      IF(LAS.EQ.3) GO TO 46
      J4M=IABS(JJ3-1)+1
      J4P=JJ3+1
      DO 40 JJ4=J4M,J4P
      JI4=4*(JJ4-1)
      J2=JI3
      M1=0
      M2=0
      J=JI4
      M=0
      CALL KLEGO
      CK1=AKG
      ABCK=ABK*CK1
      DO 49 KK4=1,JJ4
      LK=NLA(JJ4)+KK4
      GO TO (50,51,52,53,54),KK4
   50 M1=-4
      M2=4
      CALL KLEGO
      CK2=AKG
      ALK41(LK)=ALK41(LK)+ABCK*CK1*BK1*AK1
      ALK42(LK)=ALK42(LK)+ABCK*(CK1*(BK1*AK2+AK4*BK2)+CK2*
     *(AK4*BK6+BK7*AK1/2.D0))
      ALK43(LK)=ALK43(LK)+ABCK*CK2*(AK3*BK9/4.D0+BK7*AK2/2.D0)
      GO TO 49
   51 M1=0
      M=4
      CALL KLEGO
      CK3=AKG
      M1=4
      M2=0
      CALL KLEGO
      CK4=AKG
      M1=-4
      M2=8
      CALL KLEGO
      CK5=AKG
      ALK41(LK)=ALK41(LK)+ABCK*(CK3*(AK4*BK6*SQ2+
     *BK7*AK1/SQ2)+CK4*BK1*AK1/SQ2)
      ALK42(LK)=ALK42(LK)+ABCK*(CK3*(AK3*BK9/SQ32+BK7*AK2/SQ2)+
     *CK5*(AK3*BK5/SQ32+AK4*BK3/SQ2)+CK4*(BK1*AK2/SQ2+AK4*BK2/SQ2))
      GO TO 49
   52 M1=0
      M=8
      CALL KLEGO
      CK6=AKG
      M1=4
      M2=4
      CALL KLEGO
      CK7=AKG
      M1=-4
      M2=12
      CALL KLEGO
      CK8=AKG
      ALK41(LK)=ALK41(LK)+ABCK*(CK6*(AK3*BK5/2.D0+AK4*BK3)+
     *CK7*(AK4*BK6+BK7*AK1/2.D0))
      ALK42(LK)=ALK42(LK)+ABCK*(CK7*(AK3*BK9/4.D0+BK7*AK2/2.D0)+
     *CK8*AK3*BK8/4.)
      GO TO 49
   53 M1=0
      M=12
      CALL KLEGO
      CK11=AKG
      M1=4
      M2=8
      CALL KLEGO
      CK10=AKG
      ALK41(LK)=ALK41(LK)+ABCK*(CK11*AK3*BK8/SQ32+
     *CK10*(AK3*BK5/SQ32+AK4*BK3/SQ2))
      GO TO 49
   54 M2=12
      M=16
      CALL KLEGO
      CK9=AKG
      ALK41(LK)=ALK41(LK)+ABCK*CK9*AK3*BK8/4.D0
   49 CONTINUE
   40 CONTINUE
   46 CONTINUE
   15 CONTINUE
   21 IF(MESHA.EQ.1) GO TO 60
      B0LK1(4)=A40*BET4
      B0LK1(5)=A42*BET4
      B0LK1(6)=A44*BET4
      IF(LAS.EQ.1) GO TO 60
      DO 61 JJ2=1,5
      JI2=4*(JJ2-1)
      J1=8
      J2=8
      J=JI2
      M1=0
      M2=0
      M=0
      CALL KLEGO
      CK1=AKG
      DO 62 KK2=1,JJ2
      LK=NLA(JJ2)+KK2
      GO TO (63,64,65,66,67),KK2
   63 M1=4
      M2=-4
      CALL KLEGO
      CK2=AKG
      M1=8
      M2=-8
      CALL KLEGO
      CK3=AKG
      B0LK2(LK)=CK1*(A40*A40*CK1+2.D0*A42*A42*CK2+
     *2.D0*A44*A44*CK3)*C95*BET42
      GO TO 62
   64 M1=0
      M2=4
      M=4
      CALL KLEGO
      CK4=AKG
      M1=-4
      M2=8
      CALL KLEGO
      CK5=AKG
      B0LK2(LK)=CK1*2.D0*(A40*A42*CK4+A42*A44*CK5)*C95*BET42
      GO TO 62
   65 M1=4
      M2=4
      M=8
      CALL KLEGO
      CK6=AKG
      M1=0
      M2=8
      CALL KLEGO
      CK7=AKG
      B0LK2(LK)=CK1*(A42*A42*CK6+2.D0*A40*A44*CK7)*C95*BET42
      GO TO 62
   66 M1=4
      M=12
      CALL KLEGO
      CK8=AKG
      B0LK2(LK)=CK1*2.D0*A42*A44*CK8*C95*BET42
      GO TO 62
   67 M1=8
      M=16
      CALL KLEGO
      CK9=AKG
      B0LK2(LK)=CK1*A44*A44*CK9*C95*BET42
   62 CONTINUE
      IF(JJ2.LT.2.OR.JJ2.GT.4) GO TO 61
      J1=4
      J2=8
      M1=0
      M2=0
      M=0
      CALL KLEGO
      BK1=AKG
      DO 68 KK2=1,JJ2
      LK=NLA(JJ2)+KK2
      GO TO (69,70,71,72),KK2
   69 M1=4
      M2=-4
      CALL KLEGO
      BK2=AKG
      B1LK21(LK)=2.D0*BK1*A40*BK1*CQ95*BET4
      B1LK22(LK)=2.D0*BK1*A42*BK2*CQ95*BET4*SQ2
      GO TO 68
   70 M1=0
      M2=4
      M=4
      CALL KLEGO
      BK3=AKG
      M1=4
      M2=0
      CALL KLEGO
      BK4=AKG
      M1=-4
      M2=8
      CALL KLEGO
      BK5=AKG
      B1LK21(LK)=2.D0*BK1*A42*BK3*CQ95*BET4
      B1LK22(LK)=2.D0*BK1*(A40*BK4+A44*BK5)*CQ95*BET4/SQ2
      GO TO 68
   71 M1=0
      M2=8
      M=8
      CALL KLEGO
      BK6=AKG
      M1=4
      M2=4
      CALL KLEGO
      BK7=AKG
      B1LK21(LK)=2.D0*BK1*A44*BK6*CQ95*BET4
      B1LK22(LK)=2.D0*BK1*A42*BK7*CQ95*BET4/SQ2
      GO TO 68
   72 M2=8
      M=12
      CALL KLEGO
      BK8=AKG
      B1LK22(LK)=2.D0*BK1*A44*BK8*CQ95*BET4/SQ2
   68 CONTINUE
   61 CONTINUE
      IF(LAS.EQ.2) GO TO 60
      DO 73 JJ2=1,3
      JMI=IABS(JJ2-3)+1
      JMA=JJ2+2
      JI2=4*(JJ2-1)
      J1=4
      J2=4
      J=JI2
      M1=0
      M2=0
      M=0
      CALL KLEGO
      AK1=AKG
      M1=4
      M2=-4
      CALL KLEGO
      AK2=AKG
      M1=0
      M2=4
      M=4
      CALL KLEGO
      AK3=AKG
      M1=4
      M=8
      CALL KLEGO
      AK4=AKG
      DO 74 JJ3=JMI,JMA
      JI3=4*(JJ3-1)
      J1=JI2
      J2=8
      J=JI3
      M1=0
      M2=0
      M=0
      CALL KLEGO
      DK1=AKG
      DO 75 KK3=1,JJ3
      LK=NLA(JJ3)+KK3
      GO TO (76,77,78,79,80),KK3
   76 M1=4
      M2=-4
      CALL KLEGO
      DK2=AKG
      M1=8
      M2=-8
      CALL KLEGO
      DK3=AKG
      B2LK31(LK)=B2LK31(LK)+3.D0*AK1*DK1*AK1*A40*DK1*CQ95*BET4
      B2LK32(LK)=B2LK32(LK)+3.D0*AK1*DK1*(AK2*A40*DK1+
     *AK4*A44*DK3)*CQ95*BET4
      B2LK33(LK)=B2LK33(LK)+3.D0*AK1*DK1*AK3*A42*DK2*CQ95*BET4*SQ2*2.D0
      GO TO 75
   77 M1=4
      M2=0
      M=4
      CALL KLEGO
      DK4=AKG
      M1=0
      M2=4
      CALL KLEGO
      DK5=AKG
      M1=8
      M2=-4
      CALL KLEGO
      DK6=AKG
      M1=-4
      M2=8
      CALL KLEGO
      DK7=AKG
      B2LK31(LK)=B2LK31(LK)+3.D0*AK1*DK1*AK1*A42*DK5*CQ95*BET4
      B2LK32(LK)=B2LK32(LK)+3.D0*AK1*DK1*(AK2*A42*DK5+
     *AK4*A42*DK6/2.)*CQ95*BET4
      B2LK33(LK)=B2LK33(LK)+3.D0*AK1*DK1*(AK3*A40*DK4+
     *AK3*A44*DK7)*CQ95*BET4*SQ2
      GO TO 75
   78 M1=8
      M2=0
      M=8
      CALL KLEGO
      DK8=AKG
      M1=4
      M2=4
      CALL KLEGO
      DK9=AKG
      M1=0
      M2=8
      CALL KLEGO
      DK10=AKG
      B2LK31(LK)=B2LK31(LK)+3.D0*AK1*DK1*AK1*A44*DK10*CQ95*BET4
      B2LK32(LK)=B2LK32(LK)+3.D0*AK1*DK1*(AK3*A44*DK10+
     *AK4*A40*DK8/2.D0)*CQ95*BET4
      B2LK33(LK)=B2LK33(LK)+3.D0*AK1*DK1*AK3*A42*DK9*CQ95*BET4*SQ2
      GO TO 75
   79 M1=8
      M2=4
      M=12
      CALL KLEGO
      DK11=AKG
      M1=4
      M2=8
      CALL KLEGO
      DK12=AKG
      B2LK32(LK)=B2LK32(LK)+3.D0*AK1*DK1*AK4*A42*DK11*CQ95*BET4/2.D0
      B2LK33(LK)=B2LK33(LK)+3.D0*AK1*DK1*AK3*A44*DK12*CQ95*BET4*SQ2
      GO TO 75
   80 M1=8
      M=16
      CALL KLEGO
      DK13=AKG
      B2LK32(LK)=B2LK32(LK)+3.D0*AK1*DK1*AK4*A44*DK13*CQ95*BET4/2.D0
   75 CONTINUE
   74 CONTINUE
   73 CONTINUE
   60 IF(MESHO.EQ.0) GO TO 2
      GO TO (1,6),MESHO
    1 A30=1.D0
      A32=0.D0
      GO TO 8
    6 A30=COS(ETO)
      A32=SIN(ETO)/C1
    8 B0O1(2)=A30
      B0O1(3)=A32
      IF(LAS.EQ.1) GO TO 2
      DO 9 J32=1,4
      IF(J32.EQ.4) GO TO 10
      JI32=4*(J32-1)+2
      J1=4
      J2=6
      J=JI32
      M1=0
      M2=0
      M=0
      CALL KLEGO
      OK1=AKG
      DO 20 K32=1,J32
      LK=NLA(J32)+K32
      GO TO (22,23,24),K32
   22 M1=4
      M2=-4
      CALL KLEGO
      OK2=AKG
      B1O21(LK)=2.D0*OK1**2*A30*CSQ75
      B1O22(LK)=OK1*OK2*C1*A32*CSQ75*2.D0
      GO TO 20
   23 M2=0
      M=4
      CALL KLEGO
      OK3=AKG
      M1=0
      M2=4
      CALL KLEGO
      OK4=AKG
      B1O21(LK)=2.D0*OK1*OK4*A32*CSQ75
      B1O22(LK)=OK1*OK3*C1*A30*CSQ75
      GO TO 20
   24 M1=4
      M=8
      CALL KLEGO
      OK5=AKG
      B1O22(LK)=OK1*OK5*C1*A32*CSQ75
   20 CONTINUE
   10 IF(MESHA.EQ.1) GO TO 55
      J1=6
      J1=8
      J=JI32
      M1=0
      M2=0
      M=0
      CALL KLEGO
      OL1=AKG
      DO 26 K34=1,J32
      LK=NLA(J32)+K34
      GO TO (27,28,29,30),K34
   27 M1=4
      M2=-4
      CALL KLEGO
      OL2=AKG
      B0O2(LK)=2.D0*OL1*(OL1*A30*A40+2.D0*OL2*A32*A42)*C6325*BET4
      GO TO 26
   28 M2=0
      M=4
      CALL KLEGO
      OL3=AKG
      M1=0
      M2=4
      CALL KLEGO
      OL5=AKG
      M1=-4
      M2=8
      CALL KLEGO
      OL4=AKG
      B0O2(LK)=2.D0*OL1*(OL3*A32*A40+OL5*A30*A42+OL4*A32*A44)*C6325
     **BET4
      GO TO 26
   29 M1=4
      M2=4
      M=8
      CALL KLEGO
      OL6=AKG
      M1=0
      M2=8
      CALL KLEGO
      OL7=AKG
      B0O2(LK)=2.D0*OL1*(OL6*A32*A42+OL7*A30*A44)*C6325*BET4
      GO TO 26
   30 M2=4
      M=12
      CALL KLEGO
      OL8=AKG
      B0O2(LK)=2.D0*OL1*(OL8*A32*A44)*C6325*BET4
   26 CONTINUE
   55 JI33=4*(J32-1)
      J1=6
      J2=6
      J=JI33
      M1=0
      M2=0
      M=0
      CALL KLEGO
      OM1=AKG
      DO 31 K33=1,J32
      LK=NLA(J32)+K33
      GO TO (32,33,34,31),K33
   32 M1=4
      M2=-4
      CALL KLEGO
      OM2=AKG
      BO2(LK)=OM1*(A30**2*OM1+2.D0*OM2*A32**2)*C75
      GO TO 31
   33 M1=0
      M2=4
      M=4
      CALL KLEGO
      OM3=AKG
      BO2(LK)=OM1*A30*A32*2.D0*OM3*C75
      GO TO 31
   34 M1=4
      M=8
      CALL KLEGO
      OM4=AKG
      BO2(LK)=OM1*A32**2*OM4*C75
   31 CONTINUE
      IF(LAS.EQ.2) GO TO 9
      IF(J32.EQ.4) GO TO 2
      J1=4
      J2=4
      J=4*(J32-1)
      M1=0
      M2=0
      M=0
      CALL KLEGO
      ON1=AKG
      M1=4
      M2=-4
      CALL KLEGO
      ON2=AKG
      M1=0
      M2=4
      M=4
      CALL KLEGO
      ON3=AKG
      M1=4
      M=8
      CALL KLEGO
      ON4=AKG
      JMI=1
      JMA=J32+1
      IF(J32.EQ.1) JMI=2
      DO 35 J223=JMI,JMA
      JI223=4*(J223-1)+2
      J1=4*(J32-1)
      J2=6
      J=JI223
      M1=0
      M2=0
      M=0
      CALL KLEGO
      OM1=AKG
      DO 41 K223=1,J223
      LK=NLA(J223)+K223
      GO TO (42,43,44,48),K223
   42 M1=4
      M2=-4
      CALL KLEGO
      OM2=AKG
      B2O31(LK)=B2O31(LK)+ON1*OM1*ON1*OM1*A30*3.D0*CSQ75
      B2O32(LK)=B2O32(LK)+ON1*OM1*A30*ON2*OM1*3.D0*CSQ75
      B2O33(LK)=B2O33(LK)+ON1*OM1*A32*ON3*OM2*C1*3.D0*CSQ75
      GO TO 41
   43 M1=4
      M2=0
      M=4
      CALL KLEGO
      OM3=AKG
      M1=0
      M2=4
      CALL KLEGO
      OM4=AKG
      M1=8
      M2=-4
      CALL KLEGO
      OM5=AKG
      B2O31(LK)=B2O31(LK)+ON1*OM1*ON1*OM4*A32*3.D0*CSQ75
      B2O32(LK)=B2O32(LK)+ON1*OM1*(ON2*OM4+ON4*OM5/2.D0)*A32*3.D0*CSQ75
      B2O33(LK)=B2O33(LK)+ON1*OM1*A30*ON3*OM3/C1*3.D0*CSQ75
      GO TO 41
   44 M2=0
      M=8
      CALL KLEGO
      OM6=AKG
      M1=4
      M2=4
      CALL KLEGO
      OM7=AKG
      B2O31(LK)=B2O31(LK)+ON1*OM1*ON1*OM7*A32*3.D0*CSQ75
      B2O32(LK)=B2O32(LK)+ON1*OM1*(ON2*OM7*A32+ON4*OM6*A30/2.D0)*3.D0*
     *CSQ75
      GO TO 41
   48 M1=8
      M=12
      CALL KLEGO
      OM8=AKG
      B2O32(LK)=B2O32(LK)+ON1*OM1*ON4*OM8*A32*1.5D0*CSQ75
   41 CONTINUE
   35 CONTINUE
    9 CONTINUE
    2 CONTINUE
      RETURN
      END
C     *******************************************************
      SUBROUTINE KNDIT
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
C     BELOW CARD IS NO NECESSARY IF MEMORY IS MORE THAN 32Mb
      DOUBLE PRECISION PVC,PRC
      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM1.FOR'
      INCLUDE 'PRIVCOM7.FOR'  
      INCLUDE 'PRIVCOM8.FOR'
      INCLUDE 'PRIVCOM9.FOR' 
      INCLUDE 'PRIVCOM10.FOR' 
      INCLUDE 'PRIVCOM12.FOR'
      INCLUDE 'PRIVCOM13.FOR'
      INCLUDE 'PRIVCOM20.FOR'
      DIMENSION NLA(5)
      DATA NLA/0,1,3,6,10/
      
      sum2=sum(WEIGH**2)
      defnul2=sum(effdef*effdef)
      
      
      NSPIP=IDINT(ASP*2.D0+0.01)
      SQ2=DSQRT(2.D0)
C     SQ32=2.828427D0
      SQ32=SQ2**3
      C95=1.8D0
      CQ95=1.34164D0
      IILL=NCLL*NCLL*2*9
      NCLLX2=NCLL*NCLL
      NCLL4=NCLLX2*4
      DO 4 I=1,IILL
      IF(I.LE.NCLLX2) CVNRPN(I)=0.D0
      IF(I.LE.NCLL4) CVNR(I)=0.D0
    4 CVNC(I)=0.D0
      LAS2=LAS
      ICLL=NCLL*9*LAS
      LAS8=LAS*9
          LASC=1
          LAS8=9
          IF(LAS2.GE.2) LASC=2
          IF(LAS2.GE.2) LAS8=18
          IF(MEPOT.EQ.1) GO TO 10
      ICLL=NCLL*LAS8
      NCLA=NCLL*LAS2
      DO 20 K=1,NCLL
      K1=(K-1)*NCLA
      KI1=(K-1)*LAS2
      K1PN=(K-1)*NCLL
      KC1=ICLL*(K-1)
      KIC1=LAS8*(K-1)
      NU1=NNJ1(K)
      JO1=JO(NU1)
      LN1=LNJ1(K)*2
      JN1=JNJ1(K)
      NPI1=(-1)**NNO(NU1)
      JO12=JO1/2
      KMI1=1-(-1)**JO12*NPI1
      KMA1=2*(JO12/2)
      NK1=(KMA1-KMI1+2)/2
      DO 20 KK=K,NCLL
      K2=K1+(KK-1)*LAS2
      K2PN=K1PN+KK
      K2NP=(KK-1)*NCLL+K
      KI2=KI1+(KK-1)*NCLA
      KC2=KC1+(KK-1)*LAS8
      KIC2=KIC1+(KK-1)*ICLL
      NU2=NNJ1(KK)
      JO2=JO(NU2)
      LN2=LNJ1(KK)*2
      JN2=JNJ1(KK)
      JO22=JO2/2
      NPI2=(-1)**NNO(NU2)
      KMI2=1-(-1)**JO22*NPI2
      KMA2=2*(JO22/2)
      NK2=(KMA2-KMI2+2)/2
      AA=(-1)**((JSS-JO2-NSPIP+LN1+LN2)/2+(LN2-LN1)/4)*SQRT((LN1+1.D0)*
     *(LN2+1.D0)*(JN1+1.D0)*(JN2+1.D0)*(JO1+1.D0)/12.56663708D0)
      NPCV=1
      IF(LN2-LN1.LT.0) NPCV=-1
      FOK1=FOK(NNO(NU1)+1,NNO(NU2)+1,1)
      FOK2=FOK(NNO(NU1)+1,NNO(NU2)+1,2)
      IF(BET3.NE.0.) FOK1=FOK1*BET3
      IF(BET3.NE.0.) FOK2=FOK2*BET3**2
      IF(BET3.EQ.0.) FOK1=FOK1*AMUO
      IF(BET3.EQ.0.) FOK2=FOK2*AMUO**2
      CG1=TRIG(NU1,NU2,1)
      SG1=TRIG(NU1,NU2,2)
      CG2=TRIG(NU1,NU2,3)
      SG2=TRIG(NU1,NU2,4)
      CSG11=TRIG(NU1,NU2,5)
      CG3=TRIG(NU1,NU2,6)
      SG3=TRIG(NU1,NU2,7)
      CSG21=TRIG(NU1,NU2,8)
      CSG12=TRIG(NU1,NU2,9)
      CG4=TRIG(NU1,NU2,10)
      SG4=TRIG(NU1,NU2,11)
      CSG22=TRIG(NU1,NU2,12)
      CSG31=TRIG(NU1,NU2,13)
      CSG13=TRIG(NU1,NU2,14)
      FOIJ1=FOV(NU1,NU2,1)*BET0
      IF(LAS2.GE.2) FOIJ2=FOV(NU1,NU2,2)*BET0**2
      IF(LAS2.GE.3) FOIJ3=FOV(NU1,NU2,3)*BET0**3
      DO 22 L=1,LAS2
      AAA=AA*FOV(NU1,NU2,L)*BET0**L
      LL2=K2+L
      LI2=KI2+L
      LLC2=KC2+L
      LIC2=KIC2+L
      SU1=0.
      DO 23 LA=1,5
C     NEXT TWO CARDS JUST TO OVERCOME THE ERROR OR  MICROSOFT FPS COMPILER
C     v.1.00, YOU CAN DELETE THEM USING OTHER COMPILERS
      IF(LA.GT.5 .AND. MEPRI.LT.98) PRINT 9911, LA
 9911 FORMAT (I3)
      IF(L.EQ.1.AND.LA.NE.2) GO TO 83
      IF(L.EQ.2.AND.LA.GE.4) GO TO 83
      IF(L.EQ.3.AND.LA.GE.5) GO TO 83
      GO TO 84
   83 IF(MESHA.LT.1) GO TO 3
      IF(L.EQ.1.AND.LA.NE.3) GO TO 3
      GO TO 84
    3 IF(MESHO.EQ.0) GO TO 23
      IF(L.EQ.1.AND.LA.NE.2) GO TO 23
      IF(LA.GT.4) GO TO 23
   84 LAM=4*(LA-1)
      SU2=0.
      DO 25 N=1,LA
      ALK=0.
      LK=NLA(LA)+N
      IF(NPI1.NE.NPI2) GO TO 11
      GO TO (27,28,29,30),L
   27 GO TO (31,32),N
   31 ALK=ALK11(LK)*CG1
      IF(MESHA.LE.1) GO TO 33
      ALK=ALK+B0LK1(LK)/FOIJ1
      GO TO 33
   32 ALK=ALK11(LK)*SG1
      IF(MESHA.LE.1) GO TO 33
      ALK=ALK+B0LK1(LK)/FOIJ1
      GO TO 33
   28 GO TO (34,35,48),N
   34 ALK=ALK21(LK)*CG2+ALK22(LK)*SG2
      IF(MESHA.LE.1) GO TO 33
      ALK=ALK+B0LK2(LK)/FOIJ2+(B1LK21(LK)*CG1+B1LK22(LK)*SG1)*
     *FOIJ1/FOIJ2
      IF(MESHO.EQ.0) GO TO 33
      IF(MEHAO.EQ.2) GO TO 61
      ALK=ALK+BO2(LK)*FOK2/FOIJ2
      GO TO 33
   61 ALK=ALK+BO2(LK)*FOK2
      GO TO 33
   35 ALK=ALK21(LK)*CSG11
      IF(MESHA.LE.1) GO TO 33
      ALK=ALK+B0LK2(LK)/FOIJ2+(B1LK21(LK)*CG1+B1LK22(LK)*SG1)*
     *FOIJ1/FOIJ2
      IF(MESHO.EQ.0) GO TO 33
      IF(MEHAO.EQ.2) GO TO 62
      ALK=ALK+BO2(LK)*FOK2/FOIJ2
      GO TO 33
   62 ALK=ALK+BO2(LK)*FOK2
      GO TO 33
   48 ALK=ALK21(LK)*SG2
      IF(MESHA.LE.1) GO TO 33
      ALK=ALK+B0LK2(LK)/FOIJ2+(B1LK21(LK)*CG1+B1LK22(LK)*SG1)*
     *FOIJ1/FOIJ2
      IF(MESHO.EQ.0) GO TO 33
      IF(MEHAO.EQ.2) GO TO 63
      ALK=ALK+BO2(LK)*FOK2/FOIJ2
      GO TO 33
   63 ALK=ALK+BO2(LK)*FOK2
      GO TO 33
   29 GO TO (41,42,43,44),N
   41 ALK=ALK31(LK)*CG3+ALK32(LK)*CSG12
      IF(MESHA.LE.1) GO TO 33
      ALK=ALK+(B2LK31(LK)*CG2+B2LK32(LK)*SG2+B2LK33(LK)*CSG11)*
     *FOIJ2/FOIJ3
      GO TO 33
   42 ALK=ALK31(LK)*CSG21+ALK32(LK)*SG3
      IF(MESHA.LE.1) GO TO 33
      ALK=ALK+(B2LK31(LK)*CG2+B2LK32(LK)*SG2+B2LK33(LK)*CSG11)*
     *FOIJ2/FOIJ3
      GO TO 33
   43 ALK=ALK31(LK)*CSG12
      IF(MESHA.LE.1) GO TO 33
      ALK=ALK+(B2LK31(LK)*CG2+B2LK32(LK)*SG2+B2LK33(LK)*CSG11)*
     *FOIJ2/FOIJ3
      GO TO 33
   44 ALK=ALK31(LK)*SG3
      IF(MESHA.LE.1) GO TO 33
      ALK=ALK+(B2LK31(LK)*CG2+B2LK32(LK)*SG2+B2LK33(LK)*CSG11)*
     *FOIJ2/FOIJ3
   30 GO TO (55,56,57,58,59),N
   55 ALK=ALK41(LK)*CG4+ALK42(LK)*CSG22+ALK43(LK)*SG4
      GO TO 33
   56 ALK=ALK41(LK)*CSG31+ALK42(LK)*CSG13
      GO TO 33
   57 ALK=ALK41(LK)*CSG22+ALK42(LK)*SG4
      GO TO 33
   58 ALK=ALK41(LK)*CSG13
      GO TO 33
   59 ALK=ALK41(LK)*SG4
      GO TO 33
   11 LAM=4*(LA-1)+2
      IF(MESHO.EQ.0.OR.MEHAO.EQ.0) GO TO 33
      GO TO (16,17,18,23),L
   16 GO TO (5,7),N
    5 IF(MEHAO.EQ.2) GO TO 64
      ALK=B0O1(LK)*FOK1/FOIJ1
      GO TO 33
   64 ALK=B0O1(LK)*FOK1
      GO TO 33
    7 IF(MEHAO.EQ.2) GO TO 65
      ALK=B0O1(LK)*FOK1/FOIJ1
      GO TO 33
   65 ALK=B0O1(LK)*FOK1
      GO TO 33
   17 GO TO (36,37,38,39),N
   36 IF(MEHAO.EQ.2) GO TO 66
      ALK=((B1O21(LK)*CG1+B1O22(LK)*SG1)*FOIJ1+
     *B0O2(LK))/FOIJ2*FOK1
      GO TO 33
   66 ALK=((B1O21(LK)*CG1+B1O22(LK)*SG1)*FOIJ2+
     *B0O2(LK)*FOIJ1)/FOIJ2*FOK1                                        
      GO TO 33
   37 IF(MEHAO.EQ.2) GO TO 67
      ALK=((B1O21(LK)*CG1+B1O22(LK)*SG1)*FOIJ1+
     *B0O2(LK))/FOIJ2*FOK1                                              
      GO TO 33
   67 ALK=((B1O21(LK)*CG1+B1O22(LK)*SG1)*FOIJ2+
     *B0O2(LK)*FOIJ1)/FOIJ2*FOK1                                        
      GO TO 33
   38 IF(MEHAO.EQ.2) GO TO 68
      ALK=((B1O21(LK)*CG1+B1O22(LK)*SG1)*FOIJ1+
     *B0O2(LK))/FOIJ2*FOK1                                              
      GO TO 33
   68 ALK=((B1O21(LK)*CG1+B1O22(LK)*SG1)*FOIJ2+
     *B0O2(LK)*FOIJ1)/FOIJ2*FOK1                                        
      GO TO 33
   39 IF(MEHAO.EQ.2) GO TO 69
      ALK=((B1O21(LK)*CG1+B1O22(LK)*SG1)*FOIJ1+
     *B0O2(LK))/FOIJ2*FOK1                                              
      GO TO 33
   69 ALK=((B1O21(LK)*CG1+B1O22(LK)*SG1)*FOIJ2+
     *B0O2(LK)*FOIJ1)/FOIJ2*FOK1                                        
      GO TO 33
   18 GO TO (50,51,52,53),N
   50 IF(MEHAO.EQ.2) GO TO 70
      ALK=(B2O31(LK)*CG2+B2O32(LK)*SG2+B2O33(LK)*CSG11*2.D0)*
     *FOK1*FOIJ2/FOIJ3
      GO TO 33
   70 ALK=(B2O31(LK)*CG2+B2O32(LK)*SG2+B2O33(LK)*CSG11*2.D0)*
     *FOK1
      GO TO 33
   51 IF(MEHAO.EQ.2) GO TO 71
      ALK=(B2O31(LK)*CG2+B2O32(LK)*SG2+B2O33(LK)*CSG11*2.D0)*
     *FOK1*FOIJ2/FOIJ3
      GO TO 33
   71 ALK=(B2O31(LK)*CG2+B2O32(LK)*SG2+B2O33(LK)*CSG11*2.D0)*
     *FOK1
   52 IF(MEHAO.EQ.2) GO TO 72
      ALK=(B2O31(LK)*CG2+B2O32(LK)*SG2+B2O33(LK)*CSG11*2.D0)*
     *FOK1*FOIJ2/FOIJ3
      GO TO 33
   72 ALK=(B2O31(LK)*CG2+B2O32(LK)*SG2+B2O33(LK)*CSG11*2.D0)*
     *FOK1
      GO TO 33
   53 IF(MEHAO.EQ.2) GO TO 73
      ALK=(B2O31(LK)*CG2+B2O32(LK)*SG2+B2O33(LK)*CSG11*2.D0)*
     *FOK1*FOIJ2/FOIJ3
      GO TO 33
   73 ALK=(B2O31(LK)*CG2+B2O32(LK)*SG2+B2O33(LK)*CSG11*2.D0)*
     *FOK1
   33 IF(ALK.EQ.0.D0) GO TO 25
      SU=0.
      J1=LN1
      J2=LN2
      J=LAM
      M1=0
      M2=0
      M=0
      CALL KLEGO
      AB=AKG
      JA=JN1
      JB=JO1
      JE=JN2
      JD=JO2
      JC=JSS
      JF=LAM
      CALL RACAH
      AB=AB*W
      JA=LN1
      JB=JN1
      JE=LN2
      JD=JN2
      JC=NSPIP
      JF=LAM
      CALL RACAH
      AB=AB*W
      DO 24 LK1=1,NK1
      KS1=2*(KMI1+2*(LK1-1))
      AIT1=AIT(NU1,LK1)
      IF(KS1.EQ.0) AIT1=AIT1/SQ2
      DO 24 LK2=1,NK2
      KS2=2*(KMI2+2*(LK2-1))
      AIT2=AIT(NU2,LK2)
      IF(KS2.EQ.0) AIT2=AIT2/SQ2
      KA=4*(N-1)
      J1=JO1
      J2=LAM
      J=JO2
      M1=KS1
      M2=KA
      M=KS2
      CALL KLEGO
      CK1=AKG
      M1=-KS1
      CALL KLEGO
      CLK=CK1+(-1)**JO12*AKG*NPI1
      IF(KA.EQ.0) GO TO 26
      M1=KS1
      M2=-KA
      CALL KLEGO
      CLK=CLK+AKG
   26 SU=SU+CLK*AIT1*AIT2
   24 CONTINUE
      SU2=SU2+SU*ALK
   25 CONTINUE
      IF(L.EQ.2) SU2=SU2*1.4104730D0/SQRT(LAM+1.D0)
      IF(L.EQ.3) SU2=SU2*0.8897014D0/SQRT(LAM+1.D0)
      IF(L.EQ.4) SU2=SU2*0.5612105D0/SQRT(LAM+1.D0)
      SU1=SU1+SU2*AB
      IF(L.GT.2) GO TO 23
      LLLC2=LLC2+(LAM/2)*LASC
      LLIC2=LIC2+(LAM/2)*LASC
      CVNC(LLLC2)=SU2*AB*AAA
C     CVNC(LLLC2)=0.
      IF(KK.NE.K) CVNC(LLIC2)=CVNC(LLLC2)
      IF(NPI1.NE.NPI2) CVNC(LLLC2)=CVNC(LLLC2)*NPCV
      IF(NPI1.NE.NPI2)CVNC(LLIC2)=-CVNC(LLLC2)
C      if(lam/2.eq.0) write (21,9999)jo1,jo2,l,kc2,k,kk,llc2,cvnc(lllc2)
C9999 format(7i6,e12.4)
   23 CONTINUE
      CVNR(LL2)=SU1*AAA
      IF(KK.NE.K) CVNR(LI2)=CVNR(LL2)
      IF(NPI1.NE.NPI2) CVNR(LL2)=CVNR(LL2)*NPCV
      IF(NPI1.NE.NPI2)CVNR(LI2)=-CVNR(LL2)
      
      
 
   22 CONTINUE
   20 CONTINUE


      DO 115 K=1,NCLL
      K1PN=(K-1)*NCLL
      NU=NNJ1(K)
      JO1=JO(NU)
      LN1=LNJ1(K)*2
      NPI1=(-1)**NNO(NU)
      DO 115 KK=1,NCLL
      K2PN=K1PN+KK
      NU1=NNJ1(KK)
      JO2=JO(NU1)
      LN2=LNJ1(KK)*2
      NPI2=(-1)**NNO(NU1)
      CVNRPN(K2PN)=0.D0

      IF(JO1.EQ.JO2.AND.NCA(NU).NE.NCA(NU1).AND.NTU(NU).EQ.NTU(NU1)
     *.AND.NNB(NU).EQ.NNB(NU1).AND.NNG(NU).EQ.NNG(NU1)
     *.AND.NNO(NU).EQ.NNO(NU1))GO TO 121
      GO TO 115
  121 CVNRPN(K2PN)=1.D0
      NPCV=1
C     IF(LN2-LN1.LT.0) NPCV=-1
C     IF(NPI1.NE.NPI2) CVNRPN(K2PN)=CVNRPN(K2PN)*NPCV
C     IF(NPI1.NE.NPI2)CVNRPN(K2PN)=-CVNRPN(K2PN)
  
  
  115 CONTINUE




      GO TO 47
   10 JC=JSS
      LAS2=LAS/2
      NCLA=NCLL*LAS2
      BTGS=BET(2)
      AVOL=1.D0
      DO 1 K=1,NCLL
      K1=(K-1)*NCLA
      K1PN=(K-1)*NCLL
      NU=NNJ1(K)
      KO1=KO(NU)
      JO1=JO(NU)
      NNB1=NNB(NU)
      LN1=LNJ1(K)
      J1=JNJ1(K)
      BT1=BETB(NU)
      NUMB1=NUMB(NU)
      NCA1=NCA(NU)
      IF(NUMB1.EQ.NUMB(1) )BT1=BET(2)
      ALFA1=1.D0
      IF(JO(1)/2*2.NE.JO(1)) ALFA1=AIGS(NU)
      JN1=J1
      JA=J1
      JB=JO1
      DO 2 KK=1,NCLL
      K2=K1+(KK-1)*LAS2
      K2PN=K1PN+KK
      IF(NPD.EQ.0) GO TO 2
      LN2=LNJ1(KK)
      J2=JNJ1(KK)
      JN2=J2
      JE=J2
      CC=SQRT((J1+1)*(J2+1)/12.566371D0)
      NU1=NNJ1(KK)
      NNB2=NNB(NU1)
      BT2=BETB(NU1)
      NUMB2=NUMB(NU1)
      NCA2=NCA(NU1)
      IF(NUMB2.EQ.NUMB(1) )BT2=BET(2)
      ALFA2=1.D0
      IF(JO(1)/2*2.NE.JO(1)) ALFA2=AIGS(NU1)
      B12BGS=SQRT(BT1*BT2/BTGS**2)
      B12BGS=BT1*BT2/BTGS**2
      IF(NUMB1.NE.NUMB2) B12BGS=B12BGS*ALFA1*ALFA2
      KO2=KO(NU1)                         
                                
      JO2=JO(NU1)

      KO1=KO(NU)
C     COUPLING DUE TO VOLUME CONSERVATION IN CASE OF QUADRUPOLE AND OCTUPOLE SHAPE OSCILLATIONS
C     GAMMA OSCILLATIONS DO NOT CHANGE NUCLEAR VOLUME

      BTGS2=BTGS**2
      IF(NNB1.NE.NNB2) BTGS2=0.D00
       CVNRV0(K2PN)=0.D0


C      PRINT "(3I3,6E11.3)", TID,NU,NU1,EFFDEF(NU,NU1,1:6)
      
C         PAUSE 222
      











 
C     IF(MEVOL.GE.1.AND.MEDEF.EQ.1.AND.LN1.EQ.LN2.
C    *AND.J1.EQ.J2.AND.KO1.EQ.KO2.AND.JO1.EQ.JO2.AND.NCA1.EQ.NCA2)  
C    *         CVNRV0(K2PN)=-(EFFDEF(NU,NU1,1)*BTGS*2.D0+BTGS2+                           
C    *          EFFDEF(NU,NU1,6)/BTGS)*0.282094791773878D0/AVOL
     
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
C       IF(MEVOL.GE.1.AND.MEDEF.EQ.1.AND.K.EQ.KK) 
       !IF(MEVOL.GE.1.AND.K.EQ.KK)  
      IF(MERAD.EQ.0.AND.NU.EQ.NU1) THEN
          STATCOR=BTGS!+BET3**2/BTGS!!! CHANGE TO BTGS !!!!
      ELSE
          STATCOR=0.0
      END IF 
      
      !IF(NU.EQ.NU1)THEN
      !    BET32AVG=BET3**2
      !ELSE
      !    BET32AVG=0.0
      !ENDIF
      
      
      IF(MEVOL.GE.1.!AND.LN1.EQ.LN2.
     *AND.J1.EQ.J2.AND.KO1.EQ.KO2.AND.JO1.EQ.JO2.AND.NCA1.EQ.NCA2)  
     *         CVNRV0(K2PN)=-(EFFDEF(NU,NU1,1)*BTGS*2.D0+STATCOR+!BTGS2+                           
     *          (EFFDEF(NU,NU1,6)!-BET32AVG
     * )/BTGS+EFFDEF(NU,NU1,7)*BTGS
     *          )*0.282094791773878D0/AVOL*0.282094791773878D0                               
     
      !!!if(NU.ne.1) CVNRV0(K2PN)=0
      ! IF(MEDEF.EQ.0) CVNRV0(K2PN)=0.0  ! temporary workout for odd


      JD=JO2
      JPS=JSS-1-JO2+J1+J2+LN2-LN1
      CCC=CC*(-1)**(JPS/2)
      
!      write (21,'(15A4)')'J','lam','N1','I1','J1','L1','pi1','K1',
!     *'N2','I2','J2','L2','pi2','K2', 'CVNR'
      
    8 DO 9 L=1,LAS2
    
      LL2=K2+L 

      M1=-1
      M2=1
      LAM=2*L
      IF(NPO(NU).NE.NPO(NU1))LAM=LAM+1

      
      J=LAM*2
      M=0
      CALL KLEGO
      BKG=AKG
      
C      SUMMING BY K AND K' BY WEIGHTS AIT
     
      CVNR(LL2)=0.D0
      NKI1=NKIWE(NU)
      NKI2=NKIWE(NU1)
      IF(MEAXI.EQ.0) THEN
          NKI1=1
          NKI2=1
          WEIGH(1,NU)=1.D0
          WEIGH(1,NU1)=1.D0
          KWE(1,NU)=KO1
          KWE(1,NU1)=KO2
      END IF
      DO 333 K1WE=1,NKI1
      KO1=KWE(K1WE,NU)
      WE1=WEIGH(K1WE,NU)
      DO 444 K2WE=1,NKI2
      KO2=KWE(K2WE,NU1)
      WE2=WEIGH(K2WE,NU1)
   
 
      IF(NUMB1.EQ.NUMB2) GO TO 125
        
      IF(L.GT.1) GO TO 140
      IF(KO1.EQ.KO2) GO TO 125
   
      IF(JO(1)/2*2.NE.JO(1)) GO TO 200

C     COUPLING FOR EVEN-EVEN NUCLIDES
C     NON-AXIAL BAND LEVELS
      !GO TO 145
      
C     LAST JOSES'S FORMULAS
      MUBAND=IABS(KO2-KO1)
      IF(MUBAND.GT.4) GO TO 140     

      
      J1=JO2
      J2=2*LAM
      J=JO1
      M1=KO2
      M2=MUBAND
      M=KO1
      CALL KLEGO
      AKGS=AKG
      M1=-KO2
      CALL KLEGO
      AKGS=AKGS+(-1)**(JO2/2)*NPO(NU1)*AKG
      M=-KO1
      CALL KLEGO
      AKGS=AKGS+(-1)**((JO1+JO2)/2)*NPO(NU)*NPO(NU1)*AKG 
      M1=KO2
      CALL KLEGO
      AKGS=AKGS+(-1)**(JO1/2)*NPO(NU)*AKG
      
      AKG=AKGS*SQRT(J1+1.D0)!/SQ2!2.D0

      IF(KO1.EQ.0.OR.KO2.EQ.0)THEN
          AKG=AKG/SQ2
      ENDIF
      
      GO TO 129
 
      
C     END OF NEW JOSE'S FORMULA FOR EVEN-EVEN CASE


      
      


  145 CONTINUE


C     OLD JOSE FORMULA  FOR EVEN-EVEN CASE
      IF(KO1.EQ.4.AND.KO2.EQ.0) GO TO 131
      IF(KO1.EQ.0.AND.KO2.EQ.4) GO TO 130
  131 J1=JO1
      J2=2*LAM
      J=JO2
      M1=4
      M2=-4
      M=0
      CALL KLEGO
      AKG=AKG*SQRT(J1+1.D0)*(-1.D0)**(J1/2)       
 
      IF(KO1.EQ.0.OR.KO2.EQ.0)THEN
          AKG=AKG/SQ2
      ENDIF      
      
      GO TO 129
      
         
      
  130 J1=JO2
      J2=2*LAM
      J=JO1
      M1=4
      M2=-4
      M=0
      CALL KLEGO
      AKG=AKG*SQRT(J1+1.D0)
         
      IF(KO1.EQ.0.OR.KO2.EQ.0)THEN
          AKG=AKG/SQ2
      ENDIF
      
      GO TO 129   
C     END JOSES'S FORMULAS FOR EVEN-EVEN CASE

  200 CONTINUE
C     COUPLING FOR ODD NUCLIDES, NON GS STATE LEVELS WITH GS ONES 

C      IF(KO1.EQ.KO(1).OR.KO2.EQ.KO(1))  GO TO 201
      IF(NUMB1.EQ.NUMB(1).OR. NUMB2.EQ.NUMB(1)) GO TO 201
      
C     GO TO 140
  201 CONTINUE
      !MUBAND=IABS(KO(1)-KO1)
      !IF(KO(1)-KO1.EQ.0) MUBAND=IABS(KO(1)-KO2)
C     IF(MUBAND.EQ.2.OR.MUBAND.GT.6) GO TO 140
      !MUBAND=4
      MUBAND=IABS(KO2-KO1)
      IF(MUBAND.GT.4) GO TO 140     
      !IF(MUBAND.EQ.0) PAUSE 11
      J1=JO2
      J2=2*LAM
      J=JO1
      M1=KO2
      M2=MUBAND
      M=KO1
      CALL KLEGO
      AKGS=AKG
      M1=-KO2
      CALL KLEGO
C     BUG RCN 19 June 2014
C     AKGS=AKGS+(-1)**((J02-1)/2)*AKG
      AKGS=AKGS+(-1)**((JO2-1)/2)*NPO(NU1)*AKG
      M=-KO1
      CALL KLEGO
      AKGS=AKGS+(-1)**((JO1+JO2)/2-1)*NPO(NU)*NPO(NU1)*AKG 
      M1=KO2
      CALL KLEGO
      AKGS=AKGS+(-1)**((JO1-1)/2)*NPO(NU)*AKG
      
      AKG=AKGS*SQRT(J1+1.D0)!/2.D0  CHECK 1+(-1)^(lambda)*pi*pi'=1!!!!
      GO TO 129
 
      
C     END OF JOSE'S FORMULA ODD CASE
 
      
 

C    AXIAL ROTATOR COUPLING, WITH NORMALIZATION ACCOUNTING 
C    EFFECTIVE DEFORMATIONS OF OTHER BANDS 
      
      
 125  J1=JO2
      J2=2*LAM
      M1=KO(NU1)
      M2=0
      M=KO1
      J=JO1
      CALL KLEGO
      AKG=AKG*SQRT(J1+1.D0)
 129  JF=LAM*2
      CALL RACAH
      CVNN1=CCC*BKG*AKG*W
      SCALE=B12BGS**L
      IF(MEDEF.EQ.1) SCALE=0.D0       !!!!!!!!!
        if(numb1.ne.numb(1).and.numb2.ne.numb(1)) SCALE=0.d0
      IF(NUMB1.EQ.NUMB2) SCALE=1.D0
      IF(MEDEF.EQ.0.AND.MEVOL.EQ.0) GO TO 777
c      IF(MEAXI.EQ.0) GO TO 777
c      IF(MEDEF.EQ.0) GO TO 777

      !!!IF(JO(1)/2*2.NE.JO(1)) GO TO 777
  
      
       IF(L.GT.1 ) GO TO 88 
      !! NUMB is correct quantum number instead of KO 
      IF(MEDEF.EQ.2.AND.NUMB1.EQ.NUMB2) SCALE=DCOS(GAM0)
      IF(MEDEF.EQ.2.AND.IABS(KO1-KO2).EQ.4.AND.NPO(NU).EQ.NPO(NU1))
     *       SCALE=DSIN(GAM0)/SQ2
      IF(MEDEF.EQ.2) GO TO 777
    
       !! replace effective deformations of K=2 and gamma bands to ones 
       !! of rigid asymmetric rotor since in SRM they are is not correct
       IF(MEDEF.EQ.1.AND.NUMB1.EQ.NUMB2) THEN
           IF(KO1.EQ.KO2)THEN
               SCALE=!1.D0 + EFFDEF(NU,NU1,8) +
!    *             (EFFDEF(NU,NU1,1)-EFFDEF(NU,NU1,2))/AVOL
     *             DCOS(GAM0) + EFFDEF(NU,NU1,8)+EFFDEF(NU,NU1,1)
           ELSEIF(IABS(KO1-KO2).EQ.4)THEN
               SCALE=DSIN(GAM0)/SQ2!EFFDEF(NU,NU1,3)
           ELSE    
               SCALE=0.0  
           ENDIF    
       ENDIF
           

     
       IF(MEDEF.EQ.0.AND.MEVOL.GE.1.AND.NUMB(NNJ1(K)).EQ.NUMB(1).AND.
     *       NUMB(NNJ1(KK)).EQ.NUMB(1)) SCALE=1.D0 + EFFDEF(NU,NU1,8)+
     *       (EFFDEF(NU,NU1,1)-EFFDEF(NU,NU1,2))/AVOL  
            
      IF(NUMB1.EQ.NUMB2) GO TO 777
      IF(L.GT.1 ) pause
       
      !! replace effective deformations of K=2 and gamma bands to ones 
      !! of rigid asymmetric rotor since in SRM they are is not correct
      IF (MEDEF.EQ.1.AND.NUMB1.NE.NUMB2) THEN
       IF (KO1.EQ.KO2.AND.NPO(NU).EQ.NPO(NU1)) THEN
             SCALE=!(EFFDEF(NU,NU1,1)-EFFDEF(NU,NU1,2))/AVOL             
     *         DCOS(GAM0)-1.0+EFFDEF(NU,NU1,1)
          ELSE IF (IABS(KO1-KO2).EQ.4.AND.NPO(NU).EQ.NPO(NU1)) THEN
             SCALE=DSIN(GAM0)/SQ2!EFFDEF(NU,NU1,3)!/BTGS!(1.D0+BTGS**2)/AVOL
          ELSE IF (KO1.EQ.KO2.AND.NPO(NU).NE.NPO(NU1)) THEN
             SCALE=EFFDEF(NU,NU1,4)/BTGS/AVOL
C             PRINT *, SCALE,EFFDEF(NU,NU1,4),BTGS,AVOL,WE1,WE2,CVNN1
C             PAUSE 4444
          ELSE IF (IABS(KO1-KO2).EQ.4.AND.NPO(NU).NE.NPO(NU1)) THEN
             SCALE=EFFDEF(NU,NU1,5)/BTGS/AVOL
          ELSE 
             SCALE=0.D0
          END IF   
      END IF
  
      !! compensation of CoM movement
      IF (MEVOL.EQ.2.AND.MEDEF.EQ.1!.AND.DIPOLECORRECTION.EQ.1.AND.EVENEVEN
     *   .AND.NUMB1.NE.NUMB2.AND.KO1.EQ.KO2.AND.NPO(NU).NE.NPO(NU1))THEN
          ! Part of A-factor
          J1=JNJ1(K)
          J2=JNJ1(KK)
          M1=-1
          M2=1
          J=1*2 !!2*lam
          M=0
          CALL KLEGO
          AKGS=AKG      
          ! RME
          J1=JO2
          J2=2*1!! 2*lam
          M1=KO(NU1)
          M2=0
          M=KO1
          J=JO1
          CALL KLEGO
          AKG=AKG*SQRT(J1+1.D0)
          ! Racah
          JF=1*2 !! 2*lam
          CALL RACAH
          ! Scale = -[beta3]eff*sqrt(3/(7pi))*(9/2/sqrt(5)beta2+2beta4)
          !     CVNN1=CCC*BKG*AKG*W          
          CVNR(LL2)=CVNR(LL2)-CCC*AKG*AKGS*W*WE1*WE2*EFFDEF(NU,NU1,4)*
     *                (0.7433*BTGS+0.7387*BET(4))
      END IF
      !! end CoM compensation

  777 CONTINUE    
C  88  CVNR(LL2)=CVNR(LL2)*SCALE
            
   88 CVNR(LL2)=CVNR(LL2)+CVNN1*SCALE*WE1*WE2
      



C     NECESSARY TO DIVIDE BY SIN(GAMMA)/SQRT(2.)!!! for actinides
            
C 144  IF(KO1.EQ.4.AND.KO2.EQ.4) CVNR(LL2)=CVNR(LL2)/0.1D0                 !!!!!!

      GO TO 141
  140 CVNR(LL2)=0.D0
  141 J1=JN1
      J2=JN2
      
!      write(21,'(14I4,E12.4)') JSS,LAM,
!     * NU,JO1,JN1,LN1,NPO(NU),KO(NU),
!     * NU1,JO2,JN2,LN2,NPO(NU1),KO(NU1), 
!     * CVNR(LL2)
      
  444 CONTINUE
  333 CONTINUE 
    9 CONTINUE
    
     
 
  124 IF(JO1.EQ.JO2.AND.NCA(NU).NE.NCA(NU1)) GO TO 135
      GO TO 2
  135 IF(NUMB1.EQ.NUMB(1). AND. NUMB2.EQ.NUMB(1)) GO TO 101
      GO TO 2

  101 CC=SQRT((J1+1)*(J2+1)/12.566371D0)
      NU1=NNJ1(KK)
      JO2=JO(NU1)
      JD=JO2
      JPS=JSS-1-JO2+J1+J2+LN2-LN1
      CCC=CC*(-1)**(JPS/2)
   21 M1=-1
      M2=1
      LAM=0
      J=0
      M=0
      CALL KLEGO
      BKG=AKG
      J1=JO2
      J2=2*LAM
      M1=KO(NU1)
      M2=0
      M=KO1
      J=JO1
      CALL KLEGO
      AKG=AKG*SQRT(J1+1.D0)
      JF=0
      CALL RACAH
      CVNRPN(K2PN)=CCC*BKG*AKG*W
       IF(MEDEF.EQ.1.AND.L.EQ.1)CVNRPN(K2PN)=CVNRPN(K2PN)*
     *  (1.D0+EFFDEF(NU,NU1,8)+(EFFDEF(NU,NU1,1)-EFFDEF(NU,NU1,2))/AVOL)
      J1=JN1
      J2=JN2
C     WRITE(21,'(1X,8I8,1X,6(f8.4,1X))') NCLL,K,KK,KL, JO1,JO2,NCA(NU)
C     *,NCA(NU1)
C     *,CVNRPN(K2PN)
 
    2 CONTINUE
    1 CONTINUE
   47 CONTINUE
      RETURN
      END
C     *****************************************************************
      SUBROUTINE PREQU
C     *****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 

C---------------------------------
C     These commons are also used and initialized in ABCT 
      INCLUDE 'PRIVCOM10.FOR'
      INCLUDE 'PRIVCOM12.FOR'
      INCLUDE 'PRIVCOM9.FOR'
      INCLUDE 'PRIVCOM13.FOR'
      INCLUDE 'PRIVCOM1.FOR'
      INCLUDE 'PRIVCOM.FOR'         
C     INCLUDE 'PRIVCOM15.FOR'
      INCLUDE 'PRIVCOM6.FOR'
C--------------------------------- 
C
      INCLUDE 'PRIVCOM8.FOR'

      SI=SIN(GAM0+GSHAPE)
      CO=COS(GAM0+GSHAPE)
      CG1=CO
      SG1=SI
      SG2=SI*SI
      CG2=CO*CO
      CSG11=SI*CO
      CG3=CG1*CG2
      SG3=SG1*SG2
      CSG21=CG2*SG1
      CSG12=SG2*CG1
      CG4=CG2*CG2
      SG4=SG2*SG2
      CSG22=CG2*SG2
      CSG31=CG3*SG1
      CSG13=SG3*CG1
      FOLAR=1.
      DO 1 I=1,NUR
    1 JU(I)=JO(I)/2
      CALL SHEM
      IF(MESHO.EQ.0) GO TO 4
      NOM=0
      DO 5 I=1,NUR
      IF(NOM.LT.NNO(I)) NOM=NNO(I)
    5 CONTINUE
      NOM1=NOM+1
      LA=1
      IF(LAS.GE.2) LA=2
      DO 6 IO1=1,NOM1
      ANU1=ANO(IO1)
      DO 6 IO2=IO1,NOM1
      ANU2=ANO(IO2)
      DO 6 NNT=1,LA
      FOLAR=BET3**NNT
      IF(MEHAO.EQ.0) GO TO 7
      IF(MEHAO.LT.2) GO TO 9
      FOLAR=0.
      EBM=EXP(-(BET3/AMUO)**2)
      IF(NNT.EQ.1.AND.IO1.NE.IO2) FOLAR=BET3/SQRT(1.D0-EBM**2)
      IF(NNT.EQ.2.AND.IO1.EQ.IO2) FOLAR=AMUO**2/2.D0+BET3**2/
     *(1.-(-1)**IO1*EBM)
      IF(MEHAO.GE.2) GO TO 7
    9 X1=-BET3/AMUO
      X2=X1
      P1=BET3
      P2=BET3
      IF(BET3.EQ.0.) CALL OVLAO
      IF(BET3.NE.0.) CALL OVLAB
    7 IF(BET3.NE.0.) FOLAR=FOLAR/BET3**NNT
      IF(BET3.EQ.0.) FOLAR=FOLAR/AMUO**NNT
      IF(NNT.EQ.1.AND.IO1.EQ.IO2) FOLAR=0.D0
      IF(NNT.EQ.2.AND.IO1.NE.IO2) FOLAR=0.D0
      FOK(IO1,IO2,NNT)=FOLAR
      IF(IO1.NE.IO2) FOK(IO2,IO1,NNT)=FOK(IO1,IO2,NNT)
      IF(MEPRI.LT.98) THEN
        PRINT 8,IO1,IO2,NNT,FOLAR,ANU1,ANU2
        WRITE(21,88)IO1,IO2,NNT,FOLAR,ANU1,ANU2
      ENDIF
    8 FORMAT (3I4,3E20.7)
   88 FORMAT (1X,'IO1=',I3,1X,'IO2=',I3,1X,'NNT=',I3,2X,'FOLAR=',D15.7,
     *1X,'ANU1=',D15.7,2X,'ANU2=',D15.7)
    6 CONTINUE
    4 DO 2 JU1=1,NUR
      KG1=(3-(-1)**NNO(JU1))/2
      DG1=0.
      IF(KG1.EQ.2) DG1=GAMDE
      ANG1=ANG(NNG(JU1)+1,KG1)
      CD1=CD(NNG(JU1)+1,KG1)
      X1=XI(JU1)
      P1=PT(JU1)
      ANU1=ANB(JU1)
      DO 2 JU2=JU1,NUR
      KG2=(3-(-1)**NNO(JU2))/2
      DG2=0.
      IF(KG2.EQ.2) DG2=GAMDE
      ANG2=ANG(NNG(JU2)+1,KG2)
      CD2=CD(NNG(JU2)+1,KG2)
      X2=XI(JU2)
      P2=PT(JU2)
      ANU2=ANB(JU2)
      IF(MEHAM.GT.4) CALL OVLAGE
      TRIG(JU1,JU2,1)=CG1
      IF(JU1.NE.JU2) TRIG(JU2,JU1,1)=CG1
      TRIG(JU1,JU2,2)=SG1
      IF(JU1.NE.JU2) TRIG(JU2,JU1,2)=SG1
      TRIG(JU1,JU2,3)=CG2
      IF(JU1.NE.JU2) TRIG(JU2,JU1,3)=CG2
      TRIG(JU1,JU2,4)=SG2
      IF(JU1.NE.JU2) TRIG(JU2,JU1,4)=SG2
      TRIG(JU1,JU2,5)=CSG11
      IF(JU1.NE.JU2) TRIG(JU2,JU1,5)=CSG11
      TRIG(JU1,JU2,6)=CG3
      IF(JU1.NE.JU2) TRIG(JU2,JU1,6)=CG3
      TRIG(JU1,JU2,7)=SG3
      IF(JU1.NE.JU2) TRIG(JU2,JU1,7)=SG3
      TRIG(JU1,JU2,8)=CSG21
      IF(JU1.NE.JU2) TRIG(JU2,JU1,8)=CSG21
      TRIG(JU1,JU2,9)=CSG12
      IF(JU1.NE.JU2) TRIG(JU2,JU1,9)=CSG12
      TRIG(JU1,JU2,10)=CG4
      IF(JU1.NE.JU2) TRIG(JU2,JU1,10)=CG4
      TRIG(JU1,JU2,11)=SG4
      IF(JU1.NE.JU2) TRIG(JU2,JU1,11)=SG4
      TRIG(JU1,JU2,12)=CSG22
      IF(JU1.NE.JU2) TRIG(JU2,JU1,12)=CSG22
      TRIG(JU1,JU2,13)=CSG31
      IF(JU1.NE.JU2) TRIG(JU2,JU1,13)=CSG31
      TRIG(JU1,JU2,14)=CSG13
      IF(JU1.NE.JU2) TRIG(JU2,JU1,14)=CSG13
      FOLAR=1.
      DO 3 NNT=1,4     !MAY BE MORE THAN LAS
      IF(MEHAM.NE.4) CALL OVLAB
      FOV(JU1,JU2,NNT)=FOLAR
      IF(JU1.NE.JU2) FOV(JU2,JU1,NNT)=FOV(JU1,JU2,NNT)
      IF(MEPRI.LT.98) THEN
        PRINT 10,JU1,JU2,NNT,FOV(JU1,JU2,NNT),ANU1,ANU2
        WRITE(21,89)JU1,JU2,NNT,FOV(JU1,JU2,NNT),ANU1,ANU2
      ENDIF
   10 FORMAT(3I5,3E17.7)
   89 FORMAT(1X,'JU1=',I2,1X,'JU2=',I2,1X,'NNT=',I2,1X,'FOV(JU1,JU2,
     *NNT)=',D15.7,2X,'ANU1=',D15.7,2X,'ANU2=',D15.7)
    3 CONTINUE
    2 CONTINUE
      RETURN
      END
C     *******************************************************
C     END of knditd
C     *******************************************************
