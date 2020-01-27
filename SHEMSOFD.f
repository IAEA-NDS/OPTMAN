C     *******************************************************
C     START of shemsofd
C     *******************************************************
C     *****************************************************************
      SUBROUTINE FUDNU
C     *****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM9.FOR'
      SQ2=SQRT(2.d0)
      YY=Y
      Z2=Y**2
      AN2=ANU/2.D0
      EX4=DEXP(-Z2/4.D0)
      GNP1=DGAMMA(ANU+1.D0)
      M=1
       PI=4.D0*DATAN(1.D0)
       PIN=PI*ANU
       PIN2=PI*AN2
      SQPIN=0.564189583548D0
      IF(DABS(YY).GT.5.AND.DABS(YY).GT.ANU) GO TO 1
      CPN2=DCOS(PIN2)
      SPN2=DSIN(PIN2)
      GN2P12=DGAMMA(AN2+0.5D0)
      GN2P1=DGAMMA(AN2+1.D0)
      C2N2=(2.D0)**AN2
      FM=1.D0
      VM=1.D0
      DF=FM
      DV=VM
    2 Z2M=Z2/M
      FM=FM*Z2M/(2*M-1)*(M-1-AN2)
      VM=VM*Z2M/(2*M+1)*(M-0.5D0-AN2)
      DF=DF+FM
      DV=DV+VM
      IF((DABS(FM)+DABS(VM))/(DABS(DF)+DABS(DV)).LT.1.D-6)GO TO 5
      M=M+1
      GO TO 2
    5 DN=C2N2*SQPIN*EX4*(CPN2*GN2P12*DF+SQ2*Y*SPN2*GN2P1*DV)
      DNV=EX4/C2N2*(-SPN2/GN2P1*DF+SQ2*Y*CPN2/GN2P12*DV)
      GO TO 6
    1 YN=Z2**AN2
      YN1=Z2**(AN2+0.5D0)
      CPN=DCOS(PIN)
      SPN=DSIN(PIN)
      SQP2=SQ2*SQPIN
      FM=1.D0
      DF=FM
    4 Z2M=Z2*(2*M)
      VM=(2*M+ANU)*(2*M-1+ANU)/Z2M
      IF(DABS(VM).GT.1.D0)GO TO 3
      FM=FM*VM
      DF=DF+FM
      IF(DABS(FM).LT.1.D-6)GO TO 3
      M=M+1
      GO TO 4
    3 M=1
      FM=1
      DV=FM
    9 Z2M=Z2*(2*M)
      VM=-(ANU-2*M+1)*(ANU-2*M+2)/Z2M
      IF(DABS(VM).GT.1.D0)GO TO 8
      FM=FM*VM
      DV=DV+FM
      IF(DABS(FM).LT.1.D-6)GO TO 8
      M=M+1
      GO TO 9
    8 IF(Y) 10,10,11
   11 DN=EX4*YN*DV
      DNV=SQP2/EX4/YN1*DF
      GO TO 6
   10 DN=-SQP2*SPN*GNP1/EX4/YN1*DF+CPN*EX4*YN*DV
      DNV=-SPN*EX4/GNP1*YN*DV-CPN*SQP2/EX4/YN1*DF
    6 RETURN
      END
C     *****************************************************************
      SUBROUTINE MATAM
C     *****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM9.FOR'
      KMI=1-(-1)**(IS+NO)
      KMA=2*(IS/2)
      NK=(KMA-KMI+2)/2
      CALL INERMO
      A=AMO
      B=BMO
      C=CMO
      C4=C/4.D0
      AB8=(A+B)/8.
      AB16=(A-B)/16.
      DO 1 J=1,NK
      DO 1 I=1,J
      KS=KMI+2*(I-1)
      IF(I.EQ.J) GO TO 2
      IF(I+1.EQ.J) GO TO 3
    5 AM(I,J)=0.D0
      GO TO 1
    2 AM(I,J)=AB8*(IS*(IS+1.D0)-KS*KS)+C4*KS*KS
      GO TO 1
    3 AB16IK=AB16*DSQRT((IS+KS+2.D0)*(IS-KS-1.D0)*(IS+KS+1.D0)*(IS-KS))
      IF(KS.EQ.0) GO TO 4
      AM(I,J)=AB16IK
      GO TO 1
    4 AM(I,J)=AB16IK*(1+(-1)**(IS+NO))/DSQRT(2.d0)
    1 CONTINUE
      CALL VECNO
      RETURN
      END
C     *****************************************************************
      SUBROUTINE VECNO
C     *****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM9.FOR'
      DO 1 I=1,NK
      DO 1 J=1,NK
      IF(I.EQ.J) GO TO 2
      TM(I,J)=0.D0
      TM1(I,J)=0.D0
      GO TO 18
    2 TM(I,J)=1.D0
      TM1(I,J)=1.D0
   18 IF(I.GT.J) GO TO 1
      AM1(I,J)=AM(I,J)
    1 CONTINUE
      IF(NK.EQ.1) GO TO 17
    3 IM=1
      JM=2
      A=DABS(AM(1,2))
      IF(NK.EQ.2) GO TO 20
      DO 4 J=3,NK
      J1=J-1
      DO 4 I=1,J1
      B=DABS(AM(I,J))
      IF(A.GE.B) GO TO 4
      A=B
      IM=I
      JM=J
    4 CONTINUE
   20 CONTINUE
      AMII=AM(IM,IM)
      AMJJ=AM(JM,JM)
      AMIJ=AM(IM,JM)
      AL=DATAN(2.D0*AMIJ/(AMII-AMJJ))/2.D0
      ACO=DCOS(AL)
      ASI=DSIN(AL)
      DO 5 L=1,NK
      DO 5 K=1,L
      IF(K.EQ.IM) GO TO 6
      IF(K.EQ.JM) GO TO 7
      IF(L.EQ.IM) GO TO 8
      IF(L.EQ.JM) GO TO 9
      GO TO 5
    6 IF(L.EQ.IM) GO TO 10
      IF(L.EQ.JM) GO TO 11
      AM1(K,L)=ACO*AM(IM,L)+ASI*AM(JM,L)
      GO TO 5
    7 IF(L.EQ.JM) GO TO 12
      AM1(K,L)=-ASI*AM(IM,L)+ACO*AM(JM,L)
      GO TO 5
    8 AM1(K,L)=ACO*AM(K,IM)+ASI*AM(K,JM)
      GO TO 5
    9 AM1(K,L)=-ASI*AM(K,IM)+ACO*AM(K,JM)
      GO TO 5
   10 AM1(K,L)=ACO*(ACO*AMII+ASI*AMIJ)+ASI*(ACO*AMIJ+ASI*AMJJ)
      GO TO 5
   11 AM1(K,L)=0.D0
      GO TO 5
   12 AM1(K,L)=-ASI*(-ASI*AMII+ACO*AMIJ)+ACO*(-ASI*AMIJ+ACO*AMJJ)
    5 CONTINUE
      DO 13 K=1,NK
      DO 13 L=1,NK
      IF(K.EQ.IM) GO TO 14
      IF(K.EQ.JM) GO TO 15
      GO TO 13
   14 TM1(K,L)=ACO*TM(IM,L)-ASI*TM(JM,L)
      GO TO 13
   15 TM1(K,L)=ASI*TM(IM,L)+ACO*TM(JM,L)
   13 CONTINUE
      E2A=0.D0
      DO 16 K=1,NK
      DO 16 L=1,NK
      TM(K,L)=TM1(K,L)
      IF(K.GT.L) GO TO 16
      IF(K.LT.L) E2A=E2A+AM1(K,L)**2
      AM(K,L)=AM1(K,L)
   16 CONTINUE
      IF(E2A.LT.1.D-8) GO TO 17
      GO TO 3
   17 CONTINUE
      IF(IS.EQ.0) GO TO 25
      DO 24 I=1,NK
      I1=I+1
      EM=AM(I,I)
      IM=I
      IF(I1.GT.NK) GO TO 23
      DO 21 J=I1,NK
      E=AM(J,J)
      IF(EM.LE.E) GO TO 21
      EM=E
      IM=J
   21 CONTINUE
   23 DO 22 J=1,NK
      T1=TM(J,IM)
      TM(J,IM)=TM(J,I)
      TM(J,I)=T1
   22 CONTINUE
      EIN(IS+1,I)=EM
      AM(IM,IM)=AM(I,I)
      AM(I,I)=EM
   24 CONTINUE
   25 CONTINUE
      IF(IS.EQ.0) EIN(1,1)=0.D0
      IF(IS.EQ.0) TM(1,1)=1.D0
      RETURN
      END

C     *****************************************************************
      BLOCK DATA
C     *****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM11.FOR'
      DATA B0/-22.11505d0,6.64844d0,.769968d0,-.622406d0,.107521d0,
     * -.007364d0/, B1/-17.87240d0,4.66396d0,2.94376d0,-2.03488d0,
     * .51202d0,-.051336d0/, C0/-2.32285d0,2.65365d0,-1.077629d0,
     * .173691d0,-.007545d0/, C1/-.630875d0,1.198105d0,-.727493d0,
     * .164266d0,-.0091289d0/, AN/-2.d0,-2.2d0,-2.4d0,-2.6d0/
      DATA B2/-15.007288d0,5.405502d0,1.117805d0,-.831844d0,.110628d0/,
     *B3/-12.206497d0,4.48481d0,1.679626d0,-1.226475d0,.190465d0/,
     *C2/.82866d0,-.587769d0,-.005974d0,.0581198d0,-.00321186d0/,
     *C3/1.191310d0,-1.284630d0,.368034d0,-.00188923d0,.000426352d0/
      END
C     ***************************************************************
      SUBROUTINE EIT12
C     *****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM9.FOR'
      INCLUDE 'PRIVCOM10.FOR'
      INCLUDE 'PRIVCOM12.FOR'
      C1=DSQRT(2.D0)
      PI3=4.D0/3.D0*DATAN(1.D0)      
      IF (MEHAM.EQ.5.OR.MEHAM.EQ.7) GO TO 2
      AMOO=AMO
      BMOO=BMO
      CMOO=CMO
      DEGAM=0.001D0*GAM
      GAM=GAM0+DEGAM
      CALL INERMO
      AMOP=AMO
      BMOP=BMO
      CMOP=CMO
      GAM=GAM0-DEGAM
      CALL INERMO
      AMOM=AMO
      BMOM=BMO
      CMOM=CMO
      GAM=GAM0
      A1D=(AMOP-AMOM)/DEGAM/2.D0
      B1D=(BMOP-BMOM)/DEGAM/2.D0
      C1D=(CMOP-CMOM)/DEGAM/2.D0
      A2D=(AMOP+AMOM-2.D0*AMOO)/DEGAM/DEGAM
      B2D=(BMOP+BMOM-2.D0*BMOO)/DEGAM/DEGAM
      C2D=(CMOP+CMOM-2.D0*CMOO)/DEGAM/DEGAM
      AB1DM=A1D-B1D
      AB1DP=A1D+B1D
      AB2DM=A2D-B2D
      AB2DP=A2D+B2D
      NGD=2
      IF(MEHAM.EQ.5.OR.MEHAM.EQ.7) NGD=0
      SUM3=0.D0
      NGU=NROOT
      NR1=NROOT+1
      INGP=NGU+NGD+1
      IF(INGP.GT.4) INGP=4
      INGM=NGU-NGD+1
      IF(INGM.LT.1) INGM=1
      ANG1=ANG(NROOT+1,(3-(-1)**NO)/2)
      CD1=CD(NROOT+1,(3-(-1)**NO)/2)
      DO 10 INGA=1,3!INGM,INGP
      NGA=INGA-1
      ANG2=ANG(INGA,(3-(-1)**NO)/2)
      CD2=CD(INGA,(3-(-1)**NO)/2)
      DO 10 ITA=1,NK
      SUM1=0.D0
      SUM2=0.D0
      DO 1 I=1,NK
      AK=TM(I,ITAU)
      AKA=TM(I,ITA)
      KA=KMI+(I-1)*2
      SUM1=SUM1+AK*AKA*KA*KA
      D=1.D0
      DD=0.D0
      IF(KA.EQ.0) DD=1.D0
      IF(KA.EQ.0) D=DSQRT(2.D0)
      IF(I.EQ.NK) GO TO 1
      AK1=TM(I+1,ITAU)
      AKA1=TM(I+1,ITA)
      SUM2=SUM2+(AK1*AKA+AK*AKA1)*(1.D0+(-1)**IS*DD)/D*
     *SQRT((IS+KA+2.D0)*(IS+KA+1.D0)*(IS-KA-1.D0)*(IS-KA))
    1 CONTINUE
      KA=KMI+(ITAU-1)*2
      KAITA=KMI+(ITA-1)*2
      IF(ITA.EQ.ITAU.AND.NGU.EQ.NGA)GO TO 5
      IF(ITA.NE.ITAU.OR.NGU.NE.NGA)GO TO 4
      GO TO 10
    5 EPIT1=AB1DP/8.D0*IS*(IS+1)+(0.25D0*C1D-AB1DP/8.D0)*SUM1+
     *AB1DM/16.D0*SUM2
      EPIT1=EPIT1*FOG(1,NR1,INGA)
      EPIT2=AB2DP/8.D0*IS*(IS+1)+(0.25D0*C2D-AB2DP/8.D0)*SUM1+
     *AB2DM/16.D0*SUM2
      EPIT2=EPIT2*FOG(2,NR1,INGA)/2.
      GO TO 10
    4 SU1M=(0.25*C1D-AB1DP/8.)*SUM1+!*(1+(-1)**(KA/2+KAITA/2))/2+ !symmetric if K+K1=0,4,8...
     *AB1DM/16.*SUM2!*(1-(-1)**(KA/2+KAITA/2))/2                 !antisymmetric if K+K1=2,6,10... 
      IF(ITAU.EQ.ITA) SU1M=SU1M+AB1DP/8.D0*IS*(IS+1.D0)
      SU1M=SU1M*FOG(1,NR1,INGA)
      SU2M=(0.25*C2D-AB2DP/8.)*SUM1+!*(1+(-1)**(KA/2+KAITA/2))/2+
     *AB2DM/16.*SUM2!*(1-(-1)**(KA/2+KAITA/2))/2
      IF(ITAU.EQ.ITA) SU2M=SU2M+AB2DP/8.*IS*(IS+1.)
      SU2M=SU2M*FOG(2,NR1,INGA)/2.D0
      DELE=2.D0/AMG0**2*(ANG1-ANG2)+EIN(IS+1,ITAU)-EIN(IS+1,ITA)
      IF (NR1.EQ.INGA)  SUM3=SUM3+(SU1M+SU2M)**2/DELE 
   10 CONTINUE
      EPIT12=SUM3
      GO TO 3
    2 EPIT1=0.D0
      EPIT2=0.D0
      EPIT12=0.D0
    3 RETURN
      END
C     *****************************************************************
      SUBROUTINE ANUDF
C     *****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM9.FOR'
      INCLUDE 'PRIVCOM11.FOR'
      ANU=NROOT
      NRP1=NROOT+1
      IF(XIT+5) 7,7,8
    7 IF(XIT.LT.-10.D0) ANU=ANU
      IF(XIT.LT.-10.D0)GO TO 9
      Z2=2.D0*XIT*XIT
      ANU=ANU+DEXP(-Z2/2.D0)*Z2**(ANU+0.5D0)*0.39894228D0
     **(1.D0-(ANU*ANU+ANU+1.D0)/Z2)/DGAMMA(1.D0+ANU)
      GO TO 9
    8 Y=5.D0+XIT
      Y2=Y*Y
      Y3=Y*Y2
      Y4=Y*Y3
      IF(XIT-AN(NRP1)) 10,10,11
   10 Y5=Y*Y4
      GO TO (2,12,22,32,35,35,35),NRP1
    2 ANU=ANU+Y*DEXP(B0(1)+B0(2)*Y+B0(3)*Y2+B0(4)*Y3+
     *B0(5)*Y4+B0(6)*Y5)
      GO TO 9
   12 ANU=ANU+ Y*DEXP(B1(1)+B1(2)*Y+B1(3)*Y2+B1(4)*Y3+
     *B1(5)*Y4+B1(6)*Y5)
      GO TO 9
   22 ANU=ANU+Y*DEXP(B2(1)+B2(2)*Y+B2(3)*Y2+B2(4)*Y3+
     *B2(5)*Y4)
      GO TO 9
   32 ANU=ANU+Y*DEXP(B3(1)+B3(2)*Y+B3(3)*Y2+B3(4)*Y3+
     *B3(5)*Y4)
      GO TO 9
   11 GO TO (3,13,23,33,35,35,35),NRP1
    3 ANU=ANU+C0(1)+C0(2)*Y+C0(3)*Y2+C0(4)*Y3+C0(5)*Y4
      GO TO 9
   13 ANU=ANU+C1(1)+C1(2)*Y+C1(3)*Y2+C1(4)*Y3+C1(5)*Y4
      GO TO 9
   23 ANU=ANU+C2(1)+C2(2)*Y+C2(3)*Y2+C2(4)*Y3+C2(5)*Y4
      GO TO 9
   33 ANU=ANU+C3(1)+C3(2)*Y+C3(3)*Y2+C3(4)*Y3+C3(5)*Y4
      GO TO 9
   35 PRINT 36
      WRITE(21,36)
   36 FORMAT(10X,'THIS CASE IS NOT REALIZED NROOT>=4')
    9 RETURN
      END
C     *****************************************************************
      SUBROUTINE DETX12
C     *****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM9.FOR'
      Y=XIT
      CALL FUDNU
      UX1=DN
      VX1=DNV
      Y=XIT1
      CALL FUDNU
      UX2=DN
      VX2=DNV
      DET=UX1*VX2-UX2*VX1
      RETURN
      END
C     *****************************************************************
      SUBROUTINE ANDET0
C     *****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM9.FOR'
      ANU=0.D0
      NR=0
    5 STEP=0.1D0
      CALL DETX12
      A=DET
C     IF(A.EQ.0.) GO TO 2
      IF(DABS(A).LE.1.d-10) GO TO 2
    3 ANU=ANU+STEP
      CALL DETX12
      IF(A*DET)1,2,3
    1 ANU=ANU-STEP
      STEP=STEP/5.D0
      IF(STEP.LT.1.D-6) GO TO 2
      GO TO 3
    2 NR=NR+1
      IF(NR.GT.NROOT) GO TO 4
      ANU=ANU+0.01D0
      GO TO 5
    4 CDV=-DN/DNV
      RETURN
      END
C     *****************************************************************
      SUBROUTINE OVLAG
C     *****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM9.FOR'
      INCLUDE 'PRIVCOM12.FOR'
      C1=DSQRT(2.D0)
      G3=3.D0*GAM0
      IF(NNTG.EQ.11) CBG2=-27.D0/4.D0*AMG0**3*COS(G3)/SIN(G3)**3*C1
      YGH=-GSHAPE/AMG0*C1

      PI3=4.D0/3.D0*DATAN(1.D0)      
      NSTEP=16.D0*(ANG1+ANG2+4.D0)/AMG0
      STEP=PI3/AMG0/NSTEP*C1
      SUM2=0.D0
      S1M2=0.D0
      S2M2=0.D0
      SUM1=0.D0
      S1M1=0.D0
      S2M1=0.D0
      Y0=-GAM0/AMG0*C1
      NS2=NSTEP-1
      DO 1 I=1,NS2
      ANU=ANG1
      Y=Y0+STEP*I
      X=(Y+YGH)**NNTG
      IF(NNTG.EQ.11)X=(Y+YGH+CBG2)**2
      IF(NNTG.EQ.15)X=DSIN(AMG0/C1*Y+GAM0-GSHAPE)
      IF(NNTG.EQ.16)X=(1-DCOS(AMG0/C1*Y+GAM0-GSHAPE))
      IF(NNTG.EQ.17) THEN
            GAM=AMG0/C1*Y+GAM0-GSHAPE
            CALL INERMO !DSIN(AMG0/C1*Y+GAM0-GSHAPE)
            X=CMO
      END IF
      CALL FUDNU
      DN13=DN+CD1*DNV
      ANU=ANG2
      CALL FUDNU
      DN23=DN+CD2*DNV
      IF(I/2*2.EQ.I) GO TO 3
      SUM1=SUM1+X*DN13*DN23
      S1M1=S1M1+DN13**2
      S2M1=S2M1+DN23**2
      GO TO 1
    3 SUM2=SUM2+X*DN13*DN23
      S1M2=S1M2+DN13**2
      S2M2=S2M2+DN23**2
    1 CONTINUE
      FNR1=S1M2*2.D0+S1M1*4.D0
      FNR2=S2M2*2.D0+S2M1*4.D0
      FOLAG=SUM2*2.D0+SUM1*4.D0
      FOLAG=FOLAG/DSQRT(FNR1*FNR2)*(AMG0/C1)**NNTG
      IF(NNTG.EQ.11) FOLAG=FOLAG/(AMG0/C1)**NNTG
      IF(NNTG.EQ.15) FOLAG=FOLAG/(AMG0/C1)**(NNTG)
      IF(NNTG.EQ.16) FOLAG=FOLAG/(AMG0/C1)**(NNTG)
      IF(NNTG.EQ.17) FOLAG=FOLAG/(AMG0/C1)**(NNTG)
      RETURN
      END
C     *****************************************************************
      SUBROUTINE KLEGO
C     *****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM1.FOR'
      INCLUDE 'PRIVCOM4.FOR'      
      IF(J1.LT.IABS(M1)) GO TO 1
      IF(J2.LT.IABS(M2)) GO TO 1
      IF(J.LT.IABS(M)) GO TO 1
      IF(M1+M2-M)1,2,1
    1 AKG=0.D0
      GO TO 14
    2 IF(IABS(J1-J2)-J)3,3,1
    3 IF(J1+J2-J)1,4,4
    4 NF=J1+J2-J
      IF(NF/2*2.NE.NF) GO TO 1
      K=J1+J2+J
      IF(M1.EQ.0.AND.M.EQ.0.AND.K/4*4.NE.K) GO TO 1
      FLN=A(NF+2)
      CL=FLN
      NF=J1-J2+J
      FLN=A(NF+2)
      CL=CL+FLN
      NF=J+J2-J1
      FLN=A(NF+2)
      CL=CL+FLN
      NF=J1+J2+J+2
      FLN=A(NF+2)
      CL=CL-FLN
      NF=J1+M1
      FLN=A(NF+2)
      CL=CL+FLN
      NF=J1-M1
      FLN=A(NF+2)
      CL=CL+FLN
      NF=J2+M2
      FLN=A(NF+2)
      CL=CL+FLN
      NF=J2-M2
      FLN=A(NF+2)
      CL=CL+FLN
      NF=J+M
      FLN=A(NF+2)
      CL=CL+FLN
      NF=J-M
      FLN=A(NF+2)
      CL=0.5*(CL+FLN)
      NF1=J1+J2-J
      NF2=J1-M1
      NF3=J2+M2
      NF4=J-J2+M1
      NF5=J-J1-M2
      NB=NF1
      NM=-NF4
      IF(NF2-NB)5,6,6
    5 NB=NF2
    6 IF(NF3-NB)7,8,8
    7 NB=NF3
    8 IF(-NF5.GT.NM) NM=-NF5
      IF(NM.LT.0) NM=0
      NM=NM+2
      NB=NB+2
      AKG=0.D0
      IF(NB.LT.NM) GO TO 14
      DO 13 I1=NM,NB,2
      C1=1.D0
      N=I1-2
      NF=N
      FLN=A(NF+2)
      CL1=CL-FLN
      NF=NF1-N
      FLN=A(NF+2)
      CL1=CL1-FLN
      NF=NF2-N
      FLN=A(NF+2)
      CL1=CL1-FLN
      NF=NF3-N
      FLN=A(NF+2)
      CL1=CL1-FLN
      NF=NF4+N
      FLN=A(NF+2)
      CL1=CL1-FLN
      NF=NF5+N
      FLN=A(NF+2)
      CL1=CL1-FLN
      IF(N/4*4.NE.N) C1=-1.D0
   13 AKG=AKG+C1*DEXP(CL1)
      AKG=AKG*DSQRT(J+1.D0)
   14 RETURN
      END
C     *****************************************************************
      SUBROUTINE TRLAG
C     *****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM9.FOR'
      INCLUDE 'PRIVCOM12.FOR'      
      PI3=4.D0/3.D0*DATAN(1.D0)    
      C1=DSQRT(2.D0)
      NSTEP=4.D0*(ANG1+ANG2+4.D0)
      STEP=PI3/AMG0/NSTEP*C1
      SCM2=0.D0
      SSM2=0.D0
      S1M2=0.D0
      S2M2=0.D0
      SCM1=0.D0
      SSM1=0.D0
      S1M1=0.D0
      S2M1=0.D0
      SCM22=0.D0
      SSM22=0.D0
      SCM12=0.D0
      SSM12=0.D0
      SCM13=0.D0
      SCM23=0.D0
      Y0=-C1*GAM0/AMG0
      DY1=-C1*DG1/AMG0
      DY2=-C1*DG2/AMG0
      NS2=NSTEP-1
      DO 1 I=1,NS2
      ANU=ANG1
      Y=Y0+STEP*I
      SI=DSIN(GAM0+AMG0*Y/C1)
      CO=DCOS(GAM0+AMG0*Y/C1)
      SI2=DSIN(2*(GAM0+AMG0*Y/C1))
      CO2=DCOS(2*(GAM0+AMG0*Y/C1))
      CO3=DCOS(3*(GAM0+AMG0*Y/C1))
      Y=Y0+STEP*I+DY1
      CALL FUDNU
      DN13=DN+CD1*DNV
      ANU=ANG2
      Y=Y0+STEP*I+DY2
      CALL FUDNU
      DN23=DN+CD2*DNV
      IF(I/2*2.EQ.I) GO TO 3
      SSM1=SSM1+SI*DN13*DN23
      SCM1=SCM1+CO*DN13*DN23
      SSM12=SSM12+SI2*DN13*DN23
      SCM13=SCM13+CO3*DN13*DN23
      SCM12=SCM12+CO2*DN13*DN23
      S1M1=S1M1+DN13**2
      S2M1=S2M1+DN23**2
      GO TO 1
    3 SSM2=SSM2+SI*DN13*DN23
      SCM2=SCM2+CO*DN13*DN23
      SSM22=SSM22+SI2*DN13*DN23
      SCM22=SCM22+CO2*DN13*DN23
      SCM23=SCM23+CO3*DN13*DN23
      S1M2=S1M2+DN13**2
      S2M2=S2M2+DN23**2
    1 CONTINUE
      FNR1=S1M2*2.D0+S1M1*4.D0
      FNR2=S2M2*2.D0+S2M1*4.D0
      FOLAS=SSM2*2.D0+SSM1*4.D0
      FOLAC=SCM2*2.D0+SCM1*4.D0
      FOLS2=SSM22*2.D0+SSM12*4.D0
      FOLC2=SCM22*2.D0+SCM12*4.D0
      FOLC3=SCM23*2.D0+SCM13*4.D0
      DF=DSQRT(FNR1*FNR2)
      FOLAS=FOLAS/DF
      FOLAC=FOLAC/DF
      FOLS2=FOLS2/DF
      FOLC2=FOLC2/DF
      FOLC3=FOLC3/DF
      RETURN
      END
C     *****************************************************************
      SUBROUTINE OVLAGE
C     *****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INCLUDE 'PRIVCOM9.FOR'
      INCLUDE 'PRIVCOM12.FOR'
      INCLUDE 'PRIVCOM13.FOR'    
c     COMMON/ENA/EN,EL(40),BET(10),NUR,NMAX,NPD,LAS
      

      PI3=4.D0/3.D0*DATAN(1.D0)      
      C1=DSQRT(2.D0)
      NSTEP=4.*(ANG1+ANG2+4.)
      STEP=PI3/AMG0/NSTEP*C1
      SUF11=0.D0
      SUF12=0.D0
      SUF21=0.D0
      SUF22=0.D0
      C11=0.D0
      S11=0.D0
      C12=0.D0
      S12=0.D0
      C21=0.D0
      C22=0.D0
      S21=0.D0
      S22=0.D0
      CS111=0.D0
      CS112=0.D0
      C31=0.D0
      C32=0.D0
      S31=0.D0
      S32=0.D0
      CS121=0.D0
      CS122=0.D0
      CS211=0.D0
      CS212=0.D0
      C41=0.D0
      C42=0.D0
      S41=0.D0
      S42=0.D0
      CS221=0.D0
      CS222=0.D0
      CS311=0.D0
      CS312=0.D0
      CS131=0.D0
      CS132=0.D0
      Y0=-C1*GAM0/AMG0
      DY1=-C1*DG1/AMG0
      DY2=-C1*DG2/AMG0
      NS2=NSTEP-1
      DO 1 I=1,NS2
      Y=Y0+STEP*I
      SI=DSIN(GAM0+GSHAPE+AMG0*Y/C1)
      CO=DCOS(GAM0+GSHAPE+AMG0*Y/C1)
      ANU=ANG1
      Y=Y0+STEP*I+DY1
      CALL FUDNU
      DN1 =DN+CD1*DNV
      ANU=ANG2
      Y=Y0+STEP*I+DY2
      CALL FUDNU
      DN2 =DN+CD2*DNV
      DN12=DN1*DN2
      IF(I/2*2.EQ.I) GO TO 3
      SUF11=SUF11+DN1*DN1
      SUF21=SUF21+DN2*DN2
      CODN=CO*DN12
      SIDN=SI*DN12
      C11=C11+CODN
      S11=S11+SIDN
      IF(LAS.EQ.1) GO TO 1
      CO2DN=CO*CODN
      SI2DN=SI*SIDN
      C21=C21+CO2DN
      S21=S21+SI2DN
      CS111=CS111+CO*SIDN
      IF(LAS.EQ.2) GO TO 1
      CO3DN=CO*CO2DN
      SI3DN=SI*SI2DN
      C31=C31+CO3DN
      S31=S31+SI3DN
      CS121=CS121+CO*SI2DN
      CS211=CS211+SI*CO2DN
      IF(LAS.EQ.3) GO TO 1
      C41=C41+CO*CO3DN
      S41=S41+SI*SI3DN
      CS221=CS221+CO2DN*SI2DN/DN12
      CS311=CS311+SI*CO3DN
      CS131=CS131+CO*SI3DN
      GO TO 1
    3 SUF12=SUF12+DN1*DN1
      SUF22=SUF22+DN2*DN2
      CODN=CO*DN12
      SIDN=SI*DN12
      C12=C12+CODN
      S12=S12+SIDN
      IF(LAS.EQ.1) GO TO 1
      CO2DN=CO*CODN
      SI2DN=SI*SIDN
      C22=C22+CO2DN
      S22=S22+SI2DN
      CS112=CS112+CO*SIDN
      IF(LAS.EQ.2) GO TO 1
      CO3DN=CO*CO2DN
      SI3DN=SI*SI2DN
      C32=C32+CO3DN
      S32=S32+SI3DN
      CS122=CS122+CO*SI2DN
      CS212=CS212+SI*CO2DN
      IF(LAS.EQ.3) GO TO 1
      C42=C42+CO*CO3DN
      S42=S42+SI*SI3DN
      CS222=CS222+CO2DN*SI2DN/DN12
      CS312=CS312+SI*CO3DN
      CS132=CS132+CO*SI3DN
    1 CONTINUE
      SUF1=SUF11*4.D0+SUF12*2.D0
      SUF2=SUF21*4.D0+SUF22*2.D0
      ANOR=DSQRT(SUF1*SUF2)
      CG1=(C11*4.D0+C12*2.D0)/ANOR
      SG1=(S11*4.D0+S12*2.D0)/ANOR
      IF(LAS.EQ.1) GO TO 4
      CG2=(C21*4.D0+C22*2.D0)/ANOR
      SG2=(S21*4.D0+S22*2.D0)/ANOR
      CSG11=(CS111*4.D0+CS112*2.D0)/ANOR
      IF(LAS.EQ.2) GO TO 4
      CG3=(C31*4.D0+C32*2.D0)/ANOR
      SG3=(S31*4.D0+S32*2.D0)/ANOR
      CSG21=(CS211*4.D0+CS212*2.D0)/ANOR
      CSG12=(CS121*4.D0+CS122*2.D0)/ANOR
      IF(LAS.EQ.3) GO TO 4
      CG4=(C41*4.D0+C42*2.D0)/ANOR
      SG4=(S41*4.D0+S42*2.D0)/ANOR
      CSG22=(CS221*4.D0+CS222*2.D0)/ANOR
      CSG31=(CS311*4.D0+CS312*2.D0)/ANOR
      CSG13=(CS131*4.D0+CS132*2.D0)/ANOR
    4 CONTINUE
      RETURN
      END
C     *****************************************************************
      SUBROUTINE OVLAB
C     *****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM9.FOR'
      AH=0.2D0
      C1=DSQRT(2.D0)
      C=X2*P1/P2/X1
      IF(C.GT.1.) GO TO 1
      STEP1=AH
      STEP2=AH*C
      GO TO 2
    1 STEP2=AH
      STEP1=AH/C
    2 STEP=-STEP1*P1/X1/C1
      N=1
      SQX2=C1*X2
      SQX1=C1*X1
      ANU=ANU1
      Y=SQX1+STEP1
      CALL FUDNU
      DN12=DN
      ANU=ANU2
      Y=SQX2+STEP2
      CALL FUDNU
      DN22=DN
      SUM2=0.D0
      S1M2=0.D0
      S2M2=0.D0
      FUYN=STEP**NNT
      IF(NNT.EQ.10) FUYN=(STEP-1.D0)/STEP/STEP
      IF(NNT.EQ.11) FUYN=(STEP**3-1.D0)/STEP/STEP
      IF(NNT.EQ.12) FUYN=STEP-1.D0
      IF(NNT.EQ.13) FUYN=1.D0/STEP/STEP
      SUM1=FUYN*DN12*DN22
      S1M1=DN12*DN12
      S2M1=DN22*DN22
    6 N=N+1
      YY=STEP*N
      FUYN=YY**NNT
      IF(NNT.EQ.10) FUYN=(YY-1.D0)/YY/YY
      IF(NNT.EQ.11) FUYN=(YY**3-1.D0)/YY/YY
      IF(NNT.EQ.12) FUYN=YY-1.
      IF(NNT.EQ.13) FUYN=1./YY/YY
      ANU=ANU1
      Y=SQX1+N*STEP1
      CALL FUDNU
      DN13=DN
      ANU=ANU2
      Y=SQX2+N*STEP2
      CALL FUDNU
      DN23=DN
      IF (N/2*2.EQ.N) GO TO 3
      SUM1=SUM1+FUYN*DN13*DN23
      S1M1=S1M1+DN13*DN13
      S2M1=S2M1+DN23*DN23
      GO TO 4
    3 SUM2=SUM2+FUYN*DN13*DN23
      S1M2=S1M2+DN13*DN13
      S2M2=S2M2+DN23*DN23
    4 IF(DABS(DN13).LT.1.D-4*DSQRT(S1M1).AND.DABS(DN23).
     *LT.1.D-4*DSQRT(S2M1)) GO TO 5
      GO TO 6
    5 FNR1=S1M2*2.D0+S1M1*4.D0
      FNR2=S2M2*2.D0+S2M1*4.D0
      FOLAR=SUM2*2.D0+SUM1*4.D0
      FOLAR=FOLAR/DSQRT(FNR1*FNR2)
      RETURN
      END
C     *****************************************************************
      SUBROUTINE OVLAO
C     *****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM9.FOR'
      AH=0.2D0
      C1=DSQRT(2.D0)
      STEPP=AH
      STEP=STEPP*AMUO/C1
      N=1
      ANU=ANU1
      Y=STEPP
      CALL FUDNU
      DN12=DN
      ANU=ANU2
      CALL FUDNU
      DN22=DN
      SUM2=0.D0
      S1M2=0.D0
      S2M2=0.D0
      FUYN=STEP**NNT
      SUM1=FUYN*DN12*DN22
      S1M1=DN12*DN12
      S2M1=DN22*DN22
    6 N=N+1
      YY=STEP*N
      FUYN=YY**NNT
      ANU=ANU1
      Y=N*STEPP
      CALL FUDNU
      DN13=DN
      ANU=ANU2
      CALL FUDNU
      DN23=DN
      IF (N/2*2.EQ.N) GO TO 3
      SUM1=SUM1+FUYN*DN13*DN23
      S1M1=S1M1+DN13*DN13
      S2M1=S2M1+DN23*DN23
      GO TO 4
    3 SUM2=SUM2+FUYN*DN13*DN23
      S1M2=S1M2+DN13*DN13
      S2M2=S2M2+DN23*DN23
    4 IF(DABS(DN13).LT.1.D-4*DSQRT(S1M1).AND.DABS(DN23).
     *LT.1.D-4*DSQRT(S2M1)) GO TO 5
      GO TO 6
    5 FNR1=S1M2*2.D0+S1M1*4.D0
      FNR2=S2M2*2.D0+S2M1*4.D0
      FOLAR=SUM2*2.D0+SUM1*4.D0
      FOLAR=FOLAR/DSQRT(FNR1*FNR2)
      RETURN
      END
C     *****************************************************************
      SUBROUTINE INERMO
C     *****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM9.FOR'
      INCLUDE 'PRIVCOM10.FOR'
      INCLUDE 'PRIVCOM12.FOR'     
      PI23=2.0943951024D0
      C1=DSQRT(2.D0)
      C2=3.16227766017D0
      C3=2.64575131106D0
      ARG1=GAM-PI23
      ARG2=GAM-2.D0*PI23
      ARG3=GAM
      ARG1=GAM-PI23+GSHAPE
      ARG2=GAM-2.*PI23+GSHAPE
      ARG3=GAM+GSHAPE
      AMO=DSIN(ARG1)**2
      BMO=DSIN(ARG2)**2
      CMO=DSIN(ARG3)**2
      PBE22=1./PBET2**2
      IF(MESHA-2) 1,2,3
    2 AMO=AMO+BB42*2.5D0*PBE22
      BMO=BMO+BB42*2.5D0*PBE22
      GO TO 1
    3 C7=7.D0/12.D0
      C5=5.D0/12.D0
      CQ7=DSQRT(C7)
      CQ5=DSQRT(C5)
      IF(MESHA.EQ.4) CSDG=DCOS(DELG)
      IF(MESHA.EQ.4) GO TO 4
      GAMG=GAM
      CSDG=CQ7*DCOS(3.D0*GAM)
    4 SSDG=DSQRT(1.D0-CSDG*CSDG)
      CSGG=DCOS(GAMG)
      SSGG=DSIN(GAMG)
      A40=CQ7*CSDG+CQ5*SSDG*CSGG
      A42=-SSDG*SSGG/C1
      A44=(CQ5*CSDG-CQ7*SSDG*CSGG)/C1
      AB1=2.5D0*A40*A40+4.*A42*A42+A44*A44
      AB2=1.5D0*C2*A40*A42+C3*A42*A44
      AMO=AMO+BB42*(AB1+AB2)*PBE22
      BMO=BMO+BB42*(AB1-AB2)*PBE22
      CMO=CMO+2.D0*BB42*(A42*A42+4.D0*A44*A44)*PBE22
    1 PBET32=PBET3**2*PBE22
      IF(MESHO-1) 5,6,7
    6 AMO=AMO+BB32*1.5D0*PBET32
      BMO=BMO+BB32*1.5D0*PBET32
      GO TO 5
    7 C30=DSQRT(30.D0)
      A30=DCOS(ETO)
      A32=DSIN(ETO)/C1
      AB1=1.5D0*A30*A30+2.D0*A32*A32
      AB2=0.5D0*C30*A30*A32
      AMO=AMO+BB32*(AB1+AB2)*PBET32
      BMO=BMO+BB32*(AB1-AB2)*PBET32
      CMO=CMO+2.D0*BB32*A32*A32*PBET32
    5 AMO=1.D0/AMO
      BMO=1.D0/BMO
      CMO=1.D0/CMO
      RETURN
      END
C     *************************************************************
      DOUBLE PRECISION FUNCTION DGAMMA(XX)
C     *************************************************************
      DOUBLE PRECISION COF(6),STP,HALF,ONE,FPF,X,TMP,SER,XX
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=TMP**(X+HALF)*DEXP(-TMP)
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      DGAMMA=TMP*STP*SER
      RETURN
      END
C     *****************************************************************
      SUBROUTINE TRANS
C     *****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INCLUDE 'PRIVCOM1.FOR'
      INCLUDE 'PRIVCOM8.FOR'
      INCLUDE 'PRIVCOM9.FOR'
      INCLUDE 'PRIVCOM10.FOR'
      INCLUDE 'PRIVCOM12.FOR' 
      INCLUDE 'PRIVCOM13.FOR'   
  
         IF(MEPRI.LT.98) PRINT 10
          WRITE(21,10)
          IF(MEHAM.EQ.5 .AND. MEPRI.LT.98)  PRINT 13
          IF(MEHAM.EQ.5)  WRITE(21,13)
          IF(MEHAM.EQ.3 .AND. MEPRI.LT.98)  PRINT 14
          IF(MEHAM.EQ.3)  WRITE(21,14)
          IF(MEPRI.LT.98) PRINT 15
          WRITE(21,15)
          IF(MEPRI.LT.98) PRINT 11
          WRITE(21,11)
   10     FORMAT(30X,'E2-TRANSITION PROBABILITIES B(E2:I->F)')
   11     FORMAT(1X,120(1H-))
   13     FORMAT(20X,'( CALCULATIONS INVOLVE  5-PARAM. MODEL')

   14     FORMAT(20X,'( CALCULATIONS INVOLVE DAVYDOV-CHABAN MODEL')
   15     FORMAT(20X,'TAKING INTO ACCOUNT SQUARE A(LAM,MU) TERMS )')
   16     FORMAT (30X,'E3-TRANSITION PROBABILITIES B(E3;I->F)')
  110     FORMAT(30X,' E4-TRANSITION PROBABILITIES B(E4:I->F)')
      C1=DSQRT(2.D0)
      PI=4.D0*DATAN(1.D0)
      SQPI=DSQRT(PI)
      SQ5=DSQRT(5.D0)
      CQ57=DSQRT(5.D0/7.D0)
      CQ59=DSQRT(5.D0)/3.D0
      CQ79=DSQRT(7.D0)/3.D0
      C00=5.D0*SQ5/21.D0/SQPI
      IF(MESHA-2) 21,22,23
C     LAMBDA4 AXIAL DEFORMATION
   22 A40=1.D0
      A42=0.D0
      A44=0.D0
      GO TO 21
C     LAMBDA4 NON-AXIAL DEFORMATION
   23 C712=7.D0/12.D0
      C512=5.D0/12.D0
      CQ712=DSQRT(C712)
      CQ512=DSQRT(C512)
      IF(MESHA.EQ.4) CSDG=COS(DELG)
      IF(MESHA.EQ.4) GO TO 24
      GAMG=GAM
      CSDG=CQ712*DCOS(3.D0*GAM)
   24 SSDG=DSQRT(1.D0-CSDG*CSDG)
      CSGG=DCOS(GAMG)
      SSGG=DSIN(GAMG)
      A40=CQ712*CSDG+CQ512*SSDG*CSGG
      A42=-SSDG*SSGG/C1
      A44=(CQ512*CSDG-CQ712*SSDG*CSGG)/C1
   21 IF(MESHO-1) 25,26,27
C     LAMBDA3 AXIAL DEFORMATION
   26 A30=1.D0
      A32=0.D0
      GO TO 25
C     LAMBDA3 NON-AXIAL DEFORMATION
   27 A30=DCOS(ETO)
      A32=DSIN(ETO)/C1
   25 SIGO=DSIN(GAM0+GSHAPE)
      COGO=DCOS(GAM0+GSHAPE)
      FOLAC=DCOS(GAM0+GSHAPE)
      FOLAS=DSIN(GAM0+GSHAPE)
      FOLC2=DCOS(2.D0*(GAM0+GSHAPE))
      FOLS2=DSIN(2.D0*(GAM0+GSHAPE))
      FLC2=DCOS(GAM0+GSHAPE)**2
      FLS2=DSIN(GAM0+GSHAPE)**2
      FOLC3=DCOS(3.D0*(GAM0+GSHAPE))
      C3G0=DCOS(3.D0*(GAM0+GSHAPE))
      FOKS=BET3
      FOKS1=BET3*BET3
C
C     VECTOR COEFFICIENTS FOR E2-,E3- ˆ E4-TRANSITIONS
C
      QQQ00=-0.5345224838D0
      QQQ22=-QQQ00
      GQQ00=0.7171371656D0
      GQQ20=0.1195228609D0
      GQQ22=0.4629100500D0
      GQQ24=1.D0
      GQG00=-0.5096471915D0
      GQG20=0.5921565255D0
      GQG22=-0.2038588766D0
      GQG02=GQG20
      GQG42=0.3302891295D0
      GQG44=0.7135060680D0
      GQG24=GQG42
      QGG00=CQ59*GQG00
      QGG20=CQ59*GQG22
      QGG02=CQ59*GQG02
      QGG42=CQ59*GQG24
      QGG40=CQ59*GQG44
      QQG00=CQ59*GQQ00
      QQG20=CQ59*GQQ22
      QQG22=QQG20
      QQG02=CQ59*GQQ20
      QQG42=CQ59*GQQ24
      GOO00=-0.4834937784D0
      GOO20=0.5640760748D0
      GOO02=0.1395726315D0
      GOO24=-0.6741998625D0
      GOO22=GOO02
      GGG00=0.4022911408D0
      GGG20=-0.245844586D0
      GGG40=0.3128931094D0
      GGG22=GGG20
      GGG42=0.560968194D0
      GGG44=GGG40
      GGG24=GGG42
      OQO00=-0.5163977795D0
      OQO20=0.5773502692D0
      OQO02=OQO20
      OQO22=0.D0
      OOG00=-CQ79*GOO00
      OOG20=-CQ79*GOO20
      OOG22=-CQ79*GOO22
      OOG02=-CQ79*GOO20
      OOG42=-CQ79*GOO24
      QOO00=-CQ57*OQO00
      QOO02=-CQ57*OQO02
      QOO20=-CQ57*OQO22
C
C     CONSTANTS FOR E2-TRANSITIONS
C     Q20=3.*Z*e*RR*2*BETTA20/SQRT(5*PI)
C     CONST=5*Q20**2/16/PI   RESULTS TO BE MULTIPLIED BY
C
      QCG1A=6.D0*BET4/SQPI
      QSG1A=QCG1A*C1
      QC2GA=SQ5/SQPI*BET0
      QG0A=3.D0*BET4**2/SQ5/SQPI/BET0
      QO0A=7.D0/SQPI/SQ5/BET0
C     QO0A=7.D0/SQPI/SQ5/BET0*BET3**2
      QSG1B=QCG1A/C1
      QCG1B=QCG1A
      QS2GB=QC2GA/C1
      QG0B=QG0A*2.D0
      QO0B=QO0A*2.D0
      AC2=1.D0+QCG1A*QQG00**2*A40
      AS2=QSG1A*QQG00*QQG20*A42
      AC22=QC2GA*QQQ00**2
      AQG2=QG0A*QGG00*(2.D0*QGG40*A44**2+2.D0*QGG20*A42**2+
     1QGG00*A40**2)
      AQO2=QO0A*QOO00*(QOO00*A30**2+2.D0*QOO20*A32**2)
      BS2=1.D0/C1+QSG1B*QQG00*(QQG02*A40+QQG42*A44)
      BC2=QCG1B*QQG00*QQG22*A42
      BS22=QS2GB*QQQ00*QQQ22
      BQG2=QG0B*QGG00*(QGG02*A40*A42+QGG42*A42*A44)
      BQO2=QO0B*QOO00*QOO02*A30*A32
C
C     CONSTANTS FOR E3-TRANSITIONS
C     Q30=3.*Z*e*RR*3*BETTA30/SQRT(7*PI)
C     CONST=7*Q30**2/16/PI RESULTS TO BE MULTIPLIED BY
C
      OC1=2.5D0*SQ5/SQPI*BET0
      OC2=7.5D0*BET4/SQPI
      AC3=OC1*OQO00**2*DCOS(ETO)
      AS3=OC1*OQO00*OQO20*DSIN(ETO)*0.5D0
      AC00=DCOS(ETO)
      AOG3=OC2*OOG00*(OOG00*DCOS(ETO)*A40+C1*OOG20*DSIN(ETO)*A42)
      BC3=OC1*OQO00*OQO22*DSIN(ETO)/C1
      BS3=OC1*OQO00*OQO02*DCOS(ETO)/C1
      BS00=DSIN(ETO)/C1
      BOG3=OC2*OOG00*(OOG22*DCOS(ETO)*A42+DSIN(ETO)/C1*
     1(OOG02*A40+OOG42*A44))
C
C     CONSTANTS FOR E4-TRANSITIONS
C     Q40=3.*Z*e*RR*4*BETTA40/SQRT(9*PI)
C     CONST=9*Q40**2/16/PI  RESULTS TO BE MULTIPLIED BY
C
      A1=2.5D0/SQPI*BET0**2/BET4
      A2=3.D0*SQ5/SQPI*BET0
      A3=4.5D0/SQPI*BET4
      A4=3.5D0/SQPI/BET4
C     A4=3.5D0/SQPI/BET4*BET3**2
      BB1=A1/C1
      BB2=A2
      BB3=A3*2.D0
      BB4=A4*2.D0
      G1=A1/2.D0
      G2=A2
      G3=A3
      G4=A4
      AC4=A2*GQG00**2*A40
      AS4=A2*GQG00*GQG20*C1*A42
      AC24=A1*GQQ00**2
      AS24=A1*GQQ00*GQQ20
      AGG4=A3*GGG00*(GGG00*A40**2+2.D0*GGG20*A44**2+
     12.D0*GGG40*A44**2)
      AGO4=A4*GOO00*(GOO00*A30**2+2.D0*GOO20*A32**2)
      BC4=BB2*GQG00*GQG22*A42
      BS4=BB2*GQG00*(GQG02*A40+GQG42*A44)/C1
      BS24=BB1*GQQ00*GQQ22
      BGG4=BB3*GGG00*(GGG22*A40*A42+GGG42*A44*A42)
      BGO4=BB4*GOO00*GOO22*A30*A32
      CC4=G2*GQG00**2*A44
      CS4=G2*GQG00*GQG24*A42/C1
      CS24=G1*GQG00*GQG24
      CGG4=G3*GGG00*(2.D0*GGG44*A40*A44+GGG24*A42**2)
      CGO4=G4*GOO00*GOO24*A32**2
C     CALCULATIONS OF E2-TRANSITION PROBABILITIES
      J2=4
      DO 1 I=1,NUR
      JI=JU(I)
      PO1=(-1)**NNO(I)
      DG1=0.
      IF(PO1.NE.1) DG1=GAMDE
      KG1=(3-PO1)/2
      NG1I=NNG(I)+1
      ANG1=ANG(NG1I,KG1)
      NO1I=NNO(I)+1
      AONU1=ANO(NO1I)
      CD1=CD(NG1I,KG1)
      J1=2*JI
      KMI=1-(-1)**(JI+NNO(I))
      NKI=(2*(JI/2)-KMI+2)/2
      DO 1 JJ=1,NUR
      JF=JU(JJ)
      PO2=(-1)**NNO(JJ)
      DG2=0.
      IF(PO2.NE.1) DG2=GAMDE
      KG2=(3-PO2)/2
      NG1F=NNG(JJ)+1
      ANG2=ANG(NG1F,KG2)
      NO1F=NNO(JJ)+1
      AONU2=ANO(NO1F)
      CD2=CD(NG1F,KG2)
      J=2*JF
      JIF=IABS(JI-JF)
      IF(JIF.GT.2) GO TO 111
      IF(JI.EQ.0.AND.JF.EQ.0) GO TO 111
      IF(PO1.NE.PO2) GO TO 111
      IF(MEHAO.EQ.0) GO TO 45
      IF(MEHAO.LT.2) GO TO 51
      EBM=DEXP(-(BET3/AMUO)**2)
      FOKS1=AMUO**2/2.D0+BET3**2/(1.D0-(-1)**NO1I*EBM)
      IF(MEHAO.GE.2) GO TO 45
   51 NNT=2
      X1=-BET3/AMUO
      X2=X1
      P1=BET3
      P2=BET3
      ANU1=AONU1
      ANU2=AONU2
      IF(BET3.EQ.0.D0) CALL OVLAO
      IF(BET3.NE.0.D0) CALL OVLAB
      FOKS1=FOLAR
   45 IF(MEHAM.GE.5)  CALL TRLAG
      X1= XI(I)
      X2= XI(JJ)
      P1=PT(I)
      P2=PT(JJ)
      ANU1=ANB(I)
      ANU2=ANB(JJ)
      FOLAR2=1.D0
      IF(MEHAM.EQ.4) GO TO 12
      NNT=3
      CALL OVLAB
      FOLAR2=FOLAR
   12 NNT=2
      IF(MEHAM.NE.4) CALL OVLAB
      IF(MEHAM.EQ.4) FOLAR=1.
      FOLAR1=FOLAR
      IF(MEHAO.EQ.2) FOKS1=FOKS1*FOLAR1*BET0**2
      NNT=1
      IF(MEHAM.NE.4) CALL OVLAB
      IF(MEHAM.EQ.4) FOLAR=1.
      NNTG=1
      IF(MEHAM.GT.4)CALL OVLAG
      KMF=1-(-1)**(JF+NNO(JJ))
      NKF=(2*(JF/2)-KMF+2)/2
      SUM=0.D0
      DO 2 JKI=1,NKI
      KI=KMI+2*(JKI-1)
      DI=1.D0
      IF(KI.EQ.0) DI=C1
      DO 2 JKF=1,NKF
      KF=KMF+2*(JKF-1)
      SUM1=0.D0
      DF=1.D0
      IF(KF.EQ.0) DF=C1
      M=2*KF
      IF(KF.NE.KI) GO TO 3
      M1=2*KI
      M2=0
      CALL KLEGO
      T1=AKG
      IF(KF.NE.0) GO TO 4
      M1=-2*KI
      CALL KLEGO
      T1=T1+AKG*(-1)**(JI+NNO(I))
    4 SUM1=SUM1 +T1*((AC2*FOLAC+AS2*FOLAS)*FOLAR+
     1AC22*FOLC2*FOLAR1+AQG2+AQO2*FOKS1)
    3 T1=0.
      IF((KI+2).NE.KF) GO TO 5
      M1=2*KI
      M2=4
      CALL KLEGO
      T1=AKG
    5 IF((KI-2).NE.KF) GO TO 6
      M1=2*KI
      M2=-4
      CALL KLEGO
      T1=T1+AKG
    6 IF((2-KI).NE.KF) GO TO 7
      M1=-2*KI
      M2=4
      CALL KLEGO
      T1=T1+AKG*(-1)**(JI+NNO(I))
    7 SUM1=SUM1 +T1*((BC2*FOLAC+BS2*FOLAS)*FOLAR+BS22*FOLS2*FOLAR1+
     1BQG2+BQO2*FOKS1)
      SUM=SUM+SUM1*AIT(I,JKI)*AIT(JJ,JKF)/DI/DF
    2 CONTINUE
      IF(I.NE.JJ) GO TO 20
      Q1=DSQRT(JI*(2.D0*JI-1.D0)
     * /((JI+1.D0)*(2.D0*JI+1.)*(2.D0*JI+3.D0)))
      Q1=Q1*SUM*DSQRT(2.D0*JI+1.D0)
      IF(MEPRI.LT.98) 
     *  PRINT 17,ES(I),JI,NTU(I),NNB(I),NNG(I),NNO(I),Q1
      WRITE(21,17)ES(I),JI,NTU(I),NNB(I),NNG(I),NNO(I),Q1
   17 FORMAT(10X,'AVERAGE VALUE OF LAMBDA2 MOMENT FOR STATE  E=',
     *D11.4,'JI=',I2,'TAU=',I1,'NB=',I1,'NG=',I1,'NO=',I1,3X,
     *'DEVIDED BY Q20=',D11.4)
      GO TO 1
   20 BE2(I,JJ)=SUM*SUM
      IF(MEPRI.LT.98) PRINT 11
      WRITE(21,11)
      IF(MEPRI.LT.98) 
     *  PRINT 8, ES(I),JI,NTU(I),NNB(I),NNG(I),NNO(I),
     *         ES(JJ),JF,NTU(JJ),NNB(JJ),NNG(JJ),NNO(JJ),BE2(I,JJ)
      WRITE(21,8) ES(I),JI,NTU(I),NNB(I),NNG(I),NNO(I),
     *         ES(JJ),JF,NTU(JJ),NNB(JJ),NNG(JJ),NNO(JJ),BE2(I,JJ)
      IF(MEPRI.LT.98) PRINT 11
      WRITE(21,11)
      IF(MEPRI.LT.98) 
     *  PRINT 9,FOLAC,FOLAS,FOLAR,COGO,SIGO,FOKS1,FOLAG
      WRITE(21,9)FOLAC,FOLAS,FOLAR,COGO,SIGO,FOKS1,FOLAG
    9 FORMAT(5X,'I(COS)=',D11.4,'  I(SIN)=',D11.4,
     *'  FOLAR=',D11.4,'  COS(G0)=',D11.4,'  SIN(G0)=',
     *D11.4,'  FOKS1=', D11.4,' FOLAG=',D11.4//)
    8 FORMAT(5X,'B(E2:(',D11.4,')JI=',I2,'TAU=',I1,'NB=',
     *I1,'NG=',I1,'NO=',I1,'--->(',D11.4,')JF=',I2,'TAU=',I1,'NB=',
     *I1,'NG=',I1,'NO=',I1,')=',D11.4)
C     LAMBDA0 TRANSITION PROBABILITIES
  111 IF(JIF.NE.0) GO TO 1
      IF(I.EQ.JJ) GO TO 1
      IF(MEHAM.EQ.4) GO TO 1
      BE0(I,JJ)=(BET0*BET0*FOLAR1+C00*BET0**3*FOLAR2*FOLC3)**2
      IF(MEPRI.LT.98) PRINT 11
      WRITE(21,11)
      IF(MEPRI.LT.98) 
     *  PRINT 81,ES(I),JI,NTU(I),NNB(I),NNG(I),NNO(I),
     *         ES(JJ),JF,NTU(JJ),NNB(JJ),NNG(JJ),NNO(JJ),BE0(I,JJ)
      WRITE(21,81)ES(I),JI,NTU(I),NNB(I),NNG(I),NNO(I),
     *         ES(JJ),JF,NTU(JJ),NNB(JJ),NNG(JJ),NNO(JJ),BE0(I,JJ)
   81 FORMAT(5X,'B(E0:(',D11.4,')JI=',I2,'TAU=',I1,'NB=',
     *I1,'NG=',I1,'NO=',I1,'--->(',D11.4,')JF=',I2,'TAU=',I1,'NB=',
     *I1,'NG=',I1,'NO=',I1,')=',D11.4)
      IF(MEPRI.LT.98) PRINT 91,FOLC3,FOLAR2,FOLAR1,C3G0
      WRITE(21,91)FOLC3,FOLAR2,FOLAR1,C3G0
   91 FORMAT(5X,'I(COS(3*GAM0))=',D11.4,' FOLAR2 =',D11.4,
     *' FOLAR1=',D11.4,'  COS(3*G0)=',D11.4/)
    1 CONTINUE
      IF(MESHA.LT.2) GO TO 100
          IF(MEPRI.LT.98) PRINT 110
          WRITE(21,110)
          IF(MEPRI.LT.98) PRINT 11
          WRITE(21,11)
C     CALCULATIONS OF E4- TRANSITION PROBABILITIES
      J2=8
      DO 19 I=1,NUR
      JI=JU(I)
      PO1=(-1)**NNO(I)
      DG1=0.
      IF(PO1.NE.1) DG1=GAMDE
      KG1=(3-PO1)/2
      NG1I=NNG(I)+1
      ANG1=ANG(NG1I,KG1)
      NO1I=NNO(I)+1
      AONU1=ANO(NO1I)
      CD1=CD(NG1I,KG1)
      J1=2*JI
      KMI=1-(-1)**(JI+NNO(I))
      NKI=(2*(JI/2)-KMI+2)/2
      DO 19 JJ=1,NUR
      JF=JU(JJ)
      PO2=(-1)**NNO(JJ)
      DG2=0.
      IF(PO2.NE.1) DG2=GAMDE
      KG2=(3-PO2)/2
      NG1F=NNG(JJ)+1
      ANG2=ANG(NG1F,KG2)
      NO1F=NNO(JJ)+1
      AONU2=ANO(NO1F)
      CD2=CD(NG1F,KG2)
      J=2*JF
      JIF=IABS(JI-JF)
      IF(JIF.GT.4) GO TO 19
      IF(JI.EQ.0.AND.JF.EQ.0) GO TO 19
      IF(PO1.NE.PO2) GO TO 19
      IF(MEHAO.EQ.0) GO TO 46
      IF(MEHAO.LT.2) GO TO 52
      EBM=DEXP(-(BET3/AMUO)**2)
      FOKS1=AMUO**2/2.D0+BET3**2/(1.D0-(-1)**NO1I*EBM)
      IF(MEHAO.GE.2) GO TO 46
   52 NNT=2
      X1=-BET3/AMUO
      X2=X1
      P1=BET3
      P2=BET3
      ANU1=AONU1
      ANU2=AONU2
      IF(BET3.EQ.0.D0) CALL OVLAO
      IF(BET3.NE.0.D0) CALL OVLAB
      FOKS1=FOLAR
   46 IF(MEHAM.GE.5)  CALL TRLAG
      FLC2=(1.D0+FOLC2)/2.D0
      FLS2=(1.D0-FOLC2)/2.D0
      X1= XI(I)
      X2= XI(JJ)
      P1=PT(I)
      P2=PT(JJ)
      ANU1=ANB(I)
      ANU2=ANB(JJ)
      NNT=2
      IF(MEHAM.NE.4) CALL OVLAB
      IF(MEHAM.EQ.4) FOLAR=1.D0
      FOLAR1=FOLAR
      IF(MEHAO.EQ.2) FOKS1=FOKS1*FOLAR1*BET0**2
      NNT=1
      IF(MEHAM.NE.4) CALL OVLAB
      IF(MEHAM.EQ.4) FOLAR=1.
      KMF=1-(-1)**(JF+NNO(JJ))
      NKF=(2*(JF/2)-KMF+2)/2
      SUM=0.D0
      DO 120 JKI=1,NKI
      KI=KMI+2*(JKI-1)
      DI=1.D0
      IF(KI.EQ.0) DI=C1
      DO 120 JKF=1,NKF
      KF=KMF+2*(JKF-1)
      SUM1=0.D0
      DF=1.D0
      IF(KF.EQ.0) DF=C1
      M=2*KF
      IF(KF.NE.KI) GO TO 30
      M1=2*KI
      M2=0
      CALL KLEGO
      T1=AKG
      IF(KF.NE.0) GO TO 40
      M1=-2*KI
      CALL KLEGO
      T1=T1+AKG*(-1)**(JI+NNO(I))
   40 SUM1=SUM1 +T1*(A40+(AC4*FOLAC+AS4*FOLAS)*FOLAR+
     1(AC24*FLC2+AS24*FLS2)*FOLAR1+AGG4+AGO4*FOKS1)
   30 T1=0.
      IF((KI+2).NE.KF) GO TO 50
      M1=2*KI
      M2=4
      CALL KLEGO
      T1=AKG
   50 IF((KI-2).NE.KF) GO TO 60
      M1=2*KI
      M2=-4
      CALL KLEGO
      T1=T1+AKG
   60 IF((2-KI).NE.KF) GO TO 70
      M1=-2*KI
      M2=4
      CALL KLEGO
      T1=T1+AKG*(-1)**(JI+NNO(I))
   70 SUM1=SUM1 +T1*(A42+(BC4*FOLAC+BS4*FOLAS)*FOLAR+
     1BS24*FOLS2*FOLAR1+BGG4+BGO4*FOKS1)
      T1=0.
      IF((KI+4).NE.KF) GO TO 80
      M1=2*KI
      M2=8
      CALL KLEGO
      T1=T1+AKG
   80 IF((KI-4).NE.KF) GO TO 90
      M1=2*KI
      M2=-8
      CALL KLEGO
      T1=T1+AKG
   90 IF((4-KI).NE.KF) GO TO 95
      M1=-2*KI
      M2=8
      CALL KLEGO
      T1=T1+(-1)**JI*AKG
   95 SUM1=SUM1+T1*(A44+(CC4*FOLAC+CS4*FOLAS)*FOLAR+CS24*FLS2*FOLAR1+
     1CGG4+CGO4*FOKS1)
      SUM=SUM+SUM1*AIT(I,JKI)*AIT(JJ,JKF)/DI/DF
  120 CONTINUE
      IF(I.NE.JJ) GO TO 130
      Q4=SUM
      IF(MEPRI.LT.98) 
     * PRINT 117,ES(I),JI,NTU(I),NNB(I),NNG(I),NNO(I),Q4
      WRITE(21,117)ES(I),JI,NTU(I),NNB(I),NNG(I),NNO(I),Q4
  117 FORMAT(10X,'AVERAGE VALUE OF LAMBDA4 MOMENT FOT STATE  E=',
     *D11.4,'JI=',I2,'TAU=',I1,'NB=',I1,'NG=',I1,'NO=',I1,3X,
     *'DEVIDED BY Q40=',D11.4)
      GO TO 19
  130 BE4(I,JJ)=SUM*SUM
      IF(MEPRI.LT.98) PRINT 11
      WRITE (21,11)
      IF(MEPRI.LT.98) PRINT 88,ES(I),JI,NTU(I),NNB(I),NNG(I),NNO(I),
     *         ES(JJ),JF,NTU(JJ),NNB(JJ),NNG(JJ),NNO(JJ),BE4(I,JJ)
      WRITE(21,88)ES(I),JI,NTU(I),NNB(I),NNG(I),NNO(I),
     *         ES(JJ),JF,NTU(JJ),NNB(JJ),NNG(JJ),NNO(JJ),BE4(I,JJ)
      IF(MEPRI.LT.98) PRINT 11
      WRITE(21,11)
      IF(MEPRI.LT.98) PRINT 99,FLC2,FLS2,FOLS2,FOLAR1,FOKS1
      WRITE(21,99)FLC2,FLS2,FOLS2,FOLAR1,FOKS1
   99 FORMAT(5X,'I(CS**2)',D11.4,'I(SS**2)=',D11.4,
     *'I(SS(2G0))',D11.4,' FOLAR1=',D11.4,' FOKS1=',D11.4//)
   88 FORMAT(5X,'B(E4:(',D11.4,')JI=',I2,'TAU=',I1,'NB=',
     *I1,'NG=',I1,'NO=',I1,'--->(',D11.4,')JF=',I2,'TAU=',I1,'NB=',
     *I1,'NG=',I1,'NO=',I1,')=',D11.4)
   19 CONTINUE
      IF(MEHAO.EQ.0) GO TO 100
      IF(MEPRI.LT.98) PRINT 16
      WRITE(21,16)
C     CALCULATIONS OF E3- TRANSITION PROBABILITIES
      J2=6
      DO 42 I=1,NUR
      PO1=(-1)**NNO(I)
      DG1=0.D0
      IF(PO1.NE.1) DG1=GAMDE
      KG1=(3-PO1)/2
      JI=JU(I)
      NG1I=NNG(I)+1
      ANG1=ANG(NG1I,KG1)
      NO1I=NNO(I)+1
      AONU1=ANO(NO1I)
      CD1=CD(NG1I,KG1)
      J1=2*JI
      KMI=1-(-1)**(JI+NNO(I))
      NKI=(2*(JI/2)-KMI+2)/2
      DO 42 JJ=1,NUR
      JF=JU(JJ)
      PO2=(-1)**NNO(JJ)
      DG2=0.D0
      IF(PO2.NE.1) DG2=GAMDE
      KG2=(3-PO2)/2
      NG1F=NNG(JJ)+1
      ANG2=ANG(NG1F,KG2)
      NO1F=NNO(JJ)+1
      AONU2=ANO(NO1F)
      CD2=CD(NG1F,KG2)
      J=2*JF
      JIF=IABS(JI-JF)
      IF(JIF.GT.3) GO TO 42
      IF(PO1.EQ.PO2) GO TO 42
      IF(JI.EQ.0.AND.JF.EQ.0) GO TO 42
      IF(MEHAO.EQ.0) GO TO 47
      IF(MEHAO.LT.2) GO TO 53
      EBM=DEXP(-(BET3/AMUO)**2)
      FOKS=BET3/DSQRT(1.D0-EBM**2)
      IF(MEHAO.GE.2) GO TO 47
   53 NNT=1
      X1=-BET3/AMUO
      X2=X1
      P1=BET3
      P2=BET3
      ANU1=AONU1
      ANU2=AONU2
      IF(BET3.EQ.0.D0) CALL OVLAO
      IF(BET3.NE.0.D0) CALL OVLAB
      FOKS=FOLAR
   47 FOKS=FOKS/BET3
      IF(MEHAO.EQ.2) FOKS=FOKS/BET0
      IF(MEHAM.GE.5)  CALL TRLAG
      NNT=2
      X1= XI(I)
      X2= XI(JJ)
      P1=PT(I)
      P2=PT(JJ)
      ANU1=ANB(I)
      ANU2=ANB(JJ)
      IF(MEHAM.NE.4) CALL OVLAB
      IF(MEHAM.EQ.4) FOLAR=1.
      FOLAR1=FOLAR
      NNT=1
      IF(MEHAM.NE.4) CALL OVLAB
      IF(MEHAM.EQ.4) FOLAR=1.
      KMF=1-(-1)**(JF+NNO(JJ))
      NKF=(2*(JF/2)-KMF+2)/2
      SUM=0.
      DO 31 JKI=1,NKI
      KI=KMI+2*(JKI-1)
      DI=1.
      IF(KI.EQ.0) DI=C1
      DO 31 JKF=1,NKF
      KF=KMF+2*(JKF-1)
      SUM1=0.D0
      DF=1.D0
      IF(KF.EQ.0) DF=C1
      M=2*KF
      IF(KF.NE.KI) GO TO 32
      M1=2*KI
      M2=0
      CALL KLEGO
      T1=AKG
      IF(KF.NE.0) GO TO 33
      M1=-2*KI
      CALL KLEGO
      T1=T1+AKG*(-1)**(JI+NNO(I))
   33 IF(MEHAO.EQ.2) SUM1=SUM1+T1*((AC3*FOLAC+AS3*FOLAS)
     1*FOLAR1+(AC00+AOG3)*FOLAR)*FOKS*BET0
      IF(MEHAO.EQ.2) GO TO 32
      SUM1=SUM1 +T1*(AC00+(AC3*FOLAC+AS3*FOLAS)*FOLAR+
     1AOG3)*FOKS
   32 T1=0.
      IF((KI+2).NE.KF) GO TO 34
      M1=2*KI
      M2=4
      CALL KLEGO
      T1=AKG
      T1=T1+AKG*(-1)**(JI+NNO(I))
   34 IF((KI-2).NE.KF) GO TO 35
      M1=2*KI
      M2=-4
      CALL KLEGO
      T1=T1+AKG
   35 IF((2-KI).NE.KF) GO TO 36
      M1=-2*KI
      M2=4
      CALL KLEGO
      T1=T1+AKG*(-1)**(JI+NNO(I))
   36 IF(MEHAO.EQ.2) SUM1=SUM1+T1*((BC3*FOLAC+BS3*FOLAS)
     1*FOLAR1+(BS00+AOG3)*FOLAR)*FOKS*BET0
      IF(MEHAO.EQ.2) GO TO 43
      SUM1=SUM1 +T1*(BS00+(BC3*FOLAC+BS3*FOLAS)*FOLAR+
     1BOG3)*FOKS
   43 SUM=SUM+SUM1*AIT(I,JKI)*AIT(JJ,JKF)/DI/DF
   31 CONTINUE
      IF(I.NE.JJ) GO TO 37
      Q3= SUM
      IF(MEPRI.LT.98) PRINT 38,ES(I),JI,NTU(I),NNB(I),NNG(I),NNO(I),Q3
      WRITE(21,38)ES(I),JI,NTU(I),NNB(I),NNG(I),NNO(I),Q3
   38 FORMAT(10X,'AVERAGE VALUE OF LAMBDA3 MOMENT FOT STATE  E=',
     *D11.4,'JI=',I2,'TAU=',I1,'NB=',I1,'NG=',I1,'NO=',I1,3X,
     *'DEVIDED BY Q30=',D11.4)
   37 BE3(I,JJ)=SUM*SUM
      IF(MEPRI.LT.98) PRINT 11
      WRITE(21,11)
      IF(MEPRI.LT.98) PRINT 39, ES(I),JI,NTU(I),NNB(I),NNG(I),NNO(I),
     *         ES(JJ),JF,NTU(JJ),NNB(JJ),NNG(JJ),NNO(JJ),BE3(I,JJ)
      WRITE(21,39) ES(I),JI,NTU(I),NNB(I),NNG(I),NNO(I),
     *         ES(JJ),JF,NTU(JJ),NNB(JJ),NNG(JJ),NNO(JJ),BE3(I,JJ)
   39 FORMAT(5X,'B(E3:(',D11.4,')JI=',I2,'TAU=',I1,'NB=',
     *I1,'NG=',I1,'NO=',I1,'--->(',D11.4,')JF=',I2,'TAU=',I1,'NB=',
     *I1,'NG=',I1,'NO=',I1,')=',D11.4)
      IF(MEPRI.LT.98) PRINT 83,FOKS
      WRITE(21,83) FOKS
   83 FORMAT(5X,'FOKS=',D11.4/)
   42 CONTINUE
  100 STOP
      END
C     *****************************************************************
      SUBROUTINE SHEM
C     *****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM8.FOR'
      INCLUDE 'PRIVCOM9.FOR'
      INCLUDE 'PRIVCOM10.FOR'
      INCLUDE 'PRIVCOM12.FOR' 
      INCLUDE 'PRIVCOM13.FOR'   
c     COMMON/ENA/EN,EL(40),BET(10),NUR,NMAX,NPD,LAS
      PI=4.D0*DATAN(1.D0)
      SQPI=DSQRT(PI)
      AMB4=AMB0**4
      AMG4=AMG0**4
      PBET2=1.D0
      PBET3=1.D0
C     PI3=1.0471975512D0
      PI3=PI/3.D0      
      C1=DSQRT(2.D0)
      AMO2=AMUO**2/2.D0
      IF(MEHAO.EQ.0) GO TO 31
      NOM=0
      DO 30 I=1,NUR
      IF(NOM.LT.NNO(I)) NOM=NNO(I)
   30 CONTINUE
      NOM1=NOM+1
      XIT=-BET3/AMUO
      X1=XIT
      X2=XIT
      P1=BET3
      P2=BET3
      DO 32 I=1,NOM1
      NROOT=I-1
      IF(BET3.EQ.0.D0) ANO(I)=2*NROOT+1.
      IF(BET3.EQ.0.D0) GO TO 2
      CALL ANUDF
      ANO(I)=ANU
    2 IF(MEHAO.EQ.1) GO TO 45
      EBM=EXP(-XIT**2)
      FOL=EBM*(-1.D0/XIT/SQPI*(1-(-1)**I)+(-1)**I)-ERFC2(-XIT)
      FOLB(I)=FOL/(1-(-1)**I*EBM)
      GO TO 32
   45 ANU1=ANU
      ANU2=ANU
      IF(BET3.EQ.0.D0) GO TO 40
      NNT=1
      CALL OVLAB
      FOLB(I)=(FOLAR-BET3)/BET3
      GO TO 32
   40 FOLB(I)=(2.D0*ANO(I)+1.D0)*AMO2
   32 CONTINUE
   31 NGD=2
      IF(MEHAM.EQ.5.OR.MEHAM.EQ.7) NGD=0
      NGM=0
      IF(MEHAM.LE.4) GO TO 26
      DO 23 I=1,NUR
      IF(NGM.LT.NNG(I)) NGM=NNG(I)
   23 CONTINUE
      GAMM=GAM0
      NGMP=NGM+NGD+1
      DO 46 INO=1,2
      DO 24 I=1,3!NGMP
      NROOT=I-1
      XIT=-C1*GAM0/AMG0
      XIT1=C1*(PI3-GAM0)/AMG0
      IF(XIT1.GT.7.D0) GO TO 21
      CALL ANDET0
      ANG(I,INO)=ANU
      CD(I,INO)=CDV
      GO TO 24
   21 XIT=-GAM0/AMG0
      CALL ANUDF
      ANG(I,INO)=ANU
      CD(I,INO)=0.D0
   24 CONTINUE
      GAM0=GAM0+GAMDE
   46 CONTINUE
      GAM0=GAMM
      NGM1=NGM+1
      DO 25 I=1,NGM1
      JP1=I+NGD
      ANG1=ANG(I,1)
      CD1=CD(I,1)
      DO 25 J=I,JP1
      ANG2=ANG(J,1)
      CD2=CD(J,1)
      !IF(MEHAM.EQ.5.OR.MEHAM.EQ.7) GO TO 22
      !gamma-corrections are enabled if upper line is commented
      NNTG=1
      CALL OVLAG
      FOG(1,I,J)=FOLAG
      FOG(1,J,I)=FOG(1,I,J)
      NNTG=2
      CALL OVLAG
      FOG(2,I,J)=FOLAG
      FOG(2,J,I)=FOG(2,I,J)
   22 CONTINUE
      IF(I.NE.J) GO TO 25
      IF(MEHAM.NE.7) GO TO 25
      NNTG=11
      CALL OVLAG
      FOLGB(I)=FOLAG
   25 CONTINUE
      GAM=GAM0
      IF(BET3.EQ.0.) PBET3=0.D0
      IF(MESHA.EQ.1) GO TO 37
      BB420=BB42
      DEBET=0.001D0*BB420
      BB42=BB420+DEBET
      CALL INERMO
      AMOP=AMO
      BMOP=BMO
      CMOP=CMO
      BB42=BB420-DEBET
      CALL INERMO
      AMOM=AMO
      BMOM=BMO
      CMOM=CMO
      BB42=BB420
      A1D=(AMOP-AMOM)/DEBET/2.D0
      B1D=(BMOP-BMOM)/DEBET/2.D0
      C1D=(CMOP-CMOM)/DEBET/2.D0
      AB1DM=A1D-B1D
      AB1DP=A1D+B1D
   37 IF(MESHO.EQ.0) GO TO 26
      BB320=BB32
      IF(BET3.EQ.0.D0) PBET3=0.001D0/BB32*GAM**2
      IF(BET3.EQ.0.D0) GO TO 41
      DEBET=0.001*BB320
      BB32=BB320+DEBET
   41 CALL INERMO
      AMOP=AMO
      BMOP=BMO
      CMOP=CMO
      IF(BET3.EQ.0.D0) PBET3=0.D0
      IF(BET3.EQ.0.D0) GO TO 42
      BB32=BB320-DEBET
   42 CALL INERMO
      AMOM=AMO
      BMOM=BMO
      CMOM=CMO
      BB32=BB320
      IF(BET3.EQ.0.D0) DEBET=0.0005D0*GAM**2
      AOD=(AMOP-AMOM)/DEBET/2.D0
      BOD=(BMOP-BMOM)/DEBET/2.D0
      COD=(CMOP-CMOM)/DEBET/2.D0
      ABODM=AOD-BOD
      ABODP=AOD+BOD
      IF(BET3.EQ.0.) PBET3=0.D0
      GAMM=GAM0
   26 DO 1 I=1,NUR
      IS=JU(I)
      ITAU=NTU(I)
      NG=NNG(I)
      NB=NNB(I)
      NO=NNO(I)
      NROOT=NG
      GAM=GAM0
      IF((-1)**NO.NE.1) GAM=GAM+GAMDE
      CALL MATAM
      EP0(I)=EIN(IS+1,ITAU)
      IF(MEHAM.LE.4) GO TO 27
      CALL EIT12
      GIT(I)=GAM
      EP1(I)=EPIT1
      EP2(I)=EPIT2
      EP12(I)=EPIT12
   27 DO 20 J=1,NK
      AIT(I,J)=TM(J,ITAU)
   20 CONTINUE
      IF(MEHAM.EQ.4) EGB(I)=EP0(I)
      IF(MEHAM.EQ.4) GO TO 1
      IF(MEHAM.EQ.3) GO TO 29
      ANU=ANG(NG+1,(3-(-1)**NO)/2)
      EPG(I)=(ANU+0.5D0)*2./AMG0**2+EPIT1+EPIT2+EPIT12
      EE=EP0(I)+EPG(I)-EPG(1)
   29 IF(MEHAM.EQ.3) EE=EP0(I)
      IF(MEHAO.EQ.0) GO TO 43
      EDO=DPAR*(-1)**(NO+1)
      IF(I.EQ.1) EDO1=EDO
      IF(MEHAO.EQ.2) EE=EE+EDO-EDO1
   43 AA=AMB4*EE
      PIT0=1.D0
      DPI=100.D0
      IF(AA.EQ.0.D0) GO TO 5
    7 PIT=PIT0+DPI
      IF((PIT-1.D0)*PIT**3-AA) 8,6,9
    8 PIT0=PIT
      GO TO 7
    9 DPI=DPI/5.
      IF(DPI.LE.1.D-8) GO TO 6
      GO TO 7
    5 PIT=1.
    6 CONTINUE
      PT(I)=PIT
      XIT=-PIT/AMB0*(4.D0-3.D0/PIT)**0.25D0
      XI(I)=XIT
      NROOT=NB
      AP=AMB0/PIT
      AP2=AP*AP
      AP6=AP2*AP2*AP2
      CALL ANUDF
      ANB(I)=ANU
      EPIBO1=0.D0
      P1=PIT
      P2=PIT
      X1=XIT
      X2=XIT
      ANU1=ANU
      ANU2=ANU
      IF(MESHA.EQ.1.AND.MESHO.EQ.0) GO TO 33
      SUM1=0.D0
      SUM2=0.D0
      DO 11 J=1,NK
      AK=TM(J,ITAU)
      AKA=TM(J,ITAU)
      KA=KMI+(J-1)*2
      SUM1=SUM1+AK*AKA*KA*KA
      D=1.D0
      DD=0.D0
      IF(KA.EQ.0) DD=1.D0
      IF(KA.EQ.0) D=DSQRT(2.D0)
      IF(J.EQ.NK) GO TO 11
      AK1=TM(J+1,ITAU)
      AKA1=TM(J+1,ITAU)
      SUM2=SUM2+(AK1*AKA+AK*AKA1)*(1.D0+(-1)**(IS+NO)*DD)/D*
     *SQRT((IS+KA+2.D0)*(IS+KA+1.D0)*(IS-KA-1.D0)*(IS-KA))
   11 CONTINUE
      IF(MESHO.EQ.0) GO TO 34
      EPIOO=ABODP/8.D0*IS*(IS+1)+(0.25D0*COD-ABODP/8.D0)*SUM1+
     *ABODM/16.D0*SUM2
      IF(MEHAO.LT.1) GO TO 35
      NNT=13
      CALL OVLAB
      EPIBO1=HW*BB32*AMB0**2*FOLAR*EPIOO*FOLB(NO+1)
   35 NNT=10
      CALL OVLAB
      IF(BET3.EQ.0.D0) GO TO 33
      EPIBO1=-HW*BB32*FOLAR*AMB0**2*EPIOO+EPIBO1
   33 IF(MESHA.EQ.1) EPIB1=0.D0
      IF(MESHA.EQ.1) GO TO 28
      IF(MESHO.GT.0) GO TO 36
   34 NNT=10
      CALL OVLAB
   36 EPIB1=AB1DP/8.D0*IS*(IS+1)+(0.25D0*C1D-AB1DP/8.D0)*SUM1+
     *AB1DM/16.D0*SUM2
      EPIB1=EPIB1*FOLAR*AMB0**2*BB42*HW
   28 EPB(I)=EPIB1
      EGB(I)=HW*((ANU+0.5D0)*DSQRT(4.-3./PIT)+0.5*AP2*EE+0.5*AP6*EE*EE)
     *-EPIB1+EPIBO1
      IF(MEHAO.EQ.1) EGB(I)=EGB(I)+HWO*(ANO(NO+1)+0.5D0)
      IF(MEHAO.EQ.3) EGB(I)=EGB(I)+EDO
      IF(MEHAM.NE.7) GO TO 1
      G3=3.D0*GAM0
      CBG1=1.D0/AMG0**4+81.D0/4.D0*(1.D0/DSIN(G3)**2+3.D0*DCOS(G3)**2/
     *SIN(G3)**4)
      NNT=11
      X1=XIT
      P1=PIT
      X2=XIT
      P2=PIT
      ANU1=ANU
      ANU2=ANU
      CALL OVLAB
      EBG1=HW*AMB0**2/4.D0*BET0**2*CBG1*FOLAR*FOLGB(NG+1)*AMG0**2
      EGB(I)=EGB(I)+EBG1
    1 CONTINUE
      EGB1=EGB(1)
      DO 10 I=1,NUR
   10 EGB(I)=EGB(I)-EGB1
      RETURN
      END                        
C **********************************************************************
      DOUBLE PRECISION FUNCTION ERFC2(XX)
C **********************************************************************
C     Abramowitz & Stegun, p.299, Eq. 7.1.26 (Error < 1.5E-7)
C
      DOUBLE PRECISION P,A1,A2,A3,A4,A5,P1,R1,R2,R3,R4,R5,XX
      DATA P,A1,A2,A3,A4,A5/0.3275911d0,0.254829592d0,-0.284496736d0,
     *1.421413741d0,-1.453152027d0,1.061405429d0/
      P1=DEXP(-XX*XX)
      R1=1.D0/(1.D0+P*XX)
      R2=R1*R1
      R3=R2*R1
      R4=R3*R1
      R5=R4*R1
      ERFC2=(A5*R5+A4*R4+A3*R3+A2*R2+A1*R1)*P1
      RETURN
      END  
C     ****************************************************************** 
      SUBROUTINE OVLOPT
C     ****************************************************************** 
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PRIVCOM8.FOR'
      INCLUDE 'PRIVCOM9.FOR' !!TOO MUCH!
      INCLUDE 'PRIVCOM10.FOR'
      INCLUDE 'PRIVCOM12.FOR'
      INCLUDE 'PRIVCOM13.FOR'
      INCLUDE 'PRIVCOM20.FOR'

 
      CALL SHEM

      
      DO NU=1,NUR
          IS=JU(NU)
          NO=NNO(NU)
          ITAU=NTU(NU)
          KMI=1-(-1)**(IS+NO)
          KMA=2*(IS/2)
          NKIWE(NU)=(KMA-KMI+2)/2
          SUMK2=0
          WEIGH(:,NU)=0.0
          KWE(:,NU)=0
          DO KK=1,NKIWE(NU)
               WEIGH(KK,NU)=AIT(NU,KK)
               SUMK2=SUMK2+WEIGH(KK,NU)**2
               KWE(KK,NU)=(KMI+(KK-1)*2)*2
          END DO
          SUMK2=SUM(WEIGH(:,NU)**2)
          IF (DABS(SUMK2-1.0).GT.1D-10) then
              PRINT "(A,d12.4)", 
     *         '!WARNING! NONUNIT SUM OF Aitau^2: ',SUMK2 
              PAUSE 123
          END IF
      END DO

          DO II=1,NUR
            DO JJ=II,NUR
c             |BETA> inputs
              X1= XI(II)
              X2= XI(JJ)
              P1=PT(II)
              P2=PT(JJ)
              ANU1=ANB(II)
              ANU2=ANB(JJ)      
c             |GAMMA> inputs
              PO1=(-1)**NNO(II)
              DG1=0.
              IF(PO1.NE.1) DG1=GAMDE
              KG1=(3-PO1)/2
              NG1I=NNG(II)+1
              ANG1=ANG(NG1I,KG1)
              CD1=CD(NG1I,KG1)  
              !----
              PO2=(-1)**NNO(JJ)
              DG2=0.
              IF(PO2.NE.1) DG2=GAMDE
              KG2=(3-PO2)/2
              NG2I=NNG(JJ)+1
              ANG2=ANG(NG2I,KG2)
              CD2=CD(NG2I,KG2)                
c             matrix elements calculations
              NNT=1
              CALL OVLAB
              FOLAB1=FOLAR ! <b1|b/b0|b2>
              NNT=2
              CALL OVLAB
              FOLAB2=FOLAR ! <b1|(b/b0)^2|b2>
              NNT=0
              CALL OVLAB
              FOLAB0=FOLAR ! <b1|b2>
              NNTG=1
              CALL OVLAG
              FOLAG1=FOLAG ! <g1|g-g0|g2>
              NNTG=2
              CALL OVLAG
              FOLAG2=FOLAG ! <g1|(g-g0)^2|g2>
              !NNTG=0
              !CALL OVLAG
              !FOLAG0=FOLAG ! <g1|g2>
              NNTG=15
              CALL OVLAG
              FOLAG15=FOLAG ! <g1|sin(g)|g2>
              NNTG=16
              CALL OVLAG
              FOLAG16=FOLAG ! <g1|1-cos(g)|g2>
              !NNTG=17
              !CALL OVLAG
              !FOLAG17=FOLAG*FOLAB0 ! <g1|1/sin^2(g)|g2>
              
              IF (JU(II).EQ.JU(JJ).AND.NNB(II).EQ.NNB(JJ)
     *         .AND.NNO(II).EQ.NNO(JJ)) THEN
                  IS=JU(II)
                  NO=NNO(II)
                  KMI=1-(-1)**(IS+NO)
                  KMA=2*(IS/2)
                  NK=(KMA-KMI+2)/2
                  ITAU=NTU(II)
                  ITAU2=NTU(JJ)
                  NROOT=NNG(II)
                  NROOT2=NNG(JJ)
                  !CALL GAMIT12
                  FOLAG20=0!GAMMAT
                  FOLAG22=0!GAMMAT2
              ELSE
                  FOLAG20=0
                  FOLAG22=0
              END IF
              
              FOLAG20=FOLAG20+FOG(1,NROOT+1,NROOT2+1)+GAM0
              FOLAG22=FOLAG22+FOG(2,NROOT+1,NROOT2+1)**2
     *                +2*GAM0*FOG(2,NROOT+1,NROOT2+1)+GAM0**2
  
              IF (NNG(II).EQ.NNG(JJ)) THEN
                   FOLAG0=1.D0
              ELSE 
                   FOLAG0=0.D0              
              END IF
              
              IS=JU(II)
              NO=NNO(II)
              IS2=JU(JJ)
              NO2=NNO(JJ)
              ITAU=NTU(II)
              ITAU2=NTU(JJ)
              KMI=1-(-1)**(IS+NO)
              KMI2=1-(-1)**(IS2+NO2)              
              KAA=KMI+(ITAU-1)*2
              KAA2=KMI2+(ITAU2-1)*2
              
              !FOLAG1=FOLAG1*(1+(-1)**(KAA/2+KAA2/2))/2
              !FOLAG2=FOLAG2*(1+(-1)**(KAA/2+KAA2/2))/2
              !IF (ITAU.EQ.ITAU2) THEN
                DTT2=1.D0
              !ELSE
              !  DTT2=0.D0
              !END IF
              
              IF (PO1.eq.PO2) THEN
                  DB3=1.D0
              ELSE
                  DB3=0.D0
              END IF
  
               DBETEFF=(FOLAB1-FOLAB0)*FOLAG0*DB3
               GAMEFF=DTT2*(FOLAG1+GAM0*FOLAG0)*FOLAB0*DB3 !FOLAG15*FOLAB0!
               GAM2EFF=DTT2*(FOLAG2+2*GAM0*FOLAG1+GAM0**2*FOLAG0)*FOLAB0
     *                   *DB3!FOLAG16*FOLAB0!
               BET3EFF=
     *          BET3/DSQRT(1-DEXP(-2*BET3**2/AMUO**2))*
     *          FOLAB0*FOLAG0*
     *           (1-DB3)
              

               BET2SQ=(FOLAB2-2*FOLAB1+FOLAB0)*FOLAG0*DB3
               
               BET3SQ=(AMUO**2/2.D0+
     *            BET3**2/(1+PO1*DEXP(-BET3**2/AMUO**2)))*
     *                FOLAB0*FOLAG0*
     *                DB3
               
               !IF(II.ne.JJ)BET3SQ=0.D0
       
      EFFDEF(II,JJ,1:8) = (/DBETEFF,GAM2EFF/2,GAMEFF/1.41421356237309D0,
     *  BET3EFF*DCOS(ETO),BET3EFF*DSIN(ETO)/1.41421356237309D0,
     *  BET3SQ,BET2SQ,FOLAB0-1.d0/)
    
              END DO
          END DO
       END
C     ****************************************************************** 
      SUBROUTINE OVLOUT
C     ****************************************************************** 
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MAT/GAM,IS,NO
c      COMMON/INRM/AMO,BMO,CMO,BB42,GAMG,DELG
      COMMON/VEC/TM(16,16),AM(16,16),EIN(31,16),NK,KMI,ITAU
      COMMON/SHEMM/ES(40),JU(40),NTU(40),NNB(40),NNG(40),NNO(40),
     *NPI(40)
      COMMON/ENA/EN,EL(40),BET(10),NUR,NMAX,NPD,LAS
      COMMON/SHEM1/HW,AMG0,AMB0,GAM0,BET0,BET4,GAMDE,GSHAPE
     */SHEM2/GIT(40),EP0(40),EP1(40),EP2(40),EP12(40),ANG(40,2),EPB(40),
     *EPG(40),EGB(40),AIT(40,16),PT(40),FOG(2,8,8),XI(40),ANB(40),
     *CD(40,2)
      COMMON/OVLG/ANG1,ANG2,FOLAG,FOLAC,FOLAS,CD1,CD2,DG1,DG2,NNTG
     */MENU/MEJOB,MEPOT,MEHAM,MECHA,MEPRI,MESOL,MESHA,MESHO,MEHAO
     *,MEAPP,MEVOL,MEREL,MECUL,MERZZ,MERRR,MEDIS,MERIP
c     */FUDN/Y,DN,ANU,DNV,CDV
     */OVLB/X1,P1,ANU1,X2,P2,ANU2,FOLAR,NNT
c      COMMON/KG/J1,J2,M1,M2,J,M,AKG
c     */TR/B2(21,21),B4(21,21),B3(21,21)/NP/NPJ(25)
c     */OVLG1/FOLC2,FOLS2,FOLC3/TR1/B0(21,21)
     */SHAMO/BET3,ETO,AMUO,HWO,BB32,DPAR
c     */AA/ANO(40)
      COMMON/FILENAME/OUTNAME
      COMMON/GAMCOR/GAMMAT,GAMMAT2,ITAU2,NROOT2,NGAM!GAMMAT=<ITAU,NROOT|d^nT/dg^n*g^NGAM|ITAU2,NROOT2>
      
      COMMON/EFFD/EFFDEF(40,40,7)
      
      CHARACTER*256 OUTNAME
      CHARACTER*256 DEFOUTFILE

          IPOS=INDEX(OUTNAME,'.out',.false.)-1
          IF (IPOS.GT.0) OUTNAME=OUTNAME(:IPOS)
          DEFOUTFILE=trim(OUTNAME)//'.def'
c          print "(A100)" ,OUTFILENAME
          OPEN(UNIT=25,FILE=DEFOUTFILE,STATUS='UNKNOWN')
       
          IF(MEHAM.EQ.5)  PRINT 15
c          IF(MEHAM.EQ.5)  WRITE(21,15)
          IF(MEHAM.EQ.3)  PRINT 16
          IF(MEHAM.EQ.3)  WRITE(21,16)
          PRINT 14
c          WRITE(21,14)
          PRINT 10
c          WRITE(21,10)
   10     FORMAT(1X,80(1H-))
   11     FORMAT(30X,'<i|gamma|f> values')
   12     FORMAT(30X,'<i|gamma^2|f> values')
   13     FORMAT(30X,'<i|beta2|f> values')
   14     FORMAT(20X,'MATRIX ELEMENTS FOR DEFORMATIONS')
   15     FORMAT(20X,'CALCULATIONS INVOLVE  5-PARAM. MODEL')
   16     FORMAT(20X,'CALCULATIONS INVOLVE DAVYDOV-CHABAN MODEL')

         WRITE(*,"(A6,A7,A14,A14,A14,A14,A14,A14
     *            ,A14,A14,A14,A14,A14,A14)")
     *         'I','F','<Ib|b/b0|Fb>','<Ib|Fb>','<Ig|dg|Fg>',
     *         '<Ig|(dg)^2|Fg>','<Ig|Fg>'
     *         ,'db_eff/b0','sin(g)_eff','(1-cos(g))_eff','b3_eff'
     *         ,'<g_corr>','<g^2_corr>','<b_3^2>'
         WRITE(21,"(A6,A7,A14,A14,A14,A14,A14,A14
     *            ,A14,A14,A14,A14,A14,A14)")
     *         'I','F','<Ib|b/b0|Fb>','<Ib|Fb>','<Ig|dg|Fg>',
     *         '<Ig|(dg)^2|Fg>','<Ig|Fg>'
     *         ,'db_eff/b0','sin(g)_eff','(1-cos(g))_eff','b3_eff'
     *         ,'<g_corr>','<g^2_corr>','<b_3^2>'
         WRITE(25,"(A14,A14,A14,A14,A14,A14)")
     *         'db_eff','g20_eff','g22_eff'
     *         ,'b30_eff','b32_eff','<b_3^2>'

          CALL OVLOPT
         
          DO II=1,NUR
            DO JJ=II,NUR
c             |BETA> inputs
              X1= XI(II)
              X2= XI(JJ)
              P1=PT(II)
              P2=PT(JJ)
              ANU1=ANB(II)
              ANU2=ANB(JJ)      
c             |GAMMA> inputs
              PO1=(-1)**NNO(II)
              DG1=0.
              IF(PO1.NE.1) DG1=GAMDE
              KG1=(3-PO1)/2
              NG1I=NNG(II)+1
              ANG1=ANG(NG1I,KG1)
              CD1=CD(NG1I,KG1)  
              !----
              PO2=(-1)**NNO(JJ)
              DG2=0.
              IF(PO2.NE.1) DG2=GAMDE
              KG2=(3-PO2)/2
              NG2I=NNG(JJ)+1
              ANG2=ANG(NG2I,KG2)
              CD2=CD(NG2I,KG2)                
c             matrix elements calculations
              NNT=1
              CALL OVLAB
              FOLAB1=FOLAR ! <b1|b/b0|b2>
              NNT=2
              CALL OVLAB
              FOLAB2=FOLAR ! <b1|(b/b0)^2|b2>
              NNT=0
              CALL OVLAB
              FOLAB0=FOLAR ! <b1|b2>
              NNTG=1
              CALL OVLAG
              FOLAG1=FOLAG ! <g1|g-g0|g2>
              NNTG=2
              CALL OVLAG
              FOLAG2=FOLAG ! <g1|(g-g0)^2|g2>
              !NNTG=0
              !CALL OVLAG
              !FOLAG0=FOLAG ! <g1|g2>
              NNTG=15
              CALL OVLAG
              FOLAG15=FOLAG ! <g1|sin(g)|g2>
              NNTG=16
              CALL OVLAG
              FOLAG16=FOLAG ! <g1|1-cos(g)|g2>
              !NNTG=17
              !CALL OVLAG
              !FOLAG17=FOLAG*FOLAB0 ! <g1|1/sin^2(g)|g2>
              
              IF (JU(II).EQ.JU(JJ).AND.NNB(II).EQ.NNB(JJ)
     *         .AND.NNO(II).EQ.NNO(JJ)) THEN
                  IS=JU(II)
                  NO=NNO(II)
                  KMI=1-(-1)**(IS+NO)
                  KMA=2*(IS/2)
                  NK=(KMA-KMI+2)/2
                  ITAU=NTU(II)
                  ITAU2=NTU(JJ)
                  NROOT=NNG(II)
                  NROOT2=NNG(JJ)
                  !CALL GAMIT12
                  FOLAG20=0!GAMMAT
                  FOLAG22=0!GAMMAT2
              ELSE
                  FOLAG20=0
                  FOLAG22=0
              END IF
              
              FOLAG20=FOLAG20+FOG(1,NROOT+1,NROOT2+1)+GAM0
              FOLAG22=FOLAG22+FOG(2,NROOT+1,NROOT2+1)**2
     *                +2*GAM0*FOG(2,NROOT+1,NROOT2+1)+GAM0**2
  
              IF (NNG(II).EQ.NNG(JJ)) THEN
                   FOLAG0=1.D0
              ELSE 
                   FOLAG0=0.D0              
              END IF
              
              IS=JU(II)
              NO=NNO(II)
              IS2=JU(JJ)
              NO2=NNO(JJ)
              ITAU=NTU(II)
              ITAU2=NTU(JJ)
              KMI=1-(-1)**(IS+NO)
              KMI2=1-(-1)**(IS2+NO2)              
              KAA=KMI+(ITAU-1)*2
              KAA2=KMI2+(ITAU2-1)*2
              
              !FOLAG1=FOLAG1*(1+(-1)**(KAA/2+KAA2/2))/2
              !FOLAG2=FOLAG2*(1+(-1)**(KAA/2+KAA2/2))/2
              IF (ITAU.EQ.ITAU2) THEN
                DTT2=1.D0
              ELSE
                DTT2=0.D0
              END IF
              
              IF (PO1.eq.PO2) THEN
                  DB3=1.D0
              ELSE
                  DB3=0.D0
              END IF
  
               DBETEFF=(FOLAB1-FOLAB0)*FOLAG0*DB3
               GAMEFF=DTT2*(FOLAG1+GAM0*FOLAG0)*FOLAB0*DB3 !FOLAG15*FOLAB0!
               GAM2EFF=DTT2*(FOLAG2+2*GAM0*FOLAG1+GAM0**2*FOLAG0)*FOLAB0
     *                   *DB3!FOLAG16*FOLAB0!
               BET3EFF=
     *          BET3/DSQRT(1-DEXP(-2*BET3**2/AMUO**2))*
     *          FOLAB0*FOLAG0*
     *           (1-DB3)
              

               BET2SQ=(FOLAB2-2*FOLAB1+FOLAB0)*FOLAG0*DB3
               
               BET3SQ=(AMUO**2/2.D0+
     *            BET3**2/(1+PO1*DEXP(-BET3**2/AMUO**2)))*
     *                FOLAB0*FOLAG0*
     *                DB3
               
               !IF(II.ne.JJ)BET3SQ=0.D0
       WRITE(*,
     * "(I2,I1,I1,I1,I1,I3,I1,I1,I1,I1,
     *  12D14.6)")
     *  JU(II),NTU(II),NNB(II),NNG(II),NNO(II),
     *  JU(JJ),NTU(JJ),NNB(JJ),NNG(JJ),NNO(JJ),
     *  FOLAB1,FOLAB0,FOLAG1,FOLAG2,FOLAG0,
     *  DBETEFF,GAMEFF,GAM2EFF,BET3EFF,FOLAG20,FOLAG22,BET3SQ
     
       WRITE(21,
     *  "(I2,I1,I1,I1,I1,I3,I1,I1,I1,I1,
     *   12D14.6)")
     *  JU(II),NTU(II),NNB(II),NNG(II),NNO(II),
     *  JU(JJ),NTU(JJ),NNB(JJ),NNG(JJ),NNO(JJ),
     *  FOLAB1,FOLAB0,FOLAG1,FOLAG2,FOLAG0,
     *  DBETEFF,GAMEFF,GAM2EFF,BET3EFF,FOLAG20,FOLAG22,BET3SQ
      
      
      
        WRITE(25,"(6D14.6)") (EFFDEF(II,JJ,1:6))
                 
c     *  DBETEFF,GAM2EFF/2,GAMEFF/1.41421356237309D0,
c     *  BET3EFF*DCOS(ETO),BET3EFF*DSIN(ETO)/1.41421356237309D0,BET3SQ
     
              END DO
          END DO

      STOP
      END
C     ****************************************************************** 
      
C     *******************************************************
C     END of shemsofd
C     *******************************************************
