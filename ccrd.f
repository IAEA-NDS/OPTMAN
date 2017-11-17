C     *******************************************************
C     Start of ccrd
C     *******************************************************
      SUBROUTINE POTET
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM3.FOR'
      X=(R-RR*DE)/AR
      IF(X.GT.23) GO TO 7
      EX1=EXP(X)
      ARC1=-1.D0/(1.D0+EX1)
      GO TO 3
    7 EX1=EXP(-X)
      ARC1=-EX1
    3 IF(WC.EQ.0.D0) GO TO 5
      X=(R-RC*DE)/AC
      XX=((R-RW*DE)/AW)**2
      IF(X.GT.23) GO TO 4
      EX1=EXP(X)
      BRC1=-ALF/(1.D0+EX1)-(1.D0-ALF)*EXP(-XX)
      GO TO 5
    4 EX1=EXP(-X)
      BRC1=-ALF*EX1-(1.D0-ALF)*EXP(-XX)
    5 X=(R-RD*DE)/AD
      IF(X.GT.23) GO TO 6
      EX2=EXP(X)
      BRC2=-4.D0*EX2/(1.D0+EX2)/(1.D0+EX2)
      GO TO 2
    6 EX2=EXP(-X)
      BRC2=-4.D0*EX2
    2 CONTINUE
      VP=ARC1*VR+BRC2*VD+BRC1*VRDC
      WP=BRC1*WC+BRC2*WD
      RETURN
      END
C     *******************************************************
      SUBROUTINE ASFUT
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM10.FOR'
      INCLUDE 'PRIVCOM13.FOR' 
      INCLUDE 'PRIVCOM14.FOR'  
      AMI=939.56536D0
       IF(MECHA.EQ.1) AMI=938.272029D0
       
       IF(MECHA.EQ.1.AND.IDINT(ASP+0.1D0).EQ.1)AMI=1875.612859D0 
       IF(MECHA.EQ.1.AND.IDINT(ASP+0.1D0).EQ.1)ANEU=2.013553212712D0 
       
      REL=(EN+AMI)/AMI
      IF(MEREL.EQ.0) REL=1.D+0
      RELPOT=1.d0
      IF(REL.NE.1.) RELPOT=2.D+0*(EN+AMI)/(2.D+0*AMI+EN)

      EIRELM=SQRT(REL**2*AT**2+2.D+0*REL*ANEU*AT+ANEU**2)
      EIRELT=AT*SQRT(REL**2*ANEU**2+2.D+0*REL*ANEU*AT+AT**2)
      DEREL=SQRT(2.D+0*REL*ANEU*AT+AT**2+ANEU**2)
      WW=4.8257984E-2*EIRELM*EIRELT/(ANEU*EIRELM+EIRELT)/DEREL*(ANEU/
     *1.008664924)
      IF(MEREL.EQ.1.OR.MEREL.EQ.3) WW=WW*RELPOT


      LLMA1=LLMA+1
      LMA2=0.5252113*SQRT(VR*WW)*RR+4.D+0
      LMA1=1.65*(WNK(1)*RK+1.83)+4.D+0
      IF(LMA1.GT.LLMA1) LMA1=LLMA1
      DO 1 I=1,NUR
      COPH(I,1)=0.
      X=WNK(I)*(RK-STEP)
      
      ETA=CONZ/WNK(I)
      IF(MECHA.EQ.1.AND.NCA(1).NE.NCA(I)) ETA=0.D0
       

      IF(ETA.GT.400.D0) THEN
        IF (MEPRI.LT.98) PRINT 999, I,ETA
        WRITE(21,999) I,ETA
        pause 'WARNING, SEE OUTPUT FILE !!!'
      ENDIF
     
  999 FORMAT(9X,'WARNING! INCIDENT ENERGY IS TOO CLOSE TO',I3,3X,'LEVEL'
     */10X,'IT MAY CAUSE ERROR IN CALCULATED VALUES!!!!!!'//
     *14X,'YOU CAN MOVE ENERGY BY 0.01MEV,',4X,'ETA=',E12.7,2E12.5)

      IF(ETA.NE.0.D0) CALL COPHA
      COPH(I,1)=COPH(40,1)
      IF(ETA.EQ.0.D0) COPH(I,1)=0.D0
      IF(I.GT.NMAX) GO TO 2
      CALL BENEC
      GO TO 3
    2 LMA1=LMA2+1
      IF(LMA1.GT.LLMA1) LMA1=LLMA1
       WSTEP=0.D0
      CALL BESIM
    3 DO 4 K=1,LMA1
      IF(K.GE.LMA1) GO TO 8
      COPH(I,K+1)=COPH(I,K)+ATAN(ETA/K)
    8 FBR1(I,K)=FBR(K)
      FBI1(I,K)=FBI(K)
      FNR1(I,K)=FNR(K)
    4 FNI1(I,K)=FNI(K)
      X=WNK(I)*(RK+STEP)
      IF(I.GT.NMAX) GO TO 5
      CALL BENEC
      GO TO 6
    5 WSTEP=WNK(I)*2.D0*STEP
      CALL BESIM
    6 DO 7 K=1,LMA1      
      IF(K.GE.LMA1) GO TO 9
      COPH(I,K+1)=COPH(I,K)+ATAN(ETA/K)
    9 FBR2(I,K)=FBR(K)
      FBI2(I,K)=FBI(K)
      FNR2(I,K)=FNR(K)
    7 FNI2(I,K)=FNI(K)
    1 CONTINUE
      COPH0=COPH(1,1)
      DO 10 I=1,20
      DO 10 K=1,90
   10 COPH(I,K)=COPH(I,K)-COPH0
      RETURN
      END
C     *******************************************************
      SUBROUTINE BENEC
C     *******************************************************
C***  FBR, FNR - F AND G ASSYMPTOTIC FUNCTIONS FOR OPENED CHANNELS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM14.FOR'  
      DIMENSION FC(91),FCP(91),GC(91),GCP(91)
      IF(ETA.EQ.0.D0) GO TO 3
      MINL=0
      MAXL=LMA1-1
      RHS=X
      ETS=ETA
      ACCUR=10.D-14
c      STEPC=100
      STEPC=999.0
      LN=1
      NUMBR=3
c      CALL RCWFN(RHS,ETS,MINL,MAXL,FC,FCP,GC,GCP,ACCUR,STEPC)
      CALL RCWF(RHS,ETS,MINL,MAXL,FC,FCP,GC,GCP,ACCUR,STEPC,NUMBR)
      GO TO 4
    3 LN=3
      FNI(1)=0.D0
      FNI(2)=0.D0
      FBI(2)=0.D0
      FBI(1)=0.D0
      FBR(1)=SIN(X)
      FNR(1)=COS(X)
      FBR(2)=SIN(X)/X-COS(X)
      FNR(2)=COS(X)/X+SIN(X)
      IF(LMA1.LT.3) GO TO 2
    4 DO 1 K=LN,LMA1
      FBI(K)=0.D0
      FNI(K)=0.D0
      IF(ETA.EQ.0.D0) GO TO 5
      FBR(K)=FC(K)
      FNR(K)=GC(K)
      GO TO 1
    5 K1=K-1
      K2=K-2
      XK=(K+K-3)/X
      FNR(K)=XK*FNR(K1)-FNR(K2)
      FBR(K)=XK*FBR(K1)-FBR(K2)
    1 CONTINUE
    2 RETURN
      END
C     *******************************************************
      SUBROUTINE BESIM
C     *******************************************************
C***  FBR, FNR - F AND G ASSYMPTOTIC FUNCTIONS FOR CLOSED CHANNELS,
C***  SUCH, THAT G+iF IS DECREASING SPHERICAL WAVE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM14.FOR'  
      IF(ETA.GT.0.D0) GO TO 4
      EX1=EXP( X)
      EX2=1./EX1
      FBR(1)=0.D0
      FBI(1)=(EX1+EX2)/2.
      FNR(1)=(EX1-EX2)/2.
      FNI(1)=0.D0
      FBR(2)=FBI(1)/X-FNR(1)
      FBI(2)=0.D0
      FNR(2)=0.D0
      FNI(2)=FBI(1)-FNR(1)/X
      IF(LMA1.LT.3) GO TO 2
      DO 1 K=3,LMA1
      K1=K-1
      K2=K-2
      XK=(K+K-3)/X
      FBR(K)=XK*FBI(K1)-FBR(K2)
      FBI(K)=-XK*FBR(K1)-FBI(K2)
      FNR(K)=XK*FNI(K1)-FNR(K2)
    1 FNI(K)=-XK*FNR(K1)-FNI(K2)
    2 RETURN
    4 X=X-WSTEP
      CALL BESIMC
      RETURN
      END
C*******************************************************************
      SUBROUTINE WITTFU
C*******************************************************************
C
C     Calculates Whittaker functions W(-ETA, L+1/2, 2X) for
C     L=0 to L=LMA1-1, and X>0., using integral presentation form
C     See: M. Abramowitz and I. Stegun " Handbook of Mathematical
C      Functions" 1964, formulas 13.1.33, 13.2.5
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM14.FOR'  
      DO 1 I=1,LMA1
      ALW=I-1.
      AUP=0.D+00
      EPS=1.D-07
      ST=2.*I
      IF(I+ETA.GT.ST) ST=I+ETA
      WIT=0.
      EPSIN=EPS
    8 ALOW=AUP
      AUP=ALOW+ST
      CALL SIMPSW
      WIT=WIT+SUMCUR
      IF(SUMCUR/WIT-EPS) 13,13,18
   18 EPSIN=EPS*WIT/SUMCUR
      GO TO 8
   13 FNI(I)=WIT*EXP(-X-ETA*LOG(2.D+00*X))/DGAMMA(1.D+00+ALW+ETA)
    1 CONTINUE
      RETURN
      END
C*******************************************************************
      SUBROUTINE SIMPSW
C***********************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM14.FOR'  
      SUMIN=0.
      NN=4
      H=(AUP-ALOW)/NN
      Y=ALOW
C     IF(Y.EQ.0.) T3=0.
      T3=0.D0
      IF(Y.NE.0.) T3=Y**(ALW+ETA)*EXP(-Y)*(1.D+00+Y/2./X)**(ALW-ETA)
      Y=AUP
C     IF(Y.EQ.0.) T3=T3
      IF(Y.NE.0.) T3=T3+Y**(ALW+ETA)*EXP(-Y)*(1.D+00+Y/2./X)**(ALW-ETA)
      Y=(ALOW+AUP)/2.
      S2=0.D0
C     IF(Y.EQ.0.) S2=0.
      IF(Y.NE.0.) S2=Y**(ALW+ETA)*EXP(-Y)*(1.D+00+Y/2./X)**(ALW-ETA)
   1  S4=0.D0
      K1=NN-1
      Y=ALOW+H
      I1=1
   2  IF(Y.EQ.0.) S4=S4
      IF(Y.NE.0.) S4=S4+Y**(ALW+ETA)*EXP(-Y)*(1.D+00+Y/2./X)**(ALW-ETA)
      IF(I1-K1)3,4,4
   3  I1=I1+2
      Y=Y+2.*H
      GO TO 2
   4  SUMCUR=H/3.*(T3+2.*S2+4.*S4)
         C1=ABS(1.-SUMIN/SUMCUR)
       IF(C1-EPSIN)6,6,5
    5 SUMIN=SUMCUR
      S2=S2+S4
      H=H*0.5
      NN=NN*2
      GO TO 1
    6 SUMCUR=SUMIN+(SUMCUR-SUMIN)*1.066666
      RETURN
      END
C     ***********************************************************************
      SUBROUTINE BESIMR
C     *******************************************************
C***  FBR, FNR - F AND G ASSYMPTOTIC FUNCTIONS FOR CLOSED
C***  CHARGED PARTICLES CHANNELS,
C***  SUCH, THAT G+iF IS DECREASING SPHERICAL WAVE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM14.FOR'  
      DIMENSION F(91),FD(91),G(91),GD(91),SIGMA(91)
      RHO=X
       L=LMA1-1
       DO 5 K=1,LMA1
       FNI(K)=0.D0
       FBR(K)=0.D0
       FBI(K)=3.D0
       FNR(K)=1.D0
    5 CONTINUE
       IF(WSTEP.EQ.0.) GO TO 8
       A=1.D0
       B=2.D0
       WS10=WSTEP/3.D+02
       DO 7 I=1,300
       RHO=RHO+WS10*(I-1)
      CALL cocl(g,gd,f,fd,sigma,eta,rho,l)
       DO 4 K=1,LMA1
       A=A+A/F(K)*FD(K)*WS10
       B=B+B/G(K)*GD(K)*WS10
C      PRINT 999, FNI(K),FBR(K),x,eta,WSTEP,K
C 999 FORMAT (5E12.4,I8)
C      PAUSE 1
    4  CONTINUE
    7 CONTINUE
       DO 9 K=1,LMA1
       FBI(K)=A+B
       FNR(K)=B-A
C      PRINT 999, FNI(K),FBR(K),x,eta,WSTEP,K
C 999 FORMAT (5E12.4,I8)
C      PAUSE 1
    9  CONTINUE
    8 RETURN
      END
C     ***********************************************************************
      SUBROUTINE BESIMC
C     ***********************************************************************
C***  TWO LINEAR  INDEPENDENT CLOSED CHANNELS COULOMB FUNCTIONS
C***  EQUAL UNITY AT K*(X-STEP) AND INTEGRATED TO K*(X+STEP) USING
C***  SECOND DERIVATIVE FROM FROM COULOMB EQUATION.
C***  LINEAR INDEPENDENT, AS EQUAL IN FIRST MATCHING POINT K*(X-STEP) AND DIFFERENT
C***  AT THE NEXT INTEGRATION  POINT
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM14.FOR'  
      RHO=X
       ETA2=ETA*2.D0
       NSTE=10
       NSTE1=NSTE-1
       DO 5 K=1,LMA1
       FNI(K)=0.D0
       FBR(K)=0.D0
       FBI(K)=1.D0
       FNR(K)=1.D0
    5 CONTINUE
       IF(WSTEP.EQ.0.) RETURN
       DO 6 K=1,LMA1
       ALL1=K*(K-1.D0)
       F0=1.D0
       G0=1.D0
       F1=F0-2.D0/NSTE
       G1=G0+5.D0/NSTE
       WSN=WSTEP/NSTE
       WSN2=WSN*WSN
       X=RHO
       X2=X*X
       DO 7 I=1,NSTE1
       X=X+WSN
       F2=2.D0*F1-F0+WSN2*(ALL1/X2+ETA2/X+1.D0)*F1
       G2=2.D0*G1-G0+WSN2*(ALL1/X2+ETA2/X+1.D0)*G1
       F0=F1
       F1=F2
       G0=G1
       G1=G2
    7  CONTINUE
       FBI(K)=F2
       FNR(K)=G2
C      PRINT 999, FBI(K),FNR(K),x,eta,WSTEP,K
C 999 FORMAT (5E12.4,I8)
C      PAUSE 1
    6  CONTINUE
      RETURN
      END
C     ***********************************************************************
C 25/07/02                                                      
      SUBROUTINE COCL(G,GD,F,FD,SIGMA,ETA,RHO,L)                        
C CLOSED CHANNEL DECREASING COULOMB FUNCTIONS                           
C INPUT VARIABLES: ETA: COULOMB PARAMETER; ETA >= 0                     
C                  RHO: |K|*R VALUE                                     
C                  L:    MAXIMUM L VALUE                                
C OUTPUT VARIABLES: SIGMA(I)=0 FOR I = 1 TO L+1                         
C       F(I):  DECREASING SOLUTION AT (ETA,RHO) FOR I = 1 TO L+1        
C       FD(I): DERIVATIVE OF F(I) FOR I = 1 TO L+1                      
C       G(I):  INCREASING SOLUTION FOR I = 1 TO L+1                     
C       GD(I): DERIVATIVE OF F(I) FOR I = 1 TO L+1                      
C THE FUNCTIONS ARE RENORMALISED TO 1 AND SUCH THAT  F*GD - G*FD = 1    
C***********************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      LOGICAL IFEQAL                                                   
      DIMENSION F(1),FD(1),G(1),GD(1),SIGMA(1),S(7),T(3)                
      DATA S /7*0.D0/                                                  
      IF (ETA.LT.0.D0.OR.RHO.LE.0.D0) GO TO 11                          
      LP1=L+1                                                           
      IF (ETA.LT.1.D-6) GO TO 6                                         
      IS=7                                                             
      M=10.D0*RHO+1.D0                                                  
      H=M                                                               
      H=RHO/H                                                           
      RHOA=10.D0*(ETA+1.D0)                                             
      IFEQAL=RHOA.LT.RHO                                                
      IF (IFEQAL) RHOA=RHO                                             
      M=RHOA/H+0.5D0                                                   
      RHOA=H*M                                                         
      IF (.NOT.IFEQAL.AND.RHOA.LT.RHO+1.5D0*H) RHOA=RHO+2.D0*H          
C EXPANSION IN POWERS OF 1/RHOA                                         
    1 C=1.D0/RHOA                                                      
      A=1.D0                                                           
      B=A                                                               
      D=0.D0                                                           
      DO 2 M=1,26                                                      
      AM=M                                                              
      A=-A*0.5D0*(ETA+AM-1.D0)*(ETA+AM)*C/AM                           
      B=B+A                                                            
    2 D=D-A*AM*C                                                       
      F(1)=1.D0                                                         
      FD(1)=D/B-1.D0-ETA/RHOA                                           
      IF (IFEQAL) GO TO 7                                              
      S(IS)=B                                                          
      IF (IS.NE.7) GO TO 3                                              
      IS=6                                                              
      RHOA=RHOA+H                                                      
      S(7)=S(7)*EXP(H-ETA*LOG(1.D0-H/RHOA))                          
      GO TO 1                                                          
C BACKWARD INTEGRATION                                                  
    3 A=2.D0+1.D0/1.2D0*H*H                                            
      B=1.D0/6.D0*H*ETA                                                
      C=1.D0-1.D0/12.D0*H*H                                            
      M1=RHOA/H-0.5D0                                                 
      M2=RHO/H-1.5D0                                                   
      AM=M1                                                            
      T(2)=B/(AM+1.D0)                                                 
      T(3)=B/AM                                                        
      JS=M1                                                            
      DO 5 IS=M2,M1                                                    
      DO 4 I=1,6                                                       
    4 S(I)=S(I+1)/S(7)                                                 
      T(1)=T(2)                                                         
      T(2)=T(3)                                                        
      AM=JS-1                                                           
      T(3)=B/AM                                                        
      S(7)=((A+10.D0*T(2))*S(6)-(C-T(1))*S(5))/(C-T(3))                
    5 JS=JS-1                                                          
      F(1)=1.D0                                                        
      FD(1)=(1.D0/60.D0*(S(1)-S(7))+0.15D0*(S(6)-S(2))+0.75D0*(S(3)-S(5)
     1))/(H*S(4))                                                       
      GO TO 7                                                           
    6 F(1)=1.D0                                                        
      FD(1)=-1.D0                                                      
C RECURRENCE FOR L > 0
C RECURRENCE FOR L >                                                   
    7 C=1.D0/RHO                                                       
      IF (L.LE.0) GO TO 9                                              
      DO 8 M=1,L                                                       
      AM=M                                                             
      A=ETA/AM                                                         
      B=A+C*AM                                                         
      F(M+1)=1.D0                                                       
    8 FD(M+1)=(A*A-1.D0)/(B-FD(M))-B                                   
    9 DO 10 M=1,LP1                                                    
      G(M)=1.D0                                                        
      GD(M)=1.D0+FD(M)                                                 
   10 SIGMA(M)=0.D0                                                    
      RETURN                                                           
   11 WRITE (6,1000) ETA,RHO                                           
      STOP                                                             
 1000 FORMAT (17H COCL  ***  ETA =,1P,D13.5,8H,   RHO=,D13.5,26H   ARGUM
     1ENT OUT OF RANGE    )                                             
      END                                                               


C     *******************************************************
c 29/05/86  ibm version                                         
      subroutine coclold(g,gd,f,fd,sigma,eta,rho,l)                     
c closed channel decreasing coulomb functions                           
c input variables: eta: coulomb parameter; eta >= 0                     
c                  rho: |k|*r value                                     
c                  l:    maximum l value                                
c output variables: sigma(i)=0 for i = 1 to l+1                         
c       f(i):  decreasing solution at (eta,rho) for i = 1 to l+1        
c       fd(i): derivative of f(i) for i = 1 to l+1                      
c       g(i):  increasing solution for i = 1 to l+1                     
c       gd(i): derivative of f(i) for i = 1 to l+1                      
c the functions are renormalised to 1 and such that  f*gd - g*fd = 1    
c***********************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      logical ifeqal                                                    
      dimension f(1),fd(1),g(1),gd(1),sigma(1),s(7),t(3)                
      data s /7*0.d0/                                                   
      if (eta.lt.0.d0.or.rho.le.0.d0) go to 11                          
      lp1=l+1                                                           
      if (eta.lt.1.d-6) go to 6                                         
      is=7                                                             
      m=10.d0*rho+3.d0                                                
      h=m                                                             
      h=rho/h                                                         
      rhoa=10.d0*(eta+1.d0)                                          
      ifeqal=rhoa.lt.rho                                              
      m=rhoa/h+0.5d0                                                  
      rhoa=h*m                                                       
      if (ifeqal.or.rhoa.lt.rho+1.5d0*h) rhoa=rho+2.d0*h             
c expansion in powers of 1/rhoa                                       
    1 c=1.d0/rhoa                                                     
      a=1.d0                                                          
      b=1.d0-c*eta                                                    
      f(1)=a                                                         
      fd(1)=b                                                        
      do 2 m=1,26                                                    
      am=m                                                           
      d=0.5d0*(eta+am-1.d0)*(eta+am)*c/am                             
      a=-a*d                                                          
      b=-b*d-a*c                                                     
      f(1)=f(1)+a                                                     
    2 fd(1)=fd(1)+b                                                   
      if (ifeqal) go to 7                                            
      s(is)=f(1)                                                      
      if (is.ne.7) go to 3                                             
      is=6                                                           
      rhoa=rhoa+h                                                    
      s(7)=s(7)*EXP(h-eta*LOG(1.d0-h/rhoa))                         
      go to 1                                                         
c backward integration                                                
    3 a=2.d0+1.d0/1.2d0*h*h                                           
      b=1.d0/6.d0*h*eta                                                
      c=1.d0-1.d0/12.d0*h*h                                          
      m1=rhoa/h-0.5d0                                               
      m2=rho/h-1.5d0                                                 
      am=m1                                                         
      t(2)=b/(am+1.d0)                                              
      t(3)=b/am                                                    
      js=m1                                                         
      do 5 is=m2,m1                                                 
      do 4 i=1,6                                                     
    4 s(i)=s(i+1)/s(7)                                              
      t(1)=t(2)                                                    
      t(2)=t(3)                                                    
      am=js-1                                                      
      t(3)=b/am                                                       
      s(7)=((a+10.d0*t(2))*s(6)-(c-t(1))*s(5))/(c-t(3))              
    5 js=js-1                                                        
      f(1)=1.d0                                                     
      fd(1)=(1.d0/60.d0*(s(1)-s(7))+0.15d0*(s(6)-s(2))+0.75d0*(s(3)-s(5)
     1))/(h*s(4))                                                       
      go to 7                                                          
    6 f(1)=1.d0                                                        
      fd(1)=-1.d0                                                      
c recurrence for l > 0                                                 
    7 c=1.d0/rho                                                       
      if (l.le.0) go to 9                                              
      do 8 m=1,l                                                       
      am=m                                                             
      a=eta/am                                                         
      b=a+c*am                                                         
      f(m+1)=1.d0                                                     
    8 fd(m+1)=(a*a-1.d0)/(b-fd(m))-b                                   
    9 do 10 m=1,lp1                                                    
      g(m)=1.d0                                                       
      gd(m)=1.d0+fd(m)                                                
   10 sigma(m)=0.d0                                                   
      return                                                           
   11 write (6,1000) eta,rho                                           
      stop                                                             
 1000 format (17h cocl  ***  eta =,1p,d13.5,8h,  rho= ,d13.5,27h    argu
     1ment out of range         )                                       
      end                                                               
C     *******************************************************
      SUBROUTINE QUANT
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 

      DOUBLE PRECISION jc,jj
      INTEGER ltlmax, itmp, LLLMAX
      CHARACTER*1 parc
      LOGICAL unformatted

      PARAMETER (LLLMAX=100)
      DIMENSION STL(LLLMAX)
      
      real*8, parameter :: haf = 0.5D0

      CHARACTER*20 fname
      COMMON/INOUT/fname

      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM7.FOR'  
      INCLUDE 'PRIVCOM10.FOR' 
      INCLUDE 'PRIVCOM13.FOR'      
      AMI=939.56536
      IF(MECHA.EQ.1) AMI=938.272029D0
       
      IF(MECHA.EQ.1.AND.IDINT(ASP+0.1D0).EQ.1) AMI=1875.612859D0 
      IF(MECHA.EQ.1.AND.IDINT(ASP+0.1D0).EQ.1) ANEU=2.013553212712D0 
       
C     NSPI=IDINT(ASP*2+0.001)
      NSPI=INT(ASP*2+0.001)
      NDEL=1
      IF(NSPI/2*2.EQ.NSPI) NDEL=0
      DO 876 L=1,100
      TL(L)=0.D0
      TRLJ(L,1)=0.D0
      TRLJ(L,2)=0.D0
      SSR(L)=0.D0
      SSI(L)=0.D0
 876  CONTINUE
      LMAX=1.65*(WNK(1)*RK+1.83)+3
      LMAX1=LMAX+1
      IF(LMAX.GT.LLMA) LMAX1=LLMA+1
      REL=(EN+AMI)/AMI
      IF(MEREL.EQ.0) REL=1.D+0
      IF(REL.NE.1.) RELPOT=2.D+0*(EN+AMI)/(2.D+0*AMI+EN)


      EIRELM=SQRT(REL**2*AT**2+2.D+0*REL*ANEU*AT+ANEU**2)
      EIRELT=AT*SQRT(REL**2*ANEU**2+2.D+0*REL*ANEU*AT+AT**2)
      DEREL=SQRT(2.D+0*REL*ANEU*AT+AT**2+ANEU**2)
      WW=4.8257984E-2*EIRELM*EIRELT/(ANEU*EIRELM+EIRELT)/DEREL*(ANEU/
     *1.008664924)
      IF(MEREL.EQ.1.OR.MEREL.EQ.3) WW=WW*RELPOT


      LMA2=0.5252113*SQRT(VR*WW)*RR+4.D+0
      IF(MEPRI.GE.10 .AND. MEPRI.LT.98) 
     *  PRINT 18,RK,WW,WNK(1),LMA2,LMAX
      IF(MEPRI.GE.10 .AND. MEPRI.LT.98) 
     *  WRITE(21,18)RK,WW,WNK(1),LMA2,LMAX
   18 FORMAT (1X,'RK=',D15.7,1X,'WW=',D15.7,1X,'WNK(1)=',D15.7,
     *1X,'LMA2=',I3,1X,'LMAX=',I3/)
      SF0=0.D0
      SF1=0.D0
      SF2=0.D0
      RS1=1.23*AT**(1./3.)+0.8
      WRR=WN(1)*RS1*RS1
      SEN=SQRT(EN)*3141.593
      P0=SEN
      P1=SEN*3.*WRR/(1.+WRR)
      P2=SEN*5.*WRR**2/(9.+3.*WRR+WRR**2)
      NSS=0
      DO 12 N=1,40
      CSN(N)=0.D0
   12 CONTINUE
      CST=0.D0
      CSTRT=0.D0
      CSRINT=0.D0
      
      
      IF(MEPRI.GE.10 .AND. MEPRI.LT.98) PRINT 19
      IF(MEPRI.GE.10 .AND. MEPRI.LT.98) WRITE(21,19)
   19 FORMAT(20X,'JSS    NPIS    LNJ1(I)     JNJ1(I)      NNJ1(I) ')
      DO 1 NOT=1,300
      NOT1=NOT-1
      NJ=0
      DO 3 NO1=1,2
      JSS=JO(1)+NDEL+2*NOT1*(-1)**NO1
      IF(JSS.LT.0) NJ=NJ+2
      IF(JSS.LT.0) GO TO 3
      NJ1=0
      DO 2 NCH=1,2
      NSS=NSS+1
      JS(NSS)=JSS
      N1=0
      N2=0
      N3=0
      NPIS=(-1)**NCH
      DO 4 I1=1,NUR
      LMAX=1.65*(WNK(1)*RK+1.83)+3
      IF(I1.GT.NMAX) LMAX=LMA2
      IF(LMAX.GT.LLMA) LMAX=LLMA
      JSO=JO(I1)
      NPIO=NPO(I1)
      CALL LOGMO
      IF(NN1.EQ.0) GO TO 4
      DO 5 I2=1,NN1
      IF(LNO(I2).GT.LMAX) GO TO 4
      IF(N1.GE.250) GO TO 4
      N1=N1+1
      IF(I1.EQ.1) N2=N2+1
      LNJ1(N1)=LNO(I2)
      JNJ1(N1)=JNO(I2)
      NNJ1(N1)=I1
      KPJ1(N1)=NPO(I1)
      LNJ(NSS,N1)=LNJ1(N1)
      JNJ(NSS,N1)=JNJ1(N1)
    5 NNJ(NSS,N1)=NNJ1(N1)
    4 CONTINUE
      IF(N2)13,14,13
   14 NSS=NSS-1
      NJ1=NJ1+1
      IF(NJ1.EQ.2) GO TO 7
      GO TO 2
   13 CONTINUE
      IF(KODMA.LE.0) GO TO 22
      N2P=N2+1
      N1M=N1-1
      DO 21 IM=N2P,N1M
      MLL=LNJ1(IM)
      IMP=IM+1
      DO 20 JM=IMP,N1
      IF(MLL.LE.LNJ1(JM)) GO TO 20
      MLL=LNJ1(JM)
      MLNJ=LNJ1(IM)
      MJNJ=JNJ1(IM)
      MNNJ=NNJ1(IM)
      MKPJ=KPJ1(IM)
      LNJ1(IM)=LNJ1(JM)
      JNJ1(IM)=JNJ1(JM)
      NNJ1(IM)=NNJ1(JM)
      KPJ1(IM)=KPJ1(JM)
      LNJ(NSS,IM)=LNJ1(JM)
      JNJ(NSS,IM)=JNJ1(JM)
      NNJ(NSS,IM)=NNJ1(JM)
      LNJ1(JM)=MLNJ
      JNJ1(JM)=MJNJ
      NNJ1(JM)=MNNJ
      KPJ1(JM)=MKPJ
      LNJ(NSS,JM)=MLNJ
      JNJ(NSS,JM)=MJNJ
      NNJ(NSS,JM)=MNNJ
   20 CONTINUE
   21 CONTINUE
   22 INCC=N2
C     NEXT TWO CARDS JUST TO OVERCOME THE ERROR OR  MICROSOFT FPS COMPILER
C     v.1.00, YOU CAN DELETE THEM USING OTHER COMPILERS
      IF (MEPRI.LT.98) PRINT 899,(WNK(II),II=1,NUR)
  899 FORMAT (5E12.5)
      IF(N2.GT.NCMA) INCC=NCMA
      INC(NSS)=INCC
      NCLL=N1
      IF (MEPRI.LT.98) PRINT 898, N1, NCMA
  898 FORMAT(/2X,'GENERATED CC NUMBER=',I3,3X,'ALLOWED CC NUMBER=',I3/)
      IF(N1.GT.NCMA) NCLL=NCMA
      DO 23 MNC=1,NCLL
      IF(NNJ1(MNC).LE.NMAX) N3=N3+1
   23 CONTINUE
      INCR=N3
      INR(NSS)=INCR
      IF(MEPRI.GE.10 .AND. MEPRI.LT.98) 
     *  PRINT 11,   (JSS,NPIS,LNJ1(NL),JNJ1(NL),NNJ1(NL),NL=1,NCLL)
      IF(MEPRI.GE.10 .AND. MEPRI.LT.98)
     *  WRITE(21,11)(JSS,NPIS,LNJ1(NL),JNJ1(NL),NNJ1(NL),NL=1,NCLL)
   11 FORMAT (15(I5,2I3,I3,I2))
      IF(LAS.EQ.0) GO TO 100
      IF(MEPOT.EQ.1.AND.NPD.EQ.0) GO TO 100
C     CALL KNDIT
 100  CALL CMATC
      CONG=(JSS+1)/(JO(1)+1.)
      CON=12.566372*CONG/WN(1)/100./(NSPI+1.D0)
      
C      PRINT 787, JSS, JO(1),NSPI, CONG,CON
C  787 FORMAT(3I5,2E12.5)
C      PAUSE
      
      CSS=CON*CSS
      IF(MEPRI.LT.98) PRINT 999, JSS,CSS
C     IF(CSS.LT.1.E-4) NJ=NJ+1

      IF(NSS.GT.1) GO TO 677
      FIR=WNK(1)*RS1
      CSEL=FIR**2
           
      FIR=FIR-ATAN(FIR)
      CSEL=CSEL+3.D0*FIR**2
      
      FIR=FIR-ATAN(3.D0*FIR/(3.D0-FIR**2))
      CSEL=CSEL+5.D0*FIR**2 
           
      FIR=FIR-ATAN(FIR*(15.D0-FIR)**2/
     *(15.D0-6.D0*FIR**2))
      CSEL=(CSEL+7.D0*FIR**2)*CON/CONG*(NSPI+1.D0) 
      
 677  CSTJPI=0.D0
      DO 8 IC=1,INCC
      CSTJPI=CSTJPI+ABS(CI(IC,IC)*CON)
 
      CSTR=0.D0
      CSNR=0.D0
      CSNRIN=0.D0
      CSRER=0.D0 
      TDR1=0.D0 

C     RCN June 2014
      ANORM = 0.D0
C     RCN June 2014
      IF(NRESN.EQ.0) GO TO 201     
     
      EBOUND=6.D0
      CRO=EXP(-2.D0*(AT/8.D0*EBOUND)) 
      CROE=EXP(2.D0*(AT/8.D0*(EN-EBOUND))) 
      ANORMM=(1.D0-CRO)/(1.D0-2.D0*CRO+CROE)
      ANORM=1.D0-ANORMM
      ANORM2=ANORM**2
c     IF(MEPRI.LT.98) PRINT 919, EN,ANORM,ANORMM,EBOUND,CRO,CROE,CONG
c     PAUSE 777
      ANORMM=1.D0
      ANORM=1.D0
         
      DO 202 I=1,NRESN
           
      IF(JCON(I).NE.JSS.OR.LON(I).NE.LNJ1(IC)
     *.OR.JMN(I).NE.JNJ1(IC))GO TO 202
      
      ERI=ERN(I)
      REL=(ABS(ERI)+AMI)/AMI
      IF(MEREL.EQ.0.OR.ERI.LE.0.) REL=1.D+0
      IF(MEREL.NE.0.AND.ERI.GT.0.) ERIC=(REL**2-1.D+0)*AMI/2.D+0
      APAN=ANEU/1.00866491600D0 
      WNKRI=0.219677*SQRT(ABS(ERIC))*AT/SQRT(AT**2+2.D+0*ANEU*REL*AT+
     *ANEU**2)*SQRT(APAN)
      WNKR=WNKRI*RS1
      WNKR1=WNK(1)*RS1
      WRRRE=WNKR**2
      RP00=WNKR1/WNKR
      IF(LON(I).EQ.0) PR0=RP00
      IF(LON(I).EQ.1) PR0=RP00*WRR/(1.+WRR)/WRRRE*(1.D0+WRRRE)
      IF(LON(I).EQ.2) PR0=RP00*WRR**2/(9.+3.*WRR+WRR**2)/
     *WRRRE**2*(9.+3.*WRRRE+WRRRE**2)
      IF(LON(I).EQ.3) PR0=RP00*WRR**3/(225.+45.*WRR+6.*WRR**2+WRR**3)/
     *WRRRE**3/(225.+45.*WRRRE+6.*WRRRE**2+WRRRE**3)
     
      IF(LON(I).EQ.0) FIR=WNKR1
      IF(LON(I).EQ.1) FIR=WNKR1-ATAN(WNKR1)
      IF(LON(I).EQ.2) FIR=WNKR1-ATAN(3.D0*WNKR1/(3.D0-WNKR1**2))
      IF(LON(I).EQ.3) FIR=WNKR1-ATAN(WNKR1*(15.D0-WNKR1)**2/
     *(15.D0-6.D0*WNKR1**2))
     
      GNNE=GNN(I)*PR0
      GTOTN=GREN(I)+GNNE
      IF(NEL(I).NE.0)GTOTN=GREN(I)+GNNE+GNN(NEL(I))*PR0
      ADEN=(EN-ERN(I))**2+GTOTN**2/4.D0
      CMRR=GNNE*GTOTN/4.D0/ADEN
      CMRI=-GNNE*(EN-ERN(I))/2.D0/ADEN
      
        

      CR(IC,IC)=ANORM*CR(IC,IC)+CMRR*ANORMM
      CI(IC,IC)=ANORM*CI(IC,IC)+CMRI*ANORMM
      GO TO 979 
      CR(IC,IC)= ANORM*CR(IC,IC)
      CI(IC,IC)= ANORM*CI(IC,IC)
      
      
       
      CSTR=CSTR+GNNE*GTOTN/ADEN*CON/4.D0
      CSNR=CSNR+GNNE*GNNE/ADEN*CON/4.D0
      CSNRIN=CSNRIN-(GNNE*GTOTN*SIN(FIR)**2-(EN-ERN(I))*GNNE*
     *SIN(2.D0*FIR))/ADEN*CON/2.D0
      CSRER=CSRER+GNNE*GREN(I)/ADEN*CON/4.D0
      TDR1=CSRER*4.D0/CON*ANORMM
       
  979 IF(MEPRI.LT.98) THEN
        PRINT 919,CSTR,CSNR,CSRER,CSNRIN
        PRINT 919,CON,GNNE,GTOTN,ADEN
        PRINT 919, CMRI, CMRR,CI(IC,IC),CR(IC,IC)
      ENDIF
  919 FORMAT (5E12.5)
  202 CONTINUE    
  201 CONTINUE           
      
      CSTRT=CSTRT+CSTR
      CSRINT=CSRINT+CSNRIN
c   8 nspi
CST+CI(IC,IC)*CON+CSTR+CSNRIN
    8 CST=CST+CI(IC,IC)*CON 
      

      IF(MECHA.EQ.0) EPC=3.E-4
      IF(MECHA.NE.0) EPC=3.E-9
      IF(CSS.LT.EPC.AND.CSTJPI.LT.EPC) NJ=NJ+1
      IF(MEPRI.LT.98) PRINT 123,CST
  123 FORMAT (6X,F20.10)
      DO 9 I1=1,NMAX
      CM(I1)=0.
      DO 10 IC=1,INCC
      DO 10 IR=1,INCR
      IF(I1.EQ.NNJ1(IR)) CM(I1)=CM(I1)+CR(IC,IR)**2+CI(IC,IR)**2
   10 CONTINUE
    9 CSN(I1)=CSN(I1)+CM(I1)*CON
      DO 15 IC=1,INCC
      JNUM=(3+JNJ1(IC)-2*LNJ1(IC))/2
      CI2=CI(IC,IC)
      T1=CI2
C     TL,SR,SI,SFI -FOR SPIN=1/2, CHANGE 2 by NSPI+1 IN FORMULAE
      IF(LNJ1(IC).EQ.0) SF0=SF0+CI2*CONG*2.D0/(NSPI+1.D0)
      IF(LNJ1(IC).EQ.1) SF1=SF1+CI2*CONG*2.D0/(NSPI+1.D0)
      IF(LNJ1(IC).EQ.2) SF2=SF2+CI2*CONG*2.D0/(NSPI+1.D0)
      DO 16 IR=1,INCR
      ANOIR=1.D0
      IF(IC.NE.IR) ANOIR=ANORM
      ANOIR=1.D0
      CR1=CR(IC,IR)*ANOIR
      CI1=CI(IC,IR)*ANOIR
      T1=T1-CR1*CR1-CI1*CI1
      IF(LNJ1(IC).EQ.0) SF0=SF0-(CR1*CR1+CI1*CI1)*CONG*2.D0/(NSPI+1.D0)
      IF(LNJ1(IC).EQ.1) SF1=SF1-(CR1*CR1+CI1*CI1)*CONG*2.D0/(NSPI+1.D0)
      IF(LNJ1(IC).EQ.2) SF2=SF2-(CR1*CR1+CI1*CI1)*CONG*2.D0/(NSPI+1.D0)
   16 CONTINUE
      TD1(NSS,IC)=T1*4.D0 +TDR1
      IF(MEPRI.LT.98) PRINT 919, T1, TDR1, TD1(NSS,IC)
C     PAUSE 100

      TL(LNJ1(IC)+1)=TL(LNJ1(IC)+1)+(JSS+1)*4.*T1
      TRLJ(LNJ1(IC)+1,JNUM)=TRLJ(LNJ1(IC)+1,JNUM)+(JSS+1)*4.*T1
      SSR(LNJ1(IC)+1)=SSR(LNJ1(IC)+1)+(1.-2.*CI(IC,IC))*(JSS+1)
      SSI(LNJ1(IC)+1)=SSI(LNJ1(IC)+1)+2.*CR(IC,IC)*(JSS+1)
   15 CONTINUE
      
c      write(21,*) '*** C-matrix'
c      write(21,'(A,I2,A,<INCC>(A,I2,A,I2,A,I2,A))')
c     * '*** J=',JSS,'/2,',(' l_in=', LNJ1(Iiic),', j_in=',
c     *   JNJ1(Iiic),'/2, N=', NNJ1(Iiic),';' ,Iiic=1,INCC)
c      do Iiir=1,INCR
c        write(21,'(A,I2,A,I2,A,I2,A,<INCC>(ES13.5,A,ES13.5))') 'l_out=',
c     *     LNJ1(Iiir),', j_out=',JNJ1(Iiir),'/2, N=', NNJ1(Iiir),'; ',
c     *     (CR(Iiic,Iiir),'+i*',CI(Iiic,Iiir) ,Iiic=1,INCC)
c      enddo      
c      write(21,*) '*** C-matrix end'
      
      
      IF(NSS.GE.NSMA) GO TO 7
    2 CONTINUE
C     IF(NJ.EQ.2*(NSPI+1)) GO TO 7
      IF(NJ.EQ.4) GO TO 7
      IF(NOT1.EQ.0) GO TO 1
    3 CONTINUE
    1 CONTINUE
    7 NJ=NSS
      IF(NRESN.EQ.0) GO TO 222

      IF(MEPRI.LT.98) THEN
        PRINT *,' RENORMALIZATON FROM RESONANCES'
        PRINT 919, ANORM,CST,CSN(1),CSTRT,CSRINT
      ENDIF
      CST=(CST-CSN(1))*ANORM+CSN(1)+CSTRT+CSRINT
      CSN(1)=CSN(1)+CSRINT
      ltlmax = LLLMAX
      WRITE (21,*) 
      WRITE (21,*) ' RENORMALIZATON FROM RESONANCES'
      WRITE (21,919) ANORM,CST,CSN(1),CSTRT,CSRINT
C     PAUSE 999
C 220  IF(MEJOB.EQ.2) GO TO 17
      
c      IF(NSPI.GT.1) GO TO 17
 222  CONTINUE
      
      IF(MEPRI.LT.98.AND.MEJOB.EQ.1) THEN
      PRINT 987
      WRITE(21,987)
 987  FORMAT (//1X,'ORB. MOMENT',14X,'TRANSITIONS',12X,'SR',18X,'SI'/)
      ENDIF

      DO 888 L=1,LMAX1
      LL=L-1
      TLL=TL(L)/(NSPI+1.D0)/(JO(1)+1.)/(2.*LL+1.)
C     TRLJ(L,1)=TRLJ(L,1)/2./(JO(1)+1.)/(2.*LL+1)
C     TRLJ(L,2)=TRLJ(L,2)/2./(JO(1)+1.)/(2.*LL+1)
      IF(LL.EQ.0) TRLJ(1,1)=TLL
      if(tll.gt.1.d-15) ltlmax = L
      IF(LL.NE.0) THEN
        TRLJ(L,1)=TRLJ(L,1)/2./(JO(1)+1.)/(LL)
        if(TRLJ(L,1).gt.1.e-15) ltlmax = L
      ENDIF
      TRLJ(L,2)=TRLJ(L,2)/2./(JO(1)+1.)/(LL+1)
      if(TRLJ(L,2).gt.1.e-15) ltlmax = L
      SRL=SSR(L)/(NSPI+1.D0)/(JO(1)+1.)/(2.D0*LL+1.D0)
      SIL=SSI(L)/(NSPI+1.D0)/(JO(1)+1.)/(2.D0*LL+1.D0)
      TRL(L)=TLL
      IF(MEPRI.LT.98.AND.MEJOB.EQ.1) THEN
        PRINT 999,LL,TLL,SRL,SIL
        WRITE(21,999)LL,TLL,SRL,SIL
      ENDIF
c
c  [GN] 11/20015 : Add SPRT+ESW calculations to get the following parameters
c
c .SFR0,1,2       : strenght function form SPRT and ESW
c .RRPRIME0,1,2  : effective radius
c     
c       print *,'!!! BEFORE SPRT!!!: ', LL,TLL,SRL,SIL,AT,EN
       IF(LL.eq.0)then
         SFR0=0.0
       RRPRIME0=0.0
       CALL SPRT(LL,TLL,SRL,SIL,AT,EN,SFR0,RRPRIME0)
       ELSEIF(LL.eq.1)then
         SFR1=0.0
       RRPRIME1=0.0
       CALL SPRT(LL,TLL,SRL,SIL,AT,EN,SFR1,RRPRIME1) 	 
       ELSEIF(LL.eq.2)then
         SFR2=0.0
       RRPRIME2=0.0
       CALL SPRT(LL,TLL,SRL,SIL,AT,EN,SFR2,RRPRIME2) 	 
       ENDIF
c       print *,'!!! AFTER SPRT!!!: ', SFR0,RRPRIME0,SFR1,RRPRIME1
c     *     ,SFR1,RRPRIME1
c [GN]        
 888  CONTINUE
      ltlmax = Min(ltlmax,LLLMAX-1)
 
 220  IF(MEJOB.EQ.2) GO TO 17
C
C     Formatting the output for reaction codes (e.g. EMPIRE,GNASH)
C     RCN
    
C
      IF(MEPRI.EQ.98) THEN   
C
C      EMPIRE (ecis06 formatted output)
C
#ifdef FORMATTEDOUT
      unformatted = .false.  
#else
      unformatted = .true.  
#endif          

       numbtl = 0
       DO L=1,ltlmax
        IF(NPO(1)*(-1)**(L-1).EQ.+1) THEN 
          IF(L.NE.1 .AND. TRLJ(L,1).GT.0.) numbtl = numbtl + 1
          IF(TRLJ(L,2).GT.0.) numbtl = numbtl + 1
        ENDIF
        IF(NPO(1)*(-1)**(L-1).EQ.-1) THEN 
          IF(L.NE.1 .AND. TRLJ(L,1).GT.0.) numbtl = numbtl + 1
          IF(TRLJ(L,2).GT.0.) numbtl = numbtl + 1
        ENDIF
       ENDDO
C----------------
C      from ecis06
 1006  FORMAT ('<TLJ     >',F10.2,1P,D20.8,0P,F10.2,2I5)                CAL1-427
 1007  FORMAT (1X,F9.1,4X,A1,1X,I4)                                     CAL1-428
 1008  FORMAT (1X,I2,I6,F9.1,2X,1P,D18.8,0P)                            CAL1-429
C----------------
       if (unformatted) then
           open(unit=92,file=TRIM(fname)//'.TLJ', form='unformatted')
       else
           open(unit=92,file=TRIM(fname)//'.TLJ')
       endif

C      WRITE(92,'(10H<TLJ     >,F10.2,F10.5,F10.2,2I5)') 
       IF(.NOT.unformatted)
     & WRITE(92,1006)ANEU,EN,AT,NINT(0.5*JO(1)),numbtl

       parc = '+'
       numbtl = 0

       DO L=1,ltlmax
        LL = L-1 
        IF(NPO(1)*(-1)**LL.EQ.-1) CYCLE
        IF(L.NE.1 .AND. TRLJ(L,1).GT.0) THEN
          numbtl = numbtl + 1
          xz = dble(ll) - haf
          IF(unformatted) then 
            WRITE(92)                  xz, parc,1
            WRITE(92)      1, ll, xz, TRLJ(L,1)
          ELSE
C           WRITE(92,1007)         L-1-0.5,'+',1
C           WRITE(92,1008) 1, L-1, L-1-0.5, TRLJ(L,1)
            WRITE(92,1007)             xz, parc,1
            WRITE(92,1008) 1, ll, xz, TRLJ(L,1)
          ENDIF

        ENDIF
        IF(TRLJ(L,2).GT.0) THEN
          numbtl = numbtl + 1
          xz = dble(ll) + haf
          IF(unformatted) then 
            WRITE(92)                  xz, parc,1
            WRITE(92)      1, ll, xz, TRLJ(L,2)
          ELSE
C           WRITE(92,1007)         L-1+0.5,'+',1
C           WRITE(92,1008) 1, L-1, L-1+0.5, TRLJ(L,2)
            WRITE(92,1007)             xz, parc,1
            WRITE(92,1008) 1, ll, xz, TRLJ(L,2)
          ENDIF 
        ENDIF
       ENDDO

       parc='-'

       DO L=1,ltlmax
        LL = L-1 
        IF(NPO(1)*(-1)**LL.EQ.+1) CYCLE
        IF(L.NE.1 .AND. TRLJ(L,1).GT.0.) THEN
          numbtl = numbtl + 1
          xz = dble(ll) - haf
          IF(unformatted) then 
            WRITE(92)                  xz, parc,1
            WRITE(92)      1, ll, xz, TRLJ(L,1)
          ELSE 
C           WRITE(92,1007)         L-1-0.5,'-',1
C           WRITE(92,1008) 1, L-1, L-1-0.5, TRLJ(L,1)
            WRITE(92,1007)             xz, parc,1
            WRITE(92,1008) 1, ll, xz, TRLJ(L,1)
          ENDIF
        ENDIF
        IF(TRLJ(L,2).GT.0.) THEN
          numbtl = numbtl + 1
          xz = dble(ll) + haf
          IF(unformatted) then 
            WRITE(92)                  xz, parc,1
            WRITE(92)      1, ll, xz, TRLJ(L,2)
          ELSE
C           WRITE(92,1007)         L-1+0.5,'-',1
C           WRITE(92,1008) 1, L-1, L-1+0.5, TRLJ(L,2)
            WRITE(92,1007)             xz, parc,1
            WRITE(92,1008) 1, ll, xz, TRLJ(L,2)
          ENDIF
        ENDIF
       ENDDO
C      CLOSE(92)

       Stl = 0.d0

       IF(unformatted) then 
         OPEN(45,STATUS='old',FILE=TRIM(fname)//'.tlj',
     &         ERR=1200,form='unformatted')
       ELSE
         OPEN(45,STATUS='old',FILE=TRIM(fname)//'.tlj', ERR=1200)
         READ(45,*,END=1200)    ! To skip first line <TLJs.> ..
       ENDIF 

C------JC,ParC is the channel spin and parity
C------nceq is the number of coupled equations
1100   IF(unformatted) then 
           READ(45,END = 1200)      jc, parc, nceq  ! ecis06
         ELSE
           READ(45,1007,END = 1200) jc, parc, nceq  ! ecis06
         ENDIF      
C------Loop over the number of coupled equations
       DO nc = 1, nceq
C--------Reading the coupled level number nlev, the orbital momentum L,
C--------angular momentum j and Transmission coefficient Tlj,c(JC)
C--------(nlev=1 corresponds to the ground state)

         IF(unformatted) then 
           READ(45)   nlev, l, jj, dtmp
         ELSE
           READ(45,*) nlev, l, jj, dtmp
         ENDIF

C--------Selecting only ground state
         IF (dtmp.GT.1.D-15 .AND. l.LT.LLLMAX) THEN
C-----------Averaging over particle and target spin, summing over channel spin jc
            dtmp = dtmp*(2.D0*jc + 1.D0)/(2.D0*haf + 1.D0)
            Stl(l+1) = Stl(l+1) + dtmp/DBLE(2*l + 1)/DBLE(JO(1) + 1)
         ENDIF
       ENDDO
       goto 1100
 1200  CLOSE (45)
C
      ELSE
C
C      GNASH and others
C
       WRITE(22,998)EN,LMAX,(TRL(L),L=1,LMAX1)
       WRITE(24,996)EN,LLMA,TRLJ(1,2),
     *  (TRLJ(L,1),TRLJ(L+1,1),TRLJ(L,2),TRLJ(L+1,2),L=2,LLMA,2)
C
      ENDIF

 996  FORMAT (E12.6,I3/(6E11.4))
 998  FORMAT (E12.6,I3/(6F11.8))
 999  FORMAT(1X,I5,10X,4F20.10)
C     IF(MEPRI.NE.99) PRINT 111,LNJ1(IC),CR(IC,IC),CI(IC,IC)
C 111 FORMAT (3X,I10,2F20.10)

   17 CSR=CST-CSN(1)
      
   
      SINlcc = 0.d0
      IF(NMAX.LT.2) GO TO 955
      DO 55 N=2,NMAX
   55 SINlcc = SINlcc + CSN(N)

C  55 CSR=CSR-CSN(N)
  955 CONTINUE
      SF0=SF0/P0
      SF1=SF1/P1
      SF2=SF2/P2
      
      IF(MEJOB.EQ.2) GO TO 1117
      
      IF(MEPRI.EQ.98) THEN  ! EMPIRE=.TRUE.
C
C       for EMPIRE
C
        OPEN (46,FILE = trim(fname)//'_INC.LST')
        WRITE (46,'(A5,I6,1X,D12.6)') 'LMAX:', ltlmax, EN
        DO l = 0, ltlmax
          WRITE (46,*) l, stl(l+1)
        ENDDO
        WRITE (46,'(1x,A27,2x,6(D12.6,1x))') 
     &     'EL,TOT,REAC,INEL,CC,CSFus:',
     &   1000.d0*CSN(1), 1000.d0*CST, 1000.d0*CSR, 0.d0, 
     &   1000.d0*SINLcc, 1000.d0*(CSR-SINLcc)
C       WRITE (46,'(1x,I6)') 123456 
C       DO l = 1, ltlmax
C         WRITE (46,*) l-1, SNGL(sel(l))
C       ENDDO
        CLOSE (46)

        IF(unformatted) THEN
            OPEN (45,FILE = trim(fname)//'.INC',FORM = 'UNFORMATTED')
            IF (MEREL.EQ.0) then 
              itmp = 0 
              WRITE (45) ltlmax, EN, itmp
            ELSE
              itmp = 1 
              WRITE (45) ltlmax, EN, itmp
            ENDIF
            DO l = 0, ltlmax
              WRITE (45) stl(l+1)
            ENDDO
            WRITE (45)
     &       1000.d0*CSN(1), 1000.d0*CST, 1000.d0*CSR, 0.d0, 
     &       1000.d0*SINLcc, 1000.d0*(CSR-SINLcc)          
        ELSE
            OPEN (45,FILE = trim(fname)//'.INC')
            IF (MEREL.EQ.0) then 
              itmp = 0 
              WRITE (45,'(1x,I4,2x,e13.6,2x,i4)') ltlmax, EN, itmp
            ELSE
              itmp = 1 
              WRITE (45,'(1x,I4,2x,d12.6,2x,i4)') ltlmax, EN, itmp
            ENDIF
            DO l = 0, ltlmax
              WRITE (45,*) stl(l+1)
            ENDDO
            WRITE (45,'(1x,6(d12.6,1x))') 
     &       1000.d0*CSN(1), 1000.d0*CST, 1000.d0*CSR, 0.d0, 
     &       1000.d0*SINLcc, 1000.d0*(CSR-SINLcc)                      
        ENDIF


C
C       A new flag is introduced to signal storage of the Shape elastic XS (Sel(L))
C
C       l = 123456
C       WRITE (45) l 
C       DO l = 1, ltlmax
C         WRITE (45) sel(l)
C       ENDDO
        CLOSE (45)

      
C***********************************************************************
C  Kinematics:   lab  ===>  CM                                         *
C    With relativistic kinematics, the reduced mass is replaced by     *
C    the reduced total energy                                          *
C----------------------------------------------------------------------*
C  EN     = current lab kinetic energy                                 *
C  ecms   = current  CM kinetic energy                                 *
C  AI     = incident particle rest mass (in a.m.u.)                    *
C  AT     = target   nucleus  rest mass (in a.m.u.)                    *
C  AK2    = CM wave number                                             *
C----------------------------------------------------------------------*
C  AMUmev = a.m.u. in MeV                                              *
C----------------------------------------------------------------------*
C-------CONSTANTS COMPUTED FROM THE FUNDAMENTAL CONSTANTS, ATOMIC MASS, HBAR*C
C-------AND ALPHA, AS GIVEN IN THE EUROPEAN PHYSICAL JOURNAL, PAGE 73, VOLUME
C-------15 (2000) REFERRING FOR THESE VALUES TO THE 1998 CODATA SET WHICH MAY
C-------BE FOUND AT http://physics.nist.gov/constants
C-------CM=931.494013 +/- 0.000037 MeV
C-------The above value is the one used also in the ENDF-6 manual (April 2001, 2009)
C       AMUmev = 9.31494013D+02
C       CARBUN=931.49378D0  , this is the OPTMAN mass unit
        AMUmev = 9.31494013D+02
C-------CHB=197.3269601 +/- 0.0000078 (*1.0E-9 eV*cm)
        HHBarc = 197.3269601D0

        AI = ANEU
        ck2= (2.d0*AMUmev)/(HHBarc**2)
C
C-------From lab to CM (the input quantity is EN)
C
C
        IF (MEREL.EQ.0) then 
C---------Classical    kinematics
          ecms = EN*AT/(AI + AT)
          Ak2 = ck2*AI*AT/(AI + AT)*ecms
        ELSE
C
C---------Relativistic kinematics
          ecms = AMUmev*(AI + AT)*
     &           (SQRT(1.d0 + 2.d0*EN/(AMUmev*AT*((1.d0+AI/AT)**2)))
     &           - 1.d0)
          p2 = (EN*(EN + 2.d0*AMUmev*AI))
     &           /((1.d0 + AI/AT)**2 + 2.d0*EN/(AMUmev*AT))
          Ak2 = p2/(HHBarc*HHBarc)
        ENDIF
C       10 factor converts to mbarns
        dtmp = 10.D0*4.d0*atan(1.d0)/Ak2   
        sabs  =0.d0
        DO l = 0, ltlmax
          sabs   = sabs   + Stl(l + 1)*DBLE(2*l + 1)
        ENDDO
C       write(*,*) 'Test       :',dtmp*sabs,1000*(CSR-SINLcc)

C       Renormalizing TLs
        if(dtmp.gt.0.d0  .and. sabs.gt.0.d0) then
          DO l = 0, ltlmax
            Stl(l + 1)=Stl(l + 1)/(dtmp*sabs)*1000.d0*(CSR-SINLcc)
          ENDDO
          sabs  =0.d0
          DO l = 0, ltlmax
            sabs   = sabs   + Stl(l + 1)*DBLE(2*l + 1)
          ENDDO
C         write(*,*) 'Test renorm:',dtmp*sabs,1000*(CSR-SINLcc)
        endif

      ENDIF
      
 1117 CONTINUE

      CSR=CSR-SINlcc

      RETURN
      END
C     *******************************************************
      SUBROUTINE CHLOG
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM7.FOR'
      N=0
      K1=IABS(JSO-JSS)+2
      K2=JSO+JSS+2
      DO 5 L=K1,K2,2
      N=N+1
      IF(NPIO*(-1)**((L-3)/2)*NPIS)2,1,1
    1 LNO(N)=(L-3)/2
      NS1(N)=1
      GO TO 5
    2 LNO(N)=(L-1)/2
      NS1(N)=2
    5 CONTINUE
    6 NN1=N
      RETURN
      END
C     *******************************************************
      SUBROUTINE RACAH
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM1.FOR'
      INCLUDE 'PRIVCOM4.FOR'
      
      W=0.
      NB1=JA+JB+JC
      NB2=JC+JD+JE
      NB3=JA+JE+JF
      NB4=JB+JD+JF
      NM1=JA+JB+JD+JE
      NM2=JA+JC+JD+JF
      NM3=JB+JC+JE+JF
      NM=NM1
      NB=NB1
      IF(NM2-NM)1,2,2
    1 NM=NM2
    2 IF(NM3-NM)3,4,4
    3 NM=NM3
    4 IF(NB2-NB)5,5,6
    6 NB=NB2
    5 IF(NB3-NB)7,7,8
    8 NB=NB3
    7 IF(NB4-NB)9,9,10
   10 NB=NB4
    9 IF(NB-NM)12,12,11
   12 NF=JA+JB-JC
      FLN=A(NF+2)
      CL=FLN
      NF=JA-JB+JC
      FLN=A(NF+2)
      CL=CL+FLN
      NF=JB+JC-JA
      FLN=A(NF+2)
      CL=CL+FLN
      NF=NB1+2
      FLN=A(NF+2)
      CL=CL-FLN
      NF=JC+JD-JE
      FLN=A(NF+2)
      CL=CL+FLN
      NF=JC-JD+JE
      FLN=A(NF+2)
      CL=CL+FLN
      NF=JD+JE-JC
      FLN=A(NF+2)
      CL=CL+FLN
      NF=NB2+2
      FLN=A(NF+2)
      CL=CL-FLN
      NF=JA+JE-JF
      FLN=A(NF+2)
      CL=CL+FLN
      NF=JA-JE+JF
      FLN=A(NF+2)
      CL=CL+FLN
      NF=JE+JF-JA
      FLN=A(NF+2)
      CL=CL+FLN
      NF=NB3+2
      FLN=A(NF+2)
      CL=CL-FLN
      NF=JB+JD-JF
      FLN=A(NF+2)
      CL=CL+FLN
      NF=JB-JD+JF
      FLN=A(NF+2)
      CL=CL+FLN
      NF=JD+JF-JB
      FLN=A(NF+2)
      CL=CL+FLN
      NF=NB4+2
      FLN=A(NF+2)
      CL=0.5*(CL-FLN)
      NB=NB+2
      NM=NM+2
      DO 13 K=NB,NM,2
      C1=1.
      NF=K
      FLN=A(NF+2)
      CL1=CL+FLN
      I=K-2
      NF=I-NB1
      FLN=A(NF+2)
      CL1=CL1-FLN
      NF=I-NB2
      FLN=A(NF+2)
      CL1=CL1-FLN
      NF=I-NB3
      FLN=A(NF+2)
      CL1=CL1-FLN
      NF=I-NB4
      FLN=A(NF+2)
      CL1=CL1-FLN
      NF=NM1-I
      FLN=A(NF+2)
      CL1=CL1-FLN
      NF=NM2-I
      FLN=A(NF+2)
      CL1=CL1-FLN
      NF=NM3-I
      FLN=A(NF+2)
      CL1=CL1-FLN
      IF((NM1+I)/4*4.NE.NM1+I) C1=-1.
   13 W=W+C1*EXP(CL1)
   11 RETURN
      END
C     *******************************************************
      SUBROUTINE SOSIT
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
C     BELOW CARD IS NECESSARY IF MEMORY IS LESS THAN 32Mb
      DOUBLE PRECISION VL,VSL,PV,PW,WSL,PVC,PRC
      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM2.FOR'
      INCLUDE 'PRIVCOM5.FOR'
      INCLUDE 'PRIVCOM10.FOR'
      INCLUDE 'PRIVCOM13.FOR' 
      INCLUDE 'PRIVCOM14.FOR'     

C     WRITE(21,'(1X,11I6,1X,6(f8.4,1X))') NCLL

      SPI2=ASP*(ASP+1.D0)
      NCL2=NCLL*NCLL
      IF(MEPOT.EQ.1) GO TO 27
      LAS2=LAS+1
      LAS1=LAS
      GO TO 26
   27 LAS2=(LAS+2)/2
      LAS1=LAS2-1
   26 NCLA=NCLL*LAS1
      LAS8=9
      LASC=1
      IF(LAS2.GE.3) LAS8=18
      IF(LAS2.GE.3) LASC=2
      ICLL=NCLL*LAS8
      MNULA=300*NUR*LAS2
      MLA=300*LAS2
      MNU=300*NUR
      NH=NH1
      ST=1.
      M=8
      M4=4*M
      MF=M4-1
      STEP1=STEP/M

C     K-NUMBER OF INDEPENDENT SOLUTION

      DO 1 K=1,NCLL
      NCLK=(K-1)*NCLL
      KC1=ICLL*(K-1)
      K1=(K-1)*LAS1-1
      ST2=STEP1*STEP1/240.D0
      LL1=LNJ1(K)
      NPI1=KPJ1(K)
      NUK=NNJ1(K)
      MNULAK=MNULA*(NUK-1)
      MNUK=MNU*(NUK-1)
      C=ST*(STEP1)**(LL1+1)
      C2=ST*(2.*STEP1)**(LL1+1)
      C3=ST*(3.*STEP1)**(LL1+1)

C     L-LINE OF INDEPENDENT SOLUTION

      DO 2 L=1,NCLL
      KL=NCLK+L
      KC2=KC1+LAS8*(L-1)
      K2=K1+(L-1)*NCLA
      LL=LNJ1(L)
      NUL=NNJ1(L)
      NU=NUL
      MNULKL=MNULAK+MLA*(NUL-1)
      MNUKL=MNUK+300*(NUL-1)
      NPI2=KPJ1(L)
      JJ=JNJ1(L)
      CJ1=JJ/2.
      SOP=(CJ1-LL)*(CJ1+LL+1)-SPI2
      X4(KL)=0.
      Y4(KL)=0.
      X5(KL)=0.
      Y5(KL)=0.
      VR1(KL)=0.
      VI1(KL)=0.
      UR=0.
      UI=0.

      
      
      C=CVNRPN(KL)
C     IF(NU.NE.NUN.AND.JO(NU).EQ.JO(NUN).AND. JO(NU).EQ.0) C=1.D0
      IF(NCA(NU).NE.NCA(NUK).AND.JO(NU).EQ.JO(NUK))UR=UR+PV(MNULKL+1)*C
      IF(NCA(NU).NE.NCA(NUK).AND.JO(NU).EQ.JO(NUK))UI=UI+PW(MNULKL+1)*C
      
      
      
      
C     IF(NU.NE.NUK.AND.JO(NU).EQ. JO(NUK)) UR=UR+PVV(MNUKL+1)*CVNRPN(KL)
C     IF(NU.NE.NUK.AND.JO(NU).EQ. JO(NUK)) UI=UI+PWW(MNUKL+1)*CVNRPN(KL)

      IF(MEVOL.GE.1.AND.MEPOT.EQ.1) UR=UR+PV(MNULKL+2)*CVNRV0(KL)
      IF(MEVOL.GE.1.AND.MEPOT.EQ.1) UI=UI+PW(MNULKL+2)*CVNRV0(KL)

         ACS=0.D0

      IF(LAS2.LT.2) GO TO 300
      DO 30 LLL=2,LAS2
      LP1=MNULKL+LLL
      PV1=PV(LP1)
      PW1=PW(LP1)
         IF(LLL.EQ.2) PVS=PV1
         IF(LLL.EQ.2) PWS=PW1
      LL2=K2+LLL
      A=CVNR(LL2)
      IF(NPI1.NE.NPI2) GO TO 36
      UR=UR+A*PV1
      UI=UI+A*PW1
      GO TO 41
   36 UR=UR-A*PW1
      UI=UI+A*PV1
   41 IF(MEPOT.EQ.1) GO TO 30
      IF(LLL.GT.3) GO TO 30
      LAI=3
      LAF=5
      IF(LLL.EQ.3) LAI=1
      IF(LLL.EQ.3) LAF=9
      DO 43 LA=LAI,LAF
      LALAS=LASC*(LA-1)+LLL-1
      LP1=LALAS
      PV1=PVC(LP1)
      LL2=KC2+LALAS
      A=CVNC(LL2)
          IF(LLL.EQ.3.AND.LA.EQ.1) ACS=-A*CBET0
          IF(LA.EQ.1) GO TO 43
      IF(NPI1.NE.NPI2) GO TO 42
      UR=UR+A*PV1
      GO TO 43
   42 UI=UI+A*PV1
43    CONTINUE
30    CONTINUE
      IF(MEVOL.EQ.0.OR.MEPOT.EQ.1) GO TO 300
          UR=UR+ACS*PVS
          UI=UI+ACS*PWS
300   POTR2=PV(MNULKL+1)+VSL(MNUKL+1)*SOP+VL(1)*LL*(LL+1)-WN(NU)
      IF(L-K) 3,4,3
    4 X4(KL)=ST*(3.*STEP1)**(LL+1)
      X5(KL)=ST*(4.*STEP1)**(LL+1)
      UR=UR+POTR2
      UI=UI+PW(MNULKL+1)+WSL(MNUKL+1)*SOP
      IF(LL.EQ.1) VR1(KL)=2.*ST
    3 VR2(KL)=UR*C
      VI2(KL)=UI*C
      UR=0.
      UI=0.
      



      C=CVNRPN(KL)
C     IF(NU.NE.NUN.AND.JO(NU).EQ.JO(NUN).AND. JO(NU).EQ.0) C=1.D0
      IF(NCA(NU).NE.NCA(NUK).AND.JO(NU).EQ.JO(NUK))UR=UR+PV(MNULKL+2)*C
      IF(NCA(NU).NE.NCA(NUK).AND.JO(NU).EQ.JO(NUK))UI=UI+PW(MNULKL+2)*C
      

C     IF(NU.NE.NUK.AND.JO(NU).EQ. JO(NUK)) UR=UR+PVV(MNUKL+2)*CVNRPN(KL)
C     IF(NU.NE.NUK.AND.JO(NU).EQ. JO(NUK)) UI=UI+PWW(MNUKL+2)*CVNRPN(KL)

      IF(MEVOL.GE.1.AND.MEPOT.EQ.1) UR=UR+PV(MNULKL+3)*CVNRV0(KL)
      IF(MEVOL.GE.1.AND.MEPOT.EQ.1) UI=UI+PW(MNULKL+3)*CVNRV0(KL)      
      
          ACS=0.D0
      IF(LAS2.LT.2) GO TO 310
      DO 31 LLL=2,LAS2
      LL2=LAS2+LLL+MNULKL
      PV1=PV(LL2)
      PW1=PW(LL2)
         IF(LLL.EQ.2) PVS=PV1
         IF(LLL.EQ.2) PWS=PW1
      LL2=K2+LLL
      A=CVNR(LL2)
      IF(NPI1.NE.NPI2) GO TO 37
      UR=UR+A*PV1
      UI=UI+A*PW1
      GO TO 44
   37 UR=UR-A*PW1
      UI=UI+A*PV1
   44 IF(MEPOT.EQ.1) GO TO 31
      IF(LLL.GT.3) GO TO 31
      LAI=3
      LAF=5
      IF(LLL.EQ.3) LAI=1
      IF(LLL.EQ.3) LAF=9
      DO 45 LA=LAI,LAF
      LALAS=LASC*(LA-1)+LLL-1
      LP1=LALAS+LAS8
      PV1=PVC(LP1)
      LL2=KC2+LALAS
      A=CVNC(LL2)
          IF(LLL.EQ.3.AND.LA.EQ.1) ACS=-A*CBET0
          IF(LA.EQ.1) GO TO 45
      IF(NPI1.NE.NPI2) GO TO 46
      UR=UR+A*PV1
      GO TO 45
   46 UI=UI+A*PV1
   45 CONTINUE
   31 CONTINUE
      IF(MEVOL.EQ.0.OR.MEPOT.EQ.1) GO TO 310
          UR=UR+ACS*PVS
          UI=UI+ACS*PWS
310   LL2=LAS2+1+MNULKL
      POTR2=PV(LL2)+VSL(MNUKL+2)*SOP+VL(2)*LL*(LL+1)-WN(NU)
      IF(L.EQ.K) UR=UR+POTR2
      IF(L.EQ.K) UI=UI+PW(LL2)+WSL(MNUKL+2)*SOP
      VR3(KL)=UR*C2
      VI3(KL)=UI*C2
      UR=0.
      UI=0.



      C=CVNRPN(KL)
C     IF(NU.NE.NUN.AND.JO(NU).EQ.JO(NUN).AND. JO(NU).EQ.0) C=1.D0
      IF(NCA(NU).NE.NCA(NUK).AND.JO(NU).EQ.JO(NUK))UR=UR+PV(MNULKL+3)*C
      IF(NCA(NU).NE.NCA(NUK).AND.JO(NU).EQ.JO(NUK))UI=UI+PW(MNULKL+3)*C
      

C      IF(NU.NE.NUK.AND.JO(NU).EQ. JO(NUK)) UR=UR+PVV(MNUKL+3)*CVNRPN(KL)
C     IF(NU.NE.NUK.AND.JO(NU).EQ. JO(NUK)) UI=UI+PWW(MNUKL+3)*CVNRPN(KL)

C      WRITE(21,'(1X,3I6,1X,6(f8.4,1X))') k,l,ll,PV(LL2),VSL(MNUKL+2),
C     *SOP,VL(2)

C     IF(NU.NE.NUK.AND.JO(NU).EQ. JO(NUK))
C     *     WRITE(21,'(1X,11I6,1X,9(f8.4,1X))') NCLL,NU,NUK,KL, MNUKL,
C     *JNJ1(L), JNJ1(K),LNJ1(L),LNJ1(K), JO(NU),
C     *jO(NUK),PVV(MNUKL+3),CVNRPN(KL)


      IF(MEVOL.GE.1.AND.MEPOT.EQ.1) UR=UR+PV(MNULKL+4)*CVNRV0(KL)
      IF(MEVOL.GE.1.AND.MEPOT.EQ.1) UI=UI+PW(MNULKL+4)*CVNRV0(KL)

           ACS=0.D0

      IF(LAS2.LT.2) GO TO 320
      DO 32 LLL=2,LAS2
      LL2=LAS2*2+LLL+MNULKL
      PV1=PV(LL2)
      PW1=PW(LL2)
         IF(LLL.EQ.2) PVS=PV1
         IF(LLL.EQ.2) PWS=PW1
      LL2=K2+LLL
      A=CVNR(LL2)
      IF(NPI1.NE.NPI2) GO TO 38
      UR=UR+A*PV1
      UI=UI+A*PW1
      GO TO 47
   38 UR=UR-A*PW1
      UI=UI+A*PV1
   47 IF(MEPOT.EQ.1) GO TO 32
      IF(LLL.GT.3) GO TO 32
      LAI=3
      LAF=5
      IF(LLL.EQ.3) LAI=1
      IF(LLL.EQ.3) LAF=9
      DO 48 LA=LAI,LAF
      LALAS=LASC*(LA-1)+LLL-1
      LP1=LALAS+LAS8*2
      PV1=PVC(LP1)
      LL2=KC2+LALAS
      A=CVNC(LL2)
          IF(LLL.EQ.3.AND.LA.EQ.1) ACS=-A*CBET0
          IF(LA.EQ.1) GO TO 48
      IF(NPI1.NE.NPI2) GO TO 49
      UR=UR+A*PV1
      GO TO 48
   49 UI=UI+A*PV1
   48 CONTINUE
   32 CONTINUE
      IF(MEVOL.EQ.0.OR.MEPOT.EQ.1) GO TO 320
          UR=UR+ACS*PVS
          UI=UI+ACS*PWS
320   LL2=LAS2*2+1+MNULKL
      POTR2=PV(LL2)+VSL(MNUKL+3)*SOP+VL(3)*LL*(LL+1)-WN(NU)
      IF(L.EQ.K) UR=UR+POTR2
      IF(L.EQ.K) UI=UI+PW(LL2)+WSL(MNUKL+3)*SOP
      VR4(KL)=UR*C3
      VI4(KL)=UI*C3
      P(KL)=VR1(KL)
      R(KL)=VI1(KL)
    2 CONTINUE
    1 CONTINUE
      DO 5 MM=4,MF
      MM1=(MM-1)*LAS2
      MM2=MM1+1
      MMC=LAS8*(MM-1)
C     L-LINE OF INDEPENDENT SOLUTION

      DO 15 L=1,NCLL
      LPN1=(L-1)*NCLL
      KC1=ICLL*(L-1)
      K1=(L-1)*NCLA-1
      DO 6 K=1,NCLL
      VR5E(K)=0.
    6 VI5E(K)=0.
      LL=LNJ1(L)
      JJ=JNJ1(L)
      NU=NNJ1(L)
      MNULAL=MNULA*(NU-1)
      MNUL=MNU*(NU-1)
      NPI1=KPJ1(L)
      CJ1=JJ/2.
      SOP=(CJ1-LL)*(CJ1+LL+1)-SPI2

C     N-NUMBER OF LINE OF SOLUTION COUPLED WITH L-LINE

      DO 7 N=1,NCLL
      LN=LPN1+N
      KC2=KC1+LAS8*(N-1)
      K2=K1+(N-1)*LAS1
      NPI2=KPJ1(N)
      NUN=NNJ1(N)
      MNULLN=MNULAL+MLA*(NUN-1)
      MNULN=MNUL+300*(NUN-1)
      POTR5=PV(MNULLN+MM2)+VSL(MNULN+MM)*SOP+VL(MM)*LL*(LL+1)-WN(NU)
      UR=0.
      UI=0.
      C=CVNRPN(LN)
C     IF(NU.NE.NUN.AND.JO(NU).EQ.JO(NUN).AND. JO(NU).EQ.0) C=1.D0
      IF(NCA(NU).NE.NCA(NUN).AND.JO(NU).EQ.JO(NUN))
     *UR=UR+PV(MNULLN+MM2)*C
      IF(NCA(NU).NE.NCA(NUN).AND.JO(NU).EQ.JO(NUN))
     *UI=UI+PW(MNULLN+MM2)*C

C     write (21,9797) mm,ur, ui,2.2, l,n
C 9797 format(i5, 3e12.5,2i5)

C      WRITE(21,'(1X,4I6,1X,6(f8.4,1X))')mm,l,n

C     WRITE(21,'(1X,4I6,1X,6(f8.4,1X))')mm,l,n,ll,
C     *PV(MNULLN+MM2),pw(MNULLN+MM2),VSL(MNUKL+mm),SOP,VL(mm)

C     IF(NU.NE.NUn.AND.JO(NU).EQ. JO(NUn))
C     *     WRITE(21,'(1X,12I5,1X,6(f8.4,1X))') NCLL,NU,NUn,ln, MNUln,
C     *JNJ1(L), JNJ1(n),LNJ1(L),LNJ1(n), JO(NU),
C     *jO(NUn),mm, PVV(MNUln+mm),Pww(MNUln+mm),CVNRPN(ln),WN(NU)

      IF(MEVOL.GE.1.AND.MEPOT.EQ.1) UR=UR+PV(MNULLN+MM2+1)*CVNRV0(LN)
      IF(MEVOL.GE.1.AND.MEPOT.EQ.1) UI=UI+PW(MNULLN+MM2+1)*CVNRV0(LN)
          ACS=0.D0

      IF(LAS2.LT.2) GO TO 330
      DO 33 LLL=2,LAS2
      LL2=MM1+LLL+MNULLN
      PV1=PV(LL2)
      PW1=PW(LL2)
         IF(LLL.EQ.2) PVS=PV1
         IF(LLL.EQ.2) PWS=PW1
      LL2=K2+LLL
      A=CVNR(LL2)
      IF(NPI1.NE.NPI2) GO TO 39
      UR=UR+A*PV1
      UI=UI+A*PW1
C     DIAGONAL COUPLING FOR (P,N)




      GO TO 50
   39 UR=UR-A*PW1
      UI=UI+A*PV1
   50 IF(MEPOT.EQ.1) GO TO 33
      IF(LLL.GT.3) GO TO 33
      LAI=3
      LAF=5
      IF(LLL.EQ.3) LAI=1
      IF(LLL.EQ.3) LAF=9
      DO 51 LA=LAI,LAF
      LALAS=LASC*(LA-1)+LLL-1
      LP1=LALAS+MMC
      PV1=PVC(LP1)
      LL2=KC2+LALAS
      A=CVNC(LL2)
          IF(LLL.EQ.3.AND.LA.EQ.1) ACS=-A*CBET0
          IF(LA.EQ.1) GO TO 51
      IF(NPI1.NE.NPI2) GO TO 52
      UR=UR+A*PV1
      GO TO 51
   52 UI=UI+A*PV1
   51 CONTINUE
   33 CONTINUE
      IF(MEVOL.EQ.0.OR.MEPOT.EQ.1) GO TO 330
          UR=UR+ACS*PVS
          UI=UI+ACS*PWS
330   IF(L.EQ.N) UR=UR+POTR5
      IF(L.EQ.N) UI=UI+PW(MM2+MNULLN)+WSL(MNULN+MM)*SOP


C      write (21,9797) mm,ur, ui, potr5,l,n

      DO 7 K=1,NCLL
      KN=(K-1)*NCLL+N
      X5N=X5(KN)
      Y5N=Y5(KN)
      VR5E(K)=VR5E(K)+UR*X5N-UI*Y5N
    7 VI5E(K)=VI5E(K)+UR*Y5N+UI*X5N

C     K-NUMBER OF INDEPENDENT SOLUTION

      DO 19 K=1,NCLL
      NCLK=(K-1)*NCLL
      KL=NCLK+L
      VR5L=VR5E(K)
      VI5L=VI5E(K)
      WR2=VR2(KL)
      WR3=VR3(KL)
      WR4=VR4(KL)
      WI2=VI2(KL)
      WI3=VI3(KL)
      WI4=VI4(KL)
      XW5=X5(KL)
      YW5=Y5(KL)
      X6L=2.*XW5-X4(KL)+ST2*(299.*VR5L-176.*WR4+194.*WR3-
     *96.*WR2+19.*VR1(KL))
      Y6L=2.*YW5-Y4(KL)+ST2*(299.*VI5L-176.*WI4+194.*WI3-
     *96.*WI2+19.*VI1(KL))
      VR1(KL)=WR2
      VR2(KL)=WR3
      VR3(KL)=WR4
      VR4(KL)=VR5L
      VI1(KL)=WI2
      VI2(KL)=WI3
      VI3(KL)=WI4
      VI4(KL)=VI5L
      X4(KL)=XW5
      Y4(KL)=YW5
      X6(KL)=X6L
      Y6(KL)=Y6L
   19 CONTINUE
   15 CONTINUE
      IF(MM/M*M.NE.MM) GO TO 12
      IF(MM-2*M) 8,9,10
    8 DO 901 NN=1,NCL2
      FREM(1,NN)=0.
       FIEM(1,NN)=0.
      FREM(2,NN)=X6L
       FIEM(2,NN)=Y6L
      W(NN)=VR4(NN)
  901 V(NN)=VI4(NN)
      GO TO 302
    9 DO 910 NN=1,NCL2
      FREM(3,NN)=X6L
       FIEM(3,NN)=Y6L
      S(NN)=VR4(NN)
  910 B(NN)=VI4(NN)
      GO TO 302
   10 IF(MM-3*M) 11,11,302
   11 DO 902 NN=1,NCL2
      FREM(4,NN)=X6L
       FIEM(4,NN)=Y6L
      G(NN)=VR4(NN)
      Z(NN)=VI4(NN)
      E(NN)=X4(NN)
  902 F(NN)=Y4(NN)
      GO TO 302
   12 IF(MM.NE.MF) GO TO 302
      DO 911 NN=1,NCL2
       FREM(5,NN)=X6L
       FIEM(5,NN)=Y6L
      Q(NN)=X6(NN)
  911 T(NN)=Y6(NN)
  302 DO 24 L=1,NCL2
      X5(L)=X6(L)
   24 Y5(L)=Y6(L)
    5 CONTINUE
      DO 800 NB=1,NCL2
      VR1(NB)=P(NB)
      VI1(NB)=R(NB)
      VR2(NB)=W(NB)
      VI2(NB)=V(NB)
      VR3(NB)=S(NB)
      VI3(NB)=B(NB)
      VR4(NB)=G(NB)
      VI4(NB)=Z(NB)
      X4(NB)=E(NB)
      Y4(NB)=F(NB)
      X5(NB)=Q(NB)
      Y5(NB)=T(NB)
 800  CONTINUE
      ST2=STEP*STEP/240.
       NHTT=5
      DO 13 MM=M4,NH
       NHTT=NHTT+1
      MM1=(MM-1)*LAS2
      MMC=LAS8*(MM-1)
      MM2=MM1+1
      DO 17 L=1,NCLL
      LPN1=(L-1)*NCLL
      KC1=ICLL*(L-1)
      K1=(L-1)*NCLA-1
      DO 16 K=1,NCLL
      VR5E(K)=0.
   16 VI5E(K)=0.
      LL=LNJ1(L)
      JJ=JNJ1(L)
      NU=NNJ1(L)
      MNULAL=MNULA*(NU-1)
      MNUL=MNU*(NU-1)
      NPI1=KPJ1(L)
      CJ1=JJ/2.
      SOP=(CJ1-LL)*(CJ1+LL+1)-SPI2
      DO 14 N=1,NCLL
      LN=LPN1+N
      NPI2=KPJ1(N)
      NUN=NNJ1(N)
      MNULLN=MNULAL+MLA*(NUN-1)
      MNULN=MNUL+300*(NUN-1)
      POTR5=PV(MNULLN+MM2)+VSL(MNULN+MM)*SOP+VL(MM)*LL*(LL+1)-WN(NU)
      KC2=KC1+LAS8*(N-1)
      K2=K1+(N-1)*LAS1
      UR=0.
      UI=0.

      C=CVNRPN(LN)
c      IF(NU.NE.NUN.AND.JO(NU).EQ.JO(NUN).AND. JO(NU).EQ.0) C=1.D0
      IF(NCA(NU).NE.NCA(NUN).AND.JO(NU).EQ.JO(NUN))
     *UR=UR+PV(MNULLN+MM2)*C
      IF(NCA(NU).NE.NCA(NUN).AND.JO(NU).EQ.JO(NUN))
     *UI=UI+PW(MNULLN+MM2)*C

C           write (21,9797) mm,ur, ui,2.2, l,n
C 
C      WRITE(21,'(1X,4I6,1X,6(f8.4,1X))')mm,l,n

C     WRITE(21,'(1X,4I6,1X,6(f8.4,1X))')mm,l,n,ll,
C     *PV(MNULLN+MM2),pw(MNULLN+MM2),VSL(MNUKL+mm),SOP,VL(mm)

C     IF(NU.NE.NUn.AND.JO(NU).EQ. JO(NUn))
C     *     WRITE(21,'(1X,12I5,1X,6(f8.4,1X))') NCLL,NU,NUn,ln, MNUln,
C     *JNJ1(L), JNJ1(n),LNJ1(L),LNJ1(n), JO(NU),
C     *jO(NUn),mm, PVV(MNUln+mm),Pww(MNUln+mm),CVNRPN(ln)

      IF(MEVOL.GE.1.AND.MEPOT.EQ.1) UR=UR+PV(MNULLN+MM2+1)*CVNRV0(LN)
      IF(MEVOL.GE.1.AND.MEPOT.EQ.1) UI=UI+PW(MNULLN+MM2+1)*CVNRV0(LN)
         ACS=0.D0

      IF(LAS2.LT.2) GO TO 340
      DO 34 LLL=2,LAS2
      LL2=MM1+LLL+MNULLN
      PV1=PV(LL2)
      PW1=PW(LL2)
         IF(LLL.EQ.2) PVS=PV1
         IF(LLL.EQ.2) PWS=PW1
      LL2=K2+LLL
      A=CVNR(LL2)
      IF(NPI1.NE.NPI2) GO TO 40
      UR=UR+A*PV1
      UI=UI+A*PW1

C     DIAGONAL COUPLING FOR (P,N)

      


      GO TO 53
   40 UR=UR-A*PW1
      UI=UI+A*PV1
   53 IF(MEPOT.EQ.1) GO TO 34
      IF(LLL.GT.3) GO TO 34
      LAI=3
      LAF=5
      IF(LLL.EQ.3) LAI=1
      IF(LLL.EQ.3) LAF=9
      DO 54 LA=LAI,LAF
      LALAS=LASC*(LA-1)+LLL-1
      LP1=LALAS+MMC
      PV1=PVC(LP1)
      LL2=KC2+LALAS
      A=CVNC(LL2)
          IF(LLL.EQ.3.AND.LA.EQ.1) ACS=-A*CBET0
          IF(LA.EQ.1) GO TO 54
      IF(NPI1.NE.NPI2) GO TO 55
      UR=UR+A*PV1
      GO TO 54
   55 UI=UI+A*PV1
   54 CONTINUE
   34 CONTINUE
      IF(MEVOL.EQ.0.OR.MEPOT.EQ.1) GO TO 340
          UR=UR+ACS*PVS
          UI=UI+ACS*PWS
340   IF(L.EQ.N) UR=UR+POTR5
      IF(L.EQ.N) UI=UI+PW(MM2+MNULLN)+WSL(MNULN+MM)*SOP
      DO 14 K=1,NCLL
      KN=(K-1)*NCLL+N
      X5N=X5(KN)
      Y5N=Y5(KN)
      VR5E(K)=VR5E(K)+UR*X5N-UI*Y5N
   14 VI5E(K)=VI5E(K)+UR*Y5N+UI*X5N
      DO 18 K=1,NCLL
      NCLK=(K-1)*NCLL
      KL=NCLK+L
      VR5L=VR5E(K)
      VI5L=VI5E(K)
      WR2=VR2(KL)
      WR3=VR3(KL)
      WR4=VR4(KL)
      WI2=VI2(KL)
      WI3=VI3(KL)
      WI4=VI4(KL)
      XW5=X5(KL)
      YW5=Y5(KL)
      X6L=2.*XW5-X4(KL)+ST2*(299.*VR5L-176.*WR4+194.*WR3-
     *96.*WR2+19.*VR1(KL))
      Y6L=2.*YW5-Y4(KL)+ST2*(299.*VI5L-176.*WI4+194.*WI3-
     *96.*WI2+19.*VI1(KL))
      VR1(KL)=WR2
      VR2(KL)=WR3
      VR3(KL)=WR4
      VR4(KL)=VR5L
      VI1(KL)=WI2
      VI2(KL)=WI3
      VI3(KL)=WI4
      VI4(KL)=VI5L
      X4(KL)=XW5
      Y4(KL)=YW5
      X6(KL)=X6L
      Y6(KL)=Y6L
       FREM(NHTT,KL)=X6L
       FIEM(NHTT,KL)=Y6L
      IF(MM-NH+2) 18,21,22
   21 FRF1(K,L)=X6L
      FIF1(K,L)=Y6L
      GO TO 18
   22 IF(MM-NH) 18,23,18
   23 FRF2(K,L)=X6L
      FIF2(K,L)=Y6L
   18 CONTINUE
   17 CONTINUE
      DO 25 L=1,NCL2
      X5(L)=X6(L)
   25 Y5(L)=Y6(L)
      IF(NCLL.LE.20) GO TO 13
      IF(MM.NE.M4+2) GO TO 13
      DO 20 K=1,NCLL
      NCLK=(K-1)*NCLL
      KK=NCLK+K
      LL1=LNJ1(K)
      AA=(STEP*(MM+1)*2.178281/(2.*LL+1.))**(LL1+1)/X5(KK)
      DO 20 L=1,NCLL
      KL=NCLK+L
      X4(KL)=X4(KL)*AA
      X5(KL)=X5(KL)*AA
      Y4(KL)=Y4(KL)*AA
      Y5(KL)=Y5(KL)*AA
      VR1(KL)=VR1(KL)*AA
      VR2(KL)=VR2(KL)*AA
      VR3(KL)=VR3(KL)*AA
      VR4(KL)=VR4(KL)*AA
      VI1(KL)=VI1(KL)*AA
      VI2(KL)=VI2(KL)*AA
      VI3(KL)=VI3(KL)*AA
      VI4(KL)=VI4(KL)*AA
   20 CONTINUE
   13 CONTINUE
      RETURN
      END
C     *******************************************************
      SUBROUTINE MASCT
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      REAL*16 FF1,FF2,GG1,GG2,RF1,RF2,F1I,F2I,RG1,RG2,G1I,G2I,ZNI,
     *C1,C2,D1,D2 !,BMR,BMI,ABR,ABI


      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM10.FOR'
      INCLUDE 'PRIVCOM13.FOR'          
      CSS=0.
      DO 1 I=1,NCLL
      NR=NNJ1(I)
      L=LNJ1(I)+1
      IF(NR.GT.NMAX) GO TO 5
      FF1=FBR1(NR,L)
      FF2=FBR2(NR,L)
      GG1=FNR1(NR,L)
      GG2=FNR2(NR,L)
      ZN=FF2*GG1-FF1*GG2
      GO TO 6
    5 RF1=FBR1(NR,L)
      RF2=FBR2(NR,L)
      F2I=FBI2(NR,L)
      F1I=FBI1(NR,L)
      RG1=FNR1(NR,L)
      RG2=FNR2(NR,L)
      G1I=FNI1(NR,L)
      G2I=FNI2(NR,L)
      ZNI=RF2*G1I+F2I*RG1-RF1*G2I-F1I*RG2
    6 DO 2 K=1,NCLL
      K1=(K-1)*NCLL+I
      C1=FRF1(K,I)
      C2=FRF2(K,I)
      D1=FIF1(K,I)
      D2=FIF2(K,I)
      IF(NR.GT.NMAX) GO TO 7
      BMR(K1)=(C2*FF1-C1*FF2)/ZN
      BMI(K1)=(D2*FF1-D1*FF2)/ZN
      ABR(K1)=(C2*GG1-C1*GG2)/ZN-BMI(K1)
      ABI(K1)=(D2*GG1-D1*GG2)/ZN+BMR(K1)
      GO TO 2
    7 BMR(K1)=(D2*RF1+C2*F1I-D1*RF2-C1*F2I)/ZNI
      BMI(K1)=(D2*F1I-C2*RF1-D1*F2I+C1*RF2)/ZNI
      ABR(K1)=(D2*RG1+C2*G1I-D1*RG2-C1*G2I)/ZNI-BMI(K1)
      ABI(K1)=(D2*G1I-C2*RG1-D1*G2I+C1*RG2)/ZNI+BMR(K1)
    2 CONTINUE
    1 CONTINUE
      INFOR=1
#ifdef LAPACK      
      CALL INMATLU
#else
      CALL INMAT
#endif
      DO 3 I=1,NCLL
      I1=(I-1)*NCLL
      DO 3 K=1,NCLL
      NR=NNJ1(K)
      C=SQRT(WNK(NR)/WNK(1))
      C1=SQRT(ABS(WN(NR)/WN(1)))
      CCR=0.
      CCI=0.
      DO 4 L=1,NCLL
      IL=I1+L
      L1=(L-1)*NCLL
      LK=L1+K
      AABR=ABR(IL)
      BBMR=BMR(LK)
      AABI=ABI(IL)
      BBMI=BMI(LK)
      CCR=CCR-AABR*BBMR+AABI*BBMI
      CCI=CCI-AABR*BBMI-AABI*BBMR
    4 CONTINUE
      IF(I.GT.INCC.OR.K.GT.INCR) GO TO 13
      CSS=CSS+(CCR*CCR+CCI*CCI)*C1
C     IF(MEPRI.LT.98) PRINT 111,CCR,CCI
   13 CR(I,K)=CCR*C
      CI(I,K)=CCI*C
      IF(I.GT.10) GO TO 3
      CRD(NSS,I,K)=CCR*C
      CID(NSS,I,K)=CCI*C
    3 CONTINUE
      IF(MESOL.EQ.2) GO TO 14
      DO 9 I=1,NCLL
      I1=(I-1)*NCLL
      DO 9 K=1,NCLL
      IK=I1+K
      ABR(IK)=BMR(IK)
      ABI(IK)=BMI(IK)
    9 CONTINUE
      INFOR=2
#ifdef LAPACK      
      CALL INMATLU
#else
      CALL INMAT
#endif
      DO 10 I=1,NCLL
      I1=(I-1)*NCLL
      DO 10 K=1,NCLL
      IK=I1+K
      AUR(IK)=0.
      AUI(IK)=0.
      DO 10 J=1,NCLL
      IJ=I1+J
      J1=(J-1)*NCLL
      JK=J1+K
      AUR(IK)=AUR(IK)+CR(I,J)*ABR(JK)-CI(I,J)*ABI(JK)
      AUI(IK)=AUI(IK)+CI(I,J)*ABR(JK)+CR(I,J)*ABI(JK)
   10 CONTINUE
C      DO 11 I=1,NCLL
C      I1=(I-1)*NCLL
C      DO 11 K=1,NCLL
C      IK=I1+K
C      ABR(IK)=AUR(IK)
C      ABI(IK)=AUI(IK)
C   11 CONTINUE
C      INFOR=3
C      CALL INMAT
C      DO 12 I=1,NCLL
C      I1=(I-1)*NCLL
C      DO 12 K=1,NCLL
C      IK=I1+K
C      AUR(IK)=ABR(IK)
C      AUI(IK)=ABI(IK)
C   12 CONTINUE
C     IF(MEPRI.LT.98) PRINT 967,CR(1,2),CI(1,2)
C 967 FORMAT(1X,'M1',2D20.7)
C     WRITE (16) CR,CI,LNJ1,JNJ1,NNJ1
C     IF(MEPRI.LT.98) PRINT 9998,
C    *   (LNJ1(I8),JNJ1(I8),NNJ1(I8),I8=1,NCLL),NCLL
C9999 FORMAT(10X,'MASCT')
C7999 FORMAT(10X,'MASC1')
C9998 FORMAT(3I5)
   14 RETURN
      END
C     *******************************************************
      SUBROUTINE INMAT
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     IMPLICIT none
      REAL*16 PIVR,PIVI,ELEMR,ELEMI,DETR,DETI,ABRKK,ABIKK,ABRKJ,ABIKJ,
     *AMMRIK,AMMIIK,ABRKM,ABIKM 
      REAL*16  SUMDR,SUMDI,SUMNR,SUMNI!,SUMDR2,SUMDI2,SUMNR2,SUMNI2
      REAL*16 MOD2KK,AAR,AAI,MAR,MAI,EPS!,AAR2,AAI2,MAR2,MAI2
      
c      REAL*16 AMMRKJ,AMMIKJ
      
      REAL*16 AMRIK,AMIIK,ABRIK,ABIIK,AMRKM,AMIKM
      
      
      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM10.FOR'
    
               
      DETR=1.Q0
      DETI=1.Q0
      EPS=1.0Q-18
      ONEQ=1.Q0
C      IR=1
      LG=NCLL
      LG2=LG*LG
      DO 1 I=1,LG2
      AMMR(I)=ABR(I)
    1 AMMI(I)=ABI(I)
      DO 60 K=1,LG
      K1G=(K-1)*LG
      KK=K1G+K
      ABRKK=ABR(KK)
      ABIKK=ABI(KK)
      DETR=DETR*ABRKK-DETI*ABIKK
      DETI=DETR*ABIKK+DETI*ABRKK
      MOD2KK=ABRKK**2+ABIKK**2
      PIVR=ABRKK/MOD2KK
      PIVI=-ABIKK/MOD2KK
     
      DO 20 J=1,LG
      KJ=K1G+J
      ABRKJ=ABR(KJ)
      ABIKJ=ABI(KJ)
      ABR(KJ)=ABRKJ*PIVR-ABIKJ*PIVI
   20 ABI(KJ)=ABIKJ*PIVR+ABRKJ*PIVI
      ABR(KK)=PIVR
      ABI(KK)=PIVI
      
      DO 50 I=1,LG
      I1G=(I-1)*LG
      IK=I1G+K
      IF(I-K)30,50,30
   30 ELEMR=ABR(IK)
      ELEMI=ABI(IK)
      ABR(IK)=0.Q0
      ABI(IK)=0.Q0
      
      DO 40 J=1,LG
      IJ=I1G+J
      KJ=K1G+J
      ABRKJ=ABR(KJ)
      ABIKJ=ABI(KJ)
      ABR(IJ)=ABR(IJ)-ABRKJ*ELEMR+ABIKJ*ELEMI
   40 ABI(IJ)=ABI(IJ)-ABRKJ*ELEMI-ABIKJ*ELEMR
   50 CONTINUE
   60 CONTINUE
      
      
      SUMDR=0.Q0
      SUMDI=0.Q0
      SUMNR=0.Q0
      SUMNI=0.Q0
c      SUMDR2=0.Q0
c      SUMDI2=0.Q0
c      SUMNR2=0.Q0
c      SUMNI2=0.Q0
      DO 2 I=1,LG
      I1=LG*(I-1)
      DO 3 M=1,LG
      AAR=0.Q0
      AAI=0.Q0
c      AAR2=0.Q0
c      AAI2=0.Q0
      DO 4 K=1,LG
      K1=LG*(K-1)
      IK=I1+K
      KM=K1+M
      AMMRIK=AMMR(IK)
      AMMIIK=AMMI(IK)
      ABRKM=ABR(KM)
      ABIKM=ABI(KM)
C      AMMRKM=AMMR(KM)
C      AMMIKM=AMMI(KM)
C      ABRIK=ABR(IK)
C      ABIIK=ABI(IK)
      MAR=AMMRIK*ABRKM-AMMIIK*ABIKM
      MAI=AMMIIK*ABRKM+AMMRIK*ABIKM
c      MAR2=AMMRKM*ABRIK-AMMIKM*ABIIK
c      MAI2=AMMIKM*ABRIK+AMMRKM*ABIIK
      AAR=AAR+MAR
      AAI=AAI+MAI
C      AAI2=AAI2+MAI2
C      AAR2=AAR2+MAR2
    4 CONTINUE
      IF(I.EQ.M) SUMDR=SUMDR+ABS(AAR)
      IF(I.NE.M) SUMNR=SUMNR+ABS(AAR)
      IF(I.EQ.M) SUMDI=SUMDI+ABS(AAI)
      IF(I.NE.M) SUMNI=SUMNI+ABS(AAI)
c      IF(I.EQ.M) SUMDR2=SUMDR2+ABS(AAR2)
c      IF(I.NE.M) SUMNR2=SUMNR2+ABS(AAR2)
c      IF(I.EQ.M) SUMDI2=SUMDI2+ABS(AAI2)
c      IF(I.NE.M) SUMNI2=SUMNI2+ABS(AAI2)      
    3 CONTINUE
    2 CONTINUE
      SUMDR=SUMDR/LG
c      SUMDR2=SUMDR2/LG
c     
c      IF(ABS(SUMNR2).LE.EPS.AND.ABS(SUMDR2-ONEQ).LT.EPS
c     *.AND.ABS(SUMNI2).LE.EPS.AND.ABS(SUMDI2).LE.EPS) then
c          continue
c      else
c          print *,'$$ SUMS: ',SUMDR,SUMNR,SUMDI,SUMNI,SUMDR2,
c     *     SUMNR2,SUMDI2,SUMNI2
c
c      end if
c      
c      
      IF(ABS(SUMNR).LE.EPS.AND.ABS(SUMDR-ONEQ).LT.EPS
     *.AND.ABS(SUMNI).LE.EPS.AND.ABS(SUMDI).LE.EPS) go to 19 
      

C     ITTERATIONS TO MAKE INVERTED MATRIX ACCURATE
C     
      ITER=0
   12 ITER=ITER+1
   
      
      DO 5 I=1,LG
      I1G=(I-1)*LG       
      DO 9 M=1,LG
      M1G=(M-1)*LG
      AAR=0.Q0
      AAI=0.Q0
      DO 14 K=1,LG
      K1G=(K-1)*LG
      IK=I1G+K
      KM=K1G+M
      ABRKM=ABR(KM)
      ABIKM=ABI(KM)
      AMMRIK=AMMR(IK)
      AMMIIK=AMMI(IK)
      MAR=ABRKM*AMMRIK-ABIKM*AMMIIK
      MAI=ABRKM*AMMIIK+ABIKM*AMMRIK    
      AAR=AAR+MAR
      AAI=AAI+MAI
   14 CONTINUE
      IM=I1G+M
      AMAIR(IM)=AAR
      IF(I.EQ.M) AMAIR(IM)=AAR-ONEQ
      AMAII(IM)=AAI
    9 CONTINUE
    5 CONTINUE
    
      DO 6 I=1,LG
      I1G=(I-1)*LG       
      DO 7 M=1,LG
      M1G=(M-1)*LG
      AAR=0.Q0
      AAI=0.Q0
      DO 8 K=1,LG
      K1G=(K-1)*LG
      IK=I1G+K
      KM=K1G+M
      ABRKM=AMAIR(KM)
      ABIKM=AMAII(KM)
      AMRIK=ABR(IK)
      AMIIK=ABI(IK)
      MAR=ABRKM*AMRIK-ABIKM*AMIIK
      MAI=ABRKM*AMIIK+ABIKM*AMRIK    
      AAR=AAR+MAR
      AAI=AAI+MAI
    8 CONTINUE
      IM=I1G+M
      ABRC(IM)=AAR
      ABIC(IM)=AAI
    7 CONTINUE
    6 CONTINUE
      DO 10 KK=1,LG2
      ABR(KK)=ABR(KK)-ABRC(KK) 
   10 ABI(KK)=ABI(KK)-ABIC(KK)
      SUMDR=0.Q0
      SUMDI=0.Q0
      SUMNR=0.Q0
      SUMNI=0.Q0
      DO 27 I=1,LG
      I1=LG*(I-1)
      DO 33 M=1,LG
      AAR=0.Q0
      AAI=0.Q0
      DO 31 K=1,LG
      K1=LG*(K-1)
      IK=I1+K
      KM=K1+M    
      ABRIK=ABR(IK)
      ABIIK=ABI(IK)
      AMRKM=AMMR(KM)
      AMIKM=AMMI(KM)
      MAR=ABRIK*AMRKM-ABIIK*AMIKM
      MAI=ABRIK*AMIKM+ABIIK*AMRKM    
      AAR=AAR+MAR
      AAI=AAI+MAI       
   31 CONTINUE    
      IF(I.EQ.M) SUMDR=SUMDR+ABS(AAR)
      IF(I.NE.M) SUMNR=SUMNR+ABS(AAR)
      IF(I.EQ.M) SUMDI=SUMDI+ABS(AAI)
      IF(I.NE.M) SUMNI=SUMNI+ABS(AAI)
   33 CONTINUE
   27 CONTINUE
      SUMDR=SUMDR/LG
      SUMDI=SUMDI/LG
 
      IF(ABS(SUMNR).LE.EPS.AND.ABS(SUMDR-ONEQ).LT.EPS
     *.AND.ABS(SUMNI).LE.EPS.AND.ABS(SUMDI).LE.EPS) go to 19

        IF(ITER.LE.2) GO TO 12
        
      
      IF(MEPRI.LT.98) PRINT 24,SUMNR,SUMNI,SUMDR,SUMDI,INFOR,ITER 
      
      WRITE(21,24) SUMNR,SUMNI,SUMDR,SUMDI,INFOR,ITER 
C      PRINT 24,SUMNR,SUMNI,SUMDR,SUMDI,INFOR,ITER

   24 FORMAT(10X,'WARNING! MATRIX IS POORLY INVERTED'/
     *5X,'SUM OF RE NON-DIAG. ELEM-S=',D11.5,
     *5X,'SUM OF IM NON-DIAG. ELEM-S=',D11.5,/
     *5X,'SUM OF RE DIAGONAL=',D20.10,
     *5X,'SUM OF IM DIAGONAL=',D11.5, ',INFOR=',I2,',ITER=',I2)
     
   19 RETURN
      
      END
C     *******************************************************
      SUBROUTINE INVER
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      REAL*16 AMM,AMAI,AIAMAI
      REAL*16 PIV,ELEM,SUMD,SUMN,AA,A


      DIMENSION AMM(40000),AMAI(40000),AIAMAI(40000) 
      INCLUDE 'PRIVCOM.FOR' 
      INCLUDE 'PRIVCOM10.FOR'          
C     DET=1.0
      EPS=1.0D-7
      IR=1
      LG2=LG*LG
      DO 1 I=1,LG2
    1 AMM(I)=FI2(I)
      DO 60 K=1,LG
      K1G=(K-1)*LG
      KK=K1G+K
C     DET=DET*FI2(KK)
      PIV=1.0D+00/FI2(KK)
      DO 20 J=1,LG
      KJ=K1G+J
   20 FI2(KJ)=FI2(KJ)*PIV
      FI2(KK)=PIV
      DO 50 I=1,LG
      I1G=(I-1)*LG
      IK=I1G+K
      IF(I-K)30,50,30
   30 ELEM=FI2(IK)
      FI2(IK)=0.0
      DO 40 J=1,LG
      IJ=I1G+J
      KJ=K1G+J
   40 FI2(IJ)=FI2(IJ)-ELEM*FI2(KJ)
   50 CONTINUE
   60 CONTINUE
      SUMD=0.D+00
      SUMN=0.D+00
      DO 2 I=1,LG
      I1=LG*(I-1)
      DO 3 M=1,LG
      AA=0.D+00
      DO 4 K=1,LG
      K1=LG*(K-1)
      IK=I1+K
      KM=K1+M
      A=AMM(IK)*FI2(KM)
      AA=AA+A
    4 CONTINUE
      IF(I.EQ.M) SUMD=SUMD+ABS(AA)
      IF(I.NE.M) SUMN=SUMN+ABS(AA)
    3 CONTINUE
    2 CONTINUE
      SUMD=SUMD/LG
C     IF(MEPRI.LT.98) PRINT 99,SUMD,SUMN
C  99 FORMAT(10X,'M ATRIX ',2D20.8)
      IF(SUMN.LE.EPS.AND.ABS(SUMD-1.D0).LT.EPS) GO TO 19
C
C     ITTERATIONS TO MAKE INVERTED MATRIX ACCURATE
C
      ITER=0
   12 ITER=ITER+1
      DO 5 I=1,LG
      I1=LG*(I-1)
      DO 5 M=1,LG
      AA=0.D+00
      DO 14 K=1,LG
      K1=LG*(K-1)
      IK=I1+K
      KM=K1+M
      A=AMM(IK)*FI2(KM)
      AA=AA+A
   14 CONTINUE
      IM=I1+M
      AMAI(IM)=AA
    5 CONTINUE
      DO 6 I=1,LG
      I1=LG*(I-1)
      DO 6 M=1,LG
      AA=0.D+00
      DO 7 K=1,LG
      K1=LG*(K-1)
      IK=I1+K
      KM=K1+M
      A=FI2(IK)*AMAI(KM)
      AA=AA+A
    7 CONTINUE
      IM=I1+M
      AIAMAI(IM)=AA
    6 CONTINUE
      DO 8 I=1,LG2
    8 FI2(I)=2.D+00*FI2(I)-AIAMAI(I)
      SUMD=0.D+00
      SUMN=0.D+00
      DO 9 I=1,LG
      I1=LG*(I-1)
      DO 9 M=1,LG
      AA=0.D+00
      DO 10 K=1,LG
      K1=LG*(K-1)
      IK=I1+K
      KM=K1+M
      A=AMM(IK)*FI2(KM)
      AA=AA+A
   10 CONTINUE
      IF(I.EQ.M) SUMD=SUMD+ABS(AA)
      IF(I.NE.M) SUMN=SUMN+ABS(AA)
    9 CONTINUE
      SUMD=SUMD/LG
      IF(SUMN.LE.EPS.AND.ABS(SUMD-1.D0).LT.EPS) GO TO 19
      IF(ITER.LE.2) GO TO 12
      IF(MEPRI.LT.98) PRINT 24,SUMN,SUMD,INFOR,ITER
      WRITE(21,24)SUMN,SUMD,INFOR,ITER
   24 FORMAT(10X,'WARNING! MATRIX IS POORLY INVERTED'/
     *5X,'SUM OF NON-DIAG. ELEM-S=',D11.5,
     *', DIAGONAL=',D11.5,',INFOR=',I2,',ITER=',I2)
   19 RETURN
      END
C     *******************************************************
      SUBROUTINE PLEGA
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM1.FOR'
      DO 3 M=1,MTET
      TETA=TET(M)*4.D0*DATAN(1.D0)/180.D0
      PLEG(1,M)=1.D0
      PLEG(2,M)=DCOS(TETA)
      DO 3 K=3,180
      AK=K-1.
    3 PLEG(K,M)=((2.D0*AK-1.D0)*
     * DCOS(TETA)*PLEG(K-1,M)-(AK-1.D0)*PLEG(K-2,M))/AK
      RETURN
      END
C     *******************************************************
      SUBROUTINE DISCA
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM1.FOR'
      INCLUDE 'PRIVCOM13.FOR'
      INCLUDE 'PRIVCOM14.FOR'           

      IF(IDINT(ASP+0.1D0).EQ.1) THEN
      CALL DEUTR
      RETURN
      ENDIF
     
      CALL PLEGA
      LKK=1
      ETA=CONZ/WNK(1)
      DO 2 N=1,NMAX
      DO 4 L=1,180
      IF(L.LE.90.AND.N.EQ.1) COEFR(L)=0.D0
      IF(L.LE.90.AND.N.EQ.1) COEFI(L)=0.D0
    4 COEF(N,L)=0.D0
      DO 2 M=1,MTET
    2 DISC(N,M)=0.D0
      IC=JO(1)
      A1=200.D0*WN(1)*(IC+1.D0)
      DO 1 K1=1,NJ
      JS1=JS(K1)
      N1I=INC(K1)
      N1F=INR(K1)
      DO 1 N1C=1,N1I
      L1I=LNJ(K1,N1C)
      J1I=JNJ(K1,N1C)
      CORR=DCOS(2.D0*COPH(1,L1I+1))
      COII=DSIN(2.D0*COPH(1,L1I+1))
      ACRR=CRD(K1,N1C,N1C)
      ACII=CID(K1,N1C,N1C)
      ACRN=ACRR*CORR-ACII*COII
      ACIN=ACRR*COII+ACII*CORR
      DO 8 N1R=1,N1F
      L1F=LNJ(K1,N1R)
      J1F=JNJ(K1,N1R)
      NU1=NNJ(K1,N1R)
      COR=DCOS(COPH(1,L1I+1)+COPH(NU1,L1F+1))
      COI=DSIN(COPH(1,L1I+1)+COPH(NU1,L1F+1))
      IF(KODMA.GT.0) GO TO 11
      IF(NU1.LT.NUI) GO TO 8
      IF(NU1.GT.NUF) GO TO 17
      GO TO 13
   11 IF(NU1.LT.NUI) GO TO 8
      IF(NU1.GT.NUF) GO TO 8
   13 IR=JO(NU1)
      ACR=CRD(K1,N1C,N1R)
      ACI=CID(K1,N1C,N1R)
      ACR1=ACR*COR-ACI*COI
      ACI1=ACR*COI+ACI*COR
      A=(JS1+1.D0)*DSQRT((J1I+1.D0)*(J1F+1.D0))/A1
     **(-1)**((IR-IC)/2)
      A2=A
      DO 7 K2=K1,NJ
      IF(K2.NE.K1) A2=2.D0*A
      JS2=JS(K2)
      N2I=INC(K2)
      N2F=INR(K2)
      DO 7 N2C=1,N2I
      L2I=LNJ(K2,N2C)
      J2I=JNJ(K2,N2C)
      DO 5 N2R=1,N2F
      NU2=NNJ(K2,N2R)
      IF(KODMA.GT.0) GO TO 12
      IF(NU2.LT.NU1) GO TO 5
      IF(NU2.GT.NU1) GO TO 7
      GO TO 14
   12 IF(NU2.LT.NU1) GO TO 5
      IF(NU2.GT.NU1) GO TO 5
   14 L2F=LNJ(K2,N2R)
      J2F=JNJ(K2,N2R)
      COR=DCOS(COPH(1,L2I+1)+COPH(NU1,L2F+1))
      COI=DSIN(COPH(1,L2I+1)+COPH(NU1,L2F+1))
      LN=IABS(L2F-L1F)+1
      LK=L1F+L2F+1
      LL=L1I+L2I+1
      LLL=LL+LK
      IF(LLL/2*2.NE.LLL) GO TO 5
      LLN=IABS(L1I-L2I)+1
      IF(LLN.GT.LN) LN=LLN
      IF(LL.LT.LK) LK=LL
      LIO=IABS(J1I-J2I)/2+1
      LFO=(J1I+J2I)/2+1
      LLN=IABS(J1F-J2F)/2+1
      LL=(J1F+J2F)/2+1
      IF(LLN.GT.LIO) LIO=LLN
      IF(LL.LT.LFO) LFO=LL
      LLN=IABS(JS1-JS2)/2+1
      LL=(JS1+JS2)/2+1
      IF(LLN.GT.LIO) LIO=LLN
      IF(LL.LT.LFO) LFO=LL
      IF(LIO.GT.LN) LN=LN+(LIO-LN+1)/2*2
      IF(LFO.LT.LK) LK=LK-(LK-LFO+1)/2*2
      IF(LN.GT.LK) GO TO 5
      ACR=CRD(K2,N2C,N2R)
      ACI=CID(K2,N2C,N2R)
      ACR2=ACR*COR-ACI*COI
      ACI2=ACR*COI+ACI*COR
      B=A2*(JS2+1.D0)*DSQRT((J2I+1.D0)*(J2F+1.D0))*(ACR1*ACR2+ACI1*ACI2)
     **(-1)**((J1I+J2I+J1F+J2F)/2)
      IF(LK.GT.LKK) LKK=LK
      DO 6 L=LN,LK,2
      L1=L-1
      JA=JS1
      JB=J1F
      JE=JS2
      JD=J2F
      JC=IR
      JF=L1*2
      CALL RACAH
      AA=B*W
      JB=J1I
      JD=J2I
      JC=IC
      CALL RACAH
      AA=AA*W
      J1=J1I
      J2=J2I
      M1=1
      M2=-1
      J=JF
      M=0
      CALL KLEGO
      AA=AA*AKG
      J1=J1F
      J2=J2F
      CALL KLEGO
      AA=AA*AKG
      DO 3 M=1,MTET
    3 DISC(NU1,M)=DISC(NU1,M)+PLEG(L,M)*AA
      COEF(NU1,L)=COEF(NU1,L)+AA
    6 CONTINUE
    5 CONTINUE
    7 CONTINUE
    8 CONTINUE
   17 IF(ETA.EQ.0.D0) GO TO 1
      IF(NUI.GT.1) GO TO 1
      COEFR(L1I+1)=COEFR(L1I+1)+ACRN/A1*(JS1+1.D0)/(2.D0*L1I+1.D0)/A1*2
      COEFI(L1I+1)=COEFI(L1I+1)+ACIN/A1*(JS1+1.D0)/(2.D0*L1I+1.D0)/A1*2
      DO 18 M=1,MTET
      TETA=TET(M)*4.D0*DATAN(1.D0)/180.D0
      SIT22=DSIN(TETA/2.D0)**2
      ALST2=DLOG(DSIN(TETA/2.D0))
      ARGC=2.D0*(COPH(1,1)-ETA*ALST2)
C     COULOMB AMPLITUDE * by 2K
      COULR=-ETA/SIT22*DCOS(ARGC)
      COULI=-ETA/SIT22*DSIN(ARGC)
C     Cmat*Frez*+Frex*Cmat*=2.*ACC
      AAC=(COULR*ACRN+COULI*ACIN)*(JS1+1.D0)
   18 DISC(1,M)=DISC(1,M)+PLEG(L1I+1,M)/A1*AAC
    1 CONTINUE
      IF(ETA.EQ.0.D0) GO TO 16
      IF(NUI.GT.1) GO TO 16
      DO 15 M=1,MTET
      TETA=TET(M)*4.D0*DATAN(1.D0)/180.D0
   15 DISC(1,M)=DISC(1,M)+ETA**2/A1/2.D0/DSIN(TETA/2.D0)**4*(IC+1.D0)
   16 DO 10 N=NUI,NUF
      COE=COEF(N,1)
      DO 10 L=1,LKK
      COE=1.D0
   10 COEF(N,L)=COEF(N,L)/COE/(2*L-1)
      RETURN
      END
C     ***********************************************************
      SUBROUTINE ECISS
C     ***********************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
C     BELOW CARD IS NECESSARY IF MEMORY IS LESS THAN 32Mb
      DOUBLE PRECISION WPSL,PL,PSL,PR,PI,PVC,PRC
      REAL*16 ABR1,ABI1,ANR1,ANI1,ABR2,ABI2,ANR2,ANI2,
     *CFR1,CFI1,CFR2,CFI2,CRC,CIC,ARA,AIA,ZN
      REAL*16 CT,CY,APA,BPA,AL,AM,BEN,BTN
      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM7.FOR'  
      INCLUDE 'PRIVCOM10.FOR' 
      INCLUDE 'PRIVCOM12.FOR'
      INCLUDE 'PRIVCOM13.FOR'
      INCLUDE 'PRIVCOM14.FOR'
      

      
       
           
C     NEXT TWO CARDS JUST TO OVERCOME THE ERROR OR  MICROSOFT FPS COMPILER
C     v.1.00, YOU CAN DELETE THEM USING OTHER COMPILERS
      IF(MEPRI.LT.98) PRINT 999,(WNK(II),II=1,NUR)
  999 FORMAT (5E12.5)
      SPI2=ASP*(ASP+1.D0)
      CONG=(JSS+1)/(JO(1)+1.D0)
      CON=12.566372*CONG/WN(1)/100.D0/(NSPI+1)
      IF(MEPOT.EQ.1)   AKW=1.d0+BET(2)*2.5d0
      IF(MEPOT.NE.1)   AKW=1.d0+BET0*2.5d0
      NSLU=1
      IF(MEJOB.GT.3) NSLU=2
      LLNC=INCC
      NECI=1
   24 IF(NSLU.EQ.2) AKW=1.
      IF(NSLU.EQ.2) LLNC=NCLL
      IF(MEPOT.EQ.1) GO TO 30
      LAS2=LAS+1
      LAS1=LAS
      GO TO 36
   30 LAS2=(LAS+2)/2
      LAS1=LAS2-1
   36 NCLA=NCLL*LAS1
      LAS8=9
      LASC=1
      IF(LAS2.GE.3) LAS8=18
      IF(LAS2.GE.3) LASC=2
      ICLL=NCLL*LAS8
      MNULA=300*NUR*LAS2
      MLA=300*LAS2
      MNU=300*NUR
      ST=1.d0
      ST2=STEP*STEP
      ST12=ST2/12.d0
      ST6=ST2/6.d0
      DO 1 L=1,LLNC
      NCLLD=(L-1)*NCLL+L
      LLL=(L-1)*(NCLA+LAS1)-1
      LC2=(L-1)*(ICLL+LAS8)
      LL=LNJ1(L)
      ST=1.d0/(18.D+00)**LL
      JJ=JNJ1(L)
      NU=NNJ1(L)
      MNULLL=MNULA*(NU-1)+MLA*(NU-1)
      MNULL=MNU*(NU-1)+300*(NU-1)
C       IF((EN-EL(NU)/AT*(AT+ANEU)).LT.0) ST=EXP(-5.*LL)
      CJ1=JJ/2.
      SOP=(CJ1-LL)*(CJ1+LL+1)-SPI2
      FRH(1,L)=0.d0
      FIH(1,L)=0.d0
      FRH(2,L)=ST*STEP**(LL+1)
      FIH(2,L)=0.d0
      QSR1=0.d0
      QSI1=0.d0
      IF(LL.EQ.1) QSR1=-ST*ST6
      POTR=PSL(MNULL+2)*SOP+PL(2)*LL*(LL+1)-WN(NU)+PR(MNULLL+LAS2+1)
      POTI=PI(MNULLL+LAS2+1)*AKW+WPSL(MNULL+2)*SOP
      LI=LAS2+2+MNULLL
      LF=2*LAS2+MNULLL
      IF(MEVOL.GE.1.AND.MEPOT.EQ.1)
     * POTR=POTR+PR(LI)*CVNRV0(NCLLD)
      IF(MEVOL.GE.1.AND.MEPOT.EQ.1)
     * POTI=POTI+PI(LI)*AKW*CVNRV0(NCLLD)
      ACS=0.D0
      IF(LI.GT.LF) GO TO 32
      DO 2 LA=LI,LF
      A=CVNR(LLL+LA-LAS2-MNULLL)
      POTR=POTR+PR(LA)*A
      POTI=POTI+PI(LA)*A*AKW
      IF(MEPOT.EQ.1) GO TO 2
      LLLL=LA-LI+1
      IF(LLLL.GT.2) GO TO 2
      LAI=3
      LAF=5
      IF(LLLL.EQ.2) LAI=1
      IF(LLLL.EQ.2) LAF=9
      IF(LLLL.EQ.1) PVS=PR(LA)
      IF(LLLL.EQ.1) PWS=PI(LA)*AKW
      DO 62 LAC=LAI,LAF
      LALAS=LASC*(LAC-1)+LLLL
      LP1=LALAS+LAS8
      PV1=PRC(LP1)
      LL2=LC2+LALAS
      A=CVNC(LL2)
      IF(LLLL.EQ.2.AND.LAC.EQ.1) ACS=-A*CBET0
      IF(LAC.EQ.1) GO TO 62
      POTR=POTR+A*PV1
62    CONTINUE
    2 CONTINUE
      IF(MEVOL.EQ.0.OR.MEPOT.EQ.1) GO TO 32
      POTR=POTR+ACS*PVS
      POTI=POTI+ACS*PWS
  32  QSR2=FRH(2,L)*(1.d0-ST12*POTR)
      QSI2=-FRH(2,L)*ST12*POTI
      POI=POTI+ST6*POTR*POTI
      POR=POTR+ST12*(POTR*POTR-POTI*POTI)
      UR=(QSR2*POR-QSI2*POI)*ST2
      UI=(QSR2*POI+QSI2*POR)*ST2
      DO 3 NS=3,NHI
      QSR3=2.d0*QSR2-QSR1+UR
      QSI3=2.d0*QSI2-QSI1+UI
      LIC=(NS-1)*LAS8
      LI=(NS-1)*LAS2+1+MNULLL
      POTR=PSL(MNULL+NS)*SOP+PL(NS)*LL*(LL+1)-WN(NU)+PR(LI)
      POTI=PI(LI)*AKW+WPSL(MNULL+NS)*SOP
      LI=LI+1
      LF=NS*LAS2+MNULLL
      IF(MEVOL.GE.1.AND.MEPOT.EQ.1)
     * POTR=POTR+PR(LI)*CVNRV0(NCLLD)
       IF(MEVOL.GE.1.AND.MEPOT.EQ.1)
     * POTI=POTI+PI(LI)*AKW*CVNRV0(NCLLD)
      LP=LLL-LAS2*(NS-1)
      ACS=0.D0
      IF(LI.GT.LF) GO TO 34
      DO 4 LA=LI,LF
      A=CVNR(LP+LA-MNULLL)
      POTR=POTR+PR(LA)*A
      POTI=POTI+PI(LA)*A*AKW
      IF(MEPOT.EQ.1) GO TO 4
      LLLL=LA-LI+1
      IF(LLLL.GT.2) GO TO 4
      LAI=3
      LAF=5
      IF(LLLL.EQ.2) LAI=1
      IF(LLLL.EQ.2) LAF=9
      IF(LLLL.EQ.1) PVS=PR(LA)
      IF(LLLL.EQ.1) PWS=PI(LA)*AKW
      DO 65 LAC=LAI,LAF
      LALAS=LASC*(LAC-1)+LLLL
      LP1=LALAS+LIC
      PV1=PRC(LP1)
      LL2=LC2+LALAS
      A=CVNC(LL2)
      IF(LLLL.EQ.2.AND.LAC.EQ.1) ACS=-A*CBET0
      IF(LAC.EQ.1) GO TO 65
      POTR=POTR+A*PV1
65    CONTINUE
    4 CONTINUE
      IF(MEVOL.EQ.0.OR.MEPOT.EQ.1) GO TO 34
      POTR=POTR+ACS*PVS
      POTI=POTI+ACS*PWS
   34 POR=POTR+ST12*(POTR*POTR-POTI*POTI)
      POI=POTI+ST6*POTR*POTI
      UR=(QSR3*POR-QSI3*POI)*ST2
      UI=(QSR3*POI+QSI3*POR)*ST2
      FRH(NS,L)=QSR3+UR/12.d0
      FIH(NS,L)=QSI3+UI/12.d0
C     ANS=(NS-1.)*STEP
C     IF(MEPRI.LT.98) PRINT 111,FRH(NS,L),FIH(NS,L),ANS
      QSR1=QSR2
      QSI1=QSI2
      QSR2=QSR3
      QSI2=QSI3
    3 CONTINUE
      L1=LL+1
      ABR1=FBR1(NU,L1)
      ABI1=FBI1(NU,L1)
      ANR1=FNR1(NU,L1)
      ANI1=FNI1(NU,L1)
      ABR2=FBR2(NU,L1)
      ABI2=FBI2(NU,L1)
      ANR2=FNR2(NU,L1)
      ANI2=FNI2(NU,L1)
      CFR1=FRH(NHI-2,L)
      CFI1=FIH(NHI-2,L)
      CFR2=FRH(NHI,L)
      CFI2=FIH(NHI,L)
C     IF(MEPRI.LT.98) PRINT 222,L,LL,NCLL,NU
C 222 FORMAT (2X,20I5)
C     IF(MEPRI.LT.98) PRINT 111,CFR1,CFI1,CFR2,CFI2
  111 FORMAT (2X,2E13.5,4I10)
C     IF(MEPRI.LT.98) PRINT 111,ABR1,ABI1,ANR1,ANI1
C     IF(MEPRI.LT.98) PRINT 111,ABR2,ABI2,ANR2,ANI2
      CALL MATCH
      ZN=ARA*ARA+AIA*AIA
C     PRINT 1111,ARA,AIA,CRC,CIC,LL,L
C1111 FORMAT(4E10.2,4I3)
      CRH(L)=(ARA*CRC+AIA*CIC)/ZN
      CIH(L)=(ARA*CIC-AIA*CRC)/ZN
      IF(NSLU.EQ.2) GO TO 23
      CRS(L)=CRH(L)
      CIS(L)=CIH(L)
      DO 5 NS=1,NHI
      CIC=FIH(NS,L)
      CRC=FRH(NS,L)
      FRS(NS,L)=(ARA*CRC+AIA*CIC)/ZN
      FIS(NS,L)=(ARA*CIC-AIA*CRC)/ZN
    5 CONTINUE
      GO TO 1
   23 DO 25 NS=1,NHI
      CIC=FIH(NS,L)
      CRC=FRH(NS,L)
      FRH(NS,L)=(ARA*CRC+AIA*CIC)/ZN
      FIH(NS,L)=(ARA*CIC-AIA*CRC)/ZN
   25 CONTINUE
C 444    FORMAT(2X,3E15.7,9I5)
    1 CONTINUE
      NSLU=NSLU+1
      IF(NSLU.EQ.2) GO TO 24
      CSS=0.
C     IF(MEPRI.LT.98) PRINT 333
C333  FORMAT('    SECOND PART')
      MALL=NCLL
      IF(MESOL.LE.NCLL) MALL=MESOL
      DO 6 IC=1,INCC
      IC1=(IC-1)*MALL
      DO 9 L=1,NCLL
      CRI(L)=0.d0
      CII(L)=0.d0
      CRP(L)=0.d0
      CIP(L)=0.d0
      CRT(1,L)=0.d0
      CIT(1,L)=0.d0
      NPAD(L)=1
      DO 9 NS=1,NHI
      FRI(NS,L)=0.d0
    9 FII(NS,L)=0.d0
      IF(MESOL.GT.3) GO TO 47
      CRI(IC)=CRS(IC)
      CII(IC)=CIS(IC)
      CRP(IC)=CRS(IC)
      CIP(IC)=CIS(IC)
      CRT(1,IC)=CRS(IC)
      CIT(1,IC)=CIS(IC)
      DO 10 NS=1,NHI
      FRI(NS,IC)=FRS(NS,IC)
   10 FII(NS,IC)=FIS(NS,IC)
 1199 FORMAT (2E15.5)
      GO TO 48
   47 DO 49 MLL=1,MALL
      CRI(MLL)=CR(IC,MLL)
      CII(MLL)=CI(IC,MLL)
      CRP(MLL)=CR(IC,MLL)
      CIP(MLL)=CI(IC,MLL)
      CRT(1,MLL)=CR(IC,MLL)
      CIT(1,MLL)=CI(IC,MLL)
      DO 50 NS=1,NHI
      FRI(NS,MLL)=0.d0
      FII(NS,MLL)=0.d0
      DO 50 JLL=1,MALL
      JM=(JLL-1)*MALL+MLL
      MJ=(MLL-1)*MALL+JLL
      IJ=IC1+JLL
      FRI(NS,MLL)=FRI(NS,MLL)-AUR(IJ)*FREM(NS,JM)+AUI(IJ)*FIEM(NS,JM)
   50 FII(NS,MLL)=FII(NS,MLL)-AUR(IJ)*FIEM(NS,JM)-AUI(IJ)*FREM(NS,JM)
   49 CONTINUE
   48 ITER=1
      BEN=1.D0
   21 CSI=0.d0
      ITER=ITER+1
      IF(ITER.GT.25) GO TO 26
      IF(MEPOT.NE.1) BEN=BEN*BET0
      IF(MEPOT.EQ.1) BEN=BEN*BET(2)
      DO 7 L=1,NCLL
      LLL=(L-1)*(NCLA+LAS1)-1
      
      LISO=(L-1)*NCLL
      
      NCLLD=(L-1)*NCLL+L

      LC2=(L-1)*(ICLL+LAS8)
      LL=LNJ1(L)
      ST=1./(18.D+00)**LL
      JJ=JNJ1(L)
      NU=NNJ1(L)
      MNULLL=MNULA*(NU-1)+MLA*(NU-1)
      MNULL=MNU*(NU-1)+300*(NU-1)
      MNLL=MNULA*(NU-1)
      MNL=MNU*(NU-1)
C     IF((EN-EL(NU)/AT*(AT+ANEU)).LT.0) ST=EXP(-5.*LL)
      NPI1=KPJ1(L)
C     NEXT TWO CARDS JUST TO OVERCOME THE ERROR OR  MICROSOFT FPS COMPILER
C     v.1.00, YOU CAN DELETE THEM USING OTHER COMPILERS
      IF(WNK(NU)/WNK(1).LE.0.) PRINT 99,WNK(NU),WNK(1),NU
   99 FORMAT (2E12.3,I3)
      C=SQRT(WNK(NU)/WNK(1))
      CJ1=JJ/2.d0
      SOP=(CJ1-LL)*(CJ1+LL+1)-SPI2
      FRI(1,L)=0.d0
      FII(1,L)=0.d0
      FRI(2,L)=ST*STEP**(LL+1)
      FII(2,L)=0.d0
      QSR1=0.d0
      QSI1=0.d0
      IF(LL.EQ.1) QSR1=-ST*ST6
      POTR=PSL(MNULL+2)*SOP+PL(2)*LL*(LL+1)-WN(NU)+PR(MNULLL+LAS2+1)
      POTI=PI(MNULLL+LAS2+1)+WPSL(MNULL+2)*SOP

      LI=LAS2+2+MNULLL
      LF=2*LAS2+MNULLL
      IF(MEVOL.GE.1.AND.MEPOT.EQ.1)
     * POTR=POTR+PR(LI)*CVNRV0(NCLLD)
      IF(MEVOL.GE.1.AND.MEPOT.EQ.1)
     * POTI=POTI+PI(LI)*CVNRV0(NCLLD)
        ACS=0.D0
      IF(LI.GT.LF) GO TO 38
      DO 8 LA=LI,LF
      A=CVNR(LLL+LA-LAS2-MNULLL)
      POTR=POTR+PR(LA)*A
      POTI=POTI+PI(LA)*A
      IF(MEPOT.EQ.1) GO TO 8
      LLLL=LA-LI+1
      IF(LLLL.GT.2) GO TO 8
      LAI=3
      LAF=5
      IF(LLLL.EQ.2) LAI=1
      IF(LLLL.EQ.2) LAF=9
      IF(LLLL.EQ.1) PVS=PR(LA)
      IF(LLLL.EQ.1) PWS=PI(LA)
      DO 68 LAC=LAI,LAF
      LALAS=LASC*(LAC-1)+LLLL
      LP1=LALAS+LAS8
      PV1=PRC(LP1)
      LL2=LC2+LALAS
      A=CVNC(LL2)
      IF(LLLL.EQ.2.AND.LAC.EQ.1) ACS=-A*CBET0
      IF(LAC.EQ.1) GO TO 68
      POTR=POTR+A*PV1
   68 CONTINUE
    8 CONTINUE
      IF(MEVOL.EQ.0.OR.MEPOT.EQ.1) GO TO 38
      POTR=POTR+ACS*PVS
      POTI=POTI+ACS*PWS
   38 QSR2=FRI(2,L)*(1.-ST12*POTR)
      QSI2=-FRI(2,L)*ST12*POTI
      POR=POTR+ST12*(POTR*POTR-POTI*POTI)
      POI=POTI+ST6*POTR*POTI
      UR=(QSR2*POR-QSI2*POI)*ST2
      UI=(QSR2*POI+QSI2*POR)*ST2
      WR1=0.d0
      WI1=0.d0
      WR2=0.d0
      WI2=0.d0
      K1=(L-1)*NCLA-1
      KC1=ICLL*(L-1)
      DO 11 IR=1,NCLL
C        IF(L.LE.MALL.AND.IR.GT.MALL) GO TO 11
C        IF(L.GT.MALL.AND.IR.GT.1) GO TO 11
C        IF(L.GT.MALL) GO TO 11
      NPI2=KPJ1(IR)
      NUIR=NNJ1(IR)
      LISOR=LISO+IR
      MNLLIR=MNLL+MLA*(NUIR-1)
      MNLIR=MNL+300*(NUIR-1)
      K2=K1+(IR-1)*LAS1
      KC2=KC1+LAS8*(IR-1)
      LI=LAS2+2+MNLLIR
      IF(IR.EQ.L) GO TO 11
C      LI=LAS2+2+MNLLIR
      LF=2*LAS2+MNLLIR
      POTR=0.d0
      POTI=0.d0
      ACS=0.d0
      IF(LI.GT.LF) GO TO 42
      DO 12 LA=LI,LF
      A=CVNR(K2+LA-LAS2-MNLLIR)
      IF(NPI1.NE.NPI2) GO TO 39
      POTR=POTR+PR(LA)*A
      POTI=POTI+PI(LA)*A
      GO TO 70
   39 POTR=POTR-PI(LA)*A
      POTI=POTI+PR(LA)*A
   70 IF(MEPOT.EQ.1) GO TO 12
      LLLL=LA-LI+1
      IF(LLLL.GT.2) GO TO 12
      LAI=3
      LAF=5
      IF(LLLL.EQ.2) LAI=1
      IF(LLLL.EQ.2) LAF=9
      IF(LLLL.EQ.1) PVS=PR(LA)
      IF(LLLL.EQ.1) PWS=PI(LA)
      DO 71 LAC=LAI,LAF
      LALAS=LASC*(LAC-1)+LLLL
      LP1=LALAS+LAS8
      PV1=PRC(LP1)
      LL2=KC2+LALAS
      A=CVNC(LL2)
      IF(LLLL.EQ.2.AND.LAC.EQ.1) ACS=-A*CBET0
      IF(LAC.EQ.1) GO TO 71
      IF(NPI1.NE.NPI2) GO TO 72
      POTR=POTR+A*PV1
      GO TO 71
   72 POTI=POTI+A*PV1
   71 CONTINUE
   12 CONTINUE
      
      COUPL=CVNRPN(LISOR)
C      IF(NU.NE.NUIR.AND.JO(NU).EQ. JO(NUIR).AND.JO(NU).EQ.0) COUPL=1.D0
      IF(NCA(NU).NE.NCA(NUIR).AND.JO(NU).EQ. JO(NUIR))
     *POTR=POTR+PR(LI-1)*COUPL
      IF(NCA(NU).NE.NCA(NUIR).AND.JO(NU).EQ. JO(NUIR))
     *POTI=POTI+PI(LI-1)*COUPL

      IF(L.NE.IR.AND.MEVOL.GE.1.AND.MEPOT.EQ.1)
     * POTR=POTR+PR(LI)*CVNRV0(LISOR)
      IF(L.NE.IR.AND.MEVOL.GE.1.AND.MEPOT.EQ.1)
     * POTI=POTI+PI(LI)*CVNRV0(LISOR)
     
   
     
      IF(MEVOL.EQ.0.AND.MEPOT.EQ.1) GO TO 42
      POTR=POTR+ACS*PVS
      POTI=POTI+ACS*PWS
   42 FR=FRI(2,IR)
      FI=FII(2,IR)
      WR2=WR2+POTR*FR-POTI*FI
      WI2=WI2+POTR*FI+POTI*FR
   11 CONTINUE
      DO 13 NS=3,NHI
      WR3=0.d0
      WI3=0.d0
      DO 14 IR=1,NCLL
C        IF(L.LE.MALL.AND.IR.GT.MALL) GO TO 14
C        IF(L.GT.MALL.AND.IR.GT.1) GO TO 14
C        IF(L.GT.MALL) GO TO 14
      LISOR=LISO+IR
      NPI2=KPJ1(IR)
      NUIR=NNJ1(IR)
      MNLLIR=MNLL+MLA*(NUIR-1)
      MNLIR=MNL+300*(NUIR-1)
      K2=K1+(IR-1)*LAS1
      KC2=KC1+LAS8*(IR-1)
C      IF(IR.EQ.L) GO TO 14
      LI=(NS-1)*LAS2+2+MNLLIR
      IF(IR.EQ.L) GO TO 14
      LIC=(NS-1)*LAS8
      LC2=(L-1)*(ICLL+LAS8)
      POTR=0.d0
      POTI=0.d0
      LF=NS*LAS2+MNLLIR
      LP=K2-(NS-1)*LAS2
      ACS=0.d0
      IF(LI.GT.LF) GO TO 45
      DO 15 LA=LI,LF
      A=CVNR(LP+LA-MNLLIR)
      IF(NPI1.NE.NPI2) GO TO 40
      POTR=POTR+PR(LA)*A
      POTI=POTI+PI(LA)*A
      GO TO 73
   40 POTR=POTR-PI(LA)*A
      POTI=POTI+PR(LA)*A
   73 IF(MEPOT.EQ.1) GO TO 15
      LLLL=LA-LI+1
      IF(LLLL.GT.2) GO TO 15
      LAI=3
      LAF=5
      IF(LLLL.EQ.2) LAI=1
      IF(LLLL.EQ.2) LAF=9
      IF(LLLL.EQ.1) PVS=PR(LA)
      IF(LLLL.EQ.1) PWS=PI(LA)
      DO 79 LAC=LAI,LAF
      LALAS=LASC*(LAC-1)+LLLL
      LP1=LALAS+LIC
      PV1=PRC(LP1)
      LL2=KC2+LALAS
      A=CVNC(LL2)
      IF(LLLL.EQ.2.AND.LAC.EQ.1) ACS=-A*CBET0
      IF(LAC.EQ.1) GO TO 79
      IF(NPI1.NE.NPI2) GO TO 82
      POTR=POTR+A*PV1
      GO TO 79
   82 POTI=POTI+A*PV1
   79 CONTINUE
   15 CONTINUE
            
      COUPL=CVNRPN(LISOR)
C     IF(NU.NE.NUIR.AND.JO(NU).EQ. JO(NUIR).AND.JO(NU).EQ.0) COUPL=1.D0
      IF(NCA(NU).NE.NCA(NUIR).AND.JO(NU).EQ. JO(NUIR))
     *POTR=POTR+PR(LI-1)*COUPL
      IF(NCA(NU).NE.NCA(NUIR).AND.JO(NU).EQ. JO(NUIR))
     *POTI=POTI+PI(LI-1)*COUPL
     
      IF(L.NE.IR.AND.MEVOL.GE.1.AND.MEPOT.EQ.1)
     * POTR=POTR+PR(LI)*CVNRV0(LISOR)
      IF(L.NE.IR.AND.MEVOL.GE.1.AND.MEPOT.EQ.1)
     * POTI=POTI+PI(LI)*CVNRV0(LISOR)
 
C        IF(NU.NE.NUIR.AND.JO(NU).EQ. JO(NUIR))
C     *print 567, lisor,nu,nuir,jo(nu),jo(nuir),li,PR(LI-1),Pi(LI-1),
C     *     CVNRPN(LISOR)
C  567 Format (6i4,5e12.3)
C      IF(NU.NE.NUIR.AND.JO(NU).EQ. JO(NUIR))pause

      IF(MEVOL.EQ.0.OR.MEPOT.EQ.1) GO TO 45
      POTR=POTR+ACS*PVS
      POTI=POTI+ACS*PWS
   45 FR=FRI(NS,IR)
      FI=FII(NS,IR)
      WR3=WR3+POTR*FR-POTI*FI
      WI3=WI3+POTR*FI+POTI*FR
   14 CONTINUE
      QSR3=2.*QSR2-QSR1+UR+(WR3+10.*WR2+WR1)*ST12
      QSI3=2.*QSI2-QSI1+UI+(WI3+10.*WI2+WI1)*ST12
      LI=(NS-1)*LAS2+1+MNULLL
      POTR=PSL(MNULL+NS)*SOP+PL(NS)*LL*(LL+1)-WN(NU)+PR(LI)
      POTI=PI(LI)+WPSL(MNULL+NS)*SOP
      LI=LI+1
      LF=NS*LAS2+MNULLL
      
      IF(MEVOL.GE.1.AND.MEPOT.EQ.1)
     * POTR=POTR+PR(LI)*CVNRV0(NCLLD)
      IF(MEVOL.GE.1.AND.MEPOT.EQ.1)
     * POTI=POTI+PI(LI)*CVNRV0(NCLLD)

      LP=LLL-LAS2*(NS-1)
      ACS=0.d0
      IF(LI.GT.LF) GO TO 46
      DO 16 LA=LI,LF
      A=CVNR(LP+LA-MNULLL)
      POTR=POTR+PR(LA)*A
      POTI=POTI+PI(LA)*A
      IF(MEPOT.EQ.1) GO TO 16
      LLLL=LA-LI+1
      IF(LLLL.GT.2) GO TO 16
      LAI=3
      LAF=5
      IF(LLLL.EQ.2) LAI=1
      IF(LLLL.EQ.2) LAF=9
          IF(LLLL.EQ.1) PVS=PR(LA)
          IF(LLLL.EQ.1) PWS=PI(LA)
      DO 74 LAC=LAI,LAF
      LALAS=LASC*(LAC-1)+LLLL
      LP1=LALAS+LIC
      PV1=PRC(LP1)
      LL2=LC2+LALAS
      A=CVNC(LL2)
          IF(LLLL.EQ.2.AND.LAC.EQ.1) ACS=-A*CBET0
          IF(LAC.EQ.1) GO TO 74
      POTR=POTR+A*PV1
   74 CONTINUE
   16 CONTINUE
      IF(MEVOL.EQ.0.OR.MEPOT.EQ.1) GO TO 46
          POTR=POTR+ACS*PVS
          POTI=POTI+ACS*PWS
   46 POR=POTR+ST12*(POTR*POTR-POTI*POTI)
      POI=POTI+ST6*POTR*POTI
      UR=(QSR3*POR-QSI3*POI)*ST2
      UI=(QSR3*POI+QSI3*POR)*ST2
      FRI(NS,L)=QSR3+UR/12.
      FII(NS,L)=QSI3+UI/12.
      WR1=WR2
      WI1=WI2
      WR2=WR3
      WI2=WI3
      QSR1=QSR2
      QSI1=QSI2
      QSR2=QSR3
      QSI2=QSI3
   13 CONTINUE
      L1=LL+1
      ABR1=FBR1(NU,L1)
      ABI1=FBI1(NU,L1)
      ANR1=FNR1(NU,L1)
      ANI1=FNI1(NU,L1)
      ABR2=FBR2(NU,L1)
      ABI2=FBI2(NU,L1)
      ANR2=FNR2(NU,L1)
      ANI2=FNI2(NU,L1)
      CFR1=FRI(NHI-2,L)
      CFI1=FII(NHI-2,L)
      CFR2=FRI(NHI,L)
      CFI2=FII(NHI,L)
C
C     IF(MEPRI.LT.98) THEN
C       PRINT 222,L,LL,NCLL,NU
C       PRINT 111,CFR1,CFI1,CFR2,CFI2
C       PRINT 111,ABR1,ABI1,ANR1,ANI1
C       PRINT 111,ABR2,ABI2,ANR2,ANI2
C     ENDIF  
C
      CALL MATCH
      IF(IC.EQ.L) ARA=ARA-1.d0
      CHR=CRH(L)
      CHI=CIH(L)
C     IF(MEPRI.LT.98) THEN
C       PRINT 1111,ARA,AIA,CRC,CIC,LL,L,ITER
C       PRINT 1111,ARA,AIA,CRC,CIC,CHR,CHI,L
C     ENDIF  
C1111 FORMAT(6E10.2,3I2)
      CRC=(CRC-ARA*CHR+AIA*CHI)*C
      CIC=(CIC-ARA*CHI-AIA*CHR)*C
      CRT(ITER,L)=(CRC-CRI(L))/BEN
      CIT(ITER,L)=(CIC-CII(L))/BEN
      IF(CRT(ITER,L).NE.0..AND.CIT(ITER,L).NE.0.) GO TO 28
C     IF(MEPRI.LT.98) PRINT 111,CRT(ITER,L),CIT(ITER,L),L,NSS,ITER,NCLL
      IF(ITER.GT.3) GO TO 31
      CR(IC,L)=CRC
      CI(IC,L)=CIC
  31  NPAD(L)=0
 28   CRI(L)=CRC
      CII(L)=CIC
      DO 17 NS=1,NHI
      CRC=FRI(NS,L)
      CIC=FII(NS,L)
      FR=FRH(NS,L)
      FI=FIH(NS,L)
      FRI(NS,L)=CRC-ARA*FR+AIA*FI
      FII(NS,L)=CIC-ARA*FI-AIA*FR
   17 CONTINUE
    7 CONTINUE
      IF(ITER.GE.3.AND.ITER/2*2.NE.ITER) GO TO 18
      GO TO 21
   18 SUMC=0.d0
      SUMI=0.d0
      SUMT=0.d0
      LPA=ITER/2
      LPA1=LPA+1
      MPA=LPA
C     IF(MEPRI.LT.98) PRINT 1999,ITER,MPA,LPA
C1999 FORMAT(3I6)
      DO 20 L=1,INCR
C         GO TO 27
      IF(NPAD(L).EQ.0) GO TO 27
      DO 19 IT=1,ITER
   19 CT(IT)=CRT(IT,L)
C     IF(MEPRI.LT.98) PRINT 1999,ITER,MPA,LPA
      CALL PADE
      AL=0.d0
      AM=0.d0
      DO 22 LM=1,LPA1
      IF(MEPOT.NE.1) BTN=BET0**(LM-1)
      IF(MEPOT.EQ.1) BTN=BET(2)**(LM-1)
      AL=AL+APA(LM)*BTN
   22 AM=AM+BPA(LM)*BTN
C         IF(DLOG(DABS(AL))-DLOG(DABS(AM)).GT.0.D0) GO TO 26
      CR(IC,L)=AL/AM
      DO 29 IT=1,ITER
   29 CT(IT)=CIT(IT,L)
      CALL PADE
      AL=0.d0
      AM=0.d0
      DO 33 LM=1,LPA1
      IF(MEPOT.NE.1) BTN=BET0**(LM-1)
      IF(MEPOT.EQ.1) BTN=BET(2)**(LM-1)
      AL=AL+APA(LM)*BTN
   33 AM=AM+BPA(LM)*BTN
C        IF(DLOG(DABS(AL))-DLOG(DABS(AM)).GT.0.D0) GO TO 26
      CI(IC,L)=AL/AM
   27 CRC=CR(IC,L)
      CIC=CI(IC,L)
      CSI=CSI+CRC*CRC+CIC*CIC
      SUMC=SUMC+ABS(CRP(L)-CRC)+ABS(CIP(L)-CIC)
      IF(IC.EQ.L) SUMI=SUMI+ABS(CIP(L)-CIC)
      IF(IC.EQ.L) SUMT=SUMT+CIC*CON
      IF(ABS(CRC).GT.1.D+3) GO TO 26
      IF(ABS(CIC).GT.1.D+3) GO TO 26
      CRD(NSS,IC,L)=CR(IC,L)
      CID(NSS,IC,L)=CI(IC,L)
      CRP(L)=CRC
   20 CIP(L)=CIC
      CSUR=CSI*CON
      ACON=CON
      IF(MEPRI.LT.98) PRINT 456,SUMC,SUMI,ACON,CSI,CSUR,SUMT
      
CNaN       WRITE(TID1,136) SUMC,CSI,CSUR
CNaN       IF(MEPRI.LT.98) PRINT 136, SUMC,CSI,CSUR
CNaN  136  FORMAT(3E12.5)
          
       
       
CNaN       BACKSPACE TID1
CNaN       READ (TID1,37)CCNAN,CINAN,CRNAN
CNaN   37  FORMAT(3A4)
C             
CNaN       IF(MEPRI.LT.98) PRINT 37, CCNAN,CINAN,CRNAN
C       PAUSE 33
CNaN       IF(CINAN.EQ.NAN) GO TO 26
CNaN       IF(CCNAN.EQ.NAN) GO TO 26
CNaN       IF(CRNAN.EQ.NAN) GO TO 26      
      
C      IF(SUMC /= SUMC) GO TO 26
C      IF(CSI /= CSI)   GO TO 26
C      IF(CSUR /= CSUR) GO TO 26   

       IF(SUMC.NE.SUMC) GO TO 26
       IF(CSI.NE.CSI)   GO TO 26
       IF(CSUR.NE.CSUR) GO TO 26       
      
      IF(ITER.EQ.3) GO TO 21

      EPECI=0.00005d0/CON
      IF(SUMI.GT.EPECI) GO TO 21
C     RCN: CHECK
C     A division by zero occurs if optimization allowed in MSF
C
      EPECI=0.00005D0/CON/2/SQRT(CSI)
      IF(SUMC.GT.EPECI) GO TO 21
C
C     HIGH ANGULAR MOMENTUM CASE : SCATTERING EQUAL TOTAL, TO
C     AVOID NEGATIVE TRANSMISSIONS AS ABSORPTION IS VERY SMALL.
C
      IF(SUMT.GE.0.AND.CSI.LE.CI(IC,IC)) GO TO 77
      CSID=CSI-CI(IC,IC)*CI(IC,IC)
      CI(IC,IC)=0.5D0*(1.-SQRT(1.D0-4.*CSID))
      CID(NSS,IC,IC)=CI(IC,IC)
      CSI=CI(IC,IC)
C
   77 CSS=CSS+CSI
  456 FORMAT(2X,4E8.2,2E14.7)
    6 CONTINUE
      IF(MEPRI.LT.98) PRINT 123,NSS,NCLL,ITER
  123 FORMAT(3I5,E20.7)
      RETURN
   26 NECI=0
      RETURN
      END
C     ***********************************************************
      SUBROUTINE MATCH
C     ***********************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL*16 ABR1,ABI1,ANR1,ANI1,ABR2,ABI2,ANR2,ANI2,
     *CFR1,CFI1,CFR2,CFI2,CRC,CIC,ARA,AIA 
      INCLUDE 'PRIVCOM.FOR'
      ZNR=ABR2*ANR1-ABI2*ANI1-ABR1*ANR2+ABI1*ANI2
      ZNI=ABI2*ANR1+ABR2*ANI1-ABR1*ANI2-ABI1*ANR2
      IF(ZNI.NE.0.) GO TO 1
      CRC=(CFR1*ABR2-CFR2*ABR1)/ZNR
      CIC=(CFI1*ABR2-CFI2*ABR1)/ZNR
      ARA=(CFR2*ANR1-CFI2*ABR1-CFR1*ANR2+CFI1*ABR2)/ZNR
      AIA=(CFI2*ANR1+CFR2*ABR1-CFR1*ABR2-CFI1*ANR2)/ZNR
      GO TO 2
    1 CRC=(CFR1*ABI2+CFI1*ABR2-CFR2*ABI1-CFI2*ABR1)/ZNI
      CIC=(CFI1*ABI2+CFR2*ABR1-CFR1*ABR2-CFI2*ABI1)/ZNI
      FWR1=ANR1-ABI1
      FWI1=ANI1+ABR1
      FWR2=ANR2-ABI2
      FWI2=ANI2+ABR2
      ARA=(CFR2*FWI1+CFI2*FWR1-CFR1*FWI2-CFI1*FWR2)/ZNI
      AIA=(CFI2*FWI1-CFR2*FWR1-CFI1*FWI2+CFR1*FWR2)/ZNI
    2 RETURN
      END
C     ***********************************************************
      SUBROUTINE PADE
C     ***********************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      REAL*16 AMA,CY,CT,BP,BPA,APA
      INCLUDE 'PRIVCOM.FOR'
     
C     PRINT 100,CT
C     PRINT 200,MPA,LPA,KKK
C 200 FORMAT(5I10)
      IF(MPA.LE.0) GO TO 9
      DO 1 I=1,MPA
      I1=(I-1)*MPA
      DO 1 J=1,MPA
      IJ=I1+J
      LMPA=LPA-MPA+I+J
      IF(LMPA.LT.1) GO TO 2
      AMA(IJ)=CT(LMPA)
      GO TO 1
    2 AMA(IJ)=0.D0
    1 CONTINUE
      MAM=MPA
      CALL MATIN
      DO 3 I=1,MPA
    3 CY(I)=CT(LPA+I+1)
      DO 4 K=1,MPA
      K1=(K-1)*MPA
      MK=MPA+2-K
      BP=0.D0
      DO 5 I=1,MPA
      KI=K1+I
    5 BP=BP-AMA(KI)*CY(I)
    4 BPA(MK)=BP
    9 BPA(1)=1.D0
      LPA1=LPA+1
      MPA1=MPA+1
      DO 6 L=1,LPA1
      AP=CT(L)
      LMA=L
      IF(L.GT.MPA1) LMA=MPA1
      IF(LMA.LT.2) GO TO 8
      DO 7 I=2,LMA
    7 AP=AP+BPA(I)*CT(L-I+1)
    8 APA(L)=AP
    6 CONTINUE
C     PRINT 100,APA,BPA
C 100 FORMAT (2X,6E20.7)
      RETURN
      END
C     ***********************************************************
      SUBROUTINE MATIN
C     ***********************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL*16 AMA,PIV,ELEM
      INCLUDE 'PRIVCOM.FOR'
      DO 60 K=1,MAM
      K1=(K-1)*MAM
      KK=K1+K
      PIV=1.D0/AMA(KK)
      DO 20 J=1,MAM
      KJ=K1+J
   20 AMA(KJ)=AMA(KJ)*PIV
      AMA(KK)=PIV
      DO 50 I=1,MAM
      I1=(I-1)*MAM
      IK=I1+K
      IF(I-K) 30,50,30
   30 ELEM=AMA(IK)
      AMA(IK)=0.D0
      DO 40 J=1,MAM
      IJ=I1+J
      KJ=K1+J
   40 AMA(IJ)=AMA(IJ)-ELEM*AMA(KJ)
   50 CONTINUE
   60 CONTINUE
      RETURN
      END
C     ***********************************************************
      SUBROUTINE RIPAT
C     ***********************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
C     BELOW CARD IS NECESSARY IF MEMORY IS LESS THAN 32Mb
      DOUBLE PRECISION VL,VSL,PV,PW,WSL,WPSL,PL,PSL,PR,PI,PVC,PRC
      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM3.FOR'
      INCLUDE 'PRIVCOM6.FOR' 
      INCLUDE 'PRIVCOM10.FOR' 
      INCLUDE 'PRIVCOM13.FOR'
      INCLUDE 'PRIVCOM14.FOR'

      PVV=0.D0
      PWW=0.D0
      PRR=0.D0
      PII=0.D0
      PSL=0.D0
      WPSL=0.D0
      VSL=0.D0
      WSL=0.D0
      PV =0.D0
      PW =0.D0
      PR =0.D0
      PI =0.D0

      CCDE=0.D0
      CCCOUL=CCOUL
      IF(MECUL.GE.2) CCDE=CCOUL
      IF(MECUL.GE.2) CCOUL=0.D0
      CDE=CCDE*MECHA*ZNUC/AT**(1.D0/3.D0)
      CDE12=CCDE*ZNUC/AT**(1.D0/3.D0)
 
      VRDC=0.D0
      VD=0.D0
      VDISP=0.D0
      CARBUN=931.49378D0
      AMI=939.56536D0
      IF(MECHA.EQ.1) AMI=938.272029D0
     
      IF(MECHA.EQ.1.AND.IDINT(ASP+0.1D0).EQ.1) AMI=1875.612859D0 
      IF(MECHA.EQ.1.AND.IDINT(ASP+0.1D0).EQ.1) ANEU=2.013553212712D0 
      
      DAVERN=0.D0
      DAVERP=0.D0
      EFERM=EFERMN
      EPAVER=EFERM+DAVERN
CCHA  IF(MECHA.EQ.1) EFERM=EFERMP
      IF(MECHA.EQ.1) EFERM=EFERMP+CDE 
      IF(MECHA.EQ.1) EPAVER=EFERM+DAVERP
      EINT=EN
      
c     DON'T CALCULATE CORRECTLY VOLUME INTEGRALS WITHOUT LOWER LINE!!!
      PRINT *,''
c
      
CCHA  EN=EN-EFERM-CDE
      EN=EN-EFERM
          
      IF(MEPRI.LT.98) CALL POTVOL
      EN=EINT

      RZIN=RZ
      IF(MERZZ.EQ.0) GO TO 52
      IF(MERZZ.EQ.1)   RZ=RZ*(1.D00-RZBWC*(EN-EFERM)**PDIS/
     *((EN-EFERM)**PDIS+RZWID**PDIS))

   52 RRR=RR
      APAN=ANEU/1.00866491600D0 
C
C           RR ENERGY DEPENDENCE
C
      IF(MERRR.EQ.1)   RR=RR*(1.D00-RRBWC*(EN-EFERM)**PDIS/
     *((EN-EFERM)**PDIS+RRWID**PDIS))
      DO 1 I=1,NUR
      REL=(EN+AMI)/AMI
      IF(MEREL.EQ.0) REL=1.D+0
      E=EN-EL(I)/AT*(AT+ANEU*REL)
      REL=(ABS(E)+AMI)/AMI
      IF(MEREL.EQ.0.OR.E.LE.0.) REL=1.D+0
      IF(MEREL.NE.0.AND.E.GT.0.) E=(REL**2-1.D+0)*AMI/2.D+0

      WNK(I)=0.219677D0*SQRT(ABS(E))*AT/SQRT(AT**2+2.D+0*ANEU*REL*AT+
     *ANEU**2)*SQRT(APAN)
      WN(I)=WNK(I)*WNK(I)
      IF(EN-EL(I)/AT*(AT+ANEU*REL).LE.0.) WN(I)=-WN(I)
    1 CONTINUE

      
      REL=(EN+AMI)/AMI
      IF(MEREL.EQ.0) REL=1.D+0
      IF(REL.NE.1.) RELPOT=2.D+0*(EN+AMI)/(2.D+0*AMI+EN)

      EIRELM=SQRT(REL**2*AT**2+2.D0*REL*ANEU*AT+ANEU**2)
      EIRELT=AT*SQRT(REL**2*ANEU**2+2.D0*REL*ANEU*AT+AT**2)
      DEREL=SQRT(2.D0*REL*ANEU*AT+AT**2+ANEU**2)
      WW=4.8257984D-2*EIRELM*EIRELT/(ANEU*EIRELM+EIRELT)/DEREL*APAN
      IF(MEREL.EQ.1.OR.MEREL.EQ.3) WW=WW*RELPOT
      CONZ=MECHA*ZNUC*EIRELM*EIRELT/(ANEU*EIRELM+EIRELT)/DEREL*
     *3.458814365D-2*APAN
      CONZ12=ZNUC*EIRELM*EIRELT/(ANEU*EIRELM+EIRELT)/DEREL*
     *3.458814365D-2*APAN

      CONZZ=CONZ
         
CCHA      ENEF=EN-EFERM-CDE
CCHA      EN=EN-CDE
      ENEF=EN-EFERM
      EN=EN
      
      VISO=CISO*(AT-2.D0*ZNUC)/AT
      WISO=WCISO*(AT-2.D0*ZNUC)/AT
      WVISO=WDISO*(AT-2.D0*ZNUC)/AT
      IF(MECHA.EQ.0) VISO=-VISO
      IF(MECHA.EQ.0) WISO=-WISO
      IF(MECHA.EQ.0) WVISO=-WVISO
      VRLANP=VRLA+VISO

C     WCP=Wv(Ef-Ea)=Wv(Ef+Ea)

      WCP=(WCBW+WVISO)*(EA)**PDIS/
     *(EA**PDIS+WCWID**PDIS)

C     TWO OPTIONS FOR COULOMB CORRECTIONS
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      VRDIR=((VR1+2.*VR2*ENEF+3.*VR3*ENEF**2
     *-VRLA*ALAVR*EXP(-ALAVR*ENEF)))*(1.D0+VISO/(VR0+VRLA))
     **(-1.D0)
     
      IF(MEDIS.EQ.0) GO TO 391
C    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     VHF DERIVATIVE AT EN
      EEE=EN
      EN=EN-0.005D+00
      ENCON=EN

      CALL VHFROM
      VHFM=VHF

      EN=EN+0.01D+00
      ENCON=EN

      CALL VHFROM
      VHFP=VHF

      VRDIR=(VHFM-VHFP)*100
      EN=EEE

C    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      VDISP=DOM_INT_WV(Eferm,EPAVER,WCBW+WVISO,WCWID,INT(PDIS),EN,VDISD)

      T12D = 0.D0
      IF(EA.GT.0.d0.AND.WC.NE.0.D0) THEN
       T12P=WCP*DOM_INT_T1(Eferm,EA,EN+0.01)+
     *ALFNEW*DOM_INT_T2(EFERM,EA,EN+0.01)
       T12M=WCP*DOM_INT_T1(Eferm,EA,EN-0.01)+
     *ALFNEW*DOM_INT_T2(EFERM,EA,EN-0.01)
       T12D=(T12M-T12P)*50
      ENDIF
      IF(MECUL.GE.2.OR.MECHA.EQ.0) T12D=0.D0
      VRDIRC=VDISD+T12D

      VDP=DOM_INT_WS
     * (Eferm,EPAVER,WDBW+WISO,WDWID,ALAWD,INT(PDIS),EN,VDD)

      VDCUL=VDD*MECHA*ZNUC/AT**0.333333333D0*CCOUL

      PDIS1=PDIS-1
C
      WDUL=0.
      WCUL=0.

      WCUL=0.
C
      VCULC=MECHA*ZNUC/AT**0.333333333D0*CCOUL*VRDIRC
  391 VCUL=MECHA*ZNUC/AT**0.333333333D0*CCOUL*VRDIR

      IF(MECUL.EQ.1) VCUL=MECHA*ZNUC/AT**0.333333333D0*CCOUL
C
      MUP=-1.91301D0
      IF(MECHA.NE.0) MUP=2.2928D0
      MUP=0.D0
C
      VR=(VR0+VR1*ENEF+VR2*ENEF**2+VR3*ENEF**3
     *+VRLA*EXP(-ALAVR*ENEF))*(1.D00+VISO/(VR0+VRLA))+VCUL

      IF(MEDIS.EQ.0) GO TO 392
C    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     VHF AND VR AT EN

      ENCON=EN

      CALL VHFROM
      VR=VHF+VCUL


C    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      VDIS=VDISP

      VD=VDP+VDCUL

   
      IF (MEPRI.LT.98) WRITE(21,787) 
     * EN,VRLA,ALAVR,ENEF,VISO,VCUL,VR0, VR, VD, VDP, VDCUL,
     * VDIS, VDISD,VDD,T12D

  787 FORMAT(/2X,' POTENTIAL VALUES (I) FOR INCIDENT ENERGY=',F11.6,1X,
     *'MeV:'//,
     *3X,'VRLA=',F8.3,3X,'ALAVR=',F8.3,3X,
     *                       '   Ef=',F8.3,3X,'VISO=',F8.3/
     *3X,'VCUL=',F8.3,3X,'  VR0=',F8.3,3X,
     *                       '   VR=',F8.3/
     *3X,'  VD=',F8.3,3X,'  VDP=',F8.3,3X,
     *                       'VDCUL=',F8.3/
     *3X,'VDIS=',F8.3,3X,'VDISD=',F8.3,3X,
     *                       '  VDD=',F8.3,3X,'T12D=',F8.3/)

C     VR = VR + VDIS
      VRDC=VDIS+VCULC

  392 CONTINUE


      IF(EN.LT.BNDC) GO TO 800
      WD=WD0+WD1*BNDC+(ENEF-BNDC)*WDA1+(WDBW+WISO)*EXP(-ALAWD*ENEF)
     **ENEF**PDIS/(ENEF**PDIS+WDWID**PDIS)+WDUL
      WC=WC0+WC1*BNDC+(ENEF-BNDC)*WCA1+(WCBW+WVISO)*ENEF**PDIS/
     *(ENEF**PDIS+WCWID**PDIS)+WCUL
      GO TO 801
  800 WD=WD0+WD1*ENEF+(WDBW+WISO)*EXP(-ALAWD*ENEF)
     **ENEF**PDIS/(ENEF**PDIS+WDWID**PDIS)+WDUL
      WC=WC0+WC1*ENEF+(WCBW+WVISO)*ENEF**PDIS/
     *(ENEF**PDIS+WCWID**PDIS)+WCUL




C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  801 REL=(EN+AMI)/AMI
      IF(REL.LT.1.D00) REL=1.D00
      EIRELM=SQRT(REL**2*AT**2+2.D+0*REL*ANEU*AT+ANEU**2)
      EIRELT=AT*SQRT(REL**2*ANEU**2+2.D+0*REL*ANEU*AT+AT**2)
      DEREL=SQRT(2.D+0*REL*ANEU*AT+AT**2+ANEU**2)
C     ENCON=(REL**2-1.D0)*AMI*AT*(ANEU+AT)/DEREL**2/2.D+0
      WW4=1.2064496D-2*EIRELM*EIRELT/(ANEU*EIRELM+EIRELT)/DEREL*APAN

C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      WW4=1.D0
      VSO=VS*EXP(-WW4*ALASO*(ENEF+CDE))
      IF(MEDIS.EQ.0) GO TO 491
      IF(MEDIS.EQ.2) GO TO 492

CCHA  DVSO=DOM_INT_WV(Eferm,EPAVER,WSBW,WSWID,INT(PDIS),EN+CDE,dtmp)
      DVSO=DOM_INT_WV(Eferm-CDE,EPAVER,WSBW,WSWID,INT(PDIS),EN,dtmp)
      VSO=VSO+DVSO
 
  492 IF(EA.GT.0.D0.AND.WC.NE.0.D0) THEN
      VRDC=WCP*DOM_INT_T1(Eferm,EA,EN)+ALFNEW*DOM_INT_T2(EFERM,EA,EN)
     *+VRDC
      IF(EN.GT.EFERM+EA) 
     *WC=WC+ALFNEW*(SQRT(EN)+(EFERM+EA)**1.5D0/
     *(2.D0*EN)-1.5D0*SQRT(EFERM+EA))
      ENDIF

  491 CONTINUE
      WSO=WS0+WS1*ENEF+WSBW*(ENEF+CDE)**PDIS/
     *((ENEF+CDE)**PDIS+WSWID**PDIS)
      IF(WD.LT.0.) WD=0.D0
      IF(WC.LT.0.) WC=0.D0

      EN=EINT
      IF (MEPRI.LT.98) 
     *WRITE(21,789) EN,VR,VRDC,VD,WC,WD,VSO,DVSO,WSO,VHF,VCUL,VCULC,
     *RELPOT
      IF(MEPRI.LT.98)  
     *PRINT 789,EN,VR,VRDC,VD,WC,WD,VSO,DVSO,WSO,VHF,VCUL,VCULC,RELPOT
  789 FORMAT(/2X,' POTENTIAL VALUES (II) FOR INCIDENT ENERGY=',F11.6,1X,
     *'MeV:'//,
     *3X,' VR=',F8.3,3X,'VRDC(DWv)=',F8.3,3X,
     *                   '  DWs=',F8.3,3X,' Wv=',F8.3/
     *3X,' Ws=',F8.3,3X,'      VSO=',F8.3,3X,
     *                   ' DVSO=',F8.3,3X,'WSO=',F8.3/
     *3X,'VHF=',F8.3,3X,'     VCUL=',F8.3,3X,
     *                   'VCULC=',F8.3,3X,'RELPOT=',F8.3/)

      AR=AR0+AR1*ENEF
      AW=AW0+AW1*ENEF
      AC=AC0+AC1*ENEF
      AD=AD0+AD1*ENEF
      IF (ENEF.GT.BNDC) AD=AD0+AD1*BNDC
      AS=AS0+AS1*ENEF


      RK=RR+(LOG(VR/(EN-EL(NMAX)))+10.)*AR
      
      IF(WD.EQ.0.) RKD=RK
      IF(WD.EQ.0.) GO TO 802
      RKD=RD+(LOG(WD/(EN-EL(NMAX)))+10.D0)*AD
  802 IF(RKD.GT.RK) RK=RKD
c     IF(MECHA.NE.0) RK=RK*4.D0
C     IF(MECHA.NE.0) RK=RK+ 3.0*RR
c     IF(RK.GT.30.) RK=30.
      RKC=3.D0*(MECHA*ZNUC*RR**2/EN)**0.3333333333D0
      IF(RKC.GT.RK) RK=RKC

      IF(MECHA.NE.0) RK=RK

      NH=RK/AR/0.3D0+1
      NH1=RK/AD/0.8D0+1
      NH2=RK/AC/0.8D0+1
      NH3=RK*WNK(1)/0.5+1
      IF(NH1.GT.NH) NH=NH1
      IF(NH2.GT.NH) NH=NH2
      IF(NH3.GT.NH) NH=NH3
      IF(NH.GT.270) NH=270
      LAS1=LAS+1
      LAS2=(LAS+2)/2
      IF(MEPOT.EQ.1) GO TO 39
      LAS2=LAS+1
   39 IF(MEPOT.EQ.1) GO TO 40
      GO TO 41
   40 IF(NPD.EQ.0) GO TO 41
      AN=0.
      DO 2 L=2,NPD,2
    2 AN=AN+ABS(BET(L))*L*LAS2
      NAN1=AN*25
      NAN1=2*NAN1+2
      IF(NAN1.LT.40) NAN1=40
      IF(NAN1.GT.298) NAN1=298
      NAN2=NAN1+1
      ST=1./(NAN2-1)
      DO 3 N=1,NAN2
      X=ST*(N-1)
      A=0.2820947917739D0
      B=X*A
      POL(1,N)=A
      DEF(N)=1.D0
      AVOL=1.D0
C      IF(MEDEF.EQ.1.AND.MEVOL.GE.1) AVOL=1.D0+BET(2)**2
      DO 4 K=3,LAS1,2
      KK=K-1
      C=((2.*KK-1.D0)*X*B-(KK-1.D0)*A)/KK
C     POL - SPHERICAL FUNCTIONS L-1,ANGLE
      POL((KK+2)/2,N)=C*SQRT(2.D0*KK+1.D0)
      IF(KK.LE.NPD) DEF(N)=DEF(N)+BET(KK)*C*SQRT(2.D0*KK+1.D0)/AVOL
      A=B
      B=C
      KK=KK+1
      C=((2.*KK-1.D0)*X*B-(KK-1.D0)*A)/KK
      A=B
      B=C
    4 CONTINUE
    3 CONTINUE
   41 PIM=1.996852D0
      MNULA=300*NUR*LAS2
      LALAS=9
      LASC=1
      IF(LAS2.GE.3) LALAS=18
      IF(LAS2.GE.3) LASC=2
      MLA=300*LAS2
      MNU=300*NUR
      DO 38 NU1=1,NUR
      MNULA1=MNULA*(NU1-1)
      MNU1=MNU*(NU1-1)
      DO 38 NU2=1,NUR
      MNULA2=MNULA1+MLA*(NU2-1)
      MNU2=MNU1+300*(NU2-1)
      EN=EINT-(EL(NU1)+EL(NU2))/2
C     if(nu1.ne.nu2)EN=EINT-ABS(EL(NU1)-EL(NU2))/2.
      if(nu1.ne.nu2)EN=EINT-(EL(NU1)+EL(NU2))/2.D0
      IF(MEAPP.EQ.1) EN=EINT

C     SPHERICAL PART OF THE POTENTIAL AND COUPLING IN PARTITIONS  !!!!!!!!
      
      IF(NCA(NU1).NE.NCA(NU2)) GO TO 38

      VISO=CISO*(AT-2.D0*ZNUC)/AT
      WISO=WCISO*(AT-2.D0*ZNUC)/AT
      WVISO=WDISO*(AT-2.D0*ZNUC)/AT
      IF(MECHA.EQ.0) VISO=-CISO*(AT-2.D0*ZNUC)/AT
      IF(MECHA.EQ.0) WISO=-WCISO*(AT-2.D0*ZNUC)/AT
      IF(MECHA.EQ.0) WVISO=-WDISO*(AT-2.D0*ZNUC)/AT

      IF(MECHA.EQ.1.AND.NCA(1).NE.NCA(NU1))VISO=-CISO*
     * (AT-2.D0*ZNUC-2.D0)/AT
      IF(MECHA.EQ.1.AND.NCA(1).NE.NCA(NU1))WISO=-WCISO*
     * (AT-2.D0*ZNUC-2.D0)/AT
      IF(MECHA.EQ.1.AND.NCA(1).NE.NCA(NU1))WVISO=-WDISO*
     * (AT-2.D0*ZNUC-2.D0)/AT
      
C      WCP=Wv(Ef-Ea)=Wv(Ef+Ea)

      WCP=(WCBW+WVISO)*(EA)**PDIS/
     *(EA**PDIS+WCWID**PDIS)

      IF(MECHA.EQ.1.AND.NCA(1).NE.NCA(NU1)) CONZ=0.D0
      IF(MECHA.EQ.0.AND.NCA(1).NE.NCA(NU1)) CONZ=CONZ12
      IF(MECHA.EQ.1.AND.NCA(1).EQ.NCA(NU1)) CONZ=CONZ12
      IF(MECHA.EQ.0.AND.NCA(1).NE.NCA(NU1)) CONZ=0.D0
      
      IF(MECHA.EQ.1.AND.NCA(1).NE.NCA(NU1)) EFENU12=EFERMN
      IF(MECHA.EQ.0.AND.NCA(1).NE.NCA(NU1)) EFENU12=EFERMP
      IF(MECHA.EQ.1.AND.NCA(1).EQ.NCA(NU1)) EFENU12=EFERMP
      IF(MECHA.EQ.0.AND.NCA(1).EQ.NCA(NU1)) EFENU12=EFERMN
      
      IF(MECHA.EQ.1.AND.NCA(1).NE.NCA(NU1)) ENEF=EN-EFENU12
      IF(MECHA.EQ.1.AND.NCA(1).NE.NCA(NU1)) EN=EN

CCHA      IF(MECHA.EQ.1.AND.NCA(1).EQ.NCA(NU1)) ENEF=EN-EFENU12-CDE12
CCHA      IF(MECHA.EQ.1.AND.NCA(1).EQ.NCA(NU1)) EN=EN-CDE12
      IF(MECHA.EQ.1.AND.NCA(1).EQ.NCA(NU1)) ENEF=EN-EFENU12
      IF(MECHA.EQ.1.AND.NCA(1).EQ.NCA(NU1)) EN=EN
      
CCHA     IF(MECHA.EQ.0.AND.NCA(1).NE.NCA(NU1)) ENEF=EN-EFENU12-CDE12
CCHA      IF(MECHA.EQ.0.AND.NCA(1).NE.NCA(NU1)) EN=EN-CDE12
      IF(MECHA.EQ.0.AND.NCA(1).NE.NCA(NU1)) ENEF=EN-EFENU12
      IF(MECHA.EQ.0.AND.NCA(1).NE.NCA(NU1)) EN=EN

      IF(MECHA.EQ.0.AND.NCA(1).EQ.NCA(NU1)) ENEF=EN-EFENU12
      IF(MECHA.EQ.0.AND.NCA(1).EQ.NCA(NU1)) EN=EN



C     TWO OPTIONS FOR COULOMB CORRECTIONS
C
      
      COENH=1.D0

 
  876 VRLANP=VRLA+COENH*VISO

      VRDIR=((VR1+2.*VR2*ENEF+3.*VR3*ENEF**2
     *-VRLA*ALAVR*EXP(-ALAVR*ENEF)))*(1.D00+VISO/(VR0+VRLA))
     **(-1.D0)

      IF(MEDIS.EQ.0) GO TO 394


C     VHF DERIVATIVE AT EN
      EEE=EN
      EN=EN-0.005D+00
      ENCON=EN

      CALL VHFROM
      VHFM=VHF

      EN=EN+0.01D+00
      ENCON=EN

      CALL VHFROM
      VHFP=VHF

      VRDIR=(VHFM-VHFP)*100
      EN=EEE

C     DISPERSIVE WC CONTRIBUTION TO VR

      VDISP=DOM_INT_WV(Eferm,EPAVER,WCBW+WVISO,WCWID,INT(PDIS),EN,VDISD)

      T12D=0.d0
      IF(EA.GT.0.d0.AND.WC.NE.0.D0) THEN
      T12P=WCP*DOM_INT_T1(Eferm,EA,EN+0.01D0)+
     *ALFNEW*DOM_INT_T2(EFERM,EA,EN+0.01D0)
      T12M=WCP*DOM_INT_T1(Eferm,EA,EN-0.01D0)+
     *ALFNEW*DOM_INT_T2(EFERM,EA,EN-0.01D0)
      T12D=(T12M-T12P)*50.D0
      ENDIF
      IF(MECUL.GE.2.OR.MECHA.EQ.0) T12D=0.D0
      VRDIRC=VDISD+T12D

      VCULC=MECHA*ZNUC/AT**0.333333333D0*CCOUL*VRDIRC
      WDUL=0.D0
      WCUL=0.D0
      
C     DISPERSIVE WD CONTRIBUTION TO VR

  984 CONTINUE
      
      IF(ALAWD.NE.0.D0) THEN    
        VDP=DOM_INT_WS
     * (Eferm,EPAVER,WDBW+WISO,WDWID,ALAWD,INT(PDIS),EN,VDD)
      ELSE
        VDP=DOM_INT_WV(Eferm,EPAVER,WDBW+WISO,WDWID,INT(PDIS),EN,VDD)
      ENDIF

      IF(EN.GT.WDSHI.AND.ALAWD.EQ.0.D+00) THEN 
        VDP1=VDP
        VDD1=VDD
        VDP=DOM_INT_WV(Eferm,WDSHI,WDBW+WISO,WDWID2,INT(PDIS),EN,VDD)
        VDP=VDP1-VDP
        VDD=VDD1-VDD 
      ENDIF      

      IF(EN.GT.EA) VDP=VDP+VRD*(EN-EA)
      IF(EN.GT.EA) VDD=VDD+VRD

      VDCUL=VDD*MECHA*ZNUC/AT**0.333333333*CCOUL
C
       PDIS1=PDIS-1
C
      WDUL=0.

      WCUL=0.

  394 VCUL=MECHA*ZNUC/AT**0.333333333*CCOUL*VRDIR
      IF(MECUL.EQ.1) VCUL=MECHA*ZNUC/AT**0.333333333*CCOUL

      VR=(VR0+VR1*ENEF+VR2*ENEF**2+VR3*ENEF**3
     *+VRLA*EXP(-ALAVR*ENEF))*(1.D00+VISO/(VR0+VRLA))+VCUL

      IF(MEDIS.EQ.0) GO TO 395

C    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     VHF AND VR AT EN
      ENCON=EN
      CALL VHFROM
      VR=VHF+VCUL
      VRDC=VDISP+VCULC
C     ENEF=ENCON-EFERM

C    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      VD=VDP+VDCUL

  395 CONTINUE
      ENE=ENEF
      IF(MECUL.EQ.3) ENE=ENE+CDE

      ENESHID=ENE-WDSHI+EFERM                  !!!!!!NEEDS to be CHECKED
      
      IF(ENE.LT.BNDC) GO TO 900
      WD=WD0+WD1*BNDC+(ENE-BNDC)*WDA1+(WDBW+WISO)*EXP(-ALAWD*ENE)
     **ENE**PDIS/(ENE**PDIS+WDWID**PDIS)+WDUL
      IF(ENE+EFERM.GT.WDSHI.AND.ALAWD.EQ.0.D+00)WD=WD-(WDBW+WISO)
     **ENESHID**PDIS/(ENESHID**PDIS+WDWID2**PDIS)
      WC=WC0+WC1*BNDC+(ENE-BNDC)*WCA1+(WCBW+WVISO)*ENE**PDIS/
     *(ENE**PDIS+WCWID**PDIS)+WCUL
      GO TO 901
  900 WD=WD0+WD1*ENE+(WDBW+WISO)*EXP(-ALAWD*ENE)
     **ENE**PDIS/(ENE**PDIS+WDWID**PDIS)+WDUL
      IF(ENE+EFERM.GT.WDSHI.AND.ALAWD.EQ.0.D+00)WD=WD-(WDBW+WISO)
     **ENESHID**PDIS/(ENESHID**PDIS+WDWID2**PDIS)

      WC=WC0+WC1*ENE+(WCBW+WVISO)*ENE**PDIS/
     *(ENE**PDIS+WCWID**PDIS)+WCUL

CCCCC 
      
  901 CONTINUE


C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  972 REL=(EN+AMI)/AMI
      IF(REL.LT.1.D00) REL=1.D00
      EIRELM=SQRT(REL**2*AT**2+2.D+0*REL*ANEU*AT+ANEU**2)
      EIRELT=AT*SQRT(REL**2*ANEU**2+2.D+0*REL*ANEU*AT+AT**2)
      DEREL=SQRT(2.D+0*REL*ANEU*AT+AT**2+ANEU**2)
C        ENCON=(REL**2-1.D0)*AMI*AT*(ANEU+AT)/DEREL**2/2.D+0
      WW4=1.2064496D-2*EIRELM*EIRELT/(ANEU*EIRELM+EIRELT)/DEREL*APAN

C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       WW4=1.
CCHA       VSO=VS*EXP(-WW4*ALASO*(ENEF+CDE))
       VSO=VS*EXP(-WW4*ALASO*(ENEF))

      IF(MEDIS.EQ.0) GO TO 396
      IF(MEDIS.EQ.2) GO TO 397


CCHA      DVSO=DOM_INT_WV(Eferm,EPAVER,WSBW,WSWID,INT(PDIS),EN+CDE,dtmp)
      DVSO=DOM_INT_WV(Eferm,EPAVER,WSBW,WSWID,INT(PDIS),EN,dtmp)                   ???
       VSO=VSO+DVSO

 397  IF(EA.GT.0.D0.AND.WC.NE.0.D0) THEN
      VRDC=WCP*DOM_INT_T1(Eferm,EA,EN)+ALFNEW*DOM_INT_T2(EFERM,EA,EN)
     *+VRDC
      IF(EN.GT.EFERM+EA) 
     *WC=WC+ALFNEW*(SQRT(EN)+(EFERM+EA)**1.5D0/
     *(2.D0*EN)-1.5D0*SQRT(EFERM+EA))
      ENDIF
  396 CONTINUE
 
CCHA      WSO=WS0+WS1*ENEF+WSBW*(ENEF+CDE)**PDIS/
CCHA     *((ENEF+CDE)**PDIS+WSWID**PDIS)
      IF(MECUL.EQ.3) ENESO=ENEF+CDE
       WSO=WS0+WS1*ENESO+WSBW*(ENESO)**PDIS/
     *((ENESO)**PDIS+WSWID**PDIS)    
C     WD=WD+WISO
      IF(WD.LE.0.)WD=0.
      IF(WC.LE.0.)WC=0.
      AR=AR0+AR1*ENEF
      AW=AW0+AW1*ENEF
      AC=AC0+AC1*ENEF
      AD=AD0+AD1*ENEF
      IF (ENEF.GT.BNDC) AD=AD0+AD1*BNDC
      AS=AS0+AS1*ENEF
      IF(MEPOT.EQ.1) GO TO 17
      RAR=RR/AR
      RAC=RC/AC
      RAW=RW/AW
      RAD=RD/AD
   17 M=8
      M4=4*M
      KL=0
      IN=1
      KIT=M
      NIT=1
      IK=M4-1
C
      STEP=RK/NH
      STEP1=STEP/M
C
      DO 5 K=1,2
      DO 6 I=IN,IK
      LALASI=LALAS*(I-1)
      L1=(I-1)*LAS2+MNULA2
      L1AP=(I-1)*LAS2

C     Added by RCN, June 2014
      IR=I/KIT+NIT
      LALASR=LALAS*(IR-1)
      LALASI=LALAS*(I-1)
C     Added by RCN, June 2014

      IF(K.EQ.2.AND.I.EQ.IK) GO TO 14
      R=STEP1*(I-KL)
      VSLC=-MUP*CONZ/RZ**3*PIM
      IF(R.GT.RZ) VSLC=-MUP*CONZ/R**3*PIM
      X=(R-RS)/AS
      II=MNU2+I
      IF(X.GT.23) GO TO 7
      EX1=EXP(X)
      VSLF=-PIM*EX1/AS/R/(1.D0+EX1)/(1.D0+EX1)
      GO TO 8
    7 EX1=EXP(-X)
      VSLF=-PIM*EX1/AS/R
    8 VL(I)=1.D0/R/R
      VSL(II)=VSLF*VSO*WW
      WSL(II)=VSLF*WSO*WW
      IF(MEREL.EQ.1.OR.MEREL.EQ.3) VSL(II)=VSLF*VSO*WW/RELPOT
      IF(MEREL.EQ.1.OR.MEREL.EQ.3) WSL(II)=VSLF*WSO*WW/RELPOT
      IF(I/KIT*KIT.NE.I) GO TO 15
      IR=I/KIT+NIT
      LALASR=LALAS*(IR-1)
      I1=(IR-1)*LAS2+MNULA2
      IIC=MNU2+IR
      PL(IR)=VL(I)
      PSL(IIC)=VSL(II)
      WPSL(IIC)=WSL(II)
      GO TO 15
   14 R=0.
      IR=1
      I1=MNULA2
      IIC=MNU2+1
   15 IF(AZ.GT.0.) GO TO 47
      PZ=CONZ/RZ*(3.-R*R/RZ/RZ)
      IF(R.GE.RZ) PZ=2.*CONZ/R
      GO TO 46
   47 CALL SPHEPOT
      PZ=2.*DE*CONZ
   46 IF(MEPOT.EQ.1) GO TO 18
      IF(LAS2.EQ.1) GO TO 21
      DO 20 LAC=1,9
      LALALI=LALASI+LASC*(LAC-1)
      LALALR=LALASR+LASC*(LAC-1)
      LAM=LAC-1
      IF(R.EQ.0.) RLAM=0.0
      IF(R.EQ.0.) RLAM1=0.0
      IF(R.NE.0.) RLAM=R**LAM
      IF(R.NE.0.) RLAM1=R**(LAM+1)
      PVC(LALALI+1)=6.*CONZ/(2.*LAM+1)*RLAM/RZ**(LAM+1)
      IF(R.GT.RZ) PVC(LALALI+1)=6.*CONZ/(2.*LAM+1)*RZ**LAM/RLAM1
C     PVC(LALALI+1)=0.0
      IF(LAS2.EQ.2) GO TO 20
      PVC(LALALI+2)=6.*CONZ/(2.*LAM+1)*(1.-LAM)*RLAM/RZ**(LAM+1)
C     PVC(LALALI+2)=6.*CONZ/(2.*LAM+1)*(LAM+2.)*RLAM/RZ**(LAM+1)
      IF(R.GT.RZ) PVC(LALALI+2)=6.*CONZ/(2.*LAM+1)*(LAM+2)*RZ**LAM/
     *RLAM1
C     PVC(LALALI+2)=0.0
   20 CONTINUE
   21 NT0=L1+1
      NT1=L1+2
      NT2=L1+3
      NT3=L1+4
      NT4=L1+5
      X=(R-RR)/AR
      IF(X.GT.23.) GO TO 30
      EX=EXP(X)
      EX1=1.D0+EX
      EXD=1.D0/EX1
      EXDM=1.D0-EXD
      PV(NT0)=-VR*EXD*WW+PZ
      IF(LAS.EQ.0) GO TO 31
      PV(NT1)=(PV(NT0)-PZ)*EXDM*RAR
      IF(LAS.EQ.1) GO TO 31
      Y1=1.D0-2.D0*EXD
      PV(NT2)=PV(NT1)*Y1*RAR/2.D0
      IF(LAS.EQ.2) GO TO 31
      Y2=1.D0-6.D0*EXD*EXDM
      PV(NT3)=PV(NT2)/Y1*Y2*RAR/3.D0
      IF(LAS.EQ.3) GO TO 31
      Y3=1.-EXD*(14.D0-EXD*(36.D0-24.D0*EXD))
      PV(NT4)=PV(NT3)/Y2*Y3*RAR/4.D0
      GO TO 31
   30 EX1=EXP(-X)
      PV(NT0)=-VR*EX1*WW+PZ
      IF(LAS.EQ.0) GO TO 31
C     PV(NT1)=PV(NT0)*RAR
      PV(NT1)=(PV(NT0)-PZ)*RAR
      IF(LAS.EQ.1) GO TO 31
      PV(NT2)=PV(NT1)*RAR/2.D0
      IF(LAS.EQ.2) GO TO 31
      PV(NT3)=PV(NT2)*RAR/3.D0
      IF(LAS.EQ.3) GO TO 31
      PV(NT4)=PV(NT3)*RAR/4.D0
   31 X=(R-RC)/AC
      Y=(R-RW)/AW
      XX=Y*Y
      IF(XX.GT.180.) EXX=0.
      IF(XX.LE.180.) EXX=EXP(-XX)
      IF(X.GT.23.) GO TO 32
      EX=EXP(X)
      EX1=1.D0+EX
      EXD=1.D0/EX1
      EXDM=1.D0-EXD
      WCC=-WC*ALF*EXD*WW
      WWW=-WC*(1.D0-ALF)*EXX*WW
      IF(MEREL.EQ.3) WCC=WCC/RELPOT
      IF(MEREL.EQ.3) WWW=WWW/RELPOT
      PW(NT0)=WCC+WWW
      PV(NT0)=PV(NT0)+WCC*VRDC/WC
      IF(LAS.EQ.0) GO TO 33
      WCC=WCC*EXDM*RAC
      WWW=WWW*2.*RAW*Y
      PW(NT1)=WCC+WWW
       PV(NT1)=PV(NT1)+WCC*VRDC/WC
      IF(LAS.EQ.1) GO TO 33
      YC1=1.-2.*EXD
      WCC=WCC*YC1*RAC/2.
      Y1=2.*Y*Y-1.
      WWW=WWW/Y*Y1*RAW/2.
      PW(NT2)=WCC+WWW
       PV(NT2)=PV(NT2)+WCC*VRDC/WC
      IF(LAS.EQ.2) GO TO 33
      YC2=1.-6.*EXD*EXDM
      WCC=WCC/YC1*YC2*RAC/3.
      Y2=Y**3*2.-Y*3.
      WWW=WWW/Y1*Y2*RAW*2./3.
      PW(NT3)=WCC+WWW
       PV(NT3)=PV(NT3)+WCC*VRDC/WC
      IF(LAS.EQ.3) GO TO 33
      YC3=1.-EXD*(14.-EXD*(36.-24.*EXD))
      WCC=WCC/YC2*YC3*RAC/4.
      WWW=WWW/Y2*(Y**4*4.-Y*Y*12.+3.)*RAW/4.
      PW(NT4)=WCC+WWW
       PV(NT1)=PV(NT1)+WCC*VRDC/WC
      GO TO 33
   32 EX1=EXP(-X)
      WCC=-WC*ALF*EX1*WW
      WWW=-WC*(1.-ALF)*EXX*WW
      IF(MEREL.EQ.3) WCC=WCC/RELPOT
      IF(MEREL.EQ.3) WWW=WWW/RELPOT
      PW(NT0)=WCC+WWW
      PV(NT0)=PV(NT0)+WCC*VRDC/WC
      IF(LAS.EQ.0) GO TO 33
      WCC=WCC*RAC
      WWW=WWW*2.*RAW*Y
      PW(NT1)=WCC+WWW
       PV(NT1)=PV(NT1)+WCC*VRDC/WC
      IF(LAS.EQ.1) GO TO 33
      WCC=WCC*RAC/2.
      Y1=2.*Y*Y-1.
      WWW=WWW/Y*Y1*RAW/2.
      PW(NT2)=WCC+WWW
       PV(NT2)=PV(NT2)+WCC*VRDC/WC
      IF(LAS.EQ.2) GO TO 33
      WCC=WCC*RAC/3.
      Y2=Y**3*2.-Y*3.
      WWW=WWW/Y1*Y2*RAW*2./3.
      PW(NT3)=WCC+WWW
       PV(NT3)=PV(NT3)+WCC*VRDC/WC
      IF(LAS.EQ.3) GO TO 33
      WCC=WCC*RAC/4.
      WWW=WWW/Y2*(Y**4*4.-Y*Y*12.+3.)*RAW/4.
      PW(NT4)=WCC+WWW
       PV(NT4)=PV(NT4)+WCC*VRDC/WC
   33 X=(R-RD)/AD
      IF(X.GT.23.) GO TO 34
      EX=EXP(X)
      EX1=1.D0+EX
      EXD=1.D0/EX1
      EXDM=1.D0-EXD
      WDD=-4.D0*WD*EXDM*EXD*WW
      IF(MEREL.EQ.3) WDD=WDD/RELPOT
      PW(NT0)=PW(NT0)+WDD
      PV(NT0)=PV(NT0)+WDD/WD*VD
      IF(LAS.EQ.0) GO TO 35
      Y1=1.D0-2.D0*EXD
      WDD=WDD*Y1*RAD
      PW(NT1)=PW(NT1)+WDD
       PV(NT1)=PV(NT1)+WDD/WD*VD
      IF(LAS.EQ.1) GO TO 35
      Y2=1.D0-6.D0*EXD*EXDM
      WDD=WDD/Y1*Y2*RAD/2.D0
      PW(NT2)=PW(NT2)+WDD
       PV(NT2)=PV(NT2)+WDD/WD*VD
      IF(LAS.EQ.2) GO TO 35
      Y3=1.D0-EXD*(14.D0-EXD*(36.D0-24.D0*EXD))
      WDD=WDD/Y2*Y3*RAD/3.D0
      PW(NT3)=PW(NT3)+WDD
       PV(NT3)=PV(NT3)+WDD/WD*VD
      IF(LAS.EQ.3) GO TO 35
      Y4=1.D0-EXD*(30.D0-EXD*(150.D0-EXD*(240.D0-120.D0*EXD)))
      WDD=WDD/Y3*Y4*RAD/4.D0
      PW(NT4)=PW(NT4)+WDD
       PV(NT4)=PV(NT4)+WDD/WD*VD
      GO TO 35
   34 EX1=EXP(-X)
      WDD=-4.D0*WD*EX1*WW
      IF(MEREL.EQ.3) WDD=WDD/RELPOT
      PW(NT0)=PW(NT0)+WDD
       PV(NT0)=PV(NT0)+WDD/WD*VD
      IF(LAS.EQ.0) GO TO 35
      WDD=WDD*RAD
      PW(NT1)=PW(NT1)+WDD
       PV(NT1)=PV(NT1)+WDD/WD*VD
      IF(LAS.EQ.1) GO TO 35
      WDD=WDD*RAD/2.D0
      PW(NT2)=PW(NT2)+WDD
       PV(NT2)=PV(NT2)+WDD/WD*VD
      IF(LAS.EQ.2) GO TO 35
      WDD=WDD*RAD/3.D0
      PW(NT3)=PW(NT3)+WDD
       PV(NT3)=PV(NT3)+WDD/WD*VD
      IF(LAS.EQ.3) GO TO 35
      WDD=WDD*RAD/4.
      PW(NT4)=PW(NT4)+WDD
       PV(NT4)=PV(NT4)+WDD/WD*VD
   35 IF(I/KIT*KIT.NE.I) GO TO 6
      PR(I1+1)=PV(NT0)
      PI(I1+1)=PW(NT0)
      IF(LAS.EQ.0) GO TO 6
      PR(I1+2)=PV(NT1)
      PI(I1+2)=PW(NT1)
      DO 23 LAC=1,8
      LALALI=LALASI+LASC*(LAC-1)
      LALALR=LALASR+LASC*(LAC-1)
      LAM=LAC-1
      PRC(LALALR+1)=PVC(LALALI+1)
      IF(LAS2.EQ.2) GO TO 23
      PRC(LALALR+2)=PVC(LALALI+2)
   23 CONTINUE
      IF(LAS.EQ.1) GO TO 6
      PR(I1+3)=PV(NT2)
      PI(I1+3)=PW(NT2)
      IF(LAS.EQ.2) GO TO 6
      PR(I1+4)=PV(NT3)
      PI(I1+4)=PW(NT3)
      IF(LAS.EQ.3) GO TO 6
      PR(I1+5)=PV(NT4)
      PI(I1+5)=PW(NT4)
      GO TO 6
   18 DO 9 L=1,LAS2
      
C      IF(NU1.NE.NU2.AND.L.EQ.1) GO TO 9
C      IF(NCA(NU1).NE.NCA(NU2)) GO TO 9

  109 PDC=0.
      LL2=L1+L
      LL2AP=L1AP+L
      IF(MEAPP.EQ.0) GO TO 44
      IF(NU1.EQ.NU2) GO TO 45
      GO TO 42
   45 IF(NU1.NE.1) GO TO 42
   44 DE=1.
      IF(L.EQ.1) GO TO 19
      LAC=2*(L-1)
      IF(LAC.GT.NPD) GO TO 19
      IF(R.LE.RZ) PDC=6.*CONZ/(2.*LAC+1.)*R**LAC/RZ**(LAC+1)*BET(LAC)
      IF(R.GT.RZ) PDC=6.*CONZ/(2.*LAC+1.)*RZ**LAC/R**(LAC+1)*BET(LAC)

C     HIGHER COULOMB MULTIPOLES AS PROPOSES BY SATCHLER

      IF(LAS2.EQ.1) GO TO 19
      IF(R.LE.RZ) FPVC=6.*CONZ/(2.*LAC+1)*(1.-LAC)*R**LAC/RZ**(LAC+1)
      IF(R.GT.RZ) FPVC=6.*CONZ/(2.*LAC+1)*(LAC+2)*RZ**LAC/R**(LAC+1)
      IF(LAC.EQ.2) PDC=PDC+FPVC*(0.180223752*BET(2)+0.241795536*BET(4))
     **BET(2)
      IF(LAC.EQ.4) PDC=PDC+FPVC*(0.241795536*BET(2)+0.163839774*BET(4))
     **BET(2)
      IF(LAC.EQ.6) PDC=PDC+FPVC*0.238565132*BET(2)*BET(4)
C
   19 IF(NPD.NE.0) DE=DEF(1)
      CALL POTET
      IF(NPD.EQ.0) GO TO 13
      VV=VP*POL(L,1)
      WV=WP*POL(L,1)
      N=NAN2
      DE=DEF(N)
      CALL POTET
      VV=VV+VP*POL(L,N)
      WV=WV+WP*POL(L,N)
      V4=0.
      W4=0.
      DO 10 N=2,NAN2,2
      DE=DEF(N)
      CALL POTET
      V4=V4+VP*POL(L,N)
   10 W4=W4+WP*POL(L,N)
      V2=0.
      W2=0.
      DO 11 N=3,NAN1,2
      DE=DEF(N)
      CALL POTET
      V2=V2+VP*POL(L,N)
   11 W2=W2+WP*POL(L,N)

c      write(21,1799) vr,vd,vrdc,wc,wd
c1799 format(6e12.7)
c      pause 1



        

      PV(LL2)=(VV+4.*V4+2.*V2)*ST*WW*4.1887902047864D0
      PW(LL2)=(WV+4.*W4+2.*W2)*ST*WW*4.1887903047864D0
      IF(MEREL.EQ.3) PW(LL2)=PW(LL2)/RELPOT
      IF(L.NE.1) PV(LL2)=PV(LL2)+PDC
      IF(L.NE.1) GO TO 16
      PV(LL2)=PV(LL2)*0.2820947917739D0+PZ
      PW(LL2)=PW(LL2)*0.2820947917739D0
      GO TO 16
   13 PV(LL2)=VP*WW
      IF(L.EQ.1) PV(LL2)=PV(LL2)+PZ
      PW(LL2)=WP*WW
      IF(MEREL.EQ.3) PW(LL2)=PW(LL2)/RELPOT
   16 CONTINUE
      GO TO 43
   42 PV(LL2)=PV(LL2AP)
      PW(LL2)=PW(LL2AP)
   43 IF(I/KIT*KIT.NE.I) GO TO 9
      II2=I1+L
      PR(II2)=PV(LL2)
      PI(II2)=PW(LL2)
    9 CONTINUE
    6 CONTINUE
      KL=M4-4
      IN=M4
      KIT=1
      NIT=5-M4
      IK=NH+KL+2
      STEP1=STEP
    5 CONTINUE
      NH1=NH+KL
      NHI=NH+2

C      write(21,1799) vr,vd,vrdc,wc,wd,enef,WDBW,WISO
C 1799 format(8e12.5)
C      pause 1
      
   38 CONTINUE
      EN=EINT
      RR=RRR
      RZ=RZIN
      
C      write(21,1799) vr,vd,vrdc,wc,wd
C 1799 format(6e12.5)
C      pause 1
      CCOUL=CCCOUL
      CONZ=CONZZ

C           WRITE(21,9997) RK, STEP,CONZ ,3.302222, NH1
C9997  FORMAT(4E12.5,I3)
      CALL RIPAT1
      RETURN
      END
C     ***********************************************************
      SUBROUTINE RIPAT1
C     ***********************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
C     BELOW CARD IS NECESSARY IF MEMORY IS LESS THAN 32Mb
      DOUBLE PRECISION VL,VSL,PV,PW,WSL,WPSL,PL,PSL,PR,PI,PVC,PRC
      
      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM3.FOR'
      INCLUDE 'PRIVCOM6.FOR' 
      INCLUDE 'PRIVCOM8.FOR'
      INCLUDE 'PRIVCOM10.FOR'
      INCLUDE 'PRIVCOM13.FOR' 
      INCLUDE 'PRIVCOM14.FOR'  
         

      CONZZ=CONZ 

      APAN=ANEU/1.00866491600D0 

      CCDE=0.0
      CCCOUL=CCOUL
      CCOUL=CCOUL/2.0
      IF(MECUL.GE.2) CCDE=CCOUL
      IF(MECUL.GE.2) CCOUL=0.0
      CDE=CCDE*MECHA*ZNUC/AT**(1.D0/3.D0)
      CDE12=CCDE*ZNUC/AT**(1.D0/3.D0)


      VRDC=0.0
      VD=0.0
      VDISP=0.0
      CARBUN=931.49378D0
      AMI=939.56536D0
      IF(MECHA.EQ.1) AMI=938.272029D0
      
       IF(MECHA.EQ.1.AND.IDINT(ASP+0.1D0).EQ.1) AMI=1875.612859D0 
       IF(MECHA.EQ.1.AND.IDINT(ASP+0.1D0).EQ.1) ANEU=2.013553212712D0 
      
      DAVERN=0.00
      DAVERP=0.00
      EFERM=EFERMN
      EPAVER=EFERM+DAVERN
      IF(MECHA.EQ.1) EFERM=EFERMP 

      
C           THREE CARDS FOR U238 CHARGE, DELETE THEN FOR USUAL CALCUL.
C           TWO CARDS AT THE END OF SUBROUTINE
             RZIN=RZ
                 IF(MERZZ.EQ.0) GO TO 52
                 IF(MERZZ.EQ.1)   RZ=RZ*(1.D00-RZBWC*(EN-EFERM)**PDIS/
     *((EN-EFERM)**PDIS+RZWID**PDIS))
C             RZ=9.100
C             IF(EN.GT.26.) RZ=9.100-0.31324222*(EN-26.)
C             IF(EN.GT.26.) RZ=RZ-(RZ-6.28082)/9.*(EN-26.)
C             IF(EN.GT.35.) RZ=6.28082
   52        RRR=RR
      APAN=ANEU/1.00866491600D0 
C
C           RR ENERGY DEPENDENCE
C
C     IF(MERRR.EQ.1)   RR=RR*(1.D00-RRBWC*EN**2/(EN**2+RRWID**2))
      IF(MERRR.EQ.1)   RR=RR*(1.D00-RRBWC*(EN-EFERM)**PDIS/
     *((EN-EFERM)**PDIS+RRWID**PDIS))
      
     

      

      REL=(EN+AMI)/AMI
      IF(MEREL.EQ.0) REL=1.D+0
      IF(REL.NE.1.) RELPOT=2.D+0*(EN+AMI)/(2.D+0*AMI+EN)

      
      EIRELM=SQRT(REL**2*AT**2+2.D+0*REL*ANEU*AT+ANEU**2)
      EIRELT=AT*SQRT(REL**2*ANEU**2+2.D+0*REL*ANEU*AT+AT**2)
      DEREL=SQRT(2.D+0*REL*ANEU*AT+AT**2+ANEU**2)
      WW=4.8257984E-2*EIRELM*EIRELT/(ANEU*EIRELM+EIRELT)/DEREL*APAN
      IF(MEREL.EQ.1.OR.MEREL.EQ.3) WW=WW*RELPOT

            
      





      LAS1=LAS+1
      LAS2=(LAS+2)/2
      IF(MEPOT.EQ.1) GO TO 39
      LAS2=LAS+1
   39 CONTINUE

      
   41 PIM=1.996852
      
      MNULA=300*NUR*LAS2
      LALAS=9
      LASC=1
      IF(LAS2.GE.3) LALAS=18
      IF(LAS2.GE.3) LASC=2
      MLA=300*LAS2
      MNU=300*NUR



   
      EINT=EN
      
      DO 38 NU1=1,NUR
      MNULA1=MNULA*(NU1-1)
      MNU1=MNU*(NU1-1)
      DO 38 NU2=1,NUR
      MNULA2=MNULA1+MLA*(NU2-1)
      MNU2=MNU1+300*(NU2-1)
      EN=EINT-(EL(NU1)+EL(NU2))/2
C       if(nu1.ne.nu2)EN=EINT-ABS(EL(NU1)-EL(NU2))/2.
        if(nu1.ne.nu2)EN=EINT-(EL(NU1)+EL(NU2))/2.
      IF(MEAPP.EQ.1) EN=EINT
      


C     COUPLINGIN BETWEEN LEVELS OF THE DIFFERENT PARTITIONS    !!!!!!!!
      
      
      IF(NCA(NU1).EQ.NCA(NU2)) GO TO 38
      IF(MEHAM.EQ.1) GO TO 40
      IF(NTU(NU1).NE.NTU(NU2).OR.NNB(NU1).NE.NNB(NU2)
     *.OR.NNG(NU1).NE.NNG(NU2).OR.NNO(NU1).NE.NNO(NU2))GO TO 38

      
   40 EFENU12=(EFERMN+EFERMP)/2.D0
      
CCHA      ENEF=EN-EFENU12-CDE
      
CCHA      EN=EN-CDE

      ENEF=EN-EFENU12-1./2.*CDE             !!!!!!     ALD CASE  withOTHER RIGHT     !!!!!!!
      ENEF=EN-EFENU12 
      
      EN=EN     
      
      COENH=1.D0
C     
C     IF(MECUL.GE.2)    COENH=1.D0+ALAVR*VRLA*EXP(-ALAVR*ENEF)
      

  
  
      VISO=CISO*2.D0*SQRT(AT-2.D0*ZNUC)/AT
      WISO=WCISO*2.D0*SQRT(AT-2.D0*ZNUC)/AT
      WVISO=WDISO*2.D0*SQRT(AT-2.D0*ZNUC)/AT
      VRLANP=VISO*COENH
      WCBWNP=WVISO
      WDBWNP=WISO


C      WCP=Wv(Ef-Ea)=Wv(Ef+Ea)

      WCP=(WCBWNP)*(EA)**PDIS/
     *(EA**PDIS+WCWID**PDIS)

C      write(21,1799) vcp,wcbwnp

C      pause 89

      
C     TWO OPTIONS FOR COULOMB CORRECTIONS
C


      VRDIR=(VR1+2.*VR2*ENEF+3.*VR3*ENEF**2
     *-VRLA*ALAVR*EXP(-ALAVR*ENEF))*VISO/(VR0+VRLA)   
     **(-1.D0) 
      

      IF(MEDIS.EQ.0) GO TO 394
C    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     VHF DERIVATIVE AT EN
      EEE=EN
      EN=EN-0.005D+00
      ENCON=EN

      CALL VHFROM
      VHFM=VHF

      EN=EN+0.01D+00
      ENCON=EN

      CALL VHFROM
      VHFP=VHF

      VRDIR=(VHFM-VHFP)*100
      EN=EEE



C     DISPERSIVE CONTRIBUTION TO VR

  983 VDISP=DOM_INT_WV(Eferm,EPAVER,WCBWNP,WCWID,INT(PDIS),EN,VDISD)





C     NO ISOSCALOR PART OF Wv
        
      IF(EA.GT.0.d0) THEN 
       T12P=WCP*DOM_INT_T1(Eferm,EA,EN+0.01)
       T12M=WCP*DOM_INT_T1(Eferm,EA,EN-0.01)
       T12D=(T12M-T12P)*50.
      ENDIF
      IF(MECUL.GE.2.OR.MECHA.EQ.0) T12D=0.D0



      
      VRDIRC=VDISD+T12D


   
      
      VCULC=MECHA*ZNUC/AT**0.333333333*CCOUL*VRDIRC
      WDUL=0.
      WCUL=0.




C     DISPERSIVE WD CONTRIBUTION TO VR

      

  984        IF(ALAWD.NE.0.D0) THEN    
        VDP=DOM_INT_WS
     * (Eferm,EPAVER,WDBWNP,WDWID,ALAWD,INT(PDIS),EN,VDD)
      ELSE
        VDP=DOM_INT_WV(Eferm,EPAVER,WDBWNP,WDWID,INT(PDIS),EN,VDD)
      ENDIF

      IF(EN.GT.WDSHI.AND.ALAWD.EQ.0.D+00) THEN 
        VDP1=VDP
        VDD1=VDD
        VDP=DOM_INT_WV(Eferm,WDSHI,WDBWNP,WDWID2,INT(PDIS),EN,VDD)
        VDP=VDP1-VDP
        VDD=VDD1-VDD 
      ENDIF      


       


       VDCUL=VDD*MECHA*ZNUC/AT**0.333333333*CCOUL
C
       PDIS1=PDIS-1
      
C
      WDUL=0.
C
      WCUL=0.

  394 VCUL=MECHA*ZNUC/AT**0.333333333*CCOUL*VRDIR
      IF(MECUL.EQ.1) VCUL=MECHA*ZNUC/AT**0.333333333*CCOUL
      
C   

      VR=(VR0+VR1*ENEF+VR2*ENEF**2+VR3*ENEF**3
     *+VRLA*EXP(-ALAVR*ENEF))*VISO/(VR0+VRLA)+VCUL

      IF(MEDIS.EQ.0) GO TO 395
      
C    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     VHF AND VR AT EN
      
      ENCON=EN

      CALL VHFROM
      
      
      VR=VHF+VCUL               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      VRDC=VDISP+VCULC
      
      
C         ENEF=ENCON-EFERM

C    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






       VD=VDP+VDCUL

  395 CONTINUE

      ENE=ENEF                       !!!!!!!!!!!!!!!!!!!!!!!!

CCHA      IF(MECUL.EQ.3) ENE=ENE-CDE   
      IF(MECUL.EQ.3) ENE=ENE+CDE/2.D0

     
      ENESHID=ENE-WDSHI+EFERM
      
      IF(ENE.LT.BNDC) GO TO 900
      WD=WD0+WD1*BNDC+(ENE-BNDC)*WDA1+(WDBWNP)*EXP(-ALAWD*ENE)
     **ENE**PDIS/(ENE**PDIS+WDWID**PDIS)+WDUL
      IF(ENE+EFERM.GT.WDSHI.AND.ALAWD.EQ.0.D+00)WD=WD-(WDBWNP)
     **ENESHID**PDIS/(ENESHID**PDIS+WDWID2**PDIS)
      WC=WC0+WC1*BNDC+(ENE-BNDC)*WCA1+(WCBWNP)*ENE**PDIS/
     *(ENE**PDIS+WCWID**PDIS)+WCUL
      GO TO 901
  900 WD=WD0+WD1*ENE+(WDBWNP)*EXP(-ALAWD*ENE)
     **ENE**PDIS/(ENE**PDIS+WDWID**PDIS)+WDUL
      IF(ENE+EFERM.GT.WDSHI.AND.ALAWD.EQ.0.D+00)WD=WD-(WDBWNP)
     **ENESHID**PDIS/(ENESHID**PDIS+WDWID2**PDIS)

      WC=WC0+WC1*ENE+(WCBWNP)*ENE**PDIS/
     *(ENE**PDIS+WCWID**PDIS)+WCUL





  901 CONTINUE



C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  972 REL=(EN+AMI)/AMI
      IF(REL.LT.1.D00) REL=1.D00
      EIRELM=SQRT(REL**2*AT**2+2.D+0*REL*ANEU*AT+ANEU**2)
      EIRELT=AT*SQRT(REL**2*ANEU**2+2.D+0*REL*ANEU*AT+AT**2)
      DEREL=SQRT(2.D+0*REL*ANEU*AT+AT**2+ANEU**2)
C        ENCON=(REL**2-1.D0)*AMI*AT*(ANEU+AT)/DEREL**2/2.D+0
      WW4=1.2064496D-2*EIRELM*EIRELT/(ANEU*EIRELM+EIRELT)/DEREL*APAN

C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       WW4=1.
c       VSO=VS*EXP(-WW4*ALASO*(ENEF+CDE))

      IF(MEDIS.EQ.0) GO TO 396
      IF(MEDIS.EQ.2) GO TO 397


      

C
C     NO ISOSCALAR CONTRIBUTION FROM NON SYMMETRIC IMAGINARY POSITIVE POTENTIAL Wv 
C
 397  IF(EA.GT.0.D0) VRDC=WCP*DOM_INT_T1(Eferm,EA,EN)+VRDC

  396 CONTINUE
      
      
C     WD=WD+WISO
      IF(WD.LE.0.)WD=0.
      IF(WC.LE.0.)WC=0.
      AR=AR0+AR1*ENEF
      AW=AW0+AW1*ENEF
      AC=AC0+AC1*ENEF
      AD=AD0+AD1*ENEF
      IF (ENEF.GT.BNDC) AD=AD0+AD1*BNDC
      AS=AS0+AS1*ENEF
      RAR=RR/AR
      RAC=RC/AC
      RAW=RW/AW
      RAD=RD/AD
      M=8
      M4=4*M
      KL=0
      IN=1
      KIT=M
      NIT=1
      IK=M4-1

      
C
      STEP=RK/NH
      STEP1=STEP/M
      
C
      DO 5 K=1,2
      DO 6 I=IN,IK
      


      



      II=MNU2+I
      IR=I/KIT+NIT
      I1=(IR-1)*LAS2+MNULA2
      LALASR=LALAS*(IR-1)
      IIC=MNU2+IR
      LALASI=LALAS*(I-1)
      L1=(I-1)*LAS2+MNULA2
      L1AP=(I-1)*LAS2
      IF(K.EQ.2.AND.I.EQ.IK) GO TO 14
      R=STEP1*(I-KL)
      GO TO 15
   14 R=0.D0
 
C     NON-DIAGONAL TERMS BY ISOSPIN

   15 CONTINUE
      
C      write(21,9999)vr,vd,vrdc,wc,wd,ww, 7.77777777
C 9999 format(11e12.5)
      IF(MEPOT.EQ.1) GO TO 17 
      PZ=0.D0
   21 NT0=L1+1
      NT1=L1+2
      NT2=L1+3
      NT3=L1+4
      NT4=L1+5
      X=(R-RR)/AR
      IF(X.GT.23.) GO TO 30
      EX=EXP(X)
      EX1=1.D0+EX
      EXD=1.D0/EX1
      EXDM=1.D0-EXD
      IF(JO(NU1).EQ.JO(NU2)) PV(NT0)=-VR*EXD*WW
      IF(LAS.EQ.0) GO TO 31
      PV(NT1)=(-VR*EXD*WW-PZ)*EXDM*RAR
      IF(LAS.EQ.1) GO TO 31
      Y1=1.D0-2.D0*EXD
      PV(NT2)=PV(NT1)*Y1*RAR/2.D0
      IF(LAS.EQ.2) GO TO 31
      Y2=1.D0-6.D0*EXD*EXDM
      PV(NT3)=PV(NT2)/Y1*Y2*RAR/3.D0
      IF(LAS.EQ.3) GO TO 31
      Y3=1.D0-EXD*(14.D0-EXD*(36.D0-24.D0*EXD))
      PV(NT4)=PV(NT3)/Y2*Y3*RAR/4.D0
      GO TO 31
   30 EX1=EXP(-X)
      IF(JO(NU1).EQ.JO(NU2)) PV(NT0)=-VR*EX1*WW
      IF(LAS.EQ.0) GO TO 31
      PV(NT1)=(-VR*EX1*WW-PZ)*RAR
      IF(LAS.EQ.1) GO TO 31
      PV(NT2)=PV(NT1)*RAR/2.D0
      IF(LAS.EQ.2) GO TO 31
      PV(NT3)=PV(NT2)*RAR/3.D0
      IF(LAS.EQ.3) GO TO 31
      PV(NT4)=PV(NT3)*RAR/4.D0
   31 X=(R-RC)/AC
      Y=(R-RW)/AW
      XX=Y*Y
      IF(XX.GT.180.) EXX=0.
      IF(XX.LE.180.) EXX=EXP(-XX)
      IF(X.GT.23.) GO TO 32
      EX=EXP(X)
      EX1=1.D0+EX
      EXD=1.D0/EX1
      EXDM=1.D0-EXD
      WCC=-WC*ALF*EXD*WW
      WWW=-WC*(1.-ALF)*EXX*WW
      IF(MEREL.EQ.3) WCC=WCC/RELPOT
      IF(MEREL.EQ.3) WWW=WWW/RELPOT
      IF(JO(NU1).EQ.JO(NU2)) PW(NT0)=WCC+WWW
      IF(JO(NU1).EQ.JO(NU2).AND.WC.NE.0.D0) PV(NT0)=PV(NT0)+WCC*VRDC/WC
      IF(LAS.EQ.0) GO TO 33
      WCC=WCC*EXDM*RAC
      WWW=WWW*2.D0*RAW*Y
      PW(NT1)=WCC+WWW
      IF(WC.NE.0.D0) PV(NT1)=PV(NT1)+WCC*VRDC/WC
      IF(LAS.EQ.1) GO TO 33
      YC1=1.D0-2.D0*EXD
      WCC=WCC*YC1*RAC/2.D0
      Y1=2.D0*Y*Y-1.D0
      WWW=WWW/Y*Y1*RAW/2.D0
      PW(NT2)=WCC+WWW
      IF(WC.NE.0.D0) PV(NT2)=PV(NT2)+WCC*VRDC/WC
      IF(LAS.EQ.2) GO TO 33
      YC2=1.D0-6.D0*EXD*EXDM
      WCC=WCC/YC1*YC2*RAC/3.D0
      Y2=Y**3*2.D0-Y*3.D0
      WWW=WWW/Y1*Y2*RAW*2.D0/3.D0
      PW(NT3)=WCC+WWW
      IF(WC.NE.0.D0) PV(NT3)=PV(NT3)+WCC*VRDC/WC
      IF(LAS.EQ.3) GO TO 33
      YC3=1.D0-EXD*(14.D0-EXD*(36.D0-24.D0*EXD))
      WCC=WCC/YC2*YC3*RAC/4.
      WWW=WWW/Y2*(Y**4*4.-Y*Y*12.D0+3.D0)*RAW/4.D0
      PW(NT4)=WCC+WWW
      IF(WC.NE.0.D0) PV(NT4)=PV(NT4)+WCC*VRDC/WC
   32 EX1=EXP(-X)
      WCC=-WC*ALF*EX1*WW
      WWW=-WC*(1.D0-ALF)*EXX*WW
      IF(MEREL.EQ.3) WCC=WCC/RELPOT
      IF(MEREL.EQ.3) WWW=WWW/RELPOT
      IF(JO(NU1).EQ.JO(NU2)) PW(NT0)=WCC+WWW
      IF(JO(NU1).EQ.JO(NU2).AND.WC.NE.0.D0) PV(NT0)=PV(NT0)+WCC*VRDC/WC
      IF(LAS.EQ.0) GO TO 33
      WCC=WCC*RAC
      WWW=WWW*2.D0*RAW*Y
      PW(NT1)=WCC+WWW
      IF(WC.NE.0.D0) PV(NT1)=PV(NT1)+WCC*VRDC/WC
      IF(LAS.EQ.1) GO TO 33
      WCC=WCC*RAC/2.D0
      Y1=2.D0*Y*Y-1.D0
      WWW=WWW/Y*Y1*RAW/2.D0
      PW(NT2)=WCC+WWW
      IF(WC.NE.0.D0) PV(NT2)=PV(NT2)+WCC*VRDC/WC
      IF(LAS.EQ.2) GO TO 33
      WCC=WCC*RAC/3.
      Y2=Y**3*2.D0-Y*3.D0
      WWW=WWW/Y1*Y2*RAW*2.D0/3.D0
      PW(NT3)=WCC+WWW
      IF(WC.NE.0.D0) PV(NT3)=PV(NT3)+WCC*VRDC/WC
      IF(LAS.EQ.3) GO TO 33
      WCC=WCC*RAC/4.D0
      WWW=WWW/Y2*(Y**4*4.D0-Y*Y*12.D0+3.D0)*RAW/4.D0
      PW(NT4)=WCC+WWW
      IF(WC.NE.0.D0) PV(NT4)=PV(NT4)+WCC*VRDC/WC
   33 X=(R-RD)/AD
      IF(X.GT.23.) GO TO 34
      EX=EXP(X)
      EX1=1.D0+EX
      EXD=1.D0/EX1
      EXDM=1.D0-EXD
      WDD=-4.D0*WD*EXDM*EXD*WW
      IF(MEREL.EQ.3) WDD=WDD/RELPOT
      IF(JO(NU1).EQ.JO(NU2)) PW(NT0)=PW(NT0)+WDD
      IF(JO(NU1).EQ.JO(NU2).AND.WD.NE.0.D0) PV(NT0)=PV(NT0)+WDD/WD*VD
      IF(LAS.EQ.0) GO TO 35
      Y1=1.D0-2.D0*EXD
      WDD=WDD*Y1*RAD
      PW(NT1)=PW(NT1)+WDD
      IF(WD.NE.0.D0) PV(NT1)=PV(NT1)+WDD/WD*VD
      IF(LAS.EQ.1) GO TO 35
      Y2=1.D0-6.D0*EXD*EXDM
      WDD=WDD/Y1*Y2*RAD/2.D0
      PW(NT2)=PW(NT2)+WDD
      IF(WD.NE.0.D0) PV(NT2)=PV(NT2)+WDD/WD*VD
      IF(LAS.EQ.2) GO TO 35
      Y3=1.D0-EXD*(14.D0-EXD*(36.D0-24.D0*EXD))
      WDD=WDD/Y2*Y3*RAD/3.D0
      PW(NT3)=PW(NT3)+WDD
      IF(WD.NE.0.D0) PV(NT3)=PV(NT3)+WDD/WD*VD
      IF(LAS.EQ.3) GO TO 35
      Y4=1.-EXD*(30.D0-EXD*(150.D0-EXD*(240.D0-120.D0*EXD)))
      WDD=WDD/Y3*Y4*RAD/4.D0
      PW(NT4)=PW(NT4)+WDD
      IF(WD.NE.0.D0) PV(NT4)=PV(NT4)+WDD/WD*VD
      GO TO 35
   34 EX1=EXP(-X)
      WDD=-4.D0*WD*EX1*WW
      IF(MEREL.EQ.3) WDD=WDD/RELPOT
      IF(JO(NU1).EQ.JO(NU2)) PW(NT0)=PW(NT0)+WDD
      IF(JO(NU1).EQ.JO(NU2).AND.WD.NE.0.D0) PV(NT0)=PV(NT0)+WDD/WD*VD
      IF(LAS.EQ.0) GO TO 35
      WDD=WDD*RAD
      PW(NT1)=PW(NT1)+WDD
      IF(WD.NE.0.D0) PV(NT1)=PV(NT1)+WDD/WD*VD
      IF(LAS.EQ.1) GO TO 35
      WDD=WDD*RAD/2.D0
      PW(NT2)=PW(NT2)+WDD
      IF(WD.NE.0.D0) PV(NT2)=PV(NT2)+WDD/WD*VD
      IF(LAS.EQ.2) GO TO 35
      WDD=WDD*RAD/3.D0
      PW(NT3)=PW(NT3)+WDD
      IF(WD.NE.0.D0) PV(NT2)=PV(NT2)+WDD/WD*VD
      IF(LAS.EQ.3) GO TO 35
      WDD=WDD*RAD/4.D0
      PW(NT4)=PW(NT4)+WDD
      IF(WD.NE.0.D0) PV(NT4)=PV(NT4)+WDD/WD*VD
   35 IF(I/KIT*KIT.NE.I) GO TO 6
      PR(I1+1)=PV(NT0)
      PI(I1+1)=PW(NT0)

C      print 876, i1,nt0,Nu1,Nu2, JO(NU1),JO(NU2),PR(I1+1),PI(I1+1)
C  876 FORMAT (6i7,6e12.3)
C      pause 22

      IF(LAS.EQ.0) GO TO 6
      PR(I1+2)=PV(NT1)
      PI(I1+2)=PW(NT1)
      IF(LAS.EQ.1) GO TO 6
      PR(I1+3)=PV(NT2)
      PI(I1+3)=PW(NT2)
      IF(LAS.EQ.2) GO TO 6
      PR(I1+4)=PV(NT3)
      PI(I1+4)=PW(NT3)
      IF(LAS.EQ.3) GO TO 6
      PR(I1+5)=PV(NT4)
      PI(I1+5)=PW(NT4)
      GO TO 6



      
   17 DO 99 L=1,LAS2


      
      
      PDC=0.D0
      PZ=0.D0
      LL2=L1+L
      LL2AP=L1AP+L
      
      DE=DEF(1)

      
      CALL POTET
      
      VV=VP*POL(L,1)
      WV=WP*POL(L,1)
      N=NAN2
      DE=DEF(N)
      CALL POTET
      VV=VV+VP*POL(L,N)
      WV=WV+WP*POL(L,N)
      V4=0.
      W4=0.
      DO 10 N=2,NAN2,2
      DE=DEF(N)
      CALL POTET
      V4=V4+VP*POL(L,N)
   10 W4=W4+WP*POL(L,N)
      V2=0.
      W2=0.
      DO 11 N=3,NAN1,2
      DE=DEF(N)
      CALL POTET
      V2=V2+VP*POL(L,N)
   11 W2=W2+WP*POL(L,N)

      
      IF(JO(NU1).EQ.JO(NU2).AND.L.EQ.1) GO TO 780
      IF(L.EQ.1) GO TO 99
      
        
      
C      GO TO 99                                           !!!!!!!CHECK!!!!!!!


  780 PV(LL2)=(VV+4.D0*V4+2.D0*V2)*ST*WW*4.1887902047864D0
      PW(LL2)=(WV+4.D0*W4+2.D0*W2)*ST*WW*4.1887903047864D0
      IF(MEREL.EQ.3) PW(LL2)=PW(LL2)/RELPOT

C      write(21,1234) i,l,r
C 1234 format(2i5,e12.5)
C     write(21,9999)vr,vd,vrdc,wc,wd,ww, 7.77777777
C      write(21,9999)v2,v4,vv, 7.77777777,w2,w4,wv, pv(LL2),pWW(LL2), 5.5

      IF(K.EQ.2.AND.I.EQ.IK) GO TO 114
      R=STEP1*(I-KL)
      GO TO 115
  114 R=0.D0
      IR=1

      GO TO 116

  115 IF(I/KIT*KIT.NE.I) GO TO 99
      IR=I/KIT+NIT
  116 I1=(IR-1)*LAS2+MNULA2
      II2=I1+L
      PR(II2)=PV(LL2)
      PI(II2)=PW(LL2)
      
            

   99 CONTINUE

    6 CONTINUE
      KL=M4-4
      IN=M4
      KIT=1
      NIT=5-M4
      IK=NH+KL+2
      STEP1=STEP
    5 CONTINUE
      NH1=NH+KL
      NHI=NH+2
   38 CONTINUE
      EN=EINT
      RR=RRR
      RZ=RZIN
      CCOUL=CCCOUL
      CONZ=CONZZ


C      write(21,1799) vr,vd,vrdc,wc,wd, enef,WDBWNP
C 1799 format(7e12.5)
C     pause 2

C     WRITE(21,9997) RK, STEP,CONZ ,2.2222222, NH1
C9997  FORMAT(4E12.5, I3)
      RETURN
      END
C     ***********************************************************
      SUBROUTINE CMATC
C     ***********************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM10.FOR'
       INCLUDE 'PRIVCOM13.FOR'              
      MSOLIN=MESOL
C      IF(NSS.LT.5) MESOL=2
      IF(MESOL.LE.3) GO TO 5
      NCLLT=NCLL
      INCRT=INCR
      IF(MESOL.LT.INCR) INCR=MESOL
      IF(MESOL.LT.NCLL) NCLL=MESOL
      CALL KNDIT
      CALL SOSIT
      CALL MASCT
      NCLL=NCLLT
      INCR=INCRT  
      IF(MESOL.GE.NCLLT .AND. MEPRI.LT.98) PRINT 22,NSS,NCLL
      IF(MESOL.GE.NCLLT) GO TO 4
      CALL KNDIT
      GO TO 3
    5 CALL KNDIT
      IF(MESOL-2) 1,2,3
    1 IF(NSS.LT.5) GO TO 2
      IEC=INCC*(BET(2)*40.+3)
      IF(NCLL.GT.IEC) GO TO 3
    2 CALL SOSIT
      MESOL=2
      CALL MASCT
      IF(MEPRI.LT.98) PRINT 22,NSS,NCLL
  22  FORMAT (2I5,E20.7)
      GO TO 4
    3 CALL ECISS
      IF(NECI.EQ.0) GO TO 2
    4 MESOL=MSOLIN
      RETURN
      END
C*******************************************************************************
      SUBROUTINE RCWF(RHO,ETA,MINL,MAXL,FC,FCP,GC,GCP,ACCUR,STEP,NUMBR) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      DOUBLE PRECISION K,K1,K2,K3,K4,M1,M2,M3,M4
      DIMENSION FC(1),FCP(1),GC(1),GCP(1)                               
C *** COULOMB WAVEFUNCTIONS CALCULATED AT R = RHO BY THE                
C *** CONTINUED-FRACTION METHOD OF STEED   MINL,MAXL ARE ACTUAL L-VALUES
C *** SEE BARNETT FENG STEED AND GOLDFARB COMPUTER PHYSICS COMMUN 1974  
C *** RCWF            - A MODIFICATION OF THE REAL COULOMB              
C *** WAVEFUNCTION PROGRAM RCWFN.  A.R. BARNETT.                        
C *** REF. IN COMP. PHYS. COMMUN. 11 (1976) 141                         
C  INSERT AFTER CARD 48 THE NEXT SIX CARDS                              
C *** CPC 8 (1974) 377-395                                              
C *** IF NUMBR = 1 FC(MAXL) TO FC(MINL) ARE RETURNED  ARB MARCH 1975    
C *** IF NUMBR = 2 FC AND GC ARRAYS ARE RETURNED                        
C *** IF NUMBR = 3 FC,GC,FCP,GCP ARRAYS ARE RETURNED                    
      NUM = 1                                                           
      IF(NUMBR.GE.1.AND.NUMBR.LE.3) NUM = NUMBR                         
      PACE = STEP                                                       
      ACC  = ACCUR                                                      
      IF(PACE.LT.100.0) PACE = 100.0                                    
      IF(ACC.LT.1.0E-15.OR.ACC.GT.1.0E-6) ACC = 1.0E-6                  
      R    = RHO                                                        
      KTR  = 1                                                          
      LMAX = MAXL                                                       
      LMIN1= MINL + 1                                                   
      XLL1 = FLOAT(MINL*LMIN1)                                          
      ETA2 = ETA*ETA                                                    
      TURN = ETA + SQRT(ETA2 + XLL1)                                    
      IF(R.LT.TURN.AND.ABS(ETA).GE.1.0E-6) KTR = -1                     
      KTRP = KTR                                                        
      GO TO 2                                                           
1     R    = TURN                                                       
      TF   = F                                                          
      TFP  = FP                                                         
      LMAX = MINL                                                       
      KTRP = 1                                                          
2     ETAR = ETA*R                                                      
      RHO2 =   R*R                                                      
      PL   = FLOAT(LMAX + 1)                                            
      PMX  = PL + 0.5                                                   
C *** CONTINUED FRACTION FOR FP(MAXL)/F(MAXL)  XL IS F  XLPRIME IS FP **
      FP  = ETA/PL + PL/R                                               
      DK  = ETAR*2.0                                                    
      DEL = 0.0                                                         
      D   = 0.0                                                         
      F   = 1.0                                                         
      K   = (PL*PL - PL + ETAR)*(2.0*PL - 1.0)                          
      IF(PL*PL+PL+ETAR.NE.0.0) GO TO 3                                  
      R   = R + 1.0E-6                                                  
      GO TO 2                                                           
3     H   = (PL*PL + ETA2)*(1.0 - PL*PL)*RHO2                           
      K   = K + DK + PL*PL*6.0                                          
      D   =  1.0/(D*H + K)                                              
      DEL =  DEL*(D*K - 1.0)                                            
      IF(PL.LT.PMX) DEL = -R*(PL*PL + ETA2)*(PL + 1.0)*D/PL             
      PL  = PL + 1.0                                                    
      FP  = FP + DEL                                                    
      IF(D.LT.0.0) F = -F                                               
      IF(PL.GT.20000.) GO TO 11                                         
CC      IF(ABS(DEL/FP).GE.ACC) GO TO 3                                 
C  REPLACE CARD 91 BY THE NEXT CARD                                     
      IF(ABS(DEL).GE.ABS(FP)*ACC) GO TO 3                               
      FP  = F*FP                                                        
      IF( LMAX.EQ.MINL) GO TO 5                                        
      FC (LMAX+1) = F                                                   
CC      FCP(LMAX+1) = FP                                              
C  REPLACE CARD 95 BY NEXT TWO CARDS                                   
      FPL          = FP                                                 
      IF(NUM.EQ.3) FCP(LMAX+1) = FP                                     
C *** DOWNWARD RECURSION TO MINL FOR F AND FP, ARRAYS GC,GCP ARE STORAGE
      L  = LMAX                                                         
C  INSERT NEXT CARD AFTER 97                                            
      RI = 1.0/R                                                        
      DO 4 LP  = LMIN1,LMAX                                             
      PL = FLOAT(L)                                                     
CC      GC (L+1) = ETA/PL + PL/R                                       
CC      GCP(L+1) = SQRT(ETA2 + PL*PL)/PL                               
CC      FC (L)   = (GC(L+1)*FC(L+1) + FCP(L+1))/GCP(L+1)              
CC      FCP(L)   =  GC(L+1)*FC(L)   - GCP(L+1)*FC(L+1)                 
C  REPLACE CARDS 100 - 103 BY NEXT 10 CARDS                            
      SL = ETA/PL + PL*RI                                              
      RL = SQRT(1.D0 + ETA2/(PL*PL))                                   
      FC(L) = (SL*FC(L+1) + FPL)/RL                                  
      FPL1  =  SL*FC(L) - RL*FC(L+1)                                 
      FPL   = FPL1                                                  
      IF(NUM.EQ.1) GO TO 4                                          
      GC (L+1) = RL                                                 
      IF(NUM.LE.2) GO TO 4                                           
      GCP(L+1) = SL                                                  
      FCP(L)   = FPL1                                               
4     L  = L - 1                                                    
      F  = FC (LMIN1)                                                
CC      FP = FCP(LMIN1)                                             
C  REPLACE CARD 106 BY THE NEXT CARD                                 
      FP = FPL1                                                       
5     IF(KTRP.EQ.-1) GO TO 1                                         
C *** REPEAT FOR R = TURN IF RHO LT TURN                              
C *** NOW OBTAIN P + I.Q FOR MINL FROM CONTINUED FRACTION (32)         
C *** REAL ARITHMETIC TO FACILITATE CONVERSION TO IBM USING DOUBLE PRECISION
      P  = 0.d0                                                          
      Q  = R - ETA                                                      
      PL = 0.d0                                                          
      AR = -(ETA2 + XLL1)                                               
      AI =   ETA                                                        
      BR = 2.d0*Q                                                       
      BI = 2.d0                                                          
      WI = 2.d0*ETA                                                      
      DR =   BR/(BR*BR + BI*BI)                                         
      DI =  -BI/(BR*BR + BI*BI)                                        
      DP = -(AR*DI + AI*DR)                                            
      DQ =  (AR*DR - AI*DI)                                            
6     P  =  P + DP                                                     
      Q  =  Q + DQ                                                     
      PL = PL + 2.d0                                                     
      AR = AR + PL                                                     
      AI = AI + WI                                                     
      BI = BI + 2.d0                                                    
      D  = AR*DR - AI*DI + BR                                           
      DI = AI*DR + AR*DI + BI                                           
      T  = 1.d0/(D*D + DI*DI)                                           
      DR =  T*D                                                        
      DI = -T*DI                                                        
      H  = BR*DR - BI*DI - 1.d0                                          
      K  = BI*DR + BR*DI                                               
      T  = DP*H  - DQ*K                                                 
      DQ = DP*K  + DQ*H                                                
      DP = T                                                           
      IF(PL.GT.46000.d0) GO TO 11                                         
      IF(ABS(DP)+ABS(DQ).GE.(ABS(P)+ABS(Q))*ACC) GO TO 6                
      P  = P/R                                                          
      Q  = Q/R                                                          
C *** SOLVE FOR FP,G,GP AND NORMALISE F  AT L=MINL                      
      G  = (FP - P*F)/Q                                                
      GP = P*G - Q*F                                                   
      W  = 1.D0/SQRT(FP*G - F*GP)                                       
      G  = W*G                                                          
      GP = W*GP                                                        
      IF(KTR.EQ.1) GO TO 8                                             
      F  = TF                                                          
      FP = TFP                                                          
      LMAX = MAXL                                                       
C *** RUNGE-KUTTA INTEGRATION OF G(MINL) AND GP(MINL) INWARDS FROM TURN 
C ***             SEE FOX AND MAYERS 1968 PG 202                        
      IF(RHO.LT.0.2*TURN) PACE = 999.d0                                  
      R3 = 1.d0/3.0D0                                                    
      H  = (RHO - TURN)/(PACE + 1.d0)                                    
CC      H2 = 0.5*H                                                      
C  REPLACE CARD 158 BY THE NEXT CARD                                   
      H2 = 0.5D0*H                                                      
C      I2 = IFIX(PACE + 0.001)                                          
      I2 = PACE + 0.001d0                                                 
      ETAH = ETA*H                                                      
      H2LL = H2*XLL1                                                    
      S  = (ETAH + H2LL/R  )/R   - H2                                   
7     RH2= R + H2                                                       
      T  = (ETAH + H2LL/RH2)/RH2 - H2                                   
      K1 = H2*GP                                                        
      M1 =  S*G                                                         
      K2 = H2*(GP + M1)                                                 
      M2 =  T*(G  + K1)                                                 
      K3 =  H*(GP + M2)                                                 
      M3 =  T*(G  + K2)                                                 
      M3 =     M3 + M3                                                  
      K4 = H2*(GP + M3)                                                 
      RH = R + H                                                        
      S  = (ETAH + H2LL/RH )/RH  - H2                                   
      M4 =  S*(G + K3)                                                  
      G  = G  + (K1 + K2 + K2 + K3 + K4)*R3                             
      GP = GP + (M1 + M2 + M2 + M3 + M4)*R3                             
      R  = RH                                                           
      I2 = I2 - 1                                                       
      IF(ABS(GP).GT.1.0D300) GO TO 11                                   
      IF(I2.GE.0) GO TO 7                                               
      W  = 1.d0/(FP*G - F*GP)                                            
C *** UPWARD RECURSION FROM GC(MINL) AND GCP(MINL),STORED VALUES ARE R,S
C *** RENORMALISE FC,FCP FOR EACH L-VALUE                               
CC8     GC (LMIN1) = G                                                
CC      GCP(LMIN1) = GP                                                
CC      IF(LMAX.EQ.MINL) GO TO 10                                       
C  REPLACE CARDS 185 - 187 BY THE NEXT 6 CARDS                         
CC      IF(NUM.GE.2) GC (LMIN1) = G                                    
C  REPLACE CARD A056 BY                                                
8     IF(NUM.GE.2) GC (LMIN1) = G                                      
      IF(NUM.EQ.3) GCP(LMIN1) = GP                                    
      IF(NUM.EQ.3) FCP(LMIN1) = FP*W                                  
      FC(LMIN1) = F*W                                                  
      IF(LMAX.EQ.MINL) RETURN                                          
      GPL = GP                                                        
      DO  9  L = LMIN1,LMAX                                            
CC      T        = GC(L+1)                                             
CC      GC (L+1) = (GC(L)*GC (L+1) - GCP(L))/GCP(L+1)                  
CC      GCP(L+1) =  GC(L)*GCP(L+1) - GC(L+1)*T                        
CC      FC (L+1) = W*FC (L+1)                                          
CC9     FCP(L+1) = W*FCP(L+1)                                         
C  REPLACE CARDS 189 - 193 BY THE NEXT 11 CARDS                        
      IF(NUM.EQ.1) GO TO 9                                             
      IF(NUM.EQ.2) SL = ETA/FLOAT(L) + FLOAT(L)*RI                     
      IF(NUM.EQ.3) SL = GCP(L+1)                                      
      RL = GC(L+1)                                                  
      GC(L+1) = (SL*GC(L) - GPL)/RL                                  
      GPL1    =  RL*GC(L) - SL*GC(L+1)                                
      GPL     =  GPL1                                                 
      IF(NUM.LT.3) GO TO 9                                           
      GCP(L+1) = GPL1                                                 
      FCP(L+1) = FCP(L+1)*W                                           
9     FC (L+1) = FC(L+1)*W                                            
CC      FC (LMIN1) = FC (LMIN1)*W                                       
CC      FCP(LMIN1) = FCP(LMIN1)*W                                      
CC      RETURN                                                        
CC10    FC (LMIN1) = W*F                                             
CC      FCP(LMIN1) = W*FP                                              
      RETURN                                                           
11    W  = 0.d0                                                         
      G  = 0.d0                                                         
      GP = 0.d0                                                         
      GO TO 8                                                          
      END                                                             
C     ******************************************************************
      SUBROUTINE RCWFN(RHS,ETS,MINL,MAXL,FC,FCP,GC,GCP,ACCUR,STEPC)
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      DOUBLE PRECISION K,K1,K2,K3,K4,M1,M2,M3,M4
      DIMENSION FC(1),FCP(1),GC(1),GCP(1)
      DATA GPMAX,ONE,HUNDR,C1,C2,HALF,TWO,ZER,SIX,TWOUS,FSOUS,
     *OTWO,THRN,OONE/1.0D+60,1.0D+00,0.1D+02,1.0D-12,1.0D-06,
     *0.5D+00,0.2D+01,0.0D+00,0.6D+01,0.2D+05,0.46D+05,0.2D+00,
     *0.999D+03,0.1D-02/
      PAGE=STEPC
      ACC=ACCUR
      RHO=RHS
      ETA=ETS
      IF(PAGE.LT.HUNDR) PAGE=HUNDR
      IF(ACC.LT.C1.OR.ACC.GT.C2) ACC=C2
      R=RHO
      KTR=1
      LMAX=MAXL
      LMIN1=MINL+1
      XLL1=DFLOAT(MINL*LMIN1)
      ETA2=ETA*ETA
      TURN=ETA+SQRT(ETA2+XLL1)
      IF(R.LT.TURN.AND.ABS(ETA).GE.C2) KTR=-1
      KTRP=KTR
      GO TO 2
    1 R=TURN
      TF=F
      TFP=FP
      LMAX=MINL
      KTRP=1
    2 ETAR=ETA*R
      RHO2=R*R
      PL=DFLOAT(LMAX+1)
      PMX=PL+HALF
      FP=ETA/PL+PL/R
      DK=ETAR*TWO
      DEL=ZER
      D=ZER
      F=ONE
      K=(PL*PL-PL+ETAR)*(TWO*PL-ONE)
      IF((PL*PL+PL+ETAR).NE.ZER) GO TO 3
      R=R+C2
      GO TO 2
    3 H=(PL*PL+ETA2)*(ONE-PL*PL)*RHO2
      K=K+DK+PL*PL*SIX
      D=ONE/(D*H+K)
      DEL=DEL*(D*K-ONE)
      IF(PL.LT.PMX) DEL=-R*(PL*PL+ETA2)*(PL+ONE)*D/PL
      PL=PL+ONE
      FP=FP+DEL
      IF(D.LT.ZER) F=-F
      IF(PL.GT.TWOUS) GO TO 11
      IF(ABS(DEL/FP).GE.ACC) GO TO 3
      FP=F*FP
      IF(LMAX.EQ.MINL) GO TO 5
      FC(LMAX+1)=F
      FCP(LMAX+1)=FP
      L=LMAX
      DO 4 LP=LMIN1,LMAX
      PL=DFLOAT(L)
      GC(L+1)=ETA/PL+PL/R
      GCP(L+1)=SQRT(ETA2+PL*PL)/PL
      FC(L)=(GC(L+1)*FC(L+1)+FCP(L+1))/GCP(L+1)
      FCP(L)=GC(L+1)*FC(L)-GCP(L+1)*FC(L+1)
    4 L=L-1
      F=FC(LMIN1)
      FP=FCP(LMIN1)
    5 IF(KTRP.EQ.-1) GO TO 1
      P=ZER
      Q=R-ETA
      PL=ZER
      AR=-(ETA2+XLL1)
      AI=ETA
      BR=TWO*Q
      BI=TWO
      WI=TWO*ETA
      DR=BR/(BR*BR+BI*BI)
      DI=-BI/(BR*BR+BI*BI)
      DP=-(AR*DI+AI*DR)
      DQ=AR*DR-AI*DI
    6 P=P+DP
      Q=Q+DQ
      PL=PL+TWO
      AR=AR+PL
      AI=AI+WI
      BI=BI+TWO
      D=AR*DR-AI*DI+BR
      DI=AI*DR+AR*DI+BI
      T=ONE/(D*D+DI*DI)
      DR=D*T
      DI=-T*DI
      H=BR*DR-BI*DI-ONE
      K=BI*DR+BR*DI
      T=DP*H-DQ*K
      DQ=DP*K+DQ*H
      DP=T
      IF(PL.GT.FSOUS) GO TO 11
      CNT=(ABS(DP)+ABS(DQ))-((ABS(P)+ABS(Q))*ACC)
      IF(CNT) 66,66,6
   66 P=P/R
      Q=Q/R
      G=(FP-P*F)/Q
      GP=P*G-Q*F
      W=ONE/SQRT(FP*G-F*GP)
      G=W*G
      GP=W*GP
      IF(KTR.EQ.1) GO TO 8
      F=TF
      FP=TFP
      LMAX=MAXL
      IF(RHO.LT.(OTWO*TURN)) PAGE=THRN
      R3=0.333333333333333333333D+00
      H=(RHO-TURN)/(PAGE+ONE)
      H2=HALF*H
      I2=PAGE+OONE
      ETAH=ETA*H
      H2LL=H2*XLL1
      S=(ETAH+H2LL/R)/R-H2
    7 RH2=R+H2
      T=(ETAH+H2LL/RH2)/RH2-H2
      K1=H2*GP
      M1=S*G
      K2=H2*(GP+M1)
      M2=T*(G+K1)
      K3=H*(GP+M2)
      M3=T*(G+K2)
      M3=M3+M3
      K4=H2*(GP+M3)
      RH=R+H
      S=(ETAH+H2LL/RH)/RH-H2
      M4=S*(G+K3)
      G=G+(K1+K2+K2+K3+K4)*R3
      GP=GP+(M1+M2+M2+M3+M4)*R3
      R=RH
      I2=I2-1
      GPG=GP
      IF(ABS(GPG).GT.GPMAX) GO TO 11
      IF(I2.GE.0) GO TO 7
      W=ONE/(FP*G-F*GP)
    8 GC(LMIN1)=G
      GCP(LMIN1)=GP
      IF(LMAX.EQ.MINL) GO TO 10
      DO 9 L=LMIN1,LMAX
      T=GC(L+1)
      GC(L+1)=(GC(L)*GC(L+1)-GCP(L))/GCP(L+1)
      GCP(L+1)=GC(L)*GCP(L+1)-GC(L+1)*T
      FC(L+1)=W*FC(L+1)
    9 FCP(L+1)=W*FCP(L+1)
      FCP(LMIN1)=FCP(LMIN1)*W
      FC(LMIN1)=FC(LMIN1)*W
      GO  TO 12
   10 FC(LMIN1)=W*F
      FCP(LMIN1)=W*FP
      GO TO 12
   11 W=ZER
      G=ZER
      GP=ZER
      GO TO 8
   12 RETURN
      END
C     ***********************************************************
      SUBROUTINE LOGMO
C     ***********************************************************

      INCLUDE 'PRIVCOM7.FOR'
      N=0
      KJ1=IABS(JSO-JSS)+2
      KJ2=JSO+JSS+2
      DO 1 JO2=KJ1,KJ2,2
      JO=JO2-2
      KL1=IABS(JO-NSPI)+2
      KL2=JO+NSPI+2
      DO 2 LO2=KL1,KL2,2
      LO=LO2-2
      IF(NPIO*NPIS*(-1)**(LO/2).LE.0) GO TO 2
      N=N+1
      LNO(N)=LO/2
      JNO(N)=JO
      NS1(N)=JO-LO
    2 CONTINUE
    1 CONTINUE
      NN1=N
      N1M=NN1-1
      IF(NN1.LE.1) GO TO 5
      DO 3 IM=1,N1M
      LM=LNO(IM)
      IN=IM+1
      DO 4 IL=IN,NN1
      LA=LNO(IL)
      IF(LM.LE.LA) GO TO 4
      JA=JNO(IL)
      NS1A=NS1(IL)
      JNO(IL)=JNO(IM)
      LNO(IL)=LM
      NS1(IL)=NS1(IM)
      JNO(IM)=JA
      LNO(IM)=LA
      NS1(IM)=NS1A
      LM=LA
    4 CONTINUE
    3 CONTINUE
    5 CONTINUE
      RETURN
      END
C     *****************************************************************
      SUBROUTINE KLEGO1
C     *****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      COMMON/LOFAC/A(800)
      INCLUDE 'PRIVCOM1.FOR'
C     COMMON/KG/J1,J2,M1,M2,J,M,AKG

      IF(M1+M2-M)1,2,1
    1 AKG=0.
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
      AKG=0.
      IF(NB.LT.NM) GO TO 14
      DO 13 I1=NM,NB,2
      C1=1.
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
   13 AKG=AKG+C1*EXP(CL1)
      AKG=AKG*SQRT(J+1.D0)
   14 RETURN
      END
C     *******************************************************
      SUBROUTINE COPHA
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
C     COMMON/COUL/CONZ,ETA,COPH(40,90)
      INCLUDE 'PRIVCOM14.FOR'  
      DO 1 L1=1,1000
      ARG=ETA/L1
      ALF=ATAN(ARG)
      BETE=SQRT(ETA**2+L1**2)
      CO=ALF*(L1-0.5D0)+ETA*(LOG(BETE)-1.D0)-SIN(ALF)/12.D0/BETE
      SUM7=SIN(9.D0*ALF)/1188.D0/BETE**9
      ERROR=ABS(SUM7/CO)
      IF(ERROR.LT.1D-10) GO TO 2
    1 CONTINUE
      PRINT 24,ERROR
      WRITE(21,24)ERROR
   24 FORMAT(10X,'WARNING! ERROR IN COULOMB PHASE>',D12.5)
    2 LL1=L1
      COL=CO+SIN(3.D0*ALF)/360.D0/BETE**3-SIN(5.D0*ALF)/1260.D0/BETE**5
     *+SIN(7.D0*ALF)/1680.D0/BETE**7-SUM7
      LL=LL1
      IF(LL.EQ.1) GO TO 5
    4 LL=LL-1
      COL=COL-ATAN(ETA/LL)
      IF(LL.GT.1) GO TO 4
    5 COPH(40,1)=COL
      RETURN
      END
C     *******************************************************
      SUBROUTINE SPHEPOT
C     ***********************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM3.FOR'
      DE=RZ+20.*AZ
      RGR=DE
      IF(R.GE.RGR) DE=1.D0/R
      IF(R.GE.RGR) GO TO 22
      CALL SPHER
      ANO=DE
C     write(22,33) r,pp,de
C  33 FORMAT(6E12.4)
      DE=RZ+20.*AZ
      AN=DE/AZ+2.
      NAN1=AN*10
      NAN1=2*NAN1+2
      IF(NAN1.GT.96) NAN1=96
      NAN2=NAN1+1
      ST=DE/(NAN2-1)
      VV=0.D+00
      PDZ=1./DE**2
      CALL SPHER
      PDZ=PDZ*DE
    5 VV=VV+PDZ
      V4=0.
      DO 10 N=2,NAN2,2
      DE=ST*(N-1)
      PDZ=1./DE**2
      CALL SPHER
      PDZ=PDZ*DE
   10 V4=V4+PDZ
      V2=0.
      DO 11 N=3,NAN1,2
      DE=ST*(N-1)
      PDZ=1./DE**2
      CALL SPHER
      PDZ=PDZ*DE
   11 V2=V2+PDZ
      DE=(VV+4.*V4+2.*V2)*ST/3.
      PO0=DE/ANO+1./RGR
      DE=R
      IF(R.EQ.0.) PO1=0.
      IF(R.EQ.0.) GO TO 18
      AN=DE/AZ+2.
      NAN1=AN*5
      NAN1=2*NAN1+2
      IF(NAN1.GT.48) NAN1=48
      NAN2=NAN1+1
      ST=DE/(NAN2-1)
      VV=0.D+00
      PDZ=1./DE**2
      CALL SPHER
      PDZ=PDZ*DE
      VV=VV+PDZ
      V4=0.
      DO 12 N=2,NAN2,2
      DE=ST*(N-1)
      PDZ=1./DE**2
      CALL SPHER
      PDZ=PDZ*DE
   12 V4=V4+PDZ
      V2=0.
      DO 13 N=3,NAN1,2
      DE=ST*(N-1)
      PDZ=1./DE**2
      CALL SPHER
      PDZ=PDZ*DE
   13 V2=V2+PDZ
      PO1=(VV+4.*V4+2.*V2)*ST/3.
   18 DE=PO0-PO1/ANO
   22 CONTINUE
C     PRINT 33,R,PP, DE,PO0,PO1,ANO
      RETURN
      END
C     ***********************************************************
      SUBROUTINE SPHER
C     ***********************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM3.FOR' 
      IF(DE.EQ.0.) GO TO 15
      AN=DE/AZ+2.
      NAN1=AN*10
      NAN1=2*NAN1+2
      IF(NAN1.GT.96)NAN1=96
      NAN2=NAN1+1
      ST=DE/(NAN2-1)
C     PRINT 88,ST,NAN2
C     STOP
      VV=0.D+00
      X=(DE-RZ)/AZ
      IF(X.GT.23.) GO TO 4
      EX=EXP(X)
      EX1=1.D0+EX
      EXD=1.D0/EX1
      PDZ=EXD*DE**2
      GO TO 5
    4 PDZ=EXP(-X)*DE**2
    5 VV=VV+PDZ
      V4=0.
      DO 10 N=2,NAN2,2
      Y=ST*(N-1)
      X=(Y-RZ)/AZ
      IF(X.GT.23.) GO TO 6
      EX=EXP(X)
      EX1=1.D0+EX
      EXD=1.D0/EX1
      PDZ=EXD*Y**2
   88 FORMAT(4E12.3)
      GO TO 10
    6 PDZ=EXP(-X)*Y**2
   10 V4=V4+PDZ
      V2=0.
      DO 11 N=3,NAN1,2
      Y=ST*(N-1)
      X=(Y-RZ)/AZ
      IF(X.GT.23.) GO TO 7
      EX=EXP(X)
      EX1=1.D0+EX
      EXD=1.D0/EX1
      PDZ=EXD*Y**2
      GO TO 11
    7 PDZ=EXP(-X)*Y**2
   11 V2=V2+PDZ
      DE=(VV+4.D0*V4+2.D0*V2)*ST/3.D0
   15 RETURN
      END
C     ***********************************************************
      SUBROUTINE POTVOL
C     ***********************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM6.FOR' 
      INCLUDE 'PRIVCOM10.FOR' 
      INCLUDE 'PRIVCOM13.FOR'        
 
C     CALCULATES INTEGRALS OF FORM FACTOR, IT'S FIRST DERIVATIVE
C
C     MULTIPLIED BY (-1)*(Ri/Ai) AND IT'S SECOND DERIVATIVE
C     MULTIPLIED BY (1/2)*(Ri/Ai)**2
C     INTEGRALS BY RADIUS FROM ZERO TO INFINITY. TO GET VOLUME
C     INTEGRALS MULTIPLY BY 4*PI
      CBET0=1.
      IS=2
      RRR=RR
       IF(MERRR.EQ.1)   RR=RR*(1.D00-RRBWC*EN**PDIS/(EN**2+RRWID**PDIS))
      AR=AR0+AR1*EN
      AW=AW0+AW1*EN
      AC=AC0+AC1*EN
      AD=AD0+AD1*EN
      IF (EN.GT.BNDC) AD=AD0+AD1*BNDC
      RM=RR
      RAR=RR/AR
      RAD=RD/AD
      RAW=RW/AW
      RAC=RC/AC
      IF(RW.GT.RM) RM=RW
      IF(RC.GT.RM) RM=RC
      IF(RD.GT.RM) RM=RD
      AM=AR
      IF(AW.LT.AM) AM=AW
      IF(AC.LT.AM) AM=AC
      IF(AD.LT.AM) AM=AD
      AL=AR
      IF(AW.GT.AL) AL=AW
      IF(AC.GT.AL) AL=AC
      IF(AD.GT.AL) AL=AD
      AN=RM/AM+23.*AL/AM
      NAN1=AN*15
      NAN1=2*NAN1+2
      IF(NAN1.GT.148)NAN1=148
      NAN2=NAN1+1
      ST=(RM+23.*AL)/(NAN2-1)
      A42=4.D0*ST/3.D0
C
      VIR0=0.
      VIR1=0.
      VIR2=0.
      VIR3=0.
      WIC0=0.
      WIC1=0.
      WIC2=0.
      WIC3=0.
      WID0=0.
      WID1=0.
      WID2=0.
      WID3=0.
      WIW0=0.
      WIW1=0.
      WIW2=0.
      WIW3=0.
   16 DO 10 N=IS,NAN2,2
      R=ST*(N-1)
      X=(R-RR)/AR
      IF(X.GT.23.) GO TO 30
      EX=EXP(X)
      EX1=1.+EX
      EXD=1./EX1
      EXDM=1.-EXD
      PV0=EXD
      IF(LAS.EQ.0) GO TO 31
      PV1=PV0*EXDM*RAR
      IF(LAS.EQ.1) GO TO 31
      Y1=1.-2.*EXD
      PV2=PV1*Y1*RAR/2.
      IF(LAS.EQ.2) GO TO 31
      Y2=1.-6.*EXD*EXDM
      PV3=PV2/Y1*Y2*RAR/3.
      IF(LAS.EQ.3) GO TO 31
      Y3=1.-EXD*(14.-EXD*(36.-24.*EXD))
      PV4=PV3/Y2*Y3*RAR/4.
      GO TO 31
   30 EX1=EXP(-X)
      PV0=EX1
      IF(LAS.EQ.0) GO TO 31
      PV1=PV0*RAR
      IF(LAS.EQ.1) GO TO 31
      PV2=PV1*RAR/2.
      IF(LAS.EQ.2) GO TO 31
      PV3=PV2*RAR/3.
      IF(LAS.EQ.3) GO TO 31
      PV4=PV3*RAR/4.
   31 CONTINUE
      X=(R-RC)/AC
      Y=(R-RW)/AW
      XX=Y*Y
      IF(XX.GT.180.) EXX=0.
      IF(XX.LE.180.) EXX=EXP(-XX)
      IF(X.GT.23.) GO TO 32
      EX=EXP(X)
      EX1=1.+EX
      EXD=1./EX1
      EXDM=1.-EXD
      WCC=ALF*EXD
      WWW=(1.-ALF)*EXX
      PWC0=WCC
      PWW0=WWW
      IF(LAS.EQ.0) GO TO 33
      WCC=WCC*EXDM*RAC
      WWW=WWW*2.*RAW*Y
      PWC1=WCC
      PWW1=WWW
      IF(LAS.EQ.1) GO TO 33
      YC1=1.-2.*EXD
      WCC=WCC*YC1*RAC/2.
      Y1=2.*Y*Y-1.
      WWW=WWW/Y*Y1*RAW/2.
      PWC2=WCC
      PWW2=WWW
      IF(LAS.EQ.2) GO TO 33
      YC2=1.-6.*EXD*EXDM
      WCC=WCC/YC1*YC2*RAC/3.
      Y2=Y**3*2.-Y*3.
      WWW=WWW/Y1*Y2*RAW*2./3.
      PWC3=WCC
      PWW3=WWW
      IF(LAS.EQ.3) GO TO 33
      YC3=1.-EXD*(14.-EXD*(36.-24.*EXD))
      WCC=WCC/YC2*YC3*RAC/4.
      WWW=WWW/Y2*(Y**4*4.-Y*Y*12.+3.)*RAW/4.
      PWC4=WCC
      PWW4=WWW
      GO TO 33
   32 EX1=EXP(-X)
      WCC=ALF*EX1
      WWW=(1.-ALF)*EXX
      PWC0=WCC
      PWW0=WWW
      IF(LAS.EQ.0) GO TO 33
      WCC=WCC*RAC
      WWW=WWW*2.*RAW*Y
      PWC1=WCC
      PWW1=WWW
      IF(LAS.EQ.1) GO TO 33
      WCC=WCC*RAC/2.
      Y1=2.*Y*Y-1.
      WWW=WWW/Y*Y1*RAW/2.
      PWC2=WCC
      PWW2=WWW
      IF(LAS.EQ.2) GO TO 33
      WCC=WCC*RAC/3.
      Y2=Y**3*2.-Y*3.
      WWW=WWW/Y1*Y2*RAW*2./3.
      PWC3=WCC
      PWW3=WWW
      IF(LAS.EQ.3) GO TO 33
      WCC=WCC*RAC/4.
      WWW=WWW/Y2*(Y**4*4.-Y*Y*12.+3.)*RAW/4.
      PWC4=WCC
      PWW4=WWW
   33 X=(R-RD)/AD
      IF(X.GT.23.) GO TO 34
      EX=EXP(X)
      EX1=1.+EX
      EXD=1./EX1
      EXDM=1.-EXD
      WDD=4.*EXDM*EXD
      PWD0=WDD
      IF(LAS.EQ.0) GO TO 35
      Y1=1.-2.*EXD
      WDD=WDD*Y1*RAD
      PWD1=WDD
      IF(LAS.EQ.1) GO TO 35
      Y2=1.-6.*EXD*EXDM
      WDD=WDD/Y1*Y2*RAD/2.
      PWD2=WDD
      IF(LAS.EQ.2) GO TO 35
      Y3=1.-EXD*(14.-EXD*(36.-24.*EXD))
      WDD=WDD/Y2*Y3*RAD/3.
      PWD3=WDD
      IF(LAS.EQ.3) GO TO 35
      Y4=1.-EXD*(30.-EXD*(150.-EXD*(240.-120.*EXD)))
      WDD=WDD/Y3*Y4*RAD/4.
      PWD4=WDD
      GO TO 35
   34 EX1=EXP(-X)
      WDD=4.*EX1
      PWD0=WDD
      IF(LAS.EQ.0) GO TO 35
      WDD=WDD*RAD
      PWD1=WDD
      IF(LAS.EQ.1) GO TO 35
      WDD=WDD*RAD/2.
      PWD2=WDD
      IF(LAS.EQ.2) GO TO 35
      WDD=WDD*RAD/3.
      PWD3=WDD
      IF(LAS.EQ.3) GO TO 35
      WDD=WDD*RAD/4.
      PWD4=WDD
   35 CONTINUE
      VIR0=VIR0+R**2*PV0*A42
      VIR1=VIR1+R**2*PV1*A42
      VIR2=VIR2+R**2*PV2*A42
      VIR3=VIR3+R**2*PV3*A42
      WIC0=WIC0+R**2*PWC0*A42
      WIC1=WIC1+R**2*PWC1*A42
      WIC2=WIC2+R**2*PWC2*A42
      WIC3=WIC3+R**2*PWC3*A42
      WID0=WID0+R**2*PWD0*A42
      WID1=WID1+R**2*PWD1*A42
      WID2=WID2+R**2*PWD2*A42
      WID3=WID3+R**2*PWD3*A42
      WIW0=WIW0+R**2*PWW0*A42
      WIW1=WIW1+R**2*PWW1*A42
      WIW2=WIW2+R**2*PWW2*A42
   10 WIW3=WIW3+R**2*PWW3*A42
      IS=IS+1
      NAN2=NAN1
      A42=2.*ST/3.
      IF(IS.LE.3) GO TO 16
      CBETC=0.d0
      CBETD=0.d0
      CBETW=0.d0
      IF(WIC1.NE.0.)CBETC=WIC2/WIC1
      IF(WID1.NE.0.)CBETD=WID2/WID1
      IF(WIW1.NE.0.)CBETW=WIW2/WIW1
      IF(MEVOL.LE.1) GO TO 3
      CBET0=VIR2/VIR1
    3 WRITE(21,14)VIR0,VIR1,VIR2,VIR3,CBET0
   14 FORMAT(/2X,'SPHERICAL VOLUME INTEGRALS OF REAL POTENTIAL F-FACTORS
     * AND DERIVATIVES:'/2X,'VIR0=',F8.3,2X,'VIR1=',F8.3,2X,
     *'VIR2=',F8.3,2X,'VIR3=',F8.3,2X,'CBET0=',F8.3)
      WRITE(21,15)WIC0,WIC1,WIC2,WIC3,CBETC
   15 FORMAT(2X,'SPHERICAL VOLUME INTEGRALS OF IMAGINARY (WC) F-FACTORS
     * AND DERIVATIVES:'/2X,'WIC0=',F8.3,2X,'WIC1=',F8.3,2X,
     *'WIC2=',F8.3,2X,'WIC3=',F8.3,2X,'CBETC=',F8.3)
      WRITE(21,12)WID0,WID1,WID2,WID3,CBETD
   12 FORMAT(2X,'SPHERICAL VOLUME INTEGRALS OF IMAGINARY (WD) F-FACTORS
     * AND DERIVATIVES:'/2X,'WID0=',F8.3,2X,'WID1=',F8.3,2X,
     *'WID2=',F8.3,2X,'WID3=',F8.3,2X,'CBETD=',F8.3)
      WRITE(21,11)WIW0,WIW1,WIW2,WIW3,CBETW
   11 FORMAT(2X,'SPHERICAL VOLUME INTEGRALS OF IMAGINARY (WW) F-FACTORS
     * AND DERIVATIVES:'/2X,'WIW0=',F8.3,2X,'WIW1=',F8.3,2X,
     *'WIW2=',F8.3,2X,'WIW3=',F8.3,2X,'CBETW=',F8.3)
      RR=RRR
      RETURN
      END
C     ***********************************************************
      SUBROUTINE VHFROM
C     ***********************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM.FOR'
C     WW4 - ENERGY CONVERSION FACTOR DEVIDED BY 4.
C     ALAVR - POTENTIAL NON-LOCALITY IN FM.
C     ENCON - INCIDENT ENERGY
C     EMATCH - MATCHING ENERGY
C     VAL - VALUE OF VHF REAL POTENTIAL AT MATCHING POINT,
C      USING EXPONENTIAL DECAY FORMULAE
C     DIR - DERIVATIVE OF VHT POTENTIAL AT MATCHING POINT
C     USING EXPONENTIAL DECAY FORMULAE

      APAN=ANEU/1.00866491600D0 
      IFPOT=0
      EMATCH=200.00D0

       VRLNP=VRLANP
       ALAV=ALAVR
       go to 6

       IF(IFPOT.EQ.0) GO TO 7
      IF (EMATCH.LE.ENCON) GO TO 6
      REL=(EMATCH+AMI)/AMI
        IF(REL.LT.1.D00) REL=1.D00
      EIRELM=SQRT(REL**2*AT**2+2.D+0*REL*ANEU*AT+ANEU**2)
      EIRELT=AT*SQRT(REL**2*ANEU**2+2.D+0*REL*ANEU*AT+AT**2)
      DEREL=SQRT(2.D+0*REL*ANEU*AT+AT**2+ANEU**2)
C        ENCON=(REL**2-1.D0)*AMI*AT*(ANEU+AT)/DEREL**2/2.D+0
      WW4=1.2064496D-2*EIRELM*EIRELT/(ANEU*EIRELM+EIRELT)/DEREL*APAN

       VAL=VRLANP*EXP(-ALAVR**2*(EMATCH-EFERM)*WW4)
       DIR=-VAL*ALAVR**2*WW4*(1.D+00+(EMATCH-EFERM)/(EMATCH+AMI))


C     ALAVR - MARCHED NON-LOCALITY RANGE FOR PEREY BUCK FORMULA
C     VRLANP - MATCHED VHFO FOR PEREY BUCK FORMULA

       ALAVR=SQRT(-DIR/VAL/(1.D+00+DIR+(EMATCH+VAL)/(EMATCH+AMI))/WW4)
       VRLANP=VAL*EXP(ALAVR**2*(EMATCH+VAL)*WW4)

C     MORILLIOM ROMAIM FORMULAE
C     VHF AT EN BELOW EMATCH, CALCULATED BY PEREY BUCK FORMULA,
C      ACCOUNTING RELATIVISTIC EFFECTS
C     WITH PARAMETERS MATCHED WITH EXPONENTIAL DECAY DEPENDENCE

    7 REL=(ENCON+AMI)/AMI
        IF(REL.LT.1.D00) REL=1.D00
      EIRELM=SQRT(REL**2*AT**2+2.D+0*REL*ANEU*AT+ANEU**2)
      EIRELT=AT*SQRT(REL**2*ANEU**2+2.D+0*REL*ANEU*AT+AT**2)
      DEREL=SQRT(2.D+0*REL*ANEU*AT+AT**2+ANEU**2)
C        ENCON=(REL**2-1.D0)*AMI*AT*(ANEU+AT)/DEREL**2/2.D+0
      WW4=1.2064496D-2*EIRELM*EIRELT/(ANEU*EIRELM+EIRELT)/DEREL*APAN
      CONST=VRLANP*EXP(-WW4*ALAVR**2*ENCON)
       RANGE=0.000D+00
       WW4216=WW4**2*16.D+00
      VHF=0.D+00
       STEP=45.D+00
    3  VHF=VHF+STEP
      IF(VHF.GT.5.45D+02) GO TO 4
       VALUE=CONST*EXP(-WW4*ALAVR**2*VHF+WW4216*
     *(RANGE**4*(ENCON+VHF))**2)
       IF(VALUE-VHF) 1,2,3
    1 VHF=VHF-STEP
      STEP=STEP/3.0
      IF(STEP.LT.1.D-04) GO TO 2
      GO TO 3
    4 PRINT 5,ALAVR,ENCON,VALUE,WW4
      WRITE(21,5)ALAVR,ENCON,VALUE,WW4
    5 FORMAT(10X,'ERROR! CHECK ALAVR OR VHF SIGN:NO SOLUTION FOR VHF'/
     *5X,'ALAVR=',D11.5,
     *', ENCON=',D11.5,',VALUE=',D11.5,',WW4=',D11.5)
      STOP
    2 VRLANP=VRLNP
      ALAVR=ALAV
      RETURN

C     VHF AT EN AS EXPONENTIAL DECAY, WITH ALAVR AS NON-LOCALITY RANGE,
C     WITH REDUCED MASS ACCOUNTING RELATIVISTIC EFFECTS

    6 CONTINUE

C     CALCULATIONS OF RELATIVISTIC REDUCED MASS OF SYSTEM DEVIDED BY 2*h**2

      REL=(ENCON+AMI)/AMI
        IF(REL.LT.1.D00) REL=1.D00
      EIRELM=SQRT(REL**2*AT**2+2.D+0*REL*ANEU*AT+ANEU**2)
      EIRELT=AT*SQRT(REL**2*ANEU**2+2.D+0*REL*ANEU*AT+AT**2)
      DEREL=SQRT(2.D+0*REL*ANEU*AT+AT**2+ANEU**2)
C      ENCON=(REL**2-1.D0)*AMI*AT*(ANEU+AT)/DEREL**2/2.D+0
      WW4=1.2064496D-2*EIRELM*EIRELT/(ANEU*EIRELM+EIRELT)/DEREL*APAN

C     THIS IS NON RELATIVISTIC REDUCED MASS OF SYSTEM DEVIDED BY 2*h**2
      WW4=1.2064496D-2*AT/(ANEU+AT)*APAN

C     The line below is in the spirit of Morillon and Romain formulation
C     VHF=VRLANP*EXP(-ALAVR**2*(ENCON-EFERM)*WW4)
C
C     To follow our published papers on PRC
C
      VHF=VRLANP*EXP(-ALAVR*(ENCON-EFERM))
      RETURN
      END
C     *******************************************************
      SUBROUTINE THORA(IOUT)
C     *******************************************************
C     AUTHOR: R. Capote, March 2005
C     
C     * $Date: 2007/09/03 14:20:35 $
C     * $Id: thora.f,v 1.6 2007/09/03 14:20:35 Capote Exp $
C
C     GIVES THE TIME ELAPSED SINCE THE FIRST CALL
C     Note: Elapsed time must be less than one month
C
C     (g77, g95, LAHEY, and MS FORTRAN compatible)
C
      DOUBLE PRECISION BEGTIM,ENDTIM,BEGDAY,ENDDAY,DIFTIM,DIFTI1
      INTEGER IOUT
      CHARACTER*8 DATE
      CHARACTER*10 TIME
      CHARACTER*5 ZONE
      INTEGER DT
      DIMENSION DT(8)
      LOGICAL NEVER_CALLED
      DATA NEVER_CALLED/.TRUE./
      SAVE BEGTIM,BEGDAY,NEVER_CALLED

      IF (NEVER_CALLED) then
        NEVER_CALLED = .FALSE.
        CALL DATE_AND_TIME (DATE,TIME,ZONE,DT)
        BEGDAY=DT(3)
        BEGTIM=DT(5)*3600.d0+DT(6)*60.d0+DT(7)+DT(8)/1000.d0
        WRITE(IOUT,1001) time(1:2),time(3:4),time(5:6),
     >                   DATE(7:8),DATE(5:6),DATE(1:4)
      ELSE
        CALL DATE_AND_TIME (DATE,TIME,ZONE,DT)
        ENDTIM=DT(5)*3600.d0+DT(6)*60.d0+DT(7)+DT(8)/1000.d0
        ENDDAY=DT(3)
        ENDTIM = ENDTIM +  (ENDDAY-BEGDAY)*86400.d0
        DIFTIM=(ENDTIM-BEGTIM)/60.d0
        DIFTI1=(DIFTIM-INT(DIFTIM))*60.d0
        WRITE(IOUT,1002) INT(DIFTIM),NINT(DIFTI1)
      ENDIF

      RETURN
 1001 FORMAT
     >(/10X,'Start time: ',A2,':',A2,'.',A2,' (',A2,'-',A2,'-',A4,')'/)
 1002 FORMAT(//10X,' Calculation time: ',I4,' min ',I2,' s'//)
C====================================================================
      END
C     *******************************************************
C     SUBROUTINE WRITTEN BY Rui LI AND Weili SUN FROM IAPCM CHINA
	SUBROUTINE ANPOW !analyzing power
C     *******************************************************
      IMPLICIT REAL*8(A-H,O-Z)
	!----------analyzing power---------------
c	REAL*8 FC_R(150),FC_I(150),X_R(4,40,50,50,150),
c     *       X_I(4,40,50,50,150),COEF_X(4),CSRX,CSIX,
c     *       CS_ANG_MIF(40,50,50,150),CS_ANG(40,150)	           
	REAL*8 FC_R,FC_I,X_R(4,40,50,50),
     *       X_I(4,40,50,50),COEF_X(4),CSRX,CSIX,
     *       CS_ANG_MIF(40,50,50),CS_ANG(40)	           
       INTEGER MII_INDEX,MII,MII_MAX,MIF_INDEX,MIF,MIF_MAX

     	!----analyzing power---------
      INCLUDE 'PRIVCOM.FOR'
      INCLUDE 'PRIVCOM1.FOR'
      INCLUDE 'PRIVCOM13.FOR'
      INCLUDE 'PRIVCOM14.FOR'           
 
    	!----------initialization----------
	X_R=0.D0
	X_I=0.D0
	FC_R=0.D0
	FC_I=0.D0
	DISC=0.D0
C	Px=0.
C	Pz=0.
	Ay=0.D0
C	open(unit=904,file='PxPz.txt')

        CS_ANG=0.D0
      !*************initialization***********
	CALL PLM_CAL !calculating harmonic function
      LKK=1
      ETA=CONZ/WNK(1)
      IC=JO(1)
	!-------------calculate f_c(Coulomb amplitude)-------------
	DO 627 MM=1,MTET 
	X_R=0.D0
	X_I=0.D0
	FC_R=0.D0
	FC_I=0.D0
	CS_ANG=0.D0
      CS_ANG_MIF=0.D0
      
	TETA=TET(MM)*3.1415927/180
	ALST2=DLOG(DSIN(TETA/2.D0))
      ARGC=2.D0*(COPH(1,1)-ETA*ALST2)
	temp=-ETA/(20.D0*WNK(1)*DSIN(TETA/2.D0)**2)*DCOS(ARGC)
c	FC_R(M)=temp 
	FC_R=temp 
	temp=-ETA/(20.D0*WNK(1)*DSIN(TETA/2.D0)**2)*DSIN(ARGC)
c	FC_I(M)=temp       
	FC_I=temp       
c  627	continue		
      !------------------------------------------------------------	
      DO 1 K1=1,NJ
      JS1=JS(K1)
      N1I=INC(K1)
      N1F=INR(K1)
      DO 1 N1C=1,N1I 
      L1I=LNJ(K1,N1C)
      J1I=JNJ(K1,N1C)
      DO 8 N1R=1,N1F
      L1F=LNJ(K1,N1R)              !L1F is not multiplied by 2
      J1F=JNJ(K1,N1R)
      NU1=NNJ(K1,N1R)

	MII_MAX=IC+1 !number of MII:2*I+1
	DO 614 MII_INDEX=1,MII_MAX
	MII=2*MII_INDEX-IC-2 !2*MI
      MIF_MAX=JO(NU1)+1 !number of MIF: 2*I'+1
      DO 615 MIF_INDEX=1,MIF_MAX
      MIF=2*MIF_INDEX-JO(NU1)-2 
      COEF_X=2.D0*DSQRT(3.1415927D0)/(10.D0*WNK(1))
     **DSQRT(1.D0*(2.D0*L1I+1))
      J1=2*L1I
      J2=1
      M1=0
      M2=1
      J=J1I
      M=1
      call KLEGO
      COEF_X(1)=COEF_X(1)*AKG
      COEF_X(2)=COEF_X(2)*AKG
      J1=J1I
      J2=IC
      M1=1
      M2=MII
      J=JS1
      M=MII+1
      CALL KLEGO
      COEF_X(1)=COEF_X(1)*AKG
      COEF_X(2)=COEF_X(2)*AKG
      J1=2*L1F
      J2=1
      M1=MII-MIF
      M2=1
      J=J1F
      M=1+MII-MIF
      CALL KLEGO
      COEF_X(1)=COEF_X(1)*AKG
      J1=2*L1F
      J2=1
      M1=2+MII-MIF
      M2=-1
      J=J1F
      M=1+MII-MIF
      CALL KLEGO
      COEF_X(2)=COEF_X(2)*AKG
      J1=J1F
      J2=JO(NU1) 
      M1=1+MII-MIF
      M2=MIF
      J=JS1
      M=MII+1
      CALL KLEGO
      COEF_X(1)=COEF_X(1)*AKG
      COEF_X(2)=COEF_X(2)*AKG   

      J1=2*L1I
      J2=1
      M1=0
      M2=-1
      J=J1I
      M=-1
      CALL KLEGO
      COEF_X(3)=COEF_X(3)*AKG
      COEF_X(4)=COEF_X(4)*AKG
      J1=J1I
      J2=IC
      M1=-1
      M2=MII
      J=JS1
      M=MII-1
      CALL KLEGO
      COEF_X(3)=COEF_X(3)*AKG
      COEF_X(4)=COEF_X(4)*AKG
      J1=2*L1F
      J2=1
      M1=-2+MII-MIF
      M2=1
      J=J1F
      M=-1+MII-MIF
      CALL KLEGO
      COEF_X(3)=COEF_X(3)*AKG
      J1=2*L1F
      J2=1
      M1=MII-MIF
      M2=-1
      J=J1F
      M=-1+MII-MIF
      CALL KLEGO
      COEF_X(4)=COEF_X(4)*AKG
      J1=J1F
      J2=JO(NU1) 
      M1=-1+MII-MIF
      M2=MIF
      J=JS1
      M=MII-1
      CALL KLEGO
      COEF_X(3)=COEF_X(3)*AKG
      COEF_X(4)=COEF_X(4)*AKG 

      	   CORR=DCOS(COPH(NU1,L1F+1)+COPH(1,L1I+1)) 
		   COII=DSIN(COPH(NU1,L1F+1)+COPH(1,L1I+1)) 
		   CSRX=CORR*CRD(K1,N1C,N1R)-COII*CID(K1,N1C,N1R)
		   CSIX=CORR*CID(K1,N1C,N1R)+COII*CRD(K1,N1C,N1R)
c      DO 616 M=1,MTET
	IF(NU1.EQ.1 .AND. L1I.EQ.0 .AND. J1I.EQ.1 .AND. JS1.EQ.(IC+1) 
     * .AND. L1F.EQ.0	.AND. J1F.EQ.1 .AND. MII.EQ.MIF)
     *THEN !add the Coulomb amplitude(only added once for different J1I and JS1)
	X_R(1,NU1,MII_INDEX,MIF_INDEX)=X_R(1,NU1,MII_INDEX,MIF_INDEX)+
     *FC_R
	X_I(1,NU1,MII_INDEX,MIF_INDEX)=X_I(1,NU1,MII_INDEX,MIF_INDEX)+
     *FC_I
	X_R(4,NU1,MII_INDEX,MIF_INDEX)=X_R(4,NU1,MII_INDEX,MIF_INDEX)+
     *FC_R
	X_I(4,NU1,MII_INDEX,MIF_INDEX)=X_I(4,NU1,MII_INDEX,MIF_INDEX)+
     *FC_I
	COUNT_FC=COUNT_FC+1
	ENDIF
      X_R(1,NU1,MII_INDEX,MIF_INDEX)=X_R(1,NU1,MII_INDEX,MIF_INDEX)+
     *CSRX*YLM_VAL(L1F,(MII-MIF)/2,MM)*COEF_X(1)
      X_I(1,NU1,MII_INDEX,MIF_INDEX)=X_I(1,NU1,MII_INDEX,MIF_INDEX)+
     *CSIX*YLM_VAL(L1F,(MII-MIF)/2,MM)*COEF_X(1)
      X_R(2,NU1,MII_INDEX,MIF_INDEX)=X_R(2,NU1,MII_INDEX,MIF_INDEX)+
     *CSRX*YLM_VAL(L1F,1+(MII-MIF)/2,MM)*COEF_X(2)
      X_I(2,NU1,MII_INDEX,MIF_INDEX)=X_I(2,NU1,MII_INDEX,MIF_INDEX)+
     *CSIX*YLM_VAL(L1F,1+(MII-MIF)/2,MM)*COEF_X(2)
	X_R(3,NU1,MII_INDEX,MIF_INDEX)=X_R(3,NU1,MII_INDEX,MIF_INDEX)+
     *CSRX*YLM_VAL(L1F,-1+(MII-MIF)/2,MM)*COEF_X(3)
      X_I(3,NU1,MII_INDEX,MIF_INDEX)=X_I(3,NU1,MII_INDEX,MIF_INDEX)+
     *CSIX*YLM_VAL(L1F,-1+(MII-MIF)/2,MM)*COEF_X(3)
      X_R(4,NU1,MII_INDEX,MIF_INDEX)=X_R(4,NU1,MII_INDEX,MIF_INDEX)+
     *CSRX*YLM_VAL(L1F,(MII-MIF)/2,MM)*COEF_X(4)
      X_I(4,NU1,MII_INDEX,MIF_INDEX)=X_I(4,NU1,MII_INDEX,MIF_INDEX)+
     *CSIX*YLM_VAL(L1F,(MII-MIF)/2,MM)*COEF_X(4)
c  616 CONTINUE
  615 CONTINUE  
  614 CONTINUE
    8 CONTINUE
    1 CONTINUE
            
      DO 619 N=NUI,NUF
	MII_MAX=IC+1
	DO 612 MII_INDEX=1,MII_MAX
      MIF_MAX=JO(N)+1
      DO 618 MIF_INDEX=1,MIF_MAX
c      DO 617 M=1,MTET
       CS_ANG_MIF(N,MII_INDEX,MIF_INDEX)=
     *X_R(1,N,MII_INDEX,MIF_INDEX)**2+ 
     *X_I(1,N,MII_INDEX,MIF_INDEX)**2+
     *X_R(2,N,MII_INDEX,MIF_INDEX)**2+
     *X_I(2,N,MII_INDEX,MIF_INDEX)**2+
     *X_R(3,N,MII_INDEX,MIF_INDEX)**2+ 
     *X_I(3,N,MII_INDEX,MIF_INDEX)**2+
     *X_R(4,N,MII_INDEX,MIF_INDEX)**2+
     *X_I(4,N,MII_INDEX,MIF_INDEX)**2
      CS_ANG(N)=CS_ANG(N)+CS_ANG_MIF(N,MII_INDEX,MIF_INDEX)
     */MII_MAX/2.D0 !calculate the cross-section

       
      IF(KEYAP.EQ.1) GO TO 625
      
	DISC(N,MM)=DISC(N,MM)+    !calculate Py
     *2.D0*(X_R(1,N,MII_INDEX,MIF_INDEX)*X_I(2,N,MII_INDEX,MIF_INDEX
     *   )-X_I(1,N,MII_INDEX,MIF_INDEX)*X_R(2,N,MII_INDEX,MIF_INDEX
     *   ))/MII_MAX/2.D0 +
     *2.D0*(X_R(3,N,MII_INDEX,MIF_INDEX)*X_I(4,N,MII_INDEX,MIF_INDEX
     *   )-X_I(3,N,MII_INDEX,MIF_INDEX)*X_R(4,N,MII_INDEX,MIF_INDEX
     *   ))/MII_MAX/2.D0
     
      GO TO 617
	
  625 DISC(N,MM)=DISC(N,MM)+  !calculate Ay
     *2.D0*(X_R(1,N,MII_INDEX,MIF_INDEX)*X_I(3,N,MII_INDEX,MIF_INDEX
     *   )-X_I(1,N,MII_INDEX,MIF_INDEX)*X_R(3,N,MII_INDEX,MIF_INDEX
     *   ))/MII_MAX/2.D0 +
     *2.D0*(X_R(2,N,MII_INDEX,MIF_INDEX)*X_I(4,N,MII_INDEX,MIF_INDEX
     *   )-X_I(2,N,MII_INDEX,MIF_INDEX)*X_R(4,N,MII_INDEX,MIF_INDEX
     *   ))/MII_MAX/2.D0
  617 CONTINUE
  618 CONTINUE
  612 CONTINUE
c      DO 620 M=1,MTET 
 
      !print *, '%%%', TET(MM),CS_ANG(N)
      
      DISC(N,MM)=DISC(N,MM)/CS_ANG(N)
      IF(KEYAP.EQ.1) DISC(N,MM)=-DISC(N,MM)
c  620 CONTINUE    
  619 CONTINUE 

  627	continue		      
      
      pause 333
	RETURN
      END SUBROUTINE ANPOW
C     *******************************************************
	SUBROUTINE PLM_CAL
C     *******************************************************
	IMPLICIT REAL*8(A-H,O-Z)
C	COMMON/PLM/PLM_DATA(181,181,150)
C	COMMON/DISK/TET(150),MTET
      INCLUDE 'PRIVCOM1.FOR'

	PLM_DATA(1,1,:)=1./DSQRT(2.D0)
	DO 636 L=1,180
	C=1.D0
	DO 638 INDEX_L=1,L
	C=C*DSQRT(1.D0-1.D0/(2.D0*INDEX_L))
  638 CONTINUE

	DO 637 M=1,MTET
	TETA=TET(M)*3.1415927D0/180.D0
	COS_TETA=DCOS(TETA)
	SIN_TETA=DSIN(TETA)
	COT_TETA=COS_TETA/SIN_TETA
	PLM_DATA(L+1,L+1,M)=C*(-1.D0*SIN_TETA)**L !P L L
	PLM_DATA(L+1,L,M)=PLM_DATA(L+1,L+1,M)*(-2.D0*COT_TETA)*L/
     *  DSQRT(2.D0*L) !P L L-1
	DO 639 LM=L-2,0,-1 !recursion from l-2 to 0 order 
	PLM_DATA(L+1,LM+1,M)=(-2.D0*COT_TETA*(LM+1)*PLM_DATA(L+1,LM+2,M)-
     *DSQRT(1.D0*(L+LM+2))*DSQRT(1.D0*(L-LM-1))*PLM_DATA(L+1,LM+3,M))/
     *(DSQRT(1.D0*(L+LM+1))*DSQRT(1.D0*(L-LM)))
  639 CONTINUE !LM LOOP

  637 CONTINUE !ANGULAR LOOP
      PLM_DATA(L+1,1:L+1,:)=PLM_DATA(L+1,1:L+1,:)*DSQRT(L+0.5D0)
      IFLAG=1.
      DO 635 LM=0,L
      PLM_DATA(L+1,LM+1,:)=PLM_DATA(L+1,LM+1,:)*IFLAG
      IFLAG=-IFLAG
  635 CONTINUE !LM2 LOOP
  636 CONTINUE !L LOOP
	END SUBROUTINE PLM_CAL
C     *******************************************************	
	FUNCTION YLM_VAL(Y_L,Y_M,Y_ANGU) !find the value according to l and m
C     *******************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C     COMMON/PLM/PLM_DATA(181,181,150)
      INCLUDE 'PRIVCOM1.FOR'
C
      REAL*8 YLM_VAL      
      INTEGER Y_ANGU,Y_L,Y_M,YM_TEMP !Y_ANGU is the index of angular 
      
      IF(Y_M.ge.0) THEN
      YLM_VAL=(-1)**Y_M*PLM_DATA(Y_L+1,Y_M+1,Y_ANGU)
     * /DSQRT(2.D0*3.1415927D0) ! for m>=0 case
      ENDIF
      
      IF(Y_M.lt.0) THEN
       YM_TEMP=-1*Y_M
       YLM_VAL=PLM_DATA(Y_L+1,YM_TEMP+1,Y_ANGU)
     * /DSQRT(2.D0*3.1415927D0) ! for m<0 case
      ENDIF
      
      RETURN
      END FUNCTION YLM_VAL
C     ******************************************************* 

C     *******************************************************
C     SUBROUTINE WRITTEN BY Rui LI AND Weili SUN FROM IAPCM CHINA
	SUBROUTINE DEUTR !analyzing power
C     *******************************************************
      IMPLICIT REAL*8(A-H,O-Z)
	!----deutron angular distribution1---------
	REAL*8 FC_R(150),FC_I(150),X_R(9,10,20,20,150),
     *       X_I(9,10,20,20,150),COEF_X(9),CSRX,CSIX,
     *       CS_ANG_MIF(9,20,20,150),CS_ANG(9,150)	           
       INTEGER MII_INDEX,MII,MII_MAX,MIF_INDEX,MIF,MIF_MAX


     	!----deutron angular distribution2---------
     	      INCLUDE 'PRIVCOM.FOR'
            INCLUDE 'PRIVCOM1.FOR'
            INCLUDE 'PRIVCOM13.FOR' 
            INCLUDE 'PRIVCOM14.FOR'    
 
          	!----------initialization----------
	X_R=0.D0
	X_I=0.D0
	FC_R=0.D0
	FC_I=0.D0
	DISC=0.D0


      CS_ANG=0.D0
      !*************initialization***********
	CALL PLM_CAL !calculating harmonic function
      LKK=1
      ETA=CONZ/WNK(1)
      IC=JO(1)
	!-------------calculate f_c(Coulomb amplitude)-------------
	DO 627 M=1,MTET 
	TETA=TET(M)*3.1415927/180
	ALST2=DLOG(DSIN(TETA/2.D0))
      ARGC=2.D0*(COPH(1,1)-ETA*ALST2)
	temp=-ETA/(20.D0*WNK(1)*DSIN(TETA/2.D0)**2)*DCOS(ARGC)
	FC_R(M)=temp 
	temp=-ETA/(20.D0*WNK(1)*DSIN(TETA/2.D0)**2)*DSIN(ARGC)
	FC_I(M)=temp       
  627	continue		
      !------------------------------------------------------------	
      DO 1 K1=1,NJ
      JS1=JS(K1)
      N1I=INC(K1)
      N1F=INR(K1)
      DO 1 N1C=1,N1I 
      L1I=LNJ(K1,N1C)
      J1I=JNJ(K1,N1C)
      DO 8 N1R=1,N1F
      L1F=LNJ(K1,N1R)              !L1F is not multiplied by 2
      J1F=JNJ(K1,N1R)
      NU1=NNJ(K1,N1R)

	MII_MAX=IC+1 !number of MII:2*I+1
	DO 614 MII_INDEX=1,MII_MAX
	MII=2*MII_INDEX-IC-2 !2*MI
      MIF_MAX=JO(NU1)+1 !number of MIF: 2*I'+1
      DO 615 MIF_INDEX=1,MIF_MAX
      MIF=2*MIF_INDEX-JO(NU1)-2 
      COEF_X=2.D0*DSQRT(3.1415927D0)/(10.D0*WNK(1))
     **DSQRT(1.D0*(2.D0*L1I+1))
      J1=2*L1I
      J2=2
      M1=0
      M2=2
      J=J1I
      M=2
      call KLEGO
      COEF_X(1)=COEF_X(1)*AKG
      COEF_X(2)=COEF_X(2)*AKG
      COEF_X(3)=COEF_X(3)*AKG
      J1=J1I
      J2=IC
      M1=2
      M2=MII
      J=JS1
      M=MII+2
      CALL KLEGO
      COEF_X(1)=COEF_X(1)*AKG
      COEF_X(2)=COEF_X(2)*AKG
      COEF_X(3)=COEF_X(3)*AKG
      J1=2*L1F
      J2=2
      M1=MII-MIF
      M2=2
      J=J1F
      M=2+MII-MIF
      CALL KLEGO
      COEF_X(1)=COEF_X(1)*AKG
C     BUG RCN 19 June 2014
C     J1=2*LIF
      J1=2*L1F
      J2=2
      M1=2+MII-MIF
      M2=0
      J=J1F
      M=2+MII-MIF
      CALL KLEGO
      COEF_X(2)=COEF_X(2)*AKG
      J1=2*L1F
      J2=1
      M1=4+MII-MIF
      M2=-2
      J=J1F
      M=2+MII-MIF
      CALL KLEGO
      COEF_X(3)=COEF_X(3)*AKG
      J1=J1F
      J2=JO(NU1) 
      M1=2+MII-MIF
      M2=MIF
      J=JS1
      M=MII+2
      CALL KLEGO
      COEF_X(1)=COEF_X(1)*AKG
      COEF_X(2)=COEF_X(2)*AKG 
      COEF_X(3)=COEF_X(3)*AKG  
c----------------------------
      J1=2*L1I
      J2=2
      M1=0
      M2=0
      J=J1I
      M=0
      CALL KLEGO
      COEF_X(4)=COEF_X(4)*AKG
      COEF_X(5)=COEF_X(5)*AKG
      COEF_X(6)=COEF_X(6)*AKG
      J1=J1I
      J2=IC
      M1=0
      M2=MII
      J=JS1
      M=MII
      CALL KLEGO
      COEF_X(4)=COEF_X(4)*AKG
      COEF_X(5)=COEF_X(5)*AKG
      COEF_X(6)=COEF_X(6)*AKG
      J1=2*L1F
      J2=2
      M1=-2+MII-MIF
      M2=2
      J=J1F
      M=MII-MIF
      CALL KLEGO
      COEF_X(4)=COEF_X(4)*AKG
      J1=2*L1F
      J2=2
      M1=MII-MIF
      M2=0
      J=J1F
      M=MII-MIF
      CALL KLEGO
      COEF_X(5)=COEF_X(5)*AKG
      J1=2*L1F
      J2=2
      M1=2+MII-MIF
      M2=-2
      J=J1F
      M=MII-MIF
      CALL KLEGO
      COEF_X(6)=COEF_X(6)*AKG
      J1=J1F
      J2=JO(NU1) 
      M1=MII-MIF
      M2=MIF
      J=JS1
      M=MII
      CALL KLEGO
      COEF_X(4)=COEF_X(4)*AKG
      COEF_X(5)=COEF_X(5)*AKG 
      COEF_X(6)=COEF_X(6)*AKG
c----------------------------
      J1=2*L1I
      J2=2
      M1=0
      M2=-2
      J=J1I
      M=-2
      CALL KLEGO
      COEF_X(7)=COEF_X(7)*AKG
      COEF_X(8)=COEF_X(8)*AKG
      COEF_X(9)=COEF_X(9)*AKG
      J1=J1I
      J2=IC
      M1=-2
      M2=MII
      J=JS1
      M=MII-2
      CALL KLEGO
      COEF_X(7)=COEF_X(7)*AKG
      COEF_X(8)=COEF_X(8)*AKG
      COEF_X(9)=COEF_X(9)*AKG
      J1=2*L1F
      J2=2
      M1=-4+MII-MIF
      M2=2
      J=J1F
      M=-2+MII-MIF
      CALL KLEGO
      COEF_X(7)=COEF_X(7)*AKG
      J1=2*L1F
      J2=2
      M1=-2+MII-MIF
      M2=0
      J=J1F
      M=-2+MII-MIF
      CALL KLEGO
      COEF_X(8)=COEF_X(8)*AKG
      J1=2*L1F
      J2=2
      M1=MII-MIF
      M2=-2
      J=J1F
      M=-2+MII-MIF
      CALL KLEGO
      COEF_X(9)=COEF_X(9)*AKG
      J1=J1F
      J2=JO(NU1) 
      M1=-2+MII-MIF
      M2=MIF
      J=JS1
      M=-2+MII
      CALL KLEGO
      COEF_X(7)=COEF_X(7)*AKG
      COEF_X(8)=COEF_X(8)*AKG 
      COEF_X(9)=COEF_X(9)*AKG 

      CORR=DCOS(COPH(NU1,L1F+1)+COPH(1,L1I+1)) 
	COII=DSIN(COPH(NU1,L1F+1)+COPH(1,L1I+1)) 
	CSRX=CORR*CRD(K1,N1C,N1R)-COII*CID(K1,N1C,N1R)
	CSIX=CORR*CID(K1,N1C,N1R)+COII*CRD(K1,N1C,N1R)

      DO 616 M=1,MTET
	IF(NU1.EQ.1 .AND. L1I.EQ.0 .AND. J1I.EQ.2   ! .AND. JS1.EQ.(IC+2) 
     * .AND. L1F.EQ.0	.AND. J1F.EQ.2 .AND. MII.EQ.MIF)
     *THEN !add the Coulomb amplitude(only added once for different J1I and JS1)
	X_R(1,NU1,MII_INDEX,MIF_INDEX,M)=X_R(1,NU1,MII_INDEX,MIF_INDEX,M)+
     *FC_R(M)
	X_I(1,NU1,MII_INDEX,MIF_INDEX,M)=X_I(1,NU1,MII_INDEX,MIF_INDEX,M)+
     *FC_I(M)
	X_R(5,NU1,MII_INDEX,MIF_INDEX,M)=X_R(5,NU1,MII_INDEX,MIF_INDEX,M)+
     *FC_R(M)
	X_I(5,NU1,MII_INDEX,MIF_INDEX,M)=X_I(5,NU1,MII_INDEX,MIF_INDEX,M)+
     *FC_I(M)
     	X_R(9,NU1,MII_INDEX,MIF_INDEX,M)=X_R(9,NU1,MII_INDEX,MIF_INDEX,M)+
     *FC_R(M)
	X_I(9,NU1,MII_INDEX,MIF_INDEX,M)=X_I(9,NU1,MII_INDEX,MIF_INDEX,M)+
     *FC_I(M)
	COUNT_FC=COUNT_FC+1
	ENDIF
      X_R(1,NU1,MII_INDEX,MIF_INDEX,M)=X_R(1,NU1,MII_INDEX,MIF_INDEX,M)+
     *CSRX*YLM_VAL(L1F,(MII-MIF)/2,M)*COEF_X(1)
      X_I(1,NU1,MII_INDEX,MIF_INDEX,M)=X_I(1,NU1,MII_INDEX,MIF_INDEX,M)+
     *CSIX*YLM_VAL(L1F,(MII-MIF)/2,M)*COEF_X(1)
      X_R(2,NU1,MII_INDEX,MIF_INDEX,M)=X_R(2,NU1,MII_INDEX,MIF_INDEX,M)+
     *CSRX*YLM_VAL(L1F,1+(MII-MIF)/2,M)*COEF_X(2)
      X_I(2,NU1,MII_INDEX,MIF_INDEX,M)=X_I(2,NU1,MII_INDEX,MIF_INDEX,M)+
     *CSIX*YLM_VAL(L1F,1+(MII-MIF)/2,M)*COEF_X(2)
	X_R(3,NU1,MII_INDEX,MIF_INDEX,M)=X_R(3,NU1,MII_INDEX,MIF_INDEX,M)+
     *CSRX*YLM_VAL(L1F,2+(MII-MIF)/2,M)*COEF_X(3)
      X_I(3,NU1,MII_INDEX,MIF_INDEX,M)=X_I(3,NU1,MII_INDEX,MIF_INDEX,M)+
     *CSIX*YLM_VAL(L1F,2+(MII-MIF)/2,M)*COEF_X(3)
     
      X_R(4,NU1,MII_INDEX,MIF_INDEX,M)=X_R(4,NU1,MII_INDEX,MIF_INDEX,M)+
     *CSRX*YLM_VAL(L1F,-1+(MII-MIF)/2,M)*COEF_X(4)
      X_I(4,NU1,MII_INDEX,MIF_INDEX,M)=X_I(4,NU1,MII_INDEX,MIF_INDEX,M)+
     *CSIX*YLM_VAL(L1F,-1+(MII-MIF)/2,M)*COEF_X(4)
      X_R(5,NU1,MII_INDEX,MIF_INDEX,M)=X_R(5,NU1,MII_INDEX,MIF_INDEX,M)+
     *CSRX*YLM_VAL(L1F,(MII-MIF)/2,M)*COEF_X(5)
      X_I(5,NU1,MII_INDEX,MIF_INDEX,M)=X_I(5,NU1,MII_INDEX,MIF_INDEX,M)+
     *CSIX*YLM_VAL(L1F,(MII-MIF)/2,M)*COEF_X(5)   
      X_R(6,NU1,MII_INDEX,MIF_INDEX,M)=X_R(6,NU1,MII_INDEX,MIF_INDEX,M)+
     *CSRX*YLM_VAL(L1F,1+(MII-MIF)/2,M)*COEF_X(6)
      X_I(6,NU1,MII_INDEX,MIF_INDEX,M)=X_I(6,NU1,MII_INDEX,MIF_INDEX,M)+
     *CSIX*YLM_VAL(L1F,1+(MII-MIF)/2,M)*COEF_X(6)      

      X_R(7,NU1,MII_INDEX,MIF_INDEX,M)=X_R(7,NU1,MII_INDEX,MIF_INDEX,M)+
     *CSRX*YLM_VAL(L1F,-2+(MII-MIF)/2,M)*COEF_X(7)
      X_I(7,NU1,MII_INDEX,MIF_INDEX,M)=X_I(7,NU1,MII_INDEX,MIF_INDEX,M)+
     *CSIX*YLM_VAL(L1F,-2+(MII-MIF)/2,M)*COEF_X(7)
      X_R(8,NU1,MII_INDEX,MIF_INDEX,M)=X_R(8,NU1,MII_INDEX,MIF_INDEX,M)+
     *CSRX*YLM_VAL(L1F,-1+(MII-MIF)/2,M)*COEF_X(8)
      X_I(8,NU1,MII_INDEX,MIF_INDEX,M)=X_I(8,NU1,MII_INDEX,MIF_INDEX,M)+
     *CSIX*YLM_VAL(L1F,-1+(MII-MIF)/2,M)*COEF_X(8)   
      X_R(9,NU1,MII_INDEX,MIF_INDEX,M)=X_R(9,NU1,MII_INDEX,MIF_INDEX,M)+
     *CSRX*YLM_VAL(L1F,(MII-MIF)/2,M)*COEF_X(9)
      X_I(9,NU1,MII_INDEX,MIF_INDEX,M)=X_I(9,NU1,MII_INDEX,MIF_INDEX,M)+
     *CSIX*YLM_VAL(L1F,(MII-MIF)/2,M)*COEF_X(9)      
  616 CONTINUE
  615 CONTINUE  
  614 CONTINUE
    8 CONTINUE
    1 CONTINUE
    
      DO 619 N=NUI,NUF
	MII_MAX=IC+1
	DO 612 MII_INDEX=1,MII_MAX
      MIF_MAX=JO(N)+1
      DO 618 MIF_INDEX=1,MIF_MAX
      DO 617 M=1,MTET
       CS_ANG_MIF(N,MII_INDEX,MIF_INDEX,M)=
     *X_R(1,N,MII_INDEX,MIF_INDEX,M)**2+ 
     *X_I(1,N,MII_INDEX,MIF_INDEX,M)**2+
     *X_R(2,N,MII_INDEX,MIF_INDEX,M)**2+
     *X_I(2,N,MII_INDEX,MIF_INDEX,M)**2+
     *X_R(3,N,MII_INDEX,MIF_INDEX,M)**2+ 
     *X_I(3,N,MII_INDEX,MIF_INDEX,M)**2+
     *X_R(4,N,MII_INDEX,MIF_INDEX,M)**2+
     *X_I(4,N,MII_INDEX,MIF_INDEX,M)**2+
     *X_R(5,N,MII_INDEX,MIF_INDEX,M)**2+
     *X_I(5,N,MII_INDEX,MIF_INDEX,M)**2+
     *X_R(6,N,MII_INDEX,MIF_INDEX,M)**2+
     *X_I(6,N,MII_INDEX,MIF_INDEX,M)**2+     
     *X_R(7,N,MII_INDEX,MIF_INDEX,M)**2+
     *X_I(7,N,MII_INDEX,MIF_INDEX,M)**2+ 
     *X_R(8,N,MII_INDEX,MIF_INDEX,M)**2+
     *X_I(8,N,MII_INDEX,MIF_INDEX,M)**2+
     *X_R(9,N,MII_INDEX,MIF_INDEX,M)**2+
     *X_I(9,N,MII_INDEX,MIF_INDEX,M)**2     
      CS_ANG(N,M)=CS_ANG(N,M)+CS_ANG_MIF(N,MII_INDEX,MIF_INDEX,M)
     */MII_MAX/3.D0 !calculate the cross-section
  617 CONTINUE
  618 CONTINUE
  612 CONTINUE
      DO 620 M=1,MTET 
      DISC(N,M)=CS_ANG(N,M)
  620 CONTINUE    
  619 CONTINUE 


      RETURN
      END SUBROUTINE DEUTR
C     *******************************************************



C     END of ccrd
C     *******************************************************
