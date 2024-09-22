C     *******************************************************
C     START of privcom1  
C     *******************************************************
C==========================================================================
C     PRIVATE COMMONS FOR EACH PARALLEL THREAD TO BE USED AS INCLUDE BLOCK


C     FROM SUBROUTINE RACAH  
      COMMON/RACB/JA,JB,JC,JD,JE,JF,W

C     FROM SUBROUTINE PLEGA       
      COMMON/DISK/TET(150),MTET
      COMMON/DISCAN/DISC(40,150),PLEG(180,150),COEF(40,180),LKK
      
      COMMON/KGxx/J1,J2,M1,M2,J,M,AKG
      
C     FROM SUBROUTINE ANPOW              
      COMMON/POL1/Px(40,150),Pz(40,150),Ay(40,150)

C     FROM SUBROUTINE PLM_CAL              
 	COMMON/PLM/PLM_DATA(181,181,150)

CCC   /DISK/ is input      
!$OMP THREADPRIVATE(/DISK/,/PLM/)
!$OMP THREADPRIVATE(/RACB/,/DISCAN/,/KGxx/,/POL1/)



C     *******************************************************
C     END of privcom1
C     *******************************************************
