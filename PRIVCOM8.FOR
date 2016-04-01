C     *******************************************************
C     START of privcom8  
C     *******************************************************
C==========================================================================
C     PRIVATE COMMONS FOR EACH PARALLEL THREAD TO BE USED AS INCLUDE BLOCK


C     FROM SUBROUTINE SHEMM 
      COMMON/SHEMM/ES(40),JU(40),NTU(40),NNB(40),NNG(40),NNO(40),NPI(40)  

!$OMP THREADPRIVATE(/SHEMM/)

C     *******************************************************
C     END of privcom8
C     *******************************************************
