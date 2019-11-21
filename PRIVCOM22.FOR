C     *******************************************************
C     START of privcom22  
C     *******************************************************
C==========================================================================
C     PRIVATE COMMONS FOR EACH PARALLEL THREAD TO BE USED AS INCLUDE BLOCK
C     Nucleus number

      COMMON/NIND/IIS
!$OMP THREADPRIVATE(/NIND/)          
     

C     *******************************************************
C     END of pricom22
C     *******************************************************
