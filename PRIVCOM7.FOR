C     *******************************************************
C     START of privcom7  
C     *******************************************************
C==========================================================================
C     PRIVATE COMMONS FOR EACH PARALLEL THREAD TO BE USED AS INCLUDE BLOCK


C     FROM SUBROURINE LOGMO
      COMMON/LNP/JSS,NPIS,NPIO,JSO,NN1,LNO(180),NS1(180),JNO(180),NSPI
    
     
      
!$OMP THREADPRIVATE(/LNP/)


      

C     *******************************************************
C     END of pricom7
C     *******************************************************
