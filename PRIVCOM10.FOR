C     *******************************************************
C     START of privcom10  
C     *******************************************************
C==========================================================================
C     PRIVATE COMMONS FOR EACH PARALLEL THREAD TO BE USED AS INCLUDE BLOCK


C     FROM SUBROUTINE KNCOE   
      COMMON/MENU/MEJOB,MEPOT,MEHAM,MECHA,MEPRI,MESOL,MESHA,MESHO,MEHAO
     *,MEAPP,MEVOL,MEREL,MECUL,MERZZ,MERRR,MEDIS,MERIP,MEDEF,MEAXI,MERAD

!$OMP THREADPRIVATE(/MENU/)
     
C     *******************************************************
C     END of privcom10
C     *******************************************************
