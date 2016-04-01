C     *******************************************************
C     START of privcom15  
C     *******************************************************
C==========================================================================
C     PRIVATE COMMONS FOR EACH PARALLEL THREAD TO BE USED AS INCLUDE BLOCK

      COMMON /ENB/EE(50),MCHAE(50)

      COMMON/INP0/ELC(40),BETBC(40),JOC(40),NP0C(40),KOC(40),NUMBC(40),
     *NCAC(40),NTUC(40),NNBC(40),NNGC(40),NNOC(40),NPOC(40)

      COMMON/INP1/CAVR,CARR,CAAR,CARD,CAAC,ATI

!$OMP THREADPRIVATE(/ENB/,/INP0/,/INP1/)

C     *******************************************************
C     END of privcom15
C     *******************************************************
