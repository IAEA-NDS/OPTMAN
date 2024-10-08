C      VERSION 24 (Sep 2024, PARALLEL CALCULATIONS USING OPENMP
C      INSTRUCTIONS, LANE CONSISTENT COULOMB CORRECTION 
C      WITH DISPERSIVE OPTICAL POTENTIAL RELATIONS, NON-AXIAL AND
C      NEGATIVE-PARITY-BAND LEVEL COUPLING, EVEN-EVEN and ODD
C      NUCLIDES CASE, ANALYZING POWERS, GLOBAL SEARCH, NATURAL 
C      FOR TOTAL AND ISOTOPE RATIO ADJUSTMENT, LOW ENERGY RESONANCES)
C      EFFECTIVE DEFORMATIONS IN EVEN-EVEN NUCLIDE MAY BE ADJUSTABLE
C      OR FROM SRM (IN THIS CASE ADDITIONAL TO RIGID ROTATOR COUPLING
C      COMING FROM DYNAMICS IS ACCOUNTED, OPTIONALY COUPLING COMMING FROM
c      NUCLEAE VOLUME CONSERVATION CAN BE ACCOUNTED)
C
C      ENABLING "OPENMP" OPTION IN COMPILER WILL PRODUCE MULTITHREADED CODE
C      DEFAULT (WITHOUT OPENMP FLAG) COMPILATION WILL PRODUCE SINGLETHREADED CODE
C      
C      TO USE LAPACK'S MATRIX INVERSION (FASTER, USING REAL*8): 
C      1) SET "DLAPACK" OPTION FOR THE COMPILER AND
C      2) ADD LAPACK LIBRARY PATH TO COMPILER AND LINKER OR SET OTHER FLAG FOR 
C         USING COMPILER'S LAPACK LIB (E.G. "MKL" FOR INTEL FORTRAN)
C      DEFAULT (WITHOUT "DLAPACK") COMPILATION WILL USE MATRIX INVERSION FROM INMAT
C      SUBROUTINE (USING REAL*16)
C
C      ENABLED COMPILER PREPOCESSING (FPP OR CPP) IS ALWAYS REQUIRED 
C     
C      FOR EMPIRE COMPILE WITH FLAG "/DEMPMODE" 
C      FOR FORMATTED PRINT OF *.TLJ COMPILE WITH FLAG "/DFORMATTEDOUT"
C      PLEASE USE THESE FLAGS INSTEAD OF CHANGING SOURCES      
C
C      FULL DOUBLE PRECISION (REAL*8 = DOUBLE PRECISION, COMPLEX*16 = DOUBLE COMPLEX)
C      To allow automatic extension to quadruple precision (REAL*16, COMPLEX*32)
C
C      AUTHORS:
C
C      EFREM SOUKHOVITSKI - NON-AXIAL SOFT-ROTATOR NUCLEAR MODEL THEORY, MAIN
C                           CC COMPUTATIONAL ALGORITHMS AND CODING, LANE 
C                           CONSISTENCY AND (P,N) CS, PARALLELIZATION
C                           E-MAIL: esukhov@sosny.bas-net.by
C
C      SATOSHI CHIBA      - THEORY DEVELOPMENT in COOPERATION WITH
C                           E. SOUKHOVITSKI, E-MAIL: chiba.satoshi@jaea.go.jp
C
C      ROBERTO CAPOTE     - DISPERSIVE OPTICAL MODEL POTENTIAL RELATIONS,
C                           LANE CONSISTENCY AND (P,N) CS, RIPL INTERFACE,
C                           PARALLELIZATION  
C                           E-MAIL: R.CapoteNoy@iaea.org
C
C      JOSE M. QUESADA    - DISPERSIVE OPTICAL MODEL POTENTIAL RELATIONS,
C                           LANE CONSISTENCY,E-MAIL: quesada@us.es
C      
C      DMITRY MARTYANOV   - THEORY DEVELOPMENT in COOPERATION WITH
C                           E. SOUKHOVITSKI, CODE DEVELOPEMENT
C                           E-MAIL: dmart@sosny.bas-net.by
C     
C
C      MAIN REFERENCES:   1. E.SH. SOUKHOVITSKII, S. CHIBA, R. CAPOTE, JOSE M.
C                            QUESADA, S. KUNIEDA and G.B. MOROGOVSKII, TECHNICAL 
C                            REPORT,JAEA-DATA/CODE 2008-025, JAPAN ATOMIC ENERGY
C                            AGENCY, 2008.
C                         2. E.SH. SOUKHOVITSKII, S. CHIBA, O.IWAMOTO, K.SHIBATA
C                            T. FUKAHORI and G.B. MOROGOVSKII, TECHNICAL REPORT,
C                            JAERI-DATA/CODE 2005-002, JAPAN ATOMIC ENERGY
C                            INSTITUTE, 2005.
C                         3. E.SH. SOUKHOVITSKII, R. CAPOTE, J.M. QUESADA,
C                            S. CHIBA, PHYS. REV. C72, 024604 (2005).
C                         4. J.M. QUESADA, R. CAPOTE, E.Sh. SOUKHOVITSKI and
C                            S. CHIBA, PHYS. REV. C76, 057602 (2007).
C                         5. MORE DETAILS, MANUAL,CODE'S SOURCE FILES AND INPUTS
C                            WITH VARIOUS CC OPTICAL POTENTIALS CAN BE FOUND ON: 
C                            http://www-nds.iaea.org/RIPL-3/
C
C
C  ****************************************************************
      PROGRAM OPTMAN12
C     SUBROUTINE OPTMAN12(fname)
C  *************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      logical f_ex
      character*20 answ
      INCLUDE 'PRIVCOM20.FOR'

      INCLUDE 'PRIVCOM10.FOR'
      INCLUDE 'PRIVCOM17D.FOR'  
      INCLUDE 'PRIVCOM21.FOR'
      
      CHARACTER*80 TITLE

      REAL*16 LFA(400),dtmp  
C     REAL*8  LFA(400),dtmp  

      DOUBLE PRECISION A
      INCLUDE 'PRIVCOM4.FOR'
c      COMMON/LOFAC/A(800)

      INTEGER NTHREADS, TID, narg
!$    INTEGER OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
C     INTEGER omp_set_num_threads

C----------------------------------------------------------------------------
C----------------------------------------------------------------------------
C     FACTORIAL CALCULATION AVOIDING A LONG DATA STATEMENT (common /LOFAC/A)
C            but keeping the same precision
C   (a long data statement was producing errors/warnings with some compilers)
C
      A = 0.D0
      LFA(1) = 0
      DO i = 2, 400
        dtmp = i 
        LFA(i) = LOG(dtmp) + LFA(i - 1)
      ENDDO
      DO j = 6, 800
        if(mod(j,2).eq.0) A(j)=LFA(j/2 - 1)  
C       if(j.lt.15) write(*,*) j,A(j) 
      ENDDO
C
C     FORMER BLOCK DATA VALUES
C     DATA A/0.d0,0.d0,0.d0,0.d0,
C    *0.,.693147180559945309417232D+00,0.,.179175946922805500081248D+01,
C    *0.,.317805383034794561964694D+01,0.,.478749174278204599424770D+01,
C    *0.,.657925121201010099506018D+01,0.,.852516136106541430016553D+01,
C    *0.,.106046029027452502284172D+02,0.,.128018274800814696112077D+02,

C----------------------------------------------------------------------------
C----------------------------------------------------------------------------
C
C     Logical variable EMPIRE 
C     EMPIRE = .TRUE.  -> OPTMAN used within the EMPIRE system
C     EMPIRE = .FALSE. -> OPTMAN used in stand-alone mode
C
#ifdef EMPMODE
      EMPIRE = .true.  
#else
      EMPIRE = .false.  
#endif
      IF (EMPIRE) THEN 
C--------------------- EMPIRE related i/o changes ----------------------
C       Input filename fixed to OPTMAN.INP for EMPIRE
        open(unit=20,file='OPTMAN.INP',STATUS='OLD')
C       Output root filename fixed to ecis06 for EMPIRE
        fname='ecis06'  
C       Output filename fixed to OPTMAN.OUT for EMPIRE
        open(unit=21,file='OPTMAN.OUT')
        WRITE(21,'(5x,A)')
     *  '***********************************************'
        WRITE(21,'(5x,A)')
     *  '*      CODE OPTMAN VERSION 24 ( SEP 2024)     *'
!$      WRITE(21,'(5x,A)')
!$   *  '*      OPENMP version for parallel execution  *'
#ifdef LAPACK      
        WRITE(21,'(5x,A)')
     *  '*   USE LAPACK LIBRARY FOR MATRIX INVERSION   *'
#endif        
        WRITE(21,'(5x,A)')
     *  '*                                             *'
        WRITE(21,'(5x,A)')
     *  '*  DISPERSIVE RELATIONS AND LANE CONSISTENCY  *'
        WRITE(21,'(5x,A)')
     *  '*      LANE CONSISTENT COULOMB CORRECTION     *'
        WRITE(21,'(5x,A)')
     *  '*            GLOBAL POTENTIAL SEARCH          *'
        WRITE(21,'(5x,A)')
     *  '*  OTHER NON-AXIAL BANDS LEVELS COUPLING      *'
         WRITE(21,'(5x,A)')
     *  '*     OPTION USING AXIAL ROTATIONAL MODEL     *'
        WRITE(21,'(5x,A)')
     *  '*    POTENTIAL MULTIPOLES, ANALYZING POWERS   *'
        WRITE(21,'(5x,A)')
     *  '*---------------------------------------------*'
        WRITE(21,'(5x,A)')
     *  '*    COMPATIBLE WITH THE EMPIRE-3.2 SYSTEM    *'
        WRITE(21,'(5x,A)')
     *  '***********************************************'

      ELSE
C--------------------- FOR NORMAL OPERATION (NOT EMPIRE) ---------------
        WRITE(*,'(A)')' ***********************************************'
        WRITE(*,'(A)')' *      CODE OPTMAN VERSION 24 ( SEP 2024)     *'
!$      WRITE(*,'(A)')' *      OPENMP version for parallel execution  *'
#ifdef LAPACK      
        WRITE(*,'(A)')' *   USE LAPACK LIBRARY FOR MATRIX INVERSION   *'
#endif         
        WRITE(*,'(A)')' *                                             *'
        WRITE(*,'(A)')' *  DISPERSIVE RELATIONS AND LANE CONSISTENCY  *'
        WRITE(*,'(A)')' *      LANE CONSISTENT COULOMB CORRECTION     *'
        WRITE(*,'(A)')' *            GLOBAL POTENTIAL SEARCH          *'
        WRITE(*,'(A)')' *     OTHER NON-AXIAL BANDS LEVELS COUPLING   *'
        WRITE(*,'(A)')' *    OPTION USING AXIAL ROTATIONAL MODEL      *'
        WRITE(*,'(A)')' *   POTENTIAL MULTIPOLES, ANALYZING POWERS    *'
        WRITE(*,'(A)')' ***********************************************'
C
        narg = 1
        CALL getarg(narg,fname)
        if(trim(fname(1:1)).eq.'') then
          WRITE(*,'(1X,A40)') 'INPUT FILE NAME (without extension) ? =>'
          READ(*,'(A20)') fname
        endif
        if(fname(1:1).eq.'') fname='OPTMAN'

C       fname='rigidtest'
C       fname='globalFe'
C       fname='TESTCAL'

        inquire(file=TRIM(fname)//'.INP',exist=f_ex)
        if(f_ex) then
          open(unit=20,file=TRIM(fname)//'.INP',STATUS='OLD')
        else 
          print *,'Input file not found.'
          stop
        end if 
C
C        open(unit=21,file=TRIM(fname)//'.OUT',STATUS='NEW')
        inquire(file=TRIM(fname)//'.OUT',exist=f_ex)
        if(f_ex) then
          print *,'Output file already exists. Enter [yes] to'//
     *            ' overwrite or any other to cancel.'
          read(*,'(A20)') answ
          if(trim(answ).ne.'yes') stop
        end if 
        
        open(unit=21,file=TRIM(fname)//'.OUT')

        WRITE(21,'(5x,A)')
     *  '***********************************************'
        WRITE(21,'(5x,A)')
     *  '*     CODE OPTMAN VERSION 24 ( SEP 2024)      *'
!$      WRITE(21,'(5x,A)')
!$   *  '*    OPENMP version for parallel execution    *'
#ifdef LAPACK      
        WRITE(21,'(5x,A)')
     *  '*   USE LAPACK LIBRARY FOR MATRIX INVERSION   *'
#endif                
        WRITE(21,'(5x,A)')
     *  '*                                             *'
        WRITE(21,'(5x,A)')
     *  '*  DISPERSIVE RELATIONS AND LANE CONSISTENCY  *'
        WRITE(21,'(5x,A)')
     *  '*      LANE CONSISTENT COULOMB CORRECTION     *'
        WRITE(21,'(5x,A)')
     *  '*            GLOBAL POTENTIAL SEARCH          *'
        WRITE(21,'(5x,A)')
     *  '*  OTHER NON-AXIAL BANDS LEVELS COUPLING      *'
         WRITE(21,'(5x,A)')
     *  '*     OPTION USING AXIAL ROTATIONAL MODEL     *'
        WRITE(21,'(5x,A)')
     *  '*    POTENTIAL MULTIPOLES, ANALYZING POWERS   *'
        WRITE(21,'(5x,A)')
     *  '*---------------------------------------------*'
        WRITE(21,'(5x,A)')
     *  '*    COMPATIBLE WITH THE EMPIRE-3.2 SYSTEM    *'
        WRITE(21,'(5x,A)')
     *  '***********************************************'
      
      ENDIF

      CALL THORA(21)

            READ (20,4) TITLE
            WRITE(21,'(7x,A)') trim(TITLE)
            READ(20,1)MEJOB,MEPOT,MEHAM,MEPRI,MESOL,MESHA,MESHO,
     *                MEHAO,MEAPP,MEVOL,MEREL,MECUL,MERZZ,MERRR,
     *                MEDIS,MERIP,MEDEF,MEAXI,MERAD

C     FOR EMPIRE OPERATION MEPRI IS SET TO 98
      IF(EMPIRE) MEPRI=98    
C     TO MINIMIZE PRINTING DURING PARALLEL CALCULATIONS SET MEPRI to 99

    1 FORMAT(20I2)
    4 FORMAT(A80)
C#1   MEJOP 1-JOB STANDARD*2-JOB WITH POTENTIAL OPTIMIZATION.
C#2   MEPOT 1-POT-AL OF ROT.MODEL  YL0* 2-POT-AL EXPANDED BY BETTA
C#3   MEHAM 1-RM* 2-VM* 3-DCHM* 6-5PARM* 4-FDM* 5-5PAR0M* 7-COUPL.GB
C#4   MEPRI OUTPUT MANAGEMENT * 0-MINIMUM (98 = EMPIRE I.E. NO SCREEN OUTPUT AND LIMITED FILE OUTPUT
C           OUTPUT MANAGEMENT * 98- FOR EMPIRE I.E. NO SCREEN OUTPUT AND LIMITED FILE OUTPUT
C           OUTPUT MANAGEMENT * 99- FOR PARALLEL CALCS.I.E. NO SCREEN OUTPUT AND LIMITED FILE OUTPUT
C#5   MESOL 1-OPTIMIZED 2-EXACT >3-ITERATION METHODS
C           3-ZERO APPROX..-SPHERICAL OPTICAL MODEL
C          >3-ZERO APPROX. HAS THIS NUMBER OF COUPLED EQ.
C#6   MESHA 1-QUADR*2-+HEXAD.AXIAL.* 3-+HEXAD.N.AXIAL. DEF BY  GAM.
C     *4-COMMON CASE.
C#7   MESHO 0-NO *1-AXIAL.*2-NON-AXIAL OCTUPOLE DEFORMATION
C#8   MEHAO 0-NO *1-CONSIDER. OF OCT OSC. * 2-SIMMET. OCTUPOLE OSC. ScaLED BY \beta^2
C           3-2-SIMMET. OCTUPOLE OSC. NOT SCALED BY \beta^2
C#9   MEAPP 0-EXECT SOLUTION; *1-QUICK SOLUTION WITHOUT LEVEL'S  POTENTIAL DEPENDANCE
C           2- The most common case
C#10  MEVOL 0-STANDARD SOLUTION; *=1- ACCOUNT OF NUCLEAR VOLUME; *=2 - ALSO ACCOUNT 
C           CENTER OF MASS IMMOBILITY
C     CONSERVATION
C#11  MEREL 0-STANDARD SOLUTION;
C           1-ACCOUNT OF RELATIVISTIC KINEM AND POTEN. DEPENDENCE
C           2-ACCOUNT OF RELATIVISTIC KINEMATICS
C           3-ACCOUNT OF RELATIVISTIC KINEM AND REAL POTEN. DEPENDENCE
C#12  MECUL 0-COULOMB CORRECTION PROPORTIONAL TO DERIVATIVE OF REAL PORENTIAL
C           1-CONSTANT
C           2-LANE CONSISTANT COULOMB CORRECTION NUCLEAR POTENTIAL ENERGY FOR PROTONS EQUAL TO
C           INCIDENT ENERGY - CDE, APPLIED TO BOTH REAL AND UMAGINARY POTENTIAL
C           3-LANE CONSISTANT COULOMB CORRECTION NUCLEAR POTENTIAL ENERGY FOR PROTONS EQUAL TO
C           INCIDENT ENERGY - CDE, APPLIED TO REAL POTENTIAL ONLY
C#13  MERZZ 0-CHARGE RADIUS -CONSTANT
C           1-ENERGY DEPENDENT
C#14  MERRR 0-REAL RADIUS IS ENERGY INDEPENDENT
C           1-REAL RADIUS IS ENERGY DEPENDENT
C#15  MEDIS 0-WITHOUT ACCOUNT OF DISPERSION RELATIONS BETWEEN REAL AND IMAGINARY POTENTIALS
C           1-ACCOUNT OF DISPERSION RELATIONS BETWEEN REAL AND IMAGINARY POTENTIALS
C           2-ACCOUNT OF DISPERSION RELATIONS BETWEEN REAL AND IMAGINARY POTENTIALS,
C           EXCLUDING SPIN ORBIT POTENTIAL
C#16  MERIP 0-IN ABCT READS POTENTIAL ONE TIME WITH ANALITICAL DEPENDENCIES FOR ALL ENERGIES
C           1-ABCT READS POTENTIAL BLOCKS FOR EACH ENERGY TO USE RIPL COMPILED INPUTS
C#17  MEDEF 0-DEFORMATIONS OF NON-GS-BAND LEVELS ARE INPUTTED AND CAN BE ADJUSTED
C           1-DEFORMATIONS OF NON-GS-BAND LEVELS ARE CALCULATED FROM SOFT-ROTATOR MODEL
C           2-DEFORMATIONS OF GS-BAND LEVELS AND NON-AXIAL BAND ARE FROM RIGID-ROTATOR MODEL      
C@18  MEAXI 0-USES AXIAL WEIGHTS (NOT AFFECT "EFFECTIVE" DEFORMATIONS)
C           1-USES NON-AXIAL WEIGHTS
C#19   MERAD 0-TREAT R_i as CONSTANT, BUT IF MEVOL>=1 ADD STATIC TERM (BTGS2) TO VOLUME CONSERVATION  //+B_30^2 not ready yet
C            1-R_i CORRECTION DUE TO STATIC DEFORMATIONS (CONSTANT NUCLEAR DENSITY)      
C            2-NO RADIUS CORRECTION OF ANY KIND      
C
c      NTHREADS = omp_set_num_threads(16)
!$OMP PARALLEL PRIVATE(TID) 
      TID = 0
      NTHREADS = 1
!$    TID = OMP_GET_THREAD_NUM()
      IF(TID.eq.0.and.MEPRI.NE.98) then
!$      NTHREADS = OMP_GET_NUM_THREADS()
        PRINT *
        PRINT *, 'Number of threads =', NTHREADS
        PRINT *
      ENDIF
!$OMP END PARALLEL

      IF(MEJOB.NE.2) THEN
C
C        OMP CALCULATION
C
        IF(.NOT.EMPIRE) then
C         TRANSME ARRAY FOR TRANSITIONS TL
          open(unit=22,file='TRANSME')
C         GNASH ARRAY FOR TRANSITIONS TLJ, J=L+ - 1/2
          open(unit=24,file='GNASH')
          open(unit=23,file='CR-SECT')
          open(unit=25,file='ANG-DIST')
          open(unit=26,file='ANG-POL')
          open(unit=27,file='ANGDIS-yw')
          open(unit=327,file='ANAPOW-yw')
          open(unit=328,file='POLARIZATION-yw')
        ENDIF

        CALL ABCT

        IF(.NOT.EMPIRE) then
          close(22)
          close(24)
          close(23)
          close(25)
          close(26)
          close(27)
          close(327)
          close(328)
        ENDIF 

      ELSE
C
C        OMP FITTING
C
         CALL DATET

      ENDIF

      CLOSE(20)
      CLOSE(21)

C
C     RETURN
C
      END
C     *******************************************************
      SUBROUTINE ABCT
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      CHARACTER*4 ITEXT,IMOD,IPOT,MPOT,IFOR1,IFOR2,
     *IFOR3,ISMOD,IFOR4,IHMOD,IFOR5 
     
      
      INCLUDE 'PRIVCOM21.FOR'
      INCLUDE 'PRIVCOM4.FOR'
      
      INCLUDE 'PRIVCOM20.FOR'

C     Data read in OPTMAND
      INCLUDE 'PRIVCOM10.FOR'
C     FROM SUBROUTINE KNCOE   
C     COMMON/MENU/MEJOB,MEPOT,MEHAM,MECHA,MEPRI,MESOL,MESHA,MESHO,MEHAO
C    *,MEAPP,MEVOL,MEREL,MECUL,MERZZ,MERRR,MEDIS,MERIP,MEDEF,MEAXI

C
      INCLUDE 'PRIVCOM12.FOR'
C     FROM SUBROUTINE OVLAG
C     COMMON/SHEM1/HW,AMG0,AMB0,GAM0,BET0,BET4,GAMDE,GSHAPE 
C
      INCLUDE 'PRIVCOM9.FOR'
C     FROM SUBROUTINE KNCOE   
C     COMMON/SHAMO/BET3,ETO,AMUO,HWO,BB32,DPAR
C     FROM SUBROUTINE KNDIT      
C     COMMON/INRMi/BB42,GAMG,DELG
C     FROM SUBROUTINE KNDIT      
C     COMMON/ALFAi/AGSIC(40) 
C
      INCLUDE 'PRIVCOM13.FOR'
C     FROM SUBROUTINE OVLAG
C     COMMON/ENAi/BET(10),NUR,NPD,LAS 
C
      INCLUDE 'PRIVCOM1.FOR'
C     FROM SUBROUTINE PLEGA       
C     COMMON/DISK/TET(150),MTET
C
      INCLUDE 'PRIVCOM.FOR'         
C      COMMON/RAD/RR,RC,RD,RW,RS,AR,AC,AW,AD,AS,RZ
C     (a)
C     FROM SUBROUTINE QUANT
C      COMMON/RESONI/ERN(10),GNN(10),GREN(10),LON(10),JMN(10),JCON(10),
C     *NEL(10),NRESN
C     FROM SUBROUTINE RIPAT 
C      COMMON/POT2/AR1,AC1,AW1,AD1,AS1
C     */POT3/BNDC,WDA1,WCA1,CCOUL,CISO,WCISO,WS0,WS1
C      COMMON/DISPE2/VRD,WDSHI,WDWID2,ALFNEW 
C     */POTD/ALAVR,VRLA,WCBW,WCWID,WDBW,WDWID,ALAWD,EFERMN,EFERMP,ALASO
C     *,PDIS,WSBW,WSWID,RRBWC,RRWID,RZBWC,RZWID
C     (b)
C     FROM SUBROUTINE POTET
C      COMMON/RADi/ALF,AT,ANEU,ZNUC,ASP,AZ
C     FROM SUBROUTINE POTET
C      COMMON/DISPEi/EA,WDISO
C     FROM SUBROUTINE QUANT
C      COMMON/CSBi/NST
C     FROM SUBROUTINE ASFUT
C      COMMON/NCLMAi/LLMA,NCMA,NSMA
C     */NCLMA/KODMA

      INCLUDE 'PRIVCOM15.FOR'
C     COMMON /ENB/EE(50),MCHAE(50)
C     COMMON/INP0/ELC(40),BETBC(40),JOC(40),NP0C(40),KOC(40),NUMBC(40),
C    *NCAC(40),NTUC(40),NNBC(40),NNGC(40),NNOC(40),NPOC(40)
C     COMMON/INP1/CAVR,CARR,CAAR,CARD,CAAC,ATI
C
      INCLUDE 'PRIVCOM6.FOR'
C     FROM SUBROUTINE RIPAT          
C     COMMON/POT1/VR3,VR2,VR1,VR0,WC1,WC0,WD1,WD0,VS,AC0,AR0,AW0,AD0,AS0
C                                                    .        .
      INCLUDE 'PRIVCOM14.FOR'
      INCLUDE 'PRIVCOM16D.FOR'

C
      INTEGER TID, IIparal

      DIMENSION ITEXT(3),IMOD(7),IPOT(7),MPOT(2),IFOR1(2),IFOR2(9),
     *IFOR3(3),ISMOD(2),IFOR4(5),IHMOD(3),IFOR5(4)

      DATA ITEXT,IMOD/4HHAMI,4HLTON,4HIAN ,4H RV ,4H VM ,4HDCHM,
     *4H FDM,4H5PA0,4H 5PM,4HCLGB/,
     *IPOT,MPOT/4HPOTE,4HNTIA,4HL EX,4HPAND,4HED  ,4HBY  ,
     *4H    ,4HYL0 ,4HBET0/,IFOR1,IFOR2,IFOR3/4HWITH,4H AC.,
     *4HAXIA,4HL HE,4HXADE,4HCAPO,4HLE D,4HEFOR,4HMATI,4HONS ,
     *4H    ,4H    ,4H|NON,4H NON/,
     *ISMOD,IFOR4,IHMOD,IFOR5/4H    ,4H NON,4HAXIA,4HL OC,4HTUP0,
     *4HLE  ,4H    ,4HRID.,4HSOFT,4HSOFT,4H    ,4HDEFO,4HRMAT,4HIONS /


      II1=1
      IF(MEHAM.GT.1.OR.MEDEF.GT.0.OR.MEAXI.EQ.1.OR.MEVOL.GT.0) THEN     
      READ(20,2)HWIS(II1),AMB0IS(II1),AMG0IS(II1),
     *GAM0IS(II1),BET0IS(II1),BET4IS(II1),BB42IS(II1),GAMGIS(II1),
     *DELGIS(II1),BET3IS(II1),ETOIS(II1),AMUOIS(II1),HWOIS(II1),
     *BB32IS(II1),GAMDIS(II1),DPARIS(II1),GSHAEIS(II1)
      IF(MEPRI.LT.98) PRINT *, "Hamiltonian parameters are read"          
      END IF

 
C=======================================================================
C     LLMA-MAXIMUM MOMENTUM L
C     NCMA-MAXIMUM NUMBER OF COUPLED EQ.
C     NSMA-NUMBER OF SYSTEMS WITH  J AND PARITY
C     KODMA-SWITCH: 0-COUPLED STATES ARE ORDERED ONE BY ONE, NO MORE
C                   THAN NCMA
C                 1:-COUPLED STATES ARE ORDERED BY GROWING MOMENTUM L
C                    NO MORE THAN NCMA
C=======================================================================

           READ(20,211)NUR,NST,NPD,LAS,MTET,LLMA,NCMA,NSMA,KODMA
           IF(MEPRI.LT.98) PRINT *, "NUR etc... are read"
           
  211 FORMAT(20I3)
      IF(LLMA.EQ.0.OR.LLMA.GT.89) LLMA=89
      IF(NCMA.EQ.0.OR.NCMA.GT.200) NCMA=200
      IF(NSMA.EQ.0.OR.NSMA.GT.180) NSMA=180
            if(NST.LT.0) THEN
                 OPEN(99,FILE='OMINPUT.INP')
                 READ(99,*) NST
                 READ(99,*,END=212)(EE(I),I=1,NST)
  212            CLOSE(99)
              WRITE(*,*) NST
              DO I=1,NST
               WRITE(*,*) EE(I)
              MCHAE(I)=0
              ENDDO
            ELSE
             READ(20,2)(EE(I),I=1,NST)
             IF(MEPRI.LT.98) PRINT *, "EE(I) ... are read"
             READ(20,1)(MCHAE(I),I=1,NST)                          
             IF(MEPRI.LT.98) PRINT *, "MCHAE(I) ... are read"
            ENDIF
      IF(MTET.EQ.0) GO TO 13
           READ(20,2)(TET(I),I=1,MTET)
           IF(MEPRI.LT.98) PRINT *, "Angles are read"
   13 IF(MEPOT.GT.1) GO TO 16
       
           READ(20,3)(ELC(I),JOC(I),NPOC(I),KOC(I),NCAC(I),
     *     NUMBC(I),BETBC(I),AGSIC(I),NTUC(I),NNBC(I),NNGC(I),
     *                 NNOC(I),I=1,NUR)
           IF(MEPRI.LT.98) PRINT *, "ELC(I) etc... are read"
           
        GO TO 117
   16   READ(20,43)(ELC(I),JOC(I),NPOC(I), NTUC(I),NNBC(I),NNGC(I),
     *                 NNOC(I),NCAC(I),I=1,NUR)
        IF(MEPRI.LT.98) PRINT *, "ELC(I) etc... are read"
C====================================================================
C     VR=VR0+VR1*EN+VR2*EN*EN      AR=AR0+AR1*EN
C===================================================================
C                WD=WD0+WD1*EN     AD=AD0+AD1*EN
C     EN<BNDC    WC=WC0+WC1*EN     AC=AC0+AC1*EN
C ====================================================================
C                WD=WD0+WD1*BNDC+(EN-BNDC)*WDA1
C     EN>BNDC    WC=WC0+WC1*BNDC+(EN-BNDC)*WCA1
C                AD=AD0+AD1+BNDC
C====================================================================
        
  117 READ(20,211)NRESN
      IF(MEPRI.LT.98) PRINT *, "NRESN is read"
      
      IF(NRESN.EQ.0) GO TO 17      
           READ(20,213)(ERN(I),GNN(I),GREN(I),
     *         LON(I),JMN(I),JCON(I),NEL(10),I=1,NRESN)     
           IF(MEPRI.LT.98) PRINT *, "ERN(I) etc... are read"
  213      FORMAT(3E12.6,4I3)  

   17      READ(20,2)ANEU,ASP,AT,ZNUC,EFERMN,EFERMP
           IF(MEPRI.LT.98) PRINT *, "ANEU etc... are read"
           READ(20,2)VR0,VR1,VR2,VR3,VRLA,ALAVR,
     *               WD0,WD1,WDA1,WDBW,WDWID,ALAWD,
     *               WC0,WC1,WCA1,WCBW,WCWID,BNDC,
     *               VS,ALASO,WS0,WS1,WSBW,WSWID,
     *               RR,RRBWC,RRWID,PDIS,AR0,AR1,
     *               RD,AD0,AD1,RC,AC0,AC1,
     *               RW,AW0,AW1,RS,AS0,AS1,
     *               RZ,RZBWC,RZWID,AZ,CCOUL,ALF,
     *               CISO,WCISO,WDISO,EA,WDSHI,WDWID2,
     *               ALFNEW,VRD,CAVR,CARR,CAAR,CARD,
     *               CAAC,ATI,CAWD,CAWDW
           IF(MEPRI.LT.98) PRINT *, "Potential parameters are read"
     
      
C      !!!!!!!!COMPARING WITH DATET YOU NEED TO INPUT ATI - REFERENCE NUCLEI NUMBER IN ABCT !!!!      
      IF(MEPRI.LT.98) PRINT 500,ASP,AT
      WRITE(21,500)ASP,AT
  500 FORMAT( 7X,'INTERACTION OF PARTICLE, HAVING SPIN =',F5.2/19X,
     *'WITH NUCLEI',2X,'A=',F12.7/20X,'COUPLED CHANNELS METHOD')
      MESHH=MESHA-1
      IF(MEREL.EQ.0 .AND. MEPRI.LT.98) PRINT 134
      IF(MEREL.EQ.0) WRITE(21,134)
  134 FORMAT(22X,'NEWTON KINEMATICS')
      IF(MEREL.EQ.1 .AND. MEPRI.LT.98) PRINT 135
      IF(MEREL.EQ.1) WRITE(21,135)
  135 FORMAT(5X,'RELATIVISTIC KINEMATICS AND POTENTIAL DEPENDENCE')
      IF(MEREL.EQ.2 .AND. MEPRI.LT.98) PRINT 136
      IF(MEREL.EQ.2) WRITE(21,136)
  136 FORMAT(20X,'RELATIVISTIC KINEMATICS')
      IF(MEREL.EQ.3 .AND. MEPRI.LT.98) PRINT 137
      IF(MEREL.EQ.3) WRITE(21,137)
  137 FORMAT(3X,'RELATIVISTIC KINEMATICS AND REAL POTENTIAL DEPENDENCE')
C
      IF(MEDIS.EQ.0 .AND. MEPRI.LT.98) PRINT 184
      IF(MEDIS.EQ.0) WRITE(21,184)
  184 FORMAT(6X,'OPTICAL POTENTIAL WITHOUT DISPERSIVE RELATIONSHIPS')
      IF(MEDIS.GE.1 .AND. MEPRI.LT.98) PRINT 185
      IF(MEDIS.GE.1) WRITE(21,185)
  185 FORMAT(6X,'OPTICAL POTENTIAL WITH THE DISPERSIVE RELATIONSHIPS')
C
      IF(MECUL.EQ.0 .AND. MEPRI.LT.98) PRINT 154
      IF(MECUL.EQ.0) WRITE(21,154)
  154 FORMAT(5X,'COULOMB CORRECTION PROPORTIONAL REAL POTENTIAL DER-VE')
      IF(MECUL.EQ.1 .AND. MEPRI.LT.98) PRINT 155
      IF(MECUL.EQ.1) WRITE(21,155)
  155 FORMAT(15X,' COULOMB CORRECTION IS CONSTANT')
      IF(MECUL.EQ.2 .AND. MEPRI.LT.98) PRINT 156
      IF(MECUL.EQ.2) WRITE(21,156)
  156 FORMAT(/7X,' LANE CONSISTENT, EFFECTIVE PROTON ENERGY = E-CME,'/
     *13X,'BOTH FOR REAL AND IMAGINARY POTENTIALS'/)
      IF(MECUL.EQ.3 .AND. MEPRI.LT.98) PRINT 157
      IF(MECUL.EQ.3) WRITE(21,157)
  157 FORMAT(/7X,' LANE CONSISTENT, EFFECTIVE PROTON ENERGY = E-CME,'/
     *20X,' FOR REAL POTENTIAL ONLY'/)
C
      IF(MERZZ.EQ.0 .AND. MEPRI.LT.98) PRINT 164
      IF(MERZZ.EQ.0) WRITE(21,164)
  164 FORMAT(22X,'CHARGE RADIUS IS CONSTANT')
      IF(MERZZ.EQ.1 .AND. MEPRI.LT.98) PRINT 165
      IF(MERZZ.EQ.1) WRITE(21,165)
  165 FORMAT(15X,' CHARGE RADIUS IS ENERGY DEPENDENT')
C
      IF(MERRR.EQ.0 .AND. MEPRI.LT.98) PRINT 174
      IF(MERRR.EQ.0) WRITE(21,174)
  174 FORMAT(22X,'REAL RADIUS IS CONSTANT')
      IF(MERRR.EQ.1 .AND. MEPRI.LT.98) PRINT 175
      IF(MERRR.EQ.1) WRITE(21,175)
  175 FORMAT(15X,' REAL RADIUS IS ENERGY DEPENDENT')
C
      IF(MEDEF.EQ.1.AND.MEPOT.EQ.1) WRITE(21,195)
  195 FORMAT(/5X,'SOFT-ROTATOR NUCLEAR MODEL IS INVOLVED FOR
     * CALCULATIONS OF "EFFECTIVE"'/
     *15X,'DEFORMATIONS (COUPLING OF) WITH NON-GS BAND LEVELS'/)
      
      IF(MEDEF.EQ.2.AND.MEPOT.EQ.1) WRITE(21,201)
  201 FORMAT(/15X,'RIGID-ROTATOR NON-AXIAL NUCLEAR MODEL'/)
      
      IF(MEDEF.EQ.0) WRITE(21,196)
  196 FORMAT(/10X,'CALCULATIONS WITH INPUTTED "EFFECTIVE"'/
     *10X,'DEFORMATIONS (COUPLING OF)FOR NON-GS BAND LEVELS'/)
C     
      IF(MEAXI.EQ.0) WRITE(21,197)
  197 FORMAT(/20X,'CALCULATIONS FOR AXIAL CASE'/)
      IF(MEAXI.EQ.1) WRITE(21,198)
  198 FORMAT(/20X,'CALCULATIONS FOR NON-AXIAL CASE'/)
C     
      IF(MEVOL.EQ.0) WRITE(21,199)
  199 FORMAT(/20X,'CALCULATIONS WITHOUT VOLUME CONSERVATION'/)
      IF(MEVOL.EQ.1) WRITE(21,200)
  200 FORMAT(/16X,'CALCULATIONS WITH VOLUME CONSERVATION'/)  
      IF(MEVOL.EQ.2) WRITE(21,"(A80)")'CALCULATIONS WITH VOLUME 
     * CONSERVATION AND CoM IMMOBILITY'

      
      IF(MERAD.EQ.0.AND.MEVOL.GE.1) WRITE(21,202)
  202 FORMAT(/20X,'CALCULATIONS WITH STATIC PART OF VOLUME
     * CONSERVATION TERM'/)
      IF(MERAD.EQ.1) WRITE(21,203)
  203 FORMAT(/16X,'CALCULATIONS WITH RADIUS CORRECTION'/)        
      IF(MERAD.EQ.2) WRITE(21,204)
  204 FORMAT(/16X,'CALCULATIONS WITHOUT RADIUS CORRECTION'/)        
      
C
      IF(MESHA.GT.1 .AND. MEPRI.LT.98) PRINT 51,IFOR1,IFOR3(MESHH),IFOR2
      IF(MESHA.GT.1) WRITE(21,51)IFOR1,IFOR3(MESHH),IFOR2
   51 FORMAT(10X,14A4)
      IF(MEPRI.LT.98) PRINT 50,ITEXT,IMOD(MEHAM),IPOT,MPOT(MEPOT)
      WRITE (21,50)ITEXT,IMOD(MEHAM),IPOT,MPOT(MEPOT)
   50 FORMAT(/10X,4A4,6X,8A4)
      MEHA1=MEHAO+1
      IF(MESHO.GT.0 .AND. MEPRI.LT.98) 
     >  PRINT 45, IFOR1,ISMOD(MESHO),IFOR4,IHMOD(MEHA1),IFOR5
      IF(MESHO.GT.0) WRITE(21,45) IFOR1,ISMOD(MESHO),IFOR4,IHMOD(MEHAO),
     *IFOR5
   45 FORMAT(10X,13A4)
      IF(MEPRI.LT.98) PRINT 100,NUR,NPD,LAS
      WRITE(21,100)NUR,NPD,LAS
      IF(MEHAM.GT.1) GO TO 18
      IF(MEPRI.LT.98) PRINT 501,(I,ELC(I),JOC(I),NPOC(I),KOC(I),NCAC(I),
     *NUMBC(I),BETBC(I),AGSIC(I),NTUC(I),NNBC(I),NNGC(I),NNOC(I),
     *I=1,NUR)
      WRITE(21,501)(I,ELC(I),JOC(I),NPOC(I),KOC(I),NCAC(I),NUMBC(I),
     *BETBC(I),AGSIC(I),NTUC(I),NNBC(I),NNGC(I),NNOC(I),I=1,NUR)
      GO TO 19
   18 IF(MEPRI.LT.98) PRINT 20,(I,ELC(I),JOC(I),NTUC(I),NNBC(I),NNGC(I),
     *NNOC(I),NPOC(I),NCAC(I),I=1,NUR)
      WRITE(21,20)(I,ELC(I),JOC(I),NTUC(I),NNBC(I),NNGC(I),NNOC(I),
     *NPOC(I),NCAC(I),I=1,NUR)
      IF(MEPRI.LT.98) PRINT 21,HWIS(1),AMB0IS(1),AMG0IS(1),GAM0IS(1),
     * BET0IS(1),BET4IS(1),BB42IS(1),GAMGIS(1),DELGIS(1),BET3IS(1),
     * ETOIS(1),AMUOIS(1),HWOIS(1),BB32IS(1),GAMDIS(1),DPARIS(1),
     * GSHAEIS(1)
      WRITE(21,21)HWIS(1),AMB0IS(1),AMG0IS(1),GAM0IS(1),
     * BET0IS(1),BET4IS(1),BB42IS(1),GAMGIS(1),DELGIS(1),BET3IS(1),
     * ETOIS(1),AMUOIS(1),HWOIS(1),BB32IS(1),GAMDIS(1),DPARIS(1),
     * GSHAEIS(1)
   21 FORMAT(/22X,'PARAMETERS OF HAMILTONIAN '/5X,'HW=',F12.5,3X,
     *'AMB0=',F 8.5,3X,'AMG0=',F 8.5,3X,'GAM0=',F 8.5,3X,
     *'BET0=',F 8.5/
     *5X,'BET4=',F10.5,3X,'BB42=',F8.5,3X,'GAMG=',F8.5,3X,
     *'DELG=',F8.5/
     *5X,'BET3=',F10.5,3X,'ETO=',F9.5,3X,'AMUO=',F8.5,3X,
     *'HWO=',F8.5,4X,'BB32=',F8.5,3X/
     *5X,'GAMDE=',F9.5,3X,'DPAR=',F8.4,3X,'GSHAPE=',F8.5//)      
   20 FORMAT(//16X,'ENERGY ',4X,'LEVEL''S SPIN*2',4X,'NTU  ',
     *6X,'NNB  ', 6X,'NNG', 9X,'NNO',9X,'NPO',9X,'NCA'//
     *(1X,I4,8X,E14.7,7I11))
  501 FORMAT(//16X,'ENERGY',5X,'LEVEL''S SPIN*2',3X,'PARITY',10X,
     *'BAND*2',10X,'NCA',8X,'NUMB',9X,'BETB',11X,'ALFA(I)-> GS',
     *13X,'NTU',12X,'NNB',12X,'NNG',12X,'NNO'//
     *(1X,I4,6X,E12.5,I11,I14,I15,I15,I11,2E19.5,4I15))
 100  FORMAT( /15X,'NUMBER OF COUPLED LEVELS=',I3,5X,'NPD =',I2/14X,
     *'NUMBER OF TERMS IN POTENTIAL EXPANSION= ',2X,I2)
  19  IF(MEPRI.LT.98) PRINT 90
      WRITE(21,90)
  90  FORMAT(/15X,'POTENTIAL PARAMETERS V(R)')
      IF(MEPRI.LT.98) PRINT 91,VR0,VR1,VR2,RR,AR0,AR1,WD0,WD1,VR3,RD,
     *AD0,AD1,WC0,WC1,RC,AC0,AC1,RW,AW0,AW1,VS,RS,AS0,AS1,ALF,ANEU,RZ,
     *AZ,BNDC,WDA1,WCA1,CCOUL,CISO,WCISO,WS0,WS1,VRLA,ALAVR,
     *WCBW,WCWID,WDBW,WDWID,ALAWD,EFERMN,EFERMP,ALASO,PDIS,
     *WSBW,WSWID,RRBWC,RRWID,RZBWC,RZWID,EA,WDISO,WDSHI,WDWID2,
     *ALFNEW,VRD,CAVR,CARR,CAAR,CARD,CAAC
      WRITE(21,91)VR0,VR1,VR2,RR,AR0,AR1,WD0,WD1,VR3,RD,AD0,AD1,
     *WC0,WC1,RC,AC0,AC1,RW,AW0,AW1,VS,RS,AS0,AS1,ALF,ANEU,RZ,
     *AZ,BNDC,WDA1,WCA1,CCOUL,CISO,WCISO,WS0,WS1,VRLA,ALAVR,
     *WCBW,WCWID,WDBW,WDWID,ALAWD,EFERMN,EFERMP,ALASO,PDIS,
     *WSBW,WSWID,RRBWC,RRWID,RZBWC,RZWID,EA,WDISO,WDSHI,WDWID2,
     *ALFNEW,VRD,CAVR,CARR,CAAR,CARD,CAAC

   91 FORMAT(/1X,'VR0=',F7.3,5X,'VR1=',F7.4,5X,'VR2=',F10.7,2X,
     *'RR=',F7.4,5X,'AR0=',F7.4,5X,'AR1=',F7.4
     */1X,'WD0=',F7.4,5X,'WD1=',F7.4,5X,'VR3=',F10.7,2X,
     *'RD=',F7.4,5X,'AD0=',F7.4,5X,'AD1=',F7.4
     */1X,'WC0=',F7.4,5X,'WC1=',F7.4,21X,
     *'RC=',F7.4,5X,'AC0=',F7.4,5X,'AC1=',F7.4
     */49X,'RW=',F7.4,5X,'AW0=',F7.4,5X,'AW1=',F7.4
     */1X,'VSO=',F7.4,37X,'RS=',F7.4,5X,'AS0=',F7.4,5X,'AS1=',F7.4
     */1X,'ALF=',F7.4,5X,'ANEU=',F7.4,20X,'RZ=',F7.4,5X,'AZ0=',F7.4,
     */1X,'BNDC=',F7.2,4X,'WDA1=',F7.4,4X,'WCA1=',F7.4,4X,'CCOUL=',F7.4
     *,5X,'CISO=',F7.3,4X,'WCISO=',F7.3
     */1X,'WS0=',F7.4,5X,'WS1=',F7.4,5X,'VRLA=',F7.4
     *,4X,'ALAVR=',F8.5,4X,'WCBW=',F7.4,4X,'WCWID=',F7.4,/1X,'WDBW='
     *,F7.4,4X,'WDWID=',F7.4,3X,'ALAWD=',F7.4
     *,3X,'EFERMN=',F7.3,4X,'EFERMP=',F7.3,2X,'ALASO=',F7.4,
     */1X,'PDIS=',F7.4,4X,'WSBW=',F7.4,4X,'WSWID=',F7.2,3X,'RRBWC=',F7.4
     *,5X,'RRWID=',F6.2,4X,'RZBWC=',F7.4,
     */1X,'RZWID=',F7.4,3X,'EA=',F9.5,4X,'WDISO=',F7.3,
     *3X,'WDSHI=',F7.2,5X,'WDWID2=',F7.2,2X,'ALFNEW=',F6.3,
     */1X,'VRD=',F8.3,4X,'CAVR=',F8.5,3X,'CARR=',F9.6,
     *2X,'CAAR=',F9.6,4X,'CARD=',F9.6,2X,'CAAC=',F9.6/)

      IF(MEPRI.LT.98) PRINT 133,ZNUC,ATI
      WRITE(21,133) ZNUC,ATI
 133  FORMAT(/10X,'NUCLEUS CHARGE = ',F7.4,5x,
     *            'REFERENCE NUCLEUS   MASS = ',F7.3/)

      RCORR=1.d0
      
      IF(MEPOT.GT.1) GO TO 8
           IF(NPD.EQ.0) GO TO 8
           READ(20,2)(BET(I),I=2,NPD,2)
           IF(MEPRI.LT.98) PRINT *, "BET(I) are read"
      IF(MEPRI.LT.98) PRINT 96,(I,BET(I),I=2,NPD,2)
      WRITE(21,96)(I,BET(I),I=2,NPD,2)
  96  FORMAT(6X,'NPD',5X,'DEFORMATION PARAMETER VALUES'/
     *(6X,I2,13X,F7.4))
      
      BET2SUM=0.d0
      DO I=2,NPD,2
         !BET(I)=BETIS(IIIS,I)
         BET2SUM=BET2SUM+BET(I)**2
         END DO
      !BET2SUM=BET2SUM+BET3IS(1)**2

      IF(MERAD.EQ.1) RCORR=1.d0-BET2SUM*7.9577471546d-2 ! 1-bet2sum/(4*pi)      
      
      
    1 FORMAT(36I2)
    2 FORMAT(6E12.7)
    3 FORMAT(E12.7,5I2,2E12.7,4I2)
   43 FORMAT(E12.7,7I2)
    8 ASQ=AT**(1.D0/3.D0)
      VRLA=VRLA+CAVR*(AT-ATI)
      RR=(RR*RCORR+CARR*(AT-ATI))*ASQ
      RC=RC*RCORR*ASQ
      RD=(RD*RCORR+CARD*(AT-ATI))*ASQ
      !RD=RD*RCORR*ASQ
      RW=RW*RCORR*ASQ
      RS=RS*RCORR*ASQ
      RZ=RZ*RCORR*ASQ       
      AR0=AR0+CAAR*(AT-ATI)
      AC0=AC0+CAAC*(AT-ATI)
            
      WDBW=WDBW+CAWD*(AT-ATI)
      WDWID=WDWID+CAWDW*(AT-ATI)

!$OMP PARALLEL PRIVATE(IIparal,TID) 
!$OMP*  COPYIN(/MENU/)                                  ! PRIVCOM10
!$OMP*  COPYIN(/SHEM1/)                                 ! PRIVCOM12 
!$OMP*  COPYIN(/SHAMO/,/INRMi/,/ALFAi/,/INRM/,/ALFA/)   ! PRIVCOM9 
!$OMP*  COPYIN(/ENA/,/ENAi/,/ENAa/)                     ! PRIVCOM13 
!$OMP*  COPYIN(/DISK/,/RACB/,/DISCAN/,/POL1/)           ! PRIVCOM1 
!$OMP*  COPYIN(/ENB/,/INP0/,/INP1/)                     ! PRIVCOM15 
!$OMP*  COPYIN(/POT1/)                                  ! PRIVCOM6 
!$OMP*  COPYIN(/COUL/)                                  ! PRIVCOM14 
C
!$OMP*  COPYIN(/RESONI/,/POT2/,/POT3/,/POTD/,/DISPE2/)  ! PRIVCOM (a) 
!$OMP*  COPYIN(/RADi/,/DISPEi/,/CSBi/,/NCLMAi/)         ! PRIVCOM (b) 
!$OMP*  COPYIN(/NCLMA/)                                 ! PRIVCOM 
!$OMP*  COPYIN(/RAD/,/QNB/)                             ! PRIVCOM 

ccc!$OMP*  COPYIN(/LOFAC/)                                ! 
ccc!$OMP*  COPYIN(/INOUT/)                                !                        

C
!$OMP DO SCHEDULE(DYNAMIC,1)
C     This is a parallel loop
      DO IIparal=NST,1,-1
        CALL ABCTpar(IIparal)  
      ENDDO     
!$OMP ENDDO 
!$OMP END PARALLEL
      IF(MEPRI.NE.98) PRINT *, 'Exit from the parallel region'
      CALL THORA(21)

      RETURN
      END
C     *******************************************************
      SUBROUTINE DATET
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 

      CHARACTER*4 ITEXT,IMOD,IPOT,MPOT,IFOR1,IFOR2,
     *IFOR3,ISMOD,IFOR4,IHMOD,IFOR5
     
C     These common is used FOR initialization CCOULii <-> CCOUL
      INCLUDE 'PRIVCOM18D.FOR'
  
     
           
     
      INCLUDE 'PRIVCOM10.FOR'
      INCLUDE 'PRIVCOM12.FOR'
      
      INCLUDE 'PRIVCOM9.FOR'
      INCLUDE 'PRIVCOM13.FOR'
      INCLUDE 'PRIVCOM1.FOR'
      INCLUDE 'PRIVCOM.FOR'         
      INCLUDE 'PRIVCOM15.FOR'
      INCLUDE 'PRIVCOM6.FOR'
C
      INCLUDE 'PRIVCOM8.FOR'
C     COMMON/SHEMM/ES(40),JU(40),NTU(40),NNB(40),NNG(40),NNO(40),NPI(40)  

      INCLUDE 'PRIVCOM14.FOR'

      INCLUDE 'PRIVCOM17D.FOR'    ! not THREADPRIVATE
      INCLUDE 'PRIVCOM16D.FOR'    ! not THREADPRIVATE
     
      DIMENSION ITEXT(3),IMOD(7),IPOT(7),MPOT(2),IFOR1(2),IFOR2(9),
     *IFOR3(3),ISMOD(2),IFOR4(5),IHMOD(3),IFOR5(4), IBANDS(10)

      DATA ITEXT,IMOD/4HHAMI,4HLTON,4HIAN ,4H RV ,4H VM ,4HDCHM,
     *4H FDM,4H5PA0,4H 5PM,4HCLGB/,
     *IPOT,MPOT/4HPOTE,4HNTIA,4HL EX,4HPAND,4HED B,4HY    ,
     *4H    ,4HYL0 ,4HBET0/,IFOR1,IFOR2,IFOR3/4HWITH,4H AC.,
     *4HAXIA,4HL HE,4HXADE,4HCAPO,4HLE D,4HEFOR,4HMATI,4HONS ,
     *4H    ,4H    ,4H|NON,4H NON/,
     *ISMOD,IFOR4,IHMOD,IFOR5/4H    ,4H NON,4HAXIA,4HL OC,4HTUPO,
     *4HLE  ,4H    ,4HRID.,4HSOFT,4HSOFT,4H    ,4HDEFO,4HRMAT,4HIONS /
                READ(20,211)MENUC,MEBET,MEIIS,MERES,MELEV
C=======================================================================
C     
C     MENUC-NUMBER OF ADJUSTED ISOTOPES
C     MEBET-NUMBER OF ISOTOPE DEFORMATIONS OF WHICH ARE ADJUSTED, IF MEBET=0 THEN DEFORMATIONS 
C           OF ALL ISOTOPES ARE ADJUSTED (BETA2,BETA4,BETA6 FOR NOW)
C     MEIIS-NUMBER OF ISOTOPE FOR WHICH RESONSNCES ARE TO BE ADJUSTED
C     MERES-NUMBER OF RESONANCE FOR A CHOZEN ISOTOPE THAT IS TO BE ADJUSTED
C     MELEV-NUMBER OF LEVEL OF ISOTOPE WITH NUMBER MEBET DEFORMATION FOR WHICH IS TO BE ADJUSTED
C     LLMA-MAXIMUM MOMENTUM L
C     NCMA-MAXIMUM NUMBER OF COUPLED EQ.
C     NSMA-NUMBER OF SYSTEMS WITH  J AND PARITY
C     KODMA-SWITCH: 0-COUPLED STATES ARE ORDERED ONE BY ONE, NO MORE
C                   THAN NCMA
C                 1:-COUPLED STATES ARE ORDERED BY GROWING MOMENTUM L
C                    NO MORE THAN NCMA
C=======================================================================
           READ(20,211)NUR,NST,NPD,LAS,MTET,LLMA,NCMA,NSMA,KODMA
      PRINT *, "NUR etc... are read"
  211 FORMAT(20I3)
  234      FORMAT(6I3) 
      IF(LLMA.EQ.0.OR.LLMA.GT.89) LLMA=89
      IF(NCMA.EQ.0.OR.NCMA.GT.200) NCMA=200
      IF(NSMA.EQ.0.OR.NSMA.GT.180) NSMA=180
           READ(20,234)(NSTIS(I),NURIS(I),MESOIS(I),NRES(I),
     *     MEDEIS(I),MEAXIS(I), I=1, MENUC)
      PRINT *, "NSTIS etc... are read"
           READ(20,2) (WEI(I),I=1, MENUC)      
      PRINT *, "WEI etc... are read"
      DO 600 IIS=1,MENUC
      NUR=NURIS(IIS)
      NST=NSTIS(IIS)
      NRESN=NRES(IIS)
      MEDEF=MEDEIS(IIS)
      MEAXI=MEAXIS(IIS)
  
      IF(MEHAM.GT.1.OR.MEDEF.GT.0.OR.MEAXI.EQ.1.OR.MEVOL.GT.0) THEN     
          READ(20,2)HWIS(IIS),AMB0IS(IIS),AMG0IS(IIS),
     *GAM0IS(IIS),BET0IS(IIS),BET4IS(IIS),BB42IS(IIS),GAMGIS(IIS),
     *DELGIS(IIS),BET3IS(IIS),ETOIS(IIS),AMUOIS(IIS),HWOIS(IIS),
     *BB32IS(IIS),GAMDIS(IIS),DPARIS(IIS),GSHAEIS(IIS)
      PRINT *, "Ham params are read"          
      END IF     
      
           READ(20,2)(EEIS(IIS,I),I=1,NST)
      PRINT *, "EEIS etc... are read"
           READ(20,1)(MCHAIS(IIS,I),I=1,NST)
      PRINT *, "MCHAIS etc... are read"           
           IF(MEPOT.GT.1) GO TO 36
           READ(20,3)(ELIS(IIS,I),JOIS(IIS,I),NPOIS(IIS,I),
     *     KOIS(IIS,I),NCAIS(IIS,I),NUMBIS(IIS,I),BETBIS(IIS,I),
     *     AIGSIS(IIS,I),NTUIS(IIS,I),
     *     NNBIS(IIS,I),NNGIS(IIS,I),NNOIS(IIS,I),I=1,NUR)
           PRINT *, "ELIS etc... are read"  
           GO TO 37
   36      READ(20,43)(ELIS(IIS,I),JOIS(IIS,I),NPOIS(IIS,I),NTUIS(IIS,I)
     *     ,NNBIS(IIS,I),NNGIS(IIS,I),NNOIS(IIS,I),NCAIS(IIS,I),I=1,NUR)
           PRINT *, "ELIS etc... are read"  
C====================================================================
C     VR=VR0+VR1*EN+VR2*EN*EN      AR=AR0+AR1*EN
C===================================================================
C                WD=WD0+WD1*EN     AD=AD0+AD1*EN
C     EN<BNDC    WC=WC0+WC1*EN     AC=AC0+AC1*EN
C ====================================================================
C                WD=WD0+WD1*BNDC+(EN-BNDC)*WDA1
C     EN>BNDC    WC=WC0+WC1*BNDC+(EN-BNDC)*WCA1
C                AD=AD0+AD1+BNDC
C====================================================================
   37      IF(NRESN.EQ.0) GO TO 212      
           READ(20,213)(ERIS(IIS,I),GNIS(IIS,I),GREIS(IIS,I)
     *      ,LOIS(IIS,I),JMIS(IIS,I),JCOIS(IIS,I),NELA(IIS,I),I=1,NRESN) 
           PRINT *, "ERIS etc... are read"  
       
  213      FORMAT(3E12.6,4I3)  
  212      READ(20,2)ANEU,ASP,ATIS(IIS),ZNUCIS(IIS),EFISN(IIS),
     *     EFISP(IIS)
           PRINT *, "ANEU etc... are read"  
           
  600 CONTINUE 
        
         READ(20,2)VR0,VR1,VR2,VR3,VRLA,ALAVR,
     *               WD0,WD1,WDA1,WDBW,WDWID,ALAWD,
     *               WC0,WC1,WCA1,WCBW,WCWID,BNDC,
     *               VS,ALASO,WS0,WS1,WSBW,WSWID,
     *               RR,RRBWC,RRWID,PDIS,AR0,AR1,
     *               RD,AD0,AD1,RC,AC0,AC1,
     *               RW,AW0,AW1,RS,AS0,AS1,
     *               RZ,RZBWC,RZWID,AZ,CCOULii,ALF,
     *               CISO,WCISO,WDISO,EA,WDSHI,WDWID2,
     *               ALFNEW,VRD,CAVR,CARR,CAAR,CARD,
     *               CAAC,ABASE,CAWD,CAWDW
         PRINT *, "Potential params are read"
      IF(MEPRI.LT.98) 
     * PRINT 500,   ASP,(NINT(ATIS(I)),NINT(ZNUCIS(I)),I=1,MENUC)

      WRITE(21,500) ASP,(NINT(ATIS(I)),NINT(ZNUCIS(I)),I=1,MENUC)
  500 FORMAT( 7X,'INTERACTION OF PARTICLE HAVING SPIN =',F5.2/7X,
     *'WITH NUCLEI:',10(' A=',I3,' Z=',I3,' ;')/20X,
     *  'COUPLED CHANNELS METHOD')

      MESHH=MESHA-1
      IF(MEREL.EQ.0 .AND. MEPRI.LT.98) PRINT 134
      IF(MEREL.EQ.0) WRITE(21,134)
  134 FORMAT(22X,'NEWTON KINEMATICS')
      IF(MEREL.EQ.1 .AND. MEPRI.LT.98) PRINT 135
      IF(MEREL.EQ.1) WRITE(21,135)
  135 FORMAT(5X,'RELATIVISTIC KINEMATICS AND POTENTIAL DEPENDENCE')
      IF(MEREL.EQ.2 .AND. MEPRI.LT.98) PRINT 136
      IF(MEREL.EQ.2) WRITE(21,136)
  136 FORMAT(20X,'RELATIVISTIC KINEMATICS')
      IF(MEREL.EQ.3 .AND. MEPRI.LT.98) PRINT 137
      IF(MEREL.EQ.3) WRITE(21,137)
  137 FORMAT(3X,'RELATIVISTIC KINEMATICS AND REAL POTENTIAL DEPENDENCE')
C
      IF(MEDIS.EQ.0 .AND. MEPRI.LT.98) PRINT 184
      IF(MEDIS.EQ.0) WRITE(21,184)
  184 FORMAT(6X,'OPTICAL POTENTIAL WITHOUT DISPERSIVE RELATIONSHIPS')
      IF(MEDIS.GE.1 .AND. MEPRI.LT.98) PRINT 185
      IF(MEDIS.GE.1) WRITE(21,185)
  185 FORMAT(6X,'OPTICAL POTENTIAL WITH THE DISPERSIVE RELATIONSHIPS')
C
      IF(MECUL.EQ.0 .AND. MEPRI.LT.98) PRINT 154
      IF(MECUL.EQ.0) WRITE(21,154)
  154 FORMAT(5X,'COULOMB CORRECTION PROPORTIONAL REAL POTENTIAL DER-VE')
      IF(MECUL.EQ.1 .AND. MEPRI.LT.98) PRINT 155
      IF(MECUL.EQ.1) WRITE(21,155)
  155 FORMAT(15X,' COULOMB CORRECTION IS CONSTANT')
      IF(MECUL.EQ.2 .AND. MEPRI.LT.98) PRINT 156
      IF(MECUL.EQ.2) WRITE(21,156)
  156 FORMAT(/7X,' LANE CONSISTENT, EFFECTIVE PROTON ENERGY = E-CME,'/
     *13X,'BOTH FOR REAL AND IMAGINARY POTENTIALS'/)
      IF(MECUL.EQ.3 .AND. MEPRI.LT.98) PRINT 157
      IF(MECUL.EQ.3) WRITE(21,157)
  157 FORMAT(/7X,' LANE CONSISTENT, EFFECTIVE PROTON ENERGY = E-CME,'/
     *20X,' FOR REAL POTENTIAL ONLY'/)
C
      IF(MERZZ.EQ.0 .AND. MEPRI.LT.98) PRINT 164
      IF(MERZZ.EQ.0) WRITE(21,164)
  164 FORMAT(22X,'CHARGE RADIUS IS CONSTANT')
      IF(MERZZ.EQ.1 .AND. MEPRI.LT.98) PRINT 165
      IF(MERZZ.EQ.1) WRITE(21,165)
  165 FORMAT(15X,' CHARGE RADIUS IS ENERGY DEPENDENT')

C
      IF(MERRR.EQ.0 .AND. MEPRI.LT.98) PRINT 174
      IF(MERRR.EQ.0) WRITE(21,174)
  174 FORMAT(22X,'REAL RADIUS IS CONSTANT')
      IF(MERRR.EQ.1 .AND. MEPRI.LT.98) PRINT 175
      IF(MERRR.EQ.1) WRITE(21,175)
  175 FORMAT(15X,' REAL RADIUS IS ENERGY DEPENDENT')
C
      IF(MEDEF.EQ.1.AND.MEPOT.EQ.1) WRITE(21,195)
  195 FORMAT(/5X,'SOFT-ROTATOR NUCLEAR MODEL IS INVOLVED FOR
     * CALCULATIONS OF "EFFECTIVE"'/
     *15X,'DEFORMATIONS (COUPLING OF) WITH NON-GS BAND LEVELS'/)
      
      IF(MEDEF.EQ.2.AND.MEPOT.EQ.1) WRITE(21,201)
  201 FORMAT(/15X,'RIGID-ROTATOR NON-AXIAL NUCLEAR MODEL'/)
      
      IF(MEDEF.EQ.0) WRITE(21,196)
  196 FORMAT(/10X,'CALCULATIONS WITH INPUTTED "EFFECTIVE"'/
     *10X,'DEFORMATIONS (COUPLING OF)FOR NON-GS BAND LEVELS'/)
C     
      IF(MEAXI.EQ.0) WRITE(21,197)
  197 FORMAT(/20X,'CALCULATIONS FOR AXIAL CASE'/)
      IF(MEAXI.EQ.1) WRITE(21,198)
  198 FORMAT(/20X,'CALCULATIONS FOR NON-AXIAL CASE'/)
  
C     
      IF(MEVOL.EQ.0) WRITE(21,199)
  199 FORMAT(/20X,'CALCULATIONS WITHOUT VOLUME CONSERVATION'/)
      IF(MEVOL.EQ.1) WRITE(21,200)
  200 FORMAT(/16X,'CALCULATIONS WITH VOLUME CONSERVATION'/)  
      IF(MEVOL.EQ.2) WRITE(21,"(A80)")'CALCULATIONS WITH VOLUME 
     * CONSERVATION AND CoM IMMOBILITY'      
      
      IF(MERAD.EQ.0.AND.MEVOL.GE.1) WRITE(21,202)
  202 FORMAT(/20X,'CALCULATIONS WITH STATIC PART OF VOLUME
     * CONSERVATION TERM'/)
      IF(MERAD.EQ.1) WRITE(21,203)
  203 FORMAT(/16X,'CALCULATIONS WITH RADIUS CORRECTION'/)        
      IF(MERAD.EQ.2) WRITE(21,204)
  204 FORMAT(/16X,'CALCULATIONS WITHOUT RADIUS CORRECTION'/)        
      
C
      IF(MESHA.GT.1 .AND. MEPRI.LT.98) PRINT 51,IFOR1,IFOR3(MESHH),IFOR2
      IF(MESHA.GT.1) WRITE(21,51)IFOR1,IFOR3(MESHH),IFOR2
   51 FORMAT(10X,14A4)
      IF(MEPRI.LT.98) PRINT 50,ITEXT,IMOD(MEHAM),IPOT,MPOT(MEPOT)
      WRITE(21,50)ITEXT,IMOD(MEHAM),IPOT,MPOT(MEPOT)
   50 FORMAT(/10X,4A4,6X,8A4)
      MEHA1=MEHAO+1
      IF(MESHO.GT.0 .AND. MEPRI.LT.98) PRINT 45, IFOR1,ISMOD(MESHO),
     *  IFOR4,IHMOD(MEHA1),IFOR5
      IF(MESHO.GT.0) WRITE(21,45) IFOR1,ISMOD(MESHO),IFOR4,IHMOD(MEHAO),
     *IFOR5
   45 FORMAT(10X,13A4)


      DO 607 IIS=1,MENUC
      NUR=NURIS(IIS)
      NST=NSTIS(IIS)
      MESOL=MESOIS(IIS)
      MEDEF=MEDEIS(IIS)
      MEAXI=MEAXIS(IIS)
   
      IF(MEPOT.GT.1) GO TO 38
      DO 601 I=1, NUR
      EL(I)=ELIS(IIS,I)
      JO(I)=JOIS(IIS,I)
      NPO(I)=NPOIS(IIS,I)
      KO(I)=KOIS(IIS,I)
      NCA(I)=NCAIS(IIS,I)
      NUMB(I)=NUMBIS(IIS,I)
      AIGS(I)=AIGSIS(IIS,I)
      NTU(I)=NTUIS(IIS,I) 
      NNB(I)=NNBIS(IIS,I)
      NNG(I)=NNGIS(IIS,I)
      NNO(I)=NNOIS(IIS,I)
      ES(I)=EL(I) 
      JU(I)=JO(I)/2
      NPI(I)=NPO(I)
  601 BETB(I)=BETBIS(IIS,I)

      IF(MEPOT.EQ.1.AND.MEDEF.EQ.0.AND.MEAXI.EQ.0) GO TO 876

      HW=HWIS(IIS)
      AMB0=AMB0IS(IIS)
      AMG0=AMG0IS(IIS)
      GAM0=GAM0IS(IIS)
      BET0=BET0IS(IIS)
      BET4=BET4IS(IIS)
      BB42=BB42IS(IIS)
      GAMG=GAMGIS(IIS)
      DELG=DELGIS(IIS)
      BET3=BET3IS(IIS)
      ETO=ETOIS(IIS)
      AMUO=AMUOIS(IIS)
      HWO=HWOIS(IIS)
      BB32=BB32IS(IIS)
      GAMDE=GAMDIS(IIS)
      DPAR=DPARIS(IIS)
      GSHAPE=GSHAEIS(IIS)
      
 
  
      
  876 continue
       
      
                   
       GO TO 456
     
  456 IF(MEPRI.LT.98) PRINT 100,NUR,NPD,LAS
      WRITE(21,100)NUR,NPD,LAS
      
      IF(MEPRI.LT.98) PRINT 501,(I,EL(I),JO(I),NPO(I),KO(I),NCA(I)
     *,NUMB(I),BETB(I),
     *AIGS(I),NTU(I),NNB(I),NNG(I),NNO(I),I=1,NUR)
      WRITE(21,501)(I,EL(I),JO(I),NPO(I),KO(I),NCA(I),NUMB(I),BETB(I),
     *AIGS(I),NTU(I),NNB(I),NNG(I),NNO(I),I=1,NUR)
      IF(MEDEF.GT.0) GO TO 388
      GO TO 39
   38 DO 602 I=1, NUR
      EL(I)=ELIS(IIS,I)
      JO(I)=JOIS(IIS,I)
      NTU(I)=NTUIS(IIS,I)
      NNB(I)=NNBIS(IIS,I)
      NNG(I)=NNGIS(IIS,I)
      NNO(I)=NNOIS(IIS,I)      
      NPO(I)=NPOIS(IIS,I)
      NCA(I)=NCAIS(IIS,I)
      NUMB(I)=NUMBIS(IIS,I)
  602 BETB(I)=BETBIS(IIS,I)

  388 HW=HWIS(IIS)
      AMB0=AMB0IS(IIS)
      AMG0=AMG0IS(IIS)
      GAM0=GAM0IS(IIS)
      BET0=BET0IS(IIS)
      BET4=BET4IS(IIS)
      BB42=BB42IS(IIS)
      GAMG=GAMGIS(IIS)
      DELG=DELGIS(IIS)
      BET3=BET3IS(IIS)
      ETO=ETOIS(IIS)
      AMUO=AMUOIS(IIS)
      HWO=HWOIS(IIS)
      BB32=BB32IS(IIS)
      GAMDE=GAMDIS(IIS)
      DPAR=DPARIS(IIS)
      GSHAPE=GSHAEIS(IIS)    
      IF(MEPOT.EQ.1) GO TO 389
      IF(MEPRI.LT.98)PRINT 40,NINT(ATIS(IIS)),NINT(ZNUCIS(IIS)), 
     *   (I,EL(I),JO(I),NTU(I),NNB(I),NNG(I),NNO(I),NPO(I),NCA(I)
     *,I=1,NUR)
      WRITE(21,40) NINT(ATIS(IIS)),NINT(ZNUCIS(IIS)),
     *   (I,EL(I),JO(I),NTU(I),NNB(I),NNG(I),NNO(I),NPO(I),
     *NCA(I),I=1,NUR)
   40 FORMAT(//26X,'TARGET NUCLEUS A=',I3,' Z=',I3//
     *16X,'ENERGY ',4X,'LEVEL''S SPIN*2',4X,'NTU  ',
     *6X,'NNB  ', 6X,'NNG', 9X,'NNO',9X,'NPO',9X,'NCA'//
     *(1X,I4,8X,E14.7,7I11))

  389 IF(MEPRI.LT.98)PRINT 41,HW,AMB0,AMG0,GAM0,BET0,BET4,BB42,GAMG,DELG
     *,BET3,ETO,AMUO,HWO,BB32,GAMDE,DPAR,GSHAPE
      WRITE(21,41)HW,AMB0,AMG0,GAM0,BET0,BET4,BB42,GAMG,DELG
     *,BET3,ETO,AMUO,HWO,BB32,GAMDE,DPAR,GSHAPE
   41 FORMAT(/22X,'PARAMETERS OF HAMILTONIAN '/5X,'HW=',F12.5,3X,
     *'AMB0=',F 8.5,3X,'AMG0=',F 8.5,3X,'GAM0=',F 8.5,3X,
     *'BET0=',F 8.5/
     *5X,'BET4=',F10.5,3X,'BB42=',F8.5,3X,'GAMG=',F8.5,3X,
     *'DELG=',F8.5/
     *5X,'BET3=',F10.5,3X,'ETO=',F9.5,3X,'AMUO=',F8.5,3X,
     *'HWO=',F8.5,4X,'BB32=',F8.5,3X/
     *5X,'GAMDE=',F9.5,3X,'DPAR=',F8.4,3X,'GSHAPE=',F8.5//)

      IF(MEPOT.EQ.1) GO TO 39
      IF(MEHAM.GT.2) CALL PREQU
  501 FORMAT(//16X,'ENERGY',5X,'LEVEL''S SPIN*2',3X,'PARITY',10X,
     *'BAND*2',10X,'NCA',8X,'NUMB',9X,'BETB',11X,'ALFA(I)-> GS',
     *13X,'NTU',12X,'NNB',12X,'NNG',12X,'NNO'//
     *(1X,I4,6X,E12.5,I11,I14,I15,I15,I11,2E19.5,4I15))
 100  FORMAT( /15X,'NUMBER OF COUPLED LEVELS=',I3,5X,'NPD =',I2/14X,
     *'NUMBER OF TERMS IN POTENTIAL EXPANSION= ',2X,I2)
  39  EFERMNi = EFISN(IIS)
      EFERMPi = EFISP(IIS)  
      IF(MEPRI.LT.98) PRINT 90
      WRITE (21,90)
  90  FORMAT(/15X,'POTENTIAL PARAMETERS V(R)')
      IF(MEPRI.LT.98)
     *PRINT 91,VR0,VR1,VR2,RR,AR0,AR1,WD0,WD1,VR3,RD,AD0,AD1,
     *WC0,WC1,RC,AC0,AC1,RW,AW0,AW1,VS,RS,AS0,AS1,ALF,ANEU,RZ,
     *AZ,BNDC,WDA1,WCA1,CCOULii,CISO,WCISO,WS0,WS1,VRLA,ALAVR,
     *WCBW,WCWID,WDBW,WDWID,ALAWD,EFERMN,EFERMP,ALASO,PDIS,
     *WSBW,WDWID,RRBWC,RRWID,RZBWC,RZWID,EA,WDISO,WDSHI,WDWID2,
     *ALFNEW,VRD,CAVR,CARR,CAAR,CARD,CAAC
      IF(MEPRI.LT.98)
     *WRITE(21,91)VR0,VR1,VR2,RR,AR0,AR1,WD0,WD1,VR3,RD,AD0,AD1,
     *WC0,WC1,RC,AC0,AC1,RW,AW0,AW1,VS,RS,AS0,AS1,ALF,ANEU,RZ,
     *AZ,BNDC,WDA1,WCA1,CCOULii,CISO,WCISO,WS0,WS1,VRLA,ALAVR,
     *WCBW,WCWID,WDBW,WDWID,ALAWD,EFERMN,EFERMP,ALASO,PDIS,
     *WSBW,WSWID,RRBWC,RRWID,RZBWC,RZWID,EA,WDISO,WDSHI,WDWID2,
     *ALFNEW,VRD,CAVR,CARR,CAAR,CARD,CAAC
   91 FORMAT(/1X,'VR0=',F7.3,5X,'VR1=',F7.4,5X,'VR2=',F10.7,2X,
     *'RR=',F7.4,5X,'AR0=',F7.4,5X,'AR1=',F7.4
     */1X,'WD0=',F7.4,5X,'WD1=',F7.4,5X,'VR3=',F10.7,2X,
     *'RD=',F7.4,5X,'AD0=',F7.4,5X,'AD1=',F7.4
     */1X,'WC0=',F7.4,5X,'WC1=',F7.4,21X,
     *'RC=',F7.4,5X,'AC0=',F7.4,5X,'AC1=',F7.4
     */49X,'RW=',F7.4,5X,'AW0=',F7.4,5X,'AW1=',F7.4
     */1X,'VSO=',F7.4,37X,'RS=',F7.4,5X,'AS0=',F7.4,5X,'AS1=',F7.4
     */1X,'ALF=',F7.4,5X,'ANEU=',F7.4,20X,'RZ=',F7.4,5X,'AZ0=',F7.4,
     */1X,'BNDC=',F7.2,4X,'WDA1=',F7.4,4X,'WCA1=',F7.4,4X,'CCOUL=',F7.4
     *,5X,'CISO=',F7.3,4X,'WCISO=',F7.3
     */1X,'WS0=',F7.4,5X,'WS1=',F7.4,5X,'VRLA=',F7.4
     *,4X,'ALAVR=',F8.5,4X,'WCBW=',F7.4,4X,'WCWID=',F7.4,/1X,'WDBW='
     *,F7.4,4X,'WDWID=',F7.4,3X,'ALAWD=',F7.4
     *,3X,'EFERMN=',F7.3,4X,'EFERMP=',F7.3,2X,'ALASO=',F7.4,
     */1X,'PDIS=',F7.4,4X,'WSBW=',F7.4,4X,'WSWID=',F7.2,3X,'RRBWC=',F7.4
     *,5X,'RRWID=',F6.2,4X,'RZBWC=',F7.4,
     */1X,'RZWID=',F7.4,3X,'EA=',F9.5,4X,'WDISO=',F7.3,
     *3X,'WDSHI=',F7.2,5X,'WDWID2=',F7.2,2X,'ALFNEW=',F6.3,
     */1X,'VRD=',F8.3,4X,'CAVR=',F8.5,3X,'CARR=',F9.6,
     *2X,'CAAR=',F9.6,4X,'CARD=',F9.6,2X,'CAAC=',F9.6)
C     IF(MEPRI.LT.98) PRINT 133,(NINT(ZNUCIS(I)), I=1,MENUC)
C     WRITE(21,133) (NINT(ZNUCIS(I)), I=1,MENUC)
C133  FORMAT(/130X,'NUCLEI(I) CHARGE = ',10(I3,1x)/)
      IF(MEPOT.GT.1) GO TO 607
      IF(NPD.EQ.0) GO TO 607
      READ(20,2)(BETIS(IIS,I),I=2,NPD,2)
      PRINT *, "BETIS etc... are read" 
      IF(MEPRI.LT.98) PRINT 96,(I,BETIS(IIS,I),I=2,NPD,2)
      WRITE(21,96)(I,BETIS(IIS,I),I=2,NPD,2)
   96  FORMAT(6X,'NPD',5X,'DEFORMATION PARAMETER VALUES'/
     *(6X,I2,13X,F7.4))
  607 CONTINUE
C
    1 FORMAT(36I2)
    2 FORMAT(6E12.7)
    3 FORMAT(E12.7,5I2,2E12.7,4I2)
   43 FORMAT(E12.7,7I2)
    6 FORMAT(2E12.7,2I2)
    8 CONTINUE
C

C     Storing READ values INTO new VARIABLES (RR -> RRi)
      INCLUDE 'PRIVCOM19.FOR' 

      RRG=RRi
      RCG=RCi
      RDG=RDi
      RWG=RWi
      RSG=RSi
      RZG=RZi
      VRG=VRLAi
      ARG=AR0i
      ACG=AC0i
      
      WDBWG=WDBWi
      WDWIDG=WDWIDi

      READ(20,112)(NPJ(I),I=1,77)
      PRINT *, "NPJ(I) etc... are read"
      IF(MEPRI.LT.98) PRINT 99
      WRITE (21,99)

      WRITE (21,99)
   99 FORMAT(/10X,'PARAMETERS ADJUSTED'/)
      IF(MEPRI.LT.98) PRINT 111,(NPJ(I),I=1,77)
      WRITE(21,111)(NPJ(I),I=1,77)
  111 FORMAT(1X,6I2)
  112 FORMAT(6I2)
      DO 4 IIS=1,MENUC
      NST=NSTIS(IIS)
      NUR=NURIS(IIS)
      MESOL=MESOIS(IIS)
      MEDEF=MEDEIS(IIS)
      MEAXI=MEAXIS(IIS)
      DO 4 I=1,NST
      IF(MEPRI.LT.98) PRINT 114, EEIS(IIS,I)
      WRITE (21,114) EEIS(IIS,I)
  114 FORMAT(//6X,'EXPERIMENTAL DATA FOR ENERGY=',F10.6,1X,'MeV'/)
           READ(20,1) NT(IIS,I),NR(IIS,I),NGN(IIS,I),NGD(IIS,I),
     *     NSF1(IIS,I),NSF2(IIS,I),NRAT(IIS,I),NNAT(IIS,I)
      IF(MEPRI.LT.98) PRINT 111,NT(IIS,I),NR(IIS,I),NGN(IIS,I),
     *     NGD(IIS,I),NSF1(IIS,I),NSF2(IIS,I),NRAT(IIS,I),NNAT(IIS,I)
      WRITE(21,111)NT(IIS,I),NR(IIS,I),NGN(IIS,I),NGD(IIS,I),
     *NSF1(IIS,I),NSF2(IIS,I),NRAT(IIS,I),NNAT(IIS,I)
           READ(20,2)STE(IIS,I),DST(IIS,I),SRE(IIS,I),DSR(IIS,I),
     *     RATIO(IIS,I),DRAT(IIS,I)
           IF(NNAT(IIS,I).NE.0) READ(20,2) CSNAT(IIS,I),DCSNAT(IIS,I)
      IF(MEPRI.LT.98) PRINT 222,STE(IIS,I),DST(IIS,I),SRE(IIS,I),
     *  DSR(IIS,I),RATIO(IIS,I),DRAT(IIS,I),CSNAT(IIS,I),DCSNAT(IIS,I) 
      WRITE(21,222)STE(IIS,I),DST(IIS,I),SRE(IIS,I),DSR(IIS,I),
     *RATIO(IIS,I),DRAT(IIS,I), CSNAT(IIS,I),DCSNAT(IIS,I)
  222 FORMAT(1X,6E12.7)
      IF(NSF1(IIS,I).EQ.0.AND.NSF2(IIS,I).EQ.0) GO TO 9
           READ(20,2)SE1(IIS,I),DS1(IIS,I),SE2(IIS,I),DS2(IIS,I)
      IF(MEPRI.LT.98) PRINT 222,SE1(IIS,I),DS1(IIS,I),SE2(IIS,I),
     *                          DS2(IIS,I)
      WRITE(21,222)SE1(IIS,I),DS1(IIS,I),SE2(IIS,I),DS2(IIS,I)
    9 NG=NGN(IIS,I)
      IF(NG.EQ.0) GO TO 5
  987 FORMAT (2E12.6,2I3)
           READ(20,6)(SNE(IIS,I,K),DSN(IIS,I,K),NIN(IIS,I,K),
     *     NFN(IIS,I,K),K=1,NG)
      IF(MEPRI.LT.98) PRINT 666,(SNE(IIS,I,K),DSN(IIS,I,K),
     *     NIN(IIS,I,K),NFN(IIS,I,K),K=1,NG)
      WRITE(21,666)(SNE(IIS,I,K),DSN(IIS,I,K),NIN(IIS,I,K),
     *     NFN(IIS,I,K),K=1,NG)
  666 FORMAT(1X,2E12.7,2I3)
    5 NG=NGD(IIS,I)
      IF(NG.EQ.0) GO TO 4
           READ(20,1)(NID(IIS,I,K),NFD(IIS,I,K),MTD(IIS,I,K),K=1,NG)
      IF(MEPRI.LT.98) PRINT 111,(NID(IIS,I,K),NFD(IIS,I,K),
     *                           MTD(IIS,I,K),K=1,NG) 
      WRITE(21,111)(NID(IIS,I,K),NFD(IIS,I,K),MTD(IIS,I,K),K=1,NG)
      DO 7 K=1,NG
      M=MTD(IIS,I,K)
           READ(20,2)(TED(IIS,I,K,L),SNGD(IIS,I,K,L),
     *     DSD(IIS,I,K,L),L=1,M)
      IF(MEPRI.LT.98) PRINT 223,(TED(IIS,I,K,L),SNGD(IIS,I,K,L),
     *                           DSD(IIS,I,K,L),L=1,M)
      WRITE(21,223)(TED(IIS,I,K,L),SNGD(IIS,I,K,L),DSD(IIS,I,K,L),L=1,M)
  223 FORMAT(1X,6E12.7)    
    7 CONTINUE
    4 CONTINUE
      KEV=0
      IF(NPJ(1).NE.1) GO TO 11
      KEV=KEV+1
      XAD(KEV)=VR0
   11 IF(NPJ(2).NE.1) GO TO 12
      KEV=KEV+1
      XAD(KEV)=VR1
   12 IF(NPJ(3).NE.1) GO TO 13
      KEV=KEV+1
      XAD(KEV)=VR2
   13 IF(NPJ(4).NE.1) GO TO 14
      KEV=KEV+1
      XAD(KEV)=VR3
   14 IF(NPJ(5).NE.1) GO TO 15
      KEV=KEV+1
      XAD(KEV)=VRG
   15 IF(NPJ(6).NE.1) GO TO 16
      KEV=KEV+1
      XAD(KEV)=ALAVR
   16 IF(NPJ(7).NE.1) GO TO 17
      KEV=KEV+1
      XAD(KEV)=WD0
   17 IF(NPJ(8).NE.1) GO TO 18
      KEV=KEV+1
      XAD(KEV)=WD1
   18 IF(NPJ(9).NE.1) GO TO 19
      KEV=KEV+1
      XAD(KEV)=WDA1
   19 IF(NPJ(10).NE.1) GO TO 20
      KEV=KEV+1
      XAD(KEV)=WDBWG
   20 IF(NPJ(11).NE.1) GO TO 21
      KEV=KEV+1
      XAD(KEV)=WDWIDG
   21 IF(NPJ(12).NE.1) GO TO 22
      KEV=KEV+1
      XAD(KEV)=ALAWD
   22 IF(NPJ(13).NE.1) GO TO 23
      KEV=KEV+1
      XAD(KEV)=WC0
   23 IF(NPJ(14).NE.1) GO TO 24
      KEV=KEV+1
      XAD(KEV)=WC1
   24 IF(NPJ(15).NE.1) GO TO 25
      KEV=KEV+1
      XAD(KEV)=WCA1
   25 IF(NPJ(16).NE.1) GO TO 26
      KEV=KEV+1
      XAD(KEV)=WCBW
   26 IF(NPJ(17).NE.1) GO TO 27
      KEV=KEV+1
      XAD(KEV)=WCWID
   27 IF(NPJ(18).NE.1) GO TO 28
      KEV=KEV+1
      XAD(KEV)=BNDC
   28 IF(NPJ(19).NE.1) GO TO 29
      KEV=KEV+1
      XAD(KEV)=VS
   29 IF(NPJ(20).NE.1) GO TO 30
      KEV=KEV+1
      XAD(KEV)=ALASO
   30 IF(NPJ(21).NE.1) GO TO 31
      KEV=KEV+1
      XAD(KEV)=WS0
   31 IF(NPJ(22).NE.1) GO TO 32
      KEV=KEV+1
      XAD(KEV)=WS1
   32 IF(NPJ(23).NE.1) GO TO 33
      KEV=KEV+1
      XAD(KEV)=WSBW
   33 IF(NPJ(24).NE.1) GO TO 53
      KEV=KEV+1
      XAD(KEV)=WSWID
   53 IF(NPJ(25).NE.1) GO TO 54
      KEV=KEV+1
      XAD(KEV)=RRG
   54 IF(NPJ(26).NE.1) GO TO 55
      KEV=KEV+1
      XAD(KEV)=RRBWC
   55 IF(NPJ(27).NE.1) GO TO 56
      KEV=KEV+1
      XAD(KEV)=RRWID
   56 IF(NPJ(28).NE.1) GO TO 57
      KEV=KEV+1
      XAD(KEV)=PDIS
   57 IF(NPJ(29).NE.1) GO TO 58
      KEV=KEV+1
      XAD(KEV)=ARG
   58 IF(NPJ(30).NE.1) GO TO 59
      KEV=KEV+1
      XAD(KEV)=AR1
   59 IF(NPJ(31).NE.1) GO TO 60
      KEV=KEV+1
      XAD(KEV)=RDG
   60 IF(NPJ(32).NE.1) GO TO 61
      KEV=KEV+1
      XAD(KEV)=AD0
   61 IF(NPJ(33).NE.1) GO TO 62
      KEV=KEV+1
      XAD(KEV)=AD1
   62 IF(NPJ(34).NE.1) GO TO 63
      KEV=KEV+1
      XAD(KEV)=RCG
   63 IF(NPJ(35).NE.1) GO TO 64
      KEV=KEV+1
      XAD(KEV)=ACG
   64 IF(NPJ(36).NE.1) GO TO 65
      KEV=KEV+1
      XAD(KEV)=AC1
   65 IF(NPJ(37).NE.1) GO TO 66
      KEV=KEV+1
      XAD(KEV)=RWG
   66 IF(NPJ(38).NE.1) GO TO 67
      KEV=KEV+1
      XAD(KEV)=AW0
   67 IF(NPJ(39).NE.1) GO TO 68
      KEV=KEV+1
      XAD(KEV)=AW1
   68 IF(NPJ(40).NE.1) GO TO 69
      KEV=KEV+1
      XAD(KEV)=RSG
   69 IF(NPJ(41).NE.1) GO TO 70
      KEV=KEV+1
      XAD(KEV)=AS0
   70 IF(NPJ(42).NE.1) GO TO 71
      KEV=KEV+1
      XAD(KEV)=AS1
   71 IF(NPJ(43).NE.1) GO TO 72
      KEV=KEV+1
      XAD(KEV)=RZG
   72 IF(NPJ(44).NE.1) GO TO 73
      KEV=KEV+1
      XAD(KEV)=RZBWC
   73 IF(NPJ(45).NE.1) GO TO 74
      KEV=KEV+1
      XAD(KEV)=RZWID
   74 IF(NPJ(46).NE.1) GO TO 75
      KEV=KEV+1
      XAD(KEV)=AZ
   75 IF(NPJ(47).NE.1) GO TO 76
      KEV=KEV+1
      XAD(KEV)=CCOULii
   76 IF(NPJ(48).NE.1) GO TO 77
      KEV=KEV+1
      XAD(KEV)=ALF
   77 IF(NPJ(49).NE.1) GO TO 78
      KEV=KEV+1
      XAD(KEV)=CISO
   78 IF(NPJ(50).NE.1) GO TO 79
      KEV=KEV+1
      XAD(KEV)=WCISO
   79 IF(NPJ(51).NE.1) GO TO 80
      KEV=KEV+1
      XAD(KEV)=WDISO
   80 IF(NPJ(52).NE.1) GO TO 81
      KEV=KEV+1
      XAD(KEV)=EA
   81 IF(NPJ(53).NE.1) GO TO 82
      KEV=KEV+1
      XAD(KEV)=WDSHI
   82 IF(NPJ(54).NE.1) GO TO 83
      KEV=KEV+1
      XAD(KEV)=WDWID2
   83 IF(NPJ(55).NE.1) GO TO 84
      KEV=KEV+1
      XAD(KEV)=ALFNEW
   84 IF(NPJ(56).NE.1) GO TO 85
      KEV=KEV+1
      XAD(KEV)=VRD
   85 IF(NPJ(57).NE.1) GO TO 86
      KEV=KEV+1
      XAD(KEV)=BET0IS(MEBET)
   86 IF(NPJ(58).NE.1) GO TO 87
      KEV=KEV+1
      XAD(KEV)=BET3IS(MEBET)
   87 IF(NPJ(59).NE.1) GO TO 88
      KEV=KEV+1
      XAD(KEV)=BET4IS(MEBET)
   88 IF(NPJ(60).NE.1) GO TO 899
        if(MEBET.eq.0) then
          XAD(KEV+1:KEV+MENUC)=BETIS(1:MENUC,2)
          KEV=KEV+MENUC
        else
          KEV=KEV+1
          XAD(KEV)=BETIS(MEBET,2)
        endif
  899 IF(NPJ(61).NE.1) GO TO 901
        if(MEBET.eq.0) then
          XAD(KEV+1:KEV+MENUC)=BETIS(1:MENUC,4)
          KEV=KEV+MENUC
        else
          KEV=KEV+1
          XAD(KEV)=BETIS(MEBET,4)
        endif
  901 IF(NPJ(62).NE.1) GO TO 910
        if(MEBET.eq.0) then
          XAD(KEV+1:KEV+MENUC)=BETIS(1:MENUC,6)
          KEV=KEV+MENUC
        else
          KEV=KEV+1
          XAD(KEV)=BETIS(MEBET,6)
        endif
  910 IF(NPJ(63).NE.1) GO TO 92
      KEV=KEV+1
      XAD(KEV)=AMUOIS(MEBET)
   92 IF(NPJ(64).NE.1) GO TO 750
      KEV=KEV+1
      XAD(KEV)=AMG0IS(MEBET)
  750 IF(NPJ(65).NE.1) GO TO 751
      KEV=KEV+1
      XAD(KEV)=CAVR  
  751 IF(NPJ(66).NE.1) GO TO 752
      KEV=KEV+1
      XAD(KEV)=CARR      
  752 IF(NPJ(67).NE.1) GO TO 753
      KEV=KEV+1
      XAD(KEV)=CAAR
  753 IF(NPJ(68).NE.1) GO TO 754
      KEV=KEV+1
      XAD(KEV)=CARD 
  754 IF(NPJ(69).NE.1) GO TO 755
      KEV=KEV+1
      XAD(KEV)=CAAC
  755 IF(NPJ(70).NE.1) GO TO 756
      KEV=KEV+1
      XAD(KEV)=ERIS(MEIIS,MERES)
  756 IF(NPJ(71).NE.1) GO TO 757
      KEV=KEV+1
      XAD(KEV)= ABS(GNIS(MEIIS,MERES)) 
  757 IF(NPJ(72).NE.1) GO TO 758
      KEV=KEV+1
      XAD(KEV)=ABS(GREIS(MEIIS,MERES))
  758 IF(NPJ(73).NE.1) GO TO 759
        if(MELEV.eq.0.and.MEBET.ne.0)then
          IBANDS=-1
          IIIB=0
          do III=1,NURIS(MEBET)
            if(.not.ANY(IBANDS.eq.NUMBIS(MEBET,III))
     *             .and.NUMBIS(MEBET,III).ne.NUMBIS(MEBET,1)) then
              IIIB=IIIB+1
              IBANDS(IIIB)=NUMBIS(MEBET,III)
              KEV=KEV+1
              XAD(KEV)=ABS(BETBIS(MEBET,III)) 
            endif
          enddo 
        elseif(MELEV.ne.0.and.MEBET.ne.0)then
          KEV=KEV+1
          XAD(KEV)=ABS(BETBIS(MEBET,MELEV)) 
        else
          print *,'Not supported yet.'  
          stop
        endif
  759 IF(NPJ(74).NE.1) GO TO 760
      KEV=KEV+1
      XAD(KEV)=ABS(GAM0IS(MEBET))    
  760 IF(NPJ(75).NE.1) GO TO 761
      KEV=KEV+1
      XAD(KEV)=ABS(AMB0IS(MEBET))
  761 IF(NPJ(76).NE.1) GO TO 762
      KEV=KEV+1
      XAD(KEV)=CAWD
  762 IF(NPJ(77).NE.1) GO TO 93
      KEV=KEV+1
      XAD(KEV)=CAWDW
   93 NV=KEV       

      READ(20,2)(EP(K),K=1,NV)
      IF(MEPRI.LT.98) PRINT 222,(EP(K),K=1,NV)
      WRITE(21,222)(EP(K),K=1,NV)
           READ(20,2)FU

      CALL SEART
      RETURN
      END
C     *******************************************************
      SUBROUTINE SEART
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM10.FOR'  
      INCLUDE 'PRIVCOM17D.FOR'  
C     COMMON/OPT/XAD(25),GR(25),XAD1(25),XAD2(25),EP(25),EPSGR(25),NV
C     COMMON/OPB/C,GRR(25),FM,EPS1,NRL
      !DOUBLE PRECISION SCUR(25),SPREV(25),GRPREV(25),OMEGA
      
      !SCUR=0.0
      !SPREV=0.0
      !GRPREV=0.0
      !GR=0.0
      
      NCC=1
      NI=0
      DO 1 I=1,NV
      XAD1(I)=XAD(I)
    1 XAD2(I)=XAD(I)+EP(I)*2.D0
      IF(FU.EQ.0.) CALL XISQT
      FM=FU
   10 IF(MEPRI.LT.98) PRINT 101,NI,FM
      WRITE(21,101)NI,FM
      IF(MEPRI.LT.98) PRINT 102,(I,XAD(I),I=1,NV)
      WRITE(21,102)(I,XAD(I),I=1,NV)
  101 FORMAT(1X,'NI=',I5,14X,'FM=',D30.15)
  102 FORMAT(/1X,4(1X,' NV',4X, 6X,' X ', 6X)/
     *       (1X,4(I4,4X,G15.9)) )
    6 NXX=0
      N1C=0
      FU=FM
      DO 21 I=1,NV
      XAD(I)=XAD1(I)
      IF(DABS(XAD2(I)-XAD1(I)).LE.DABS(EP(I))) NXX=NXX+1
   21 CONTINUE
      IF(NXX.EQ.NV)GO TO 9
      
      !GRPREV=GR
     
      CALL DEFGT
      
      !SPREV=SCUR
      !IF(SUM(GRPREV*GRPREV).GT.0) THEN
      !   OMEGA=MAX(0.0,SUM((GR-GRPREV)*GR)/SUM(GRPREV*GRPREV))
      !ELSE
      !    OMEGA=0.0
      !ENDIF   
      !SCUR=GR+OMEGA*SPREV

      
      LL=0
      NNK=0
      DO 13 I=1,NV
   13 XAD2(I)=XAD1(I)
      NI=NI+1
      C=0.
      DO 7 I=1,NV
      GRR(I)=GR(I)!SCUR(I)!
    7 C=C+GRR(I)**2*(300*EP(I))**2
      !C=SUM(SCUR*GR*(300*EP)**2)
      EPS1=FM/C
    8 NX=0
      DO 30 I=1,NV
      EPSGR(I)=EPS1*GRR(I)*(300*EP(I))**2
      IF(DABS(EPSGR(I)).GT.0.3D0*DABS(XAD1(I))) GO TO 15
      IF(DABS(EPSGR(I)).LE.DABS(EP(I))) NX=NX+1
   30 CONTINUE
      IF(NNK.EQ.0) NX=0
   17 DO 2 I=1,NV
      HS=EPS1*GRR(I)*(300*EP(I))**2
      HS1=XAD1(I)
      XAD(I)=HS1-HS
    2 CONTINUE
      CALL XISQT
      
CNaN       WRITE(30,36) FU
C      PRINT 36, FU
C      PRINT *, 'Outside the parallel region: FU=',SNGL(FU)
CNaN   36  FORMAT(E12.5)
          
CNaN       BACKSPACE 30
CNaN       READ (30,37)RNAN
CNaN   37  FORMAT(A4)
C      PRINT 37, RNAN,NAN 
CNaN       IF(RNAN.EQ.NAN) GO TO 4
       IF(FU.NE.FU) GO TO 4
      
      IF(FU-FM) 3,4,4
    3 DO 5 I=1,NV
    5 XAD1(I)=XAD(I)
      LL=LL+1
      IF(LL.GT.3)EPS1=EPS1*5.D0
      N1C=1
      NCC=0
      FM=FU
      IF(MEPRI.LT.98) PRINT 101,NI,FM
      WRITE(21,101)NI,FM
      IF(MEPRI.LT.98) PRINT 102,(I,XAD(I),I=1,NV)
      WRITE(21,102)(I,XAD(I),I=1,NV)
      GO TO 17
    4 IF(NX.EQ.NV) GO TO 6
      IF(N1C.EQ.1)NCC=NCC+1
      IF(NCC-2)15,14,14
   14 EPS1=-EPS1
      NCC=0
      GO TO 17
   15 EPS1=EPS1/5.D0
      NNK=1
      LL=0
      NCC=1
      GO TO 8
    9 CONTINUE
    
      RETURN
      END
C     *******************************************************
      SUBROUTINE DEFGT
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      INCLUDE 'PRIVCOM17D.FOR'
C      COMMON/OPT/XAD(25),GR(25),XAD1(25),XAD2(25),EP(25),EPSGR(25),NV
C      COMMON/OPB/C,GRR(25),FM,EPS1,NRL
C      COMMON/CHISQC/FU
      DIMENSION TEMPAD(35)
      F1=FU
      DO I=1,NV
      DL=EP(I)
      TEMPAD(I)=XAD(I)
      XAD(I)=XAD(I)+DL
C     write(*,*) 'XAD(i)=',i,XAD(i) 
      CALL XISQT
C     write(*,*) 'XAD(i)=',i,XAD(i) 
      GR(I)=(FU-F1)/DL
      !XAD(I)=XAD(I)-DL
      XAD(I)=TEMPAD(I)
      ENDDO
      RETURN
      END
C     *******************************************************
C     *******************************************************
      SUBROUTINE XISQT
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      DIMENSION XPRN(35), IBANDS(10),DBANDS(10)
      
C     These common is used FOR initialization CCOULii <-> CCOUL
      INCLUDE 'PRIVCOM18D.FOR'

      INCLUDE 'PRIVCOM10.FOR'
      INCLUDE 'PRIVCOM12.FOR'
      INCLUDE 'PRIVCOM9.FOR'
      INCLUDE 'PRIVCOM13.FOR'
      INCLUDE 'PRIVCOM1.FOR'
      INCLUDE 'PRIVCOM.FOR'         
      INCLUDE 'PRIVCOM15.FOR'
      INCLUDE 'PRIVCOM6.FOR'
C
      INCLUDE 'PRIVCOM8.FOR'
C     COMMON/SHEMM/ES(40),JU(40),NTU(40),NNB(40),NNG(40),NNO(40),NPI(40)  

      INCLUDE 'PRIVCOM14.FOR'

      INCLUDE 'PRIVCOM17D.FOR'    ! not THREADPRIVATE
      INCLUDE 'PRIVCOM16D.FOR'    ! not THREADPRIVATE

      INCLUDE 'PRIVCOM22.FOR'   

      INTEGER IIparal
        
C     APRN=1.D0/ATIS(1)**(1.D0/3.D0)
C     APRN=1.D0 

C     Restoring READ values for other nuclei in the loop after PARALLEL execution 
      INCLUDE 'PRIVCOM18.FOR' 
            
      KEV=0
      IF(NPJ(1).NE.1) GO TO 11
      KEV=KEV+1
      VR0=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   11 IF(NPJ(2).NE.1) GO TO 12
      KEV=KEV+1
      VR1=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   12 IF(NPJ(3).NE.1) GO TO 13
      KEV=KEV+1
      VR2=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   13 IF(NPJ(4).NE.1) GO TO 14
      KEV=KEV+1
      VR3=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   14 IF(NPJ(5).NE.1) GO TO 15
      KEV=KEV+1
      VRG=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
   15 IF(NPJ(6).NE.1) GO TO 16
      KEV=KEV+1
      ALAVR=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   16 IF(NPJ(7).NE.1) GO TO 17
      KEV=KEV+1
      WD0=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
   17 IF(NPJ(8).NE.1) GO TO 18
      KEV=KEV+1
      WD1=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
   18 IF(NPJ(9).NE.1) GO TO 19
      KEV=KEV+1
      WDA1=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
   19 IF(NPJ(10).NE.1) GO TO 20
      KEV=KEV+1
      WDBWG=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
   20 IF(NPJ(11).NE.1) GO TO 21
      KEV=KEV+1
      WDWIDG=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
   21 IF(NPJ(12).NE.1) GO TO 22
      KEV=KEV+1
      ALAWD=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
   22 IF(NPJ(13).NE.1) GO TO 23
      KEV=KEV+1
      WC0=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
   23 IF(NPJ(14).NE.1) GO TO 24
      KEV=KEV+1
      WC1=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
   24 IF(NPJ(15).NE.1) GO TO 25
      KEV=KEV+1
      WCA1=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
   25 IF(NPJ(16).NE.1) GO TO 26
      KEV=KEV+1
      WCBW=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   26 IF(NPJ(17).NE.1) GO TO 27
      KEV=KEV+1
      WCWID=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   27 IF(NPJ(18).NE.1) GO TO 28
      KEV=KEV+1
      BNDC=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   28 IF(NPJ(19).NE.1) GO TO 34
      KEV=KEV+1
      VS=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   34 IF(NPJ(20).NE.1) GO TO 35
      KEV=KEV+1
      ALASO=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   35 IF(NPJ(21).NE.1) GO TO 36
      KEV=KEV+1
      WS0=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   36 IF(NPJ(22).NE.1) GO TO 37
      KEV=KEV+1
      WS1=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   37 IF(NPJ(23).NE.1) GO TO 38
      KEV=KEV+1
      WSBW=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   38 IF(NPJ(24).NE.1) GO TO 52
      KEV=KEV+1
      WSWID=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
   52 IF(NPJ(25).NE.1) GO TO 53
      KEV=KEV+1
      RRG=ABS(XAD(KEV))
C     XPRN(KEV)=XAD(KEV)*APRN
      XPRN(KEV)=XAD(KEV)
   53 IF(NPJ(26).NE.1) GO TO 54
      KEV=KEV+1
      RRBWC=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   54 IF(NPJ(27).NE.1) GO TO 55
      KEV=KEV+1
      RRWID=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   55 IF(NPJ(28).NE.1) GO TO 56
      KEV=KEV+1
      PDIS=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   56 IF(NPJ(29).NE.1) GO TO 57
      KEV=KEV+1
      ARG=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
   57 IF(NPJ(30).NE.1) GO TO 58
      KEV=KEV+1
      AR1=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   58 IF(NPJ(31).NE.1) GO TO 59
      KEV=KEV+1
      RDG=XAD(KEV)
C     XPRN(KEV)=XAD(KEV)*APRN
      XPRN(KEV)=XAD(KEV)
   59 IF(NPJ(32).NE.1) GO TO 60
      KEV=KEV+1
      AD0=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
   60 IF(NPJ(33).NE.1) GO TO 61
      KEV=KEV+1
      AD1=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
   61 IF(NPJ(34).NE.1) GO TO 62
      KEV=KEV+1
      RCG=ABS(XAD(KEV))
C     XPRN(KEV)=XAD(KEV)*APRN
      XPRN(KEV)=XAD(KEV)
C      rc=rr
   62 IF(NPJ(35).NE.1) GO TO 63
      KEV=KEV+1
      ACG=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
C      ac0=ar0
   63 IF(NPJ(36).NE.1) GO TO 64
      KEV=KEV+1
      AC1=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   64 IF(NPJ(37).NE.1) GO TO 65
      KEV=KEV+1
      RWG=XAD(KEV)
C     XPRN(KEV)=XAD(KEV)*APRN
      XPRN(KEV)=XAD(KEV)
   65 IF(NPJ(38).NE.1) GO TO 66
      KEV=KEV+1
      AW0=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
   66 IF(NPJ(39).NE.1) GO TO 67
      KEV=KEV+1
      AW1=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   67 IF(NPJ(40).NE.1) GO TO 68
      KEV=KEV+1
      RSG=XAD(KEV)
C     XPRN(KEV)=XAD(KEV)*APRN
      XPRN(KEV)=XAD(KEV)
   68 IF(NPJ(41).NE.1) GO TO 69
      KEV=KEV+1
      AS0=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   69 IF(NPJ(42).NE.1) GO TO 70
      KEV=KEV+1
      AS1=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   70 IF(NPJ(43).NE.1) GO TO 71
      KEV=KEV+1
      RZG=XAD(KEV)
C     XPRN(KEV)=XAD(KEV)*APRN
      XPRN(KEV)=XAD(KEV)
   71 IF(NPJ(44).NE.1) GO TO 72
      KEV=KEV+1
      RZBWC=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   72 IF(NPJ(45).NE.1) GO TO 73
      KEV=KEV+1
      RZWID=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   73 IF(NPJ(46).NE.1) GO TO 74
      KEV=KEV+1
      AZ=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
   74 IF(NPJ(47).NE.1) GO TO 75
      KEV=KEV+1
      CCOULii=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   75 IF(NPJ(48).NE.1) GO TO 76
      KEV=KEV+1
      ALF=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   76 IF(NPJ(49).NE.1) GO TO 77
      KEV=KEV+1
      CISO=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   77 IF(NPJ(50).NE.1) GO TO 78
      KEV=KEV+1
      WCISO=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
   78 IF(NPJ(51).NE.1) GO TO 79
      KEV=KEV+1
      WDISO=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
   79 IF(NPJ(52).NE.1) GO TO 80
      KEV=KEV+1
      EA=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
   80 IF(NPJ(53).NE.1) GO TO 81
      KEV=KEV+1
      WDSHI=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   81 IF(NPJ(54).NE.1) GO TO 82
      KEV=KEV+1
      WDWID2=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   82 IF(NPJ(55).NE.1) GO TO 83
      KEV=KEV+1
      ALFNEW=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
   83 IF(NPJ(56).NE.1) GO TO 84
      KEV=KEV+1
      VRD=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
   84 IF(NPJ(57).NE.1) GO TO 85
      KEV=KEV+1
      BET0IS(MEBET)=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   85 IF(NPJ(58).NE.1) GO TO 86
      KEV=KEV+1
      CBM=BET3IS(MEBET)/AMUOIS(MEBET)
      BET3IS(MEBET)=XAD(KEV)
      AMUOIS(MEBET)=BET3IS(MEBET)/CBM
      XPRN(KEV)=XAD(KEV)
   86 IF(NPJ(59).NE.1) GO TO 87
      KEV=KEV+1
      BET4IS(MEBET)=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
   87 IF(NPJ(60).NE.1) GO TO 88
        if(MEBET.eq.0) then
          BETIS(1:MENUC,2)=XAD(KEV+1:KEV+MENUC)
          XPRN(KEV+1:KEV+MENUC)=XAD(KEV+1:KEV+MENUC)
          KEV=KEV+MENUC
        else
          KEV=KEV+1
          BETIS(MEBET,2)=XAD(KEV)
          XPRN(KEV)=XAD(KEV)
        endif
   88 IF(NPJ(61).NE.1) GO TO 899
        if(MEBET.eq.0) then
          BETIS(1:MENUC,4)=XAD(KEV+1:KEV+MENUC)
          XPRN(KEV+1:KEV+MENUC)=XAD(KEV+1:KEV+MENUC)
          KEV=KEV+MENUC
        else
          KEV=KEV+1
          BETIS(MEBET,4)=XAD(KEV)
          XPRN(KEV)=XAD(KEV)
        endif
  899 IF(NPJ(62).NE.1) GO TO 901
        if(MEBET.eq.0) then
          BETIS(1:MENUC,6)=XAD(KEV+1:KEV+MENUC)
          XPRN(KEV+1:KEV+MENUC)=XAD(KEV+1:KEV+MENUC)
          KEV=KEV+MENUC
        else
          KEV=KEV+1
          BETIS(MEBET,6)=XAD(KEV)
          XPRN(KEV)=XAD(KEV)
        endif
  901 IF(NPJ(63).NE.1) GO TO 910
      KEV=KEV+1
      CMB=AMUOIS(MEBET)**2*BB32
      CBM=BET3IS(MEBET)/AMUOIS(MEBET)
      AMUOIS(MEBET)=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
      IF(BET3.EQ.0.) BB32=CMB/AMUOIS(MEBET)**2
      IF(BET3.NE.0.) BET3IS(MEBET)=AMUOIS(MEBET)*CBM
  910 IF(NPJ(64).NE.1) GO TO 751
      KEV=KEV+1
      AMG0IS(MEBET)=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
  751 IF(NPJ(65).NE.1) GO TO 752
      KEV=KEV+1
      CAVR=XAD(KEV)
      XPRN(KEV)=XAD(KEV)  
  752 IF(NPJ(66).NE.1) GO TO 753
      KEV=KEV+1
      CARR=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
  753 IF(NPJ(67).NE.1) GO TO 754
      KEV=KEV+1
      CAAR=XAD(KEV)
      XPRN(KEV)=XAD(KEV)
  754 IF(NPJ(68).NE.1) GO TO 755
      KEV=KEV+1
      CARD=XAD(KEV)
      XPRN(KEV)=XAD(KEV)  
  755 IF(NPJ(69).NE.1) GO TO 756
      KEV=KEV+1
      CAAC=XAD(KEV)
      XPRN(KEV)=XAD(KEV)        
  756 IF(NPJ(70).NE.1) GO TO 757
      KEV=KEV+1
      ERIS(MEIIS,MERES)=XAD(KEV)
      IF(NELA(MEIIS,MERES).NE.0) ERIS(MEIIS,NELA(MEIIS,MERES))=XAD(KEV)
      XPRN(KEV)=XAD(KEV)  
  757 IF(NPJ(71).NE.1) GO TO 758
      KEV=KEV+1
      GNIS(MEIIS,MERES)=ABS(XAD(KEV))
      IF(NELA(MEIIS,MERES).NE.0) GNIS(MEIIS,NELA(MEIIS,MERES))=
     *ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)        
  758 IF(NPJ(72).NE.1) GO TO 759
      KEV=KEV+1
      GREIS(MEIIS,MERES)=ABS(XAD(KEV))
      IF(NELA(MEIIS,MERES).NE.0) GREIS(MEIIS,NELA(MEIIS,MERES))=
     *ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)
  759 IF(NPJ(73).NE.1) GO TO 760
        if(MELEV.eq.0.and.MEBET.ne.0)then
          IBANDS=-1
          DBANDS=0.0
          IIIB=1
          IBANDS(1)=NUMBIS(MEBET,1)
          do III=1,NURIS(MEBET)
            if(NUMBIS(MEBET,III).ne.NUMBIS(MEBET,1))then
              do IIII=1,IIIB
                if(NUMBIS(MEBET,III).eq.IBANDS(IIII))then
                  BETBIS(MEBET,III)=DBANDS(IIII)
                  exit
                endif
              enddo
              if(IIII.eq.IIIB+1)then
                IIIB=IIIB+1
                IBANDS(IIIB)=NUMBIS(MEBET,III)
                KEV=KEV+1
                DBANDS(IIIB)=ABS(XAD(KEV))
                BETBIS(MEBET,III)=DBANDS(IIIB)
                XPRN(KEV)=XAD(KEV)              
              endif
            endif 
          enddo
        elseif(MELEV.ne.0.and.MEBET.ne.0)then
          KEV=KEV+1
          !BETBIS(MEBET,MELEV)=ABS(XAD(KEV))
          DO III=1, NURIS(MEBET)
             IF(NUMBIS(MEBET,III).eq.NUMBIS(MEBET,MELEV)) THEN
                 BETBIS(MEBET,III)=ABS(XAD(KEV))
             ENDIF
          ENDDO
          XPRN(KEV)=XAD(KEV)
        else
          print *,'Not supported yet.'  
          stop
        endif      

      
  760 IF(NPJ(74).NE.1) GO TO 761
      KEV=KEV+1
      GAM0IS(MEBET)=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)      
  761 IF(NPJ(75).NE.1) GO TO 762
      KEV=KEV+1
      AMB0IS(MEBET)=ABS(XAD(KEV))
      XPRN(KEV)=XAD(KEV)        
      
  762 IF(NPJ(76).NE.1) GO TO 763
      KEV=KEV+1
      CAWD=XAD(KEV)
      XPRN(KEV)=XAD(KEV)        
  763 IF(NPJ(77).NE.1) GO TO 93
      KEV=KEV+1
      CAWDW=XAD(KEV)
      XPRN(KEV)=XAD(KEV)        
      
   93 FU=0.D0
c      IF(NPJ(58).EQ.1.OR.NPJ(63).EQ.1.OR.NPJ(64).EQ.1.OR.NPJ(74).EQ.1.
c     *OR.NPJ(75).EQ.1) CALL PREQU  
      NNTT=0
      NNTM=0
      FUM=0.D0

      KODMAi = KODMA
      NPDi   = NPD
      LASi   = LAS
      ASPin  = ASP

      DO 600 IIS=1, MENUC

      maxit = NSTIS(IIS)
C
C     Restoring READ values for other nuclei in the loop after PARALLEL execution 
      KODMA = KODMAi
      NPD   = NPDi
      LAS   = LASi
      ATI   = 1
      ASP   = ASPin
      CCOUL=CCOULii
C     PRINT*,CCOULii,CCOUL,IIPARAL,IIS      

!$OMP PARALLEL DEFAULT (PRIVATE) SHARED (maxit) 
!$OMP*  COPYIN(/MENU/)                                ! PRIVCOM10
!$OMP*  COPYIN(/ENAa/)                                ! PRIVCOM13  
!$OMP*  COPYIN(/INP1/)                                ! PRIVCOM15 
!$OMP*  COPYIN(/POT1/)                                ! PRIVCOM6 
C
!$OMP*  COPYIN(/RESONI/,/POT2/,/POT3/,/POTD/,/DISPE2/)! PRIVCOM (a) 
!$OMP*  COPYIN(/RADi/,/DISPEi/,/NCLMAi/)              ! PRIVCOM (b) 
!$OMP*  COPYIN(/NCLMA/)                               ! PRIVCOM 
!$OMP*  COPYIN(/RAD/,/QNB/)                           ! PRIVCOM 
!$OMP*  COPYIN(/NIND/)                                ! psssing IIS 
C
!$OMP DO SCHEDULE(DYNAMIC,1)
C     This is a parallel loop
      DO IIparal= maxit,1,-1
      CALL DATETpar(IIparal)
      ENDDO     
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

      NNTI=NNTT-NNTM
      NNTI=NNTT-NNTM
      FI=FU-FUM
      FUM=FU
      NNTM=NNTT
      FI=FI/NNTI

      IF(MEPRI.LT.98) PRINT 160,NINT(ATIS(IIS)),NINT(ZNUCIS(IIS)),
     * FI,FU/NNTT!ATIS(IIS), FI

      WRITE (*,122) NINT(ATIS(IIS)), FI, NNTI, FU/NNTT, NNTT
  122 FORMAT(1X,/,1X,
     * 'Outside parallel reg. for A=',I3,' FI=',D13.6,
     * ' NNTI=',I4,' FU=',D13.6,' NNTT=',I4/)
      
      WRITE (21,160) NINT(ATIS(IIS)),NINT(ZNUCIS(IIS)), FI,FU/NNTT
  160 FORMAT(/1X,'NUCLEUS MASS IS=',I3, 3X,'CHARGE IS=',I3,
     *        5X,'FU FOR NUCLEUS IS=',D14.8, ' CUMULATIVE FU=',D14.8/)
         
  600 CONTINUE 

      FU=FU/DBLE(NNTT)

      IF(KEV.LE.4) THEN
        IF(MEPRI.LT.98) PRINT 1380,(I,XAD(I),I=1,KEV)
        WRITE(21,1380)(I,XAD(I),I=1,KEV)
        PRINT 1380,(I,XPRN(I),I=1,KEV)
 1380   FORMAT(1X,4(1X,'KEV',4X,6X,' X ',6X)/
     *        (1X,4(I4,4X,G15.9)) )
        WRITE(21,1380)(I,XPRN(I),I=1,KEV)
      ELSE 
        IF(MEPRI.LT.98) PRINT 138,(I,XAD(I),I=1,KEV)
        WRITE(21,138)(I,XAD(I),I=1,KEV)
        WRITE(21,138)(I,XPRN(I),I=1,KEV)
      ENDIF 
  138 FORMAT(1X,6(1X,'KEV',4X,6X,' X ',6X)/
     *           (1X,6(I4,4X,G15.9)) )
C
      PRINT 139,FU
      WRITE (21,139) FU
  139 FORMAT(/1X,'FU=',E14.7,' starting next iteration'/)
      IF(NPJ(63).EQ.1.OR.NPJ(58).EQ.1 .AND. MEPRI.LT.98)
     *            PRINT 140,AMUOIS(MEBET)
      IF(NPJ(63).EQ.1.OR.NPJ(58).EQ.1) WRITE (21,140) AMUOIS(MEBET)
  140 FORMAT (40X,'AMUO=',F20.12/)
C
      CALL THORA(21)
      RETURN
      END

      CHARACTER*1 FUNCTION cpar(integ)
      INTEGER integ
      cpar = '*'
      if(integ.eq.-1) cpar = '-'
      if(integ.eq.+1) cpar = '+'
      RETURN
      END
     
     
C     *******************************************************
C     END of optmand
C     *******************************************************
