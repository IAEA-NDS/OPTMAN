c-----------------------------------------------------------------------------
c subroutine SPRT
c G. NOGUERE - CEA Cadarache (11/2015)
c R-Matrix parameters, hard sphere approximation, direct contribution
c
c Calculate the following parameters:
c .STRt      : "total" strenght function form SPRT and ESW
c .STRc      : compound contribution
c .STRd      : direct contribution
c .PHIomc    : phase shift from omc
c .PHIsqw    : phase shift from ESW
c .RPRIME    : effective radius 
c .RINFINITY : distant level parameter
c
c-----------------------------------------------------------------------------

	subroutine SPRT(lo,TLL,SRL,SIL,AT,EN,STRc,RPRIME)
	implicit none

        double precision  TLL,SRL,SIL,AT,EN
        double precision  STRt,STRc,STRd,PHIomc,PHIsqw,RPRIME,RINFINITY

        double precision  ll,A,Elab
        double precision  Tlcoef(0:100)
        double precision  ReS(0:100),ImS(0:100),ReC(0:100),ImC(0:100)
        integer           lo
        integer           nbLmax
        integer           i
        double precision  ac,k,Rfrohner
        double precision  P(0:100),Shift(0:100),Phi(0:100)
        double precision  Rinf(0:100)
        double precision  Str(0:100),Scomp(0:100),Sdir(0:100)
        double precision  pole(0:100),pcomp(0:100),pdir(0:100)
        double precision  Tcoef(0:100),Xspot(0:100),Xstot(0:100)
        double precision  deltaL(0:100),RprimeL(0:100),PhiNew(0:100)
        double precision  mod
	
c-----------------------------------------------------------------------------
c Initialize parameter
c-----------------------------------------------------------------------------

	Tlcoef(lo)=TLL
        ReS(lo)=SRL
        ImS(lo)=SIL
        A=AT
        Elab=EN*1.0e06

        ll=FLOAT(lo)
        nbLmax=lo-1
	
c-----------------------------------------------------------------------------
c Transmission min
c-----------------------------------------------------------------------------

c       if(TLL.le.0.0000001)return

c-----------------------------------------------------------------------------
c Energy max
c-----------------------------------------------------------------------------

       if(EN.gt.1.0)return

c-----------------------------------------------------------------------------
c Collision matrix elements 	
c-----------------------------------------------------------------------------
	
        ImC(lo)=(1.0-ReS(lo))/2.0
        ReC(lo)=ImS(lo)/2.0

c-----------------------------------------------------------------------------
c Channel radius 	
c-----------------------------------------------------------------------------

	call CR(ac,A)
	
c	write(6,'("  - Channel Radius (fm)",26x," = ",f9.4)')ac
c	write(6,*)
	
c-----------------------------------------------------------------------------
c Wave number 	
c-----------------------------------------------------------------------------
	 
        call WaveNumber(k,A,Elab)
 
c-----------------------------------------------------------------------------
c Penetrability factor and Phase shift
c Potential scattering phase shift 	
c-----------------------------------------------------------------------------

	call PenF(P,Shift,ll,k,ac)
	call PhiF(P,Shift,Phi,ll)

c-----------------------------------------------------------------------------
c R-Matrix parameters 	
c-----------------------------------------------------------------------------
	
c	write(6,'("  R-matrix elements from generalized SPRT ...")')

	call Rinfiny(Rinf(lo),P(lo),Phi(lo),ReC(lo),ImC(lo))
	call polefunction(pole(lo),P(lo),Phi(lo),ReC(lo),ImC(lo))
	call pole2Sfunction(Str(lo),pole(lo),P(0),Elab)

c	write(6,'("  - Strenght function ",27x," = ",f8.4)')Str(lo)
c	write(6,'("  - Distant parameter ",27x," = ",f8.4)')Rinf(lo)
c	write(6,*)' '

c-----------------------------------------------------------------------------
c Cross sections 	
c-----------------------------------------------------------------------------
	
c	write(6,'("  Cross sections ...")')
	
        call shape(XSpot(lo),P(lo),Phi(lo),Rinf(lo),pole(lo),k)
        call TransCoef(Tcoef(lo),P(lo),Rinf(lo),pole(lo))
        call XStotal(XStot(lo),P(lo),Phi(lo),Rinf(lo),pole(lo),k)

c	write(6,'("  - Shape cross section ",25x," = ",E15.7)')Xspot(lo)/1000.0
c	write(6,'("  - Transmission coef   ",25x," = ",E15.7)')Tcoef(lo)
c	write(6,'("  - Total cross section ",25x," = ",E15.7)')Xstot(lo)/1000.0
c	write(6,*)' '
	

c-----------------------------------------------------------------------------
c Phase shift from OMC	
c
c Initial value : Radius from Frohner formula or channel radius ac
c-----------------------------------------------------------------------------

c	write(6,'("  Hard sphere phase shift from OMC ...")')

	Rfrohner=ac*(1.0-(2.0*ll+1.0)*Rinf(lo))**(1.0/(2.0*ll+1.0))
	call Phase(ReS(lo),ImS(lo),deltaL(lo))
	call FindRprime(deltaL(lo),ll,k,Rfrohner,RprimeL(lo),PhiNew)
c	call FindRprime(deltaL(lo),ll,k,ac,RprimeL(lo),PhiNew)

c	write(6,'("  - Phase shift from ac    ",22x," = ",E15.7)')Phi(lo)
c	write(6,'("  - Phase shift from OMC   ",22x," = ",E15.7)')deltaL(lo)
c	write(6,'("  - Phase shift from Rprime",22x," = ",E15.7)')PhiNew(lo)
c	write(6,*)' '

c-----------------------------------------------------------------------------
c Rprime from OMC	
c-----------------------------------------------------------------------------

c	write(6,'("  Hard sphere radius from OMC  ...")')

c	write(6,'("  - Rprime from OMC        ",22x," = ",E15.7)')RprimeL(lo)
c       write(6,'("  - Rprime from Frohner    ",22x," = ",E15.7)')Rfrohner
c	write(6,*)' '

c-----------------------------------------------------------------------------
c R-matrix parameters with hard sphere corrections	
c-----------------------------------------------------------------------------

c	write(6,'("  R-matrix elements with hard sphere correction  ...")')

	call PenF(P,Shift,ll,k,RprimeL(lo))
	call PhiF(P,Shift,Phi,ll)
	call Rinfiny(Rinf(lo),P(lo),Phi(lo),ReC(lo),ImC(lo))
	call polefunction(pole(lo),P(lo),Phi(lo),ReC(lo),ImC(lo))
	call pole2Sfunction(Str(lo),pole(lo),P(0),Elab)

c	write(6,'("  - Strenght function ",27x," = ",f8.4)')Str(lo)
c	write(6,'("  - Distant parameter ",27x," = ",f8.4)')Rinf(lo)
c	write(6,*)' '

c-----------------------------------------------------------------------------
c Breakdown strength function
c
c Tcoef : total transmission calculated from scattering matrix elements
c Tlcoef: transmission calculated by Optman 
c	
c-----------------------------------------------------------------------------

        call correctionPole(Tcoef(lo),Tlcoef(lo),
     +  pole(lo),pcomp(lo),pdir(lo),P(lo),Rinf(lo),Elab)
	call pole2Sfunction(Scomp(lo),pcomp(lo),P(0),Elab)
	call pole2Sfunction(Sdir(lo),pdir(lo),P(0),Elab)

c	write(6,'("  Pole function contributions  ...")')

c	write(6,'("  - Pole function total        ",18x," = ",f8.4)')pole(lo)
c	write(6,'("  - Pole function compound     ",18x," = ",f8.4)')pcomp(lo)
c	write(6,'("  - Pole function direct       ",18x," = ",f8.4)')pdir(lo)
c	write(6,*)' '

c	write(6,'("  Strenght function contributions  ...")')

c	write(6,'("  - Strenght function total    ",18x," = ",f8.4)')Str(lo)
c	write(6,'("  - Strenght function compound ",18x," = ",f8.4)')Scomp(lo)
c	write(6,'("  - Strenght function direct   ",18x," = ",f8.4)')Sdir(lo)
c	write(6,*)' '

c Asymptotic formula	
c-----------------------------------------------------------------------------

c	write(6,'("  Asymptotic formula ...")')
	
c        mod=ReC(lo)**2+ImC(lo)**2
c	call Asymptoticpolefunction(pole(lo),P(lo),Phi(lo),mod,ImC(lo))
c	call pole2Sfunction(Str(lo),pole(lo),P(0),Elab)

c	write(6,'("  - Asymptotic Strenght function ",16x," = ",f8.4)')Str(lo)	
 	
	
c Final Value
c
c STRt      : "total" strenght function form SPRT and ESW
c STRc      : compound contribution
c STRd      : direct contribution
c PHIomc    : phase shift from omc
c PHIsqw    : phase shift from ESW
c RPRIME    : effective radius 
c RINFINITY : distant level parameter
c 	
c-----------------------------------------------------------------------------
      
        STRt=Str(lo)*1.0e-4
        STRc=Scomp(lo)*1.0e-4
        STRd=Sdir(lo)*1.0e-4
        PHIomc=deltaL(lo)
        PHIsqw=PhiNew(lo)
        RPRIME=RprimeL(lo)
        RINFINITY=Rinf(lo)

	end 
 
c-----------------------------------------------------------------------------
c Subroutine channel radius 
c 
c use ENDF-6 convention	
c-----------------------------------------------------------------------------
	subroutine CR(ac,A)
	implicit none
	double precision :: ac,A,Mn
		
	Mn=1.008665
        ac=1.23*A**(1.0/3.0)+0.8
	
	return
	end

c-----------------------------------------------------------------------------
c Subroutine Wave number 
c 
c k is given in eV	
c-----------------------------------------------------------------------------
	subroutine WaveNumber(k,A,ELab)
	implicit none
	double precision :: k,A,Elab,cst
	
 	cst=2.1968E-04
        k=cst*A*sqrt(Elab)/(A+1.0)

	return
	end

c-----------------------------------------------------------------------------
c Subroutine penetration factor and level shift actor	
c-----------------------------------------------------------------------------
	subroutine PenF(P,Shift,ll,k,ac)
	implicit none
	double precision    ::  P(0:100),Shift(0:100)
	double precision    ::  k,ac
	double precision    ::  ii
	double precision    ::  num,den
	double precision    ::  ll
        double precision    ::  Jl,Nl
	integer             ::  i
	
        P(0)=k*ac
	Shift(0)=0.0
	
	ii=0.0
	
c  Iterative calculations for  l>0
c-----------------------------------------------------------------------------

	do i=1,Int(ll)

	 ii=ii+1.0

	 num=P(i-1)*P(0)**2
	 den=P(i-1)**2+(ii - Shift(i-1))**2
	 P(i)=num/den

	 num=(ii - Shift(i-1))*P(0)**2
	 den=P(i-1)**2+(ii - Shift(i-1))**2
	 Shift(i)=num/den-ii

	enddo
	
c  Correct calculations for l=1,2,3
c-----------------------------------------------------------------------------

	if(ll.ge.1.0.and.ll.le.3)then
         i=INT(ll)
         call Bessel(P(0),Jl,i)
         call Neumann(P(0),Nl,i)
         P(i)=1.0/(P(0)*Jl**2+P(0)*Nl**2)
        endif

c  Test for  l=3
c-----------------------------------------------------------------------------
c	if(Int(ll).eq.3)then
c	 num=P(0)**7
c	 den=225.0+45.0*P(0)**2+6.0*P(0)**4+P(0)**6
c	 P(ll)=num/den
c	endif
	
        return
	end


c-----------------------------------------------------------------------------
c Subroutine potential phase shift factor
c Value for l=0 to 4 :	
c	Phi(0)=r
c	Phi(1)=r-ATan(r)
c	Phi(2)=r-ATan(3.0*r/(3.0-r**2))
c	Phi(3)=r-ATan(r*(15.0-r**2)/(15.0-6*r**2))
c	Phi(4)=r-ATan(r*(105.0-10.0*r**2)/(105.0-45.0*r**2+r**4))
c-----------------------------------------------------------------------------
	subroutine PhiF(P,Shift,Phi,ll)
	implicit none
	double precision    ::  P(0:100),Shift(0:100),Phi(0:100)
	double precision    ::  ii
	double precision    ::  ll
        double precision    ::  Jl,Nl
	integer             ::  i
	
	Phi(0)=P(0)
	
	ii=0.0
	
c  Iterative calculatioons for  l>0
c-----------------------------------------------------------------------------

	do i=1,Int(ll)
	 ii=ii+1.0
	 Phi(i)=Phi(i-1) - ATan(P(i-1)/(ii - Shift(i-1)))
	enddo

c  Correct calculations for l=1,2,3
c-----------------------------------------------------------------------------

	if(ll.ge.1.0.and.ll.le.3)then
         i=INT(ll)
         call Bessel(Phi(0),Jl,i)
         call Neumann(Phi(0),Nl,i)
         Phi(i)=-1.0*Atan(Jl/Nl)
        endif
	
c  Test for  l=3
c-----------------------------------------------------------------------------
c	if(Int(ll).eq.3)then
c	Phi(ll)=P(0)-Atan(P(0)*(15.0-P(0)**2)/(15.0-6.0*P(0)**2))	
c	endif 
	
        return
	end

c-----------------------------------------------------------------------------
c Subroutine Bessel
c 
c Give Spherical Bessel of the 1st kind
c-----------------------------------------------------------------------------

	 subroutine Bessel(rho,Jl,i)
	 implicit none
	 double precision    ::    rho,Jl,C,S
	 integer             ::    i	 
	  
	 C=cos(rho)
	 S=sin(rho)
          
         if(i.eq.1)Jl=S/rho**2-C/rho
         if(i.eq.2)Jl=(3.0/rho**3-1.0/rho)*S-3.0*C/rho**2
         if(i.eq.3)Jl=(15.0/rho**4-6.0/rho**2)*S-(15.0/rho**3-1.0/rho)*C
		 
	 return
	 end

c-----------------------------------------------------------------------------
c Subroutine Neumann
c 
c Give Neumann function 
c-----------------------------------------------------------------------------

	 subroutine Neumann(rho,Nl,i)
	 implicit none
	 double precision    ::   rho,Nl,C,S
	 integer             ::   i	 
	  
	 C=cos(rho)
	 S=sin(rho)
          
         if(i.eq.1)Nl=-C/rho**2-S/rho
         if(i.eq.2)Nl=(-3.0/rho**3+1.0/rho)*C-3.0*S/rho**2
         if(i.eq.3)Nl=(-15./rho**4+6.0/rho**2)*C-(15./rho**3-1./rho)*S
		 
	 return
	 end

c-----------------------------------------------------------------------------
c Subroutine Rinfinity	
c----------------------------------------------------------------------- ------
	subroutine Rinfiny(Rinf,P,Phi,a,b)
	implicit none	
	double precision :: Rinf,P,Phi
	double precision :: a,b,C,S
	double precision :: num,den
	
	
	C=cos(2.0*Phi)
	S=sin(2.0*Phi)
			
	num=2.0*b*S-2.0*a*C-S
        den=P*(2.0*b-2.0*b*b-C+2.0*b*C+2.0*a*S-1.0-2.0*a*a)
	Rinf=num/den
	
	return
	end
		
c-----------------------------------------------------------------------------
c Subroutine polefunction	
c----------------------------------------------------------------------- ------
	subroutine polefunction(pole,P,Phi,a,b)
	implicit none	
	double precision :: Pole,P,Phi
	double precision :: a,b,S,C
	double precision :: num,den
	double precision :: Pi
	
	Pi=3.1415926
	
	C=cos(2.0*Phi)
	S=sin(2.0*Phi)
		
	num=2.0*(a*a-b+b*b)
	den=P*Pi*(2.0*b-2.0*b*b-C+2.0*b*C+2.0*a*S-1.0-2.0*a*a)
	pole=num/den
		
	return
	end

c-----------------------------------------------------------------------------
c Subroutine pole2Sfunction	
c----------------------------------------------------------------------- ------
	subroutine pole2Sfunction(S,pole,P,x)
	implicit none
	double precision :: S,pole,P,x
	
        S=2.0*pole*P/sqrt(x)
        S=S*1.0E+04
		
	return
	end

c-----------------------------------------------------------------------------
c Subroutine Shapelastic
c----------------------------------------------------------------------- ------
	subroutine shape(XSpot,P,Phi,R,pole,k)
	implicit none
	double precision ::  XSpot,P,Phi,R,pole,k
	double precision ::  cst,C,S,Num
	double precision ::  Pi
	
	Pi=3.1415926
	
	C=cos(Phi)
	S=sin(Phi)

	Num=P*P*(R*R+Pi*Pi*pole*pole)
	cst=10.0*4.0*Pi/k**2
	XSpot=cst*(Num*C*C+S*(S-2.0*P*R*C))/(1.0+2.0*P*Pi*pole+Num)
	
	return
	end
	
	
c-----------------------------------------------------------------------------
c Subroutine Transmission coef
c----------------------------------------------------------------------- ------
       subroutine TransCoef(T,P,R,pole)
       implicit none
       double precision ::  T,P,R,pole
       double precision ::  Num,Den
       double precision ::  Pi
       
       Pi=3.1415926
	       
       num=4.0*P*Pi*pole
       den=1.0+2.0*P*Pi*pole+P*P*(R*R+Pi*Pi*pole*pole)
       T=num/den
       
       return
       end
		
c-----------------------------------------------------------------------------
c Subroutine Xs total
c----------------------------------------------------------------------- ------
	subroutine XStotal(XS,P,Phi,R,pole,k)
	implicit none
	double precision :: XS,P,Phi,R,pole,k
	double precision :: XSpot,Tc
	double precision :: cst
	double precision :: Pi
	
	Pi=3.1415926
	
	cst=10.0*Pi/k**2
	
	call shape(XSpot,P,Phi,R,pole,k)
	call TransCoef(Tc,P,R,pole)
		
	XS=cst*Tc+XSpot
   
	return
	end

c-----------------------------------------------------------------------------
c Subroutine Sfunction2pole	
c-----------------------------------------------------------------------------
	subroutine Sfunction2pole(S,pole,P,x)
	implicit none
	double precision    :: S,pole,P,x
	
	pole=S*sqrt(x)/(2.0*P)
	pole=pole*1.0E-4
	
	return
	end

c-----------------------------------------------------------------------------
c Subroutine Asymtoticpolefunction	
c----------------------------------------------------------------------- ------
	subroutine Asymptoticpolefunction(pole,P,Phi,a,b)
	implicit none	
	double precision ::  Pole,P,Phi
	double precision ::  a,b
	double precision ::  Pi
	
	Pi=3.1415926
	
	pole=(b-a)/(P*Pi)
		
	return
	end

c-----------------------------------------------------------------------------
c Subroutine Phase
c-----------------------------------------------------------------------------

	 subroutine Phase(Re,Im,delta)
	 implicit none
	 double precision  ::    Im,Re,delta
	 double precision  ::    ModS
	 double precision  ::    deltaC,deltaS
	 double precision  ::    sign
	 	 
         ModS=Re*Re+Im*Im
	 deltaC=ACOS(Re/sqrt(ModS)) 
	 deltaS=ASIN(Im/sqrt(ModS)) 
       	 sign=1.0
	 if(deltaC.le.0.0)sign=-1.0
	 if(deltaS.le.0.0)sign=-1.0
	 delta=-0.5*sign*abs(deltaC)
	 
	 return
	 end

c-----------------------------------------------------------------------------
c Subroutine FindRprime
c-----------------------------------------------------------------------------
	subroutine FindRprime(delta,ll,k,ac,Rprime,phi)
	implicit none
	double precision k,Ratio,a,b,ll,perturb
	double precision Rprime,ac,Rprimep
	double precision P(0:100),Shift(0:100)
	double precision delta,Phi(0:100),Phip(0:100)
	integer          lo
	integer          i 
		
	lo=INT(ll)
	Phi(lo)=0.0
        Phip(lo)=0.0
	i=0

c Initial radius
c-----------------------------------------------------------------------------

	Rprime=ac     

c Perturbation parameter
c-----------------------------------------------------------------------------

	perturb=0.01
		
 100    continue
        i=i+1

c 1st calculation
c-----------------------------------------------------------------------------
             
	call PenF(P,Shift,ll,k,Rprime)
	call PhiF(P,Shift,Phi,ll)
	
        if(delta.le.0.0)return
	if(phi(lo).le.0.0)return

	Ratio=abs(delta/phi(lo))

	if(i.gt.100000)return
	if(Ratio.ge.0.99999.and.Ratio.le.1.00001)return
	
c 2nd calculation 
c-----------------------------------------------------------------------------

	Rprimep=Rprime+Rprime*perturb

	call PenF(P,Shift,ll,k,Rprimep)
	call PhiF(P,Shift,Phip,ll)
        
        a=(phip(lo)-phi(lo))/(Rprime*perturb)
	b=phi(lo)-a*Rprime

	Rprime=(delta-b)/a
			
	goto 100
	 
        return
	end
	
c-----------------------------------------------------------------------------
c Subroutine correctionPole
c 
c Take into account direct reaction
c p0: before correction
c p1: after correction
c-----------------------------------------------------------------------------

	 subroutine correctionPole(Ttot,Tc,ptot,pcomp,pdir,P,R,Elab)
	 implicit none
	 double precision    P,R,Elab,Pi,Ttot,Tc,ptot,pcomp,pdir
	 
	 Pi=3.1415926
	
	 pdir=(Ttot-Tc)*((1.0+Pi*P*(ptot))**2+(P*R)**2)/(4.0*Pi*P)
	 
         pcomp=ptot-pdir	         

	 return
	 end
