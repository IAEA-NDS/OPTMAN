! variation on rosenbrock function
subroutine r(m, n, x, fvec)
  implicit none
  integer m, n
  double precision x(n), fvec(m)
  
  x->XAD(i)
  
  call XISQ
  
  fvec = (/ x(1), 100*(x(2) - x(1)**2) /)
end subroutine r

! jacobian of rosenbrock function
subroutine j(m, n, x, fjac)
  implicit none
  integer m, n
  double precision x(n), fjac(m, n)
  fjac = 0
end subroutine j

! second direction derivative function
subroutine Avv(m, n, x, v, acc)
  implicit none
  integer m, n
  double precision x(n), v(n), acc(m)
  acc = 0
end subroutine Avv

! Dummy callback routine
subroutine callback(m,n,x,v,a,fvec,fjac,acc,lam,dtd,fvec_new,accepted,info)
  implicit none
  integer m, n, accepted, info
  double precision x(n), v(n), a(n), fvec(m), fjac(m,n), acc(m), lam, dtd(n,n)
  double precision fvec_new(m)
end subroutine callback 
    
subroutine levMarSearch 
  implicit none
  ! dimensions of problem
  m = 2
  n = 2

  ! starting point
  x = (/ -160d0, 100d0 /)

  analytic_jac = .false.!.true.
  analytic_Avv = .false.!.true.

  ! Irrelevant since we use analytic derivatives
  center_diff = .true.!.false. 
  h1 = 1.0d-8
  h2 = 1.0d-1

  ! damping matrix
  dtd = reshape( (/ 1d0, 0d0, 0d0, 1d0 /),  (/ 2, 2 /) )

  ! dynamically update dtd
  damp_mode = 1 

  ! stopping criterion
  maxiter = 1000
  maxfev = -1
  maxjev = -1
  maxaev = -1
  maxlam = -1.0d0
  minlam = -1.0d0
  ! this is the preferred stopping criterion for most problems in which m > n
  ! I recommend using 0.001.
  ! since m = n in this case we do not use it (indicated by a negative value)
  artol = -1.0d0 
  ! instead, we only use the Cgoal as a stopping criterion
  Cgoal = 1.0d-8
  gtol = -1.0d0
  xtol = -1.0d0
  xrtol = -1.0d0
  ftol = -1.0d0
  frtol = -1.0d0

  ! setup printing of algorithm status
  print_level = 1
  print_unit = 6

  ! flags to determine algorithm
  imethod = 11
  iaccel = 0
  ibold = 0
  ibroyden = 0
  initialfactor = 100d0
  factoraccept = 2d0
  factorreject = 3d0
  avmax = 0.75d0
  
  call geodesiclm(r, j, Avv, x, fvec, fjac, n, m, callback, info, &
       & analytic_jac, analytic_Avv, & 
       & center_diff, h1, h2, dtd, damp_mode, &
       & niters, nfev, njev, naev, &
       & maxiter, maxfev, maxjev, maxaev, maxlam, minlam, &
       & artol, Cgoal, gtol, xtol, xrtol, ftol, frtol, &
       & converged, print_level, print_unit, &
       & imethod, iaccel, ibold, ibroyden, &
       & initialfactor, factoraccept, factorreject, avmax)

 
end subroutine levMarSearch

    
    