program scf
  implicit none
  double precision :: a, G, K, gamma, deltac,xl,xr, C
  double precision, allocatable :: rho(:), phi(:),x(:),sumx(:),sumrho(:), &
                                   grav(:),sumgrav(:), preciphi(:),enth(:)
  integer :: n,i,R,conv 
  double precision :: h,polytropic_n,enth_max,c_prev


!!!!Define grid size  
  polytropic_n = 1
  n=501                      !!!Change
  gamma=   1+1/polytropic_n !6/5.0!
!!!!Define grid boundery
  xl=0d0
  xr=1
  
!!!!Stepsize    
  h=(xr-xl)/(n-1)!!!!!!Include boundery points in array


!!!!Allocate arrays  
  allocate (rho(n))
  allocate (phi(n))
  allocate (x(n))
  allocate (sumx(n))
  allocate (sumrho(n))
  allocate (grav(n))
  allocate (sumgrav(n))
  allocate (preciphi(n))
  allocate (enth(n))
  
  

!!!!Assume rho  
  do i=1,n
    rho(i)= (1-x(i))
  end do  

!!!!Iterate Till Convergence  
  deltac=1
  conv=0
  c=1
  
  do while (deltac.gt.1d-4)	
    conv=conv+1
    c_prev=c 

!!Poisson Solve
    call poisson_solve(rho,phi)	
!    call mass_solve(rho,phi)
!!Get rho	
    
    c=phi(n)
    
    do i=1,n
      enth(i)= c -phi(i)
    end do 
    
    enth_max=enth(1)
      
    do i=1,n
      rho(i)= (enth(i)/enth_max)**polytropic_n   
    end do
    
!Normalize
	
    do i=2,n
      rho(i)=rho(i)/rho(1)
    end do  	  
     rho(1)=1.0
   
    deltac=abs((c-c_prev)/c_prev)
    print*, "Iteration", conv
    print*, "deltac",deltac
    
  end do

!  print*,"rho"
!  do i=1,n  
!    print*, rho(i)
!  end do 
  
  


  
  print*, "After", conv, "iterations!!"
  print*, "gamma=",gamma
  print*,  "K=", K
  
    open(unit=14,file='rho')	 !!!Change
    do i=1,n
      write(14,*) rho(i)
    enddo
   close(14)
  print*, "File rho printed."      !!!Change
  
  
end program scf


subroutine poisson_solve(rho,phi)
  implicit none
  double precision :: rho(501), phi(501)             !!!Change
  double precision :: a, G, K, gamma, deltac,xl,xr, C
  double precision, allocatable :: x(:),dphi(:)
  integer :: n,i,test, count
  double precision :: h  

  
  n=501       !!!Change

  xl=0d0
  xr=1
  
!!!!Stepsize    
  h=(xr-xl)/(n-1)!!!!!!Include boundery points in array


!!!!Allocate arrays  

  allocate (x(n))
  allocate (dphi(n))
 
  do i=1,n
    dphi(i)=0
  end do
  
  do i=1,n
    x(i)=(i-1)*h
  end do


  !!!!First Poisson solved, trial phi
  phi(1)=0
  phi(2)=0
  phi(n)=0


  do i=1,n-1
    phi(n)=phi(n)+h*rho(i)/((n-i)*h)
  end do
  
!  print*,phi(n)
 
  do i=3,n-1
!    phi(i+1)=(4*3.14*h**2*x(i)**2*rho(i)+(x(i)-h/2)**2*(phi(i)-phi(i-1))+(x(i)+h/2)**2*phi(i))/(x(i)+h/2)**2
     phi(i)= phi(n)/(1-h)*x(i)
  end do
  

  
  do i=1,n
    dphi(i)=phi(i)
  end do
  
  !!!!Poisson Relaxation
  test=1
  count=0
  do while (test==1)
    
    do i=2,n-1
      phi(i)=(phi(i+1)*(x(i)+h/2)**2+phi(i-1)*(x(i)-h/2)**2-4*3.14*h**2*x(i)**2*rho(i))/((x(i)+h/2)**2+(x(i)-h/2)**2)
    end do
    
    phi(1)=phi(2)
    
    test=0
    do i=1,n

      dphi(i)=abs(dphi(i)-phi(i))
      if (dphi(i).gt.1d-10) then
        test=1
      end if 
    end do  

    do i=1,n
      dphi(i)=phi(i)
    end do  
  

  count =count +1
  end do
!  print*, "Internal Iterations",count
  
!  print*,"phi"
!  do i=1,n  
!    print*, phi(i)
!  end do 
!  print*,count

end subroutine poisson_solve







subroutine goisson_solve(rho,phi)
  implicit none
  double precision :: a, G, K, gamma, deltac,xl,xr, C
  double precision, allocatable :: x(:),dphi(:)
  integer :: n,i,test, count
  double precision :: h  
  double precision :: rho(2001), phi(2001)        !!!Change

  
  n=2001       !!!Change

  xl=0d0
  xr=1
  
!!!!Stepsize    
  h=(xr-xl)/(n-1)!!!!!!Include boundery points in array


!!!!Allocate arrays  

  allocate (x(n))
  allocate (dphi(n))
  
  do i=1,n
    dphi(i)=0
  end do
  
  do i=1,n
    x(i)=xl+(i-1)*h
  end do
  
  
  !!!!First Poisson solved, trial phi
  phi(1)=0
  phi(2)=0
  phi(n)=0
  
  
  do i=1,n-1
    phi(n)=phi(n)+h*rho(i)/((n-i)*h)
  end do
  
!  print*,phi(n)
  
  do i=3,n-1
!    phi(i+1)=(4*3.14*h**2*x(i)**2*rho(i)+(x(i)-h/2)**2*(phi(i)-phi(i-1))+(x(i)+h/2)**2*phi(i))/(x(i)+h/2)**2
     phi(i)= phi(n)/(1-h)*x(i)
  end do
  
  
  
  do i=1,n
    dphi(i)=phi(i)
  end do
  
  !!!!Poisson Relaxation
  test=1
  count=0
  do while (test==1)
    
    do i=2,n-1
      phi(i)=(phi(i+1)*(x(i)+h/2)**2+phi(i-1)*(x(i)-h/2)**2-4*3.14*h**2*x(i)**2*rho(i))/((x(i)+h/2)**2+(x(i)-h/2)**2)
    end do
    
    phi(1)=phi(2)
    
    test=0
    do i=1,n

      dphi(i)=abs(dphi(i)-phi(i))
      if (dphi(i).gt.1d-6) then
        test=1
      end if 
    end do  

    do i=1,n
      dphi(i)=phi(i)
    end do  
  

  count =count +1
  end do
  !print*, "Internal Iterations",count
  
!  print*,"phi"
!  do i=1,n  
!    print*, phi(i)
!  end do 
!  print*,count

end subroutine goisson_solve
  



subroutine mass_solve(rho,phi)
  implicit none
  double precision :: a, Pi, dR
  double precision :: rho(2001), phi(2001), MR(2001), gR(2001), x(2001)    !!!Change
  integer :: Ns, i
  
  
  Ns=2001          !!!Change
  
  a=1
  dR=a/(Ns-1)
  Pi=3.14
  
  do i=1, Ns
    x(i)=(i-1)*dR
  end do

  do i=1, Ns
    MR(i)=4*Pi*dR*Rho(i)*x(i)**2
  enddo
  do i=2, Ns
    MR(i)=MR(i-1)+MR(i)
  enddo

  gR(1)=0
  do i=2, Ns
    gR(i)= Mr(i)/x(i)**2
  enddo
  
  
  do i=1, Ns
    phi(i)=-dR*gR(i)
  enddo
  do i=2, Ns
    phi(i)=phi(i-1)+phi(i)
  enddo
  
  
end subroutine mass_solve  
  
  
  
