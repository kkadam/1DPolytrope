program ana
  implicit none
  double precision :: pi
  double precision, allocatable :: rho(:), x(:)
  integer :: n,i
  
  
  n=101
  pi =3.14159265359
  
  allocate (rho(n))
  allocate (x(n))
  
  do i=1,n
    x(i)=i*pi/n
    rho(i)=sin(x(i))/x(i)
  enddo
  
  
  open(unit=14,file='ana_n=1')	
    do i=1,n
      write(14,*) rho(i)
    enddo
   close(14)
  print*, "File ana_n=1 printed."

end program ana
