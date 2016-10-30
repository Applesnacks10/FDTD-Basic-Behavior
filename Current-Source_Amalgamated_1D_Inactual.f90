implicit none

!
!~~~ fundamental constants [all numbers are in SI units]~~~!
!
double precision, parameter :: pi=3.1415926535897932384626433832795D0,c=299792458.0D0
double precision, parameter :: mu0=4.0D-7*pi,eps0=1.0D0/(c*c*mu0)
double precision, parameter :: ev_to_radsec=2.0*pi*2.4180e14

!
!~~~ number of grid points & time steps ~~~!
!
integer, parameter :: Nt= 1000

integer, parameter :: Ny=501,N_loc=Ny-1 !N_loc must equal Ny-1 for 1 proc
double precision, parameter :: y0=-250E-9,yM=250E-9

!
!~~~ Spatial and Temporal steps; Spatial Indexing ~~~!
!
double precision, parameter :: dy=(yM-y0)/(Ny-1)
double precision, parameter :: dt=dy/(2.0*c)
double precision y(N_loc),yM2(N_loc)

!
!~~~ eps for host media; abbreviated coefficients ~~~!
!
double precision, parameter :: eps_delectric=1.0
double precision, parameter :: dt_eps0=dt/eps0,dt_mu0=dt/mu0

!
!~~~ EM field components ~~~!
!
double precision Ex(N_loc),Hz(N_loc)
double precision Ex_inc(N_loc),Hz_inc(N_loc)

!
!~~~ Field Input ~~~!
!
integer, parameter :: js = N_loc/2 !Place in the midpoint in order to confirm dual-propagation
double precision, parameter :: tau=0.36d-15,E0=1.0,omega=ev_to_radsec*3.0
double precision aBH(4)
double precision Jx(Nt)

!
!~~~ Loop Indices; time ~~~!
!
integer i,ii,j,jj,n,nn,k,a,b
double precision t

!
!~~~ Grid Return ~~~!
!
double precision, parameter:: y_return1 = y0, y_return2 = yM
integer, parameter :: Nreturn = 5
integer n_return(Nreturn)
integer j_return1, j_return2
logical GR
 character(len = 9), parameter :: prefix = 'Snapshot-'
 character(len = 4), parameter :: suffix = '.dat'
 character(len = 2), parameter :: str_Ex = 'Ex', str_Hz = 'Hz'
 character(len = 50) filename, str_n

!
!~~~ MPI part ~~~!
!
integer, parameter :: myrank = 0
integer j_glob

!~~~ grid ~~~!

do j=1,N_loc
 j_glob=myrank*N_loc+j
 y(j)=y0+dy*(j_glob-1)
 yM2(j)=y0+dy*(j_glob-1)+dy/2.0
enddo

!~~~ Grid Return ~~~!


 if( (y(1)<=y_return2).and.(y(N_loc)>=y_return1) )then
  GR = .true. !Processor has at least one grid-point within the return-zone
 else
  GR = .false.
 endif
 
 if(GR)then !Assign local return boundaries
  do j = 1,N_loc
   if( y(j) >= y_return1 )then !lower y-bound
    j_return1 = j
    exit
   endif
  enddo
  do j = 1,N_loc
   if( y(N_loc-j+1) <= y_return2 )then !upper y-bound
    j_return2 = N_loc-j+1
    exit
   endif
  enddo
 endif
 
! j_return1 = 1
! j_return2 = N_loc 

 n_return(1) = 50
 n_return(2) = 100
 n_return(3) = 200
 n_return(4) = 400
 n_return(5) = 800

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
   !~~~ Initialize to Zero ~~~!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

Ex=0.0
Hz=0.0

Ex_inc=0.0
Hz_inc=0.0

Jx = 0.0

!~~~ Source ~~~!
do n=1,Nt
 t = dt*dble(n)
 if(t <= tau)then
  Jx(n)=-dy/dt_eps0*E0*cos(omega*t)*( &
                  aBH(1)+ &
		  aBH(2)*cos(2.0*pi*t/tau)+ &
		  aBH(3)*cos(2.0*pi*2.0*t/tau)+ &
		  aBH(4)*cos(2.0*pi*3.0*t/tau))
 endif
enddo

do n=1,Nt
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Hz ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

!~~~ total ~~~!
do j=1,N_loc-1
 
 Hz(j)=Hz(j)+dt_mu0*(Ex(j+1)-Ex(j))/dy
 
enddo

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Ex ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

!~~~ total ~~~!  
do j=2,N_loc
 
 Ex(j)=Ex(j)+dt_eps0*(Hz(j)-Hz(j-1))/dy
 if(j == js)then !add current source
  Ex(j) = Ex(j) - dt_eps0/dy*(Jx(n))
 endif
  
enddo

!--------------------------------------------------------------------------!
!~~~~~~~~~~~~~~~~~~~~~~~~==========================~~~~~~~~~~~~~~~~~~~~~~~~!
!--------------------------------------------------------------------------!
!~~~~~~~~~~~~~~~~~~~~~~~~        Grid Return         ~~~~~~~~~~~~~~~~~~~~~~!
!--------------------------------------------------------------------------!
!~~~~~~~~~~~~~~~~~~~~~~~~==========================~~~~~~~~~~~~~~~~~~~~~~~~!
!--------------------------------------------------------------------------!

if(Nreturn > 0.and.GR)then 
 do a = 1,Nreturn
  if(n == n_return(a))then
   
   write(str_n,*) n
   filename = str_Ex//trim(adjustl(str_n))//suffix
!   filename = prefix//filename
   open(file=trim(adjustl(filename)),position = 'append',unit=a+10)
    do j = j_return1,j_return2
     write(a+10,*) Ex(j)
    enddo
   close(unit=a+10)
  
!   if(a == 1)then 
!    open(file='Return_Grid_y0.dat',position = 'append',unit=90)
!     do j = j_return1,j_return2
!      write(90,*) y(j)
!     enddo
!    close(unit=90)
!   

  endif
 enddo
endif !GR

enddo !Nt

end !main

