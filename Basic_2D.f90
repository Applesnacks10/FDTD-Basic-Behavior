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
integer, parameter :: Nt= 5000


integer, parameter :: Ny=101,N_loc=Ny-1 !N_loc must equal Ny-1 for 1 proc
double precision, parameter :: y0=-50E-9,yM=50E-9

integer, parameter :: Nx=101
double precision, parameter :: x0=-50E-9,xM=50E-9

!
!~~~ Spatial and Temporal steps; Spatial Indexing ~~~!
!
double precision, parameter :: dx=(xM-x0)/(Nx-1),dy=(yM-y0)/(Ny-1)
double precision, parameter :: dt=dy/(2.0*c)
double precision x(Nx),xM2(Nx-1),y(N_loc),yM2(N_loc)

!
!~~~ eps for host media; abbreviated coefficients ~~~!
!
double precision, parameter :: eps_delectric=1.0
double precision, parameter :: dt_eps0=dt/eps0,dt_mu0=dt/mu0
!double precision !Insert den_ here -- Note: replaced den_ with 1/d_

!
!~~~ EM field components; Field Input ~~~!
!
double precision Ex(Nx-1,N_loc),Ey(Nx,N_loc),Hz(Nx-1,N_loc)
double precision Ex_inc(N_loc),Hz_inc(N_loc)
double precision pulse(Nt)
double precision aBH(4)
double precision, parameter :: tau=0.12d-15,E0=1.0,omega=ev_to_radsec*3.0
integer, parameter :: js = -1, is = -1 ! No pulse

!
!~~~ Loop Indices; time ~~~!
!
integer i,ii,j,jj,n,nn,k,a,b
double precision t

!
!~~~ Grid Return ~~~!
!
double precision, parameter:: y_return1 = y0, y_return2 = yM, x_return1 = x0, x_return2 = xM!-1E-9
integer, parameter :: Nreturn = 5           !xM or xM-1E-9, depending on x-size of field array.
integer n_return(Nreturn)
integer i_return1, i_return2, j_return1, j_return2
logical GR
 character(len = 9), parameter :: prefix = 'Snapshot-'
 character(len = 4), parameter :: suffix = '.dat'
 character(len = 2), parameter :: str_Ex = 'Ex', str_Hz = 'Hz', str_Ey = 'Ey'
 character(len = 50) filename, str_n

!
!~~~ MPI part ~~~!
!
integer, parameter :: myrank = 0
integer j_glob

!~~~ grid ~~~!
do i=1,Nx
 x(i)=x0+dx*(i-1)
enddo
do i=1,(Nx-1)
 xM2(i)=x0+dx*(i-1)+dx/2.0
enddo
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
  
  do i = 1,Nx
   if( x(i) >= x_return1 )then !lower x-bound
    i_return1 = i
    exit
   endif
  enddo
  do i = 1,Nx
   if( x(Nx-i+1) <= x_return2 )then !upper x-bound
    i_return2 = Nx-i+1
    exit
   endif
  enddo
 endif
  
! j_return1 = 1
! j_return2 = N_loc
! i_return1 = 1
! i_return2 = Nx-1 
! itag6 = N_loc+2
! requisite_send = 1
! requisite_recv = 0
 n_return(1) = 200*2
 n_return(2) = 400*2
 n_return(3) = 800*2
 n_return(4) = 1200*2
 n_return(5) = 1600*2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
   !~~~ Initialize to Zero ~~~!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

Ex=0.0
Ey=0.0
Hz=0.0

Ex_inc=0.0
Hz_inc=0.0

pulse = 0.0

!   filename = 'Snapshot-Ex8.dat'
!   open(file=filename,position = 'append',unit=a+10)
!    do j = j_return1,j_return2
!     write(a+10,*) Ex(i_return1:i_return2,j)
!    enddo
!   close(unit=a+10)
  
!    open(file='Return_Grid_y.dat',position = 'append',unit=90)
!     do j = j_return1,j_return2
!      write(90,*) (y(j), i = i_return1,i_return2)
!     enddo
!    close(unit=90)
!   
!    open(file='Return_Grid_x.dat',position = 'append',unit=91)
!     do j = j_return1,j_return2
!       write(91,*) (x(i), i = i_return1,i_return2) 
!     enddo
!    close(unit=91)

!~~~ Pulse ~~~!

aBH(1)=0.353222222
aBH(2)=-0.488
aBH(3)=0.145
aBH(4)=-0.010222222

do n=1,Nt
 t=dt*dble(n)
 if(t<=tau)then
   pulse(n)=E0*cos(omega*t)*( &
                 aBH(1)+ &
				 aBH(2)*cos(2.0*pi*t/tau)+ &
				 aBH(3)*cos(2.0*pi*2.0*t/tau)+ &
				 aBH(4)*cos(2.0*pi*3.0*t/tau))
  else
   pulse(n)=0.0
 endif
enddo

do n=1,Nt
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Hz ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

!~~~ total ~~~!
do i=1,Nx-1
 do j=1,N_loc-1
 
  Hz(i,j)=Hz(i,j)+dt_mu0*((Ey(i,j)-Ey(i+1,j))/dx + &
			               (Ex(i,j+1)-Ex(i,j))/dy)
 enddo
enddo

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Ex ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

!~~~ total ~~~! 
do i=1,Nx-1  
 do j=2,N_loc
 
  Ex(i,j)=Ex(i,j)+dt_eps0*(Hz(i,j)-Hz(i,j-1))/dy
  if(j == js)then !js = -1
   Ex(i,j) = Ex(i,j) + pulse(n)
  endif
  
 enddo
enddo

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Ey ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

i=1
do j=1,N_loc
 
 Ey(i,j)=Ey(i,j)+dt_eps0*(Hz(Nx-1,j)-Hz(i,j))/dx !Is this the periodic BC??

enddo
 
do i=2,Nx-1
 do j=1,N_loc
  
  Ey(i,j)=Ey(i,j)+dt_eps0*(Hz(i-1,j)-Hz(i,j))/dx
  if(i == is.and.j=/N_loc)then ! laser pulse. At N_loc, Hz is not updated.
   Ey(i,j) = Ey(i,j) + pulse(n)
  endif
 
 enddo
enddo

i=Nx
do j=1,N_loc
 
 Ey(i,j)=Ey(i,j)+dt_eps0*(Hz(i-1,j)-Hz(1,j))/dx !Nx is the Correspondent to i = 1

enddo

! Ey(1) = Ey(Nx). i = 1 and i = Nx are interchangable 
! (picture a grid wrapped into a cylinder, where the endpoints become coincident)
! ((the coincidence of the endpoints is an analog to the reciprocity of index)) 


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
   filename = prefix//str_Ey//trim(adjustl(str_n))//suffix
   open(file=trim(adjustl(filename)),position = 'append',unit=a+10)
    do j = j_return1,j_return2
     write(a+10,*) Ey(i_return1:i_return2,j)
    enddo
   close(unit=a+10)
  
!   if(a == 1)then 
!    open(file='Return_Grid_y0.dat',position = 'append',unit=90)
!     do j = j_return1,j_return2
!      write(90,*) (y(j), i = i_return1,i_return2)
!     enddo
!    close(unit=90)
!   
!    open(file='Return_Grid_x0.dat',position = 'append',unit=91)
!     do j = j_return1,j_return2
!       write(91,*) (x(i), i = i_return1,i_return2) 
!     enddo
!    close(unit=91)
!   endif   

  endif
 enddo
endif !GR

enddo !Nt

end !main

