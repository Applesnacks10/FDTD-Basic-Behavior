program Convergence_Vacuum
implicit none
include 'mpif.h'

integer, parameter :: Nr = 1
integer, parameter, dimension(Nr) :: res_array = (/1/)
integer, parameter, dimension(2) :: pml_add = (/0,1/)
double precision :: Convergence(Nr,2), Rel_error(Nr)
integer :: a,b !loop variables

!-------------------------------------------------------------------
!------------------------ Global Parameters ------------------------
!-------------------------------------------------------------------

double precision, parameter :: length_add = 9.0E-9 
!
!~~~ fundamental constants [all numbers are in SI units]~~~!
!
double precision, parameter :: pi=3.1415926535897932384626433832795D0,c=299792458.0D0
double precision, parameter :: mu0=4.0D-7*pi,eps0=1.0D0/(c*c*mu0)
double precision, parameter :: h=1.054571628D-34
double complex, parameter :: Im=(0.0,1.0)
double precision, parameter :: ev_to_radsec=2.0*pi*2.4180e14

double precision, parameter :: eps_delectric=1.0
!
!~~~ Source ~~~!
!

double precision, parameter :: H0 = 1.0, wavelength = length_add
double precision, parameter :: omega = 2*pi*c/wavelength
double precision, parameter :: x_source = 0.0, y_source = 0.0
double precision, parameter :: x_detect = x_source, y_detect = y_source

! Drude Specification

!double precision, parameter :: eps_r=8.926,omegaD=ev_to_radsec*11.585,GammaD=ev_to_radsec*0.203
!double precision, parameter :: A1=(2.0-GammaD*dt)/(2.0+GammaD*dt),A2=eps0*omegaD*omegaD*dt/(2.0+GammaD*dt)
!double precision, parameter :: C1=(eps_r*eps0/dt-0.5*A2)/(eps_r*eps0/dt+0.5*A2)
!double precision, parameter :: C3=1.0/(eps_r*eps0/dt+0.5*A2)
!double precision, parameter :: C4=0.5*(A1+1.0)/(eps_r*eps0/dt+0.5*A2)

!CPML Scaling
integer, parameter :: m = 3, ma = 1
double precision, parameter :: alphaCPML = 0.05, kappaCPML = 5.0

!
!~~~ Indices ~~~!
!
integer i,ii,j,jj,n,nn,k
double precision t

!
!~~~ MPI part ~~~!
!
integer ierr,nprocs,myrank,j_glob,mfdtd,n1
integer :: istatus(MPI_STATUS_SIZE)
integer itag,ireq,itag1,itag2,itag3,itag4,itag5,itag6

!-------------------------------------------------------------------
!----------------- Minimum Resolution Variables --------------------
!-------------------------------------------------------------------

 double precision, parameter :: y0_min = -100E-9, yM_min = 100E-9
 double precision, parameter :: x0_min = -320E-9, xM_min = 320E-9

 double precision, parameter ::	               &
      dx_max = 1.0E-9 , dy_max = dx_max 
            
 double precision, parameter ::                &
      dt_max = dx_max/(2.0*c)

   INTEGER, parameter :: &
      Nt_min = 10000,                     & 
      Nx_min = 200 + 1,                  & 
      Ny_min = 640 + 1,                  &
      N_loc_min = 20
      
   INTEGER, parameter :: &
      npml_min = 9 + 1
      
!-------------------------------------------------------------------
!----------------- Maximum Resolution Variables --------------------
!-------------------------------------------------------------------

 double precision, parameter :: y0_max = y0_min - length_add, yM_max = yM_min + length_add
 double precision, parameter :: x0_max = x0_min - length_add, xM_max = xM_min + length_add 
      
  double precision, parameter :: &
      dx_min = (dx_max)/res_array(Nr), dy_min =(dy_max)/res_array(Nr)
            
 double precision, parameter :: &
      dt_min = (dt_max)/res_array(Nr)

   INTEGER, parameter :: &
      Nt_max = res_array(Nr)*Nt_min,                  &                
      Nx_max = res_array(Nr)*(Nx_min-1) + 2*length_add/dx_min + 1,                  &
      Ny_max = res_array(Nr)*(Ny_min-1) + 2*length_add/dx_min + 1,                  &
      N_loc_max = res_array(Nr)*N_loc_min + length_add/dx_min !Interior processors inherit (npml_add) unused y-spaces

      
   INTEGER, parameter :: &
      npml_max = 9*res_array(Nr) + length_add/dx_min + 1

!-------------------------------------------------------------------
!------------------------ Field and CPML Arrays --------------------
!-------------------------------------------------------------------

double precision :: Hz(Nx_max-1, N_loc_max), Hz_get(Nx_max-1), Hz_send(Nx_max-1)

double precision :: Ex(Nx_max-1, N_loc_max), Ex_get(Nx_max-1), Ex_send(Nx_max-1)

double precision :: Ey(Nx_max, N_loc_max)

!Physical Grid

double precision :: x(Nx_max), xM2(Nx_max-1), y(N_loc_max), yM2(N_loc_max)

!Source

double precision :: pulse(Nt_max)

!PML

double precision :: psi_Hzx_1(npml_max-1,N_loc_max)

double precision :: psi_Hzx_2(npml_max-1,N_loc_max)

double precision :: psi_Eyx_1(npml_max,N_loc_max)

double precision :: psi_Eyx_2(npml_max,N_loc_max)

double precision :: psi_Hzy_1(Nx_max-1,npml_max-1)

double precision :: psi_Hzy_2(Nx_max-1,npml_max-1)                          

double precision :: psi_Exy_1(Nx_max-1,npml_max)

double precision :: psi_Exy_2(Nx_max-1,npml_max)                        

double precision, dimension(npml_max-1) :: bh_x, ch_x, alphah_x, sigh_x, kappah_x

double precision, dimension(npml_max) :: be_x, ce_x, alphae_x, sige_x, kappae_x

double precision :: bh_y(npml_max-1), ch_y(npml_max-1), alphah_y(npml_max-1), sigh_y(npml_max-1), kappah_y(npml_max-1)

double precision :: be_y(npml_max), ce_y(npml_max), alphae_y(npml_max), sige_y(npml_max), kappae_y(npml_max)

!denominators for the update equations

double precision :: den_ex(Nx_max), den_hx(Nx_max)

double precision :: den_ey(N_loc_max), den_hy(N_loc_max)

!-------------------------!
 call MPI_INIT(ierr)
 call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
 call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      
!-----------------------------------------------------------------------
!--------------------------  Convergence Loop --------------------------
!-----------------------------------------------------------------------

 Convergence = 0.0

do a = 1,Nr
 do b = 1,2
 
  Convergence(a,b) = Vacuum_CPML()
  
 enddo! 2 cpml lengths
 
 Rel_error(a) = abs((Convergence(a,2) - Convergence(a,1))/Convergence(a,1))
 

 if(Rel_error(a) /= 0.0)then 
  open(file = 'Relative Errors.dat', position = 'append', unit = 28)
    write(28,*) res_array(a), Rel_error(a)
  close(unit = 28)
 endif
 
enddo! Nr resolutions

 call MPI_FINALIZE(ierr)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!------------------------ Internal Function ----------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

 contains
 
function Vacuum_CPML() result(P_sum)

   double precision :: P_sum

!-------------------------------------------------------------------
!----------------------- Variable Declaration ----------------------
!-------------------------------------------------------------------

   double precision :: dx, dy, x0, y0, xM, yM
   
   double precision :: dt

   integer :: Nt, Nx, Ny, N_loc
   
   integer :: npml
   
!!  ..................................
!!  Convergence Detection Zone
!   integer :: i_start, i_end, & 
!              j_start, j_end, 
!   integer ::                                    
!      isource, jsource
   double precision :: sigmaCPML
   double precision :: dt_eps0, dt_mu0
   
!
!~~~ Grid Return ~~~!
!
double precision, parameter :: y_return1 = y_source - 10*dx_max, y_return2 = y_source + 10*dx_max
double precision, parameter :: x_return1 = x_source - 10*dx_max, x_return2 = x_source + 10*dx_max 
integer, parameter :: Nreturn = 2
integer n_return(Nreturn)
integer i_return1, i_return2, j_return1, j_return2
logical GR
 character(len = 9), parameter :: prefix = 'Snapshot-'
 character(len = 4), parameter :: suffix = '.dat'
 character(len = 2), parameter :: str_Ex = 'Ex', str_Hz = 'Hz', str_Ey = 'Ey'
 character(len = 50) filename, str_n

!-------------------------------------------------------------------
!----------------------- Variable Assignment -----------------------
!-------------------------------------------------------------------

 dx = dx_max/res_array(a)
 dy = dy_max/res_array(a)
 
 x0 = x0_min - pml_add(b)*length_add/dx
 xM = xM_min + pml_add(b)*length_add/dx
 
 y0 = y0_min - pml_add(b)*length_add/dx
 yM = yM_min + pml_add(b)*length_add/dx
 
 dt = dx/(2.0*c)
 
 dt_eps0 = dt/eps0
 dt_mu0 = dt/mu0

 npml = res_array(a)*(npml_min-1) + pml_add(b)*length_add/dx + 1
 
 Nt = res_array(a)*Nt_min
 Nx = res_array(a)*(Nx_min-1) + 2*pml_add(b)*length_add/dx + 1
 Ny = res_array(a)*(Ny_min-1) + 2*pml_add(b)*length_add/dy + 1
 N_loc = res_array(a)*(N_loc_min)

if(myrank == 0 .or. myrank == nprocs-1)then !PML Processors take on all of the extra load since npml_add < nprocs
 N_loc = res_array(a)*N_loc_min + pml_add(b)*length_add/dx
endif
 
 sigmaCPML=0.8*(m+1)/(dx*(mu0/eps0*eps_delectric)**0.5)

! isource = (Nx-1)/2 + 1
! jsource = (Ny-1)/2 + 1
!
! i_start = isource
! j_start = jsource
! i_end   = i_start   
! j_end   = j_start

!
!~~~ Specify Source ~~~!
!

pulse=0.0
do n=1,Nt
 t=dt*dble(n-1)
 
 pulse(n) = H0*sin(omega*t)

enddo

!
!~~~ physical grid ~~~!
!

do i=1,Nx
 x(i)=x0+dx*(i-1)
enddo

do i=1,(Nx-1)
 xM2(i)=x0+dx*(i-1)+dx/2.0
enddo

do j=1,N_loc

 if(myrank == 0)then
  j_glob = j
 endif

 if(myrank > 0 .and. myrank < nprocs-1)then
  j_glob=myrank*N_loc + pml_add(b)*length_add/dx + j
 endif

 if(myrank == nprocs-1)then
  j_glob = myrank*(N_loc - pml_add(b)*length_add/dx) + pml_add(b)*length_add/dx + j
 endif

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

 nn = 24
 n_return(1) = Nt/100
 n_return(2) = Nt

!!FBx=.false.
!!FBy=.false.

!~~~ structure ~~~!
!do i=1,Nx-1
! do j=1,N_loc
!  if( &
!    ((y(j)>z1).and.(y(j)<z2).and.(xM2(i)<(-R)).and.(xM2(i)>(-slit_length/2.0))).or. &
!    ((y(j)>z1).and.(y(j)<z2).and.(xM2(i)>R).and.(xM2(i)<(slit_length/2.0))) &
!     )then
!    !!FBx(i,j)=.true.
!   else
!    !!FBx(i,j)=.false.
!  endif
! enddo
!enddo
!
!do i=1,Nx
! do j=1,N_loc
!  if( &
!    ((yM2(j)>z1).and.(yM2(j)<z2).and.(x(i)<(-R)).and.(x(i)>(-slit_length/2.0))).or. &
!    ((yM2(j)>z1).and.(yM2(j)<z2).and.(x(i)>R).and.(x(i)<(slit_length/2.0))) &
!     )then
!    !!FBy(i,j)=.true.
!   else
!    !!FBy(i,j)=.false.
!  endif
! enddo
!enddo

!!FBy = .false.
!!FBx = .false.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         !~~~ CPML vectors ~~~!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

!~~~ set CPML vectors ~~~!
do i=1,npml
 sige_x(i)=sigmaCPML*((npml-i)/(npml-1.0))**m
 alphae_x(i)=alphaCPML*((i-1.0)/(npml-1.0))**ma
 kappae_x(i)=1.0+(kappaCPML-1.0)*((npml-i)/(npml-1.0))**m
 be_x(i)=exp(-(sige_x(i)/kappae_x(i)+alphae_x(i))*dt/eps0)
 if( &
    (sige_x(i)==0.0).and. &
    (alphae_x(i)==0.0).and. & 
    (i==npml) &
   )then
   ce_x(i)=0.0
  else
   ce_x(i)=sige_x(i)*(be_x(i)-1.0)/(sige_x(i)+kappae_x(i)*alphae_x(i))/ kappae_x(i)
 endif
enddo

do i=1,npml-1
 sigh_x(i)=sigmaCPML*((npml-i-0.5)/(npml-1.0))**m
 alphah_x(i)=alphaCPML*((i-0.5)/(npml-1.0))**ma
 kappah_x(i)=1.0+(kappaCPML-1.0)*((npml-i-0.5)/(npml-1.0))**m
 bh_x(i)=exp(-(sigh_x(i)/kappah_x(i)+alphah_x(i))*dt/eps0)
 ch_x(i)=sigh_x(i)*(bh_x(i)-1.0)/(sigh_x(i)+kappah_x(i)*alphah_x(i))/kappah_x(i)
enddo

do j=1,npml
 sige_y(j)=sigmaCPML*((npml-j)/(npml-1.0))**m
 alphae_y(j)=alphaCPML*((j-1)/(npml-1.0))**ma
 kappae_y(j)=1.0+(kappaCPML-1.0)*((npml-j)/(npml-1.0))**m
 be_y(j)=exp(-(sige_y(j)/kappae_y(j)+alphae_y(j))*dt/eps0)
 if( &
    (sige_y(j)==0.0).and.&
    (alphae_y(j)==0.0).and. &
    (j==npml) &
   )then
   ce_y(j)=0.0
  else
   ce_y(j)=sige_y(j)*(be_y(j)-1.0)/(sige_y(j)+kappae_y(j)*alphae_y(j))/kappae_y(j)
 endif
enddo
   
do j=1,npml-1
 sigh_y(j)=sigmaCPML*((npml-j-0.5)/(npml-1.0))**m
 alphah_y(j)=alphaCPML*((j-0.5)/(npml-1.0))**ma
 kappah_y(j)=1.0+(kappaCPML-1.0)*((npml-j-0.5)/(npml-1.0))**m
 bh_y(j)=exp(-(sigh_y(j)/kappah_y(j)+alphah_y(j))*dt/eps0)
 ch_y(j)=sigh_y(j)*(bh_y(j)-1.0)/(sigh_y(j)+kappah_y(j)*alphah_y(j))/kappah_y(j)
enddo

den_hy=1.0/dy
if(myrank==0)then
  do j=1,N_loc
   if(j<=(npml-1))then
    den_hy(j)=1.0/(kappah_y(j)*dy)
   endif
  enddo
 elseif(myrank==(nprocs-1))then
  jj=npml-1
  do j=1,(N_loc-1)
   if(j>=(N_loc+1-npml))then
     den_hy(j)=1.0/(kappah_y(jj)*dy)
     jj=jj-1
   endif
  enddo
endif

den_ey=1.0/dy
if(myrank==0)then
  do j=1,N_loc
   if(j<=npml)then
    den_ey(j)=1.0/(kappae_y(j)*dy)
   endif
  enddo
 elseif(myrank==(nprocs-1))then
  jj=npml
  do j=1,(N_loc-1)
   if(j>=(N_loc+1-npml))then
     den_ey(j)=1.0/(kappae_y(jj)*dy)
     jj=jj-1
   endif
  enddo
endif

ii=npml-1
do i=1,Nx-1
 if(i<=(npml-1))then
   den_hx(i)=1.0/(kappah_x(i)*dx)
  elseif(i>=(Nx+1-npml))then
   den_hx(i)=1.0/(kappah_x(ii)*dx)
   ii=ii-1
  else
   den_hx(i)=1.0/dx
 endif
enddo

ii=npml
do i=1,Nx-1
 if(i<=npml)then
   den_ex(i)=1.0/(kappae_x(i)*dx)
  elseif (i>=(Nx+1-npml))then
   den_ex(i)=1.0/(kappae_x(ii)*dx)
   ii=ii-1
  else
   den_ex(i)=1.0/dx
 endif
enddo

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         !~~~ end of CPML ~~~!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

Ex=0.0
Ey=0.0
Hz=0.0

!PDx=0.0
!PDy=0.0

psi_Hzy_1=0.0
psi_Exy_1=0.0
psi_Hzy_2=0.0
psi_Exy_2=0.0
psi_Hzx_1=0.0
psi_Eyx_1=0.0
psi_Hzx_2=0.0
psi_Eyx_2=0.0

Ex_get=0.0
Ex_send=0.0
Hz_get=0.0
Hz_send=0.0

  if(myrank == 0)then
   write(*,*)"res: ", res_array(a)
   write(*,*)"pml_add: ", pml_add(b)
   write(*,*)"begin time-stepping"
  endif

do n=1,Nt
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Hz ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
 do j=1,(nprocs-1) !Ex send/recieve
  itag=j
  itag1=j+1
  if(myrank==j)then
    do i=1,(Nx-1)
     Ex_send(i)=Ex(i,1)
    enddo    
    call MPI_SEND(Ex_send,(Nx-1),MPI_doUBLE_PRECISION,(j-1),itag,MPI_COMM_WORLD,ierr)
   elseif(myrank==(j-1))then
    call MPI_RECV(Ex_get,(Nx-1),MPI_doUBLE_PRECISION,j,itag,MPI_COMM_WORLD,istatus,ierr)
  endif!send_recv
 enddo!nprocs

if(myrank==0)then !rank=0, here PML in y-direction is only applied to the bottom part
 do i=1,Nx-1
  do j=1,N_loc-1
   Hz(i,j)=Hz(i,j)+dt_mu0*((Ey(i,j)-Ey(i+1,j))*den_hx(i)+ &
			               (Ex(i,j+1)-Ex(i,j))*den_hy(j))
  enddo
    
  j=N_loc
   Hz(i,j)=Hz(i,j)+dt_mu0*((Ey(i,j)-Ey(i+1,j))*den_hx(i)+ &
			               (Ex_get(i)-Ex(i,j))*den_hy(j))
 enddo

 do j=1,N_loc
!  PML for left Hz, x-direction
  do i=1,npml-1
   psi_Hzx_1(i,j)=bh_x(i)*psi_Hzx_1(i,j)+ch_x(i)*(Ey(i,j)-Ey(i+1,j))/dx
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzx_1(i,j)
  enddo
!  PML for right Hz, x-direction
  ii=npml-1
  do i=Nx+1-npml,Nx-1
   psi_Hzx_2(ii,j)=bh_x(ii)*psi_Hzx_2(ii,j)+ch_x(ii)*(Ey(i,j)-Ey(i+1,j))/dx
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzx_2(ii,j)
   ii=ii-1
  enddo
 enddo

 do i=1,Nx-1
 !  PML for bottom Hz [bottom only here! since myrank=0], y-direction
  do j=1,npml-1
   psi_Hzy_1(i,j)=bh_y(j)*psi_Hzy_1(i,j)+ch_y(j)*(Ex(i,j+1)-Ex(i,j))/dy
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzy_1(i,j)
  enddo
 enddo
 
endif !myrank == 0

if((myrank>0).and.(myrank<(nprocs-1)))then !no PML for y-direction here
 do i=1,Nx-1
  do j=1,N_loc-1
   Hz(i,j)=Hz(i,j)+dt_mu0*((Ey(i,j)-Ey(i+1,j))*den_hx(i)+ &
			               (Ex(i,j+1)-Ex(i,j))*den_hy(j))
  enddo
    
  j=N_loc
   Hz(i,j)=Hz(i,j)+dt_mu0*((Ey(i,j)-Ey(i+1,j))*den_hx(i)+ &
			               (Ex_get(i)-Ex(i,j))*den_hy(j))
 enddo
 
!
!~~~ Source and Detection~~~!
!

if(myrank == (nprocs)/2)then
! do j = 1,N_loc
!  do i = 1,Nx
!   if(x(i) == x_source .and. y(j) == y_source)then
   i = (Nx-1)/2
   j = N_loc
    Hz(i,j) = Hz(i,j) + pulse(n)
!   endif
!   if(x(i) == x_detect .and. y(j) == y_detect)then
    P_sum = P_sum + (Hz(i,j) + Hz(i-1,j) + Hz(i,j-1) + Hz(i-1,j-1) )/4.0
!   endif
!  enddo
! enddo
endif

 
 do j=1,N_loc
!  PML for left Hz, x-direction
  do i=1,npml-1
   psi_Hzx_1(i,j)=bh_x(i)*psi_Hzx_1(i,j)+ch_x(i)*(Ey(i,j)-Ey(i+1,j))/dx
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzx_1(i,j)
  enddo
!  PML for right Hz, x-direction
  ii=npml-1
  do i=Nx+1-npml,Nx-1
   psi_Hzx_2(ii,j)=bh_x(ii)*psi_Hzx_2(ii,j)+ch_x(ii)*(Ey(i,j)-Ey(i+1,j))/dx
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzx_2(ii,j)
   ii=ii-1
  enddo
 enddo

endif ! 0 < myrank < nprocs-1

if(myrank==(nprocs-1))then !rank=(nprocs-1), here PML in y-direction is only applied to the top part
 do i=1,Nx-1
  do j=1,N_loc-1
   Hz(i,j)=Hz(i,j)+dt_mu0*((Ey(i,j)-Ey(i+1,j))*den_hx(i)+ &
			               (Ex(i,j+1)-Ex(i,j))*den_hy(j))
  enddo
 enddo
 
 do j=1,N_loc
!  PML for left Hz, x-direction
  do i=1,npml-1
   psi_Hzx_1(i,j)=bh_x(i)*psi_Hzx_1(i,j)+ch_x(i)*(Ey(i,j)-Ey(i+1,j))/dx
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzx_1(i,j)
  enddo
!  PML for right Hz, x-direction
  ii=npml-1
  do i=Nx+1-npml,Nx-1
   psi_Hzx_2(ii,j)=bh_x(ii)*psi_Hzx_2(ii,j)+ch_x(ii)*(Ey(i,j)-Ey(i+1,j))/dx
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzx_2(ii,j)
   ii=ii-1
  enddo
 enddo

!  PML for top Hz [top only here! since myrank=(nrpocs-1)], y-direction
 do i=1,Nx-1    
  jj=npml-1
  do j=N_loc+1-npml,N_loc-1
   psi_Hzy_2(i,jj)=bh_y(jj)*psi_Hzy_2(i,jj)+ch_y(jj)*(Ex(i,j+1)-Ex(i,j))/dy
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzy_2(i,jj)
   jj=jj-1
  enddo
 enddo

endif ! myrank == nprocs-1

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Ex ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
 do j=1,(nprocs-1) !Hz send/recieve for Ex simple updates
  itag=j
  itag1=j+1
  if(myrank==(j-1))then
    do i=1,Nx-1
     Hz_send(i)=Hz(i,N_loc)
    enddo
    call MPI_SEND(Hz_send,(Nx-1),MPI_doUBLE_PRECISION,j,itag,MPI_COMM_WORLD,ierr)
   elseif(myrank==j)then
    call MPI_RECV(Hz_get,(Nx-1),MPI_doUBLE_PRECISION,(j-1),itag,MPI_COMM_WORLD,istatus,ierr)
  endif!send_recv
 enddo!nprocs
 
if(myrank==0)then !rank=0, here PML in y-direction is only applied to the bottom part
 do i=1,Nx-1
  do j=2,N_loc
   Ex(i,j)=Ex(i,j)+dt_eps0*(Hz(i,j)-Hz(i,j-1))*den_ey(j)
  enddo
 enddo
   
!  PML for bottom Ex [bottom only here! since myrank=0], y-direction
 do i=1,Nx-1
  do j=2,npml
   psi_Exy_1(i,j)=be_y(j)*psi_Exy_1(i,j)+ce_y(j)*(Hz(i,j)-Hz(i,j-1))/dy
   Ex(i,j)=Ex(i,j)+dt_eps0*psi_Exy_1(i,j)
  enddo
 enddo

endif ! myrank == 0

if((myrank>0).and.(myrank<(nprocs-1)))then !no PML for y-direction here
 do i=1,Nx-1
  j=1
!   if(!FBx(i,j))then
!     tmpE=C1*Ex(i,j)+C3*(Hz(i,j)-Hz_get(i))*den_ey(j)-C4*PDx(i,j)
!     PDx(i,j)=A1*PDx(i,j)+A2*(tmpE+Ex(i,j))
!     Ex(i,j)=tmpE
!    else
     Ex(i,j)=Ex(i,j)+dt_eps0*(Hz(i,j)-Hz_get(i))*den_ey(j)
!   endif
  
  do j=2,N_loc
!   if(!FBx(i,j))then
!     tmpE=C1*Ex(i,j)+C3*(Hz(i,j)-Hz(i,j-1))*den_ey(j)-C4*PDx(i,j)
!     PDx(i,j)=A1*PDx(i,j)+A2*(tmpE+Ex(i,j))
!     Ex(i,j)=tmpE
!    else
     Ex(i,j)=Ex(i,j)+dt_eps0*(Hz(i,j)-Hz(i,j-1))*den_ey(j)
!   endif
  enddo
 enddo

endif ! 0 < myrank < nprocs-1

if(myrank==(nprocs-1))then !rank=(nprocs-1), here PML in y-direction is only applied to the top part
 j=1
  do i=1,Nx-1
   Ex(i,j)=Ex(i,j)+dt_eps0*(Hz(i,j)-Hz_get(i))*den_ey(j)
  enddo
   
 do i=1,Nx-1
  do j=2,(N_loc-1)
   Ex(i,j)=Ex(i,j)+dt_eps0*(Hz(i,j)-Hz(i,j-1))*den_ey(j)
  enddo
 enddo

!  PML for top Ex [top only here! since myrank=(nrpocs-1)], y-direction
 do i=1,Nx-1
  jj=npml
  do j=N_loc+1-npml,N_loc-1
   psi_Exy_2(i,jj)=be_y(jj)*psi_Exy_2(i,jj)+ce_y(jj)*(Hz(i,j)-Hz(i,(j-1)))/dy
   Ex(i,j)=Ex(i,j)+dt_eps0*psi_Exy_2(i,jj)
   jj=jj-1
  enddo
 enddo

endif ! myrank == nprocs-1

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Ey ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
if((myrank>=0).and.(myrank<(nprocs-1)))then
! do i=2,Nx-1
!  do j=1,N_loc
!   if(!FBy(i,j))then
!     tmpE=C1*Ey(i,j)+C3*(Hz(i-1,j)-Hz(i,j))*den_ex(i)-C4*PDy(i,j)
!     PDy(i,j)=A1*PDy(i,j)+A2*(tmpE+Ey(i,j))
!     Ey(i,j)=tmpE
!	else
!	 Ey(i,j)=Ey(i,j)+dt_eps0*(Hz(i-1,j)-Hz(i,j))*den_ex(i)
!   endif
!  enddo
! enddo

 do i=2,Nx-1
  do j=1,N_loc
   Ey(i,j)=Ey(i,j)+dt_eps0*((Hz(i-1,j)-Hz(i,j))*den_ex(i))
  enddo
 enddo
 
 do j=1,N_loc
!  PML for bottom Ey, x-direction
  do i=2,npml
   psi_Eyx_1(i,j)=be_x(i)*psi_Eyx_1(i,j)+ce_x(i)*(Hz(i-1,j)-Hz(i,j))/dx
   Ey(i,j)=Ey(i,j)+dt_eps0*psi_Eyx_1(i,j)
  enddo
!  PML for top Ey, x-direction
 ii=npml
 do i=Nx+1-npml,Nx-1
  psi_Eyx_2(ii,j)=be_x(ii)*psi_Eyx_2(ii,j)+ce_x(ii)*(Hz(i-1,j)-Hz(i,j))/dx
  Ey(i,j)=Ey(i,j)+dt_eps0*psi_Eyx_2(ii,j)
  ii=ii-1
  enddo
 enddo

endif ! 0 <= myrank < nprocs-1

if(myrank==(nprocs-1))then
 do i=2,Nx-1
  do j=1,N_loc-1
   Ey(i,j)=Ey(i,j)+dt_eps0*((Hz(i-1,j)-Hz(i,j))*den_ex(i))
  enddo
 enddo

 do j=1,N_loc-1
!  PML for bottom Ey, x-direction
  do i=2,npml
   psi_Eyx_1(i,j)=be_x(i)*psi_Eyx_1(i,j)+ce_x(i)*(Hz(i-1,j)-Hz(i,j))/dx
   Ey(i,j)=Ey(i,j)+dt_eps0*psi_Eyx_1(i,j)
  enddo
!  PML for top Ey, x-direction
  ii=npml
  do i=Nx+1-npml,Nx-1
   psi_Eyx_2(ii,j)=be_x(ii)*psi_Eyx_2(ii,j)+ce_x(ii)*(Hz(i-1,j)-Hz(i,j))/dx
   Ey(i,j)=Ey(i,j)+dt_eps0*psi_Eyx_2(ii,j)
   ii=ii-1
  enddo
 enddo
endif

!--------------------------------------------------------------------------!
!~~~~~~~~~~~~~~~~~~~~~~~~==========================~~~~~~~~~~~~~~~~~~~~~~~~!
!--------------------------------------------------------------------------!
!~~~~~~~~~~~~~~~~~~~~~~~~        Grid Return         ~~~~~~~~~~~~~~~~~~~~~~!
!--------------------------------------------------------------------------!
!~~~~~~~~~~~~~~~~~~~~~~~~==========================~~~~~~~~~~~~~~~~~~~~~~~~!
!--------------------------------------------------------------------------!

!if( b == 1 .and. (a == 1 .or. a == 3))then
!
!if(Nreturn > 0.and.GR)then
! do k = 1,Nreturn 
!  if(n == n_return(k))then
!  
!   nn = nn + 1
!   write(str_n,*) n
!   
!   if(a == 1)then
!    filename = str_Hz//trim(adjustl(str_n))//'a=1'//suffix
!   elseif(a == 3)then
!    filename = str_Hz//trim(adjustl(str_n))//'a=3'//suffix
!   endif
!   
!   !filename = prefix//filename
!   open(file=trim(adjustl(filename)),position = 'append',unit=nn)
!    do j = j_return1,j_return2
!     write(nn,*) Hz(i_return1:i_return2,j)
!    enddo
!   close(unit=nn)
!  
!   
!   write(str_n,*) n
!   
!   if(a == 1)then
!    filename = str_Ex//trim(adjustl(str_n))//'a=1'//suffix
!   elseif(a == 3)then
!    filename = str_Ex//trim(adjustl(str_n))//'a=3'//suffix
!   endif
!   
!   !filename = prefix//filename
!   open(file=trim(adjustl(filename)),position = 'append',unit=nn*3)
!    do j = j_return1,j_return2
!     write(nn*3,*) Ex(i_return1:i_return2,j)
!    enddo
!   close(unit=nn*3) 
!
!  endif
! enddo
!endif
!
!endif !GR

if(myrank == 0.or.myrank == (nprocs)/2.or.myrank == nprocs-1)then
 if( b == 1 .and. a == 1 )then
  if(n == 100)then
   nn = 30 + myrank
   write(str_n,*) myrank
   
   filename = str_Hz//trim(adjustl(str_n))//suffix
   open(file = trim(adjustl(filename)), unit = nn)
    write(nn,*) Hz
   close(unit = nn)
   
   filename = str_Ey//trim(adjustl(str_n))//suffix
   open(file = trim(adjustl(filename)), unit = nn*4)
    write(nn*4,*) Ey
   close(unit = nn*4)
   
  endif 
 endif
endif !GR

enddo !Nt
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  END TIME STEP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!.:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:.

if(myrank == 0)then
 write(*,*) "P_sum = ", P_sum
 WRITE(*,*)"done time-stepping"
endif
    
end function Vacuum_CPML

end !Main
