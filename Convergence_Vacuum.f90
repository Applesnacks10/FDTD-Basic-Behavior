program Convergence_fdtd3D_CPML
implicit none
include 'mpif.h'

integer, parameter :: Nr = 3
integer, parameter, dimension(Nr) :: res_array = (/1,2,4/)
integer, parameter, dimension(2) :: pml_add = (/0,1/)
double precision :: Convergence(Nr,2), Rel_error(Nr)
integer :: a,b !loop variables

!-------------------------------------------------------------------
!------------------------ Global Parameters ------------------------
!-------------------------------------------------------------------

   double precision, parameter :: length_add = 1.0E-2 
!  ..................................
!  Input Fundamental Constants (MKS units)
   double precision, PARAMETER ::                            &
      pi = 3.14159265358979, C = 2.99792458E8, &
      muO = 4.0 * pi * 1.0E-7, epsO = 1.0/(C*C*muO)

!  ..................................
!  Specify Material Relative Permittivity and Conductivity
   double precision, PARAMETER::                      &
      epsR = 1.0, sigM1 = 0.0   ! free space
      
!  ..................................
!  Specify the Impulsive Source (See Equation 7.134)
   double precision, PARAMETER ::                                        &
      tw = 53.0E-12, tO = 4.0*tw  

!  ..................................
!  Specify the CPML Order and Other Parameters
   INTEGER, PARAMETER ::                        & 
      m = 3, ma = 1 
   double precision, PARAMETER  ::                     &
!      sig_x_max = 0.75 * (0.8*(m+1)/(dx*(muO/epsO*epsR)**0.5)),   &
!      sig_y_max = 0.75 * (0.8*(m+1)/(dy*(muO/epsO*epsR)**0.5)),   &
!      sig_z_max = 0.75 * (0.8*(m+1)/(dz*(muO/epsO*epsR)**0.5)),   &
      sig_x_max = 6.370604950428188  ,&
      sig_y_max = sig_x_max  ,&
      sig_z_max = sig_x_max  ,&
      alpha_x_max = 0.24,   &
      alpha_y_max = alpha_x_max, alpha_z_max = alpha_x_max, &
      kappa_x_max = 15.0, &
      kappa_y_max = kappa_x_max, kappa_z_max = kappa_x_max
      
   INTEGER ::                                                &
	i,j,ii,jj,k,kk,n

   double precision  ::                                                   &
      source, P1, P2
   
   double precision :: DA, DB

!-------------------------------------------------------------------
!----------------- Minimum Resolution Variables --------------------
!-------------------------------------------------------------------

 double precision, parameter :: &
      dx_max = , dy_max = 
            
 double precision, parameter :: &
      dt_max =

   INTEGER, parameter :: &
      nmax_min =,                     & 
      Nx_min =  + 1,                  & 
      Ny_min =  + 1,                  &
      N_loc_min = ceiling(real((Ny_min - 1))/real(nprocs))
      
!  ..................................
!  Convergence Detection Zone
   INTEGER, parameter :: &
      i_start_min = (Nx_min-1)/2 - 2, i_end_min = (Nx_min-1)/2 + 2,                  &
      j_start_min = (N_loc_min-1)/2 - 4, j_end_min = (N_loc_min-1)/2 + 4,                  &

   INTEGER, parameter :: &
      istart_min = (Nx_min-1)/2-1,  iend_min = (Nx_min-1)/2,                    &
      jstart_min = (N_loc_min-1)/2-10, jend_min = (N_loc_min-1)/2,                    &                  &
      isource_min = istart_min, jsource_min = jstart_min
      
   INTEGER, parameter :: &
      npml_min = 10 + 1
      
!-------------------------------------------------------------------
!----------------- Maximum Resolution Variables --------------------
!-------------------------------------------------------------------
      
  double precision, parameter :: &
      dx_min = (dx_max)/res_array(Nr), dy_min =(dy_max)/res_array(Nr)
            
 double precision, parameter :: &
      dt_min = (dt_max)/res_array(Nr)

   INTEGER, parameter :: &
      nmax_max = res_array(Nr)*nmax_min,                  &                
      Nx_max = res_array(Nr)*(Nx_min-1) + length_add/dx_min + 1,                  &
      Ny_max = res_array(Nr)*(Ny_min-1) + length_add/dy_min + 1,                  &
      N_loc_max = ceiling(real((Ny_max - 1))/real(nprocs))
      
!  ..................................
!  Convergence Detection Zone
   INTEGER, parameter :: &
      i_start_max = (Nx_max-1)/2 - 2*res_array(Nr), i_end_max = (Nx_max-1)/2 + 2*res_array(Nr),                  &
      j_start_max = (N_loc_max-1)/2 - 4*res_array(Nr), j_end_max = (N_loc_max-1)/2 + 4*res_array(Nr),                  &

   INTEGER, parameter :: &
      istart_max = (Nx_max-1)/2-1*res_array(Nr),  iend_max = (Nx_max-1)/2,                    &
      jstart_max = (Ny_max-1)/2-10*res_array(Nr), jend_max = (Ny_max-1)/2,                    &
      isource_max = istart_max, jsource_max = jstart_max
      
   INTEGER, parameter :: &
      npml_max = 10*res_array(Nr) + length_add/dx_min + 1

!-------------------------------------------------------------------
!------------------------ Field and CPML Arrays --------------------
!-------------------------------------------------------------------

double precision :: Hz(Nx_max-1, N_loc_max)

double precision :: Ex(Nx_max-1, N_loc_max)

double precision :: Ey(Nx_max, N_loc_max)

!PML

double precision :: psi_Hzx_1(npml_max-1,N_loc_max)

double precision :: psi_Hzx_2(npml_max-1,N_loc_max)

double precision :: psi_Eyx_1(npml_max,N_loc_max)

double precision :: psi_Eyx_2(npml_max,N_loc_max)

double precision :: psi_Hzy_1(Nx_max-1,npml_max-1)

double precision :: psi_Hzy_2(Nx_max-1,npml_max-1)                          

double precision :: psi_Exy_1(Nx_max-1,npml_max)

double precision :: psi_Exy_2(Nx_max-1,npml_max)                        

double precision :: bh_x(npml_max-1), ch_x(npml_max-1), alphah_x(npml_max-1), sigh_x(npml_max-1), kappah_x(npml_max-1)

double precision :: be_x(npml_max), ce_x(npml_max), alphae_x(npml_max), sige_x(npml_max), kappae_x(npml_max)

double precision :: bh_y(npml_max-1), ch_y(npml_max-1), alphah_y(npml_max-1), sigh_y(npml_max-1), kappah_y(npml_max-1)

double precision :: be_y(npml_max), ce_y(npml_max), alphae_y(npml_max), sige_y(npml_max), kappae_y(npml_max)

!denominators for the update equations

double precision :: den_ex(Nx_max-1), den_hx(Nx_max-1)

double precision :: den_ey(N_loc_max), den_hy(N_loc_max)
      
!-----------------------------------------------------------------------
!--------------------------  Convergence Loop --------------------------
!-----------------------------------------------------------------------

 Convergence = 0.0

do a = 1,Nr
 do b = 1,2
 
  Convergence(a,b) = fdtd3D_CPML()
  
 enddo! 2 cpml lengths
 
 Rel_error(a) = abs((Convergence(a,2) - Convergence(a,1))/Convergence(a,1))
 
 open(file = 'Relative Errors.dat', position = 'append', unit = 40)
   write(40,*) res_array(a), Rel_error(a)
 close(unit = 40)
 
enddo! Nr resolutions

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!------------------------ Internal Function ----------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

 contains
 
function Vacuum_CPML() result(P_sum_fdtd)

   double precision :: P_sum_fdtd
   double precision :: Px, Py, Pz, P_sum

!-------------------------------------------------------------------
!----------------------- Variable Declaration ----------------------
!-------------------------------------------------------------------


!  ..................................
!  Specify Grid Cell Size in Each Direction and Calculate the 
!  Resulting Courant-Stable Time Step
   double precision ::                                        &
      dx, dy
!  ..................................
!  Specify Number of Time Steps and Grid Size Parameters
   INTEGER ::                                     &
      nmax, Nx, Ny, N_loc
      
!  ..................................
!  Convergence Detection Zone
   integer :: i_start, i_end, &
              j_start, j_end, &
     
      
   double precision ::             &
       dt

!  ..................................
!  Specify the PEC Plate Boundaries and the Source/Recording Points
   INTEGER ::                                    &
      istart, iend, jstart, jend,                       &
      isource, jsource
!  ..................................
!  Specify the CPML Thickness in Each Direction (Value of Zero 
!  Corresponds to No PML, and the Grid is Terminated with a PEC)
   INTEGER ::
      npml


!-------------------------------------------------------------------
!----------------------- Variable Assignment -----------------------
!-------------------------------------------------------------------

 dx = dx_max/res_array(a)
 dy = dy_max/res_array(a)
 
 nmax = res_array(a)*nmax_min
 Nx = res_array(a)*(Nx_min-1) + pml_add(b)*length_add/dx + 1
 Ny = res_array(a)*(Ny_min-1) + pml_add(b)*length_add/dy + 1 !Beware Truncated Length
 N_loc = res_array(a)*(N_loc_min) + ceiling(real(pml_add(b)*length_add)/real(dy)/real(nprocs))

 i_start = (Nx-1)/2 - 2*res_array(a)
 j_start = (N_loc)/2 - 4*res_array(a)
 i_end   = (Nx-1)/2 + 2*res_array(a)   
 j_end   = (N_loc)/2 + 4*res_array(a)

 dt = 1.906574869531006E-12/res_array(a)
 
 istart = (Nx-1)/2 - 1*res_array(a)
 iend   = (Nx-1)/2 + 1*res_array(a)
 jstart = (N_loc)/2 - 10*res_array(a)
 jend   = (N_loc)/2 + 10*res_array(a)
 
 isource = istart
 jsource = jstart
 
 npml = res_array(a)*(npml_min-1) + pml_add(b)*length_add/dx + 1
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  INITIALIZE VARIABLES
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   P1 = 0.0
   P2 = 0.0
   
   Hz(:,:,:) = 0.0
   Ex(:,:,:) = 0.0
   Ey(:,:,:) = 0.0
   sig(:,:,:) = sigM1
   eps(:,:,:) = epsR*epsO

   psi_Exy_1(:,:,:) = 0.0
   psi_Exy_2(:,:,:) = 0.0
   psi_Exz_1(:,:,:) = 0.0
   psi_Exz_2(:,:,:) = 0.0
   psi_Eyx_1(:,:,:) = 0.0
   psi_Eyx_2(:,:,:) = 0.0
   psi_Eyz_1(:,:,:) = 0.0
   psi_Eyz_2(:,:,:) = 0.0
   psi_Hzy_1(:,:,:) = 0.0
   psi_Hzy_2(:,:,:) = 0.0
   psi_Hzx_1(:,:,:) = 0.0
   psi_Hzx_2(:,:,:) = 0.0
   
   P_sum_fdtd = 0.0

   write(*,*)"res: ", res_array(a)
   write(*,*)"pml_add: ", pml_add(b)
  !defect check
   if(istart < 1.or. iend > Nx-1)then
    write(*,*)"Error -- i(end/start) reference not within range [1,Nx-1]"
   endif
   if(jstart < 1.or. jend > N_loc)then
    write(*,*)"Error -- j(end/start) reference not within range [1,N_loc]"
   endif
   
!
!~~~ fundamental constants [all numbers are in SI units]~~~!
!
double precision, parameter :: pi=3.1415926535897932384626433832795D0,c=299792458.0D0
double precision, parameter :: mu0=4.0D-7*pi,eps0=1.0D0/(c*c*mu0)
double precision, parameter :: h=1.054571628D-34
double complex, parameter :: Im=(0.0,1.0)
double precision, parameter :: ev_to_radsec=2.0*pi*2.4180e14
!
!~~~ number of grid points & time steps ~~~!
!
integer, parameter :: Nt=200000,N_w=400
double precision, parameter :: omega_min=ev_to_radsec*1.5,omega_max=ev_to_radsec*4.0

integer, parameter :: Ny=1281,N_loc=40
double precision, parameter :: y0=-640.0D-9,yM=640.0D-9

integer, parameter :: Nx=301
double precision, parameter :: x0=-150.0e-9,xM=150.0e-9
!
!~~~ scattered field zone ~~~!
!
integer, parameter :: i0=26,i1=276  !<--- +/- 125nm
integer, parameter :: mj0=1,j0=11   !<--- - 590nm
integer, parameter :: mj1=28,j1=21  !<--- + 500nm

integer, parameter :: ms=29,js=21

double precision, parameter :: tau=0.36d-15,E0=1.0,omega=ev_to_radsec*3.0
double precision aBH(4)
double precision pulse(Nt)
double precision tmp1,tmp2,omega_P(N_w),SN(N_w,2)
!
!~~~ detection ~~~!
!
double precision tmp

integer, parameter :: iw1=i0,iw2=i1
integer, parameter :: mwR=30,jwR=31  !<--- + 590nm
integer, parameter :: mwT=1,jwT=21   !<--- - 580nm

double precision Ex_temp(iw1:iw2,N_w,2),Hz_temp(iw1:iw2,N_w,2)
double precision Ex_temp_inc(iw1:iw2,N_w,2),Hz_temp_inc(iw1:iw2,N_w,2)
double precision Ex_w,Hz_w,av_x1,av_x2,av_y
double complex sum_int,cTMP1,cTMP2
double precision P_sct(N_w),P_inc(N_w)
!
!~~~ ph_zsical grid ~~~!
!
double precision, parameter :: dx=(xM-x0)/(Nx-1),dy=(yM-y0)/(Ny-1)
double precision x(Nx),xM2(Nx-1),y(N_loc),yM2(N_loc)
!
!~~~ dielectric constant of the host media (=1 for vacuum) ~~~!
!
double precision, parameter :: eps_delectric=1.0
!
!~~~ time step ~~~!
!
double precision, parameter :: dt=dy/(2.0*c)
!
!~~~ CPML ~~~!
!
integer, parameter :: npml=19,m=3,ma=1 
double precision sigmaCPML,alphaCPML,kappaCPML
double precision psi_Hzy_1(Nx-1,npml-1),psi_Exy_1(Nx-1,npml)                              
double precision psi_Hzy_2(Nx-1,npml-1),psi_Exy_2(Nx-1,npml)
double precision be_y(npml),ce_y(npml),alphae_y(npml),sige_y(npml),kappae_y(npml)
double precision bh_y(npml-1),ch_y(npml-1),alphah_y(npml-1),sigh_y(npml-1),kappah_y(npml-1)
double precision den_ex(Nx),den_hx(Nx),den_ey(N_loc),den_hy(N_loc)

double precision psi_Hzx_1(npml-1,N_loc),psi_Eyx_1(npml,N_loc)
double precision psi_Hzx_2(npml-1,N_loc),psi_Eyx_2(npml,N_loc)
double precision be_x(npml),ce_x(npml),alphae_x(npml),sige_x(npml),kappae_x(npml)
double precision bh_x(npml-1),ch_x(npml-1),alphah_x(npml-1),sigh_x(npml-1),kappah_x(npml-1)

double precision psi_Hzy_1_inc(npml-1),psi_Exy_1_inc(npml)                              
double precision psi_Hzy_2_inc(npml-1),psi_Exy_2_inc(npml)
!
!~~~ Drude model for Ag ~~~!
!
double precision, parameter :: eps_r=8.926,omegaD=ev_to_radsec*11.585,GammaD=ev_to_radsec*0.203
double precision, parameter :: A1=(2.0-GammaD*dt)/(2.0+GammaD*dt),A2=eps0*omegaD*omegaD*dt/(2.0+GammaD*dt)
double precision, parameter :: C1=(eps_r*eps0/dt-0.5*A2)/(eps_r*eps0/dt+0.5*A2)
double precision, parameter :: C3=1.0/(eps_r*eps0/dt+0.5*A2)
double precision, parameter :: C4=0.5*(A1+1.0)/(eps_r*eps0/dt+0.5*A2)

double precision tmpE

double precision PDx(Nx-1,N_loc),PDy(Nx,N_loc)

logical FBx(Nx-1,N_loc),FBy(Nx,N_loc)

double precision, parameter :: R=83.2525D-9
double precision, parameter :: z1=-75.2525D-9,z2=75.2525d-9
double precision, parameter :: slit_length=200.2525d-9 !should be < 2*x(i1)

!
!~~~ EM field components ~~~!
!
double precision Ex(Nx-1,N_loc),Ey(Nx,N_loc),Hz(Nx-1,N_loc)
double precision Ex_inc(N_loc),Hz_inc(N_loc)
!
!~~~ other parameters and variables ~~~!
!
integer i,ii,j,jj,n,nn,k
double precision t
double precision, parameter :: dt_eps0=dt/eps0,dt_mu0=dt/mu0
!
!~~~ MPI part ~~~!
!
integer ierr,nprocs,myrank,j_glob,mfdtd,n1
integer :: istatus(MPI_STATUS_SIZE)
integer itag,ireq,itag1,itag2,itag3,itag4,itag5,itag6
double precision Ex_get(Nx-1),Ex_send(Nx-1)
double precision Hz_get(Nx-1),Hz_send(Nx-1)
double precision Ex_send_inc,Ex_get_inc
double precision Hz_send_inc,Hz_get_inc

!-------------------------!
 call MPI_INIT(ierr)
 call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
 call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

do nn=1,N_w
 omega_P(nn)=omega_min+(omega_max-omega_min)*(nn-1)/(N_w-1)
enddo

aBH(1)=0.353222222
aBH(2)=-0.488
aBH(3)=0.145
aBH(4)=-0.010222222

pulse=0.0
SN=0.0
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

 do nn=1,N_w
  tmp1=sin(omega_P(nn)*t)
  tmp2=cos(omega_P(nn)*t)
  SN(nn,1)=SN(nn,1)+pulse(n)*tmp1
  SN(nn,2)=SN(nn,2)+pulse(n)*tmp2
 enddo
enddo

do nn=1,N_w
 tmp1=sqrt(SN(nn,1)**2+SN(nn,2)**2)
 SN(nn,2)=-atan2(SN(nn,1),SN(nn,2))
 SN(nn,1)=tmp1
enddo

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
enddo
do j=1,N_loc
 j_glob=myrank*N_loc+j
 yM2(j)=y0+dy*(j_glob-1)+dy/2.0
enddo

FBx=.false.
FBy=.false.

!~~~ structure ~~~!
do i=1,Nx-1
 do j=1,N_loc
  if( &
    ((y(j)>z1).and.(y(j)<z2).and.(xM2(i)<(-R)).and.(xM2(i)>(-slit_length/2.0))).or. &
    ((y(j)>z1).and.(y(j)<z2).and.(xM2(i)>R).and.(xM2(i)<(slit_length/2.0))) &
     )then
    FBx(i,j)=.true.
   else
    FBx(i,j)=.false.
  endif
 enddo
enddo

do i=1,Nx
 do j=1,N_loc
  if( &
    ((yM2(j)>z1).and.(yM2(j)<z2).and.(x(i)<(-R)).and.(x(i)>(-slit_length/2.0))).or. &
    ((yM2(j)>z1).and.(yM2(j)<z2).and.(x(i)>R).and.(x(i)<(slit_length/2.0))) &
     )then
    FBy(i,j)=.true.
   else
    FBy(i,j)=.false.
  endif
 enddo
enddo

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         !~~~ CPML vectors ~~~!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
sigmaCPML=0.8*(m+1)/(dx*(mu0/eps0*eps_delectric)**0.5)
alphaCPML=0.05
kappaCPML=5.0
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
 sigh_y(i)=sigmaCPML*((npml-i-0.5)/(npml-1.0))**m
 alphah_y(i)=alphaCPML*((i-0.5)/(npml-1.0))**ma
 kappah_y(i)=1.0+(kappaCPML-1.0)*((npml-i-0.5)/(npml-1.0))**m
 bh_y(i)=exp(-(sigh_y(i)/kappah_y(i)+alphah_y(i))*dt/eps0)
 ch_y(i)=sigh_y(i)*(bh_y(i)-1.0)/(sigh_y(i)+kappah_y(i)*alphah_y(i))/kappah_y(i)
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

den_h_z=1.0/dy
if(myrank==0)then
  do j=1,N_loc
   if(j<=(npml-1))then
    den_h_z(j)=1.0/(kappah_y(j)*dy)
   endif
  enddo
 elseif(myrank==(nprocs-1))then
  jj=npml-1
  do j=1,(N_loc-1)
   if(j>=(N_loc+1-npml))then
     den_h_z(j)=1.0/(kappaFLAG(jj)*dy)
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
   den_h_z(i)=1.0/(kappah_y(i)*dx)
  elseif(i>=(Nx+1-npml))then
   den_h_z(i)=1.0/(kappaFLAG(ii)*dx)
   ii=ii-1
  else
   den_h_z(i)=1.0/dx
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

PDx=0.0
PDy=0.0

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

Ex_inc=0.0
Hz_inc=0.0

Ex_get_inc=0.0
Ex_send_inc=0.0
Hz_get_inc=0.0
Hz_send_inc=0.0

psi_Hzy_1_inc=0.0
psi_Exy_1_inc=0.0
psi_Hzy_2_inc=0.0
psi_Exy_2_inc=0.0

Ex_temp=0.0
Hz_temp=0.0

Ex_temp_inc=0.0
Hz_temp_inc=0.0

if(myrank==mwT)then
 call cpu_time(cpu1)
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
    Ex_send_inc=Ex_inc(1)
    call MPI_SEND(Ex_send,(Nx-1),MPI_doUBLE_PRECISION,(j-1),itag,MPI_COMM_WORLD,ierr)
    call MPI_SEND(Ex_send_inc,1,MPI_doUBLE_PRECISION,(j-1),itag1,MPI_COMM_WORLD,ierr)
   elseif(myrank==(j-1))then
    call MPI_RECV(Ex_get,(Nx-1),MPI_doUBLE_PRECISION,j,itag,MPI_COMM_WORLD,istatus,ierr)
    call MPI_RECV(Ex_get_inc,1,MPI_doUBLE_PRECISION,j,itag1,MPI_COMM_WORLD,istatus,ierr)
  endif!send_recv
 enddo!nprocs

if(myrank==0)then !rank=0, here PML in y-direction is only applied to the bottom part
 do i=1,Nx-1
  do j=1,N_loc-1
   Hz(i,j)=Hz(i,j)+dt_mu0*((Ey(i,j)-Ey(i+1,j))*den_h_z(i)+ &
			               (Ex(i,j+1)-Ex(i,j))*den_h_z(j))
  enddo
    
  j=N_loc
   Hz(i,j)=Hz(i,j)+dt_mu0*((Ey(i,j)-Ey(i+1,j))*den_h_z(i)+ &
			               (Ex_get(i)-Ex(i,j))*den_h_z(j))
 enddo

 do j=1,N_loc
!  PML for left Hz, x-direction
  do i=1,npml-1
   psi_Hzx_1(i,j)=bh_y(i)*psi_Hzx_1(i,j)+ch_y(i)*(Ey(i,j)-Ey(i+1,j))/dx
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzx_1(i,j)
  enddo
!  PML for right Hz, x-direction
  ii=npml-1
  do i=Nx+1-npml,Nx-1
   psi_Hzx_2(ii,j)=bFLAG(ii)*psi_Hzx_2(ii,j)+cFLAG(ii)*(Ey(i,j)-Ey(i+1,j))/dx
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

!~~~~ incident ~~~~!
 do j=1,N_loc-1
  Hz_inc(j)=Hz_inc(j)+dt_mu0*(Ex_inc(j+1)-Ex_inc(j))*den_h_z(j)
 enddo
    
 j=N_loc
  Hz_inc(j)=Hz_inc(j)+dt_mu0*(Ex_get_inc-Ex_inc(j))*den_h_z(j)

!  PML for bottom Hz [bottom only here! since myrank=0], y-direction
 do j=1,npml-1
  psi_Hzy_1_inc(j)=bh_y(j)*psi_Hzy_1_inc(j)+ch_y(j)*(Ex_inc(j+1)-Ex_inc(j))/dy
  Hz_inc(j)=Hz_inc(j)+dt_mu0*psi_Hzy_1_inc(j)
 enddo
endif

if((myrank>0).and.(myrank<(nprocs-1)))then !no PML for y-direction here
 do i=1,Nx-1
  do j=1,N_loc-1
   Hz(i,j)=Hz(i,j)+dt_mu0*((Ey(i,j)-Ey(i+1,j))*den_h_z(i)+ &
			               (Ex(i,j+1)-Ex(i,j))*den_h_z(j))
  enddo
    
  j=N_loc
   Hz(i,j)=Hz(i,j)+dt_mu0*((Ey(i,j)-Ey(i+1,j))*den_h_z(i)+ &
			               (Ex_get(i)-Ex(i,j))*den_h_z(j))
 enddo

 do j=1,N_loc
!  PML for left Hz, x-direction
  do i=1,npml-1
   psi_Hzx_1(i,j)=bh_y(i)*psi_Hzx_1(i,j)+ch_y(i)*(Ey(i,j)-Ey(i+1,j))/dx
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzx_1(i,j)
  enddo
!  PML for right Hz, x-direction
  ii=npml-1
  do i=Nx+1-npml,Nx-1
   psi_Hzx_2(ii,j)=bFLAG(ii)*psi_Hzx_2(ii,j)+cFLAG(ii)*(Ey(i,j)-Ey(i+1,j))/dx
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzx_2(ii,j)
   ii=ii-1
  enddo
 enddo

!~~~~ incident ~~~~!
 do j=1,N_loc-1
  Hz_inc(j)=Hz_inc(j)+dt_mu0*(Ex_inc(j+1)-Ex_inc(j))*den_h_z(j)
 enddo
    
 j=N_loc
  Hz_inc(j)=Hz_inc(j)+dt_mu0*(Ex_get_inc-Ex_inc(j))*den_h_z(j)

! scattered/total field updates
 if(myrank==mj0)then
  do i=i0,i1-1
   Hz(i,j0-1)=Hz(i,j0-1)-dt_mu0*Ex_inc(j0)/dy
  enddo
 endif
 
 if(myrank==mj1)then
  do i=i0,i1-1
   Hz(i,j1)=Hz(i,j1)+dt_mu0*Ex_inc(j1)/dy
  enddo
 endif
endif

if(myrank==(nprocs-1))then !rank=(nprocs-1), here PML in y-direction is only applied to the top part
 do i=1,Nx-1
  do j=1,N_loc-1
   Hz(i,j)=Hz(i,j)+dt_mu0*((Ey(i,j)-Ey(i+1,j))*den_h_z(i)+ &
			               (Ex(i,j+1)-Ex(i,j))*den_h_z(j))
  enddo
 enddo
 
 do j=1,N_loc
!  PML for left Hz, x-direction
  do i=1,npml-1
   psi_Hzx_1(i,j)=bh_y(i)*psi_Hzx_1(i,j)+ch_y(i)*(Ey(i,j)-Ey(i+1,j))/dx
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzx_1(i,j)
  enddo
!  PML for right Hz, x-direction
  ii=npml-1
  do i=Nx+1-npml,Nx-1
   psi_Hzx_2(ii,j)=bFLAG(ii)*psi_Hzx_2(ii,j)+cFLAG(ii)*(Ey(i,j)-Ey(i+1,j))/dx
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzx_2(ii,j)
   ii=ii-1
  enddo
 enddo

!  PML for top Hz [top only here! since myrank=(nrpocs-1)], y-direction
 do i=1,Nx-1    
  jj=npml-1
  do j=N_loc+1-npml,N_loc-1
   psi_Hzy_2(i,jj)=bFLAG(jj)*psi_Hzy_2(i,jj)+cFLAG(jj)*(Ex(i,j+1)-Ex(i,j))/dy
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzy_2(i,jj)
   jj=jj-1
  enddo
 enddo

!~~~~ incident ~~~~!
 do j=1,N_loc-1
  Hz_inc(j)=Hz_inc(j)+dt_mu0*(Ex_inc(j+1)-Ex_inc(j))*den_h_z(j)
 enddo

!  PML for top Hz [top only here! since myrank=(nrpocs-1)], y-direction
 jj=npml-1
 do j=N_loc+1-npml,N_loc-1
  psi_Hzy_2_inc(jj)=bFLAG(jj)*psi_Hzy_2_inc(jj)+cFLAG(jj)*(Ex_inc(j+1)-Ex_inc(j))/dy
  Hz_inc(j)=Hz_inc(j)+dt_mu0*psi_Hzy_2_inc(jj)
  jj=jj-1
 enddo
endif

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
    Hz_send_inc=Hz_inc(N_loc)
    call MPI_SEND(Hz_send,(Nx-1),MPI_doUBLE_PRECISION,j,itag,MPI_COMM_WORLD,ierr)
    call MPI_SEND(Hz_send_inc,1,MPI_doUBLE_PRECISION,j,itag1,MPI_COMM_WORLD,ierr)
   elseif(myrank==j)then
    call MPI_RECV(Hz_get,(Nx-1),MPI_doUBLE_PRECISION,(j-1),itag,MPI_COMM_WORLD,istatus,ierr)
    call MPI_RECV(Hz_get_inc,1,MPI_doUBLE_PRECISION,(j-1),itag1,MPI_COMM_WORLD,istatus,ierr)
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

!~~~ incident ~~~!
 do j=2,N_loc
  Ex_inc(j)=Ex_inc(j)+dt_eps0*(Hz_inc(j)-Hz_inc(j-1))*den_ey(j)
 enddo
   
!  PML for bottom Ex [bottom only here! since myrank=0], y-direction
 do j=2,npml
  psi_Exy_1_inc(j)=be_y(j)*psi_Exy_1_inc(j)+ce_y(j)*(Hz_inc(j)-Hz_inc(j-1))/dy
  Ex_inc(j)=Ex_inc(j)+dt_eps0*psi_Exy_1_inc(j)
 enddo
endif

if((myrank>0).and.(myrank<(nprocs-1)))then !no PML for y-direction here
 do i=1,Nx-1
  j=1
   if(FBx(i,j))then
     tmpE=C1*Ex(i,j)+C3*(Hz(i,j)-Hz_get(i))*den_ey(j)-C4*PDx(i,j)
     PDx(i,j)=A1*PDx(i,j)+A2*(tmpE+Ex(i,j))
     Ex(i,j)=tmpE
    else
     Ex(i,j)=Ex(i,j)+dt_eps0*(Hz(i,j)-Hz_get(i))*den_ey(j)
   endif
  
  do j=2,N_loc
   if(FBx(i,j))then
     tmpE=C1*Ex(i,j)+C3*(Hz(i,j)-Hz(i,j-1))*den_ey(j)-C4*PDx(i,j)
     PDx(i,j)=A1*PDx(i,j)+A2*(tmpE+Ex(i,j))
     Ex(i,j)=tmpE
    else
     Ex(i,j)=Ex(i,j)+dt_eps0*(Hz(i,j)-Hz(i,j-1))*den_ey(j)
   endif
  enddo
 enddo

!~~~ incident ~~~!
 j=1
  Ex_inc(j)=Ex_inc(j)+dt_eps0*(Hz_inc(j)-Hz_get_inc)*den_ey(j)
   
 do j=2,N_loc
  if((myrank==ms).and.(j==js))then !laser pulse
    Ex_inc(j)=Ex_inc(j)+dt_eps0*(Hz_inc(j)-Hz_inc(j-1))/dy+ &
 	                    pulse(n)
   else
	Ex_inc(j)=Ex_inc(j)+dt_eps0*(Hz_inc(j)-Hz_inc(j-1))*den_ey(j)
  endif
 enddo

! scattered/total field updates
 if(myrank==mj0)then
  do i=i0,i1-1
   Ex(i,j0)=Ex(i,j0)-dt_eps0*Hz_inc(j0-1)/dy
  enddo
 endif

 if(myrank==mj1)then
  do i=i0,i1-1
   Ex(i,j1)=Ex(i,j1)+dt_eps0*Hz_inc(j1)/dy
  enddo
 endif
endif

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

!~~~ incident ~~~!
 j=1
  Ex_inc(j)=Ex_inc(j)+dt_eps0*(Hz_inc(j)-Hz_get_inc)*den_ey(j)
   
 do j=2,N_loc-1
  Ex_inc(j)=Ex_inc(j)+dt_eps0*(Hz_inc(j)-Hz_inc(j-1))*den_ey(j)
 enddo

!  PML for top Ex [top only here! since myrank=(nrpocs-1)], y-direction
 jj=npml
 do j=N_loc+1-npml,N_loc-1
  psi_Exy_2_inc(jj)=be_y(jj)*psi_Exy_2_inc(jj)+ce_y(jj)*(Hz_inc(j)-Hz_inc(j-1))/dy
  Ex_inc(j)=Ex_inc(j)+dt_eps0*psi_Exy_2_inc(jj)
  jj=jj-1
 enddo
endif

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Ey ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
if((myrank>=0).and.(myrank<(nprocs-1)))then
 do i=2,Nx-1
  do j=1,N_loc
   if(FBy(i,j))then
     tmpE=C1*Ey(i,j)+C3*(Hz(i-1,j)-Hz(i,j))*den_ex(i)-C4*PDy(i,j)
     PDy(i,j)=A1*PDy(i,j)+A2*(tmpE+Ey(i,j))
     Ey(i,j)=tmpE
	else
	 Ey(i,j)=Ey(i,j)+dt_eps0*(Hz(i-1,j)-Hz(i,j))*den_ex(i)
   endif
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

! scattered/total field updates
 if(myrank==mj0)then
  do j=j0,N_loc
   Ey(i0,j)=Ey(i0,j)+dt_eps0*Hz_inc(j)/dy
   Ey(i1,j)=Ey(i1,j)-dt_eps0*Hz_inc(j)/dy
  enddo
 endif

 if((myrank>mj0).and.(myrank<mj1))then
  do j=1,N_loc
   Ey(i0,j)=Ey(i0,j)+dt_eps0*Hz_inc(j)/dy
   Ey(i1,j)=Ey(i1,j)-dt_eps0*Hz_inc(j)/dy
  enddo
 endif

 if(myrank==mj1)then
  do j=1,j1-1
   Ey(i0,j)=Ey(i0,j)+dt_eps0*Hz_inc(j)/dy
   Ey(i1,j)=Ey(i1,j)-dt_eps0*Hz_inc(j)/dy
  enddo
 endif
endif

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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  END TIME STEP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!.:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:.
    WRITE(*,*)"done time-stepping"
    
end function fdtd3D_CPML

end !Main
