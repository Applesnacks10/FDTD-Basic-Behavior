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


integer, parameter :: Ny=201,N_loc=Ny-1 !N_loc must equal Ny-1 for 1 proc
double precision, parameter :: y0=-100E-9,yM=100E-9

integer, parameter :: Nx=11
double precision, parameter :: x0=-5E-9,xM=5E-9

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

!
!~~~ EM field components ~~~!
!
double precision Ex(Nx-1,N_loc),Ey(Nx,N_loc),Hz(Nx-1,N_loc)
double precision Ex_inc(N_loc),Hz_inc(N_loc)

!
!~~~ Field Input ~~~!
!
integer, parameter :: js = N_loc/2, is = -1
double precision aBH(4)
double precision, parameter :: tau=0.36d-15,E0=1.0,omega=ev_to_radsec*3.0
double precision Jx(Nt),Jy(Nt)

!
!~~~ Drude model for Ag ~~~!
!
double precision, parameter :: eps_r=8.926,omegaD=ev_to_radsec*11.585,GammaD=ev_to_radsec*0.203
double precision, parameter :: A1=(2.0-GammaD*dt)/(2.0+GammaD*dt),A2=eps0*omegaD*omegaD*dt/(2.0+GammaD*dt)
double precision, parameter :: C1=(eps_r*eps0/dt-0.5*A2)/(eps_r*eps0/dt+0.5*A2)
double precision, parameter :: C3=1.0/(eps_r*eps0/dt+0.5*A2)
double precision, parameter :: C4=0.5*(A1+1.0)/(eps_r*eps0/dt+0.5*A2)
double precision PDy(Nx,N_loc), PDx(Nx-1,N_loc)
double precision tmpE

!
!~~~ Scattering Presence ~~~!
!

logical :: FBx(Nx-1,N_loc), FBy(Nx,N_loc)
double precision, parameter :: z1 = -5*dy, z2 = 5*dy

!
!~~~ Loop Indices; time ~~~!
!
integer i,ii,j,jj,n,nn,k,a,b
double precision t

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
!~~~ Grid Return ~~~!
!
double precision, parameter:: y_return1 = y0, y_return2 = yM, x_return1 = x0, x_return2 = xM-dx
integer, parameter :: Nreturn = 5
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

!~~~ Scatterer Presence ~~~!

do i = 1,Nx-1
 do j = 1,N_loc
  if(y(j) >= z1.and.y(j) <= z2)then
   FBx(i,j) = .true.
  else
   FBx(i,j) = .false.
  endif
 enddo
enddo

do i = 1,Nx
 do j = 1,N_loc
  if(yM2(j) >= z1.and.yM2(j) <= z2)then
   FBy(i,j) = .true.
  else
   FBy(i,j) = .false.
  endif
 enddo
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

 n_return(1) = 100*2
 n_return(2) = 200*2
 n_return(3) = 300*2
 n_return(4) = 400*2
 n_return(5) = 500*2
 
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
jj=npml-1
do j=1,N_loc
 if(j<=(npml-1))then
  den_hy(j)=1.0/(kappah_y(j)*dy)
 elseif(j >= (N_loc-(npml-1)))then
  den_hy(j)=1.0/(kappah_y(jj)*dy)
  jj=jj-1
 endif
enddo

den_ey=1.0/dy
jj=npml
do j=1,N_loc
 if(j<=npml)then
  den_ey(j)=1.0/(kappae_y(j)*dy)
 elseif(j >= (N_loc-(npml-1)))then
  den_ey(j)=1.0/(kappae_y(jj)*dy)
  jj=jj-1
 endif
enddo

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


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
   !~~~ Initialize to Zero ~~~!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

Ex=0.0
Ey=0.0
Hz=0.0

Ex_inc=0.0
Hz_inc=0.0

Jx = 0.0
Jy = 0.0


!~~~ Source ~~~!

aBH(1)=0.353222222
aBH(2)=-0.488
aBH(3)=0.145
aBH(4)=-0.010222222

do n=1,Nt
 t=dt*dble(n)
 if(t<=tau)then
  Jx(n)= -dy/dt_eps0*E0*cos(omega*t)*( &
                  aBH(1)+ &
		  aBH(2)*cos(2.0*pi*t/tau)+ &
		  aBH(3)*cos(2.0*pi*2.0*t/tau)+ &
		  aBH(4)*cos(2.0*pi*3.0*t/tau))
  Jy(n) = Jx(n)*dx/dy
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
 
  Hz(i,j)=Hz(i,j)+dt_mu0*((Ey(i,j)-Ey(i+1,j))*den_hx(i) + &
			               (Ex(i,j+1)-Ex(i,j))*den_hy(j))
  
 enddo
enddo

			              
do j=1,N_loc-1
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

!  PML for top Hz [top only here! since myrank=(nrpocs-1)], y-direction
do i=1,Nx-1    
 jj=npml-1
 do j=N_loc+1-npml,N_loc-1
  psi_Hzy_2(i,jj)=bh_y(jj)*psi_Hzy_2(i,jj)+ch_y(jj)*(Ex(i,j+1)-Ex(i,j))/dy
  Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzy_2(i,jj)
  jj=jj-1
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
  
  if(FBx(i,j))then !Drude update
   tmpE=C1*Ex(i,j)+C3*(Hz(i,j)-Hz(i,j-1))/dy-C4*PDx(i,j)
   PDx(i,j)=A1*PDx(i,j)+A2*(tmpE+Ex(i,j))
   Ex(i,j)=tmpE
  if(j == js)then
   Ex(i,j) = Ex(i,j) - dt_eps0/dy*Jx(n) !add current source AFTER Drude polarization current update.
  endif
  else !Vacuum update
   Ex(i,j)=Ex(i,j)+dt_eps0*(Hz(i,j)-Hz(i,j-1))/dy
   if(j == js)then !add current source
    Ex(i,j) = Ex(i,j) - dt_eps0/dy*Jx(n)
   endif
  endif
  
 enddo
enddo

!  PML for bottom Ex [bottom only here! since myrank=0], y-direction
do i=1,Nx-1
 do j=2,npml
  psi_Exy_1(i,j)=be_y(j)*psi_Exy_1(i,j)+ce_y(j)*(Hz(i,j)-Hz(i,j-1))/dy
  Ex(i,j)=Ex(i,j)+dt_eps0*psi_Exy_1(i,j)
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


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Ey ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
 
do i=2,Nx-1
 do j=1,N_loc
  
  if(FBy(i,j))then !Drude update
   tmpE=C1*Ey(i,j)+C3*(Hz(i-1,j)-Hz(i,j))/dx-C4*PDy(i,j)
   PDy(i,j)=A1*PDy(i,j)+A2*(tmpE+Ey(i,j))
   Ey(i,j)=tmpE
   if(i == is.and.j/=N_loc)then
    Ey(i,j) = Ey(i,j) - dt_eps0/dx*Jy(n) !add current source AFTER Drude current update
   endif
  else !Vacuum update
   Ey(i,j)=Ey(i,j)+dt_eps0*(Hz(i-1,j)-Hz(i,j))/dx
   if(i == is.and.j/=N_loc)then !add current source
    Ey(i,j) = Ey(i,j) - dt_eps0/dx*Jy(n)
   endif
  endif
 
 enddo
enddo

do j=1,N_loc
!  PML for left Ey, x-direction
 do i=2,npml
  psi_Eyx_1(i,j)=be_x(i)*psi_Eyx_1(i,j)+ce_x(i)*(Hz(i-1,j)-Hz(i,j))/dx
  Ey(i,j)=Ey(i,j)+dt_eps0*psi_Eyx_1(i,j)
 enddo
!  PML for right Ey, x-direction
 ii=npml
 do i=Nx+1-npml,Nx-1
  psi_Eyx_2(ii,j)=be_x(ii)*psi_Eyx_2(ii,j)+ce_x(ii)*(Hz(i-1,j)-Hz(i,j))/dx
  Ey(i,j)=Ey(i,j)+dt_eps0*psi_Eyx_2(ii,j)
  ii=ii-1
 enddo
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
   filename = str_Ey//trim(adjustl(str_n))//suffix
   !filename = prefix//filename
   open(file=trim(adjustl(filename)),position = 'append',unit=a+10)
    do j = j_return1,j_return2
     write(a+10,*) Ey(i_return1:i_return2,j)
    enddo
   close(unit=a+10)
  
   
   write(str_n,*) n
   filename = str_Ex//trim(adjustl(str_n))//suffix
   !filename = prefix//filename
   open(file=trim(adjustl(filename)),position = 'append',unit=a+20)
    do j = j_return1,j_return2
     write(a+20,*) Ex(i_return1:i_return2,j)
    enddo
   close(unit=a+20)
  
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
