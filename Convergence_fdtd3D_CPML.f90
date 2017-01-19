program Convergence_fdtd3D_CPML
implicit none

integer, parameter :: Nr = 2
integer, parameter, dimension(Nr) :: res_array = (/1,2/)
integer, parameter, dimension(2) :: pml_add = (/0,1/)
double precision :: Convergence(Nr,2), Rel_error(Nr)
integer :: a,b !loop variables
double precision :: Px, Py, Pz, P_sum
 Convergence = 0.0

do a = 1,Nr
 do b = 1,2
 
  Convergence(a,b) = fdtd3D_CPML()
  
 enddo! 2 cpml lengths
 
 Rel_error(a) = abs((Convergence(a,2) - Convergence(a,1))/Convergence(a,1))
 
 open(file = 'Relative Errors.dat', unit = 40)
   write(40,*) res_array(a), Rel_error(a)
 close(unit = 40)
 
enddo! Nr resolutions

!Transmit Relative Errors

! open(file = 'Relative Errors', unit = 40)
!  do a = 1,Nr 
!   write(40,*) res_array(a), Rel_error(a)
!  enddo
! close(unit = 40)

!-----------------------------------------------------------------------
!------------------------ Internal Function ----------------------------
!-----------------------------------------------------------------------
 contains
 
function fdtd3D_CPML() result(P_sum_fdtd)

   double precision :: P_sum_fdtd

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

   REAL  ::                                                   &
      source, P1, P2

!-------------------------------------------------------------------
!---------------------- Start of Un-Parametrized Variables ---------
!-------------------------------------------------------------------


!  ..................................
!  Specify Grid Cell Size in Each Direction and Calculate the 
!  Resulting Courant-Stable Time Step
   double precision ::                                        &
      dx, dy, dz

!  ..................................
!  Specify Number of Time Steps and Grid Size Parameters
   INTEGER ::                                     &
      nmax, Imax, Jmax, Kmax
      
!  ..................................
!  Convergence Detection Zone
   integer :: i_start, i_end, &
              j_start, j_end, &
              k_start, k_end
     
      
   double precision ::             &
       dt

!  ..................................
!  Specify the PEC Plate Boundaries and the Source/Recording Points
   INTEGER ::                                    &
      istart, iend, jstart,   &
      jend, kstart, kend,    &
      isource, jsource, ksource

!  ..................................
!  Specify the CPML Thickness in Each Direction (Value of Zero 
!  Corresponds to No PML, and the Grid is Terminated with a PEC)
   INTEGER ::                         &
      nxPML_1, nxPML_2, nyPML_1,      &
      nyPML_2, nzPML_1, nzPML_2
      
   REAL :: DA, DB

!-------------------------------------------------------------------
!---------------------- Start of Allocatable Arrays ----------------
!-------------------------------------------------------------------


!     TM components
   REAL, allocatable, DIMENSION(:,:,:)  ::                      &
      Ez, CA, CB, sig, eps

   REAL, allocatable, DIMENSION(:,:,:)  ::                      &
      Hy

   REAL, allocatable, DIMENSION(:,:,:)  ::                      &
      Hx

!     TE components
   REAL, allocatable, DIMENSION(:,:,:)  ::                      &
      Hz

   REAL, allocatable, DIMENSION(:,:,:)  ::                      &
      Ex

   REAL, allocatable, DIMENSION(:,:,:)  ::                      &
      Ey

!  PML
   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Ezx_1

   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Ezx_2

   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Hyx_1

   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Hyx_2

   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Ezy_1                               

   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Ezy_2

   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Hxy_1                               

   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Hxy_2

   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Hxz_1

   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Hxz_2

   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Hyz_1

   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Hyz_2

   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Exz_1

   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Exz_2

   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Eyz_1

   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Eyz_2

   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Hzx_1
   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Eyx_1

   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Hzx_2
   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Eyx_2

   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Hzy_1                               
   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Exy_1                               

   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Hzy_2
   REAL , allocatable, DIMENSION(:,:,:) ::                       &
      psi_Exy_2

   REAL , allocatable, DIMENSION(:) ::                       &
      be_x_1, ce_x_1, alphae_x_PML_1, sige_x_PML_1, kappae_x_PML_1
   REAL , allocatable, DIMENSION(:) ::                       &
      bh_x_1, ch_x_1, alphah_x_PML_1, sigh_x_PML_1, kappah_x_PML_1

   REAL , allocatable, DIMENSION(:) ::                       &
      be_x_2, ce_x_2, alphae_x_PML_2, sige_x_PML_2, kappae_x_PML_2
   REAL , allocatable, DIMENSION(:) ::                       &
      bh_x_2, ch_x_2, alphah_x_PML_2, sigh_x_PML_2, kappah_x_PML_2

   REAL , allocatable, DIMENSION(:) ::                       &
      be_y_1, ce_y_1, alphae_y_PML_1, sige_y_PML_1, kappae_y_PML_1
   REAL , allocatable, DIMENSION(:) ::                       &
      bh_y_1, ch_y_1, alphah_y_PML_1, sigh_y_PML_1, kappah_y_PML_1

   REAL , allocatable, DIMENSION(:) ::                       &
      be_y_2, ce_y_2, alphae_y_PML_2, sige_y_PML_2, kappae_y_PML_2
   REAL , allocatable, DIMENSION(:) ::                       &
      bh_y_2, ch_y_2, alphah_y_PML_2, sigh_y_PML_2, kappah_y_PML_2

   REAL , allocatable, DIMENSION(:) ::                       &
      be_z_1, ce_z_1, alphae_z_PML_1, sige_z_PML_1, kappae_z_PML_1
   REAL , allocatable, DIMENSION(:) ::                       &
      bh_z_1, ch_z_1, alphah_z_PML_1, sigh_z_PML_1, kappah_z_PML_1

   REAL , allocatable, DIMENSION(:) ::                       &
      be_z_2, ce_z_2, alphae_z_PML_2, sige_z_PML_2, kappae_z_PML_2
   REAL , allocatable, DIMENSION(:) ::                       &
      bh_z_2, ch_z_2, alphah_z_PML_2, sigh_z_PML_2, kappah_z_PML_2

!     denominators for the update equations
   REAL, allocatable, DIMENSION(:)  ::                      &
      den_ex, den_hx

   REAL, allocatable, DIMENSION(:)  ::                      &
      den_ey, den_hy

   REAL, allocatable, DIMENSION(:)  ::                      &
      den_ez, den_hz
      
      
      
!-------------------------------------------------------------------
!---------------------- Start of Variable Assignment ---------------
!-------------------------------------------------------------------

 dx = 1.0D-3/res_array(a)
 dy = 1.0D-3/res_array(a)
 dz = 1.0D-3/res_array(a)
 
 nmax = res_array(a)*2100
 Imax = res_array(a)*51+pml_add(b)*length_add/dx
 Jmax = res_array(a)*127+pml_add(b)*length_add/dy
 Kmax = res_array(a)*27+pml_add(b)*length_add/dy

 i_start = ((Imax-1)/2 - 2)*res_array(a)
 j_start = ((Jmax-1)/2 - 4)*res_array(a)
 k_start = ((Kmax-1)/2 - 1)*res_array(a)
 i_end = ((Imax-1)/2 + 2)*res_array(a)    
 j_end = ((Jmax-1)/2 + 4)*res_array(a)
 k_end = ((Kmax-1)/2 + 1)*res_array(a)
     
 dt = 1.906574869531006E-12/res_array(a)
 
 istart = ((Imax-1)/2-11)*res_array(a)
 iend = istart+24*res_array(a)
 jstart = Jmax/2-49*res_array(a)
 jend = jstart + 99*res_array(a)
 kstart = Kmax/2
 kend = kstart
 
 isource = istart
 jsource = jstart
 ksource = kstart
 
 nxPML_1 = res_array(a)*11+pml_add(b)*length_add/dx
 nxPML_2 = nxPML_1
 nyPML_1 = nxPML_1
 nyPML_2 = nxPML_1
 nzPML_1 = nxPML_1
 nzPML_2 = nxPML_2
 
!-------------------------------------------------------------------
!---------------------- Start of Array Allocation ------------------
!-------------------------------------------------------------------
allocate(Ez(Imax, Jmax, Kmax-1), CA(Imax, Jmax, Kmax-1), CB(Imax, Jmax, Kmax-1), sig(Imax, Jmax, Kmax-1), eps(Imax, Jmax, Kmax-1))

allocate(Hy(Imax-1, Jmax, Kmax-1))

allocate(Hx(Imax,Jmax-1, Kmax-1))

allocate(Hz(Imax-1, Jmax-1, Kmax))

allocate(Ex(Imax-1, Jmax, Kmax))

allocate(Ey(Imax,Jmax-1, Kmax))

!PML

allocate(psi_Ezx_1(nxPML_1,Jmax,Kmax))

allocate(psi_Ezx_2(nxPML_2,Jmax,Kmax))

allocate(psi_Hyx_1(nxPML_1-1,Jmax,Kmax))

allocate(psi_Hyx_2(nxPML_2-1,Jmax,Kmax))

allocate(psi_Ezy_1(Imax,nyPML_1,Kmax))                               

allocate(psi_Ezy_2(Imax,nyPML_2,Kmax))

allocate(psi_Hxy_1(Imax,nyPML_1-1,Kmax))                               

allocate(psi_Hxy_2(Imax,nyPML_2-1,Kmax))

allocate(psi_Hxz_1(Imax,Jmax-1,nzPML_1-1))

allocate(psi_Hxz_2(Imax,Jmax-1,nzPML_2-1))

allocate(psi_Hyz_1(Imax-1,Jmax,nzPML_1-1))

allocate(psi_Hyz_2(Imax-1,Jmax,nzPML_2-1))

allocate(psi_Exz_1(Imax-1,Jmax,nzPML_1))

allocate(psi_Exz_2(Imax-1,Jmax,nzPML_2))

allocate(psi_Eyz_1(Imax,Jmax-1,nzPML_1))

allocate(psi_Eyz_2(Imax,Jmax-1,nzPML_2))

allocate(psi_Hzx_1(nxPML_1-1,Jmax-1,Kmax))

allocate(psi_Eyx_1(nxPML_1,Jmax-1,Kmax))

allocate(psi_Hzx_2(nxPML_2-1,Jmax-1,Kmax))

allocate(psi_Eyx_2(nxPML_2,Jmax-1,Kmax))

allocate(psi_Hzy_1(Imax-1,nyPML_1-1,Kmax))                              

allocate(psi_Exy_1(Imax-1,nyPML_1,Kmax))                               

allocate(psi_Hzy_2(Imax-1,nyPML_2-1,Kmax))

allocate(psi_Exy_2(Imax-1,nyPML_2,Kmax))

allocate(be_x_1(nxPML_1), ce_x_1(nxPML_1), alphae_x_PML_1(nxPML_1), sige_x_PML_1(nxPML_1), kappae_x_PML_1(nxPML_1))

allocate(bh_x_1(nxPML_1-1), ch_x_1(nxPML_1-1), alphah_x_PML_1(nxPML_1-1), sigh_x_PML_1(nxPML_1-1), kappah_x_PML_1(nxPML_1-1))

allocate(be_x_2(nxPML_2), ce_x_2(nxPML_2), alphae_x_PML_2(nxPML_2), sige_x_PML_2(nxPML_2), kappae_x_PML_2(nxPML_2))

allocate(bh_x_2(nxPML_2-1), ch_x_2(nxPML_2-1), alphah_x_PML_2(nxPML_2-1), sigh_x_PML_2(nxPML_2-1), kappah_x_PML_2(nxPML_2-1))

allocate(be_y_1(nyPML_1), ce_y_1(nyPML_1), alphae_y_PML_1(nyPML_1), sige_y_PML_1(nyPML_1), kappae_y_PML_1(nyPML_1))

allocate(bh_y_1(nyPML_1-1), ch_y_1(nyPML_1-1), alphah_y_PML_1(nyPML_1-1), sigh_y_PML_1(nyPML_1-1), kappah_y_PML_1(nyPML_1-1))

allocate(be_y_2(nyPML_2), ce_y_2(nyPML_2), alphae_y_PML_2(nyPML_2), sige_y_PML_2(nyPML_2), kappae_y_PML_2(nyPML_2))

allocate(bh_y_2(nyPML_2-1), ch_y_2(nyPML_2-1), alphah_y_PML_2(nyPML_2-1), sigh_y_PML_2(nyPML_2-1), kappah_y_PML_2(nyPML_2-1))

allocate(be_z_1(nzPML_1), ce_z_1(nzPML_1), alphae_z_PML_1(nzPML_1), sige_z_PML_1(nzPML_1), kappae_z_PML_1(nzPML_1))

allocate(bh_z_1(nzPML_1-1), ch_z_1(nzPML_1-1), alphah_z_PML_1(nzPML_1-1), sigh_z_PML_1(nzPML_1-1), kappah_z_PML_1(nzPML_1-1))

allocate(be_z_2(nzPML_2), ce_z_2(nzPML_2), alphae_z_PML_2(nzPML_2), sige_z_PML_2(nzPML_2), kappae_z_PML_2(nzPML_2))

allocate(bh_z_2(nzPML_2-1), ch_z_2(nzPML_2-1), alphah_z_PML_2(nzPML_2-1), sigh_z_PML_2(nzPML_2-1), kappah_z_PML_2(nzPML_2-1))

!denominators for the update equations

allocate(den_ex(Imax-1), den_hx(Imax-1))

allocate(den_ey(Jmax-1), den_hy(Jmax-1))

allocate(den_ez(Kmax-1), den_hz(Kmax-1))

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  INITIALIZE VARIABLES
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   P1 = 0.0
   P2 = 0.0

   Ez(:,:,:) = 0.0
   Hz(:,:,:) = 0.0
   Ex(:,:,:) = 0.0
   Hx(:,:,:) = 0.0
   Ey(:,:,:) = 0.0
   Hy(:,:,:) = 0.0
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
   psi_Ezy_1(:,:,:) = 0.0
   psi_Ezy_2(:,:,:) = 0.0
   psi_Ezx_1(:,:,:) = 0.0
   psi_Ezx_2(:,:,:) = 0.0
   psi_Hxy_1(:,:,:) = 0.0
   psi_Hxy_2(:,:,:) = 0.0
   psi_Hxz_1(:,:,:) = 0.0
   psi_Hxz_2(:,:,:) = 0.0
   psi_Hyx_1(:,:,:) = 0.0
   psi_Hyx_2(:,:,:) = 0.0
   psi_Hyz_1(:,:,:) = 0.0
   psi_Hyz_2(:,:,:) = 0.0
   psi_Hzy_1(:,:,:) = 0.0
   psi_Hzy_2(:,:,:) = 0.0
   psi_Hzx_1(:,:,:) = 0.0
   psi_Hzx_2(:,:,:) = 0.0
   
   P_sum_fdtd = 0.0

   write(*,*)"res: ", res_array(a)
   write(*,*)"pml_add: ", pml_add(b)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  SET CPML PARAMETERS IN EACH DIRECTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   DO i = 1,nxPML_1
      sige_x_PML_1(i) = sig_x_max * ( (nxPML_1 - i) / (nxPML_1 - 1.0) )**m
      alphae_x_PML_1(i) = alpha_x_max*((i-1.0)/(nxPML_1-1.0))**ma
      kappae_x_PML_1(i) = 1.0+(kappa_x_max-1.0)*   &
                                 ((nxPML_1 - i) / (nxPML_1 - 1.0))**m
      be_x_1(i) = EXP(-(sige_x_PML_1(i) / kappae_x_PML_1(i) +   &
                                 alphae_x_PML_1(i))*dt/epsO)
      if ((sige_x_PML_1(i) == 0.0) .and.        &
         (alphae_x_PML_1(i) == 0.0) .and. (i == nxPML_1)) then
         ce_x_1(i) = 0.0
      else
         ce_x_1(i) = sige_x_PML_1(i)*(be_x_1(i)-1.0)/       &
            (sige_x_PML_1(i)+kappae_x_PML_1(i)*alphae_x_PML_1(i)) &
            / kappae_x_PML_1(i)
      endif
   ENDDO
   DO i = 1,nxPML_1-1
      sigh_x_PML_1(i) = sig_x_max * ( (nxPML_1 - i - 0.5)/(nxPML_1-1.0))**m
      alphah_x_PML_1(i) = alpha_x_max*((i-0.5)/(nxPML_1-1.0))**ma
      kappah_x_PML_1(i) = 1.0+(kappa_x_max-1.0)*   &
                            ((nxPML_1 - i - 0.5) / (nxPML_1 - 1.0))**m
      bh_x_1(i) = EXP(-(sigh_x_PML_1(i) / kappah_x_PML_1(i) +   &
                                 alphah_x_PML_1(i))*dt/epsO)
      ch_x_1(i) = sigh_x_PML_1(i)*(bh_x_1(i)-1.0)/      &
                  (sigh_x_PML_1(i)+kappah_x_PML_1(i)*alphah_x_PML_1(i)) &
                  / kappah_x_PML_1(i)
   ENDDO

   DO i = 1,nxPML_2
      sige_x_PML_2(i) = sig_x_max * ( (nxPML_2 - i) / (nxPML_2 - 1.0) )**m
      alphae_x_PML_2(i) = alpha_x_max*((i-1.0)/(nxPML_2-1.0))**ma
      kappae_x_PML_2(i) = 1.0+(kappa_x_max-1.0)*   &
                                 ((nxPML_2 - i) / (nxPML_2 - 1.0))**m
      be_x_2(i) = EXP(-(sige_x_PML_2(i) / kappae_x_PML_2(i) +   &
                                 alphae_x_PML_2(i))*dt/epsO)
      if ((sige_x_PML_2(i) == 0.0) .and.        &
         (alphae_x_PML_2(i) == 0.0) .and. (i == nxPML_2)) then
         ce_x_2(i) = 0.0
      else
         ce_x_2(i) = sige_x_PML_2(i)*(be_x_2(i)-1.0)/       &
            (sige_x_PML_2(i)+kappae_x_PML_2(i)*alphae_x_PML_2(i)) &
            / kappae_x_PML_2(i)
      endif
   ENDDO
   DO i = 1,nxPML_2-1
      sigh_x_PML_2(i) = sig_x_max * ( (nxPML_2 - i - 0.5)/(nxPML_2-1.0))**m
      alphah_x_PML_2(i) = alpha_x_max*((i-0.5)/(nxPML_2-1.0))**ma
      kappah_x_PML_2(i) = 1.0+(kappa_x_max-1.0)*   &
                            ((nxPML_2 - i - 0.5) / (nxPML_2 - 1.0))**m
      bh_x_2(i) = EXP(-(sigh_x_PML_2(i) / kappah_x_PML_2(i) +   &
                                 alphah_x_PML_2(i))*dt/epsO)
      ch_x_2(i) = sigh_x_PML_2(i)*(bh_x_2(i)-1.0)/      &
                  (sigh_x_PML_2(i)+kappah_x_PML_2(i)*alphah_x_PML_2(i)) &
                  / kappah_x_PML_2(i)
   ENDDO

   DO j = 1,nyPML_1
      sige_y_PML_1(j) = sig_y_max * ( (nyPML_1 - j ) / (nyPML_1 - 1.0) )**m
      alphae_y_PML_1(j) = alpha_y_max*((j-1)/(nyPML_1-1.0))**ma
      kappae_y_PML_1(j) = 1.0+(kappa_y_max-1.0)*   &
                                 ((nyPML_1 - j) / (nyPML_1 - 1.0))**m
      be_y_1(j) = EXP(-(sige_y_PML_1(j) / kappae_y_PML_1(j) +   &
                                 alphae_y_PML_1(j))*dt/epsO)
      if ((sige_y_PML_1(j) == 0.0) .and.        &
         (alphae_y_PML_1(j) == 0.0) .and. (j == nyPML_1)) then
         ce_y_1(j) = 0.0
      else
         ce_y_1(j) = sige_y_PML_1(j)*(be_y_1(j)-1.0)/       &
            (sige_y_PML_1(j)+kappae_y_PML_1(j)*alphae_y_PML_1(j)) &
            / kappae_y_PML_1(j)
      endif
   ENDDO
   DO j = 1,nyPML_1-1
      sigh_y_PML_1(j) = sig_y_max * ( (nyPML_1 - j - 0.5)/(nyPML_1-1.0))**m
      alphah_y_PML_1(j) = alpha_y_max*((j-0.5)/(nyPML_1-1.0))**ma
      kappah_y_PML_1(j) = 1.0+(kappa_y_max-1.0)*   &
                            ((nyPML_1 - j - 0.5) / (nyPML_1 - 1.0))**m
      bh_y_1(j) = EXP(-(sigh_y_PML_1(j) / kappah_y_PML_1(j) +   &
                                 alphah_y_PML_1(j))*dt/epsO)
      ch_y_1(j) = sigh_y_PML_1(j)*(bh_y_1(j)-1.0)/      &
                  (sigh_y_PML_1(j)+kappah_y_PML_1(j)*alphah_y_PML_1(j)) &
                  / kappah_y_PML_1(j)
   ENDDO
   DO j = 1,nyPML_2
      sige_y_PML_2(j) = sig_y_max * ( (nyPML_2 - j ) / (nyPML_2 - 1.0) )**m
      alphae_y_PML_2(j) = alpha_y_max*((j-1)/(nyPML_2-1.0))**ma
      kappae_y_PML_2(j) = 1.0+(kappa_y_max-1.0)*   &
                                 ((nyPML_2 - j) / (nyPML_2 - 1.0))**m
      be_y_2(j) = EXP(-(sige_y_PML_2(j) / kappae_y_PML_2(j) +   &
                                 alphae_y_PML_2(j))*dt/epsO)
      if ((sige_y_PML_2(j) == 0.0) .and.        &
         (alphae_y_PML_2(j) == 0.0) .and. (j == nyPML_2)) then
         ce_y_2(j) = 0.0
      else
         ce_y_2(j) = sige_y_PML_2(j)*(be_y_2(j)-1.0)/       &
            (sige_y_PML_2(j)+kappae_y_PML_2(j)*alphae_y_PML_2(j)) &
            / kappae_y_PML_2(j)
      endif
   ENDDO
   DO j = 1,nyPML_2-1
      sigh_y_PML_2(j) = sig_y_max * ( (nyPML_2 - j - 0.5)/(nyPML_2-1.0))**m
      alphah_y_PML_2(j) = alpha_y_max*((j-0.5)/(nyPML_2-1.0))**ma
      kappah_y_PML_2(j) = 1.0+(kappa_y_max-1.0)*   &
                            ((nyPML_2 - j - 0.5) / (nyPML_2 - 1.0))**m
      bh_y_2(j) = EXP(-(sigh_y_PML_2(j) / kappah_y_PML_2(j) +   &
                                 alphah_y_PML_2(j))*dt/epsO)
      ch_y_2(j) = sigh_y_PML_2(j)*(bh_y_2(j)-1.0)/      &
                  (sigh_y_PML_2(j)+kappah_y_PML_2(j)*alphah_y_PML_2(j)) &
                  / kappah_y_PML_2(j)
   ENDDO

   DO k = 1,nzPML_1
      sige_z_PML_1(k) = sig_z_max * ( (nzPML_1 - k ) / (nzPML_1 - 1.0) )**m
      alphae_z_PML_1(k) = alpha_z_max*((k-1)/(nzPML_1-1.0))**ma
      kappae_z_PML_1(k) = 1.0+(kappa_z_max-1.0)*   &
                                 ((nzPML_1 - k) / (nzPML_1 - 1.0))**m
      be_z_1(k) = EXP(-(sige_z_PML_1(k) / kappae_z_PML_1(k) +   &
                                 alphae_z_PML_1(k))*dt/epsO)
      if ((sige_z_PML_1(k) == 0.0) .and.        &
         (alphae_z_PML_1(k) == 0.0) .and. (k == nzPML_1)) then
         ce_z_1(k) = 0.0
      else
         ce_z_1(k) = sige_z_PML_1(k)*(be_z_1(k)-1.0)/       &
            (sige_z_PML_1(k)+kappae_z_PML_1(k)*alphae_z_PML_1(k)) &
            / kappae_z_PML_1(k)
      endif
   ENDDO
   DO k = 1,nzPML_1-1
      sigh_z_PML_1(k) = sig_z_max * ( (nzPML_1 - k - 0.5)/(nzPML_1-1.0))**m
      alphah_z_PML_1(k) = alpha_z_max*((k-0.5)/(nzPML_1-1.0))**ma
      kappah_z_PML_1(k) = 1.0+(kappa_z_max-1.0)*   &
                            ((nzPML_1 - k - 0.5) / (nzPML_1 - 1.0))**m
      bh_z_1(k) = EXP(-(sigh_z_PML_1(k) / kappah_z_PML_1(k) +   &
                                 alphah_z_PML_1(k))*dt/epsO)
      ch_z_1(k) = sigh_z_PML_1(k)*(bh_z_1(k)-1.0)/      &
                  (sigh_z_PML_1(k)+kappah_z_PML_1(k)*alphah_z_PML_1(k)) &
                  / kappah_z_PML_1(k)
   ENDDO

   DO k = 1,nzPML_2
      sige_z_PML_2(k) = sig_z_max * ( (nzPML_2 - k ) / (nzPML_2 - 1.0) )**m
      alphae_z_PML_2(k) = alpha_z_max*((k-1)/(nzPML_2-1.0))**ma
      kappae_z_PML_2(k) = 1.0+(kappa_z_max-1.0)*   &
                                 ((nzPML_2 - k) / (nzPML_2 - 1.0))**m
      be_z_2(k) = EXP(-(sige_z_PML_2(k) / kappae_z_PML_2(k) +   &
                                 alphae_z_PML_2(k))*dt/epsO)
      if ((sige_z_PML_2(k) == 0.0) .and.        &
         (alphae_z_PML_2(k) == 0.0) .and. (k == nzPML_2)) then
         ce_z_2(k) = 0.0
      else
         ce_z_2(k) = sige_z_PML_2(k)*(be_z_2(k)-1.0)/       &
            (sige_z_PML_2(k)+kappae_z_PML_2(k)*alphae_z_PML_2(k)) &
            / kappae_z_PML_2(k)
      endif
   ENDDO
   DO k = 1,nzPML_2-1
      sigh_z_PML_2(k) = sig_z_max * ( (nzPML_2 - k - 0.5)/(nzPML_2-1.0))**m
      alphah_z_PML_2(k) = alpha_z_max*((k-0.5)/(nzPML_2-1.0))**ma
      kappah_z_PML_2(k) = 1.0+(kappa_z_max-1.0)*   &
                            ((nzPML_2 - k - 0.5) / (nzPML_2 - 1.0))**m
      bh_z_2(k) = EXP(-(sigh_z_PML_2(k) / kappah_z_PML_2(k) +   &
                                 alphah_z_PML_2(k))*dt/epsO)
      ch_z_2(k) = sigh_z_PML_2(k)*(bh_z_2(k)-1.0)/      &
                  (sigh_z_PML_2(k)+kappah_z_PML_2(k)*alphah_z_PML_2(k)) &
                  / kappah_z_PML_2(k)
   ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  FILL IN UPDATING COEFFICIENTS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   DA = 1.0
   DB = (dt/(muO)) 

   DO i = 1,Imax
      DO j = 1,Jmax
         DO k = 1,Kmax-1
            CA(i,j,k) = (1.0 - sig(i,j,k)*dt / (2.0*eps(i,j,k))) /        &
               (1.0 + sig(i,j,k) * dt / (2.0*eps(i,j,k)))
            CB(i,j,k) = (dt/(eps(i,j,k))) /                            &
               (1.0 + sig(i,j,k)*dt / (2.0*eps(i,j,k)))
            ENDDO
	ENDDO
   ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  FILL IN DENOMINATORS FOR FIELD UPDATES
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ii = nxPML_2-1
   DO i = 1,Imax-1
      if (i <= nxPML_1-1) then
         den_hx(i) = 1.0/(kappah_x_PML_1(i)*dx)
      elseif (i >= Imax+1-nxPML_2) then
         den_hx(i) = 1.0/(kappah_x_PML_2(ii)*dx)
         ii = ii-1
      else
         den_hx(i) = 1.0/dx
      endif
   ENDDO
   jj = nyPML_2-1
   DO j = 1,Jmax-1
      if (j <= nyPML_1-1) then
         den_hy(j) = 1.0/(kappah_y_PML_1(j)*dy)
      elseif (j >= Jmax+1-nyPML_2) then
         den_hy(j) = 1.0/(kappah_y_PML_2(jj)*dy)
         jj = jj-1
      else
         den_hy(j) = 1.0/dy
      endif
   ENDDO
   kk = nzPML_2-1
   DO k = 1,Kmax-1
      if (k <= nzPML_1-1) then
         den_hz(k) = 1.0/(kappah_z_PML_1(k)*dz)
      elseif (k >= Kmax+1-nzPML_2) then
         den_hz(k) = 1.0/(kappah_z_PML_2(kk)*dz)
         kk = kk - 1
      else
         den_hz(k) = 1.0/dz
      endif
   ENDDO
   ii = nxPML_2
   DO i = 1,Imax-1
      if (i <= nxPML_1) then
         den_ex(i) = 1.0/(kappae_x_PML_1(i)*dx)
      elseif (i >= Imax+1-nxPML_2) then
         den_ex(i) = 1.0/(kappae_x_PML_2(ii)*dx)
         ii = ii-1
      else
         den_ex(i) = 1.0/dx
      endif
   ENDDO
   jj = nyPML_2
   DO j = 1,Jmax-1
      if (j <= nyPML_1) then
         den_ey(j) = 1.0/(kappae_y_PML_1(j)*dy)
      elseif (j >= Jmax+1-nyPML_2) then
         den_ey(j) = 1.0/(kappae_y_PML_2(jj)*dy)
         jj = jj-1
      else
         den_ey(j) = 1.0/dy
      endif
   ENDDO
   kk = nzPML_2
   DO k = 1,Kmax-1
      if (k <= nzPML_1) then
         den_ez(k) = 1.0/(kappae_z_PML_1(k)*dz)
      elseif (k >= Kmax+1-nzPML_2) then
         den_ez(k) = 1.0/(kappae_z_PML_2(kk)*dz)
         kk = kk - 1
      else
         den_ez(k) = 1.0/dz
      endif
   ENDDO

!.:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  BEGIN TIME STEP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  write(*,*)"begin time-stepping"
  DO n = 1,nmax
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  UPDATE Hx
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   DO k = 1,Kmax-1
      DO i = 1,Imax-1
	   DO j = 1,Jmax-1
	      Hx(i,j,k) = DA * Hx(i,j,k) + DB *       &
			( (Ez(i,j,k) - Ez(i,j+1,k))*den_hy(j)  +    &
			  (Ey(i,j,k+1) - Ey(i,j,k))*den_hz(k) )
	   ENDDO
	ENDDO 
      DO i = 1,Imax-1
!.....................................................................
!  PML for bottom Hx, j-direction
!.....................................................................
         DO j = 1,nyPML_1-1
  	      psi_Hxy_1(i,j,k) = bh_y_1(j)*psi_Hxy_1(i,j,k)                 &
	 			+ ch_y_1(j) *(Ez(i,j,k) - Ez(i,j+1,k))/dy
            Hx(i,j,k) = Hx(i,j,k) + DB*psi_Hxy_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Hx, j-direction
!.....................................................................
         jj = nyPML_2-1
         DO j = Jmax+1-nyPML_2,Jmax-1
  	      psi_Hxy_2(i,jj,k) = bh_y_2(jj)*psi_Hxy_2(i,jj,k)       &
	 			+ ch_y_2(jj) *(Ez(i,j,k) -       &
                                Ez(i,(j+1),k))/dy
            Hx(i,j,k) = Hx(i,j,k) + DB*psi_Hxy_2(i,jj,k)
            jj = jj-1
         ENDDO
	ENDDO
   ENDDO
   DO i = 1,Imax-1
      DO j = 1,Jmax-1
!.....................................................................
!  PML for bottom Hx, k-direction
!.....................................................................
         DO k = 1,nzPML_1-1
  	      psi_Hxz_1(i,j,k) = bh_z_1(k)*psi_Hxz_1(i,j,k)            &
	 			+ ch_z_1(k) *(Ey(i,j,k+1) - Ey(i,j,k))/dz
            Hx(i,j,k) = Hx(i,j,k) + DB*psi_Hxz_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Hx, k-direction
!.....................................................................
         kk = nzPML_2-1
         DO k = Kmax+1-nzPML_2,Kmax-1
  	      psi_Hxz_2(i,j,kk) = bh_z_2(kk)*psi_Hxz_2(i,j,kk)             &
	 			+ ch_z_2(kk) *(Ey(i,j,k+1) -       &
                                Ey(i,j,k))/dz
            Hx(i,j,k) = Hx(i,j,k) + DB*psi_Hxz_2(i,j,kk)
           kk = kk-1
         ENDDO
	ENDDO
   ENDDO

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  UPDATE Hy
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   DO k = 1,Kmax-1
      DO i = 1,Imax-1
	   DO j = 1,Jmax-1
            Hy(i,j,k) = DA * Hy(i,j,k) + DB *             &
			( (Ez(i+1,j,k) - Ez(i,j,k))*den_hx(i) +      &
			  (Ex(i,j,k) - Ex(i,j,k+1))*den_hz(k) )
	   ENDDO 
      ENDDO
      DO j = 1,Jmax-1
!.....................................................................
!  PML for bottom Hy, i-direction
!.....................................................................
         DO i = 1,nxPML_1-1
	      psi_Hyx_1(i,j,k) = bh_x_1(i)*psi_Hyx_1(i,j,k)                   &
				+ ch_x_1(i)*(Ez(i+1,j,k) - Ez(i,j,k))/dx
	      Hy(i,j,k) = Hy(i,j,k) + DB*psi_Hyx_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Hy, i-direction
!.....................................................................
         ii = nxPML_2-1
         DO i = Imax+1-nxPML_2,Imax-1
	      psi_Hyx_2(ii,j,k) = bh_x_2(ii)*psi_Hyx_2(ii,j,k)     &
				+ ch_x_2(ii)*(Ez(i+1,j,k) -    &
                                Ez(i,j,k))/dx
	      Hy(i,j,k) = Hy(i,j,k) + DB*psi_Hyx_2(ii,j,k)
            ii = ii-1
	   ENDDO
     ENDDO
   ENDDO
   DO i = 1,Imax-1
      DO j = 1,Jmax-1
!.....................................................................
!  PML for bottom Hy, k-direction
!.....................................................................
         DO k = 1,nzPML_1-1
	      psi_Hyz_1(i,j,k) = bh_z_1(k)*psi_Hyz_1(i,j,k)          &
				+ ch_z_1(k)*(Ex(i,j,k) - Ex(i,j,k+1))/dz
	      Hy(i,j,k) = Hy(i,j,k) + DB*psi_Hyz_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Hy, k-direction
!.....................................................................
         kk = nzPML_2-1
         DO k = Kmax+1-nzPML_2,Kmax-1
	    psi_Hyz_2(i,j,kk) = bh_z_2(kk)*psi_Hyz_2(i,j,kk)               &
				+ ch_z_2(kk)*(Ex(i,j,k) -    &
                                Ex(i,j,k+1))/dz
	    Hy(i,j,k) = Hy(i,j,k) + DB*psi_Hyz_2(i,j,kk)
            kk = kk-1
         ENDDO
     ENDDO
   ENDDO

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  UPDATE Hz
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   DO k = 2,Kmax-1
      DO i = 1,Imax-1
        DO j = 1,Jmax-1
            Hz(i,j,k) = DA * Hz(i,j,k) + DB       &
                  * ((Ey(i,j,k) - Ey(i+1,j,k))*den_hx(i) +        &
			    (Ex(i,j+1,k) - Ex(i,j,k))*den_hy(j))
	   ENDDO
      ENDDO
      DO j = 1,Jmax-1
!.....................................................................
!  PML for bottom Hz, x-direction
!.....................................................................
         DO i = 1,nxPML_1-1
   	      psi_Hzx_1(i,j,k) = bh_x_1(i)*psi_Hzx_1(i,j,k)                 &
	 			+ ch_x_1(i) *(Ey(i,j,k) - Ey(i+1,j,k))/dx
	      Hz(i,j,k) = Hz(i,j,k) + DB*psi_Hzx_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Hz, x-direction
!.....................................................................
         ii = nxPML_2-1
         DO i = Imax+1-nxPML_2,Imax-1
   	      psi_Hzx_2(ii,j,k) = bh_x_2(ii)*psi_Hzx_2(ii,j,k)            &
	 			+ ch_x_2(ii) *(Ey(i,j,k) -       &
                                Ey(i+1,j,k))/dx
	      Hz(i,j,k) = Hz(i,j,k) + DB*psi_Hzx_2(ii,j,k)
            ii = ii-1
	   ENDDO
      ENDDO
      DO i = 1,Imax-1
!.....................................................................
!  PML for bottom Hz, y-direction
!.....................................................................
         DO j = 1,nyPML_1-1
            psi_Hzy_1(i,j,k) = bh_y_1(j)*psi_Hzy_1(i,j,k)                   &
				+ ch_y_1(j)*(Ex(i,j+1,k) - Ex(i,j,k))/dy
	      Hz(i,j,k) = Hz(i,j,k) + DB*psi_Hzy_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Hz, y-direction
!.....................................................................
         jj = nyPML_2-1
         DO j = Jmax+1-nyPML_2,Jmax-1
            psi_Hzy_2(i,jj,k) = bh_y_2(jj)*psi_Hzy_2(i,jj,k)               &
				+ ch_y_2(jj)*(Ex(i,j+1,k) -    &
                                Ex(i,j,k))/dy
	      Hz(i,j,k) = Hz(i,j,k) + DB*psi_Hzy_2(i,jj,k)
            jj = jj-1
	   ENDDO
      ENDDO
   ENDDO

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  UPDATE Ex
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   DO k = 2,Kmax-1
      DO i = 1,Imax-1
	   DO j = 2,Jmax-1
              IF (i >= istart .and. i <= iend .and. j >= jstart .and.  &
                   j <= jend .and. k >= kstart .and. k <= kend) THEN
                 Ex(i,j,k) = 0.0
              ELSE
	         Ex(i,j,k) = CA(i,j,k) * Ex(i,j,k) + CB(i,j,k) *       &
  			( (Hz(i,j,k) - Hz(i,j-1,k))*den_ey(j)  +    &
			  (Hy(i,j,k-1) - Hy(i,j,k))*den_ez(k) )
              ENDIF
	   ENDDO
	ENDDO 
      DO i = 1,Imax-1
!.....................................................................
!  PML for bottom Ex, j-direction
!.....................................................................
         DO j = 2,nyPML_1
  	      psi_Exy_1(i,j,k) = be_y_1(j)*psi_Exy_1(i,j,k)                 &
	 			+ ce_y_1(j) *(Hz(i,j,k) - Hz(i,j-1,k))/dy
            Ex(i,j,k) = Ex(i,j,k) + CB(i,j,k)*psi_Exy_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Ex, j-direction
!.....................................................................
         jj = nyPML_2
         DO j = Jmax+1-nyPML_2,Jmax-1
  	      psi_Exy_2(i,jj,k) = be_y_2(jj)*psi_Exy_2(i,jj,k)       &
	 			+ ce_y_2(jj) *(Hz(i,j,k) -       &
                                Hz(i,(j-1),k))/dy
            Ex(i,j,k) = Ex(i,j,k) + CB(i,j,k)*psi_Exy_2(i,jj,k)
            jj = jj-1
         ENDDO
	ENDDO
   ENDDO
   DO i = 1,Imax-1
      DO j = 2,Jmax-1
!.....................................................................
!  PML for bottom Ex, k-direction
!.....................................................................
         DO k = 2,nzPML_1
  	      psi_Exz_1(i,j,k) = be_z_1(k)*psi_Exz_1(i,j,k)                 &
	 			+ ce_z_1(k) *(Hy(i,j,k-1) - Hy(i,j,k))/dz
            Ex(i,j,k) = Ex(i,j,k) + CB(i,j,k)*psi_Exz_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Ex, k-direction
!.....................................................................
         kk = nzPML_2
         DO k = Kmax+1-nzPML_2,Kmax-1
  	      psi_Exz_2(i,j,kk) = be_z_2(kk)*psi_Exz_2(i,j,kk)             &
	 			+ ce_z_2(kk) *(Hy(i,j,k-1) -       &
                                Hy(i,j,k))/dz
            Ex(i,j,k) = Ex(i,j,k) + CB(i,j,k)*psi_Exz_2(i,j,kk)
            kk = kk-1
         ENDDO
	ENDDO
   ENDDO

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  UPDATE Ey
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   DO k = 2,Kmax-1
      DO i = 2,Imax-1
	   DO j = 1,Jmax-1
              IF (i >= istart .and. i <= iend .and. j >= jstart .and. &
                   j <= jend .and. k >= kstart .and. k <= kend) THEN
                 Ey(i,j,k) = 0.0
              ELSE
                 Ey(i,j,k) = CA(i,j,k) * Ey(i,j,k) + CB(i,j,k) *    &
			( (Hz(i-1,j,k) - Hz(i,j,k))*den_ex(i) +         &
			  (Hx(i,j,k) - Hx(i,j,k-1))*den_ez(k) )
              ENDIF
	   ENDDO 
      ENDDO
      DO j = 1,Jmax-1
!.....................................................................
!  PML for bottom Ey, i-direction
!.....................................................................
         DO i = 2,nxPML_1
	      psi_Eyx_1(i,j,k) = be_x_1(i)*psi_Eyx_1(i,j,k)                   &
				+ ce_x_1(i)*(Hz(i-1,j,k) - Hz(i,j,k))/dx
	      Ey(i,j,k) = Ey(i,j,k) + CB(i,j,k)*psi_Eyx_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Ey, i-direction
!.....................................................................
         ii = nxPML_2
         DO i = Imax+1-nxPML_2,Imax-1
	      psi_Eyx_2(ii,j,k) = be_x_2(ii)*psi_Eyx_2(ii,j,k)              &
				+ ce_x_2(ii)*(Hz(i-1,j,k) -    &
                               Hz(i,j,k))/dx
	      Ey(i,j,k) = Ey(i,j,k) + CB(i,j,k)*psi_Eyx_2(ii,j,k)
            ii = ii-1
	   ENDDO
     ENDDO
   ENDDO
   DO i = 2,Imax-1
      DO j = 1,Jmax-1
!.....................................................................
!  PML for bottom Ey, k-direction
!.....................................................................
         DO k = 2,nzPML_1
	      psi_Eyz_1(i,j,k) = be_z_1(k)*psi_Eyz_1(i,j,k)                   &
				+ ce_z_1(k)*(Hx(i,j,k) - Hx(i,j,k-1))/dz
	      Ey(i,j,k) = Ey(i,j,k) + CB(i,j,k)*psi_Eyz_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Ey, k-direction
!.....................................................................
         kk = nzPML_2
         DO k = Kmax+1-nzPML_2,Kmax-1
	    psi_Eyz_2(i,j,kk) = be_z_2(kk)*psi_Eyz_2(i,j,kk)               &
				+ ce_z_2(kk)*(Hx(i,j,k) -    &
                                Hx(i,j,k-1))/dz
	    Ey(i,j,k) = Ey(i,j,k) + CB(i,j,k)*psi_Eyz_2(i,j,kk)
            kk = kk-1
         ENDDO
     ENDDO
   ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  UPDATE Ez
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   DO k = 1,Kmax-1
      DO i = 2,Imax-1
         DO j = 2,Jmax-1
            Ez(i,j,k) = CA(i,j,k) * Ez(i,j,k) + CB(i,j,k)       &
                  * ((Hy(i,j,k) - Hy(i-1,j,k))*den_ex(i) +        &
			    (Hx(i,j-1,k) - Hx(i,j,k))*den_ey(j))
	   ENDDO
      ENDDO
      DO j = 2,Jmax-1
!.....................................................................
!  PML for bottom Ez, x-direction
!.....................................................................
         DO i = 2,nxPML_1
   	      psi_Ezx_1(i,j,k) = be_x_1(i)*psi_Ezx_1(i,j,k)             &
	 			+ ce_x_1(i) *(Hy(i,j,k) - Hy(i-1,j,k))/dx
	      Ez(i,j,k) = Ez(i,j,k) + CB(i,j,k)*psi_Ezx_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Ez, x-direction
!.....................................................................
         ii = nxPML_2
         DO i = Imax+1-nxPML_2,Imax-1
   	      psi_Ezx_2(ii,j,k) = be_x_2(ii)*psi_Ezx_2(ii,j,k)       &
	 			+ ce_x_2(ii) *(Hy(i,j,k) -       &
                                Hy(i-1,j,k))/dx
	      Ez(i,j,k) = Ez(i,j,k) + CB(i,j,k)*psi_Ezx_2(ii,j,k)
            ii = ii-1
	   ENDDO
      ENDDO
      DO i = 2,Imax-1
!.....................................................................
!  PML for bottom Ez, y-direction
!.....................................................................
         DO j = 2,nyPML_1
            psi_Ezy_1(i,j,k) = be_y_1(j)*psi_Ezy_1(i,j,k)            &
				+ ce_y_1(j)*(Hx(i,j-1,k) - Hx(i,j,k))/dy
	      Ez(i,j,k) = Ez(i,j,k) + CB(i,j,k)*psi_Ezy_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Ez, y-direction
!.....................................................................
         jj = nyPML_2
         DO j = Jmax+1-nyPML_2,Jmax-1
            psi_Ezy_2(i,jj,k) = be_y_2(jj)*psi_Ezy_2(i,jj,k)         &
				+ ce_y_2(jj)*(Hx(i,j-1,k) -    &
                                Hx(i,j,k))/dy
	      Ez(i,j,k) = Ez(i,j,k) + CB(i,j,k)*psi_Ezy_2(i,jj,k)
            jj = jj-1
	   ENDDO
      ENDDO
   ENDDO

!-----------------------------------------------------------------------
!   SOURCE
!-----------------------------------------------------------------------
   i = isource
   j = jsource
   k = ksource
   source = -2.0*((n*dt-tO)/tw) * exp(-((n*dt-tO)/tw)**2.0)
   Ez(i,j,k) = Ez(i,j,k) - CB(i,j,k)*source

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Update P_sum_fdtd
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Px = 0.0
Py = 0.0
Pz = 0.0
P_sum = 0.0

do i = i_start,i_end !Add Poynting magnitudes in 3D
 do j = j_start,j_end
  do k = k_start,k_end
   Px = Ey(i,j,k)*Hz(i,j,k) - Ez(i,j,k)*Hy(i,j,k)
   Py = Ez(i,j,k)*Hx(i,j,k) - Ex(i,j,k)*Hz(i,j,k)
   Pz = Ex(i,j,k)*Hy(i,j,k) - Ey(i,j,k)*Hx(i,j,k)
   P_sum = P_sum + sqrt(Px**2 + Py**2 + Pz**2)
  enddo
 enddo
enddo

P_sum_fdtd = P_sum_fdtd + P_sum/res_array(a) !res is the integration equalizer

   ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  END TIME STEP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!.:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:.
    WRITE(*,*)"done time-stepping"
    
!-------------------------------------------------------------------
!---------------------- Start of Array Deallocation ------------------
!-------------------------------------------------------------------
deallocate(Ez(Imax, Jmax, Kmax-1), CA(Imax, Jmax, Kmax-1), CB(Imax, Jmax, Kmax-1), sig(Imax, Jmax, Kmax-1), eps(Imax, Jmax, Kmax-1))

deallocate(Hy(Imax-1, Jmax, Kmax-1))

deallocate(Hx(Imax,Jmax-1, Kmax-1))

deallocate(Hz(Imax-1, Jmax-1, Kmax))

deallocate(Ex(Imax-1, Jmax, Kmax))

deallocate(Ey(Imax,Jmax-1, Kmax))

!PML

deallocate(psi_Ezx_1(nxPML_1,Jmax,Kmax))

deallocate(psi_Ezx_2(nxPML_2,Jmax,Kmax))

deallocate(psi_Hyx_1(nxPML_1-1,Jmax,Kmax))

deallocate(psi_Hyx_2(nxPML_2-1,Jmax,Kmax))

deallocate(psi_Ezy_1(Imax,nyPML_1,Kmax))                               

deallocate(psi_Ezy_2(Imax,nyPML_2,Kmax))

deallocate(psi_Hxy_1(Imax,nyPML_1-1,Kmax))                               

deallocate(psi_Hxy_2(Imax,nyPML_2-1,Kmax))

deallocate(psi_Hxz_1(Imax,Jmax-1,nzPML_1-1))

deallocate(psi_Hxz_2(Imax,Jmax-1,nzPML_2-1))

deallocate(psi_Hyz_1(Imax-1,Jmax,nzPML_1-1))

deallocate(psi_Hyz_2(Imax-1,Jmax,nzPML_2-1))

deallocate(psi_Exz_1(Imax-1,Jmax,nzPML_1))

deallocate(psi_Exz_2(Imax-1,Jmax,nzPML_2))

deallocate(psi_Eyz_1(Imax,Jmax-1,nzPML_1))

deallocate(psi_Eyz_2(Imax,Jmax-1,nzPML_2))

deallocate(psi_Hzx_1(nxPML_1-1,Jmax-1,Kmax))

deallocate(psi_Eyx_1(nxPML_1,Jmax-1,Kmax))

deallocate(psi_Hzx_2(nxPML_2-1,Jmax-1,Kmax))

deallocate(psi_Eyx_2(nxPML_2,Jmax-1,Kmax))

deallocate(psi_Hzy_1(Imax-1,nyPML_1-1,Kmax))                              

deallocate(psi_Exy_1(Imax-1,nyPML_1,Kmax))                               

deallocate(psi_Hzy_2(Imax-1,nyPML_2-1,Kmax))

deallocate(psi_Exy_2(Imax-1,nyPML_2,Kmax))

deallocate(be_x_1(nxPML_1), ce_x_1(nxPML_1), alphae_x_PML_1(nxPML_1), sige_x_PML_1(nxPML_1), kappae_x_PML_1(nxPML_1))

deallocate(bh_x_1(nxPML_1-1), ch_x_1(nxPML_1-1), alphah_x_PML_1(nxPML_1-1), sigh_x_PML_1(nxPML_1-1), kappah_x_PML_1(nxPML_1-1))

deallocate(be_x_2(nxPML_2), ce_x_2(nxPML_2), alphae_x_PML_2(nxPML_2), sige_x_PML_2(nxPML_2), kappae_x_PML_2(nxPML_2))

deallocate(bh_x_2(nxPML_2-1), ch_x_2(nxPML_2-1), alphah_x_PML_2(nxPML_2-1), sigh_x_PML_2(nxPML_2-1), kappah_x_PML_2(nxPML_2-1))

deallocate(be_y_1(nyPML_1), ce_y_1(nyPML_1), alphae_y_PML_1(nyPML_1), sige_y_PML_1(nyPML_1), kappae_y_PML_1(nyPML_1))

deallocate(bh_y_1(nyPML_1-1), ch_y_1(nyPML_1-1), alphah_y_PML_1(nyPML_1-1), sigh_y_PML_1(nyPML_1-1), kappah_y_PML_1(nyPML_1-1))

deallocate(be_y_2(nyPML_2), ce_y_2(nyPML_2), alphae_y_PML_2(nyPML_2), sige_y_PML_2(nyPML_2), kappae_y_PML_2(nyPML_2))

deallocate(bh_y_2(nyPML_2-1), ch_y_2(nyPML_2-1), alphah_y_PML_2(nyPML_2-1), sigh_y_PML_2(nyPML_2-1), kappah_y_PML_2(nyPML_2-1))

deallocate(be_z_1(nzPML_1), ce_z_1(nzPML_1), alphae_z_PML_1(nzPML_1), sige_z_PML_1(nzPML_1), kappae_z_PML_1(nzPML_1))

deallocate(bh_z_1(nzPML_1-1), ch_z_1(nzPML_1-1), alphah_z_PML_1(nzPML_1-1), sigh_z_PML_1(nzPML_1-1), kappah_z_PML_1(nzPML_1-1))

deallocate(be_z_2(nzPML_2), ce_z_2(nzPML_2), alphae_z_PML_2(nzPML_2), sige_z_PML_2(nzPML_2), kappae_z_PML_2(nzPML_2))

deallocate(bh_z_2(nzPML_2-1), ch_z_2(nzPML_2-1), alphah_z_PML_2(nzPML_2-1), sigh_z_PML_2(nzPML_2-1), kappah_z_PML_2(nzPML_2-1))

!denominators for the update equations

deallocate(den_ex(Imax-1), den_hx(Imax-1))

deallocate(den_ey(Jmax-1), den_hy(Jmax-1))

deallocate(den_ez(Kmax-1), den_hz(Kmax-1))


end function fdtd3D_CPML

end !Main
