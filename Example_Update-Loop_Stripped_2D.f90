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

 do j = 1,N_loc
  do i = 1,Nx
   if(x(i) == x_source .and. y(j) == y_source)then
    Hz(i,j) = Hz(i,j) + pulse(n)
   endif
   if(x(i) == x_detect .and. y(j) == y_detect)then
    P_sum = P_sum + (Hz(i,j) + Hz(i-1,j) + Hz(i,j-1) + Hz(i-1,j-1) )/4.0
   endif
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

enddo !Nt
