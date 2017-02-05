function Convergence_Collect(D, Ex, Ey, Ez, Hx, Hy, Hz, &
                             i_start, i_end, j_start, j_end, k_start, k_end) result(P_sum)
                             
integer, intent(in) :: D, i_start, i_end, j_start, j_end, k_start, k_end

if(D == 3)then
double precision, intent(in) :: Ex(:,:,:), Ey(:,:,:), Ez(:,:,:), Hx(:,:,:), Hy(:,:,:), Hz(:,:,:)
elseif(D == 2)then
double precision, intent(in) :: Ex(:,:), Ey(:,:), Ez(:,:), Hx(:,:), Hy(:,:), Hz(:,:)
elseif(D == 1)then
double precision, intent(in) :: Ex(:), Ey(:), Ez(:), Hx(:), Hy(:), Hz(:)
endif

double precision, intent(out) :: P_sum
double precision :: Px, Py, Pz
integer :: i,j,k

P_sum = 0.0
Px = 0.0
Py = 0.0
Pz = 0.0

if(D == 3)then !add Poynting magnitudes in 3D
 do i = i_start,i_end
  do j = j_start,j_end
   do k = k_start,k_end
    Px = Ey(i,j,k)*Hz(i,j,k) - Ez(i,j,k)*Hy(i,j,k)
    Py = Ez(i,j,k)*Hx(i,j,k) - Ex(i,j,k)*Hz(i,j,k)
    Pz = Ex(i,j,k)*Hy(i,j,k) - Ey(i,j,k)*Hx(i,j,k)
    P_sum = P_sum + sqrt(Px**2 + Py**2 + Pz**2)
   enddo
  enddo
 enddo
endif

if(D == 2)then !add Poynting magnitudes in 2D
 do i = i_start,i_end
  do j = j_start,j_end
   Px = Ey(i,j)*Hz(i,j) - Ez(i,j)*Hy(i,j)
   Py = Ez(i,j)*Hx(i,j) - Ex(i,j)*Hz(i,j)
   Pz = Ex(i,j)*Hy(i,j) - Ey(i,j)*Hx(i,j)
   P_sum = P_sum + sqrt(Px**2 + Py**2 + Pz**2) 
  enddo
 enddo
endif

if(D == 1)then !add Poynting magnitudes in 1D
 do i = i_start,i_end
  Px = Ey(i)*Hz(i) - Ez(i)*Hy(i)
  Py = Ez(i)*Hx(i) - Ex(i)*Hz(i)
  Pz = Ex(i)*Hy(i) - Ey(i)*Hx(i)
  P_sum = P_sum + sqrt(Px**2 + Py**2 + Pz**2)
 enddo
endif
 
if( D <= 0.or.D => 4 )then !returns a marked P_sum if dimension is non-Euclidean
 P_sum = -1.0
endif

end function Convergence_Collect
