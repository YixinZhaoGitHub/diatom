module aggcal
    implicit none
integer, parameter :: dp=kind(0.d0) ! double precision

contains
subroutine getagg(ngrid,D1, D2, D3, D4, D5, Dd)
! Arguments declarations
integer :: j
integer, intent(in):: ngrid
real(dp),intent(inout) :: D1(ngrid)
real(dp),intent(inout) :: D2(ngrid)
real(dp),intent(inout) :: D3(ngrid)
real(dp),intent(inout) :: D4(ngrid)
real(dp),intent(inout) :: D5(ngrid)
real(dp),intent(inout) :: Dd(ngrid)

  do j = 1, ngrid
      D1(j)= 200*D1(j)
      D2(j)= 200*D2(j)
      D3(j)= 200*D3(j)
      D4(j)= 200*D4(j)
      D5(j)= 200*D5(j)
      Dd(j)= 200*Dd(j)
  end do

end subroutine getagg
end module aggcal


