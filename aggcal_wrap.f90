module aggcal_wrap
 use iso_c_binding, only: c_double, c_int
 use aggcal

  implicit none

contains

  subroutine f_getagg(ngrid,D1, D2, D3, D4, D5, Dd) bind(c)

    integer(c_int), intent(in), value:: ngrid
    real(c_double),  intent(inout):: D1(ngrid),D2(ngrid),D3(ngrid),D4(ngrid),D5(ngrid),Dd(ngrid)


    call getagg(ngrid,D1, D2, D3, D4, D5, Dd)

  end subroutine f_getagg

end module aggcal_wrap

