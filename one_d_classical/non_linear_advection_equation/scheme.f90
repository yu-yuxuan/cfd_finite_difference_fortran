! This file contains multiple schemes I need to use

module scheme

  implicit none
  real(kind=8)::pi
  save

contains
  !============one_d_advection_burgers:FTBS============
  subroutine one_d_adv_burgers_FTBS(dt,t,dx,x,u0,u)
    implicit none
    real(kind=8),intent(in)::dt,t,dx,x
    integer::ispace,itime
    real(kind=8),intent(in),dimension(:)::u0
    real(kind=8),intent(out),dimension(size(u0))::u
    real(kind=8),dimension(size(u0))::ut_last

    u=u0
    print *, nint(t/dt)
    print *, nint(x/dx)
    Do itime=1,nint(t/dt)+1
       ut_last=u
       Do ispace=2,nint(x/dx)
          u(ispace)=ut_last(ispace)-ut_last(ispace)*dt/dx*(ut_last(ispace)-ut_last(ispace-1))
       enddo
       u(1)=ut_last(1)-ut_last(1)*dt/dx*(ut_last(1)-ut_last(nint(x/dx)))
       u(nint(x/dx)+1)=u(1)
    enddo
  end subroutine one_d_adv_burgers_FTBS

end module scheme
