! This file contains multiple schemes I need to use

module scheme

  implicit none
  real(kind=8)::pi
  save

contains

  !============two dimension linear convection--around 1============


  subroutine two_d_linear_adv_FTBS(c,dt,t,dx,x,dy,y,u0,u)


    implicit none
    real(kind=8),intent(in)::c,dt,dx,t,x,dy,y
    integer::itime,ispace,jspace
    real(kind=8),dimension(:,:),intent(in)::u0
    real(kind=8),dimension(size(u0,1),size(u0,2)),intent(out)::u
    real(kind=8),dimension(size(u0,1),size(u0,2))::ut_last
    u=u0
    print *, nint(t/dt)
    print *, nint(y/dy),nint(x/dx),c*dt/dx,c*dt/dy
    Do itime=1,nint(t/dt)+1
       ut_last=u
       Do jspace=2,nint(y/dy)

          Do ispace=2,nint(x/dx)

             u(jspace, ispace)=ut_last(jspace,ispace)-c*dt/dx*(ut_last(jspace,ispace)-ut_last(jspace,ispace-1)) &
                  -c*dt/dy*(ut_last(jspace,ispace)-ut_last(jspace-1,ispace))

          enddo
       enddo
       u(jspace, ispace)=ut_last(jspace,ispace)-c*dt/dx*(ut_last(jspace,ispace)-ut_last(jspace,ispace-1)) &
            -c*dt/dy*(ut_last(jspace,ispace)-ut_last(jspace-1,ispace))


    enddo

    !============test============


    ! u(1,1)=4
    ! u(3,2)=7
    ! print *, u(:,1)-2*u(:,2)
  end subroutine two_d_linear_adv_FTBS



end module scheme
