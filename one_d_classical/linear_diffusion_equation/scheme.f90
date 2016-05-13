! This file contains multiple schemes I need to use

module scheme

  implicit none
  real(kind=8)::pi
  save

contains
  !============one_d_diffusion:FTCS============


  subroutine one_d_diff_FTCS (nu,dt,t,dx,x,u0,u)
    implicit none
    real(kind=8), intent(in) :: nu,dt,t,dx,x
    integer::ispace,itime
    real(kind=8),intent(in),dimension(:)::u0
    real(kind=8),intent(out),dimension(size(u0))::u
    real(kind=8),dimension(size(u0))::ut_last

    u=u0
    print *, nint(t/dt)
    print *, nint(x/dx),nu*dt/(dx**2)
    Do itime=1,nint(t/dt)+1
       ut_last=u
       Do ispace=2,nint(x/dx)

          u(ispace)=ut_last(ispace)+nu*dt/(dx**2)*(ut_last(ispace+1)-2*ut_last(ispace)+ut_last(ispace-1))

       enddo
       u(1)=ut_last(1)+nu*dt/(dx**2)*(ut_last(2)-2*ut_last(1)+ ut_last(nint(x/dx)))

       u(nint(x/dx)+1)=u(1)

    enddo


  end subroutine one_d_diff_FTCS


  !============one_d_diff:BTCS============


  subroutine one_d_diff_implicit (nu,dt,t,dx,x,u0,u,epsilon)
    implicit none
    real(kind=8), intent(in) :: nu,dt,t,dx,x,epsilon
    real(kind=8)::sum,s_nu
    integer::ispace,itime
    real(kind=8),intent(in),dimension(:)::u0
    real(kind=8),intent(out),dimension(size(u0))::u
    real(kind=8),dimension(size(u0))::ut_last,utemp_J
    s_nu=nu*dt/(dx**2)
    u=u0
    print *, nint(t/dt)
    print *, nint(x/dx),nu*dt/(dx**2)

    Do itime=1,nint(t/dt)+1
       ut_last=u
       Do
          sum=0
          utemp_J=u
          Do ispace=2,nint(x/dx)

             u(ispace)=ut_last(ispace)/(1+2*s_nu)+s_nu*(utemp_J(ispace+1)+utemp_J(ispace-1))/(1+2*s_nu)

             sum=(u(ispace)-utemp_J(ispace))**2+sum
          enddo
          u(1)=ut_last(1)/(1+2*s_nu)+s_nu*(utemp_J(2)+ utemp_J(nint(x/dx)))/(1+2*s_nu)

          u(nint(x/dx)+1)=u(1)

          sum=(u(1)-utemp_J(1))**2+(u(nint(x/dx)+1)-utemp_J(nint(x/dx)+1))+sum
          sum=(sum/(ispace+1))**(1.d0/2)
          if(sum<epsilon) exit
       enddo

    enddo
  end subroutine one_d_diff_implicit

end module scheme
