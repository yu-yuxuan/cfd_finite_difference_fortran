! This file contains multiple schemes I need to use

module scheme

  implicit none
  real(kind=8)::pi
  save

contains


  !============one_d_linear_advection:FTBS:boundary=constant============


  subroutine one_d_linear_adv_FTBS(c,dt,t,dx,x,u0,u)

    implicit none
    real(kind=8),intent(in)::c,dt,dx,t,x
    integer::itime,ispace
    real(kind=8),dimension(:),intent(in)::u0
    real(kind=8),dimension(size(u0)),intent(out)::u
    real(kind=8),dimension(size(u0))::ut_last
    u=u0
    ! print *, nint(t/dt)
    ! print *, nint(x/dx),c*dt/dx
    Do itime=1,nint(t/dt)+1
       ut_last=u
       Do ispace=2,nint(x/dx)

          u(ispace)=ut_last(ispace)-c*dt/dx*(ut_last(ispace)-ut_last(ispace-1))

       enddo
       u(1)=ut_last(1)-c*dt/dx*(ut_last(1)-ut_last(nint(x/dx)))
       u(nint(x/dx)+1)=u(1)
    enddo

  end subroutine one_d_linear_adv_FTBS



  !============one_d_linear_advection:FTCS:boundary=constant============


  subroutine one_d_linear_adv_FTCS(c,dt,t,dx,x,u0,u)

    implicit none
    real(kind=8),intent(in)::c,dt,dx,t,x
    integer::itime,ispace
    real(kind=8),dimension(:),intent(in)::u0
    real(kind=8),dimension(size(u0)),intent(out)::u
    real(kind=8),dimension(size(u0))::ut_last
    u=u0
    ! print *, nint(t/dt)
    ! print *, nint(x/dx),c*dt/dx
    Do itime=1,nint(t/dt)+1
       ut_last=u
       Do ispace=2,nint(x/dx)

          u(ispace)=ut_last(ispace)-c*dt/dx/2*(ut_last(ispace+1)-ut_last(ispace-1))

       enddo
       u(1)=ut_last(1)-c*dt/dx/2*(ut_last(2)-ut_last(nint(x/dx)))
       u(nint(x/dx)+1)=u(1)

    enddo

  end subroutine one_d_linear_adv_FTCS



  !============one_d_linear_advection_implicit:BTCS:boundary=constant============


  subroutine one_d_linear_adv_BTCS(c,dt,t,dx,x,u0,u,epsilon)

    implicit none
    real(kind=8),intent(in)::c,dt,dx,t,x,epsilon
    real(kind=8)::sum
    integer::itime,ispace
    real(kind=8),dimension(:),intent(in)::u0
    real(kind=8),dimension(size(u0)),intent(out)::u
    real(kind=8),dimension(size(u0))::ut_last,utemp_J
    !  real(kind=8),external::norm2
    u=u0
    ! print *, nint(t/dt)
    ! print *, nint(x/dx),c*dt/dx
    Do itime=1,nint(t/dt)+1

       !============Jacobi iteration============
       ut_last=u

       do
          sum=0
          utemp_J=u

          do ispace=2,nint(x/dx)  !get the Xn+1 by Xn

             u(ispace)=ut_last(ispace)-c*dt/2.d0/dx*(utemp_J(ispace+1)-utemp_J(ispace-1))
             sum=(u(ispace)-utemp_J(ispace))**2+sum
             !    print * ,"used",u(ispace),utemp_J(ispace)
          enddo
          u(1)=ut_last(1)-c*dt/2.d0/dx*(utemp_J(2)-utemp_J(nint(x/dx)))
          u(nint(x/dx)+1)=u(1)
          sum=(u(1)-utemp_J(1))**2+(u(nint(x/dx)+1)-utemp_J(nint(x/dx)+1))**2+sum
          !    print *, sum
          sum=(sum/(ispace+1))**(1.d0/2)  !get norm 2
          !    if ( norm2(utemp_J-u)<epsilon ) exit
          ! print *, sum,epsilon
          if ( sum<epsilon ) exit

       enddo
       ! print *, itime
    enddo

  end subroutine one_d_linear_adv_BTCS



  !============one_d_linear_advection:Lax Friedrichs:boundary=constant============


  subroutine one_d_linear_adv_Lax_Fr(c,dt,t,dx,x,u0,u)

    implicit none
    real(kind=8),intent(in)::c,dt,dx,t,x
    integer::itime,ispace
    real(kind=8),dimension(:),intent(in)::u0
    real(kind=8),dimension(size(u0)),intent(out)::u
    real(kind=8),dimension(size(u0))::ut_last
    u=u0
    print *, nint(t/dt)
    print *, nint(x/dx),c*dt/dx
    Do itime=1,nint(t/dt)+1
       ut_last=u
       Do ispace=2,nint(x/dx)

          u(ispace)=1.d0/2*(ut_last(ispace+1)+ut_last(ispace-1))-c*dt/dx/2*(ut_last(ispace+1)-ut_last(ispace-1))

       enddo
       u(1)=1.d0/2*(ut_last(2)+ut_last(nint(x/dx)))-c*dt/dx/2*(ut_last(2)-ut_last(nint(x/dx)))
       u(nint(x/dx)+1)=u(1)

    enddo

  end subroutine one_d_linear_adv_Lax_Fr

  !============one_d_linear_advection:Lax Wendroff:boundary=constant============


  subroutine one_d_linear_adv_Lax_Wen(c,dt,t,dx,x,u0,u)

    implicit none
    real(kind=8),intent(in)::c,dt,dx,t,x
    integer::itime,ispace
    real(kind=8),dimension(:),intent(in)::u0
    real(kind=8),dimension(size(u0)),intent(out)::u
    real(kind=8),dimension(size(u0))::ut_last
    u=u0
    print *, nint(t/dt)
    print *, nint(x/dx),c*dt/dx
    Do itime=1,nint(t/dt)+1
       ut_last=u
       Do ispace=2,nint(x/dx)

          u(ispace)=ut_last(ispace)-c*dt/dx/2*(ut_last(ispace+1)-ut_last(ispace-1)) &
               +(c*dt/dx)**2/2*(ut_last(ispace+1)-2*ut_last(ispace)+ut_last(ispace-1))

       enddo
       u(1)=ut_last(1)-c*dt/dx/2*(ut_last(2)-ut_last(nint(x/dx))) &
            +(c*dt/dx)**2/2*(ut_last(2)-2*ut_last(1)+ut_last(nint(x/dx)))
       u(nint(x/dx)+1)=u(1)

    enddo

  end subroutine one_d_linear_adv_Lax_Wen

  !============one_d_linear_advection Leapfrog FTBS:boundary=constant============


  subroutine one_d_linear_adv_Leapfrg_FTBS(c,dt,t,dx,x,u0,u)

    implicit none
    real(kind=8),intent(in)::c,dt,dx,t,x
    integer::itime,ispace
    real(kind=8),dimension(:),intent(in)::u0
    real(kind=8),dimension(size(u0)),intent(out)::u
    real(kind=8),dimension(size(u0))::ut_last_1,ut_last_2
    u=u0
    print *, nint(t/dt)
    print *, nint(x/dx),c*dt/dx

    ut_last_1=u
    Do ispace=2,nint(x/dx)

       u(ispace)=ut_last_1(ispace)-c*dt/dx*(ut_last_1(ispace)-ut_last_1(ispace-1))

     !  print u

    enddo
    u(1)=ut_last_1(1)-c*dt/dx*(ut_last_1(1)-ut_last_1(nint(x/dx)))
    u(nint(x/dx)+1)=u(1)
    ut_last_2=u

    Do itime=1,nint(t/dt)+1
       Do ispace=2,nint(x/dx)

          u(ispace)=ut_last_1(ispace)-c*dt/dx*(ut_last_2(ispace+1)-ut_last_2(ispace-1))
       enddo
       u(1)=ut_last_1(1)-c*dt/dx*(ut_last_2(2)-ut_last_2(nint(x/dx)))
       u(nint(x/dx)+1)=u(1)
       ut_last_1=ut_last_2
       ut_last_2=u

    enddo


  end subroutine one_d_linear_adv_Leapfrg_FTBS






end module scheme
