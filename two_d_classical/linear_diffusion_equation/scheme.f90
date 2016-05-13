! This file contains multiple schemes I need to use

module scheme

  implicit none
  real(kind=8)::pi
  save

contains



  !============two_d_diffusion:FTCS============


  subroutine two_d_diff_FTCS (nu,dt,t,dx,x,dy,y,u0,u)
    implicit none
    real(kind=8), intent(in) :: nu,dt,t,dx,x,dy,y
    integer::ispace,itime,jspace
    real(kind=8),dimension(:,:),intent(in)::u0
    real(kind=8),dimension(size(u0,1),size(u0,2)),intent(out)::u
    real(kind=8),dimension(size(u0,1),size(u0,2))::ut_last

    u=u0
    print *, nint(t/dt), nu,dt,t,dx,x,dy,y
    print *, nint(x/dx),nu,nu*dt/(dx**2),nu*dt/(dy**2)
    Do itime=1,nint(t/dt)+1
       ut_last=u
       Do jspace=2,nint(y/dy)
          Do ispace=2,nint(x/dx)

             u(jspace,ispace)=ut_last(jspace,ispace)&
                  +nu*dt/(dx**2)*(ut_last(jspace,ispace+1)-2*ut_last(jspace,ispace)+ut_last(jspace,ispace-1))&
                  +nu*dt/(dy**2)*(ut_last(jspace+1,ispace)-2*ut_last(jspace,ispace)+ut_last(jspace-1,ispace))

          enddo
       enddo
       !u(1)=ut_last(1)+nu*dt/(dx**2)*(ut_last(2)-2*ut_last(1)+ ut_last(nint(x/dx)))

       ! u(nint(x/dx)+1)=u(1)

    enddo


  end subroutine two_d_diff_FTCS


end module scheme
