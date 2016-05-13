! This gives the boundary conditions for all the algorithm

module boundary_conditions

  implicit none
  real(kind=8)::pi
  save

contains

  !     !============two_d_laplace_p0 ============
  subroutine two_d_laplace_p0(itime,ispace,jspace,dt,t,dx,x,dy,y,t_value,x_value,y_value,p0,debug)

    implicit none
    integer,intent(in)::itime,ispace,jspace
    integer::init_loop,jnit_loop,isetval
    real(kind=8),intent(in)::t,x,y,dt,dx,dy
    real(kind=8),intent(out),dimension(ispace+1)::x_value
    real(kind=8),intent(out),dimension(jspace+1)::y_value
    real(kind=8),intent(out),dimension(jspace+1,ispace+1)::p0
    real(kind=8),intent(out),dimension(itime+1)::t_value
    logical,intent(in)::debug

    character(40)::fmt,ch,filename

    do isetval=1,ispace+1
       x_value(isetval)=(isetval-1)*x/ispace
    end do

    do isetval=1,jspace+1
       y_value(isetval)=(isetval-1)*y/jspace
    end do

    do isetval=1,itime+1
       t_value(isetval)=(isetval-1)*t/itime
    end do

    do jnit_loop=1,jspace+1
       do init_loop=1,ispace+1

          p0(jnit_loop,init_loop)=0.d0

       end do
    enddo
    p0(:,1)=0  !
    p0(:,(nint(x/dx)+1))=y_value
    p0(1,:)=p0(2,:)
    p0((nint(y/dy)+1),:)=p0(nint(y/dy),:)


    !============test============





    open(unit=20,file='two_d_initial_p0.txt',status='replace')
    print *,size(p0),size(p0,1),size(p0,2),ispace+1
    write(ch,"(I10)") ispace+1
    fmt="("//trim(adjustl(ch))//"es24.14"//")" ! concatenate to make format statement
    print *, fmt, ch
    write(20,fmt) ((p0(isetval,:)),isetval=1,jspace+1)
    close(unit=20)

    open(unit=1002,file='two_d_initial_x_value.txt',status='replace')
    Do isetval=1,ispace+1
       write(1002, 302) x_value(isetval)
    enddo
302 format(es24.14)
    close(unit=1002)

    open(unit=1003,file='two_d_initial_y_value.txt',status='replace')
    Do isetval=1,jspace+1
       write(1003, 303) y_value(isetval)
    enddo
303 format(es24.14)
    close(unit=1003)


  end subroutine two_d_laplace_p0

end module boundary_conditions
