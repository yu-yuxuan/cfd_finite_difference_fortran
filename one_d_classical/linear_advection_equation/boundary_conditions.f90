! This gives the boundary conditions for all the algorithm

module boundary_conditions
  implicit none
  real(kind=8)::pi
  save

contains

  !============triangle_bound============



  subroutine triangle_bound(ispace,x,x_value,u0,debug)

    implicit none
    integer,intent(in)::ispace
    integer::init_loop,isetval
    real(kind=8),intent(in)::x
    real(kind=8),intent(out),dimension(ispace+1)::x_value,u0
    logical,intent(in)::debug

    !============The first time to set initial: dull============


    ! Do init_loop=1,nint(ispace*0.9d0/2)
    !    u0(init_loop)=0
    ! enddo

    ! do init_loop=nint(ispace*0.9d0/2)+1,nint(ispace*1.d0/2)
    !    u0(init_loop)=10*(init_loop*2/ispace-0.9)
    ! end do

    ! do init_loop=nint(ispace*1.d0/2)+1,nint(ispace*1.1d0/2)
    !    u0(init_loop)=10*(1.1-init_loop*2/ispace)
    ! end do

    ! do init_loop=nint(ispace*1.1/2)+1,ispace+1
    !    u0(init_loop)=0
    ! end do

    !============the second============

    do isetval=1,ispace+1
       x_value(isetval)=(isetval-1)*x/ispace
    end do


    do init_loop=1,ispace+1

       if ( x_value(init_loop)>=0.9d0 .and. (x_value(init_loop)<=1.d0) ) then
          u0(init_loop)=10*(x_value(init_loop)-0.9)
       elseif ((x_value(init_loop)>=1.d0) .and. (x_value(init_loop)<=1.1d0)) then
          u0(init_loop)=10*(1.1-x_value(init_loop))
       else
          u0(init_loop)=0
       end if

    end do


    !============Test============


    if (debug) then
       Do isetval=1,ispace+1
          print 101,x_value(isetval), u0(isetval)
       enddo
101    format(2es24.14)
    endif


  end subroutine triangle_bound

  !     !============Dira delta function============


  !     call init_cond(0.d0,2.d0,ispace+1,0.d0,u0)
  !     u0(51)=1


  !     !============Test============


  !     if (debug) then
  !        print *,c,dt,t,dx,x,vis*dt/(dx**2)
  !        Do isetval=1,ispace+1
  !           print 101, isetval, u0(isetval)
  !        enddo
  ! 101    format(i3,es24.14)
  !     endif

  !     !============heaviside ============
  subroutine heaviside(itime,ispace,t,x,t_value,x_value,u0,debug)

    implicit none
    integer,intent(in)::itime,ispace
    integer::init_loop,isetval
    real(kind=8),intent(in)::t,x
    real(kind=8),intent(out),dimension(ispace+1)::x_value,u0
    real(kind=8),intent(out),dimension(itime+1)::t_value
    logical,intent(in)::debug

    do isetval=1,ispace+1
       x_value(isetval)=(isetval-1)*x/ispace
    end do

    do isetval=1,itime+1
       t_value(isetval)=(isetval-1)*t/itime
    end do


    do init_loop=1,ispace+1


       u0(init_loop)=0

    end do

    u0(1)=1
    !============Test============


    if (debug) then
       Do isetval=1,ispace+1
          print 101, x_value(isetval), u0(isetval)
       enddo
101    format(2es24.14)
    endif

  end subroutine heaviside

  !     !============square wave ============
  subroutine square_wave(itime,ispace,t,x,t_value,x_value,u0,debug)

    implicit none
    integer,intent(in)::itime,ispace
    integer::init_loop,isetval
    real(kind=8),intent(in)::t,x
    real(kind=8),intent(out),dimension(ispace+1)::x_value,u0
    real(kind=8),intent(out),dimension(itime+1)::t_value
    logical,intent(in)::debug

    do isetval=1,ispace+1
       x_value(isetval)=(isetval-1)*x/ispace
    end do

    do isetval=1,itime+1
       t_value(isetval)=(isetval-1)*t/itime
    end do


    do init_loop=1,ispace+1

       if ( x_value(init_loop)>=0.5d0 .and. (x_value(init_loop)<=1.d0) ) then
          u0(init_loop)=2.d0
       else
          u0(init_loop)=1.d0
       end if

    end do


    !============Test============


    if (debug) then
       Do isetval=1,ispace+1
          print 101, x_value(isetval), u0(isetval)
       enddo
101    format(2es24.14)
    endif


  end subroutine square_wave

end module boundary_conditions
