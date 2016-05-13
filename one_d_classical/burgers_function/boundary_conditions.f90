! This gives the boundary conditions for all the algorithm

module boundary_conditions
  implicit none
  real(kind=8)::pi
  save

contains

 !     !============Burgers_initial_conditions ============
  subroutine burgers_ini_cond(itime,ispace,t,x,t_value,x_value,u0,vis,debug)

    implicit none
    integer,intent(in)::itime,ispace
    integer::init_loop,isetval
    real(kind=8),intent(in)::t,x,vis
    real(kind=8),intent(out),dimension(ispace+1)::x_value,u0
    real(kind=8),dimension(ispace+1)::phi,dphi
    real(kind=8),intent(out),dimension(itime+1)::t_value
    logical,intent(in)::debug

    do isetval=1,ispace+1
       x_value(isetval)=(isetval-1)*x/ispace
    end do

    do isetval=1,itime+1
       t_value(isetval)=(isetval-1)*t/itime
    end do

    print *, pi
    do init_loop=1,ispace+1


       phi(init_loop)=exp(-x_value(init_loop)**2/4/vis)+exp(-(x_value(init_loop)-2*pi)**2/4/vis)
       dphi(init_loop)=exp(-x_value(init_loop)**2/4/vis)*(-2*x_value(init_loop)/4/vis)+&
            exp(-(x_value(init_loop)-2*pi)**2/4/vis)*(-2*(x_value(init_loop)-2*pi)/4/vis)

       u0(init_loop)=-2*vis/phi(init_loop)*dphi(init_loop)+4

    end do



    !============Test============

    open(unit=1002,file='burgers_initial_results.txt',status='replace')
    Do isetval=1,ispace+1

       write(1002, 302) x_value(isetval),u0(isetval)
    enddo
302 format(2es24.14)
    close(unit=1002)


  end subroutine burgers_ini_cond


end module boundary_conditions
