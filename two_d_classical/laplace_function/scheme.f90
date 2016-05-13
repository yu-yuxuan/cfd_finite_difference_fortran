! This file contains multiple schemes I need to use

module scheme

  implicit none
  real(kind=8)::pi
  save

contains


   !============two_d_laplace_b============


  subroutine two_d_laplace_b(dt,t,dx,x,x_value,dy,y,y_value,p0,b)


    implicit none
    real(kind=8),intent(in)::dx,x,dy,y,dt,t
    real(kind=8)::sum
    integer::itime,ispace,jspace
    real(kind=8),intent(in),dimension(:)::x_value
    real(kind=8),intent(in),dimension(:)::y_value
    real(kind=8),dimension(:,:),intent(in)::p0
    real(kind=8),dimension(size(p0,1),size(p0,2)),intent(out)::b

    do jspace=1,nint(y/dy)+1
       do ispace=1,nint(x/dx)+1  !get the Xn+1 by Xn
          b(jspace,ispace)=0.d0
       enddo
    enddo

  end subroutine two_d_laplace_b


  !============two_d_laplace_p============


  subroutine two_d_laplace_p(dx,x,x_value,dy,y,y_value,p0,p,b,epsilon,iepsilon)


    implicit none
    real(kind=8),intent(in)::dx,x,dy,y,epsilon
    real(kind=8)::sum
    integer::itime,ispace,jspace
    integer,intent(in)::iepsilon
    real(kind=8),intent(in),dimension(:)::x_value
    real(kind=8),intent(in),dimension(:)::y_value
    real(kind=8),dimension(:,:),intent(in)::p0,b
    real(kind=8),dimension(size(p0,1),size(p0,2)),intent(inout)::p
    real(kind=8),dimension(size(p0,1),size(p0,2))::ptemp_J

    !============Jacobi iteration============
    !   p=p0
    itime=0
    do
       itime=itime+1
       sum=0
       ptemp_J=p
       do jspace=2,nint(y/dy)
          do ispace=2,nint(x/dx)  !get the Xn+1 by Xn

             p(jspace,ispace)=((ptemp_J(jspace,ispace+1)+ptemp_J(jspace,ispace-1))/(dx**2.d0)&
                  +(ptemp_J(jspace+1,ispace)+ptemp_J(jspace-1,ispace))/(dy**2.d0) -b(jspace,ispace))&
                  /(2.d0/(dx**2.d0)+2.d0/(dy**2.d0))

             !    print * ,"used",u(ispace),utemp_J(ispace)
          enddo
       enddo


       ! p(:,1)=p(:,2)  !
       ! p(:,(nint(x/dx)+1))=p(:,(nint(x/dx)))
       ! p(1,:)=p(2,:)
       ! p((nint(y/dy)+1),:)=0.d0

       p(:,1)=0  !
       p(:,(nint(x/dx)+1))=y_value
       p(1,:)=p(2,:)
       p((nint(y/dy)+1),:)=p(nint(y/dy),:)
       !    print *, sum
       !============error============
       do jspace=1,nint(y/dy)+1
          do ispace=1,nint(x/dx)+1  !get the Xn+1 by Xn


             sum=(p(jspace,ispace)-ptemp_J(jspace,ispace))**2+sum
          enddo
       enddo
       sum=(sum/(ispace+1)/(jspace+1))**(1.d0/2)  !get norm 2
       !    if ( norm2(utemp_J-u)<epsilon ) exit
       ! print *, sum,epsilon
       if ( sum<epsilon ) exit
       if (itime==iepsilon) exit
    enddo
    !      u(1)=ut_last(1)-c*dt/dx*(ut_last(1)-ut_last(nint(x/dx)))
    print *, itime


  end subroutine two_d_laplace_p

end module scheme
