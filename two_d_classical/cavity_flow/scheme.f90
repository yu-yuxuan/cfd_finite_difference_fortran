! This file contains multiple schemes I need to use

module scheme

  implicit none
  real(kind=8)::pi
  save

contains


 !============two_d_laplace_b============


  subroutine two_d_laplace_b(dt,t,dx,x,x_value,dy,y,y_value,p0,b,u,v,rho)


    implicit none
    real(kind=8),intent(in)::dx,x,dy,y,dt,t,rho
    real(kind=8)::sum
    integer::itime,ispace,jspace
    real(kind=8),intent(in),dimension(:)::x_value
    real(kind=8),intent(in),dimension(:)::y_value
    real(kind=8),dimension(:,:),intent(in)::p0
    real(kind=8),dimension(size(p0,1),size(p0,2)),intent(out)::b,u,v

    do jspace=1,nint(y/dy)+1
       do ispace=1,nint(x/dx)+1  !get the Xn+1 by Xn
          b(jspace,ispace)=0.d0
       enddo
    enddo

    !    b(nint(y/dy/4),nint(x/dx/4))=100
    !   b(nint(3.d0*y/dy/4),nint(3.d0*x/dx/4))=-100


    do jspace=2,nint(y/dy)
       do ispace=2,nint(x/dx)  !get the Xn+1 by Xn
          b(jspace,ispace)&
               =rho*(1.d0/dt*((u(jspace,ispace+1)-u(jspace,ispace-1))/2.d0/dx+(v(jspace+1,ispace)-v(jspace-1,ispace))/2.d0/dy) &
             -((u(jspace,ispace+1)-u(jspace,ispace-1))/2.d0/dx)**2.d0 &
               -2.d0*(u(jspace+1,ispace)-u(jspace-1,ispace))/2.d0/dy*(v(jspace,ispace+1)-v(jspace,ispace-1))/2.d0/dx &
               -((v(jspace+1,ispace)-v(jspace-1,ispace))/2.d0/dy)**2.d0 )
       enddo
    enddo

!============peridoic condition============


    do jspace=2,nint(y/dy)
       b(jspace,1)&
            =rho*(1.d0/dt*((u(jspace,2)-u(jspace,nint(x/dx)))/2.d0/dx+(v(jspace+1,1)-v(jspace-1,1))/2.d0/dy) &
            -((u(jspace,2)-u(jspace,nint(x/dx)))/2.d0/dx)**2.d0 &
            -2.d0*(u(jspace+1,1)-u(jspace-1,1))/2.d0/dy*(v(jspace,2)-v(jspace,nint(x/dx)))/2.d0/dx &
            -((v(jspace+1,1)-v(jspace-1,1))/2.d0/dy)**2.d0 )
    enddo

    b(:,nint(x/dx)+1)=b(:,1)

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


       p(:,1)=p(:,2)  !
       p(:,(nint(x/dx)+1))=p(:,(nint(x/dx)))
       p(1,:)=p(2,:)
       p((nint(y/dy)+1),:)=0.d0

       ! p(:,1)=0  !
       ! p(:,(nint(x/dx)+1))=y_value
       ! p(1,:)=p(2,:)
       ! p((nint(y/dy)+1),:)=p(nint(y/dy),:)
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

  !============two dimension cavity function============

  subroutine two_d_uv_cavity(nu,rho,dt,t,dx,x,dy,y,u0,u,v0,v,p)


    implicit none
    real(kind=8),intent(in)::nu,rho,dt,dx,t,x,dy,y
    integer::itime,ispace,jspace
    real(kind=8),dimension(:,:),intent(in)::u0,v0,p
    real(kind=8),dimension(size(u0,1),size(u0,2)),intent(out)::u,v
    real(kind=8),dimension(size(u0,1),size(u0,2))::ut_last,vt_last
    ut_last=u
    vt_last=v
    Do jspace=2,nint(y/dy)

       Do ispace=2,nint(x/dx)

          u(jspace, ispace)=ut_last(jspace,ispace)&
               -ut_last(jspace,ispace)*dt/dx*(ut_last(jspace,ispace)-ut_last(jspace,ispace-1)) &
               -vt_last(jspace,ispace)*dt/dy*(ut_last(jspace,ispace)-ut_last(jspace-1,ispace)) &
               -dt/rho*(p(jspace,ispace+1)-p(jspace,ispace-1))/2.d0/dx &
               +nu*dt/(dx**2.d0)*(ut_last(jspace,ispace+1)-2.d0*ut_last(jspace,ispace)+ut_last(jspace,ispace-1))&
               +nu*dt/(dy**2.d0)*(ut_last(jspace+1,ispace)-2.d0*ut_last(jspace,ispace)+ut_last(jspace-1,ispace))
                  !This is a source added, it will be delleted depending what condtions you are


       enddo
    enddo
    Do jspace=2,nint(y/dy)
       Do ispace=2,nint(x/dx)
          v(jspace, ispace)=vt_last(jspace,ispace)&
               -ut_last(jspace,ispace)*dt/dx*(vt_last(jspace,ispace)-vt_last(jspace,ispace-1)) &
               -vt_last(jspace,ispace)*dt/dy*(vt_last(jspace,ispace)-vt_last(jspace-1,ispace)) &
               -dt/rho*(p(jspace+1,ispace)-p(jspace-1,ispace))/2.d0/dy &
               +nu*dt/(dx**2.d0)*(vt_last(jspace,ispace+1)-2.d0*vt_last(jspace,ispace)+vt_last(jspace,ispace-1))&
               +nu*dt/(dy**2.d0)*(vt_last(jspace+1,ispace)-2.d0*vt_last(jspace,ispace)+vt_last(jspace-1,ispace))


       enddo
    enddo

    Do ispace=1,nint(x/dx)+1
       u(nint(y/dy)+1,ispace)=1.0d0
    enddo

    !============test============


    ! u(1,1)=4
    ! u(3,2)=7
    ! print *, u(:,1)-2*u(:,2)
  end subroutine two_d_uv_cavity



end module scheme
