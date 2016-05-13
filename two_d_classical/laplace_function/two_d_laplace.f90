  !This is for multiple solution of one dimension of CFD



Program two_d_laplace
  use boundary_conditions
  use scheme

  !============parameter settings============


  implicit none
  integer::isetval,jsetval
  integer,parameter ::ispace=40,itime=40,jspace=40,iepsilon=50000
  real(kind=8),parameter::c=1,nu=0.1,rho=1.d0,sigma=0.0009,epsilon=1.0d-8
  real(kind=8)::dt,dx,dy,t,x,y
  real(kind=8),dimension(ispace+1)::x_value
  real(kind=8),dimension(jspace+1)::y_value
  real(kind=8),dimension(itime+1)::t_value
  real(kind=8),dimension(jspace+1,ispace+1)::u0,u,v0,v,p0,p,b
  logical::debug=.True., Need_p=.True.
  character(40)::fmt,ch,filename_p
  ! Need_p=.False.
  ! Needsecond=.False.
  debug=.False.

  x=2
  y=2
  dx=x/(ispace)
  dy=y/ispace
  dt=0.001
  t=itime*dt


  !========================

  !============boundary_conditions============


  !============two_d_laplace_p0 ============
  call two_d_laplace_p0(itime,ispace,jspace,dt,t,dx,x,dy,y,t_value,x_value,y_value,p0,debug)


  !============two_d_laplace============
  call two_d_laplace_b(dt,t,dx,x,x_value,dy,y,y_value,p0,b)
  call two_d_laplace_p(dx,x,x_value,dy,y,y_value,p0,p,b,epsilon,iepsilon);filename_p="two_d_laplace_p"


!!!!  ============write to file============

  if ( Need_p ) then

     !============write to second file============
     open(unit=2001,file=trim(filename_p)//'_results.txt',status='replace')
     write(ch,"(I10)") ispace+1
     fmt="("//trim(adjustl(ch))//"es24.14"//")" ! concatenate to make format statement
     print *, fmt, ch
     write(2001,fmt) ((p(isetval,:)),isetval=1,jspace+1)
     !     print "("//trim(adjustl(ch))//"es10.1"//")" , ((b(isetval,:)),isetval=1,jspace+1)

     close(unit=2001)
  end if




end program two_d_laplace
