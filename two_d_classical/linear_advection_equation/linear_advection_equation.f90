  !This is for multiple solution of one dimension of CFD



Program two_d_linear_advection
  use boundary_conditions
  use scheme

  !============parameter settings============


  implicit none
  integer::isetval,jsetval
  integer,parameter ::ispace=40,itime=700,jspace=40
  real(kind=8),parameter::c=1
  real(kind=8)::dt,dx,dy,t,x,y
  real(kind=8),dimension(ispace+1)::x_value
  real(kind=8),dimension(jspace+1)::y_value
  real(kind=8),dimension(itime+1)::t_value
  real(kind=8),dimension(jspace+1,ispace+1)::u0,u
  logical::debug=.True.
  character(40)::fmt,ch,filename
  debug=.False.

  x=2
  y=2
  dx=x/(ispace)
  dy=y/(jspace)
  dt=0.001
  t=itime*dt



  !========================

  !============boundary_conditions============
  !     !============two_d_square_wave ============
  call two_d_square_wave(itime,ispace,jspace,t,x,y,t_value,x_value,y_value,u0,debug)
  !============two_d_linear_convection: FTBS============


  call two_d_linear_adv_FTBS(c,dt,t,dx,x,dy,y,u0,u); filename="two_d_linear_advection_FTBS"




!!!!  ============write to file============
  open(unit=2001,file=trim(filename)//'_results.txt',status='replace')
  write(ch,"(I10)") ispace+1
  fmt="("//trim(adjustl(ch))//"es24.14"//")" ! concatenate to make format statement
  print *, fmt, ch
  write(2001,fmt) ((u(isetval,:)),isetval=1,jspace+1)

  close(unit=2001)



end program two_d_linear_advection
