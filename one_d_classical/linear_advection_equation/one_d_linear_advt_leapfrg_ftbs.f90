  !This is for multiple solution of one dimension of CFD



Program one_d_linear_advt_leapfrg_ftbs
  use boundary_conditions
  use scheme

  !============parameter settings============


  implicit none
  integer::isetval
  integer,parameter ::ispace=100,itime=40
  real(kind=8),parameter::c=1,sigma=0.8
  real(kind=8)::dt,dx,t,x
  real(kind=8),dimension(ispace+1)::u0,x_value,u
  real(kind=8),dimension(itime+1)::t_value
  character(40)::ch, fmt, filename
  logical::debug=.False.
  ! logical::debug=.True.
  !========================

  !============set some global parameters============
  call initialize()
  x=4
  dx=x/(ispace)
  dt=dx*sigma/c
  t=dt*itime
  print 1010,t

1010 format(' t = ', es24.14)

  !============boundary_conditions============

  !     !============heaviside ============
  ! call  heaviside(itime,ispace,t,x,t_value,x_value,u0,debug)
  call  square_wave(itime,ispace,t,x,t_value,x_value,u0,debug)
  ! call triangle_bound(ispace,x,x_value,u0,debug)

  !============one_d_linear_advection: ============

  call one_d_linear_adv_Leapfrg_FTBS(c,dt,t,dx,x,u0,u);filename='one_d_linear_advt_leapfrg_ftbs'



  !  ============test============


  if (debug) then
     Do isetval=1,ispace+1
        print 201, x_value(isetval), u(isetval)
     enddo
201  format(2es24.10)
  endif


  !============write to file============


  open(unit=1001,file=trim(filename)//'_results.txt',status='replace')
  Do isetval=1,ispace+1

     write(1001, 301) x_value(isetval),u(isetval)
  enddo
301 format(2es24.14)
  close(unit=1001)


end program one_d_linear_advt_leapfrg_ftbs
