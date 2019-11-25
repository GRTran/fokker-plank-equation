program FokkerPlankSolution
  use ODEsolver
  use FiniteVolume1D
  implicit none

  external dlsoda

  integer, parameter :: neqn = 100

  character(len=15) :: num_nodes

  real(wp), allocatable :: z_coords(:)
  real(wp), allocatable :: temp(:)
  real(wp), allocatable :: outputs(:,:)
  real(wp) :: delta
  real(wp) :: region_start, region_end, radius, mass_calc, mass_act
  real(wp) :: start_time, end_time, time_step, time, jdum
  integer :: num_steps
  INTEGER :: itol, itask, istate, iopt, liw, lrw, jt
  real(wp) :: atol, rtol

  real(wp) :: initial_density, p0

  real(wp), allocatable, dimension( : )      ::    work
  DOUBLE PRECISION, ALLOCATABLE :: rwork(:)
    INTEGER, ALLOCATABLE :: i2work(:)
  INTEGER                                                   ::    iwork(5)
  integer :: iflag

  integer :: i, status


  region_start = 0._wp
  region_end = 50.0_wp
  radius = 0.4_wp

  start_time = 0._wp
  end_time = 100._wp
  num_steps = 1000

  allocate( work( 100+21*neqn ) )

  allocate(z_coords(neqn))
  allocate(temp(neqn))
  allocate(outputs(neqn,num_steps + 1))

  delta = (region_end - region_start) / neqn
  z_coords(1) = region_start + delta / 2._wp

  do i = 2, neqn
    z_coords(i) = z_coords(i-1) + delta
  enddo

  time_step = ( end_time - start_time ) / num_steps

  initial_density = 3.

  mass_act = initial_density * (region_end-region_start)

  outputs(:,1) = initial_density
  temp = initial_density

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! SHAMPINE GORDON SOLVE
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! open(565,file='time_data.dat')
  ! time = start_time
  ! write(565,*) time
  ! do i = 2, num_steps + 1
  !   iflag = 1
  !   call ode(grad4, neqn, temp, time, time+time_step,1d-5,1d-5,iflag,work,iwork)
  !   if (iflag /= 2) stop 'Error Shamp Gordon Solver not solved'
  !   outputs(:,i) = temp
  !   write(565,*) time
  !   ! check the mass conservation
  !   mass_calc = sum(outputs(:,i)*(delta))
  !
  !   write(*,*) time, mass_calc, mass_act!, (p / p0)
  ! enddo
  ! close(565)
  ! write(num_nodes,*) neqn
  ! open(145,file='density_data.dat')
  ! do i = 1, num_steps
  !   write(145,*) outputs(:,i)
  ! enddo
  ! close(145)
  ! open(145,file='geometric_data.dat')
  ! write(145,*) z_coords
  ! close(145)
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! END OF SHAMPINE GORDON SOLVE
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! ODEPACK SOLVE
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  outputs(:,1) = initial_density
  temp = initial_density

  ! set control parameters
  itol = 1
  rtol = 1e-5
  atol = 1e-5
  itask = 1
  istate = 1
  iopt = 0
  lrw = 22 + neqn * MAX(16, neqn + 9) + neqn**2
  ALLOCATE( rwork(  lrw ) )
  liw = 20 + neqn
  ALLOCATE( i2work(  liw ) )
  rwork = 0._wp; iwork = 0
  jdum = 2
  jt = 2

  open(565,file='time_data_odepack.dat')
  time = start_time
  write(565,*) time
  do i = 2, num_steps + 1
    ! istate=1;! iopt = 0;! jdum = 2; jt = 2

    iopt = 0
    istate = 1
    CALL dlsoda ( grad5, neqn, temp, time, time+time_step, itol, rtol, atol, itask, istate, iopt, rwork,&
                      lrw, iwork, liw, JEX, jt )
    write(*,*) istate
    ! if (istate < 0) then
    !   write(*,*) istate
    !   stop
    ! endif

    outputs(:,i) = temp
    write(565,*) time
    ! check the mass conservation
    mass_calc = sum(outputs(:,i)*(delta))

    write(*,*) i, time, mass_calc, mass_act!, (p / p0)
  enddo
  write(*,*) 'DONE'
  close(565)
  ! write(num_nodes,*) neqn
  open(145,file='density_data_odepack.dat')
  do i = 1, num_steps
    write(145,*) outputs(:,i)
  enddo
  close(145)
  open(145,file='geometric_data_odepack.dat')
  write(145,*) z_coords
  close(145)
  WRITE(*,*) size(i2work), i2work(1)

  write(*,*) allocated(i2work), allocated(rwork)
  ! deallocate(i2work, rwork, stat = status)
  ! if(status/=0) write(*,*) 'Error deallocating work'
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! END OF ODEPACK CALL
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

contains

  ! develop a static fokker-planck solution utilising a change of variable to represent the system.
  ! in this case there is no moving boundary, i.e. s(t) = 0.5
  subroutine grad2(t,y,yp)
    implicit none
    real(wp), intent(in) :: t
    real(wp), intent(inout) :: y(neqn)
    real(wp) :: yp(neqn)
    real(wp) :: adv, diff, st
    integer :: i
    real(wp) :: region_start, region_end, delta
    real(wp) :: bottom_ghost_node, top_ghost_node

    region_start = 0._wp
    region_end = 1.0_wp
    delta = (region_end - region_start) / neqn

    adv = 0.01_wp
    diff = 1._wp
    diff = diff
    st = 50._wp

    ! ghost nodes
    top_ghost_node = ( ((diff/(st*delta)) - (adv/2._wp) ) / &
                     ((diff/(st*delta)) + (adv/2._wp) ) ) * y(neqn)
    ! bottom_ghost_node = - y(1) - (delta*diff/(2._wp*adv))*(top_ghost_node - y(neqn))
    bottom_ghost_node = ( ( (diff/(st*delta)) + (adv/2._wp) ) / &
                        ((diff/(st*delta)) - (adv/2._wp) ) ) *y(1)

    ! bottom node
    yp(1) = -(adv/(st*delta))*(y(1)-y(2)) + (diff/(st**2*delta**2))*(bottom_ghost_node-2._wp*y(1)+y(2))
    ! intermediate spatial neqn
    do i = 2, size(yp) - 1
        yp(i) = -(adv/(st*delta))*(y(i)-y(i+1)) + (diff/(st**2*delta**2))*(y(i-1)-2._wp*y(i)+y(i+1))
    enddo
    ! top node
    yp(neqn) = -(adv/(st*delta))*(y(neqn)-top_ghost_node) + (diff/(st**2*delta**2))*(y(neqn-1)-2._wp*y(neqn)+top_ghost_node)

  end subroutine

  subroutine grad(t,y,yp)
    implicit none
    real(wp), intent(in) :: t
    real(wp), intent(inout) :: y(neqn)
    real(wp) :: yp(neqn)
    real(wp) :: adv, diff
    integer :: i
    real(wp) :: region_start, region_end, delta
    real(wp) :: bottom_ghost_node, top_ghost_node

    region_start = 0._wp
    region_end = 50.0_wp
    delta = (region_end - region_start) / neqn

    adv = 0.01_wp
    diff = 1._wp

    ! ghost nodes
    top_ghost_node = ( ((diff/(delta)) - (adv/2._wp) ) / &
                     ((diff/(delta)) + (adv/2._wp) ) ) * y(neqn)
    ! bottom_ghost_node = - y(1) - (delta*diff/(2._wp*adv))*(top_ghost_node - y(neqn))
    bottom_ghost_node = ( ( (adv/2._wp) + (diff/(delta)) ) / &
                        ((diff/(delta)) - (adv/2._wp) ) ) *y(1)

    ! bottom node
    yp(1) = -(adv/delta)*(y(1)-y(2)) + (diff/(delta**2))*(bottom_ghost_node-2._wp*y(1)+y(2))
    ! intermediate spatial neqn
    do i = 2, size(yp) - 1
        yp(i) = -(adv/delta)*(y(i)-y(i+1)) + (diff/(delta**2))*(y(i-1)-2._wp*y(i)+y(i+1))
    enddo
    ! top node
    yp(neqn) = -(adv/delta)*(y(neqn)-top_ghost_node) + (diff/(delta**2))*(y(neqn-1)-2._wp*y(neqn)+top_ghost_node)

    ! write(*,*) yp
  end subroutine

  subroutine grad3(t,y,yp)
    implicit none
    real(wp), intent(in) :: t
    real(wp), intent(inout) :: y(neqn)
    real(wp) :: yp(neqn)
    real(wp) :: adv, diff, st
    integer :: i
    real(wp) :: region_start, region_end, delta

    region_start = 0._wp
    region_end = 1.0_wp
    delta = (region_end - region_start) / neqn

    adv = 0.01_wp
    diff = 0.1_wp
    diff = diff
    st = 50._wp


    !! Mass conserved derivation
    yp(1) = -(adv/(st*delta))*(-y(2)) + (diff/(st**2*delta**2))*(-y(1)+y(2))

    ! write(*,*) 1, advection_a, advection_b, diffusion

    do i = 2, neqn - 1
      yp(i) = -(adv/(st*delta))*(y(i)-y(i+1)) + (diff/(st**2*delta**2))*(y(i-1)-2._wp*y(i)+y(i+1))
      write(*,*) -(adv/(st*delta)) -2._Wp*(diff/(st**2*delta**2))
      ! write(*,*) i, advection_a, advection_b, diffusion
    enddo

    !! Mass conserved derivation
    yp(neqn) = -(adv/(st*delta))*(y(neqn)) + (diff/(st**2*delta**2))*(y(neqn-1)-y(neqn))

  end subroutine

  subroutine grad4(t,y,yp)
    implicit none
    type(fv1DTemporalSpatialClass) :: fv
    real(wp), intent(in) :: t
    real(wp), intent(inout) :: y(neqn)
    real(wp) :: yp(neqn)
    real(wp) :: adv, diff, st
    integer :: i
    real(wp) :: region_start, region_end, delta
    real(wp) :: boundary
    real(wp), dimension(neqn) :: diff_arr, adv_arr, delta_arr, vol_arr

    region_start = 0._wp
    region_end = 1.0_wp
    delta = (region_end - region_start) / neqn

    adv = 0.01_wp
    diff = 0.1_wp
    diff = diff
    st = 50._wp
    boundary = (adv*st)/diff

    diff_arr = diff
    adv_arr = adv
    delta_arr = delta
    vol_arr = 0._wp

    ! call fv%initialise_mass_matrix(diff_arr/st**2, adv_arr/st, vol_arr, delta_arr, neqn, .true., (/boundary,boundary/))
    call fv%initialise_mass_matrix(diff/st**2, adv_arr/st, 0._wp, delta, neqn, .false., (/boundary,boundary/))
    ! call fv%initialise_mass_matrix(diff/st**2, adv/st, 0._wp, delta, neqn, .true., (/boundary,boundary/))

    ! call fv%print_mass_matrix()
    call fv%solve_for_gradient(y,yp)
    ! write(*,*) y
    call fv%destroy()
  end subroutine

  subroutine grad5(neq, t,y,yp)
    implicit none
    type(fv1DTemporalSpatialClass) :: fv
    real(wp), intent(in) :: t
    real(wp), intent(inout) :: y(neqn)
    real(wp) :: yp(neqn)
    real(wp) :: adv, diff, st
    integer :: i, neq
    real(wp) :: region_start, region_end, delta
    real(wp) :: boundary
    real(wp), dimension(neqn) :: diff_arr, adv_arr, delta_arr, vol_arr

    region_start = 0._wp
    region_end = 1.0_wp
    delta = (region_end - region_start) / neqn

    adv = 0.01_wp
    diff = 0.1_wp
    diff = diff
    st = 50._wp
    boundary = (adv*st)/diff

    diff_arr = diff
    adv_arr = adv
    delta_arr = delta
    vol_arr = 0._wp

    ! call fv%initialise_mass_matrix(diff_arr/st**2, adv_arr/st, vol_arr, delta_arr, neqn, .true., (/boundary,boundary/))
    call fv%initialise_mass_matrix(diff/st**2, adv_arr/st, 0._wp, delta, neqn, .false., (/boundary,boundary/))
    ! call fv%initialise_mass_matrix(diff/st**2, adv/st, 0._wp, delta, neqn, .true., (/boundary,boundary/))

    ! call fv%print_mass_matrix()
    call fv%solve_for_gradient(y,yp)
    ! write(*,*) y
    call fv%destroy()
  end subroutine
  ! subroutine analyticalExpression(
  !
  !   alpha = (2._wp*adv)/diff
  !   beta = (2._wp*lambda)/diff
  !   w = sqrt(beta-alpha**2/4)
  !   omega = w*(region_end-region_start)

  SUBROUTINE JEX(NEQ,T,Y,ML,MU,PD,NRPD)
  IMPLICIT NONE
  INTEGER, INTENT (IN) :: NEQ, ML, MU, NRPD
  DOUBLE PRECISION, INTENT (IN) :: T
  DOUBLE PRECISION, INTENT (IN) :: Y(NEQ)
  DOUBLE PRECISION, INTENT (OUT) :: PD(NRPD,NEQ)
  PD(1,1) = -.04D0
  PD(1,2) = 1.D4*Y(3)
  PD(1,3) = 1.D4*Y(2)
  PD(2,1) = .04D0
  PD(2,3) = -PD(1,3)
  PD(3,2) = 6.E7*Y(2)
  PD(2,2) = -PD(1,2) - PD(3,2)
  RETURN
END SUBROUTINE JEX



end program
