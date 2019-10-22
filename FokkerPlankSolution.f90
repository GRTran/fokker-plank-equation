program FokkerPlankSolution
  use ODEsolver
  implicit none
  integer, parameter :: wp = selected_real_kind(8)
  integer, parameter :: neqn = 100
  real(wp), parameter :: pi = 3.1415926535846264338327950228

  character(len=15) :: num_nodes

  real(wp), allocatable :: z_coords(:)
  real(wp), allocatable :: temp(:)
  real(wp), allocatable :: outputs(:,:)
  real(wp) :: delta
  real(wp) :: region_start, region_end, radius, mass_calc, mass_act
  real(wp) :: start_time, end_time, time_step, time
  integer :: num_steps

  real(wp) :: initial_density, p0

  real(wp), allocatable, dimension( : )      ::    work
  INTEGER                                                   ::    iwork(5)
  integer :: iflag

  integer :: i

  region_start = 0._wp
  region_end = 50.0_wp
  radius = 0.4_wp

  start_time = 0._wp
  end_time = 1000._wp
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

  open(565,file='time_data.dat')
  time = start_time
  write(565,*) time
  do i = 2, num_steps + 1
    iflag = 1
    call ode(grad2, neqn, temp, time, time+time_step,1d-5,1d-5,iflag,work,iwork)
    if (iflag /= 2) stop 'Error Shamp Gordon Solver not solved'
    outputs(:,i) = temp
    write(565,*) time
    ! check the mass conservation
    mass_calc = sum(outputs(:,i)*(delta))

    write(*,*) time, mass_calc, mass_act!, (p / p0)
  enddo
  close(565)
  write(num_nodes,*) neqn
  open(145,file='density_data.dat')
  do i = 1, num_steps
    write(145,*) outputs(:,i) / initial_density
  enddo
  close(145)
  open(145,file='geometric_data.dat')
  write(145,*) z_coords
  close(145)

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

  ! subroutine analyticalExpression(
  !
  !   alpha = (2._wp*adv)/diff
  !   beta = (2._wp*lambda)/diff
  !   w = sqrt(beta-alpha**2/4)
  !   omega = w*(region_end-region_start)



end program
