program FokkerPlankSolution
  use ODEsolver
  implicit none
  integer, parameter :: wp = selected_real_kind(8)
  integer, parameter :: neqn = 100
  real(wp), parameter :: pi = 3.1415926535846264338327950228

  character(len=15) :: num_nodes

  real(wp), allocatable :: z_coords(:)
  real(wp), allocatable :: p(:)
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
  region_end = 0.5_wp
  radius = 0.4_wp

  start_time = 0._wp
  end_time = 500._wp
  num_steps = 1000

  allocate( work( 100+21*neqn ) )

  allocate(z_coords(neqn))
  allocate(p(neqn))
  allocate(outputs(neqn,num_steps + 1))

  delta = (region_end - region_start) / neqn
  z_coords(1) = region_start + delta / 2._wp

  do i = 2, neqn
    z_coords(i) = z_coords(i-1) + delta
  enddo

  time_step = ( end_time - start_time ) / num_steps

  p = 1._wp/(region_end-region_start)
  p0 = (1._wp/(region_end-region_start))*sum(p)*delta
  initial_density = 2400

  mass_act = initial_density * (region_end-region_start)

  outputs(:,1) = initial_density * (p / p0)

  open(565,file='time_data.dat')
  time = start_time
  write(565,*) time
  do i = 2, num_steps + 1
    iflag = 1
    call ode(grad, neqn, p, time, time+time_step,1d-5,1d-5,iflag,work,iwork)
    if (iflag /= 2) stop 'Error Shamp Gordon Solver not solved'
    p0 = (1._wp/(region_end-region_start))*sum(p)*delta
    outputs(:,i) = initial_density * (p / p0)
    write(565,*) time

    ! check the mass conservation
    mass_calc = (1._wp/(region_end-region_start))*sum(outputs(:,i)*(delta)) * (region_end-region_start)

    write(*,*) time, mass_calc, mass_act!, (p / p0)
  enddo
  close(565)
  write(num_nodes,*) neqn
  open(145,file='density_data_'//trim(num_nodes)//'.dat')
  do i = 1, num_steps
    write(145,*) outputs(:,i) / initial_density
  enddo
  close(145)
  open(145,file='geometric_data_'//trim(num_nodes)//'.dat')
  write(145,*) z_coords/(region_end-region_start)
  close(145)

contains

  subroutine grad(t,y,yp)
    implicit none
    real(wp), intent(in) :: t
    real(wp), intent(inout) :: y(neqn)
    real(wp) :: yp(neqn)
    real(wp) :: mu, sigma
    integer :: i
    real(wp) :: region_start, region_end, delta
    real(wp) :: bottom_ghost_node, top_ghost_node

    region_start = 0._wp
    region_end = 0.5_wp
    delta = (region_end - region_start) / neqn

    mu = 0.001_wp
    sigma = 0.01_wp

    ! ghost nodes
    top_ghost_node = ( ((sigma**2/(2._wp*delta)) - (mu/2._wp) ) / &
                     ((sigma**2/(2._wp*delta)) + (mu/2._wp) ) ) * y(neqn)
    ! bottom_ghost_node = - y(1) - (delta*sigma**2/(2._wp*mu))*(top_ghost_node - y(neqn))
    bottom_ghost_node = ( ( -(mu/2._wp) + (sigma**2/(2._wp*delta)) ) / &
                        ((sigma**2/(2._wp*delta)) + (mu/2._wp) ) ) *y(1)

    ! bottom node
    yp(1) = -(mu/delta)*(y(1)-y(2)) + (sigma**2/(2._wp*delta**2))*(bottom_ghost_node-2*y(1)+y(2))
    ! intermediate spatial neqn
    do i = 2, size(yp) - 1
        yp(i) = -(mu/delta)*(y(i)-y(i+1)) + (sigma**2/(2._wp*delta**2))*(y(i-1)-2*y(i)+y(i+1))
    enddo
    ! top node
    yp(neqn) = -(mu/delta)*(y(neqn)-top_ghost_node) + (sigma**2/(2._wp*delta**2))*(y(neqn-1)-2*y(neqn)+top_ghost_node)

  end subroutine

  ! subroutine analyticalExpression(
  !
  !   alpha = (2._wp*mu)/sigma**2
  !   beta = (2._wp*lambda)/sigma**2
  !   w = sqrt(beta-alpha**2/4)
  !   omega = w*(region_end-region_start)



end program
