module FiniteVolume1D
  implicit none

  type fv1DTemporalSpatialClass
    real(wp), allocatable                   ::    lower(:)
    real(wp), allocatable                   ::    centre(:)
    real(wp), allocatable                   ::    upper(:)
    real(wp), allocatable                   ::    source(:)
    real(wp), allocatable                   ::    cell_spacings(:)
    integer :: nodes
  contains
    generic, public :: initialise_mass_matrix => initialise_1d_mass_matrix_data, initialise_1d_mass_matrix_const
    procedure, public, pass :: initialise_1d_mass_matrix_data
    procedure, public, pass :: initialise_1d_mass_matrix_const
    procedure, public, pass :: initialise_source_vector
    procedure, public, pass :: solve_for_gradient
  end type

contains

  subroutine initialise_1d_mass_matrix_data( this, diffusive_vector, advective_vector, volume_vector, cell_spacings, nodes, upwind_direction )
    class(fv1DTemporalSpatialClass), intent(out) :: this
    real(wp), intent(in) :: diffusive_vector(:) ! the diffusive vector is defined within cells
    real(wp), intent(in) :: advective_vector(:) ! the advective vector is defined on the node boundaries, excludes one of the boundaries
    real(wp), intent(in) :: volume_vector(:) ! this represents a volume amount within the system
    real(wp), intent(in) :: cell_spacings(:)
    integer, intent(in) :: nodes
    logical, intent(in) :: upwind_direction ! 1 if flow is going in direction of increasing node size, 0 otherwise

    real(wp) :: matrix_coefficients(3)
    integer :: status

    allocate( this%cell_spacings(nodes), this%lower(nodes), this%centre(nodes), this%upper(nodes), this%source(nodes), stat=status )
    if( status /= 0 ) stop 'Error allocating cell spacings in 1D finite volume matrix generation'

    this%cell_spacings = cell_spacings
    this%lower = 0._wp; this%upper = 0._wp; this%centre = 0._wp ; this%source = 0._wp
    this%nodes = nodes

    associate( d => diffusive_vector, v => advective_vector, delta => cell_spacings )

      ! lower boundary condition
      this%centre(1) = (d(1)*gamma1(1) + v(1) )/ delta(1) + volume_vector(1)
      this%upper(1) = -2._wp / (delta(2)/d(2) + delta(1)/d(1))
      if( upwind_direction .eqv. .true. ) this%upper(1) = this%upper(1) - v(1)
      this%upper(1) = this%upper(1) / delta(1)

      do i = 2, nodes - 1

        ! diffusive terms

        ! node in the direction of decreasing spatial indexing
        this%lower(i) = -2._wp / (delta(i-1)/d(i-1) + delta(i)/d(i))
        ! node in the direction of increasing spatial indexing
        this%upper(i) = -2._wp / (delta(i+1)/d(i+1) + delta(i)/d(i))
        ! centre of the cell
        this%centre(i) = volume_vector*delta(i) - (this%lower(i)+this%upper(i))

        ! advective terms
        if( upwind_direction .eqv. .true. ) then
          this%lower(i) = this%lower(i) - v(i-1)
          this%centre(i) = this%centre(i) + v(i)
        else
          this%upper(i) = this%upper(i) - v(i)
          this%centre(i) = this%centre(i) + v(i-1)
        endif

        ! divide by spacing
        this%lower(i) = this%lower(i) / delta(i)
        this%centre(i) = this%centre(i) / delta(i)
        this%upper(i) = this%upper(i) / delta(i)
      enddo
    end associate

    ! upper boundary condition
    this%centre(nodes) = (d(nodes)*gamma1(2) + v(nodes) )/ delta(nodes) + volume_vector(nodes)
    this%lower(nodes) =  -2._wp / (delta(nodes-1)/d(nodes-1) + delta(nodes)/d(nodes))
    if( upwind_direction .eqv. .true. ) this%lower(nodes-1) = this%lower(nodes-1) - v(nodes-1)
    this%lower(nodes) = this%lower(nodes) / delta(nodes)

  end subroutine

  subroutine initialise_1d_mass_matrix_const( this, diffusive, advective, volume, cell_spacings, nodes, upwind_direction, gamma1, sigma )
    class(fv1DTemporalSpatialClass), intent(out) :: this
    real(wp), intent(in) :: diffusive
    real(wp), intent(in) :: advective
    real(wp), intent(in) :: volume
    real(wp), intent(in) :: cell_spacings
    real(wp), intent(in) :: gamma1(2) ! du/dn + gamma u = sigma form of boundary
    real(wp), intent(in) :: sigma(2)
    integer, intent(in) :: nodes
    logical, intent(in) :: upwind_direction ! 1 if flow is going in direction of increasing node size, 0 otherwise

    real(wp) :: matrix_coefficients(3)
    integer :: status

    allocate( this%cell_spacings(nodes), this%lower(nodes), this%centre(nodes), this%upper(nodes), this%source(nodes), stat=status )
    if( status /= 0 ) stop 'Error allocating cell spacings in 1D finite volume matrix generation'

    this%cell_spacings = cell_spacings
    this%lower = 0._wp; this%upper = 0._wp; this%centre = 0._wp ; this%source = 0._wp
    this%nodes = nodes

    associate( d => diffusive, v => advective, delta => cell_spacings )

      ! lower boundary condition
      this%centre(1) = (d*gamma1(1) + v )/ delta + volume_vector
      this%upper(1) = -2._wp / (delta/d + delta/d)
      if( upwind_direction .eqv. .true. ) this%upper(1) = this%upper(1) - v
      this%upper(1) = this%upper(1) / delta

      do i = 2, nodes - 1

        ! diffusive terms

        ! node in the direction of decreasing spatial indexing
        this%lower(i) = -2._wp / (delta/d + delta/d)
        ! node in the direction of increasing spatial indexing
        this%upper(i) = -2._wp / (delta/d + delta/d)
        ! centre of the cell
        this%centre(i) = volume_vector - (this%lower(i)+this%upper(i))

        ! advective terms
        if( upwind_direction .eqv. .true. ) then
          this%lower(i) = this%lower(i) - v
          this%centre(i) = this%centre(i) + v
        else
          this%upper(i) = this%upper(i) - v
          this%centre(i) = this%centre(i) + v
        endif

        ! divide by spacing
        this%lower(i) = this%lower(i) / delta
        this%centre(i) = this%centre(i) / delta
        this%upper(i) = this%upper(i) / delta
      enddo
    end associate

    ! upper boundary condition
    this%centre(nodes) = (d*gamma1(2) + v )/ delta + volume_vector
    this%lower(nodes) =  -2._wp / (delta/d + delta()/d())
    if( upwind_direction .eqv. .true. ) this%lower(nodes) = this%lower(nodes) - v
    this%lower(nodes) = this%lower(nodes) / delta

  end subroutine

  subroutine initialise_source_vector( this, source_vector, cell_spacings )
    class(fv1fv1DTemporalSpatialClass), intent(inout) :: this
    real(wp), intent(in) :: source_vector(:)

    this%source = source_vector
  end subroutine

  subroutine solve_for_gradient( this, input_vector, gradient_vector )
    class(fv1DTemporalSpatialClass), intent(in) :: this
    real(wp), intent(in) :: input_vector(:)
    real(wp), intent(out) :: gradient_vector(:)

    integer :: i

    do i = 2, this%nodes - 1
      gradient_vector(i) = this%lower(i)*input_vector(i-1) + this%centre(i)*input_vector(i) + &
                        this%upper(i)*input_vector(i+1) + this%source(i)
    enddo

  end subroutine
end module
