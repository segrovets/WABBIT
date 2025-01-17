!> \brief call subroutines that read mesh and fields as initial condition from files
!! input:    - parameter array
!! output:   - light data array
!!           - heavy data array
!!           - number of active blocks (light and heavy)
!!           - time and iteration
! ********************************************************************************************

subroutine get_inicond_from_file(params, lgt_block, hvy_block, hvy_n, lgt_n, time, iteration)

    implicit none
    type (type_params), intent(in)        :: params                     !> user defined parameter structure
    integer(kind=ik), intent(inout)       :: lgt_block(:, :)            !> light data array
    real(kind=rk), intent(inout)          :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
    integer(kind=ik), intent(inout)       :: hvy_n, lgt_n               !> number of heavy and light active blocks
    real(kind=rk), intent(inout)          :: time                       !> time loop variables
    integer(kind=ik), intent(inout)       :: iteration
    real(kind=rk), dimension(3)           :: domain
    real(kind=rk)                         :: t0                         ! cpu time variables for running time calculation
    integer(kind=ik)                      :: N_files                    ! number of files to read from
    integer(kind=ik)                      :: dF, tc_length, dim, i      ! loop variable
    integer(kind=ik), dimension(3)        :: Bs
    logical                               :: periodic_BC(1:3), symmetry_BC(1:3)

    ! number of files to read from
    N_files = params%n_eqn
    ! start time
    t0 = MPI_wtime()

    if (params%rank==0) write(*,*) "Reading initial condition from file"

    ! read time, iteration, domain size and total number of blocks from first input file
    call read_attributes(params%input_files(1), lgt_n, time, iteration, domain, Bs, &
    tc_length, dim, periodic_BC=periodic_BC, symmetry_BC=symmetry_BC)

    ! print time, iteration and domain on screen
    if (params%rank==0) then
        write(*,'(80("_"))')
        write(*,'("READING: Reading from file ",A)') trim(adjustl(params%input_files(1)))
        write(*,'("time=",g12.4," iteration=", i5)') time, iteration
        write(*,'("Bs(1)=",i4,"Bs(2)=",i4,"Bs(3)=",i4," dim=", i1, " tc_length=", i2)') Bs(1),Bs(2),Bs(3), dim, tc_length
        write(*,'("Lx=",g12.4," Ly=",g12.4," Lz=",g12.4)') domain
        ! if the domain size doesn't match, proceed, but yell.
        if ((abs(params%domain_size(1)-domain(1))>1e-12_rk).or.(abs(params%domain_size(2)-domain(2))>1e-12_rk) &
            .or.(abs(params%domain_size(3)-domain(3))>1e-12_rk)) then
            write (*,'(A)') " WARNING! Domain size mismatch."
            write (*,'("in memory:   Lx=",es12.4,"Ly=",es12.4,"Lz=",es12.4)') &
                params%domain_size(1), params%domain_size(2), params%domain_size(3)
            write (*,'("but in file: Lx=",es12.4,"Ly=",es12.4,"Lz=",es12.4)') domain
            write (*,'(A)') "proceed, with fingers crossed."
        end if
    end if

    do i = 1, params%dim
        if (periodic_BC(i).neqv.params%periodic_BC(i)) then
            if (params%rank==0) write(*,*) "WARNING WARNING WARNING periodic_BC in inifile and HDF5 file are different!"
        endif

        if (symmetry_BC(i).neqv.params%symmetry_BC(i)) then
            if (params%rank==0) write(*,*) "WARNING WARNING WARNING symmetry_BC in inifile and HDF5 file are different!"
        endif
    enddo

    if (lgt_n > size(lgt_block,1)) then
        call abort(743734, 'ERROR: Not enough memory allocated for the saved field!')
    endif

    ! read treecode from first input file
    call read_mesh(params%input_files(1), params, lgt_n, hvy_n, lgt_block)

    ! read datafields from files into hvy_block array
    do dF = 1, N_files
        call read_field(params%input_files(dF), dF, params, hvy_block, hvy_n )
    end do

    ! timing
    call toc( "read_data", MPI_wtime()-t0 )
end subroutine get_inicond_from_file
