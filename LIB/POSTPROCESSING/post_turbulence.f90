!> \brief postprocessing routine for subsequent dissipation calculation from datafields ux, uy (, uz) saved in .h5 files
!-----------------------------------------------------------------------------------------------------

subroutine post_turbulence(params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_forest
    use module_mpi
    use module_operators
   ! use module_t_files

    ! for use of code that determines the domain_n vector for finding block ajacency to BC
    use module_treelib

    !for use of parameters from ini file
    use module_ini_files_parser_mpi

    implicit none

    !> parameter struct
    type (type_params), intent(inout)      :: params
    character(len=cshort)                  :: file_ux, file_uy, file_uz, operator, ini_file
    real(kind=rk)                          :: time
    integer(kind=ik)                       :: iteration, k, lgt_id, lgt_n, hvy_n, tc_length, g
    integer(kind=ik), dimension(3)         :: Bs
    character(len=2)                       :: order

    integer(kind=ik), allocatable      :: lgt_block(:, :)
    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :)
    integer(kind=ik), allocatable      :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable      :: lgt_active(:), hvy_active(:)
    integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:)
    character(len=cshort)              :: fname
    real(kind=rk), dimension(3)        :: dx, x0
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain
    integer(kind=ik)                   :: nwork
    !real(kind=rk)                      :: dissipation_sum = 0.0_rk
    real(kind=rk)                      :: nu, L_sponge
    integer(kind=2)                    :: domain_n(1:3)
    integer(kind=ik)                   :: N_mask_components = 0_ik
    integer(kind=ik)                   :: N_sponge_cells = 1_ik
    type(inifile)                      :: FILE
    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename and level for interpolation)
    call get_command_argument(1, operator)
    call get_command_argument(2, file_ux)
    ! does the user need help?
    if (file_ux=='--help' .or. file_ux=='--h') then
        if (params%rank==0) then
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " Wabbit postprocessing: turbulent energy dissipation"
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " Computes quantity from velocity files. Output is stored"
            write(*,*) " in predefined files."
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " --energy-dissipation"
            write(*,*) "./wabbit-post --energy-dissipation source_ux.h5 source_uy.h5 [source_uz.h5] [ORDER] [Re]"
            write(*,*) " Computes (3D) energy dissipation, saves in " !or 1 (2D) vorticity component
            write(*,*) " diss_*.h5"
            write(*,*) " order = 2 or 4"
            write(*,*) "-----------------------------------------------------------"
        end if
        return
    endif

    call check_file_exists(trim(file_ux))

    call get_command_argument(3, file_uy)
    call check_file_exists(trim(file_uy))

    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(file_ux, lgt_n, time, iteration, domain, Bs, tc_length, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    if (params%dim == 3) then
        call get_command_argument(4, file_uz)
        call check_file_exists(trim(file_uz))
        call get_command_argument(5, order)
        call get_command_argument(6, ini_file)
    else
        call get_command_argument(4, order)
        call get_command_argument(5, ini_file)
    end if

    ! read the file, only process 0 should create output on screen
    call set_lattice_spacing_mpi(1.0d0)
    call read_ini_file_mpi(FILE, ini_file, .true.)
    call read_param_mpi(FILE, 'ACM-new', 'nu', nu, 1e-1_rk)
    call read_param_mpi(FILE, 'Sponge', 'L_sponge', L_sponge, 1.0_rk )
    call clean_ini_file_mpi(FILE)
 

    ! decide which order for discretization and predictor is used. Note predictor
    ! is used in ghost nodes syncing
    if (order == "4") then
        params%order_discretization = "FD_4th_central_optimized"
        params%order_predictor = "multiresolution_4th"
        params%n_ghosts = 4_ik

    elseif (order == "2") then
        params%order_discretization = "FD_2nd_central"
        params%order_predictor = "multiresolution_2nd"
        params%n_ghosts = 2_ik

    else
        call abort(8765,"chosen discretization order invalid or not (yet) implemented. choose between 4 (FD_4th_central_optimized) and 2 (FD_2nd_central)")

    end if

    params%max_treelevel  = tc_length
    params%n_eqn          = params%dim
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)
    params%Bs = Bs
    allocate(params%butcher_tableau(1,1))

    allocate(params%symmetry_vector_component(1:params%n_eqn))
    params%symmetry_vector_component(1) = "x"
    params%symmetry_vector_component(2) = "y"
    if (params%dim==3) then
        params%symmetry_vector_component(3) = "z"
    endif

    Bs = params%Bs
    g  = params%n_ghosts

    ! no refinement is made in this postprocessing tool; we therefore allocate about
    ! the number of blocks in the file (and not much more than that)
    params%number_blocks = ceiling(  real(lgt_n)/real(params%number_procs) )

    nwork = 1

    ! allocate data
    call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, &
    lgt_active, hvy_active, lgt_sortednumlist, hvy_tmp=hvy_tmp, neqn_hvy_tmp=nwork)

    ! read mesh
    call read_mesh(file_ux, params, lgt_n, hvy_n, lgt_block)

    ! read actual data (velocity)
    call read_field(file_ux, 1, params, hvy_block, hvy_n)
    call read_field(file_uy, 2, params, hvy_block, hvy_n)

    if (params%dim == 3) then
        call read_field(file_uz, 3, params, hvy_block, hvy_n)
    end if

    ! create lists of active blocks (light and heavy data)
    ! update list of sorted nunmerical treecodes, used for finding blocks
    call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_ID=1)

    ! update neighbor relations
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_active, hvy_n )

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

    ! calculate energy dissipation
    do k = 1, hvy_n
        call hvy2lgt(lgt_id, hvy_active(k), params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

        ! get the surface normals (determine BC adjacency)
        call get_adjacent_boundary_surface_normal( lgt_block(lgt_id, 1:lgt_block(lgt_id,params%max_treelevel+IDX_MESH_LVL)), &
        params%domain_size, params%Bs, params%dim, domain_n )

        ! determine extent of sponge layer
        N_sponge_cells = NINT(L_sponge / dx(1))

        if (operator == "--energy-dissipation") then
            call compute_energy_dissipation( hvy_block(:,:,:,1,hvy_active(k)), &
            hvy_block(:,:,:,2,hvy_active(k)), hvy_block(:,:,:,3,hvy_active(k)),&
            dx, Bs, g, params%order_discretization, hvy_tmp(:,:,:,1,hvy_active(k)),&
            domain_n, N_sponge_cells,&
            nu, params%rank)
        else
            call abort(1812011, "operator is not --energy-dissipation")

        endif
    end do

    if (operator == "--energy-dissipation") then
        write( fname,'(a, "_", i12.12, ".h5")') 'diss', nint(time * 1.0e6_rk)

        call write_field(fname, time, iteration, 1, params, lgt_block, hvy_tmp, lgt_active, lgt_n, hvy_n, hvy_active )
        !call append_t_file('dissipation.t', (/time*1.0e6_rk, dissipation_sum /))

    endif

    call close_all_t_files()
end subroutine post_turbulence