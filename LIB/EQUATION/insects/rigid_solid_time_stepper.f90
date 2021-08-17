!-------------------------------------------------------------------------------
! This is in INSECT%STATE:
! STATE(1) : x-position of body
! STATE(2) : y-position of body
! STATE(3) : z-position of body
! STATE(4) : x-velocity of body
! STATE(5) : y-velocity of body
! STATE(6) : z-velocity of body
! STATE(7) : 1st component of body quaternion
! STATE(8) : 2nd component of body quaternion
! STATE(9) : 3rd component of body quaternion
! STATE(10) : 4th component of body quaternion
! STATE(11) : x-angular velocity of body
! STATE(12) : y-angular velocity of body
! STATE(13) : z-angular velocity of body
! STATE(14) : 1st component of left wing quaternion
! STATE(15) : 2nd component of left wing quaternion
! STATE(16) : 3rd component of left wing quaternion
! STATE(17) : 4th component of left wing quaternion
! STATE(18) : x-angular velocity of left wing
! STATE(19) : y-angular velocity of left wing
! STATE(20) : z-angular velocity of left wing
!-------------------------------------------------------------------------------
! Insect free flight dynamics.
! RHS of the ODE system.
! TASK: from STATE compute RHS
!-------------------------------------------------------------------------------

subroutine rigid_solid_rhs(time, it, state, rhs, force_g, torque_g, Insect)
    implicit none
    integer, intent(in) :: it
    real(kind=rk),intent(in) :: time, force_g(1:3), torque_g(1:3)
    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(inout) :: state(1:Insect%state_array_len)
    real(kind=rk),intent(inout) :: rhs(1:Insect%state_array_len)

    ! if Insect%use_dynamic_tree then call fractal_tree_time_stepper
    if (Insect%fractal_tree_motion) then
        call rigid_solid_rhs_fractaltree(time, it, state, rhs, force_g, torque_g, Insect)
    else
        call rigid_solid_rhs_insect(time, it, state, rhs, force_g, torque_g, Insect)
    endif
end subroutine

subroutine rigid_solid_init(time, Insect, resume_backup)
    implicit none

    real(kind=rk), intent(in) :: time
    type(diptera), intent(inout) :: Insect
    logical, intent(in) :: resume_backup

    if (Insect%fractal_tree_motion) then
        call rigid_solid_init_fractaltree(time, Insect, resume_backup)
    else
        call rigid_solid_init_insect(time, Insect, resume_backup)
    endif
end subroutine

subroutine rigid_solid_rhs_insect(time, it, state, rhs, force_g, torque_g, Insect)
    implicit none

    integer, intent(in) :: it
    real(kind=rk),intent(in) :: time, force_g(1:3), torque_g(1:3)
    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(inout) :: state(1:Insect%state_array_len)
    real(kind=rk),intent(inout) :: rhs(1:Insect%state_array_len)
    real(kind=rk) :: m, Jx, Jy,Jz,Jxy, s, Tx, Ty, Tz, T_wing_g(1:3), T_wing_w(1:3)
    real(kind=rk) :: omx,omy,omz, Jxx,Jyy,Jzz
    real(kind=rk), dimension(0:3) :: ep
    real(kind=rk), dimension(1:3) :: ROT, torque_body
    real(kind=rk), dimension(1:3,1:3) :: M_body, M_wing_l

    ! initialization
    rhs = 0.d0

    if (Insect%BodyMotion /= "free_flight") then
        call abort(900,"Insect%BodyMotion"//trim(adjustl(Insect%BodyMotion))//" but using free-flight?")
    endif

    ! copy some shortcuts (this is easier to code)
    m  = Insect%mass
    Jx = Insect%Jroll_body
    Jy = Insect%Jpitch_body
    Jz = Insect%Jyaw_body

    ! extract rotation (unit) quaternion from state vector, and create the rotation
    ! matrix from it.
    ep = Insect%STATE(7:10)
    call rotation_matrix_from_quaternion( ep, M_body )
    ! The equations of motion for the rotation are written in the body reference
    ! frame, thus the fluid torque has to be transformed to the body system
    torque_body = matmul( M_body, torque_g )

    ! extract angular velocity vector from state vector
    ! note this is actually rot_body_b
    ROT = Insect%STATE(11:13)

    ! startup conditioner (to avoid problems with impulsively started motion)
    if (Insect%startup_conditioner=="yes") then
        s = startup_conditioner(time, 0.1d0, 0.5d0)
    else
        s = 1.d0
    endif

    !-----------------------------------------------------------------------------
    ! Integrate the 13 equations of motion (written in quaternion formulation)
    ! The underlying eqns can be found in
    ! M. Maeda et al. (2012) A free-flight Simulation of Insect flapping flight
    ! J. Aero Aqua Bio-Mech(1):1,71-79
    !-----------------------------------------------------------------------------
    ! To avoid confusion with the famous unsteady corrections, the actual translation eqn reads:
    !   rho_s*V*u_dot = (rho_s-rho_f)*V*g + F_fluid
    ! but the latter term is
    !   rho_s*V*u_dot = (rho_s-rho_f)*V*g + Integral(Penal) + V*rho_f*u_dot
    ! so we can put this, as uhlmann does, on the right hand side:
    !   (rho_s-rho_f)*V*u_dot = (rho_s-rho_f)*V*g + Integral(Penal)
    ! We have thus two options regarding this, either
    !   m_corrected * u_dot = m_corrected * g + Integral(penal)       [implicit unst corrections]
    ! or
    !   m * u_dot = m_corrected * g + Integral(penal) + force_unst    [explicit unst corrections]
    ! the last equation can also be rewritten using twice the same m
    !   m * u_dot = m * g_corrected + Integral(penal) + force_unst    [explicit unst corrections]
    ! where g_corrected = (rho_s-rho_f)/rho_s * g
    !-----------------------------------------------------------------------------
    ! For the torque, the unsteady corrections are explicitly taken into account
    ! above, thus cal_unst=1; is required in the parameter file.
    !-----------------------------------------------------------------------------
    ! integrate coordinates (dx/dt = vx) Note: this is in global reference frame
    rhs(1) = Insect%STATE(4)
    rhs(2) = Insect%STATE(5)
    rhs(3) = Insect%STATE(6)
    ! integrate velocities (dvx/dt = F) Note: this is in global reference frame
    rhs(4) = s*force_g(1)/m + Insect%gravity_x
    rhs(5) = s*force_g(2)/m + Insect%gravity_y
    rhs(6) = s*force_g(3)/m + Insect%gravity
    ! integrate quaternion attitudes
    rhs(7)  = 0.5d0*(-ep(1)*ROT(1)-ep(2)*ROT(2)-ep(3)*ROT(3))
    rhs(8)  = 0.5d0*(+ep(0)*ROT(1)-ep(3)*ROT(2)+ep(2)*ROT(3))
    rhs(9)  = 0.5d0*(+ep(3)*ROT(1)+ep(0)*ROT(2)-ep(1)*ROT(3))
    rhs(10) = 0.5d0*(-ep(2)*ROT(1)+ep(1)*ROT(2)+ep(0)*ROT(3))
    ! integrate angular velocities
    rhs(11) = ( (Jy-Jz)*ROT(2)*ROT(3) + s*torque_body(1) )/Jx
    rhs(12) = ( (Jz-Jx)*ROT(3)*ROT(1) + s*torque_body(2) )/Jy
    rhs(13) = ( (Jx-Jy)*ROT(1)*ROT(2) + s*torque_body(3) )/Jz

    ! turn on or off degrees of freedom for free flight solver. The string from
    ! ini file contains 6 characters 1 or 0 that turn on/off x,y,z,yaw,pitch,roll
    ! degrees of freedom by multiplying the respective RHS by zero, keeping the
    ! value thus constant
    rhs(4) = rhs(4) * Insect%DoF_on_off(1)   ! x translation
    rhs(5) = rhs(5) * Insect%DoF_on_off(2)   ! y translation
    rhs(6) = rhs(6) * Insect%DoF_on_off(3)   ! z translation
    rhs(13) = rhs(13) * Insect%DoF_on_off(4) ! yaw rotation
    rhs(12) = rhs(12) * Insect%DoF_on_off(5) ! pitch rotation
    rhs(11) = rhs(11) * Insect%DoF_on_off(6) ! roll rotation

    if (insect%gravity_y/= 0.0) then
        call append_t_file('forces_rk.t', (/time, rhs(2), rhs(5), force_g(2)/))
    else
        call append_t_file('forces_rk.t', (/time, rhs(3), rhs(6), force_g(3)/))
    endif

    !**********************************
    !** Wing fsi model               **
    !**********************************
    ! wing fsi solver is used. we therefore have to solve 7 equations for the left
    ! wing, which are the quaternion components and the angular velocity
    ! if (Insect%wing_fsi == "yes") then
    !     ! the wing fsi model must assume resting body: if not, we have a more complicated
    !     ! multi-body problem to solve
    !     if (maxval(Insect%DoF_on_off) == 1) then
    !         call abort(1313,"you must not use wing fsi AND body fsi at the same time")
    !     endif
    !     ! for the model, we indeed need the wing inertia..
    !     if (maxval( (/Insect%Jxx,Insect%Jyy,Insect%Jzz,Insect%Jxy/) ) < 1.0d-8) then
    !         call abort(1314,"you forgot to set wing inertia tensor...")
    !     endif
    !
    !     Insect%STATE(19:20) =Insect%STATE(19:20) *0.d0
    !
    !     ! define abbreviations for simplicity
    !     Jxx = Insect%Jxx
    !     Jyy = Insect%Jyy
    !     Jzz = Insect%Jzz
    !     Jxy = Insect%Jxy
    !     ep = Insect%STATE(14:17)
    !
    !     call body_rotation_matrix( Insect, M_body )
    !     call wing_left_rotation_matrix( Insect, M_wing_l )
    !     ! get the current muscle moment
    !     call muscle_moment(time, Insect)
    !
    !     ! aero parts
    !     T_wing_g = Insect%PartIntegrals(2)%Torque(1:3)
    !     T_wing_w = matmul(M_wing_l,matmul(M_body,T_wing_g)) + Insect%torque_muscle_l_b
    !     ! T_wing_w = matmul(M_wing_l,matmul(M_body,T_wing_g)) + matmul(M_wing_l, Insect%torque_muscle_l_b)
    !
    !     ! T_wing_w = T_wing_w +  Insect%torque_muscle_l_b
    !
    !     ! total torque in wing system with startup conditioner 's'
    !     Tx = T_wing_w(1) * s
    !     Ty = T_wing_w(2) * s
    !     Tz = T_wing_w(3) * s
    !
    !     ! angular velocity components
    !     omx = Insect%STATE(18)
    !     omy = Insect%STATE(19)
    !     omz = Insect%STATE(20)
    !
    !     Insect%RHS_this(14) = (-ep(1)*omx-ep(2)*omy-ep(3)*omz)/2.d0   ! 1st quaternion component
    !     Insect%RHS_this(15) = (+ep(0)*omx-ep(3)*omy+ep(2)*omz)/2.d0   ! 2nd quaternion component
    !     Insect%RHS_this(16) = (+ep(3)*omx+ep(0)*omy-ep(1)*omz)/2.d0   ! 3rd quaternion component
    !     Insect%RHS_this(17) = (-ep(2)*omx+ep(1)*omy+ep(0)*omz)/2.d0   ! 4th quaternion component
    !
    !     ! x-component ang. velocity
    !     Insect%RHS_this(18) = Jxy*(Ty-omz*(Jxx*omx+Jxy*omy)+Jzz*omx*omy) -(Jyy*(Tx+omz*(Jxy*omx+Jyy*omy)-Jzz*omy*omz))
    !     Insect%RHS_this(18) = Insect%RHS_this(18) / (Jxy**2-Jxx*Jyy)
    !     ! y-component ang. velocity
    !     Insect%RHS_this(19) = Jxy*(Tx+omz*(Jyy*omy+Jxy*omx)-Jzz*omy*omz) -(Jxx*(Ty-omz*(Jxy*omy+Jxx*omx)+Jzz*omx*omz))
    !     Insect%RHS_this(19) = Insect%RHS_this(19) / (Jxy**2-Jxx*Jyy)
    !     ! z-component ang. velocity
    !     Insect%RHS_this(20) = (Tz+omy*(Jxx*omx+Jxy*omy)-omx*(Jxy*omx+Jyy*omy))/Jzz
    !
    !     if (root) then
    !         write(*,'("t=",f8.4,2x,3(es12.4,1x),"// ",3(es12.4,1x),"// ",3(es12.4,1x),&
    !         &" alpha=",f6.2," phi=",f6.2," theta=",f6.2)') &
    !         time, T_wing_w, Insect%torque_muscle_l_b, Insect%STATE(18:20), &
    !         rad2deg(Insect%alpha_l), rad2deg(Insect%phi_l), rad2deg(Insect%theta_l)
    !     endif
    ! endif
end subroutine

subroutine rigid_solid_init_insect(time, Insect, resume_backup)
    implicit none

    real(kind=rk), intent(in) :: time
    type(diptera), intent(inout) :: Insect
    logical, intent(in) :: resume_backup
    real(kind=rk) :: yaw,pitch,roll,a,t,p
    real(kind=rk), ALLOCATABLE :: array(:,:)
    integer :: mpicode, n_lines, n_cols, n_header, it, n_candidates
    real(kind=rk), dimension(0:3) :: ep

    Insect%time = time
    Insect%STATE = 0.d0
    n_header = 0

    if (root) write(*,'(A)') "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
    if (root) write(*,'("rigid solid init at time=",es12.4)')  Insect%time
    if (root) write(*,'(A)') "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"

    if (resume_backup) then
        ! resuming the rigid solid solver from a backup
        ! NOTE: in old versions, this read a single *.fsi_bckp file which contained the state vector.
        ! In WABBIT, we write the state vector to *.t file in every time step, and so we look for the
        ! data in this file. As simulations may fail and be resumed from a different time step, we use the
        ! last suitable entry in this file.
        if (root) then
            write(*,'(A)') "Rigid solid solver is resuming from file: we read Insect%STATE from ./insect_state_vector.t"
            write(*,*) "time=", time

            call count_lines_in_ascii_file('insect_state_vector.t', n_lines, n_header)
            call count_cols_in_ascii_file('insect_state_vector.t', n_cols, n_header)

            if (n_cols < 21) call abort(202117021, "For some reason insect_state_vector.t contains not enough columns....something is wrong?")

            allocate( array(1:n_lines, 1:n_cols) )
            call read_array_from_ascii_file('insect_state_vector.t', array, n_header)

            n_candidates = 0
            do it = 1, n_lines
                if ( abs(array(it,1)-time) <= 1.0e-6 ) n_candidates = n_candidates +1
            end do

            write(*,*) "In insect_state_vector.t we found ", n_candidates, "possible time stamps and we use the last one!"

            if (n_candidates == 0) then
                call abort(20210291, "Resuming from insect_state_vector.t was impossible as no time stamp is sufficiently close to what we need.")
            endif

            do it = n_lines, 1, -1
                if ( abs(array(it,1)-time) <= 1.0e-6 ) then
                    Insect%STATE = array(it,2:21)
                    write(*,*) "Found suitable entry in line=", it, " time=", array(it,1)
                    write(*,'("Insect%STATE=(",20(es15.8,1x),")")')  Insect%STATE
                    exit
                endif
            end do

            deallocate(array)
        endif

        call MPI_BCAST( Insect%time, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
        call MPI_BCAST( Insect%STATE, size(Insect%STATE), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )

    else ! no backup

        ! free flight solver based on quaternions. the task here is to initialize
        ! the "attitude" quaternion (Insect%quaternion) from yaw, pitch and roll
        ! angles. Note that for the free flight solver, this is the last time
        ! body yaw,pitch,roll are meaningful. from now on, quaternions are used

        ! initialization, these values are read from parameter file
        Insect%gamma = Insect%yawpitchroll_0(1)
        Insect%beta  = Insect%yawpitchroll_0(2)
        Insect%psi   = Insect%yawpitchroll_0(3)
        Insect%xc_body_g = Insect%x0
        Insect%vc_body_g = Insect%v0
        Insect%rot_body_b = 0.d0

        ! create initial value for attitude quaternion
        yaw   =  Insect%gamma / 2.d0
        pitch =  Insect%beta  / 2.d0
        roll  =  Insect%psi   / 2.d0
        Insect%STATE(7) = cos(roll)*cos(pitch)*cos(yaw) + sin(roll)*sin(pitch)*sin(yaw)
        Insect%STATE(8) = sin(roll)*cos(pitch)*cos(yaw) - cos(roll)*sin(pitch)*sin(yaw)
        Insect%STATE(9) = cos(roll)*sin(pitch)*cos(yaw) + sin(roll)*cos(pitch)*sin(yaw)
        Insect%STATE(10) = cos(roll)*cos(pitch)*sin(yaw) - sin(roll)*sin(pitch)*cos(yaw)

        Insect%STATE(1:3) = Insect%xc_body_g
        Insect%STATE(4:6) = Insect%vc_body_g
        Insect%STATE(11:13) = Insect%rot_body_b

        !**********************************
        !** Wing fsi model               **
        !**********************************
        ! if in use, inititalize the wing-fsi solver here. The tasks are similar to the
        ! free flight solver:
        ! (a) initialize the quaternion state from the initial values of the wing angles
        ! (c) initialize angular velocity and acceleration (maybe zero?)
        ! (d) if applicable, read muscle moment in body system from file
        if ( Insect%wing_fsi == "yes" ) then
            ! the intial angles are read from the parameter file. note division by 2 is
            ! because of the quaternion, it does not divide the desired angle by two.
            a = deg2rad( Insect%init_alpha_phi_theta(1) ) / 2.d0  ! alpha
            p = deg2rad( Insect%init_alpha_phi_theta(2) ) / 2.d0  ! phi
            t = deg2rad( Insect%init_alpha_phi_theta(3) ) / 2.d0  ! theta

            Insect%STATE(14) = cos(a)*cos(t)*cos(p) + sin(a)*sin(t)*sin(p)  ! 1st quaternion component
            Insect%STATE(15) = cos(a)*cos(t)*sin(p) - sin(a)*cos(p)*sin(t)  ! 2nd quaternion component
            Insect%STATE(16) =-cos(a)*sin(t)*sin(p) + cos(t)*cos(p)*sin(a)  ! 3rd quaternion component
            Insect%STATE(17) = cos(a)*cos(p)*sin(t) + sin(a)*cos(t)*sin(p)  ! 4th quaternion component

            Insect%STATE(18) = 0.d0    ! x-component ang. velocity
            Insect%STATE(19) = 0.d0    ! x-component ang. velocity
            Insect%STATE(20) = 0.d0    ! x-component ang. velocity
        endif
    endif ! of backup/no backup if

    if(root) write(*,'("Insect%STATE=(",20(es15.8,1x),")")')  Insect%STATE
end subroutine
