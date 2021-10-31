!-------------------------------------------------------------------------------
! This is in INSECT%STATE:
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
! STATE(11) : x-angular velocity of body *not sure if around the body mass center?
! STATE(12) : y-angular velocity of body
! STATE(13) : z-angular velocity of body
! STATE(14) : angle x -- around origin
! STATE(15) : angle y -- around origin
! STATE(16) : angle z -- around origin
! STATE(17) : angular velocity -- around origin
! STATE(18) : angular velocity -- around origin
! STATE(19) : angular velocity -- around origin
! STATE(20) : 
! STATE(21) :
! STATE(22) :
! STATE(23) :
! STATE(24) :
! STATE(25) : 
!-------------------------------------------------------------------------------
! Insect free flight dynamics.
! RHS of the ODE system.
! TASK: from STATE compute RHS
!-------------------------------------------------------------------------------
subroutine rigid_solid_rhs_fractaltree(time, it, state, rhs, force_g, torque_g, Insect)
    implicit none

    ! T0 DO
    ! check dimensionality
    ! check i convert back to global frame and if thats necessary or whatever
    integer, intent(in) :: it
    real(kind=rk),intent(in) :: time, force_g(1:3), torque_g(1:3)
    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(inout) :: state(1:Insect%state_array_len)
    real(kind=rk),intent(inout) :: rhs(1:Insect%state_array_len)
    real(kind=rk) :: Jx, Jy,Jz,Jxy, Tx, Ty, Tz, T_wing_g(1:3), T_wing_w(1:3)
    real(kind=rk) :: omx,omy,omz, Jxx,Jyy,Jzz
    real(kind=rk), dimension(0:3) :: ep
    real(kind=rk), dimension(1:3) :: ROT, torque_body
    real(kind=rk), dimension(1:3,1:3) :: M_body, M_wing_l
    ! 1DOF system, oscillating in radial direction
    real(kind=rk) :: m, l, a, b, s, res_coef
    ! initialization
    rhs = 0.d0

    if (Insect%BodyMotion /= "free_flight") then
        call abort(900,"Insect%BodyMotion"//trim(adjustl(Insect%BodyMotion))//" but using free-flight?")
    endif

    ! copy some shortcuts (this is easier to code)
    m  = Insect%mass
    ! l is length of pendulum arm
    l =  8
    
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

    res_coef = 10
    secondquadrant: if (Insect%STATE(1)<=0) then
        a = -1 ! restorative force coefficient
    else !first quadrant or zero
        a = 1 !swinging to positive direction force is negtive      
    endif secondquadrant
    b = 0.5_rk
    !Insect%STATE(3) = global_position(3) + l*sin(thetaX)

    ! integrate coordinates (dx/dt = vx) Note: this is in global reference frame
    rhs(1) = Insect%STATE(4) ! velocity x
    rhs(2) = Insect%STATE(5) !velocity y
    rhs(3) = -(Insect%STATE(1)-Insect%x0(1))/SQRT(l**2 - (Insect%STATE(1) - Insect%x0(1))**2) ! velocity z
    ! integrate velocities (dvx/dt = F) Note: this is in global reference frame
    !         drag                  gravity             dissipation           restorative
    rhs(4) = s*force_g(1)/(m*l) + Insect%gravity_x + b*Insect%STATE(4)*Insect%STATE(3)/((l**2)*m) + a*res_coef*((pi/2) - asin(Insect%STATE(3)/l))*Insect%STATE(3)/l
    rhs(5) = s*force_g(2)/(m*l) + Insect%gravity_y
    rhs(6) = s*force_g(3)/(m*l) + Insect%gravity   - b*Insect%STATE(6)*SQRT(1 - Insect%STATE(3)**2/l**2)/(l*m) - a*res_coef*((pi/2) - asin(Insect%STATE(3)/l))*SQRT(1 - Insect%STATE(3)**2/l**2)
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
end subroutine



subroutine rigid_solid_init_fractaltree(time, Insect, resume_backup)
    implicit none

    real(kind=rk), intent(in) :: time
    type(diptera), intent(inout) :: Insect
    logical, intent(in) :: resume_backup
    real(kind=rk) :: yaw,pitch,roll,a,t,p
    real(kind=rk), ALLOCATABLE :: array(:,:)
    integer :: mpicode, n_lines, n_cols, n_header, it, n_candidates
    real(kind=rk), dimension(0:3) :: ep
    integer :: state_length 
    
    state_length = Insect%state_array_len + 1

    Insect%time = time
    Insect%STATE = 0.d0
    n_header = 0

    if (root) write(*,'(A)') "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
    if (root) write(*,'("rigid solid init at time=",es12.4)')  Insect%time
    if (root) write(*,'(A)') "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"

    if (resume_backup) then
        ! resuming the rigid solid solver from a backup
        ! In WABBIT, we write the state vector to *.t file in every time step, and so we look for the
        ! data in this file. As simulations may fail and be resumed from a different time step, we use the
        ! last suitable entry in this file.
        if (root) then
            write(*,'(A)') "Fractal tree solver is resuming from file: we read Insect%STATE from ./insect_state_vector.t"
            write(*,*) "time=", time

            call count_lines_in_ascii_file('insect_state_vector.t', n_lines, n_header)
            call count_cols_in_ascii_file('insect_state_vector.t', n_cols, n_header)

            if (n_cols < state_length) call abort(202117021, "For some reason insect_state_vector.t contains not enough columns....something is wrong?")

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
                    Insect%STATE = array(it,2:state_length)
                    write(*,*) "Found suitable entry in line=", it, " time=", array(it,1)
                    !write(*,'("Insect%STATE=(",25(es15.8,1x),")")')  Insect%STATE 
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
            Insect%STATE(19) = 0.d0    ! y-component ang. velocity
            Insect%STATE(20) = 0.d0    ! z-component ang. velocity
        endif


        

    endif ! of backup/no backup if


    if(root) write(*,'("Insect%STATE=(",25(es15.8,1x),")")')  Insect%STATE
end subroutine
