!> \brief compute energy dissipation for time step t (for saving it on disk)
! ********************************************************************************************
subroutine compute_energy_dissipation(u, v, w, dx, Bs, g, discretization, dissipation,  nu, rank)
    use module_precision

    implicit none

    real(kind=rk), dimension(3), intent(in)        :: dx                        !> spacing of the block
    real(kind=rk), dimension(:,:,:), intent(in)    :: u, v, w                   !> local datafields
    real(kind=rk), intent(in)                      :: nu                        !> 1 / reynolds number
    real(kind=rk), dimension(:,:,:), intent(out)   :: dissipation
   ! real(kind=rk), intent(inout)                   :: dissipation_sum
    character(len=*), intent(in)                   :: discretization
    integer(kind=ik), intent(in)                   :: g                         !> grid parameters
    integer(kind=ik), dimension(3), intent(in)     :: Bs

    !> derivatives
    real(kind=rk)                                  :: u_dx, u_dy, u_dz, v_dx, v_dy, v_dz, w_dx, w_dy, w_dz
    real(kind=rk)                                  :: dx_inv, dy_inv, dz_inv   !> inverse of dx, dy, dz, re
    integer(kind=ik)                               :: ix, iy, iz                ! loop variables
    integer(kind=ik)                               :: mxx, mxy, mxz, mnx, mny, mnz, rank     !> sponge lims in # cells

    ! coefficients for Tam&Webb (4th order 1st derivative)
    real(kind=rk), parameter :: a(-3:3) = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, &
    0.02651995_rk/)
    ! coefficients for a standard centered 4th order 1st derivative
    real(kind=rk), parameter :: a_FD4(-2:2) = (/1.0_rk/12.0_rk, -2.0_rk/3.0_rk, 0.0_rk, +2.0_rk/3.0_rk, -1.0_rk/12.0_rk/)


    dissipation = 0.0_rk

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)

    if (size(u,3)>2) then ! 3D case
        dz_inv = 1.0_rk / dx(3)

        select case(discretization)
        case("FD_2nd_central")
            do ix = g+1, Bs(1)+g
                do iy = g+1, Bs(2)+g
                    do iz = g+1, Bs(3)+g
                        u_dx = (u(ix-1,iy,iz)-u(ix-1,iy,iz))*dx_inv*0.5_rk
                        u_dy = (u(ix,iy+1,iz)-u(ix,iy-1,iz))*dy_inv*0.5_rk
                        u_dz = (u(ix,iy,iz+1)-u(ix,iy,iz-1))*dz_inv*0.5_rk
                        
                        v_dx = (v(ix+1,iy,iz)-v(ix-1,iy,iz))*dx_inv*0.5_rk
                        v_dy = (v(ix,iy+1,iz)-v(ix,iy-1,iz))*dy_inv*0.5_rk
                        v_dz = (v(ix,iy,iz+1)-v(ix,iy,iz-1))*dz_inv*0.5_rk
                        
                        w_dx = (w(ix+1,iy,iz)-w(ix-1,iy,iz))*dx_inv*0.5_rk
                        w_dy = (w(ix,iy+1,iz)-w(ix,iy-1,iz))*dy_inv*0.5_rk
                        w_dz = (w(ix,iy,iz+1)-w(ix,iy,iz-1))*dz_inv*0.5_rk

                        dissipation(ix,iy,iz) = (2.d0*(u_dx**2 + v_dy**2 + w_dz**2) &
                                                + u_dy**2+u_dz**2 &
                                                + v_dx**2+v_dz**2 &
                                                + w_dx**2+w_dy**2 &
                                                + 2.d0*(v_dx*u_dy + w_dx*u_dz + w_dy*v_dz))*nu                        
                    end do
                end do
            end do

        case("FD_4th_central")
            ! Note: a_FD4(0) does NOT appear (it is zero...)
            do ix = g+1, Bs(1)+g
                do iy = g+1, Bs(2)+g
                    do iz = g+1, Bs(3)+g
 
                        u_dx = (a_FD4(-2)*u(ix-2,iy,iz) + a_FD4(-1)*u(ix-1,iy,iz) + a_FD4(+1)*u(ix+1,iy,iz) + a_FD4(+2)*u(ix+2,iy,iz))*dx_inv
                        u_dy = (a_FD4(-2)*u(ix,iy-2,iz) + a_FD4(-1)*u(ix,iy-1,iz) + a_FD4(+1)*u(ix,iy+1,iz) + a_FD4(+2)*u(ix,iy+2,iz))*dy_inv
                        u_dz = (a_FD4(-2)*u(ix,iy,iz-2) + a_FD4(-1)*u(ix,iy,iz-1) + a_FD4(+1)*u(ix,iy,iz+1) + a_FD4(+2)*u(ix,iy,iz+2))*dz_inv

                        v_dx = (a_FD4(-2)*v(ix-2,iy,iz) + a_FD4(-1)*v(ix-1,iy,iz) + a_FD4(+1)*v(ix+1,iy,iz) + a_FD4(+2)*v(ix+2,iy,iz))*dx_inv
                        v_dy = (a_FD4(-2)*v(ix,iy-2,iz) + a_FD4(-1)*v(ix,iy-1,iz) + a_FD4(+1)*v(ix,iy+1,iz) + a_FD4(+2)*v(ix,iy+2,iz))*dy_inv
                        v_dz = (a_FD4(-2)*v(ix,iy,iz-2) + a_FD4(-1)*v(ix,iy,iz-1) + a_FD4(+1)*v(ix,iy,iz+1) + a_FD4(+2)*v(ix,iy,iz+2))*dz_inv

                        w_dx = (a_FD4(-2)*w(ix-2,iy,iz) + a_FD4(-1)*w(ix-1,iy,iz) + a_FD4(+1)*w(ix+1,iy,iz) + a_FD4(+2)*w(ix+2,iy,iz))*dx_inv
                        w_dy = (a_FD4(-2)*w(ix,iy-2,iz) + a_FD4(-1)*w(ix,iy-1,iz) + a_FD4(+1)*w(ix,iy+1,iz) + a_FD4(+2)*w(ix,iy+2,iz))*dy_inv
                        w_dz = (a_FD4(-2)*w(ix,iy,iz-2) + a_FD4(-1)*w(ix,iy,iz-1) + a_FD4(+1)*w(ix,iy,iz+1) + a_FD4(+2)*w(ix,iy,iz+2))*dz_inv
                       
                        dissipation(ix,iy,iz) = (2.d0*(u_dx**2 + v_dy**2 + w_dz**2) &
                                                + u_dy**2+u_dz**2 &
                                                + v_dx**2+v_dz**2 &
                                                + w_dx**2+w_dy**2 &
                                                + 2.d0*(v_dx*u_dy + w_dx*u_dz + w_dy*v_dz))*nu
                    end do
                end do
            end do

        case("FD_4th_central_optimized")
            ! Note: a(0) does NOT appear (it is zero...)
            do ix = g+1, Bs(1)+g
                do iy = g+1, Bs(2)+g
                    do iz = g+1, Bs(3)+g

                        u_dx = (a(-3)*u(ix-3,iy,iz) + a(-2)*u(ix-2,iy,iz) &
                        + a(-1)*u(ix-1,iy,iz) + a(+1)*u(ix+1,iy,iz) + a(+2)*u(ix+2,iy,iz) &
                        + a(+3)*u(ix+3,iy,iz))*dx_inv

                        u_dy = (a(-3)*u(ix,iy-3,iz) + a(-2)*u(ix,iy-2,iz) &
                        + a(-1)*u(ix,iy-1,iz) + a(+1)*u(ix,iy+1,iz) + a(+2)*u(ix,iy+2,iz) &
                        + a(+3)*u(ix,iy+3,iz))*dy_inv

                        u_dz = (a(-3)*u(ix,iy,iz-3) + a(-2)*u(ix,iy,iz-2) &
                        + a(-1)*u(ix,iy,iz-1) + a(+1)*u(ix,iy,iz+1) + a(+2)*u(ix,iy,iz+2) &
                        + a(+3)*u(ix,iy,iz+3))*dz_inv

                        v_dx = (a(-3)*v(ix-3,iy,iz) + a(-2)*v(ix-2,iy,iz) &
                        + a(-1)*v(ix-1,iy,iz) + a(+1)*v(ix+1,iy,iz) + a(+2)*v(ix+2,iy,iz) &
                        + a(+3)*v(ix+3,iy,iz))*dx_inv

                        v_dy = (a(-3)*v(ix,iy-3,iz) + a(-2)*v(ix,iy-2,iz) &
                        + a(-1)*v(ix,iy-1,iz) + a(+1)*v(ix,iy+1,iz) + a(+2)*v(ix,iy+2,iz) &
                        + a(+3)*v(ix,iy+3,iz))*dy_inv

                        v_dz = (a(-3)*v(ix,iy,iz-3) + a(-2)*v(ix,iy,iz-2) &
                        + a(-1)*v(ix,iy,iz-1) + a(+1)*v(ix,iy,iz+1) + a(+2)*v(ix,iy,iz+2) &
                        + a(+3)*v(ix,iy,iz+3))*dz_inv

                        w_dx = (a(-3)*w(ix-3,iy,iz) + a(-2)*w(ix-2,iy,iz) &
                        + a(-1)*w(ix-1,iy,iz) + a(+1)*w(ix+1,iy,iz) + a(+2)*w(ix+2,iy,iz) &
                        + a(+3)*w(ix+3,iy,iz))*dx_inv

                        w_dy = (a(-3)*w(ix,iy-3,iz) + a(-2)*w(ix,iy-2,iz) &
                        + a(-1)*w(ix,iy-1,iz) + a(+1)*w(ix,iy+1,iz) + a(+2)*w(ix,iy+2,iz) &
                        + a(+3)*w(ix,iy+3,iz))*dy_inv

                        w_dz = (a(-3)*w(ix,iy,iz-3) + a(-2)*w(ix,iy,iz-2) &
                        + a(-1)*w(ix,iy,iz-1) + a(+1)*w(ix,iy,iz+1) + a(+2)*w(ix,iy,iz+2) &
                        + a(+3)*w(ix,iy,iz+3))*dz_inv

                        dissipation(ix,iy,iz) = (2.d0*(u_dx**2 + v_dy**2 + w_dz**2) &
                                                + u_dy**2+u_dz**2 &
                                                + v_dx**2+v_dz**2 &
                                                + w_dx**2+w_dy**2 &
                                                + 2.d0*(v_dx*u_dy + w_dx*u_dz + w_dy*v_dz))*nu      
                    end do
                end do
            end do
        case default
            write(*,*) discretization
            write(*,*) "abort command commented out"
            !call abort(240321, "ERROR: discretization method in discretization is unknown")
        end select

    else
        ! 2D case
        select case(discretization)
        case("FD_2nd_central")
            do ix = g+1, Bs(1)+g
                do iy = g+1, Bs(2)+g
                    u_dy = (u(ix,iy+1,1)-u(ix,iy-1,1))*dy_inv*0.5_rk
                    v_dx = (v(ix+1,iy,1)-v(ix-1,iy,1))*dx_inv*0.5_rk

                    dissipation(ix,iy,1) = (u_dy*u_dy + 2*u_dy*v_dx + v_dx*v_dx)*0.5_rk*nu
                end do
            end do

        case("FD_4th_central")
            do ix = g+1, Bs(1)+g
                do iy = g+1, Bs(2)+g
                    u_dy = (a_FD4(-2)*u(ix,iy-2,1) + a_FD4(-1)*u(ix,iy-1,1) + a_FD4(+1)*u(ix,iy+1,1) &
                    + a_FD4(+2)*u(ix,iy+2,1))*dy_inv
                    v_dx = (a_FD4(-2)*v(ix-2,iy,1) + a_FD4(-1)*v(ix-1,iy,1) + a_FD4(+1)*v(ix+1,iy,1) &
                    + a_FD4(+2)*v(ix+2,iy,1))*dx_inv

                    dissipation(ix,iy,1) = (u_dy*u_dy + 2*u_dy*v_dx + v_dx*v_dx)*0.5_rk*nu
                end do
            end do

        case("FD_4th_central_optimized")
            do ix = g+1, Bs(1)+g
                do iy = g+1, Bs(2)+g
                    u_dy = (a(-3)*u(ix,iy-3,1) + a(-2)*u(ix,iy-2,1) + a(-1)*u(ix,iy-1,1) + a(+1)*u(ix,iy+1,1) &
                    + a(+2)*u(ix,iy+2,1) + a(+3)*u(ix,iy+3,1))*dy_inv
                    v_dx = (a(-3)*v(ix-3,iy,1) + a(-2)*v(ix-2,iy,1) + a(-1)*v(ix-1,iy,1) + a(+1)*v(ix+1,iy,1) &
                    + a(+2)*v(ix+2,iy,1) + a(+3)*v(ix+3,iy,1))*dx_inv

                    dissipation(ix,iy,1) = (u_dy*u_dy + 2*u_dy*v_dx + v_dx*v_dx)*0.5_rk*nu
                end do
            end do
        case default
            write(*,*) discretization
            write(*,*) "abort command commented out"
            ! weird compilaiton issue here
            !call abort(2403211, "ERROR: discretization method in discretization is unknown")
        end select
    end if

end subroutine compute_energy_dissipation
