!< Multiphase calculation
module multi_phase
!< Multiphase calculation

#include "../../../debug.h"
#include "../../../error.h"
   use vartypes
   use utils,         only: alloc
   use clsvof_incomp, only: perform_clsvof_incomp => perform_multiphase
   implicit none
   integer :: imx, jmx, kmx, n_var
   private

   !< Public Members
   public :: setup_multiphase_scheme
   public :: perform_multiphase
   public :: compute_residue

   contains
         subroutine setup_multiphase_scheme(residue, F_x, F_y, F_z, control, dims)
            !< sets up the mutiphase schemes
            implicit none
            type(controltype), intent(in) :: control
            !< Control parameters
            type(extent), intent(in) :: dims
            !< extent of the 3D domain
            real(wp), dimension(:, :, :, :), allocatable, intent(out), target :: residue
            !< Store residue at each cell-center
            real(wp), dimension(:, :, :, :), allocatable, intent(out) :: F_x
            !< Store fluxes throught the I faces
            real(wp), dimension(:, :, :, :), allocatable, intent(out) :: F_y
            !< Store fluxes throught the J faces
            real(wp), dimension(:, :, :, :), allocatable, intent(out) :: F_z
            !< Store fluxes throught the K faces

            imx = dims%imx
            jmx = dims%jmx
            kmx = dims%kmx

            n_var = control%n_var

            !call setup_interpolant_scheme(dims)

            call alloc(F_x, 1, imx, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'F_x - Scheme')
            call alloc(F_y, 1, imx-1, 1, jmx, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'F_y - Scheme')
            call alloc(F_z, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'F_z - Scheme')
            call alloc(residue, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'residue - Scheme')
         end subroutine setup_multiphase_scheme

         subroutine perform_multiphase(F_surface, Ifaces, Jfaces, Kfaces, scheme, flow, bc, dims, qp, nodes, cells, epsilon, sigma, vof)
            !< Authorises the type of multiphase operation
            implicit none
            type(schemetype), intent(in) :: scheme
            !< finite-volume Schemes
            !type(flowtype), intent(in) :: flow
            !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
            type(extent), intent(in) :: dims
            !< extent of the 3D domain
            real(wp), dimension(:, :, :, :), intent(inout) :: F_x
            !< Store fluxes throught the I faces
            real(wp), dimension(:, :, :, :), intent(inout) :: F_y
            !< Store fluxes throught the J faces
            real(wp), dimension(:, :, :, :), intent(inout) :: F_z
            !< Store fluxes throught the K faces
            type(nodetype), dimension(-2:dims%imx+3,-2:dims%jmx+3,-2:dims%kmx+3), intent(in) :: nodes
            !< Grid points
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            !< Stores cell parameter: volume
            real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(inout), target :: qp
            real(wp), intent(in) :: epsilon
            !< Numerical interface width
                        real(wp), intent(in) :: sigma
            !< Surface tension at interface - fluid property
            type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
            !< Store face quantites for I faces 
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
            !< Store face quantites for J faces 
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
            !< Store face quantites for K faces 
            type(boundarytype), intent(in) :: bc
            !< boundary conditions and fixed values

            select case (trim(scheme%multiphase))
                case ("clsvof")
                  !< For the incompressible Coupled Level-Set Volume of Fluid case
                  call perform_clsvof_incomp(dims, nodes, cells, Ifaces, Jfaces, vof, qp, sigma, epsilon, flow, F_surface)
                  !call compute_fluxes_van_leer(F_x,F_y,F_z,x_qp_left,x_qp_right,y_qp_left,y_qp_right,z_qp_left,z_qp_right,Ifaces,Jfaces,Kfaces,flow,bc,dims)
                case ("clsvof_c")
                  !< For the compressible Coupled Level-Set Volume of Fluid case
                  !call compute_fluxes_ausm(F_x,F_y,F_z,x_qp_left,x_qp_right,y_qp_left,y_qp_right,z_qp_left,!z_qp_right,Ifaces,Jfaces,Kfaces,flow,bc,dims)
            
                  ! To do
                  continue
                case ("dpm")
                  !< For the Dispersed Phase Model case
                  !call compute_fluxes_ausmP(F_x,F_y,F_z,x_qp_left,x_qp_right,y_qp_left,y_qp_right,!z_qp_left,z_qp_right,Ifaces,Jfaces,Kfaces,flow,bc,dims)

                  ! To do
                  continue
                case default
                    Fatal_error
            end select
            
         end subroutine perform_multiphase

         subroutine compute_residue(residue,F_x,F_y,F_z,dims)
            
            implicit none
            type(extent), intent(in) :: dims
            !< extent of the 3D domain
            real(wp), dimension(:, :, :, :), intent(out)  :: residue
            !< Store residue at each cell-center
            real(wp), dimension(:, :, :, :), intent(in) :: F_x
            !< Store fluxes throught the I faces
            real(wp), dimension(:, :, :, :), intent(in) :: F_y
            !< Store fluxes throught the J faces
            real(wp), dimension(:, :, :, :), intent(in) :: F_z
            !< Store fluxes throught the K faces
            
            integer :: i, j, k, l

            DebugCall('compute_residue')

            do l = 1, dims%n_var
             do k = 1, dims%kmx - 1
              do j = 1, dims%jmx - 1
               do i = 1, dims%imx - 1
               ! residue(i, j, k, l) = (F_x(i+1, j, k, l) - F_x(i, j, k, l)) &
               !                     + (F_y(i, j+1, k, l) - F_y(i, j, k, l)) &
               !                     + (F_z(i, j, k+1, l) - F_z(i, j, k, l))
               end do
              end do
             end do
            end do
        
        end subroutine compute_residue