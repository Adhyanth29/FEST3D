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
         subroutine setup_multiphase_scheme(F_surface, control, dims)
            !< sets up the mutiphase schemes
            implicit none
            type(controltype), intent(in) :: control
            !< Control parameters
            type(extent), intent(in) :: dims
            !< extent of the 3D domain
            real(wp), dimension(:, :, :, :), allocatable, intent(inout) :: F_surface
            !< Stores the surface tension force with 3 dimensions

            imx = dims%imx
            jmx = dims%jmx
            kmx = dims%kmx

            n_var = control%n_var

            !call setup_interpolant_scheme(dims)

            call alloc(F_surface, 1, imx, 1, jmx, 1, kmx, 1, 3, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'Surface Tension Force - Multiphase')
         end subroutine setup_multiphase_scheme


         subroutine perform_multiphase(F_surface, delta_t, Ifaces, Jfaces, Kfaces, scheme, flow, dims, qp, nodes)
            !< Authorises the type of multiphase operation
            implicit none
            type(schemetype), intent(in) :: scheme
            !< finite-volume Schemes
            !type(flowtype), intent(in) :: flow
            !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
            type(extent), intent(in) :: dims
            !< extent of the 3D domain
            real(wp), dimension(:, :, :, :), intent(inout) :: F_surface
            !< Stores the surface tension of the domain
            type(nodetype), dimension(-2:dims%imx+3,-2:dims%jmx+3,-2:dims%kmx+3), intent(in) :: nodes
            !< Grid points
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            !< Stores cell parameter: volume
            real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(inout), target :: qp
            type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
            !< Store face quantites for I faces 
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
            !< Store face quantites for J faces 
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
            !< Store face quantites for K faces
            real(wp) , dimension(1:dims%imx-1, 1:dims%jmx-1, 1:dims%kmx-1), intent(in) :: delta_t
            !< Local time increment value at each cell center
            real(wp) :: sigma, epsilon

            select case (trim(scheme%multiphase))
                case ("clsvof")
                  !< For the incompressible Coupled Level-Set Volume of Fluid case
                  sigma = flow%sigma
                  epsilon = flow%epsilon
                  call perform_clsvof_incomp(dims, nodes, cells, Ifaces, Jfaces, sigma, epsilon, F_surface, delta_t)
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