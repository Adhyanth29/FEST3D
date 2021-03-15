   !< Incompressible Two Phase Model - CLSVOF
module clsvof_incomp
   !<
   !< Reference: Yeng-Yung Tsui, Cheng-Yen Liu & Shi-Wen Lin "Coupled 
   !< level-set and volume-of-fluid method for two-phase flow 
   !< calculations", Numerical Heat Transfer, Part B: Fundamentals,
   !< 71:2, 173-185, 2017, DOI: 10.1080/10407790.2016.1265311
   !<
   !   -- Adhyanth Giri Ajay
   !   -- First build

#include "../error.h"
#include "../debug.h"

   use vartypes
   use gradients, only : compute_gradient_G

   implicit none
   
   private
   real(wp), parameter :: pi
   !< An accurate definition of pi
   real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2) :: Fx
   !< Surface tension force to be calcualted - X component
   real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2) :: Fy
   !< Surface tension force to be calcualted - Y component
   real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2) :: Fz
   !< Surface tension force to be calcualted - Z component
   type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2) :: del_h
   !< Stores Cell size of each cell (approximated)

   pi = 4.0*atan(1.0)

   public :: pi

   contains
   
         subroutine perform_clsvof_incomp()
            !< Performs the overall computation of the CLSVOF algorithm
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            type(nodetype), dimension(-2:dims%imx+3,-2:dims%jmx+3,-2:dims%kmx+3), intent(out) :: nodes
            !< Grid points
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            !< Stores cell parameter: volume
            type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
            !< Input varaible which stores I faces' area and unit normal
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
            !< Input varaible which stores J faces' area and unit normal
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
            !< Input variable which stores K faces' area and unit normal
            real(wp), intent(in) :: del_t
            !< Time step
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(out) :: vof_n
            !< Output the next time-step of volume fraction
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: vof_o
            !< Storing the previous time-step volume fraction
            real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in), Target :: qp
            real(wp), dimension(:, :, :), pointer :: x_speed      
            !< U pointer, point to slice of qp (:,:,:,2)
            real(wp), dimension(:, :, :), pointer :: y_speed      
            !< V pointer, point to slice of qp (:,:,:,3)
            real(wp), dimension(:, :, :), pointer :: z_speed      
            !< W pointer, point to slice of qp (:,:,:,4)
            !< Store primitive variable at cell center
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: phi
            !< New Level-Set function
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: phi_init
            !< New Level-Set function
            real(wp), dimension(0:dims%imx,0:dims%jmx,0:dims%kmx), intent(out) :: grad_phi_x
            !< Stores value of level-set gradient
            real(wp), dimension(0:dims%imx,0:dims%jmx,0:dims%kmx), intent(out) :: grad_phi_y
            !< Stores value of level-set gradient
            real(wp), dimension(0:dims%imx,0:dims%jmx,0:dims%kmx), intent(out) :: grad_phi_z
            !< Stores value of level-set gradient
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(out) :: d_delta
            !< Output variable storing the Dirac Delta smooth function
            real(wp), dimension(:,:,:), allocatable, intent(out) :: K
            !< Output variable storing the gradient of curvature
            real(wp), intent(in) :: sigma
            !< Surface tension at interface - fluid property
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(out) :: H
            !< Output variable storing the Heaviside values
            real(wp), intent(in) :: epsilon
            !< Numerical interface width

            !!!!< Parameters for smoothen subroutine need to be entered

            !< Pointer allocation
            x_speed(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2) => qp(:, :, :, 2)
            y_speed(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2) => qp(:, :, :, 3)
            z_speed(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2) => qp(:, :, :, 4)

            call cell_size(cells, dims)
            !< Finding the approximate cell size
            call vof_adv(vof_n, vof_o, qp, cells, Ifaces, Jfaces, Kfaces, del_t, nodes, dims)
            !< Performs vof advection to find new timestep vof
            call level_set_coupling(phi_init, vol_n, dims)
            !< Coupled the Level-set function with new VOF value
            call level_set_advancement(phi, phi_init, grad_phi_x, grad_phi_y, grad_phi_z, &
            del_tau, cells, Ifaces, Jfaces, Kfaces, dims)
            !< Performs time advancement of Level set
            call dirac_delta(d_delta, phi, epsilon, cells, dims)
            !< Finds dirac delta function using new LS function
            call curvature(K, grad_phi_x, grad_phi_y, grad_phi_z, cells, Ifaces, Jfaces, Kfaces, dims)
            !< Finds curvature of the new level set function
            call surface_tension_force(sigma, K, d_delta, grad_phi_x, grad_phi_y, grad_phi_z, dims)
            !< Finds the surface tension force using K, Dirac Delta
            !< and the gradient of level set value
            call heaviside(H, phi, epsilon, cells, dims)
            !< computes heaviside funciton to use for smoothening
            call smoothen_G(density)
            !< Smoothens density using a heaviside formulation
            call smoothen_G(viscosity)
            !< Smoothens viscosity using a heaviside formulation
         end subroutine perform_clsvof_incomp

         subroutine cell_size(cells, dims)
            !< to find the cell size required for this module
            implicit none 
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            !< Stores cell parameter: volume
            del_h(:,:,:) = cells(:,:,:)%volume**(1.0/3.0)
         end subroutine cell_size


         subroutine vof_adv(vof_n, vof_o, qp, cells, Ifaces, Jfaces, Kfaces, del_t, nodes, dims)
            !< to account for the volume fraction advection VOF
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            real(wp), intent(in) :: del_t
            !< Time step
            type(nodetype), dimension(-2:dims%imx+3,-2:dims%jmx+3,-2:dims%kmx+3), intent(out) :: nodes
            !< Grid points
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(out) :: vof_n
            !< Output the next time-step of volume fraction
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: vof_o
            !< Storing the previous time-step volume fraction
            real(wp), dimension(:, :, :), pointer :: x_speed      
            !< U pointer, point to slice of qp (:,:,:,2)
            real(wp), dimension(:, :, :), pointer :: y_speed      
            !< V pointer, point to slice of qp (:,:,:,3)
            real(wp), dimension(:, :, :), pointer :: z_speed      
            !< W pointer, point to slice of qp (:,:,:,4)
            real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in), Target :: qp
            !< Store primitive variable at cell center
            type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
            !< Input varaible which stores I faces' area and unit normal
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
            !< Input varaible which stores J faces' area and unit normal
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
            !< Input variable which stores K faces' area and unit normal
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            !< Stores cell parameter: volume
            real(wp), dimension(:,:,:), allocatable :: s
            !< Stores the sum of the product of face velocity with wetted area
            integer :: i,j,k

            !< Pointer allocation
            x_speed(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2) => qp(:, :, :, 2)
            y_speed(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2) => qp(:, :, :, 3)
            z_speed(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2) => qp(:, :, :, 4)

            !< Calling interface reconstruction to find wetted area
            call interface_reconstruction(vof_o, cells, nodes, dims)

            !!!!!< ASSUMING FACE STATES AS CELL CENTERS (NEED TO PERFORM MUSCL SCHEME FOR PROPER VALUES)
            do k=0,dims%kmx
               do j=0,dims%jmx
                  do i=0,dims%imx
                     s(i,j,k) = x_speed(i,j,k)*Ifacewet(i,j,k)*Ifaces(i,j,k)%nx&
                                + x_speed(i+1,j,k)*Ifacewet(i+1,j,k)*Ifaces(i+1,j,k)%nx&
                                + y_speed(i,j,k)*Jfacewet(i,j,k)*Jfaces(i,j,k)%ny&
                                + y_speed(i,j+1,k)*Jfacewet(i,j+1,k)*Jfaces(i,j+1,k)%ny&
                                + z_speed(i,j,k)*Kfacewet(i,j,k)*Kfaces(i,j,k)%nz&
                                + z_speed(i,j,k+1)*Kfacewet(i,j,k+1)*Kfaces(i,j,k+1)%nz
                  end do
               end do
            end do
            vof_n(:,:,:) = vof_o(:,:,:) - del_t/cells(:,:,:)%volume*s(:,:,:)
            vof_o(:,:,:) = vof_n(:,:,:)
            
            !< Calling VoF correction for the two filling and two depletion cases
            !< Calling it twice to remove any induced filling/depletion cases
            !< from the first time
            call vof_correction(vof_n, x_speed, y_speed, z_speed, Ifaces, Jfaces, Kfaces, cells, dims)
            call vof_correction(vof_n, x_speed, y_speed, z_speed, Ifaces, Jfaces, Kfaces, cells, dims)
            !!! NEED TO ACCOUNT FOR THE INDEXES IN LOOPS

            ! This will NOT work as of now. Need to account for the bondary grad values to make sure the volume fractions are calculated accordingly for the full domain
         end subroutine vof_adv


         subroutine vof_correction(vof, x_speed, y_speed, z_speed, Ifaces, Jfaces, Kfaces, cells, dims)
            !< Corrects the vof to account for the four cases:
            !< Over-depletion, Over-filling, Under-depletion and Under-filling
            implicit none
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(inout) :: vof
            !< Output the next time-step of volume fraction
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
            !< Input varaible which stores I faces' area and unit normal
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
            !< Input varaible which stores J faces' area and unit normal
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
            !< Input variable which stores K faces' area and unit normal
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            !< Stores cell parameter: volume
            real(wp), dimension(:, :, :), pointer, intent(in) :: x_speed      
            !< U pointer, point to slice of qp (:,:,:,2) 
            real(wp), dimension(:, :, :), pointer, intent(in) :: y_speed      
            !< V pointer, point to slice of qp (:,:,:,3) 
            real(wp), dimension(:, :, :), pointer, intent(in) :: z_speed      
            !< W pointer, point to slice of qp (:,:,:,4)
            real(wp), dimension(-2:dims%imx+3,-2:dims%jmx+3,-2:dims%kmx+3) :: vof_node
            !< Stores the values of vof at the nodes
            real(wp), dimension(6) :: c, w
            real(wp) :: sum = 0, w_sum = 0
            !< Stores the correction weights for I, J, and K face directions and sum
            real(wp), dimension(:), allocatable :: c
            !< Temp calculation value plus vof node storing variable
            real(wp), dimension(:), allocatable, pointer :: vof_no
            !< Points to the values of vof_node that would be calculated
            integer :: i,j,k,m

            !< To find the vof value at the nodes using adjacent cell centers
            do k = 0,dims%kmx+1
               do j = 0,dims%jmx+1
                  do i = 0,dims%imx+1
                     !< Calculating weights using inverse volume
                     w_sum = 1.0/cells(i-1,j-1,k-1)%volume + 1.0/cells(i,j-1,k-1)%volume + &
                           1.0/cells(i-1,j,k-1)%volume + 1.0/cells(i,j,k-1)%volume + &
                           1.0/cells(i-1,j-1,k)%volume + 1.0/cells(i,j-1,k)%volume + &
                           1.0/cells(i-1,j,k)%volume + 1.0/cells(i,j,k)%volume

                     !< sum of local weight*vof / sum of local weights
                     vof_node(i,j,k) = (vof(i-1,j-1,k-1)/cells(i-1,j-1,k-1)%volume + &
                                       vof(i,j-1,k-1)/cells(i,j-1,k-1)%volume + &
                                       vof(i-1,j,k-1)/cells(i-1,j,k-1)%volume + &
                                       vof(i,j,k-1)/cells(i,j,k-1)%volume + &
                                       vof(i-1,j-1,k)/cells(i-1,j-1,k)%volume + &
                                       vof(i,j-1,k)/cells(i,j-1,k)%volume + &
                                       vof(i-1,j,k)/cells(i-1,j,k)%volume + &
                                       vof(i,j,k)/cells(i,j,k)%volume)/w_sum
                  end do
               end do
            end do

            do k = 0,dims%kmax
               do j = 0,dims%jmx
                  do i = 0,dims%imx
                     c(1) = x_speed(i,j,k)*Ifaces(i,j,k)*Ifaces(i,j,k)%nx,0
                     c(2) = x_speed(i+1,j,k)*Ifaces(i+1,j,k)*Ifaces(i+1,j,k)%nx
                     c(3) = y_speed(i,j,k)*Jfaces(i,j,k)*Jfaces(i,j,k)%ny
                     c(4) = y_speed(i,j+1,k)*Jfaces(i,j+1,k)*Jfaces(i,j+1,k)%ny
                     c(5) = z_speed(i,j,k)*Kfaces(i,j,k)*Kfaces(i,j,k)%nz
                     c(6) = z_speed(i,j,k+1)*Kfaces(i,j,k+1)*Kfaces(i,j,k+1)%nz
                     do m = 1:6
                        sum = sum + max(c(m),0)
                     end do
                     do m =1:6
                        w(m) = max(c(m),0)/sum
                     end do

                     if (vof(i,j,k) > 1.0) then
                        !< Over-filling
                        vof(i-1,j,k) = vof(i-1,j,k) + w(1)*(vof(i,j,k)-1)* &
                                       cells(i,j,k)%volume/cells(i-1,j,k)%volume
                        vof(i+1,j,k) = vof(i+1,j,k) + w(2)*(vof(i,j,k)-1)* &
                                       cells(i,j,k)%volume/cells(i+1,j,k)%volume
                        vof(i,j-1,k) = vof(i,j-1,k) + w(3)*(vof(i,j,k)-1)* &
                                       cells(i,j,k)%volume/cells(i,j-1,k)%volume
                        vof(i,j+1,k) = vof(i,j+1,k) + w(4)*(vof(i,j,k)-1)* &
                                       cells(i,j,k)%volume/cells(i,j+1,k)%volume
                        vof(i,j,k-1) = vof(i,j,k-1) + w(5)*(vof(i,j,k)-1)* &
                                       cells(i,j,k)%volume/cells(i,j,k-1)%volume
                        vof(i,j,k+1) = vof(i,j,k+1) + w(6)*(vof(i,j,k)-1)* &
                                       cells(i,j,k)%volume/cells(i,j,k+1)%volume
                        vof(i,j,k)   = 1.0
                     else if (vof(i,j,k) < 0.0) then
                        !< Over-depleting
                        vof(i-1,j,k) = vof(i-1,j,k) + w(1)*vof(i,j,k)* &
                                       cells(i,j,k)%volume/cells(i-1,j,k)%volume
                        vof(i+1,j,k) = vof(i+1,j,k) + w(2)*vof(i,j,k)* &
                                       cells(i,j,k)%volume/cells(i+1,j,k)%volume
                        vof(i,j-1,k) = vof(i,j-1,k) + w(3)*vof(i,j,k)* &
                                       cells(i,j,k)%volume/cells(i,j-1,k)%volume
                        vof(i,j+1,k) = vof(i,j+1,k) + w(4)*vof(i,j,k)* &
                                       cells(i,j,k)%volume/cells(i,j+1,k)%volume
                        vof(i,j,k-1) = vof(i,j,k-1) + w(5)*vof(i,j,k)* &
                                       cells(i,j,k)%volume/cells(i,j,k-1)%volume
                        vof(i,j,k+1) = vof(i,j,k+1) + w(6)*vof(i,j,k)* &
                                       cells(i,j,k)%volume/cells(i,j,k+1)%volume
                        vof(i,j,k)   = 0.0
                     else
                        vof_no(1) => vof_node(i,j,k)
                        vof_no(2) => vof_node(i+1,j,k)
                        vof_no(3) => vof_node(i,j+1,k)
                        vof_no(4) => vof_node(i+1,j+1,k)
                        vof_no(5) => vof_node(i,j,k+1)
                        vof_no(6) => vof_node(i+1,j,k+1)
                        vof_no(7) => vof_node(i,j+1,k+1)
                        vof_no(8) => vof_node(i+1,j+1,k+1)
                        if (vof(i,j,k) < 1.0 .and. vof_n(:)>0.5) then
                           !< Under-filling
                           vof(i-1,j,k) = vof(i-1,j,k) + w(1)*(vof(i,j,k)-1)* &
                                          cells(i,j,k)%volume/cells(i-1,j,k)%volume
                           vof(i+1,j,k) = vof(i+1,j,k) + w(2)*(vof(i,j,k)-1)* &
                                          cells(i,j,k)%volume/cells(i+1,j,k)%volume
                           vof(i,j-1,k) = vof(i,j-1,k) + w(3)*(vof(i,j,k)-1)* &
                                          cells(i,j,k)%volume/cells(i,j-1,k)%volume
                           vof(i,j+1,k) = vof(i,j+1,k) + w(4)*(vof(i,j,k)-1)* &
                                          cells(i,j,k)%volume/cells(i,j+1,k)%volume
                           vof(i,j,k-1) = vof(i,j,k-1) + w(5)*(vof(i,j,k)-1)* &
                                          cells(i,j,k)%volume/cells(i,j,k-1)%volume
                           vof(i,j,k+1) = vof(i,j,k+1) + w(6)*(vof(i,j,k)-1)* &
                                          cells(i,j,k)%volume/cells(i,j,k+1)%volume
                           vof(i,j,k)   = 1.0
                        else if (vof(i,j,k) > 0.0 .and. vof_n(:)<0.5) then
                           !< Under-depletion
                           vof(i-1,j,k) = vof(i-1,j,k) + w(1)*vof(i,j,k)* &
                                          cells(i,j,k)%volume/cells(i-1,j,k)%volume
                           vof(i+1,j,k) = vof(i+1,j,k) + w(2)*vof(i,j,k)* &
                                          cells(i,j,k)%volume/cells(i+1,j,k)%volume
                           vof(i,j-1,k) = vof(i,j-1,k) + w(3)*vof(i,j,k)* &
                                          cells(i,j,k)%volume/cells(i,j-1,k)%volume
                           vof(i,j+1,k) = vof(i,j+1,k) + w(4)*vof(i,j,k)* &
                                          cells(i,j,k)%volume/cells(i,j+1,k)%volume
                           vof(i,j,k-1) = vof(i,j,k-1) + w(5)*vof(i,j,k)* &
                                          cells(i,j,k)%volume/cells(i,j,k-1)%volume
                           vof(i,j,k+1) = vof(i,j,k+1) + w(6)*vof(i,j,k)* &
                                          cells(i,j,k)%volume/cells(i,j,k+1)%volume
                           vof(i,j,k)   = 0.0
                        end if
                     end if
                  end do
               end do
            end do

         end subroutine vof_correction   


         subroutine interface_reconstruction(vof, cells, nodes, dims)
            !< to reconstruct interface using vof 0.5
            !< Finds values at nodes and interpolates along the 
            !< face to identify where vof 0.5 occurs
            implicit none
            type(nodetype), dimension(-2:dims%imx+3,-2:dims%jmx+3,-2:dims%kmx+3), intent (in) :: nodes
            !< Stores the location of nodes
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            !< Stores cell parameter: volume
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2), intent(in) :: vof
            !< Input variables for vof at cell centers
            real(wp), dimension(-2:dims%imx+3,-2:dims%jmx+3,-2:dims%kmx+3) :: vof_node
            !< Stores the values of vof at the nodes
            type(interfacetype), dimension(:,:,:), allocatable :: inter_x
            !< Stores intercepts on X direction mesh
            type(interfacetype), dimension(:,:,:), allocatable :: inter_y
            !< Stores intercepts on Y direction mesh
            type(interfacetype), dimension(:,:,:), allocatable :: inter_z
            !< Stores intercepts on Z direction mesh
            type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(out) :: Ifacewet
            !< Input varaible which stores I faces' area and unit normal
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(out) :: Jfacewet
            !< Input varaible which stores J faces' area and unit normal
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(out) :: Kfacewet
            real(wp), dimension()
            integer :: i,j,k
            real(wp) :: w_sum
            !< Stores sum of weights of neighbouring cells
            !< Initialising to 0
            inter_x(:,:,:) = 0.0
            inter_y(:,:,:) = 0.0
            inter_z(:,:,:) = 0.0
            Ifacewet(:,:,:) = 0.0
            Jfacewet(:,:,:) = 0.0
            Kfacewet(:,:,:) = 0.0

            ! Sweep to make wetted area as 1 for all cells where VOF = 1
            do k = 0,dims%kmx+1
               do j = 0,dims%jmx+1
                  do i = 0,dims%imx+1
                     if (vof(:,:,:) == 1) then
                        Ifacewet(i,j,k) = Ifaces(i,j,k)%A
                        Ifacewet(i+1,j,k) = Ifaces(i+1,j,k)%A
                        Jfacewet(i,j,k) = Jfaces(i,j,k)%A
                        Jfacewet(i,j+1,k) = Jfaces(i,j+1,k)%A
                        Kfacewet(i,j,k) = Kfaces(i,j,k)%A
                        Kfacewet(i,j,k+1) = Kfaces(i,j,k+1)%A
                     end if
                  end do
               end do
            end do

            !< To find the vof value at the nodes using adjacent cell centers
            do k = 1,dims%kmx+1
               do j = 1,dims%jmx+1
                  do i = 1,dims%imx+1
                     !< Calculating weights using inverse volume
                     w_sum = 1.0/cells(i-1,j-1,k-1)%volume + 1.0/cells(i,j-1,k-1)%volume + &
                             1.0/cells(i-1,j,k-1)%volume + 1.0/cells(i,j,k-1)%volume + &
                             1.0/cells(i-1,j-1,k)%volume + 1.0/cells(i,j-1,k)%volume + &
                             1.0/cells(i-1,j,k)%volume + 1.0/cells(i,j,k)%volume

                     !< sum of local weight*vof / sum of local weights
                     vof_node(i,j,k) = (vof(i-1,j-1,k-1)/cells(i-1,j-1,k-1)%volume + &
                                       vof(i,j-1,k-1)/cells(i,j-1,k-1)%volume + &
                                       vof(i-1,j,k-1)/cells(i-1,j,k-1)%volume + &
                                       vof(i,j,k-1)/cells(i,j,k-1)%volume + &
                                       vof(i-1,j-1,k)/cells(i-1,j-1,k)%volume + &
                                       vof(i,j-1,k)/cells(i,j-1,k)%volume + &
                                       vof(i-1,j,k)/cells(i-1,j,k)%volume + &
                                       vof(i,j,k)/cells(i,j,k)%volume)/w_sum
                  end do
               end do
            end do
            
            !< Correction for avoiding false intercepts
            vof_node(0,:,:) = vof_node(1,:,:)
            vof_node(:,0,:) = vof_node(:,1,:)
            vof_node(:,:,0) = vof_node(:,:,1)
            vof_node(dims%imx+2,:,:) = vof_node(dims%imx+1,:,:)
            vof_node(:,dims%jmx+2,:) = vof_node(:,dims%jmx+1,:)
            vof_node(:,:,dims%kmx+2) = vof_node(:,:,dims%kmx+1)

            !< Linear interpolation to find intercept locations where vof = 0.5
            do k = 1,dims%kmx+1
               do j = 1,dims%jmx+1
                  do i = 1,dims%imx+1
                     if(vof_node(i,j,k) == 0.5) then
                        inter_x(i,j,k)%x = nodes(i,j,k)%x
                        inter_x(i,j,k)%y = nodes(i,j,k)%y
                        inter_x(i,j,k)%z = nodes(i,j,k)%z
                        inter_y(i,j,k)%x = nodes(i,j,k)%x
                        inter_y(i,j,k)%y = nodes(i,j,k)%y
                        inter_y(i,j,k)%z = nodes(i,j,k)%z
                        inter_z(i,j,k)%x = nodes(i,j,k)%x
                        inter_z(i,j,k)%y = nodes(i,j,k)%y
                        inter_z(i,j,k)%z = nodes(i,j,k)%z
                     end if
                     if((vof_node(i,j,k) < 0.5 .and. vof_node(i+1,j,k) > 0.5) .or. &
                        (vof_node(i,j,k) > 0.5 .and. vof_node(i+1,j,k) < 0.5)) then
                           !< Intercepts on face edge in I direcion
                           inter_x(i,j,k)%x = nodes(i,j,k)%x + (nodes(i+1,j,k)%x - nodes(i,j,k)%x)*&
                               (0.5 - vof_node(i,j,k))/(vof_node(i+1,j,k) - vof_node(i,j,k))
                           inter_x(i,j,k)%y = nodes(i,j,k)%y + (nodes(i+1,j,k)%y - nodes(i,j,k)%y)*&
                           (0.5 - vof_node(i,j,k))/(vof_node(i+1,j,k) - vof_node(i,j,k))
                           inter_x(i,j,k)%z = nodes(i,j,k)%z + (nodes(i+1,j,k)%z - nodes(i,j,k)%z)*&
                           (0.5 - vof_node(i,j,k))/(vof_node(i+1,j,k) - vof_node(i,j,k))
                     end if
                     if((vof_node(i,j,k) < 0.5 .and. vof_node(i,j+1,k) > 0.5) .or. &
                        (vof_node(i,j,k) > 0.5 .and. vof_node(i,j+1,k) < 0.5)) then
                           !< Intercepts on face edge in J direcion
                           inter_y(i,j,k)%x = nodes(i,j,k)%x + (nodes(i,j+1,k)%x - nodes(i,j,k)%x)*&
                           (0.5 - vof_node(i,j,k))/(vof_node(i,j+1,k) - vof_node(i,j,k))
                           inter_y(i,j,k)%y = nodes(i,j,k)%y + (nodes(i,j+1,k)%y - nodes(i,j,k)%y)*&
                                (0.5 - vof_node(i,j,k))/(vof_node(i,j+1,k) - vof_node(i,j,k))
                           inter_y(i,j,k)%z= nodes(i,j,k)%z + (nodes(i,j+1,k)%z - nodes(i,j,k)%z)*&
                           (0.5 - vof_node(i,j,k))/(vof_node(i,j+1,k) - vof_node(i,j,k))
                     end if
                     if((vof_node(i,j,k) < 0.5 .and. vof_node(i,j,k+1) > 0.5) .or. &
                        (vof_node(i,j,k) > 0.5 .and. vof_node(i,j,k+1) < 0.5)) then
                           !< Intercepts on face edge in K direcion
                           inter_z(i,j,k)%x = nodes(i,j,k)%x + (nodes(i,j,k+1)%x - nodes(i,j,k)%x)*&
                           (0.5 - vof_node(i,j,k))/(vof_node(i,j,k+1) - vof_node(i,j,k))
                           inter_z(i,j,k)%y = nodes(i,j,k)%y + (nodes(i,j,k+1)%y - nodes(i,j,k)%y)*&
                           (0.5 - vof_node(i,j,k))/(vof_node(i,j,k+1) - vof_node(i,j,k))
                           inter_z(i,j,k)%z = nodes(i,j,k)%z + (nodes(i,j,k+1)%z - nodes(i,j,k)%z)*&
                                (0.5 - vof_node(i,j,k))/(vof_node(i,j,k+1) - vof_node(i,j,k))
                     end if
                  end do
               end do
            end do

            !< Subroutine calls to find "wetted" area of each face in cell
            !< In each case the intersect points change order based on m and n
            call compute_wetted_face_area(Ifacewet, Ifaces, inter_z, inter_y, nodes, dims, 'x')
            call compute_wetted_face_area(Jfacewet, Jfaces, inter_x, inter_z, nodes, dims, 'y')
            call compute_wetted_face_area(Kfacewet, Kfaces, inter_x, inter_y, nodes, dims, 'z')

         end subroutine interface_reconstruction


         subroutine compute_wetted_face_area(A, face, inter_m, inter_n, node, dims, dir)
            !< Computes the area of the wetted surface
            !< Generalised for all directions: I, J, and K.
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            type(nodetype), dimension(-2:dims%imx+3,-2:dims%jmx+3,-2:dims%kmx+3), intent (in) :: nodes
            !< Stores the location of nodes
            type(facetype), dimension(:,:,:), allocatable, intent(inout) :: A
            !< Storing the Area of "wetted" face
            type(facetype), dimension(:,:,:), allocatable, intent(in) :: face
            !< Storing the Area and normal of the directional face
            type(interfacetype), dimension(:,:,:), allocatable, intent(in) :: inter_m
            !< Stores intercept location in the M direction
            type(interfacetype), dimension(:,:,:), allocatable, intent(in) :: inter_n
            !< Stores intercept location in the N direction
            character(len=*), intent(in) :: dir
            integer :: i,j,k
            real(wp) :: l, b, h
            !< To store height and lengths for temporary calculations

            select case(dir)
            case('x')
               !< Visualise Y from bottom to top and Z from left to right
               do k = 1,dims%kmx+1
                  do j = 1,dims%jmx+1
                     do i = 1,dims%imx+1
                        if((inter_n(i,j,k)%z /= 0.0) .and. &
                           (inter_m(i,j+1,k)%z /= 0.0)) then
                           !< slope positive - case 1
                           b = sqrt((inter_m(i,j+1,k)%z - nodes(i,j+1,k)%z)**2&
                                 + (inter_m(i,j+1,k)%y - nodes(i,j+1,k)%y)**2)
                           h = sqrt((inter_n(i,j,k)%z - nodes(i,j+1,k)%z)**2&
                                 + (inter_n(i,j,k)%y - nodes(i,j+1,k)%y)**2)
                           A(i,j,k) = 1/2*b*h
                           if (vof_node(i,j,k+1)>=0.5) then
                              A(i,j,k) = face(i,j,k)%A - A(i,j,k)
                           end if

                        else if((inter_n(i,j,k+1)%z /= 0.0) .and. &
                           (inter_m(i,j+1,k)%z /= 0.0)) then
                           !< slope negative - case 1a
                           b = sqrt((inter_m(i,j+1,k)%z - nodes(i,j+1,k+1)%z)**2&
                                 + (inter_m(i,j+1,k)%y - nodes(i,j+1,k+1)%y)**2)
                           h = sqrt((inter_n(i,j,k+1)%z - nodes(i,j+1,k+1)%z)**2&
                                 + (inter_n(i,j,k+1)%y - nodes(i,j+1,k+1)%y)**2)
                           A(i,j,k) = 1/2*b*h
                           if (vof_node(i,j,k)>=0.5) then
                              A(i,j,k) = face(i,j,k)%A - A(i,j,k)
                           end if

                        else if((inter_n(i,j,k+1)%z /= 0.0) .and. &
                           (inter_m(i,j,k)%z /= 0.0)) then
                           !< Slope positive - case 2
                           b = sqrt((inter_m(i,j,k)%z - nodes(i,j,k+1)%z)**2&
                                 + (inter_m(i,j,k)%y - nodes(i,j,k+1)%y)**2)
                           h = sqrt((inter_n(i,j,k+1)%z - nodes(i,j,k+1)%z)**2&
                                 + (inter_n(i,j,k+1)%y - nodes(i,j,k+1)%y)**2)
                           A(i,j,k) = 1/2*b*h
                           if (vof_node(i,j+1,k)>=0.5) then
                              A(i,j,k) = face(i,j,k)%A - A(i,j,k)
                           end if

                        else if((inter_n(i,j,k)%z /= 0.0) .and. &
                           (inter_m(i,j,k)%z /= 0.0)) then
                           !< Slope negative - case 2a
                           b = sqrt((inter_m(i,j,k)%z - nodes(i,j,k)%z)**2&
                                 + (inter_m(i,j,k)%y - nodes(i,j,k)%y)**2)
                           h = sqrt((inter_n(i,j,k)%z - nodes(i,j,k)%z)**2&
                                 + (inter_n(i,j,k)%y - nodes(i,j,k)%y)**2)
                           A(i,j,k) = 1/2*b*h
                           if (vof_node(i,j+1,k+1)>=0.5) then
                              A(i,j,k) = face(i,j,k)%A - A(i,j,k)
                           end if

                        else if((inter_n(i,j,k)%z /= 0.0) .and. &
                           (inter_n(i,j,k+1)%z /= 0.0)) then
                           !< Slope doesn't matter - case 3
                           h = sqrt((nodes(i,j,k+1)%z - nodes(i,j,k)%z)**2&
                                 + (nodes(i,j,k+1)%y - nodes(i,j,k)%y)**2)
                           a = sqrt((inter_n(i,j,k)%z - nodes(i,j,k)%z)**2&
                                 + (inter_n(i,j,k)%y - nodes(i,j,k)%y)**2)
                           b = sqrt((inter_n(i,j,k+1)%z - nodes(i,j,k)%z)**2&
                                 + (inter_n(i,j,k+1)%y - nodes(i,j,k+1)%y)**2)
                           A(i,j,k) = (a+b)/2*h 
                           if (vof_node(i,j+1,k)>=0.5) then
                              A(i,j,k) = face(i,j,k)%A - A(i,j,k)
                           end if

                        else if ((inter_m(i,j+1,k)%z /= 0.0) .and. &
                                 (inter_m(i,j,k)%z /= 0.0)) then
                           !< Slope somewhat matters - case 4
                           h = sqrt((nodes(i,j+1,k+1)%z - nodes(i,j,k+1)%z)**2&
                                 + (nodes(i,j+1,k+1)%y - nodes(i,j,k+1)%y)**2)
                           a = sqrt((inter_m(i,j,k)%z - nodes(i,j,k)%z)**2&
                                 + (inter_m(i,j,k)%y - nodes(i,j,k)%y)**2)
                           b = sqrt((inter_m(i,j+1,k)%z - nodes(i,j+1,k)%z)**2&
                                 + (inter_m(i,j+1,k)%y - nodes(i,j+1,k)%y)**2)
                           A(i,j,k) = (a+b)/2*h
                           if (vof_node(i,j,k)>=0.5) then
                              A(i,j,k) = face(i,j,k)%A - A(i,j,k)
                           end if
                           !! This case alone needs to be verified for the two possibilities
                        end if
                     end do
                  end do
               end do

            case('y')
               !< Visualise X from left to right and Z from bottom to top
               do k = 1,dims%kmx+1
                  do j = 1,dims%jmx+1
                     do i = 1,dims%imx+1
                        if((inter_n(i,j,k)%x /= 0.0) .and. &
                           (inter_m(i,j,k+1)%x /= 0.0)) then
                           !< slope positive - case 1
                           b = sqrt((inter_m(i,j,k+1)%x - nodes(i,j,k+1)%x)**2&
                                 + (inter_m(i,j,k+1)%z - nodes(i,j,k+1)%z)**2)
                           h = sqrt((inter_n(i,j,k)%x - nodes(i,j,k+1)%x)**2&
                                 + (inter_n(i,j,k)%z - nodes(i,j,k+1)%z)**2)
                           A(i,j,k) = 1/2*b*h
                           if (vof_node(i+1,j,k)>=0.5) then
                              A(i,j,k) = face(i,j,k)%A - A(i,j,k)
                           end if

                        else if((inter_n(i+1,j,k)%x /= 0.0) .and. &
                           (inter_m(i,j,k+1)%x /= 0.0)) then
                           !< slope negative - case 1a
                           b = sqrt((inter_m(i,j,k+1)%x - nodes(i+1,j,k+1)%x)**2&
                                 + (inter_m(i,j,k+1)%z - nodes(i+1,j,k+1)%z)**2)
                           h = sqrt((inter_n(i+1,j,k)%x - nodes(i+1,j,k+1)%x)**2&
                                 + (inter_n(i+1,j,k)%z - nodes(i+1,j,k+1)%z)**2)
                           A(i,j,k) = 1/2*b*h
                           if (vof_node(i,j,k)>=0.5) then
                              A(i,j,k) = face(i,j,k)%A - A(i,j,k)
                           end if

                        else if((inter_n(i+1,j,k)%x /= 0.0) .and. &
                           (inter_m(i,j,k)%x /= 0.0)) then
                           !< Slope positive - case 2
                           b = sqrt((inter_m(i,j,k)%x - nodes(i+1,j,k)%x)**2&
                                 + (inter_m(i,j,k)%z - nodes(i+1,j,k)%z)**2)
                           h = sqrt((inter_n(i+1,j,k)%x - nodes(i+1,j,k)%x)**2&
                                 + (inter_n(i+1,j,k)%z - nodes(i+1,j,k)%z)**2)
                           A(i,j,k) = 1/2*b*h
                           if (vof_node(i,j,k+1)>=0.5) then
                              A(i,j,k) = face(i,j,k)%A - A(i,j,k)
                           end if

                        else if((inter_n(i,j,k)%x /= 0.0) .and. &
                           (inter_m(i,j,k)%x /= 0.0)) then
                           !< Slope negative - case 2a
                           b = sqrt((inter_m(i,j,k)%x - nodes(i,j,k)%x)**2&
                                 + (inter_m(i,j,k)%z - nodes(i,j,k)%z)**2)
                           h = sqrt((inter_n(i,j,k)%x - nodes(i,j,k)%x)**2&
                                 + (inter_n(i,j,k)%z - nodes(i,j,k)%z)**2)
                           A(i,j,k) = 1/2*b*h
                           if (vof_node(i+1,j,k+1)>=0.5) then
                              A(i,j,k) = face(i,j,k)%A - A(i,j,k)
                           end if

                        else if((inter_n(i,j,k)%x /= 0.0) .and. &
                           (inter_n(i+1,j,k)%x /= 0.0)) then
                           !< Slope doesn't matter - case 3
                           h = sqrt((nodes(i+1,j,k)%x - nodes(i,j,k)%x)**2&
                                 + (nodes(i+1,j,k)%z - nodes(i,j,k)%z)**2)
                           a = sqrt((inter_n(i,j,k)%x - nodes(i,j,k)%x)**2&
                                 + (inter_n(i,j,k)%z - nodes(i,j,k)%z)**2)
                           b = sqrt((inter_n(i+1,j,k)%x - nodes(i+1,j,k)%x)**2&
                                 + (inter_n(i+1,j,k)%z - nodes(i+1,j,k)%z)**2)
                           A(i,j,k) = (a+b)/2*h 
                           if (vof_node(i,j,k+1)>=0.5) then
                              A(i,j,k) = face(i,j,k)%A - A(i,j,k)
                           end if

                        else if ((inter_m(i,j,k+1)%y /= 0.0) .and. &
                                 (inter_m(i,j,k)%y /= 0.0)) then
                           !< Slope somewhat matters - case 4
                           h = sqrt((nodes(i+1,j,k+1)%x - nodes(i+1,j,k)%x)**2&
                                 + (nodes(i+1,j,k+1)%z - nodes(i+1,j,k)%z)**2)
                           a = sqrt((inter_m(i,j,k)%x - nodes(i,j,k)%x)**2&
                                 + (inter_m(i,j,k)%z - nodes(i,j,k)%z)**2)
                           b = sqrt((inter_m(i,j,k+1)%x - nodes(i,j,k+1)%x)**2&
                                 + (inter_m(i,j,k+1)%z - nodes(i,j,k+1)%z)**2)
                           A(i,j,k) = (a+b)/2*h
                           if (vof_node(i,j,k)>=0.5) then
                              A(i,j,k) = face(i,j,k)%A - A(i,j,k)
                           end if
                           !! This case alone needs to be verified for the two possibilities
                        end if
                     end do
                  end do
               end do

            case('z')
               !< Visualise X from left to right and Y from bottom to top
               do k = 1,dims%kmx+1
                  do j = 1,dims%jmx+1
                     do i = 1,dims%imx+1
                        if((inter_n(i,j,k)%x /= 0.0) .and. &
                           (inter_m(i,j+1,k)%x /= 0.0)) then
                           !< slope positive - case 1
                           b = sqrt((inter_m(i,j+1,k)%x - nodes(i,j+1,k)%x)**2&
                                 + (inter_m(i,j+1,k)%y - nodes(i,j+1,k)%y)**2)
                           h = sqrt((inter_n(i,j,k)%x - nodes(i,j+1,k)%x)**2&
                                 + (inter_n(i,j,k)%y - nodes(i,j+1,k)%y)**2)
                           A(i,j,k) = 1/2*b*h
                           if (vof_node(i+1,j,k)>=0.5) then
                              A(i,j,k) = face(i,j,k)%A - A(i,j,k)
                           end if

                        else if((inter_n(i+1,j,k)%x /= 0.0) .and. &
                           (inter_m(i,j+1,k)%x /= 0.0)) then
                           !< slope negative - case 1a
                           b = sqrt((inter_m(i,j+1,k)%x - nodes(i+1,j+1,k)%x)**2&
                                 + (inter_m(i,j+1,k)%y - nodes(i+1,j+1,k)%y)**2)
                           h = sqrt((inter_n(i+1,j,k)%x - nodes(i+1,j+1,k)%x)**2&
                                 + (inter_n(i+1,j,k)%y - nodes(i+1,j+1,k)%y)**2)
                           A(i,j,k) = 1/2*b*h
                           if (vof_node(i,j,k)>=0.5) then
                              A(i,j,k) = face(i,j,k)%A - A(i,j,k)
                           end if

                        else if((inter_n(i+1,j,k)%x /= 0.0) .and. &
                           (inter_m(i,j,k)%x /= 0.0)) then
                           !< Slope positive - case 2
                           b = sqrt((inter_m(i,j,k)%x - nodes(i+1,j,k)%x)**2&
                                 + (inter_m(i,j,k)%y - nodes(i+1,j,k)%y)**2)
                           h = sqrt((inter_n(i+1,j,k)%x - nodes(i+1,j,k)%x)**2&
                                 + (inter_n(i+1,j,k)%y - nodes(i+1,j,k)%y)**2)
                           A(i,j,k) = 1/2*b*h
                           if (vof_node(i,j+1,k)>=0.5) then
                              A(i,j,k) = face(i,j,k)%A - A(i,j,k)
                           end if

                        else if((inter_n(i,j,k)%x /= 0.0) .and. &
                           (inter_m(i,j,k)%x /= 0.0)) then
                           !< Slope negative - case 2a
                           b = sqrt((inter_m(i,j,k)%x - nodes(i,j,k)%x)**2&
                                 + (inter_m(i,j,k)%y - nodes(i,j,k)%y)**2)
                           h = sqrt((inter_n(i,j,k)%x - nodes(i,j,k)%x)**2&
                                 + (inter_n(i,j,k)%y - nodes(i,j,k)%y)**2)
                           A(i,j,k) = 1/2*b*h
                           if (vof_node(i+1,j+1,k)>=0.5) then
                              A(i,j,k) = face(i,j,k)%A - A(i,j,k)
                           end if

                        else if((inter_n(i,j,k)%x /= 0.0) .and. &
                           (inter_n(i+1,j,k)%x /= 0.0)) then
                           !< Slope doesn't matter - case 3
                           h = sqrt((nodes(i+1,j,k)%x - nodes(i,j,k)%x)**2&
                                 + (nodes(i+1,j,k)%y - nodes(i,j,k)%y)**2)
                           a = sqrt((inter_n(i,j,k)%x - nodes(i,j,k)%x)**2&
                                 + (inter_n(i,j,k)%y - nodes(i,j,k)%y)**2)
                           b = sqrt((inter_n(i+1,j,k)%x - nodes(i+1,j,k)%x)**2&
                                 + (inter_n(i+1,j,k)%y - nodes(i+1,j,k)%y)**2)
                           A(i,j,k) = (a+b)/2*h 
                           if (vof_node(i,j+1,k)>=0.5) then
                              A(i,j,k) = face(i,j,k)%A - A(i,j,k)
                           end if

                        else if ((inter_m(i,j+1,k)%x /= 0.0) .and. &
                                 (inter_m(i,j,k)%x /= 0.0)) then
                           !< Slope somewhat matters - case 4
                           h = sqrt((nodes(i+1,j+1,k)%x - nodes(i+1,j,k)%x)**2&
                                 + (nodes(i+1,j+1,k)%y - nodes(i+1,j,k)%y)**2)
                           a = sqrt((inter_m(i,j,k)%x - nodes(i,j,k)%x)**2&
                                 + (inter_m(i,j,k)%y - nodes(i,j,k)%y)**2)
                           b = sqrt((inter_m(i,j+1,k)%x - nodes(i,j+1,k)%x)**2&
                                 + (inter_m(i,j+1,k)%y - nodes(i,j+1,k)%y)**2)
                           A(i,j,k) = (a+b)/2*h
                           if (vof_node(i,j,k)>=0.5) then
                              A(i,j,k) = face(i,j,k)%A - A(i,j,k)
                           end if
                           !! This case alone needs to be verified for the two possibilities
                        end if
                     end do
                  end do
               end do
            case DEFAULT
               print*, 'Error: wetted surface calculation'
               Fatal error
            end select

         end subroutine compute_wetted_face_area



         subroutine level_set_coupling(phi_init, vol_frac, dims)
            !< initiating the level set with the volume fraction
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(out) :: phi_init
            !< Output initial value of Level set after coupling
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: vol_frac
            !< Storing the volume fraction value in the domain

            phi_init(:,:,:) = 2*del_h(:,:,:)*(vol_frac(:,:,:)-0.5)

         end subroutine level_set_coupling


         subroutine level_set_advancement(phi, phi_init, grad_phi_x, grad_phi_y, grad_phi_z, &
                                          del_tau, cells, Ifaces, Jfaces, Kfaces, dims)
            !< acquiring the converged value of level-set in
            !< ficticious time
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(out) :: phi
            !< Outputs value of Level set after coupling
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: phi_init
            !< Storing initial value of Level set after coupling
            real(wp), intent(in) :: del_tau
            !< Ficticious time step
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2) :: sign_phi
            !< Storing the value of the sign function
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            !< Input cell quantities: cell volume
            real(wp), dimension(0:dims%imx,0:dims%jmx,0:dims%kmx), intent(out) :: grad_phi_x
            !< Stores value of level-set gradient
            real(wp), dimension(0:dims%imx,0:dims%jmx,0:dims%kmx), intent(out) :: grad_phi_y
            !< Stores value of level-set gradient
            real(wp), dimension(0:dims%imx,0:dims%jmx,0:dims%kmx), intent(out) :: grad_phi_z
            !< Stores value of level-set gradient
            type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
            !< Input varaible which stores I faces' area and unit normal
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
            !< Input varaible which stores J faces' area and unit normal
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
            !< Input varaible which stores K faces' area and unit normal
            real(wp), dimension(:,:,:), allocatable :: mag = 0
            !< Temporary variable for magnitude of gradient
            
            !< Initialiser
            phi(:,:,:) = phi_init(:,:,:)
            call sign_function(sign_phi, phi_init, dims)
            do while(mag /= 1)
               !!< grad_phi is a vector. Need to make use of qp format as shown
               !!< to extract grad_phi in vector form for other calculations
               call compute_gradient_phi(grad_phi_x, phi, phi_init, cells, &
                                       Ifaces, Jfaces, Kfaces, dims, 'x')
               call compute_gradient_phi(grad_phi_y, phi, phi_init, cells, &
                                       Ifaces, Jfaces, Kfaces, dims, 'y')
               call compute_gradient_phi(grad_phi_z, phi, phi_init, cells, &
                                       Ifaces, Jfaces, Kfaces, dims, 'z')
               mag = 1/sqrt(grad_phi_x**2 + grad_phi_y**2 + grad_phi_z**2)
               phi = phi + del_tau*(sign_phi - sign_phi*mag)
            end do
            !!!!< NEED TO INCLUDE BOUNDARY CONDITION VALUES for Phi
         end subroutine level_set_advancement
         

         subroutine sign_function(sign_phi, phi_init, dims)
            !< placing the sign based on level-set  - smoothened 
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: phi_init
            !< Storing initial value of Level set after coupling
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(out) :: sign_phi
            !< Storing the value of the sign function from
            !< the smoothening function
            sign_phi(:,:,:) = phi_init(:,:,:)/(sqrt(phi_init(:,:,:)**2 + del_h(:,:,:)**2))

         end subroutine sign_function

         subroutine compute_gradient_phi(grad, phi, phi_init, cells, Ifaces, Jfaces, Kfaces, dims, dir)
            !< Computes the gradient of the Level-Set
            !< function based on the face value logical assignment
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            character(len=*)                           , intent(in) :: dir
            !< Direction with respect to which gradients are calculated
            real(wp), dimension( 0:dims%imx  , 0:dims%jmx  , 0:dims%kmx  ), intent(out) :: grad
            !< Output variable storing the graident of phi
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: phi
            !< Input variable of phi
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: phi_init
            !< Input variable of phi initial
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            !< Input cell quantities: volume
            type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
            !< Input varaible which stores I faces' area and unit normal
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
            !< Input varaible which stores J faces' area and unit normal
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
            real(wp), dimension(0:1) :: phi_I
            !< Variable to store the face value of phi at Ifaces
            real(wp), dimension(0:1) :: phi_J
            !< Variable to store the face value of phi at Jfaces
            real(wp), dimension(0:1) :: phi_K
            !< Variable to store the face value of phi at Kfaces
            integer :: i, j, k, l

            DebugCall("compute_gradient_phi")
            grad(:,:,:) = 0.0
            select case (dir)
            case ('x')
               do k=0,dims%kmx
                  do j=0,dims%jmx
                    do i=0,dims%imx
                     !< Calculating face values of phi for current cell
                     call face_value_phi(phi_I(0), phi(i,j,k), phi(i-1,j,k), phi_init(i,j,k))
                     call face_value_phi(phi_I(1), phi(i,j,k), phi(i+1,j,k), phi_init(i,j,k))
                     call face_value_phi(phi_J(0), phi(i,j,k), phi(i,j-1,k), phi_init(i,j,k))
                     call face_value_phi(phi_J(1), phi(i,j,k), phi(i,j+1,k), phi_init(i,j,k))
                     call face_value_phi(phi_K(0), phi(i,j,k), phi(i,j,k-1), phi_init(i,j,k))
                     call face_value_phi(phi_K(1), phi(i,j,k), phi(i,j,k+1), phi_init(i,j,k))

                     !< Calculating the gradient using the face values in x-dir
                     grad(i,j,k) = (-phi_I(0)*Ifaces(i,j,k)%nx*Ifaces(i,j,k)%A &
                                   - phi_J(0)*Jfaces(i,j,k)%nx*Jfaces(i,j,k)%A &
                                   - phi*K(0)*Kfaces(i,j,k)%nx*Kfaces(i,j,k)%A &
                                   + phi_I(1)*Ifaces(i+1,j,k)%nx*Ifaces(i+1,j,k)%A &
                                   + phi_J(1)*Jfaces(i,j+1,k)%nx*Jfaces(i,j+1,k)%A &
                                   + phi_K(1)*Kfaces(i,j,k+1)%nx*Kfaces(i,j,k+1)%A &
                                   )/(cells(i,j,k)%volume)
                    end do
                  end do
               end do
            case('y')
               do k=0,dims%kmx
                  do j=0,dims%jmx
                    do i=0,dims%imx
                     !< Calculating face values of phi for current cell
                     call face_value_phi(phi_I(0), phi(i,j,k), phi(i-1,j,k), phi_init(i,j,k))
                     call face_value_phi(phi_I(1), phi(i,j,k), phi(i+1,j,k), phi_init(i,j,k))
                     call face_value_phi(phi_J(0), phi(i,j,k), phi(i,j-1,k), phi_init(i,j,k))
                     call face_value_phi(phi_J(1), phi(i,j,k), phi(i,j+1,k), phi_init(i,j,k))
                     call face_value_phi(phi_K(0), phi(i,j,k), phi(i,j,k-1), phi_init(i,j,k))
                     call face_value_phi(phi_K(1), phi(i,j,k), phi(i,j,k+1), phi_init(i,j,k))

                     !< Calculating the gradient using the face values in y-dir
                     grad(i,j,k) = (-phi_I(0)*Ifaces(i,j,k)%ny*Ifaces(i,j,k)%A &
                                   - phi_J(0)*Jfaces(i,j,k)%ny*Jfaces(i,j,k)%A &
                                   - phi*K(0)*Kfaces(i,j,k)%ny*Kfaces(i,j,k)%A &
                                   + phi_I(1)*Ifaces(i+1,j,k)%ny*Ifaces(i+1,j,k)%A &
                                   + phi_J(1)*Jfaces(i,j+1,k)%ny*Jfaces(i,j+1,k)%A &
                                   + phi_K(1)*Kfaces(i,j,k+1)%ny*Kfaces(i,j,k+1)%A &
                                   )/(cells(i,j,k)%volume)
                    end do
                  end do
               end do
            case ('z')
               do k=0,dims%kmx
                  do j=0,dims%jmx
                    do i=0,dims%imx
                     !< Calculating face values of phi for current cell
                     call face_value_phi(phi_I(0), phi(i,j,k), phi(i-1,j,k), phi_init(i,j,k))
                     call face_value_phi(phi_I(1), phi(i,j,k), phi(i+1,j,k), phi_init(i,j,k))
                     call face_value_phi(phi_J(0), phi(i,j,k), phi(i,j-1,k), phi_init(i,j,k))
                     call face_value_phi(phi_J(1), phi(i,j,k), phi(i,j+1,k), phi_init(i,j,k))
                     call face_value_phi(phi_K(0), phi(i,j,k), phi(i,j,k-1), phi_init(i,j,k))
                     call face_value_phi(phi_K(1), phi(i,j,k), phi(i,j,k+1), phi_init(i,j,k))

                     !< Calculating the gradient using the face values in z-dir
                     grad(i,j,k) = (-phi_I(0)*Ifaces(i,j,k)%nz*Ifaces(i,j,k)%A &
                                   - phi_J(0)*Jfaces(i,j,k)%nz*Jfaces(i,j,k)%A &
                                   - phi*K(0)*Kfaces(i,j,k)%nz*Kfaces(i,j,k)%A &
                                   + phi_I(1)*Ifaces(i+1,j,k)%nz*Ifaces(i+1,j,k)%A &
                                   + phi_J(1)*Jfaces(i,j+1,k)%nz*Jfaces(i,j+1,k)%A &
                                   + phi_K(1)*Kfaces(i,j,k+1)%nz*Kfaces(i,j,k+1)%A &
                                   )/(cells(i,j,k)%volume)
                    end do
                  end do
               end do
            case DEFAULT
               print *, "ERROR: gradient direction error for grad phi"
            end select
                     end subroutine compute_gradient_phi


         subroutine face_value_phi(phi_A, phi_P, phi_Ci, phi_init)
            !< Calcualting face values of phi for particular cell
            implicit none
            real(wp), intent(in) :: phi_P
            !< primary cell LS value
            real(wp), intent(in) :: phi_Ci
            !< neighbouring cell LS value
            real(wp), intent(in) :: phi_init
            !< current initial LS value to indicate the phase
            real(wp), intent(out) :: phi_A
            !< Returns face value of phi at face A

            if (phi_init > 0) then   !< For fluid phase 1
               if (phi_Ci - phi_P > 0) then
                  phi_A = phi_P
               else
                  phi_A = phi_Ci
               end if
            else if (phi_init < 0) then  !< For fluid phase 2
               if (phi_Ci - phi_P < 0) then
                  phi_A = phi_P
               else
                  phi_A = phi_Ci
               end if
            end if
         end subroutine face_value_phi


         ! subroutine level_set_face ()
         !    !< approximating face value of the level set based on
         !    !< S(\phi^0).\vec n
         !    implicit none
         ! end subroutine level_set_face


         subroutine surface_tension_force(sigma, K, d_delta, grad_phi_x, grad_phi_y, grad_phi_z, dims)
            !< obtaining surface tension force from dirac delta,
            !< curvature, and new level set function
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            real(wp), dimension(:,:,:), allocatable, intent(in) :: K
            !< Curvature of the level set - Kappa
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: d_delta
            !< Input cell quantities: cell center
            real(wp), dimension(0:dims%imx,0:dims%jmx,0:dims%kmx), intent(in) :: grad_phi_x
            !< Stores the value of gradient of level set
            real(wp), dimension(0:dims%imx,0:dims%jmx,0:dims%kmx), intent(in) :: grad_phi_y
            !< Stores the value of gradient of level set
            real(wp), dimension(0:dims%imx,0:dims%jmx,0:dims%kmx), intent(in) :: grad_phi_z
            !< Stores the value of gradient of level set
            real(wp), intent(in) :: sigma
            !< Surface tension at interface - fluid property
            integer :: i, j, k

            ! To calculate surface tension force
            do k = 0,dims%kmx
               do j = 0,dims%jmx
                  do i = 0,dims%imx
                     F_x(i,j,k) = sigma*K(i,j,k)*d_delta(i,j,k)*grad_phi_x(i,j,k)
                     F_y(i,j,k) = sigma*K(i,j,k)*d_delta(i,j,k)*grad_phi_y(i,j,k)
                     F_z(i,j,k) = sigma*K(i,j,k)*d_delta(i,j,k)*grad_phi_z(i,j,k)
                  end do
               end do
            end do

         end subroutine surface_tension_force


         subroutine curvature(K, grad_phi_x, grad_phi_y, grad_phi_z, cells, Ifaces, Jfaces, Kfaces, dims)
            !< getting curvature from level set
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            real(wp), dimension(:,:,:), allocatable, intent(out) :: K
            !< Output variable storing the gradient of curvature
            real(wp), dimension(0:dims%imx,0:dims%jmx,0:dims%kmx), intent(in) :: grad_phi_x
            !< Stores the value of gradient of level set
            real(wp), dimension(0:dims%imx,0:dims%jmx,0:dims%kmx), intent(in) :: grad_phi_y
            !< Stores the value of gradient of level set
            real(wp), dimension(0:dims%imx,0:dims%jmx,0:dims%kmx), intent(in) :: grad_phi_z
            !< Stores the value of gradient of level set
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            !< Input cell quantities: cell volume
            type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
            !< Input varaible which stores I faces' area and unit normal
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
            !< Input varaible which stores J faces' area and unit normal
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
            !< Input varaible which stores K faces' area and unit normal

            real(wp), dimension(:,:,:), allocatable :: K_x
            !< Temporary variable storing the gradient of curvature
            real(wp), dimension(:,:,:), allocatable :: K_y
            !< Temporary variable storing the gradient of curvature
            real(wp), dimension(:,:,:), allocatable :: K_z
            !< Temporary variable storing the gradient of curvature
            real(wp), dimension(:,:,:), allocatable :: mag
            !< To temporarily store the value of the magnitude of grad_phi
            
            mag = 1/sqrt(grad_phi_x**2 + grad_phi_y**2 + grad_phi_z**2)
            call compute_gradient_G(K_x, grad_phi_x/mag, cells, Ifaces, Jfaces, Kfaces, dims, 'x')
            call compute_gradient_G(K_y, grad_phi_y/mag, cells, Ifaces, Jfaces, Kfaces, dims, 'y')
            call compute_gradient_G(K_z, grad_phi_z/mag, cells, Ifaces, Jfaces, Kfaces, dims, 'z')
            K = K_x + K_y + K_z

         end subroutine curvature


         subroutine dirac_delta(d_delta, phi, epsilon, cells, dims)
            !< initialising smooth dirac delta with level set function
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(out) :: d_delta
            !< Output variable storing the Dirac Delta smooth function
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: phi
            !< New Level-Set function
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            !< Input cell quantities: cell centers
            real(wp), intent(in) :: epsilon
            !< Numerical interface width
            ! To calcualte Dirac Delta function
            d_delta(:,:,:) = 0
            if (abs(phi(:,:,:)) <= epsilon) then
               ! this is the tiny portion within the interface
               d_delta(:,:,:) = 1.0/(2.0*epsilon)*(1 + cos(pi*phi(:,:,:)/epsilon))
            end if
           
         end subroutine dirac_delta


         subroutine heaviside(H, phi, epsilon, cells, dims)
            !< forming heaviside function based on level-set
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(out) :: H
            !< Output variable storing the Heaviside values
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: phi
            !< New Level-Set function
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            !< Input cell quantities: cell centers
            real(wp), intent(in) :: epsilon
            !< Numerical interface width

            ! To calcualte heaviside function
            if (phi(:,:,:) < -1*epsilon) then 
               ! when LS is below interface limit
               H(:,:,:) = 0
            else if (phi(:,:,:) > epsilon) then
               ! when LS is above interface limit
               H(:,:,:) = 1
            else
               ! when LS is in the interface band
               H(:,:,:) = 0.5*(1.0 + phi(:,:,:)/epsilon + sin(pi*phi(:,:,:)/epslion))
            end if
                      
         end subroutine heaviside


         subroutine smoothen_G(G, G1, G2, H, cells, dims)
            !< to smoothen two scalars over interface
            !< between two fluids using Heaviside function
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            real(wp), intent(in) :: G1
            !< Input variable for scalar 1 over one half of the interface
            real(wp), intent(in) :: G2
            !< Input variable for scalar 2 over second half of the interface
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: H
            !< New Heaviside function
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(out) :: G
            !< Smoothened function G that uses G1 and G2 to form a field variable G
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            !< Input cell quantities: cell centers
            ! To Smoothen function
            G(:,:,:) = G1*(1-H(:,:,:)) + G2*H(:,:,:)
            !!!! NEED TO CHANGE G1 and G2 to matrices

         end subroutine smoothen_G

         subroutine area_of_polygon(X,Y,n,area)
            !< Calculates area of a polygon using number of points
            !< and vertices coordinates
            implicit none
            integer, intent(in) :: n
            !< Stores the number of points in the polygon
            real(wp), dimension(:), allocatable, intent(in) :: X
            real(wp), dimension(:), allocatable, intent(in) :: Y
            !< Stores X and Y coordinate data
            real(wp), intent(inout) :: area = 0
            integer :: i,j

            j = n - 1
            do i = 1,n
               area = area + (X(i) + X(j))*(Y(i)+Y(j))
               j = i
            end do
            area = abs(area/2)
         end subroutine area_of_polygon

                           ! ! ! ! subroutine setup_clsvof(control, scheme, flow, dims)
         ! ! ! !    !< allocate array memory for data communication
         ! ! ! !    implicit none
         ! ! ! !    type(controltype), intent(in) :: control
         ! ! ! !    !< Control parameters
         ! ! ! !    type(schemetype), intent(in) :: scheme
         ! ! ! !    !< finite-volume Schemes
         ! ! ! !    type(flowtype), intent(in) :: flow
         ! ! ! !    !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
         ! ! ! !    type(extent), intent(in) :: dims
         ! ! ! !    !< Extent of the domain:imx,jmx,kmx
         ! ! ! !    character(len=*), parameter :: errmsg="module: CLSVOF_incomp, subroutine setup"
         ! ! ! !    !< Error message
         
         ! ! ! !    imx = dims%imx
         ! ! ! !    jmx = dims%jmx
         ! ! ! !    kmx = dims%kmx
         ! ! ! !    n_var = control%n_var
         ! ! ! !    gm = flow%gm
         ! ! ! !    mu_ref = flow%mu_ref
         ! ! ! !    Reynolds_number = flow%Reynolds_number
         ! ! ! !    free_stream_tu = flow%tu_inf
         ! ! ! !    tk_inf = flow%tk_inf
         ! ! ! !    tkl_inf = flow%tkl_inf
         ! ! ! !    tpr = flow%tpr
         ! ! ! !    pr = flow%pr
         ! ! ! !    R_gas = flow%R_gas
         
         ! ! ! !    call alloc(delQ, 0, imx, 0, jmx, 0, kmx, 1, n_var)
         ! ! ! !    call alloc(delQstar, 0, imx, 0, jmx, 0, kmx, 1, n_var)
         
         ! ! ! !    if(mu_ref==0.0 .or. scheme%turbulence=='none') then
         ! ! ! !      call alloc(dummy, 0, imx, 0, jmx, 0, kmx)
         ! ! ! !      dummy = 0.0
         ! ! ! !    end if
         ! ! ! !    if(mu_ref==0.0)then
         ! ! ! !      mmu => dummy
         ! ! ! !    else
         ! ! ! !      mmu => mu
         ! ! ! !    end if
         ! ! ! !    if(trim(scheme%turbulence)=='none')then
         ! ! ! !      tmu => dummy
         ! ! ! !    else
         ! ! ! !      tmu => mu_t
         ! ! ! !    end if
         ! ! ! ! end subroutine setup_clsvof