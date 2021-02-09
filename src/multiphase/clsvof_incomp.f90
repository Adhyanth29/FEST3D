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
   pi = 4.0*atan(1.0)
   !< An accurate definition of pi

   contains

         subroutine vol_frac_adv(vol_frac_n, vol_frac_o, qp, cells, Ifaces, Jfaces, Kfaces, del_t, del_h, dims)
            !< to account for the volume fraction advection VOF
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            real(wp), intent(in) :: del_t
            !< Time step
            real(wp), intent(out) :: del_h
            !< Stores the value of cell size
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(out) :: vol_frac_n
            !< Output the next time-step of volume fraction
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: vol_frac_o
            !< Storing the previous time-step volume fraction
            real(wp), dimension(:, :, :), pointer :: x_speed      
            !< U pointer, point to slice of qp (:,:,:,2) 
            real(wp), dimension(:, :, :), pointer :: y_speed      
            !< V pointer, point to slice of qp (:,:,:,3) 
            real(wp), dimension(:, :, :), pointer :: z_speed      
            !< W pointer, point to slice of qp (:,:,:,4)
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2) :: grad_x
            !< To store the gradient of velocity times volume fraction in x
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2) :: grad_y
            !< To store the gradient of velocity times volume fraction in y
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2) :: grad_z
            !< To store the gradient of velocity times volume fraction in z
            type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
            !< Input varaible which stores I faces' area and unit normal
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
            !< Input varaible which stores J faces' area and unit normal
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces

            !< Pointer allocation
            x_speed(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2) => qp(:, :, :, 2)
            y_speed(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2) => qp(:, :, :, 3)
            z_speed(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2) => qp(:, :, :, 4)
      
            call compute_gradient_G(grad_x, x_speed*vol_frac_o, cells, Ifaces, Jfaces, Kfaces, dims, 'x')
            call compute_gradient_G(grad_y, y_speed*vol_frac_o, cells, Ifaces, Jfaces, Kfaces, dims, 'y')
            call compute_gradient_G(grad_z, z_speed*vol_frac_o, cells, Ifaces, Jfaces, Kfaces, dims, 'z')
            vol_frac_n(:,:,:) = vol_frac_o(:,:,:) - del_t/del_h*(grad_x(:,:,:) + grad_y(:,:,:) & 
            + grad_z(:,:,:))
            !!!
            !!!
            !!!
            ! This will NOT work as of now. Need to account for the bondary grad values to make sure the volume fractions are calculated accordingly for the full domain

            !< Should integrate this with the interface reconstruction to obtain solution
            !< Also account for the ghost cells to ensure proper boundary values
         end subroutine vol_frac_adv

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

         subroutine cell_size(del_h, cells, face, dims)
            !< to find the cell size required for this module
            implicit none 
            real(wp), intent(out) :: del_h
            !< Stores the value of cell size
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            !< Stores cell parameter: volume
            type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: face
            !< Input varaible which stores surface area of face
            del_h = cells%volume / face%A

         end subroutine cell_size      

         subroutine interface_recons()
            !< to reconstruct interface using vof 0.5
            implicit none
         end subroutine interface_recons


         subroutine level_set_coupling(phi_init, vol_frac, del_h, dims)
            !< initiating the level set with the volume fraction
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            real(wp), intent(in) :: del_h
            !< Storing value of cell size
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(out) :: phi_init
            !< Output initial value of Level set after coupling
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: vol_frac
            !< Storing the volume fraction value in the domain

            phi_init(:,:,:) = 2*del_h*(vol_frac(:,:,:)-0.5)

         end subroutine level_set_coupling

         subroutine compute_gradient_phi()
            !< Computes the gradient of the Level-Set
            !< function based on the face value mechanism
            !< in the paper
            implicit none
         end subroutine compute_gradient_phi


         subroutine level_set_advancement(phi, phi_init, grad_phi_x, grad_phi_y, grad_phi_z, sign_phi, del_t, cells, Ifaces, Jfaces, Kfaces, dims)
            !< acquiring the converged value of level-set in
            !< ficticious time
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(out) :: phi
            !< Outputs value of Level set after coupling
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: phi_init
            !< Storing initial value of Level set after coupling
            real(wp), intent(in) :: del_t
            !< Storing the value of time step
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(out) :: sign_phi
            !< Storing the value of the sign function
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            !< Input cell quantities: cell volume
            real(wp), dimension(0:dims%imx,0:dims%jmx,0:dims%kmx), intent(in) :: grad_phi_x
            !< Stores value of level-set gradient
            real(wp), dimension(0:dims%imx,0:dims%jmx,0:dims%kmx), intent(in) :: grad_phi_y
            !< Stores value of level-set gradient
            real(wp), dimension(0:dims%imx,0:dims%jmx,0:dims%kmx), intent(in) :: grad_phi_z
            !< Stores value of level-set gradient
            type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
            !< Input varaible which stores I faces' area and unit normal
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
            !< Input varaible which stores J faces' area and unit normal
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
            !< Input varaible which stores K faces' area and unit normal

            !!< grad_phi is a vector. Need to make use of qp format as shown
            !!< to extract grad_phi in vector form for other calculations


            
         end subroutine level_set_advancement

         subroutine sign_function(sign_phi, phi_init, del_h, dims)
            !< placing the sign based on level-set  - smoothened 
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: phi_init
            !< Storing initial value of Level set after coupling
            real(wp), intent(in) :: del_h
            !< Storing the value of cell size
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(out) :: sign_phi
            !< Storing the value of the sign function from
            !< the smoothening function
            sign_phi(:,:,:) = phi_init(:,:,:)/(sqrt( phi_init(:,:,:)**2 + del_h**2))

         end subroutine sign_function


         ! subroutine level_set_face ()
         !    !< approximating face value of the level set based on
         !    !< S(\phi^0).\vec n
         !    implicit none
         ! end subroutine level_set_face


         subroutine surface_tension_force(F, sigma, K, d_delta, grad_phi_x, grad_phi_y, grad_phi_z, dims)
            !< obtaining surface tension force from dirac delta,
            !< curvature, and new level set function
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(out) :: F
            !< Surface tension force to be calcualted
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
            do k = 0:dims%kmx
               do j = 0:dims%jmx
                  do i = 0:dims%imx
                     F(i,j,k) = sigma*K(i,j,k)*d_delta(i,j,k)*grad_phi(i,j,k)
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
               d_delta(:,:,:) = 1/(2*epsilon)*(1 + cos(pi*phi(:,:,:)/epsilon))
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
               H(:,:,:) = 0.5*(1 + phi(:,:,:)/epsilon + sin(pi*phi(:,:,:)/epslion))
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

         end subroutine smoothen_G