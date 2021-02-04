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


   contains

         subroutine volume_fraction_advection()
            !< to account for the volume fraction advection VOF
            implicit none
         end subroutine volume_fraction_advection

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

         subroutine cell_size(del_h, cells, face)
            !< to find the cell size required for this module
            implicit none 
            real(wp), intent(out) :: del_h
            !< Stores the value of cell size
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            !< Stores cell parameter: volume
            type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: face
            !< Input varaible which stores surface area of face
            del_h = cells%volume / face%A
         end subroutine cell_size      

         subroutine interface_recons()
            !< to reconstruct interface with the 4 filling cases
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

            integer :: i, j, k

            do k= -2:dims%kmx+2
               do j= -2:dims%jmx+2
                  do i = -2:dims%imx+2
                     phi_init(i,j,k) = 2*del_h*(vol_frac(i,j,k)-0.5)
                  end do
               end do
            end do
         end subroutine level_set_coupling


         subroutine level_set_advancement()
            !< acquiring the converged value of level-set in
            !< ficticious time
            !< also calculates gradient of level set for convergence
            implicit none
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
            integer :: i, j, k

            do k = -2:dims%kmx+2
               do j = -2:dims%jmx+2
                  do i = -2:dims%imx+2
                     sign_phi(i,j,k) = phi_init(i,j,k)/(sqrt( phi_init(i,j,k)**2 + del_h**2))
                  end do
               end do
            end do
         end subroutine sign_function


         ! subroutine level_set_face ()
         !    !< approximating face value of the level set based on
         !    !< S(\phi^0).\vec n
         !    implicit none
         ! end subroutine level_set_face


         subroutine surface_tension_force(F, sigma, K, d_delta, grad_phi, dims)
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
            real(wp), dimension(0:dims%imx,0:dims%jmx,0:dims%kmx), intent(in) :: grad_phi
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


         subroutine curvature(K, grad_phi_x, dims, dir)
            !< getting curvature from level set
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            real(wp), dimension(:,:,:), allocatable, intent(out) :: K
            !< Output variable storing the gradient of curvature
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
            integer :: i
            integer :: j
            integer :: k

            ! To calcualte Dirac Delta function
            d_delta(:,:,:) = 0
            do k = 0:dims%kmx
               do j = 0:dims%jmx
                  do i = 0:dims%imx
                     if (abs(phi(i,j,k)) <= epsilon) then
                        ! this is the tiny portion within the interface
                        d_delta(i,j,k) = 1/(2*epsilon)*(1 + cos(pi*phi(i,j,k)/epsilon))
                     end if
                  end do
               end do
            end do             
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
            integer :: i
            integer :: j
            integer :: k

            ! To calcualte heaviside function
            do k = 0:dims%kmx
               do j = 0:dims%jmx
                  do i = 0:dims%imx
                     if (phi(i,j,k) < -1*epsilon) then 
                        ! when LS is below interface limit
                        H(i,j,k) = 0
                     else if (phi(i,j,k) > epsilon) then
                        ! when LS is above interface limit
                        H(i,j,k) = 1
                     else
                        ! when LS is in the interface band
                        H(i,j,k) = 0.5*(1 + phi(i,j,k)/epsilon + sin(pi*phi/epslion))
                     end if
                  end do
               end do
            end do             
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
            integer :: i
            integer :: j
            integer :: k

            ! To Smoothen function
            do k = 0:dims%kmx
               do j = 0:dims%jmx
                  do i = 0:dims%imx
                     G(i,j,k) = G1*(1-H(i,j,k)) + G2*H(i,j,k)
                  end do
               end do
            end do
         end subroutine smoothen_G