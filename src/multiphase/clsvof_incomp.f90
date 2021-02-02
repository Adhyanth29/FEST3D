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
   real(wp), parameter :: pi = 3.14159265


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
         

         subroutine interface_recons()
            !< to reconstruct interface with the 4 filling cases
            implicit none
         end subroutine interface_recons


         subroutine level_set_coupling()
            !< initiating the level set with the volume fraction
            implicit none
         end subroutine level_set_coupling


         subroutine level_set_advancement()
            !< acquiring the converged value of level-set in
            !< ficticious time
            !< also calculates gradient of level set for convergence
            implicit none
         end subroutine level_set_advancement

         subroutine sign_function()
            !< placing the sign based on level-set  - smoothened 
            implicit none
         end subroutine sign_function


         subroutine level_set_face ()
            !< approximating face value of the level set based on
            !< S(\phi^0).\vec n
            implicit none
         end subroutine      


         subroutine surface_tension_force()
            !< obtaining surface tension force from dirac delta,
            !< curvature, and new level set function
            implicit none
         end subroutine surface_tension_force


         subroutine curvature()
            !< getting curvature from level set
            implicit none
         end subroutine curvature


         subroutine dirac_delta()
            !< initialising smooth dirac delta with level set function
            implicit none
         end subroutine dirac_delta


         subroutine heaviside(phi, epsilon, cells, dims)
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
            do k = 1:dims%kmx-1
               do j = 1:dims%jmx-1
                  do i = 1:dims%imx-1
                     if (phi(i,j,k) < -1*epsilon) then 
                        ! when LS is below interface limit
                        H(i,j,k) = 0
                     else if (phi(i,j,k) > epsilon) then
                        ! when LS is above interface limit
                        H(i,j,k) = 1
                     else
                        ! when LS is in the interface band
                        H(i,j,k) = 0.5*(1 + phi(i,j,k)/epsilon + SIN(pi*phi/epslion))
                     end if
                  end do
               end do
            end do             
         end subroutine heaviside


         subroutine smoothen_G(G1, G2, H, cells, dims)
            !< to smoothen two scalars over interface
            !< between two fluids using Heaviside function
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx/2,-2:dims%kmx+2), intent(out) :: G1
            !< Input variable for scalar 1 over one half of the interface
            real(wp), dimension(-2:dims%imx+2,dims%jmx/2+1:dims%jmx +2,-2:dims%kmx+2), intent(out) :: G2
            !< Input variable for scalar 2 over second half of the interface
            real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: H
            !< New Heaviside function
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            !< Input cell quantities: cell centers
            real(wp), intent(in) :: epsilon
            !< Numerical interface width
            integer :: i
            integer :: j
            integer :: k

            ! To Smoothen function
            do k = 1:dims%kmx-1
               do j = 1:dims%jmx-1
                  do i = 1:dims%imx-1
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
         end subroutine smoothen_G