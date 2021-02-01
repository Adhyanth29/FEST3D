   !< Incompressible Two Phase Model - CLSVOF
module clsvof_incomp
   !<
   !< Reference: Yeng-Yung Tsui, Cheng-Yen Liu & Shi-Wen Lin "Coupled 
   !< level-set and volume-of-fluid method for two-phase flow 
   !< calculations", Numerical Heat Transfer, Part B: Fundamentals,
   !< 71:2, 173-185, 2017, DOI: 10.1080/10407790.2016.1265311
   !<
#include ".../error.h"
#include "../debug.h"

   use vartypes
   use gradients, only : compute_gradient_G
   implicit none
   private


   contains

         subroutine volume_fraction_advection()
            !< to account for the volume fraction advection VOF
            implicit none
         end subroutine volume_fraction_advection


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

         subroutine heaviside()
            !< forming heaviside function based on level-set
            implicit none
         end subroutine heaviside


         subroutine smoothen_G()
            !< to smoothen a scalar G along the interface between
            implicit none 
            !< two fluids
         end subroutine smoothen_G