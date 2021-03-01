!< Multiphase calculation
module multi_phase
!< Multiphase calculation

#include "../../../debug.h"
#include "../../../error.h"
   use vartypes
   use clsvof_incomp, only: perform_clsvof_incomp => perform_multiphase
   implicit none
   private

   contains
         subroutine setup_multiphase_scheme()
            !< sets up the mutiphase schemes
            implicit none
         end subroutine setup_multiphase_scheme

         subroutine perform_multiphase()
            !< Authorises the type of multiphase operation
            implicit none
         end subroutine perform_multiphase