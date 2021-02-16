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
   real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2) :: Fx
   !< Surface tension force to be calcualted - X component
   real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2) :: Fy
   !< Surface tension force to be calcualted - Y component
   real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2) :: Fz
   !< Surface tension force to be calcualted - Z component
   type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2) :: del_h
   !< Cell_size of each cell (approximated)

   !< Public members
   public :: cell_size

   contains

         subroutine vof_adv(vof_n, vof_o, qp, cells, Ifaces, Jfaces, Kfaces, del_t, dims)
            !< to account for the volume fraction advection VOF
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            real(wp), intent(in) :: del_t
            !< Time step
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
      
            call compute_gradient_G(grad_x, x_speed*vof_o, cells, Ifaces, Jfaces, Kfaces, dims, 'x')
            call compute_gradient_G(grad_y, y_speed*vof_o, cells, Ifaces, Jfaces, Kfaces, dims, 'y')
            call compute_gradient_G(grad_z, z_speed*vof_o, cells, Ifaces, Jfaces, Kfaces, dims, 'z')
            vof_n(:,:,:) = vof_o(:,:,:) - del_t/del_h*(grad_x(:,:,:) + grad_y(:,:,:) & 
            + grad_z(:,:,:))
            !!!
            !!!
            !!!
            ! This will NOT work as of now. Need to account for the bondary grad values to make sure the volume fractions are calculated accordingly for the full domain

            !< Should integrate this with the interface reconstruction to obtain solution
            !< Also account for the ghost cells to ensure proper boundary values
         end subroutine vof_adv

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

         subroutine cell_size(cells, face, dims)
            !< to find the cell size required for this module
            implicit none 
            real(wp), intent(out) :: del_h
            !< Stores the value of cell size
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            !< Stores cell parameter: volume
            del_h(:,:,:) = cells%volume(:,:,:)**(1.0/3.0)
         end subroutine cell_size      

         subroutine interface_recons(vof, cells, nodes, dims)
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
            !< Stores intercept location in the Z direction
            type(interfacetype), dimension(:,:,:), allocatable :: inter_y
            !< Stores intercept location in the Z direction
            type(interfacetype), dimension(:,:,:), allocatable :: inter_z
            !< Stores intercept location in the Z direction
            type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(out) :: Ifacewet
            !< Input varaible which stores I faces' area and unit normal
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(out) :: Jfacewet
            !< Input varaible which stores J faces' area and unit normal
            type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(out) :: Kfaceswet
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

            if (vof(:,:,:) == 1) then
               Ifacewet(:,:,:) = Ifaces(:,:,:)%A
               Jfacewet(:,:,:) = Jfaces(:,:,:)%A
               Kfacewet(:,:,:) = Kfaces(:,:,:)%A
            end if
            !< To find the vof value at the nodes using adjacent cell centers
            do k = 0:dims%kmx+1
               do j = 0:dims%jmx+1
                  do i = 0:dims%imx+1
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
            
            !< Linear interpolation to find intercept locations where vof = 0.5
            do k = 0:dims%kmx
               do j = 0:dims%jmx
                  do i = 0:dims%imx
                     if((vof_node(i,j,k) < 0.5 .or. vof_node(i+1,j,k) < 0.5) .and. &
                        (vof_node(i,j,k) >= 0.5 .or. vof_node(i+1,j,k) >= 0.5)) then
                           !!! Write expression to signify node point x interpolation
                           inter_x(i,j,k)%x = nodes%x(i,j,k) + (nodes%x(i+1,j,k) - nodes%x(i,j,k))*&
                               (0.5 - vof_node(i,j,k))/(vof_node(i+1,j,k) - vof_node(i,j,k))
                           inter_x(i,j,k)%y = nodes%y(i,j,k)
                           inter_x(i,j,k)%z = nodes%z(i,j,k)
                     end if
                     if((vof_node(i,j,k) < 0.5 .or. vof_node(i,j+1,k) < 0.5) .and. &
                        (vof_node(i,j,k) >= 0.5 .or. vof_node(i,j+1,k) >= 0.5)) then
                           !!! Write expression to signify node point y interpolation
                           inter_y(i,j,k)%x = nodes%x(i,j,k)
                           inter_y(i,j,k)%y = nodes%y(i,j,k) + (nodes%y(i,j+1,k) - nodes%y(i,j,k))*&
                                (0.5 - vof_node(i,j,k))/(vof_node(i,j+1,k) - vof_node(i,j,k))
                           inter_y(i,j,k)%z= nodes%z(i,j,k)
                     end if
                     if((vof_node(i,j,k) < 0.5 .or. vof_node(i,j,k+1) < 0.5) .and. &
                        (vof_node(i,j,k) >= 0.5 .or. vof_node(i,j,k+1) >= 0.5)) then
                           !!! Write expression to signify node point z interpolation
                           inter_z(i,j,k)%x = nodes%x(i,j,k)
                           inter_z(i,j,k)%y = nodes%y(i,j,k)
                           inter_z(i,j,k)%z = nodes%z(i,j,k) + (nodes%z(i,j,k+1) - nodes%z(i,j,k))*&
                                (0.5 - vof_node(i,j,k))/(vof_node(i,j,k+1) - vof_node(i,j,k))
                     end if
                  end do
               end do
            end do

            !< Subroutine calls to find "wetted" area of each face in cell
            call wetted_area(Ifacewet, Ifaces, inter_z, inter_y, nodes, dims, 'x')
            call wetted_area(Jfacewet, Jfaces, inter_x, inter_z, nodes, dims, 'y')
            call wetted_area(Kfacewet, Kfaces, inter_x, inter_y, nodes, dims, 'z')

         end subroutine interface_recons

         subroutine wetted_area(A, face, inter_m, inter_n, node, dims, dir)
            !< Computes the area of the wetted surface
            implicit none
            type(extent), intent(in) :: dims
            !< Extent of domain: imx, jmx, kmx
            type(nodetype), dimension(-2:dims%imx+3,-2:dims%jmx+3,-2:dims%kmx+3), intent (in) :: nodes
            !< Stores the location of nodes
            type(facetype), dimension(:,:,:), allocatable, intent(inout) :: A
            !< Storing the Area of "wetted" face
            type(facetype), dimension(:,:,:), allocatable, intent(in) :: face
            !< Storing the Area and normal of the directional face
            type(interfacetype), dimension(:,:,:), allocatable, intent(in) :: inter_x
            !< Stores intercept location in the M direction
            type(interfacetype), dimension(:,:,:), allocatable, intent(in) :: inter_y
            !< Stores intercept location in the N direction
            character(len=*), intent(in) :: dir
            integer :: i,j,k
            real(wp) :: a, b, h
            !< To store height and lengths

            select case(dir)
            case('x')
               do k = 0:dims%kmx
                  do j = 0:dims%jmx
                     do i = 0:dims%imx
                        if((inter_y(i,j,k)%x == nodes(i,j,k)%x) .and. &
                           inter_x(i+1,j,k)%y == nodes(i+1,j,k)%y) then
                                 b = abs(inter_x(i+1,j,k)%x - nodes(i+1,j,k)%x)
                                 h = abs(inter_y(i,j,k)%y - nodes(i,j+1,k)%y)
                                 A(i,j,k) = 1/2*b*h
                        else if((inter_y(i+1,j,k)%x == nodes(i+1,j,k)%x) .and. &
                                 inter_x(i,j,k)%y == nodes(i,j,k)%y)) then
                                 b = abs(inter_x(i,j,k)%x - nodes(i+1,j,k)%x)
                                 h = abs(inter_y(i+1,j,k)%y - nodes(i+1,j,k)%y)
                                 A(i,j,k) = 1/2*b*h
                        else if((inter_y(i,j,k)%x == nodes(i,j,k)%x) .and. &
                                 inter_y(i+1,j,k)%x == nodes(i+1,j,k)%x)) then
                                 h = abs(nodes(i,j,k)%x - nodes(i+1,j,k)%x)
                                 a = abs(nodes(i,j,k)%y - inter_y(i,j,k)%y)
                                 b = abs(nodes(i+1,j,k)%y - inter_y(i+1,j,k)%y)
                              A(i,j,k) = (a+b)/2*h 
                        else if ((inter_x(i,j+1,k)%y == nodes(i,j+1,k)%y) .and. &
                                 inter_x(i,j,k)%y == nodes(i,j,k)%y)) then
                                 h = abs(nodes(i+1,j+1,k)%y - nodes(i+1,j,k)%x)
                                 a = abs(nodes(i+1,j,k)%x - inter_x(i,j,k)%x)
                                 b = abs(nodes(i+1,j+1,k)%x - inter_x(i,j+1,k)%x)
                                 A(i,j,k) = (a+b)/2*h 
                        end if
                     end do
                  end do
               end do
            case('y')
               do k = 0:dims%kmx
                  do j = 0:dims%jmx
                     do i = 0:dims%imx
                        if((inter_y(i,j,k)%x == nodes(i,j,k)%x) .and. &
                           inter_x(i+1,j,k)%y == nodes(i+1,j,k)%y) then
                                 b = abs(inter_x(i+1,j,k)%x - nodes(i+1,j,k)%x)
                                 h = abs(inter_y(i,j,k)%y - nodes(i,j+1,k)%y)
                                 A(i,j,k) = 1/2*b*h
                        else if((inter_y(i+1,j,k)%x == nodes(i+1,j,k)%x) .and. &
                                 inter_x(i,j,k)%y == nodes(i,j,k)%y)) then
                                 b = abs(inter_x(i,j,k)%x - nodes(i+1,j,k)%x)
                                 h = abs(inter_y(i+1,j,k)%y - nodes(i+1,j,k)%y)
                                 A(i,j,k) = 1/2*b*h
                        else if((inter_y(i,j,k)%x == nodes(i,j,k)%x) .and. &
                                 inter_y(i+1,j,k)%x == nodes(i+1,j,k)%x)) then
                                 h = abs(nodes(i,j,k)%x - nodes(i+1,j,k)%x)
                                 a = abs(nodes(i,j,k)%y - inter_y(i,j,k)%y)
                                 b = abs(nodes(i+1,j,k)%y - inter_y(i+1,j,k)%y)
                              A(i,j,k) = (a+b)/2*h 
                        else if ((inter_x(i,j+1,k)%y == nodes(i,j+1,k)%y) .and. &
                                 inter_x(i,j,k)%y == nodes(i,j,k)%y)) then
                                 h = abs(nodes(i+1,j+1,k)%y - nodes(i+1,j,k)%x)
                                 a = abs(nodes(i+1,j,k)%x - inter_x(i,j,k)%x)
                                 b = abs(nodes(i+1,j+1,k)%x - inter_x(i,j+1,k)%x)
                                 A(i,j,k) = (a+b)/2*h 
                        end if
                     end do
                  end do
               end do
            case('z')
               do k = 0:dims%kmx
                  do j = 0:dims%jmx
                     do i = 0:dims%imx
                        if((inter_y(i,j,k)%x == nodes(i,j,k)%x) .and. &
                           inter_x(i+1,j,k)%y == nodes(i+1,j,k)%y) then
                                 b = abs(inter_x(i+1,j,k)%x - nodes(i+1,j,k)%x)
                                 h = abs(inter_y(i,j,k)%y - nodes(i,j+1,k)%y)
                                 A(i,j,k) = 1/2*b*h
                        else if((inter_y(i+1,j,k)%x == nodes(i+1,j,k)%x) .and. &
                                 inter_x(i,j,k)%y == nodes(i,j,k)%y)) then
                                 b = abs(inter_x(i,j,k)%x - nodes(i+1,j,k)%x)
                                 h = abs(inter_y(i+1,j,k)%y - nodes(i+1,j,k)%y)
                                 A(i,j,k) = 1/2*b*h
                        else if((inter_y(i,j,k)%x == nodes(i,j,k)%x) .and. &
                                 inter_y(i+1,j,k)%x == nodes(i+1,j,k)%x)) then
                                 h = abs(nodes(i,j,k)%x - nodes(i+1,j,k)%x)
                                 a = abs(nodes(i,j,k)%y - inter_y(i,j,k)%y)
                                 b = abs(nodes(i+1,j,k)%y - inter_y(i+1,j,k)%y)
                              A(i,j,k) = (a+b)/2*h 
                        else if ((inter_x(i,j+1,k)%y == nodes(i,j+1,k)%y) .and. &
                                 inter_x(i,j,k)%y == nodes(i,j,k)%y)) then
                                 h = abs(nodes(i+1,j+1,k)%y - nodes(i+1,j,k)%x)
                                 a = abs(nodes(i+1,j,k)%x - inter_x(i,j,k)%x)
                                 b = abs(nodes(i+1,j+1,k)%x - inter_x(i,j+1,k)%x)
                                 A(i,j,k) = (a+b)/2*h 
                        end if
                     end do
                  end do
               end do



         end subroutine wetted_area


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
            integer :: i = 0
            !< Initialiser
            call sign_function(sign_phi, phi_init, dims)
            do while(mag /= 1)
               !!< grad_phi is a vector. Need to make use of qp format as shown
               !!< to extract grad_phi in vector form for other calculations
               call compute_gradient_phi(grad_phi_x, phi, phi_init, cells, &
                                       Ifaces, Jfaces, Kfaces, dims, 'x')
               call compute_gradient_phi(grad_phi_y, phi, phi_init, cells, &
                                       Ifaces, Jfaces, Kfaces, dims, dir'y')
               call compute_gradient_phi(grad_phi_z, phi, phi_init, cells, &
                                       Ifaces, Jfaces, Kfaces, dims, dir'z')
               mag = 1/sqrt(grad_phi_x**2 + grad_phi_y**2 + grad_phi_z**2)
               if (i == 0) then
                  phi = phi_init + del_tau(sign_phi - sign_phi*mag)
                  i = 1
               else
                  phi = phi + del_tau(sign_phi - sign_phi*mag)
               end if
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
            sign_phi(:,:,:) = phi_init(:,:,:)/(sqrt( phi_init(:,:,:)**2 + del_h**2))

         end subroutine sign_function


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
            do k = 0:dims%kmx
               do j = 0:dims%jmx
                  do i = 0:dims%imx
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

         end subroutine smoothen_G