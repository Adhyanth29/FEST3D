module read_vof
!< Module reads the dataset, stores 
!< final data set after program end 

   use vartypes
   ! use mpi
   ! use mapping, only : read_interface_map
    
#include "error.h"
#include "debug.h"
   private

   ! Public methods
   public :: populate_vof

   contains

      subroutine populate_vof(files, vof, dims)
         !< Read the VOF file and initialize VOF
         !-----------------------------------------------------------

         implicit none
         type(filetype), intent(in) :: files
         !< Files' name and handler
         ! type(controltype), intent(in) :: control
         ! !< Control parameters
         real(wp), dimension(:,:,:), allocatable, intent(out) :: vof
         !< To read Volume of Fluid values from file
         type(extent), intent(out) :: dims
         !< Extent of the domain:imx,jmx,kmx
         character(len=STRING_BUFFER_LENGTH) :: line
         !< store read line
         integer :: i, j, k
         integer :: ios  
         !< input/output  status
         
         DebugCall('populate_vof')

         open(files%VOF_FILE_UNIT, file=files%voffile)

         ! call populate_vof(files%VOF_FILE_UNIT, dims)

         ! ! allocate memory for storing grid points
         ! allocate(nodes(-2:dims%imx+3, -2:dims%jmx+3, -2:dims%kmx+3))

         ! !read interface mapping
         ! call read_interface_map(files, control, bc, dims)

         ! ghost grid exchange
         ! call populate_vof(files%VOF_FILE_UNIT, vof, dims)

         ! DebugCall('populate_grid_point')
         !  print *, imx, jmx, kmx

         ! Read grid points from the grid file
         do k = 0, dims%kmx+1
               do j = 0, dims%jmx+1
                  do i = 0, dims%imx+1
                     read(files%VOF_FILE_UNIT, '(A)', iostat=ios) line
                     if (ios /= 0) then
                           print *, 'Error while reading vof line.'
                           print *, 'Current grid point: ', i, j, k
                           print *, 'Current buffer length is set to: ', &
                                    STRING_BUFFER_LENGTH
                           print *, 'Exiting program.'
                           !stop
                     end if
                     !call extract_grid_point(line, i, j, k)
                     read(line, *) vof(i, j, k)
                  end do
               end do
         end do

         close(files%VOF_FILE_UNIT)

         ! ! populate ghost grid points
         ! call ghost_grid(nodes, dims)


      end subroutine populate_vof

      ! subroutine populate_vof(file_handler, dims)
      !    !< Extract the grid size from the grid file header
      !    !
      !    ! We assume that the grid could be in 1 or 2 dimensions. If
      !    ! the grid is in 1 dimension, jmx will be set to 1.
      !    ! We assume that at least one number is specified in the 
      !    ! header, i.e., the grid has atleast one dimension.
      !    !-----------------------------------------------------------

      !    implicit none
      !    integer, intent(in) :: file_handler
      !    !< (input)file handling unit
      !    character(len=STRING_BUFFER_LENGTH) :: header
      !    !< store header
      !    type(extent), intent(out) :: dims
      !    !< Extent of the domain:imx,jmx,kmx
      !    integer :: ios  ! io operation status

      !    DebugCall('populate_vof')

      !    read(file_handler, '(A)', iostat=ios) header
      !    if (ios /= 0) then
      !          print *, 'Error while reading grid file header.'
      !          print *, 'Current buffer length is set to: ', &
      !                STRING_BUFFER_LENGTH
      !          !stop
      !    end if

      !    ! Try to read constants corresponding to two dimensions.
      !    read(header, *, iostat=ios) dims%imx, dims%jmx, dims%kmx
      !    if (ios /= 0) then
      !       print*, "Not able to read dimension from the grid file"
      !       print*, "Make sure you provdie 3D grid"
      !       Fatal_error
      !    end if

      ! end subroutine populate_vof

      ! subroutine populate_vof(file_handler, vof, dims)
      !    !< Use the grid file to populate the grid points.
      !    !-----------------------------------------------------------

      !    implicit none
      !    integer, intent(in) :: file_handler
      !    !< (input)file handling unit
      !    type(extent), intent(in) :: dims
      !    !< Extent of the domain:imx,jmx,kmx
      !    real(wp), dimension(:,:,:), allocatable, intent(out) :: vof
      !    !< To read Volume of Fluid values from file
      !    character(len=STRING_BUFFER_LENGTH) :: line
      !    !< store read line
      !    integer :: i, j, k
      !    integer :: ios  
      !    !< input/output  status

      !    DebugCall('populate_grid_point')
      ! !  print *, imx, jmx, kmx

      !    ! Read grid points from the grid file
      !    do k = 1, dims%kmx
      !          do j = 1, dims%jmx
      !             do i = 1, dims%imx
      !                read(file_handler, '(A)', iostat=ios) line
      !                if (ios /= 0) then
      !                      print *, 'Error while reading grid line.'
      !                      print *, 'Current grid point: ', i, j, k
      !                      print *, 'Current buffer length is set to: ', &
      !                               STRING_BUFFER_LENGTH
      !                      print *, 'Exiting program.'
      !                      !stop
      !                end if
      !                !call extract_grid_point(line, i, j, k)
      !                read(line, *) nodes(i, j, k)%x, nodes(i, j, k)%y, nodes(i, j, k)%z
      !             end do
      !          end do
      !    end do

      ! end subroutine populate_vof
