!******************************************************************************
! parhip.F90
!
! Source of KaHIP -- Karlsruhe High Quality Partitioning.
!
! Fortran interface to the ParHIP library (interface/parhip_interface.h).
!
! The distributed graph is passed in the same format as in the C interface and
! in ParMETIS, with 0-based indices: PE p owns the global nodes
! vtxdist(p+1) .. vtxdist(p+2)-1, xadj and adjncy hold the local CSR adjacency
! with global node IDs, and part receives the block of every local node.
!
! Node weights (vwgt) and edge weights (adjwgt) are optional: leave them out
! (pass the following arguments by keyword) for an unweighted graph.
!
! comm is a Fortran MPI communicator handle, e.g. MPI_COMM_WORLD from the
! module mpi (with mpi_f08, pass MPI_COMM_WORLD%MPI_VAL).
!******************************************************************************

module parhip
  use, intrinsic :: iso_c_binding, only: c_int, c_long_long, c_double, c_bool
  implicit none
  private

  ! idxtype of the C interface (unsigned long long; values must stay below 2**63)
  integer, parameter, public :: parhip_idx = c_long_long

  ! preconfigurations (mode)
  integer(c_int), parameter, public :: PARHIP_ULTRAFASTMESH   = 0
  integer(c_int), parameter, public :: PARHIP_FASTMESH        = 1
  integer(c_int), parameter, public :: PARHIP_ECOMESH         = 2
  integer(c_int), parameter, public :: PARHIP_ULTRAFASTSOCIAL = 3
  integer(c_int), parameter, public :: PARHIP_FASTSOCIAL      = 4
  integer(c_int), parameter, public :: PARHIP_ECOSOCIAL       = 5

  public :: parhip_partition_kway

  interface

    subroutine parhip_partition_kway(vtxdist, xadj, adjncy, vwgt, adjwgt, nparts, imbalance, &
                                     suppress_output, seed, mode, edgecut, part, comm) &
                                     bind(C, name="parhip_fortran_partition_kway")
      import :: c_int, c_double, c_bool, parhip_idx
      integer(parhip_idx), intent(in)           :: vtxdist(*)
      integer(parhip_idx), intent(in)           :: xadj(*)
      integer(parhip_idx), intent(in)           :: adjncy(*)
      integer(parhip_idx), intent(in), optional :: vwgt(*)
      integer(parhip_idx), intent(in), optional :: adjwgt(*)
      integer(c_int),      intent(in)           :: nparts
      real(c_double),      intent(in)           :: imbalance
      logical(c_bool),     value                :: suppress_output
      integer(c_int),      value                :: seed
      integer(c_int),      value                :: mode
      integer(c_int),      intent(out)          :: edgecut
      integer(parhip_idx), intent(out)          :: part(*)
      integer,             intent(in)           :: comm
    end subroutine parhip_partition_kway

  end interface

end module parhip
