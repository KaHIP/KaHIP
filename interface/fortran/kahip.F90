!******************************************************************************
! kahip.F90
!
! Source of KaHIP -- Karlsruhe High Quality Partitioning.
!
! Fortran interface to the KaHIP library (interface/kaHIP_interface.h).
!
! Graphs are passed in the same compressed (CSR) format as in the C interface
! and in METIS, with 0-based indices: the neighbors of node i (i = 0..n-1) are
! adjncy(xadj(i+1)+1 : xadj(i+2)), and every entry of adjncy is a 0-based node
! ID. Outputs (part, ordering, separator) are 0-based as well.
!
! Node weights (vwgt) and edge weights (adjcwgt) are optional: leave them out
! (pass the following arguments by keyword) for an unweighted graph.
!
! The integer kind kahip_idx matches the kahip_idx type of the C library
! (32 bit, or 64 bit if KaHIP is built with -D64BITMODE=On). The function
! kahip_sizeof_idx() returns the size the library was built with.
!******************************************************************************

module kahip
  use, intrinsic :: iso_c_binding, only: c_int, c_int32_t, c_int64_t, c_double, c_bool
  implicit none
  private

#ifdef KAHIP_64BIT
  integer, parameter, public :: kahip_idx = c_int64_t
#else
  integer, parameter, public :: kahip_idx = c_int32_t
#endif

  ! preconfigurations (mode)
  integer(c_int), parameter, public :: KAHIP_FAST         = 0
  integer(c_int), parameter, public :: KAHIP_ECO          = 1
  integer(c_int), parameter, public :: KAHIP_STRONG       = 2
  integer(c_int), parameter, public :: KAHIP_FASTSOCIAL   = 3
  integer(c_int), parameter, public :: KAHIP_ECOSOCIAL    = 4
  integer(c_int), parameter, public :: KAHIP_STRONGSOCIAL = 5
  integer(c_int), parameter, public :: KAHIP_ULTRA        = 6

  ! process mapping modes
  integer(c_int), parameter, public :: KAHIP_MAPMODE_MULTISECTION = 0
  integer(c_int), parameter, public :: KAHIP_MAPMODE_BISECTION    = 1

  public :: kahip_sizeof_idx
  public :: kaffpa, kaffpa_balance, kaffpa_connected_blocks, kaffpa_balance_NE
  public :: process_mapping, node_separator, reduced_nd, edge_partitioning

  interface

    ! size of kahip_idx in the library in bytes (4 or 8)
    function kahip_sizeof_idx() bind(C, name="kahip_sizeof_idx") result(bytes)
      import :: c_int
      integer(c_int) :: bytes
    end function kahip_sizeof_idx

    ! k-way partition; part(1:n) receives the block (0..nparts-1) of every node
    subroutine kaffpa(n, vwgt, xadj, adjcwgt, adjncy, nparts, imbalance, &
                      suppress_output, seed, mode, edgecut, part) bind(C, name="kaffpa")
      import :: c_int, c_double, c_bool, kahip_idx
      integer(c_int),            intent(in)           :: n
      integer(c_int),            intent(in), optional :: vwgt(*)
      integer(kahip_idx),        intent(in)           :: xadj(*)
      integer(kahip_idx),        intent(in), optional :: adjcwgt(*)
      integer(kahip_idx),        intent(in)           :: adjncy(*)
      integer(c_int),            intent(in)           :: nparts
      real(c_double),            intent(in)           :: imbalance
      logical(c_bool),           value                :: suppress_output
      integer(c_int),            value                :: seed
      integer(c_int),            value                :: mode
      integer(kahip_idx),        intent(out)          :: edgecut
      integer(c_int),            intent(out)          :: part(*)
    end subroutine kaffpa

    ! as kaffpa, with an additional switch for perfect balance
    subroutine kaffpa_balance(n, vwgt, xadj, adjcwgt, adjncy, nparts, imbalance, &
                              perfectly_balance, suppress_output, seed, mode, edgecut, part) &
                              bind(C, name="kaffpa_balance")
      import :: c_int, c_double, c_bool, kahip_idx
      integer(c_int),            intent(in)           :: n
      integer(c_int),            intent(in), optional :: vwgt(*)
      integer(kahip_idx),        intent(in)           :: xadj(*)
      integer(kahip_idx),        intent(in), optional :: adjcwgt(*)
      integer(kahip_idx),        intent(in)           :: adjncy(*)
      integer(c_int),            intent(in)           :: nparts
      real(c_double),            intent(in)           :: imbalance
      logical(c_bool),           value                :: perfectly_balance
      logical(c_bool),           value                :: suppress_output
      integer(c_int),            value                :: seed
      integer(c_int),            value                :: mode
      integer(kahip_idx),        intent(out)          :: edgecut
      integer(c_int),            intent(out)          :: part(*)
    end subroutine kaffpa_balance

    ! as kaffpa, but every block is connected (the input graph must be connected);
    ! only KAHIP_STRONG is supported, for other modes edgecut is set to -1
    subroutine kaffpa_connected_blocks(n, vwgt, xadj, adjcwgt, adjncy, nparts, imbalance, &
                                       suppress_output, seed, mode, edgecut, part) &
                                       bind(C, name="kaffpa_connected_blocks")
      import :: c_int, c_double, c_bool, kahip_idx
      integer(c_int),            intent(in)           :: n
      integer(c_int),            intent(in), optional :: vwgt(*)
      integer(kahip_idx),        intent(in)           :: xadj(*)
      integer(kahip_idx),        intent(in), optional :: adjcwgt(*)
      integer(kahip_idx),        intent(in)           :: adjncy(*)
      integer(c_int),            intent(in)           :: nparts
      real(c_double),            intent(in)           :: imbalance
      logical(c_bool),           value                :: suppress_output
      integer(c_int),            value                :: seed
      integer(c_int),            value                :: mode
      integer(kahip_idx),        intent(out)          :: edgecut
      integer(c_int),            intent(out)          :: part(*)
    end subroutine kaffpa_connected_blocks

    ! as kaffpa, with a balance constraint on nodes and edges
    subroutine kaffpa_balance_NE(n, vwgt, xadj, adjcwgt, adjncy, nparts, imbalance, &
                                 suppress_output, seed, mode, edgecut, part) &
                                 bind(C, name="kaffpa_balance_NE")
      import :: c_int, c_double, c_bool, kahip_idx
      integer(c_int),            intent(in)           :: n
      integer(c_int),            intent(in), optional :: vwgt(*)
      integer(kahip_idx),        intent(in)           :: xadj(*)
      integer(kahip_idx),        intent(in), optional :: adjcwgt(*)
      integer(kahip_idx),        intent(in)           :: adjncy(*)
      integer(c_int),            intent(in)           :: nparts
      real(c_double),            intent(in)           :: imbalance
      logical(c_bool),           value                :: suppress_output
      integer(c_int),            value                :: seed
      integer(c_int),            value                :: mode
      integer(kahip_idx),        intent(out)          :: edgecut
      integer(c_int),            intent(out)          :: part(*)
    end subroutine kaffpa_balance_NE

    ! partition and map the blocks onto a hierarchy of processors
    ! (hierarchy_parameter and distance_parameter have hierarchy_depth entries)
    subroutine process_mapping(n, vwgt, xadj, adjcwgt, adjncy, hierarchy_parameter, &
                               distance_parameter, hierarchy_depth, mode_partitioning, &
                               mode_mapping, imbalance, suppress_output, seed, edgecut, &
                               qap, part) bind(C, name="process_mapping")
      import :: c_int, c_double, c_bool, kahip_idx
      integer(c_int),            intent(in)           :: n
      integer(c_int),            intent(in), optional :: vwgt(*)
      integer(kahip_idx),        intent(in)           :: xadj(*)
      integer(kahip_idx),        intent(in), optional :: adjcwgt(*)
      integer(kahip_idx),        intent(in)           :: adjncy(*)
      integer(c_int),            intent(in)           :: hierarchy_parameter(*)
      integer(c_int),            intent(in)           :: distance_parameter(*)
      integer(c_int),            value                :: hierarchy_depth
      integer(c_int),            value                :: mode_partitioning
      integer(c_int),            value                :: mode_mapping
      real(c_double),            intent(in)           :: imbalance
      logical(c_bool),           value                :: suppress_output
      integer(c_int),            value                :: seed
      integer(kahip_idx),        intent(out)          :: edgecut
      integer(c_int),            intent(out)          :: qap
      integer(c_int),            intent(out)          :: part(*)
    end subroutine process_mapping

    ! node separator; separator must have room for n entries, the first
    ! num_separator_vertices of them receive the separator nodes (0-based)
    subroutine node_separator(n, vwgt, xadj, adjcwgt, adjncy, nparts, imbalance, &
                              suppress_output, seed, mode, num_separator_vertices, separator) &
                              bind(C, name="kahip_fortran_node_separator")
      import :: c_int, c_double, c_bool, kahip_idx
      integer(c_int),            intent(in)           :: n
      integer(c_int),            intent(in), optional :: vwgt(*)
      integer(kahip_idx),        intent(in)           :: xadj(*)
      integer(kahip_idx),        intent(in), optional :: adjcwgt(*)
      integer(kahip_idx),        intent(in)           :: adjncy(*)
      integer(c_int),            intent(in)           :: nparts
      real(c_double),            intent(in)           :: imbalance
      logical(c_bool),           value                :: suppress_output
      integer(c_int),            value                :: seed
      integer(c_int),            value                :: mode
      integer(c_int),            intent(out)          :: num_separator_vertices
      integer(c_int),            intent(out)          :: separator(*)
    end subroutine node_separator

    ! fill-reducing ordering of an unweighted graph (reduced nested dissection);
    ! ordering(1:n) receives the new position of every node
    subroutine reduced_nd(n, xadj, adjncy, suppress_output, seed, mode, ordering) &
                          bind(C, name="reduced_nd")
      import :: c_int, c_bool, kahip_idx
      integer(c_int),            intent(in)           :: n
      integer(kahip_idx),        intent(in)           :: xadj(*)
      integer(kahip_idx),        intent(in)           :: adjncy(*)
      logical(c_bool),           value                :: suppress_output
      integer(c_int),            value                :: seed
      integer(c_int),            value                :: mode
      integer(c_int),            intent(out)          :: ordering(*)
    end subroutine reduced_nd

    ! edge partition; part(1:m) receives the block of every entry of adjncy
    ! (m = xadj(n+1)); the C default for infinity_edge_weight is 1000
    subroutine edge_partitioning(n, vwgt, xadj, adjcwgt, adjncy, nparts, imbalance, &
                                 suppress_output, seed, mode, vertexcut, part, &
                                 infinity_edge_weight) bind(C, name="edge_partitioning")
      import :: c_int, c_double, c_bool, kahip_idx
      integer(c_int),            intent(in)           :: n
      integer(c_int),            intent(in), optional :: vwgt(*)
      integer(kahip_idx),        intent(in)           :: xadj(*)
      integer(kahip_idx),        intent(in), optional :: adjcwgt(*)
      integer(kahip_idx),        intent(in)           :: adjncy(*)
      integer(c_int),            intent(in)           :: nparts
      real(c_double),            intent(in)           :: imbalance
      logical(c_bool),           value                :: suppress_output
      integer(c_int),            value                :: seed
      integer(c_int),            value                :: mode
      integer(c_int),            intent(out)          :: vertexcut
      integer(c_int),            intent(out)          :: part(*)
      integer(kahip_idx),        value                :: infinity_edge_weight
    end subroutine edge_partitioning

  end interface

end module kahip
