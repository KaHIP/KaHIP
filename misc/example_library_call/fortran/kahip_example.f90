!******************************************************************************
! kahip_example.f90
!
! Source of KaHIP -- Karlsruhe High Quality Partitioning.
!
! Calls every routine of the Fortran module kahip on a 10 x 10 grid graph and
! checks the results. Stops with a nonzero exit code if a check fails.
!******************************************************************************

program kahip_example
  use, intrinsic :: iso_c_binding, only: c_int, c_double, c_bool
  use kahip
  implicit none

  integer, parameter :: side = 10
  integer(c_int) :: n, k, i, num_sep, qap, vertexcut
  integer(kahip_idx), allocatable :: xadj(:), adjncy(:), adjcwgt(:)
  integer(c_int), allocatable :: vwgt(:), part(:), separator(:), ordering(:), edge_part(:)
  integer(kahip_idx) :: cut
  integer(c_int) :: hierarchy(2), distance(2)
  real(c_double) :: eps
  logical(c_bool), parameter :: quiet = .true.

  if (kahip_sizeof_idx() /= storage_size(0_kahip_idx) / 8) &
    call fail("kahip_idx does not match the library (64BITMODE?)")

  call grid_graph(side, n, xadj, adjncy)
  allocate(part(n), separator(n), ordering(n), edge_part(xadj(n + 1)))
  k   = 4
  eps = 0.03_c_double

  ! unweighted k-way partition (the optional weights are left out)
  call kaffpa(n, xadj=xadj, adjncy=adjncy, nparts=k, imbalance=eps, suppress_output=quiet, &
              seed=0, mode=KAHIP_STRONG, edgecut=cut, part=part)
  call check_partition("kaffpa", part, k, cut, 1_kahip_idx)
  print '(a,i0)', "kaffpa:                  edge cut ", cut

  ! weighted: every edge has weight 2, so the edge cut doubles
  allocate(vwgt(n), adjcwgt(size(adjncy)))
  vwgt    = 1
  adjcwgt = 2
  call kaffpa(n, vwgt, xadj, adjcwgt, adjncy, k, eps, quiet, 0, KAHIP_ECO, cut, part)
  call check_partition("kaffpa (weighted)", part, k, cut, 2_kahip_idx)
  print '(a,i0)', "kaffpa (weighted):       edge cut ", cut

  call kaffpa_balance(n, xadj=xadj, adjncy=adjncy, nparts=k, imbalance=0.0_c_double, &
                      perfectly_balance=.true._c_bool, suppress_output=quiet, seed=0, &
                      mode=KAHIP_ECO, edgecut=cut, part=part)
  call check_partition("kaffpa_balance", part, k, cut, 1_kahip_idx)
  do i = 0, k - 1
    if (count(part == i) /= n / k) call fail("kaffpa_balance: blocks not perfectly balanced")
  end do
  print '(a,i0)', "kaffpa_balance:          edge cut ", cut

  call kaffpa_connected_blocks(n, xadj=xadj, adjncy=adjncy, nparts=k, imbalance=eps, &
                               suppress_output=quiet, seed=0, mode=KAHIP_STRONG, edgecut=cut, part=part)
  call check_partition("kaffpa_connected_blocks", part, k, cut, 1_kahip_idx)
  do i = 0, k - 1
    if (.not. block_connected(i)) call fail("kaffpa_connected_blocks: a block is not connected")
  end do
  print '(a,i0)', "kaffpa_connected_blocks: edge cut ", cut

  call kaffpa_balance_NE(n, xadj=xadj, adjncy=adjncy, nparts=k, imbalance=eps, &
                         suppress_output=quiet, seed=0, mode=KAHIP_ECO, edgecut=cut, part=part)
  call check_partition("kaffpa_balance_NE", part, k, cut, 1_kahip_idx)
  print '(a,i0)', "kaffpa_balance_NE:       edge cut ", cut

  ! map onto 2 x 2 processors: distance 1 inside a pair, 10 between pairs
  hierarchy = [2, 2]
  distance  = [1, 10]
  call process_mapping(n, xadj=xadj, adjncy=adjncy, hierarchy_parameter=hierarchy, &
                       distance_parameter=distance, hierarchy_depth=2, mode_partitioning=KAHIP_ECO, &
                       mode_mapping=KAHIP_MAPMODE_MULTISECTION, imbalance=eps, suppress_output=quiet, &
                       seed=0, edgecut=cut, qap=qap, part=part)
  call check_partition("process_mapping", part, 4_c_int, cut, 1_kahip_idx)
  if (qap < cut) call fail("process_mapping: qap smaller than edge cut")
  print '(a,i0,a,i0)', "process_mapping:         edge cut ", cut, ", qap ", qap

  call node_separator(n, xadj=xadj, adjncy=adjncy, nparts=2_c_int, imbalance=0.2_c_double, &
                      suppress_output=quiet, seed=0, mode=KAHIP_STRONG, &
                      num_separator_vertices=num_sep, separator=separator)
  call check_separator()
  print '(a,i0)', "node_separator:          separator size ", num_sep

  call reduced_nd(n, xadj, adjncy, quiet, 0, KAHIP_ECO, ordering)
  do i = 0, n - 1
    if (count(ordering == i) /= 1) call fail("reduced_nd: ordering is not a permutation")
  end do
  print '(a)', "reduced_nd:              ordering is a permutation"

  call edge_partitioning(n, xadj=xadj, adjncy=adjncy, nparts=k, imbalance=eps, suppress_output=quiet, &
                         seed=0, mode=KAHIP_FAST, vertexcut=vertexcut, part=edge_part, &
                         infinity_edge_weight=1000_kahip_idx)
  if (any(edge_part < 0) .or. any(edge_part >= k)) call fail("edge_partitioning: block out of range")
  if (vertexcut <= 0) call fail("edge_partitioning: vertex cut should be positive")
  print '(a,i0)', "edge_partitioning:       vertex cut ", vertexcut

  print '(a)', "all Fortran interface checks passed"

contains

  ! side x side grid in 0-based CSR format
  subroutine grid_graph(side, n, xadj, adjncy)
    integer, intent(in) :: side
    integer(c_int), intent(out) :: n
    integer(kahip_idx), allocatable, intent(out) :: xadj(:), adjncy(:)
    integer :: r, c, e, u
    n = side * side
    allocate(xadj(n + 1), adjncy(4 * side * (side - 1)))
    e = 0
    xadj(1) = 0
    do r = 0, side - 1
      do c = 0, side - 1
        u = r * side + c
        if (r > 0)        call push(adjncy, e, u - side)
        if (c > 0)        call push(adjncy, e, u - 1)
        if (c < side - 1) call push(adjncy, e, u + 1)
        if (r < side - 1) call push(adjncy, e, u + side)
        xadj(u + 2) = e
      end do
    end do
  end subroutine grid_graph

  subroutine push(adjncy, e, v)
    integer(kahip_idx), intent(inout) :: adjncy(:)
    integer, intent(inout) :: e
    integer, intent(in) :: v
    e = e + 1
    adjncy(e) = v
  end subroutine push

  ! blocks in range, all used, and the reported cut equals weight * cut edges
  subroutine check_partition(name, part, k, cut, weight)
    character(*), intent(in) :: name
    integer(c_int), intent(in) :: part(:), k
    integer(kahip_idx), intent(in) :: cut, weight
    integer(kahip_idx) :: crossing, j
    integer :: u
    if (any(part < 0) .or. any(part >= k)) call fail(name // ": block out of range")
    do u = 0, k - 1
      if (count(part == u) == 0) call fail(name // ": empty block")
    end do
    crossing = 0
    do u = 1, size(part)
      do j = xadj(u) + 1, xadj(u + 1)
        if (part(u) /= part(adjncy(j) + 1)) crossing = crossing + 1
      end do
    end do
    if (cut /= weight * crossing / 2) call fail(name // ": reported edge cut is wrong")
  end subroutine check_partition

  ! breadth-first search inside block b
  logical function block_connected(b)
    integer(c_int), intent(in) :: b
    block_connected = component_size(b) == count(part == b)
  end function block_connected

  ! removing the separator must disconnect the grid
  subroutine check_separator()
    logical :: is_sep(n)
    integer :: j
    if (num_sep <= 0 .or. num_sep >= n) call fail("node_separator: bad separator size")
    if (any(separator(1:num_sep) < 0) .or. any(separator(1:num_sep) >= n)) &
      call fail("node_separator: node out of range")
    is_sep = .false.
    do j = 1, num_sep
      is_sep(separator(j) + 1) = .true.
    end do
    part = merge(1_c_int, 0_c_int, is_sep)
    if (component_size(0_c_int) == n - num_sep) &
      call fail("node_separator: removing the separator does not disconnect the graph")
  end subroutine check_separator

  ! size of the component of the first node with part == b, using only nodes with part == b
  integer function component_size(b) result(size_found)
    integer(c_int), intent(in) :: b
    integer :: queue(n), head, tail, u, v
    integer(kahip_idx) :: j
    logical :: seen(n)
    seen = .false.
    size_found = 0
    u = findloc_first(b)
    if (u == 0) return
    head = 1
    tail = 1
    queue(1) = u
    seen(u) = .true.
    do while (head <= tail)
      u = queue(head)
      head = head + 1
      do j = xadj(u) + 1, xadj(u + 1)
        v = int(adjncy(j)) + 1
        if (.not. seen(v) .and. part(v) == b) then
          seen(v) = .true.
          tail = tail + 1
          queue(tail) = v
        end if
      end do
    end do
    size_found = tail
  end function component_size

  integer function findloc_first(b)
    integer(c_int), intent(in) :: b
    integer :: u
    findloc_first = 0
    do u = 1, n
      if (part(u) == b) then
        findloc_first = u
        return
      end if
    end do
  end function findloc_first

  subroutine fail(msg)
    character(*), intent(in) :: msg
    print '(a)', "FAILED: " // msg
    error stop 1
  end subroutine fail

end program kahip_example
