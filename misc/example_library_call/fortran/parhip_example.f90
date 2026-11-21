!******************************************************************************
! parhip_example.f90
!
! Source of KaHIP -- Karlsruhe High Quality Partitioning.
!
! Partitions a 16 x 16 grid graph with the Fortran module parhip. Every PE owns
! a contiguous range of rows. Run with mpiexec -n 2 (or more PEs, up to 16).
! Checks that all blocks are used and that the reported edge cut matches the
! recomputed one; stops with a nonzero exit code if a check fails.
!******************************************************************************

program parhip_example
  use, intrinsic :: iso_c_binding, only: c_int, c_double, c_bool
  use mpi
  use parhip
  implicit none

  integer, parameter :: side = 16
  integer :: ierr, rank, nprocs, p, r, c, u, e, first_row, last_row
  integer(parhip_idx), allocatable :: vtxdist(:), xadj(:), adjncy(:), part(:), global_part(:)
  integer(parhip_idx), allocatable :: vwgt(:), adjwgt(:)
  integer, allocatable :: counts(:), displs(:)
  integer(c_int) :: k, edgecut, weighted_cut, local_ok, all_ok
  integer :: local_crossing, crossing, b, j

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  ! rows are distributed as evenly as possible; vtxdist holds global node IDs (0-based)
  allocate(vtxdist(nprocs + 1))
  do p = 0, nprocs
    vtxdist(p + 1) = int(side * ((p * side) / nprocs), parhip_idx)
  end do
  first_row = int(vtxdist(rank + 1)) / side
  last_row  = int(vtxdist(rank + 2)) / side - 1

  allocate(xadj((last_row - first_row + 1) * side + 1), adjncy(4 * (last_row - first_row + 1) * side))
  e = 0
  xadj(1) = 0
  do r = first_row, last_row
    do c = 0, side - 1
      u = r * side + c
      if (r > 0)        call add(u - side)
      if (c > 0)        call add(u - 1)
      if (c < side - 1) call add(u + 1)
      if (r < side - 1) call add(u + side)
      xadj((r - first_row) * side + c + 2) = e
    end do
  end do

  k = 4
  allocate(part(size(xadj) - 1))
  call parhip_partition_kway(vtxdist, xadj, adjncy, nparts=k, imbalance=0.03_c_double, &
                             suppress_output=.true._c_bool, seed=0, mode=PARHIP_FASTMESH, &
                             edgecut=edgecut, part=part, comm=MPI_COMM_WORLD)

  local_ok = 1
  if (any(part < 0) .or. any(part >= k)) local_ok = 0
  call check_partition(edgecut, 1)

  ! weighted: every edge has weight 2, so the reported cut is twice the number of cut edges
  allocate(vwgt(size(part)), adjwgt(size(adjncy)))
  vwgt   = 1
  adjwgt = 2
  call parhip_partition_kway(vtxdist, xadj, adjncy, vwgt, adjwgt, k, 0.03_c_double, .true._c_bool, &
                             0, PARHIP_FASTMESH, weighted_cut, part, MPI_COMM_WORLD)
  if (any(part < 0) .or. any(part >= k)) local_ok = 0
  call check_partition(weighted_cut, 2)
  call MPI_Allreduce(local_ok, all_ok, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)

  if (rank == 0) then
    if (all_ok /= 1) then
      print '(a)', "FAILED: parhip_partition_kway returned an invalid partition"
    else
      print '(a,i0)', "parhip_partition_kway:            edge cut ", edgecut
      print '(a,i0)', "parhip_partition_kway (weighted): edge cut ", weighted_cut
      print '(a)', "all Fortran ParHIP interface checks passed"
    end if
  end if
  call MPI_Finalize(ierr)
  if (all_ok /= 1) error stop 1

contains

  ! gather the whole partition, recompute the cut and check that all blocks are used
  subroutine check_partition(cut, weight)
    integer(c_int), intent(in) :: cut
    integer, intent(in) :: weight
    if (.not. allocated(global_part)) allocate(counts(nprocs), displs(nprocs), global_part(side * side))
    do p = 0, nprocs - 1
      counts(p + 1) = int(vtxdist(p + 2) - vtxdist(p + 1))
      displs(p + 1) = int(vtxdist(p + 1))
    end do
    call MPI_Allgatherv(part, size(part), MPI_INTEGER8, global_part, counts, displs, MPI_INTEGER8, &
                        MPI_COMM_WORLD, ierr)
    local_crossing = 0
    do u = 1, size(part)
      do j = int(xadj(u)) + 1, int(xadj(u + 1))
        if (part(u) /= global_part(adjncy(j) + 1)) local_crossing = local_crossing + 1
      end do
    end do
    call MPI_Allreduce(local_crossing, crossing, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (cut /= weight * (crossing / 2)) local_ok = 0
    do b = 0, k - 1
      if (count(global_part == b) == 0) local_ok = 0
    end do
  end subroutine check_partition

  subroutine add(v)
    integer, intent(in) :: v
    e = e + 1
    adjncy(e) = v
  end subroutine add

end program parhip_example
