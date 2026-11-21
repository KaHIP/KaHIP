!******************************************************************************
! parhip_example.f90
!
! Source of KaHIP -- Karlsruhe High Quality Partitioning.
!
! Partitions a 16 x 16 grid graph with the Fortran module parhip. Every PE owns
! a contiguous range of rows. Run with mpiexec -n 2 (or more PEs, up to 16).
! Stops with a nonzero exit code if a check fails.
!******************************************************************************

program parhip_example
  use, intrinsic :: iso_c_binding, only: c_int, c_double, c_bool
  use mpi
  use parhip
  implicit none

  integer, parameter :: side = 16
  integer :: ierr, rank, nprocs, p, r, c, u, e, first_row, last_row
  integer(parhip_idx), allocatable :: vtxdist(:), xadj(:), adjncy(:), part(:)
  integer(c_int) :: k, edgecut, local_ok, all_ok

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
  if (edgecut <= 0) local_ok = 0
  call MPI_Allreduce(local_ok, all_ok, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)

  if (rank == 0) then
    if (all_ok /= 1) then
      print '(a)', "FAILED: parhip_partition_kway returned an invalid partition"
    else
      print '(a,i0)', "parhip_partition_kway: edge cut ", edgecut
      print '(a)', "all Fortran ParHIP interface checks passed"
    end if
  end if
  call MPI_Finalize(ierr)
  if (all_ok /= 1) error stop 1

contains

  subroutine add(v)
    integer, intent(in) :: v
    e = e + 1
    adjncy(e) = v
  end subroutine add

end program parhip_example
