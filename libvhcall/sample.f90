program main
  use vhcall_fortran
  !use salmon_global, only: natom,kion,izatom,yn_jm
  use structures, only: s_rgrid, s_dft_system
  implicit none
  type(s_rgrid)       :: lg
  type(s_dft_system)  :: system
  !real(8)             :: system_primitive_a(3,3)
  !real(8),allocatable :: system_Rion(:,:)
  integer             :: fp
  character(60)       :: suffix
  character(30)       :: phys_quantity
  !real(8)             :: rmat(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
  real(8)             :: rmat(1:2,1:2,1:2)
  integer             :: comm_is_root
  character(1)        :: yn_jm
  integer             :: natom
  integer,allocatable :: kion(:)
  integer             :: izatom(1)

  integer(8) :: handle
  integer(8) :: sym
  integer(8) :: ca
  integer    :: ir

  lg%num(1) = 2
  lg%num(2) = 2
  lg%num(3) = 2

  lg%is(1) = 1
  lg%ie(1) = 2
  lg%is(2) = 1
  lg%ie(2) = 2
  lg%is(3) = 1
  lg%ie(3) = 2

  allocate(lg%coordinate(1,1:3))
  lg%coordinate(1,1) = 1.0d0
  lg%coordinate(1,2) = 2.0d0
  lg%coordinate(1,3) = 3.0d0

  fp = 100
  suffix = "test"
  phys_quantity = "psi"

  rmat(1,1,1) = 1.0d0
  rmat(2,1,1) = 2.0d0
  rmat(1,2,1) = 3.0d0
  rmat(2,2,1) = 4.0d0
  rmat(1,1,2) = 5.0d0
  rmat(2,1,2) = 6.0d0
  rmat(1,2,2) = 7.0d0
  rmat(2,2,2) = 8.0d0

  system%primitive_a(1,1) = 1.0d0
  system%primitive_a(2,1) = 2.0d0
  system%primitive_a(3,1) = 3.0d0
  system%primitive_a(1,2) = 4.0d0
  system%primitive_a(2,2) = 5.0d0
  system%primitive_a(3,2) = 6.0d0
  system%primitive_a(1,3) = 7.0d0
  system%primitive_a(2,3) = 8.0d0
  system%primitive_a(3,3) = 9.0d0
 
  allocate(system%Rion(1:3,1))
  system%Rion(1,1) = 1.0d0
  system%Rion(2,1) = 2.0d0
  system%Rion(3,1) = 3.0d0

  ! if( comm_is_root(nproc_id_global) ) comm_is_root = 1 else 0
  comm_is_root = 1
  yn_jm = 'n'
  natom = 1

  allocate(kion(1))
  kion(1) = 1

  izatom(1) = 1

write(*,*) "test start !!"

  handle = fvhcall_install('./libvhcall.so')
  sym    = fvhcall_find(handle,'write_cube')
  ca     = fvhcall_args_alloc()

  ir     = fvhcall_args_set(ca, fvhcall_intent_in, 1, lg%is)
  ir     = fvhcall_args_set(ca, fvhcall_intent_in, 2, lg%ie)
  ir     = fvhcall_args_set(ca, fvhcall_intent_in, 3, lg%num)
  ir     = fvhcall_args_set(ca, fvhcall_intent_in, 4, lg%coordinate(1,1))
  ir     = fvhcall_args_set(ca, fvhcall_intent_in, 5, lg%coordinate(1,2))
  ir     = fvhcall_args_set(ca, fvhcall_intent_in, 6, lg%coordinate(1,3))
  ir     = fvhcall_args_set(ca, fvhcall_intent_in, 7, fp)
  ir     = fvhcall_args_set(ca, fvhcall_intent_in, 8, suffix)
  ir     = fvhcall_args_set(ca, fvhcall_intent_in, 9, phys_quantity)
  ir     = fvhcall_args_set(ca, fvhcall_intent_in,10, rmat)
  ir     = fvhcall_args_set(ca, fvhcall_intent_in,11, system%primitive_a)
  ir     = fvhcall_args_set(ca, fvhcall_intent_in,12, system%Rion)
  ir     = fvhcall_args_set(ca, fvhcall_intent_in,13, comm_is_root)
  ir     = fvhcall_args_set(ca, fvhcall_intent_in,14, yn_jm)
  ir     = fvhcall_args_set(ca, fvhcall_intent_in,15, natom)
  ir     = fvhcall_args_set(ca, fvhcall_intent_in,16, kion)
  ir     = fvhcall_args_set(ca, fvhcall_intent_in,17, izatom)
  ir     = fvhcall_invoke_with_args(sym, ca)

  CALL fvhcall_args_free(ca)
  ir = fvhcall_uninstall(handle)

write(*,*) "test end !!"

end program
