subroutine write_cube(lg_is,lg_ie,lg_num,lg_coordinate1,lg_coordinate2,lg_coordinate3,fp,suffix,phys_quantity,&
rmat,system_primitive_a,system_Rion,comm_is_root,yn_jm,natom,kion,izatom)
  !use salmon_global, only: natom,kion,izatom,yn_jm
  !use structures, only: s_rgrid, s_dft_system
  implicit none
  !type(s_rgrid)     ,intent(in) :: lg
  integer,dimension(3),intent(in) :: lg_is,lg_ie,lg_num
  real(8)           ,intent(in) :: lg_coordinate1,lg_coordinate2,lg_coordinate3
  !type(s_dft_system),intent(in) :: system
  real(8),intent(in)            :: system_primitive_a(3,3),system_Rion(1:3,*)
  integer           ,intent(in) :: fp
  character(60)     ,intent(in) :: suffix
  character(30)     ,intent(in) :: phys_quantity
  !real(8)           ,intent(in) :: rmat(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
  real(8)           ,intent(in) :: rmat(lg_is(1):lg_ie(1),lg_is(2):lg_ie(2),lg_is(3):lg_ie(3))
  integer           ,intent(in) :: comm_is_root
  character(1)      ,intent(in) :: yn_jm
  integer           ,intent(in) :: natom
  integer           ,intent(in) :: kion(*)
  integer           ,intent(in) :: izatom(*)
  !
  character(60):: filename
  integer :: j,iatom
  integer :: ix,iy,iz
  integer :: ik
  real(8) :: daa(3,3)

  do j=1,3
     !daa(1:3,j) = system%primitive_a(1:3,j)/dble(lg%num(j))
     daa(1:3,j) = system_primitive_a(1:3,j)/dble(lg_num(j))
  end do

  !if(comm_is_root(nproc_id_global))then
  if(comm_is_root .ne. 0)then
    filename=trim(suffix)//".cube"
    open(fp,file=filename)
    if(phys_quantity=="psi")then
      write(fp,*) "Molecular Orbital"
    else if(phys_quantity=="dns")then
      write(fp,*) "Electron Density"
    else if(phys_quantity=="pbcd")then
      write(fp,*) "Positive Background Charge Density"
    else if(phys_quantity=="elf")then
      write(fp,*) "Electron Localization Function"
    else if(phys_quantity=="dnsdiff")then
      write(fp,*) "Difference of Electron Density"
    else if(phys_quantity=="exsta")then
      write(fp,*) "x Component of Static Electric Field"
    else if(phys_quantity=="eysta")then
      write(fp,*) "y Component of Static Electric Field"
    else if(phys_quantity=="ezsta")then
      write(fp,*) "z Component of Static Electric Field"
    end if
    write(fp,*) "All values here are in a.u."
    !write(fp,'(i5,3f12.6)') natom,lg%coordinate(lg%is(1),1),lg%coordinate(lg%is(2),2),lg%coordinate(lg%is(3),3)
    !write(fp,'(i5,3f12.6)') lg%num(1),daa(1:3,1)
    !write(fp,'(i5,3f12.6)') lg%num(2),daa(1:3,2)
    !write(fp,'(i5,3f12.6)') lg%num(3),daa(1:3,3)
    write(fp,'(i5,3f12.6)') natom,lg_coordinate1,lg_coordinate2,lg_coordinate3
    write(fp,'(i5,3f12.6)') lg_num(1),daa(1:3,1)
    write(fp,'(i5,3f12.6)') lg_num(2),daa(1:3,2)
    write(fp,'(i5,3f12.6)') lg_num(3),daa(1:3,3)
    if(yn_jm=='n')then
      do iatom=1,natom
        ik=Kion(iatom)
        !write(fp,'(i5,4f12.6)') izatom(ik),dble(izatom(ik)),(system%Rion(j,iatom),j=1,3)
        write(fp,'(i5,4f12.6)') izatom(ik),dble(izatom(ik)),(system_Rion(j,iatom),j=1,3)
      end do
    else
      write(fp,'(i5,4f12.6)') 1,1.0d0,0.0d0,0.0d0,0.0d0
    end if

    !do ix=lg%is(1),lg%ie(1)
    !do iy=lg%is(2),lg%ie(2)
    !  write(fp,'(6(1X,E23.15E3))', advance="yes") (rmat(ix,iy,iz),iz=lg%is(3),lg%ie(3))
    do ix=lg_is(1),lg_ie(1)
    do iy=lg_is(2),lg_ie(2)
      write(fp,'(6(1X,E23.15E3))', advance="yes") (rmat(ix,iy,iz),iz=lg_is(3),lg_ie(3))
!    do iz=lg%is(3),lg%ie(3)
!      if(mod(iz+1-lg%is(3),6)==0)then
!        write(fp,'(e13.5)', advance="yes") abs(rmat(ix,iy,iz))
!      else
!        write(fp,'(e13.5)', advance="no") abs(rmat(ix,iy,iz))
!      endif
!    end do
!    write(fp,*)
    end do
    end do
    close(fp)
  end if

end subroutine write_cube
