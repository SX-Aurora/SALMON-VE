!
!  Copyright 2019-2020 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!

subroutine zstencil_typical_seq(is_array,ie_array,is,ie,idx,idy,idz,igs,ige &
                               ,tpsi,htpsi,V_local,lap0,lapt,nabt &
                               )
  !$acc routine worker
  implicit none

  integer,intent(in) :: is_array(3),ie_array(3),is(3),ie(3)
  integer,intent(in) :: idx(is(1)-4:ie(1)+4),idy(is(2)-4:ie(2)+4),idz(is(3)-4:ie(3)+4)
  integer,intent(in) :: igs(3),ige(3)

  complex(8),intent(in)  :: tpsi   (is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  complex(8),intent(out) :: htpsi  (is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
  real(8),   intent(in)  :: V_local(is(1):ie(1),is(2):ie(2),is(3):ie(3))
  real(8),   intent(in)  :: lap0
  real(8),   intent(in)  :: lapt(12), nabt(12)

  complex(8), parameter :: zI=(0.d0,1.d0)

  integer    :: ix,iy,iz
  complex(8) :: v,w
  complex(8) :: t(8)
#ifdef __ve__
  integer, parameter :: blksz = 256
  logical, save :: first = .true.
  logical, save :: collapse = .true.
  complex(8), allocatable, save :: tmp(:,:,:,:)
  integer, save :: xs, xe, xl, ys, ye, yl, zs, ze, zl, ll, ll_

  if(first) then
    xs = lbound(tpsi,1)
    xe = ubound(tpsi,1)
    xl = xe - xs + 1
    ys = lbound(tpsi,2)
    ye = ubound(tpsi,2)
    yl = ye - ys + 1
    zs = lbound(tpsi,3)
    ze = ubound(tpsi,3)
    zl = ze - zs + 1

    allocate(tmp(xs:xe,ys:ye,zs:ze,2))

    if((ige(1)-igs(1)+1 > 256) .or. (ige(2)-igs(2)+1 > 256)) collapse = .false.
    if(is_array(1) /= is(1)) collapse = .false.
    if(ie_array(1) /= ie(1)) collapse = .false.
    if(is_array(2) /= is(2)) collapse = .false.
    if(ie_array(2) /= ie(2)) collapse = .false.
    if(is_array(1) /= xs) collapse = .false.
    if(ie_array(1) /= xe) collapse = .false.
    if(is_array(2) /= ys) collapse = .false.
    if(ie_array(2) /= ye) collapse = .false.

    first = .false.
  end if
#endif

#ifdef __INTEL_COMPILER
#if defined(__KNC__) || defined(__AVX512F__)
#   define MEM_ALIGN   64
#   define VECTOR_SIZE 4
# else
#   define MEM_ALIGN   32
#   define VECTOR_SIZE 2
# endif

!dir$ assume_aligned V_local:MEM_ALIGN
!dir$ assume_aligned tpsi   :MEM_ALIGN
!dir$ assume_aligned htpsi  :MEM_ALIGN
#endif

#define DX(dt) idx(ix+(dt)),iy,iz
#define DY(dt) ix,idy(iy+(dt)),iz
#define DZ(dt) ix,iy,idz(iz+(dt))

#ifdef __ve__
  if(.not. collapse) then
#endif
  !$acc loop collapse(3)
  do iz=igs(3),ige(3)
  do iy=igs(2),ige(2)

!dir$ assume_aligned V_local(is(1),iy,iz):MEM_ALIGN
!dir$ assume_aligned tpsi(is_array(1),iy,iz)   :MEM_ALIGN
!dir$ assume_aligned htpsi(is_array(1),iy,iz)  :MEM_ALIGN

  do ix=igs(1),ige(1)
    t(1) = tpsi(DX( 4))
    t(2) = tpsi(DX( 3))
    t(3) = tpsi(DX( 2))
    t(4) = tpsi(DX( 1))
    t(5) = tpsi(DX(-1))
    t(6) = tpsi(DX(-2))
    t(7) = tpsi(DX(-3))
    t(8) = tpsi(DX(-4))

    v=(lapt(1)*(t(4)+t(5)) &
    & +lapt(2)*(t(3)+t(6)) &
    & +lapt(3)*(t(2)+t(7)) &
    & +lapt(4)*(t(1)+t(8)))
    w=(nabt(1)*(t(4)-t(5)) &
    & +nabt(2)*(t(3)-t(6)) &
    & +nabt(3)*(t(2)-t(7)) &
    & +nabt(4)*(t(1)-t(8)))

    t(1) = tpsi(DY( 1))
    t(2) = tpsi(DY( 2))
    t(3) = tpsi(DY( 3))
    t(4) = tpsi(DY( 4))
    t(5) = tpsi(DY(-1))
    t(6) = tpsi(DY(-2))
    t(7) = tpsi(DY(-3))
    t(8) = tpsi(DY(-4))

    v=(lapt(5)*(t(1)+t(5)) &
    & +lapt(6)*(t(2)+t(6)) &
    & +lapt(7)*(t(3)+t(7)) &
    & +lapt(8)*(t(4)+t(8))) + v
    w=(nabt(5)*(t(1)-t(5)) &
    & +nabt(6)*(t(2)-t(6)) &
    & +nabt(7)*(t(3)-t(7)) &
    & +nabt(8)*(t(4)-t(8))) + w

    t(1) = tpsi(DZ( 1))
    t(2) = tpsi(DZ( 2))
    t(3) = tpsi(DZ( 3))
    t(4) = tpsi(DZ( 4))
    t(5) = tpsi(DZ(-1))
    t(6) = tpsi(DZ(-2))
    t(7) = tpsi(DZ(-3))
    t(8) = tpsi(DZ(-4))

    v=(lapt( 9)*(t(1)+t(5)) &
    & +lapt(10)*(t(2)+t(6)) &
    & +lapt(11)*(t(3)+t(7)) &
    & +lapt(12)*(t(4)+t(8))) + v
    w=(nabt( 9)*(t(1)-t(5)) &
    & +nabt(10)*(t(2)-t(6)) &
    & +nabt(11)*(t(3)-t(7)) &
    & +nabt(12)*(t(4)-t(8))) + w

    htpsi(ix,iy,iz) = V_local(ix,iy,iz)*tpsi(ix,iy,iz) &
                    + lap0*tpsi(ix,iy,iz) &
                    - 0.5d0 * v - zI * w
  end do

  end do
  end do
#ifdef __ve__
  else
! call ftrace_region_begin('x-direc')
  do ll_=1,yl*zl,blksz
!NEC$ outerloop_unroll(16)
!NEC$ shortloop
  do ix=igs(1),ige(1)
!NEC$ nointerchange
!NEC$ shortloop
  do ll=ll_,min(ll_+blksz-1,yl*zl)
    t(1) = tpsi(idx(ix+4),ys+ll-1,zs)
    t(2) = tpsi(idx(ix+3),ys+ll-1,zs)
    t(3) = tpsi(idx(ix+2),ys+ll-1,zs)
    t(4) = tpsi(idx(ix+1),ys+ll-1,zs)
    t(5) = tpsi(idx(ix-1),ys+ll-1,zs)
    t(6) = tpsi(idx(ix-2),ys+ll-1,zs)
    t(7) = tpsi(idx(ix-3),ys+ll-1,zs)
    t(8) = tpsi(idx(ix-4),ys+ll-1,zs)

    v=(lapt(1)*(t(4)+t(5)) &
    & +lapt(2)*(t(3)+t(6)) &
    & +lapt(3)*(t(2)+t(7)) &
    & +lapt(4)*(t(1)+t(8)))
    w=(nabt(1)*(t(4)-t(5)) &
    & +nabt(2)*(t(3)-t(6)) &
    & +nabt(3)*(t(2)-t(7)) &
    & +nabt(4)*(t(1)-t(8)))

    tmp(ix,ys+ll-1,zs,1) = - 0.5d0 * v - zI * w
  end do
  end do
  end do
! call ftrace_region_end('x-direc')

! call ftrace_region_begin('y-direc')
  do ll_=1,xl*zl,blksz
!NEC$ outerloop_unroll(16)
!NEC$ shortloop
  do iy=igs(2),ige(2)
!NEC$ shortloop
!NEC$ gather_reorder
!NEC$ vovertake
!NEC$ vob
!NEC$ ivdep
  do ll=ll_,min(ll_+blksz-1,xl*zl)
    ix=mod(ll-1,xl)+1
    iz=(ll-1)/xl+1
    t(1) = tpsi(DY( 1))
    t(2) = tpsi(DY( 2))
    t(3) = tpsi(DY( 3))
    t(4) = tpsi(DY( 4))
    t(5) = tpsi(DY(-1))
    t(6) = tpsi(DY(-2))
    t(7) = tpsi(DY(-3))
    t(8) = tpsi(DY(-4))

    v=(lapt(5)*(t(1)+t(5)) &
    & +lapt(6)*(t(2)+t(6)) &
    & +lapt(7)*(t(3)+t(7)) &
    & +lapt(8)*(t(4)+t(8)))
    w=(nabt(5)*(t(1)-t(5)) &
    & +nabt(6)*(t(2)-t(6)) &
    & +nabt(7)*(t(3)-t(7)) &
    & +nabt(8)*(t(4)-t(8)))

    tmp(ix,iy,iz,2) = - 0.5d0 * v - zI * w
  end do
  end do
  end do
! call ftrace_region_end('y-direc')

! call ftrace_region_begin('z-direc')
!NEC$ novector
  do ll_=1,xl*yl,blksz
!NEC$ outerloop_unroll(4)
  do iz=igs(3),ige(3)
!NEC$ shortloop
  do ll=ll_,min(ll_+blksz-1,xl*yl)
    t(1) = tpsi(xs+ll-1,ys,idz(iz+1))
    t(2) = tpsi(xs+ll-1,ys,idz(iz+2))
    t(3) = tpsi(xs+ll-1,ys,idz(iz+3))
    t(4) = tpsi(xs+ll-1,ys,idz(iz+4))
    t(5) = tpsi(xs+ll-1,ys,idz(iz-1))
    t(6) = tpsi(xs+ll-1,ys,idz(iz-2))
    t(7) = tpsi(xs+ll-1,ys,idz(iz-3))
    t(8) = tpsi(xs+ll-1,ys,idz(iz-4))

    v=(lapt( 9)*(t(1)+t(5)) &
    & +lapt(10)*(t(2)+t(6)) &
    & +lapt(11)*(t(3)+t(7)) &
    & +lapt(12)*(t(4)+t(8)))
    w=(nabt( 9)*(t(1)-t(5)) &
    & +nabt(10)*(t(2)-t(6)) &
    & +nabt(11)*(t(3)-t(7)) &
    & +nabt(12)*(t(4)-t(8)))

    htpsi(xs+ll-1,ys,iz) = V_local(xs+ll-1,ys,iz)*tpsi(xs+ll-1,ys,iz) &
                    + lap0*tpsi(xs+ll-1,ys,iz) &
                    - 0.5d0 * v - zI * w + tmp(xs+ll-1,ys,iz,1) + tmp(xs+ll-1,ys,iz,2)
  end do
  end do
  end do
! call ftrace_region_end('z-direc')
  end if
#endif
end subroutine
