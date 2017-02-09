program main
use omp_lib
implicit none
#include "dimen.f90"

character(256),parameter::workdir='/data/s2/yuyu22/testgit/'
! file name variables
character(8),parameter::zstring='0.000'
character(8),parameter::suffix='.bin'

integer(4),parameter::niter=100
character(64),parameter::idstring='_i100_'

! input switcher
logical(4)::inputdensity=.true.
character(64),parameter::infden='den'
logical(4)::inputposition=.true.
character(64),parameter::infpos='pos'
          ! inputposition will overwrite the input density
          ! inputposition is required by regrid process
! restart
logical(4)::restart=.false.
integer(4)::restartstep=100
character(10)::rstring_u='restartu'
character(10)::rstring_d='restartdef'
character(4)::sstring
! output switcher
logical(4)::flag_def=.false. ! output displacement potential at final
logical(4)::flag_gradient=.false. ! output reconstructed displacement field at final
logical(4)::flag_density=.true. ! output reconstructed density at final
logical(4)::flag_check=.true.  ! write down checkpoint after last iteration

! regrid process
logical(4),parameter::regrid_mode=.true.
integer(4),parameter::regridfrequency=10
logical(4)::regridcondition=.false.

! parameters
real(4),parameter::cmpmax=10.
real(4),parameter::dtol=3.
real(4),parameter::dtaumesh=1.
real(4),parameter::cmax=1.
real(4),parameter::dt=1.

! halo variables
integer(4),parameter::npmax=10000000 ! estimated halo number for static memo
integer(4)::np
real(4)::pos(3,npmax)
real(4)::newpos(3,npmax)

! multigrid variables
integer(4),parameter::ng1=NG
integer(4),parameter::ng2=NG
integer(4),parameter::ng3=NG
integer(4),parameter::nfluid=5
real(4)::defp(ng1,ng2,ng3), def(ng1,ng2,ng3),u(nfluid,ng1,ng2,ng3),tmp(ng1,ng2,ng3),tmp2(ng1+2,ng2,ng3)
real(4)::grad(ng1,ng2,ng3,3)
real(4)::delta(ng1,ng2,ng3)
real(4)::umean
! misc
integer(4) i,j,k
character(4)::gstring
character(512)::filename
integer(4)::iter,begin
real(8)::ompt1,ompt2
real(4)::t1,t2
real(4)::rdummy
include 'globalpa.f90'

write(*,*) ' '
write(*,*) 'BEGIN'
if (regrid_mode) then
  write(*,*) ' '
  write(*,*) 'REGRID mode ON!'
  write(*,*) 'regrid every',regridfrequency,'steps'
endif

! grid issues
write(gstring,'(i4.4)') NG
write(*,*) 'grid:',NG

compressmax=cmpmax
u=0.
def=0.
begin=1
! restart
if (restart) then
  begin=restartstep+1
  call readrestart
  if (regrid_mode) call readposition
else
  if (inputdensity) call readdensity
  if (regrid_mode) inputposition=.true.
  if (inputposition) call readposition
  if (inputposition) call voronoi(pos,np,u(1,:,:,:),NG)
endif

write(*,*) ' '
write(*,*) 'BEGIN MAIN LOOP'
! main loop
do iter=begin,niter
  call cpu_time(t1)
  ompt1=omp_get_wtime()
  write(*,*) ' '
  write(*,*) 'iter=',iter

  call calcdefp(defp,tmp,tmp2,def,u,dtol,dtaumesh,nfluid)
  write(*,*) 'CALCDEFP'
  write(*,*) 'defp sigma=',real(sqrt(sum(real(defp**2,8))/float(ng1)/float(ng2)/float(ng3)))

  def=def+defp
  write(*,*) 'DEF+DEFP'
  write(*,*) 'def sigma=',real(sqrt(sum(real(def**2,8))/float(ng1)/float(ng2)/float(ng3)))

  call relaxing(u,def,defp,cmax,dt,nfluid)
  write(*,*) 'RELAXING'
  write(*,*) 'u sigma=',real(sqrt(sum(real((u(1,:,:,:)-1.)**2,8))/float(ng1)/float(ng2)/float(ng3)))
  write(*,*) 'u min max=',minval(u(1,:,:,:)),maxval(u(1,:,:,:))
  write(*,*) 'u < 0 :',count(u(1,:,:,:).lt.0.)

  call cpu_time(t2)
  ompt2=omp_get_wtime()
  write(*,*) 'real time consumed for one iteration:',ompt2-ompt1,'seconds'
  write(*,*) 'cpu time consumed for one iteration:',t2-t1,'seconds'

  ! regrid process
  if (mod(iter,regridfrequency).eq.0) regridcondition=.true.
  if (regrid_mode .and. regridcondition) then
    call regrid(pos,np,def,NG,newpos)
    call voronoi(newpos,np,u(1,:,:,:),NG)
    regridcondition=.false.
  endif
enddo
! end of main loop
if (flag_check) call checkpoint
write(*,*) 'END MAIN LOOP'

write(*,*) ' '
write(*,*) 'def:'
write(*,*) 'avg=',real(sum(real(def,8))/float(ng1)/float(ng2)/float(ng3))
write(*,*) 'sigma=',sqrt(real(sum(real(def**2,8))/float(ng1)/float(ng2)/float(ng3)))

!begin output deformation potential field
if (flag_def) then
  filename=trim(workdir)//trim(zstring)//'def'//trim(idstring)//gstring//trim(suffix)
  write(*,*) 'writing: ',trim(filename)
  open(31,file=filename,status='replace',access='stream')
  write(31) def
  close(31)
endif
!end output deformation potential field

!begin output gradient field
if (flag_gradient) then
  call cpu_time(t1)
  call scalar2gradient(def,grad,(/ng1,ng2,ng3/))

  filename=trim(workdir)//trim(zstring)//'psix'//trim(idstring)//gstring//trim(suffix)
  write(*,*) 'writing: ',trim(filename)
  open(31,file=filename,status='replace',access='stream')
  write(31) grad(:,:,:,1)
  close(31)

  filename=trim(workdir)//trim(zstring)//'psiy'//trim(idstring)//gstring//trim(suffix)
  write(*,*) 'writing: ',trim(filename)
  open(31,file=filename,status='replace',access='stream')
  write(31) grad(:,:,:,2)
  close(31)

  filename=trim(workdir)//trim(zstring)//'psiz'//trim(idstring)//gstring//trim(suffix)
  write(*,*) 'writing: ',trim(filename)
  open(31,file=filename,status='replace',access='stream')
  write(31) grad(:,:,:,3)
  close(31)
  call cpu_time(t2)

  write(*,*) 'time consumed for defp2phi:',t2-t1,'seconds'
endif
! end output gradient field

! begin output reconstructed density field
if (flag_density) then
  call cpu_time(t1)
  call defp2delta(def,delta,(/ng1,ng2,ng3/))

  filename=trim(workdir)//trim(zstring)//'recon'//trim(idstring)//gstring//trim(suffix)
  write(*,*) 'writing: ',trim(filename)
  open(31,file=filename,status='replace',access='stream')
  write(31) delta(:,:,:)
  close(31)

  write(*,*) ' '
  call cpu_time(t2)
  write(*,*) 'time consumed for defp2delta:',t2-t1,'seconds'
  write(*,*) 'delta:'
  write(*,*) 'avg=',real(sum(real(delta(:,:,:),8))/float(ng1)/float(ng2)/float(ng3))
  write(*,*) 'sigma=',real(sqrt(sum(real(delta(:,:,:)**2,8))/float(ng1)/float(ng2)/float(ng3)))
endif
! end output reconstructed density field

contains

subroutine readrestart
implicit none
write(*,*) ' '
write(*,*) 'RESTART!'
write(sstring,'(i4.4)') restartstep
filename=trim(workdir)//trim(zstring)//trim(rstring_u)//sstring//'-'//gstring//trim(suffix)
write(*,*) 'reading: ',trim(filename)
open(31,file=filename,status='old',access='stream')
read(31) u(1,:,:,:)
close(31)
filename=trim(workdir)//trim(zstring)//trim(rstring_d)//sstring//'-'//gstring//trim(suffix)
write(*,*) 'reading: ',trim(filename)
open(31,file=filename,status='old',access='stream')
read(31) def
close(31)
write(*,*) 'def sigma=',real(sqrt(sum(real((def)**2,8))/float(ng1)/float(ng2)/float(ng3)))
write(*,*) 'u sigma=',real(sqrt(sum(real((u(1,:,:,:)-1.)**2,8))/float(ng1)/float(ng2)/float(ng3)))
endsubroutine readrestart

subroutine checkpoint
implicit none
write(*,*) ' '
write(*,*) 'CHECKPOINT!'
write(sstring,'(i4.4)') niter
filename=trim(workdir)//trim(zstring)//trim(rstring_u)//sstring//'-'//gstring//trim(suffix)
write(*,*) 'reading: ',trim(filename)
open(31,file=filename,status='replace',access='stream')
write(31) u(1,:,:,:)
close(31)
filename=trim(workdir)//trim(zstring)//trim(rstring_d)//sstring//'-'//gstring//trim(suffix)
write(*,*) 'reading: ',trim(filename)
open(31,file=filename,status='replace',access='stream')
write(31) def
close(31)
endsubroutine checkpoint

subroutine readdensity
implicit none
! input density field with mean of 1
write(*,*) ' '
call cpu_time(t1)
filename=trim(workdir)//trim(zstring)//trim(infden)//gstring//trim(suffix)
write(*,*) 'reading: ',trim(filename)
open(31,file=filename,status='old',access='stream')
read(31) u(1,:,:,:)
close(31)
call cpu_time(t2)
write(*,*) 'time consumed for reading:',t2-t1,'seconds'

write(*,*) ' '
write(*,*) 'u:'
umean=sum(real(u(1,:,:,:),8))/float(ng1)/float(ng2)/float(ng3)
!call calcmean(u(1,:,:,:),umean,ng1,ng2,ng3,'u')
if (abs(umean-1.).gt.1.e-5) then
  write(*,*) 'umean should be zero, program ended'
  stop
endif
write(*,*) 'u mean=',real(sum(real(u(1,:,:,:),8))/float(ng1)/float(ng2)/float(ng3))
write(*,*) 'u sigma=',real(sqrt(sum(real((u(1,:,:,:)-1.)**2,8))/float(ng1)/float(ng2)/float(ng3)))
write(*,*) 'u min, max=',minval(u(1,:,:,:)),maxval(u(1,:,:,:))
endsubroutine readdensity

subroutine readposition
implicit none
real(4)::box(6)
call cpu_time(t1)
write(*,*) ' '
filename=trim(workdir)//trim(zstring)//trim(infpos)//gstring//trim(suffix)
write(*,*) 'reading: ',trim(filename)
open(31,file=filename,status='old',access='stream')
read(31) np
write(*,*) 'particle number:',np
read(31) box(1:6)
read(31) pos(1:3,1:np)
close(31)
pos(1,:)=(pos(1,:)-box(1))/(box(2)-box(1))*ng1
pos(2,:)=(pos(2,:)-box(3))/(box(4)-box(3))*ng2
pos(3,:)=(pos(3,:)-box(5))/(box(6)-box(5))*ng3
where(pos.eq.0.) pos=pos+NG
endsubroutine readposition

endprogram main
