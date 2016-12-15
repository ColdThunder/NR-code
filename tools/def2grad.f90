! compile with ifort def2grad.f90 -fpp

program def2grad
implicit none
#include "../dimen.f90"

character(256),parameter::workdir='/data/s2/yuyu22/testmmh/tides30/'
character(8)::zstring='0.000'

integer(4)::step=75
character(32)::idstring='_i0075_'
character(4)::sstring                                                               

character(4)::gstring
character(8)::prefix='.bin'
integer(4),parameter::ng1=NG
integer(4),parameter::ng2=NG
integer(4),parameter::ng3=NG
real(4),allocatable::defp(:,:,:)
real(4),allocatable::grad(:,:,:,:)
character(512)::filename

real(4)::t1,t2

allocate(defp(ng1,ng2,ng3))
allocate(grad(ng1,ng2,ng3,3))
write(gstring,'(i4.4)') NG
write(sstring,'(i4.4)') step
write(*,*) 'GRID ',NG

filename=trim(workdir)//trim(zstring)//'restartdef'//sstring//'-'//gstring//'.bin'
write(*,*) 'reading: ',trim(filename)
open(31,file=filename,status='old',access='stream')
read(31) defp
close(31)

call cpu_time(t1)                                                                 
call scalar2gradient(defp,grad,(/ng1,ng2,ng3/))
filename=trim(workdir)//trim(zstring)//'psix'//trim(idstring)//gstring//trim(prefix)
write(*,*) 'writing: ',trim(filename)
open(31,file=filename,status='replace',access='stream')
write(31) grad(:,:,:,1)
close(31)
filename=trim(workdir)//trim(zstring)//'psiy'//trim(idstring)//gstring//trim(prefix)
write(*,*) 'writing: ',trim(filename)
open(31,file=filename,status='replace',access='stream')
write(31) grad(:,:,:,2)
close(31)
filename=trim(workdir)//trim(zstring)//'psiz'//trim(idstring)//gstring//trim(prefix)
write(*,*) 'writing: ',trim(filename)
open(31,file=filename,status='replace',access='stream')
write(31) grad(:,:,:,3)
close(31)
call cpu_time(t2)
write(*,*) 'time consumed for defp2phi:',t2-t1,'seconds'

deallocate(defp)
deallocate(grad)
endprogram def2grad

include '../gradient.f90'
include '../recipes12.f90'
