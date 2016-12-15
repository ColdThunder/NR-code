! compile with ifort caldisp.f90 -fpp       

program caleigen
implicit none
#include '../dimen.f90'

character(256),parameter::workdir='/data/s2/yuyu22/testmmh/tides30/'
integer(4),parameter::ng1=NG
integer(4),parameter::ng2=NG
integer(4),parameter::ng3=NG
character(8)::zstring='0.000'
character(4)::gstring
real(4),allocatable::def(:,:,:),grad(:,:,:,:)

integer(4)::step=1500
character(4)::sstring

character(512)::filename
integer(4)::i,j,k,ip,im,kp,km,jp,jm

write(gstring,'(i4.4)') NG
write(sstring,'(i4.4)') step
allocate(def(ng1,ng2,ng3))
allocate(grad(ng1,ng2,ng3,3))
! read in
filename=trim(workdir)//trim(zstring)//'restartdef'//sstring//'-'//gstring//'.bin'
write(*,*) 'reading: ',trim(filename)
open(31,file=filename,status='old',access='stream')
read(31) def
close(31)
! 
do k=1,ng3
  write(*,*) 'k=',k
  kp=mod(k,ng3)+1
  km=mod(k+ng3-2,ng3)+1
  do j=1,ng2
    jp=mod(j,ng2)+1
    jm=mod(j+ng2-2,ng2)+1
    do i=1,ng1
      ip=mod(i,ng1)+1
      im=mod(i+ng1-2,ng1)+1

      grad(i,j,k,1)=-(def(ip,j ,k )-def(im,j ,k ))/2
      grad(i,j,k,2)=-(def(i ,jp,k )-def(i ,jm,k ))/2
      grad(i,j,k,3)=-(def(i ,j ,kp)-def(i ,j ,km))/2

    enddo
  enddo
enddo

filename=trim(workdir)//trim(zstring)//'grad1-'//sstring//'-'//gstring//'.bin'
write(*,*) 'writing: ',trim(filename)
open(31,file=filename,status='replace',access='stream')
write(31) grad(:,:,:,1)
close(31)

filename=trim(workdir)//trim(zstring)//'grad2-'//sstring//'-'//gstring//'.bin'
write(*,*) 'writing: ',trim(filename)
open(31,file=filename,status='replace',access='stream')
write(31) grad(:,:,:,2)
close(31)

filename=trim(workdir)//trim(zstring)//'grad3-'//sstring//'-'//gstring//'.bin'
write(*,*) 'writing: ',trim(filename)
open(31,file=filename,status='replace',access='stream')
write(31) grad(:,:,:,3)
close(31)

deallocate(def)
deallocate(grad)

endprogram caleigen

include 'recipes08.f90'
include 'recipes11.f90'
