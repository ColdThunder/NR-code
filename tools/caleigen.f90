! compile with ifort caleigen.f90 -fpp       

program caleigen
implicit none
#include '../dimen.f90'

character(256),parameter::workdir='/data/s2/yuyu22/testmmh/tides30/'
integer(4),parameter::ng1=NG
integer(4),parameter::ng2=NG
integer(4),parameter::ng3=NG
character(8)::zstring='0.000'
character(4)::gstring
real(4),allocatable::def(:,:,:),den(:,:,:),eigen(:,:,:,:)

integer(4)::step=1500
character(4)::sstring

real(4)::phi(3,3)
real(4)::lambda(3),eigenvector(3,3)
integer(4)::nrot

character(512)::filename
integer(4)::i,j,k,ip,im,kp,km,jp,jm

write(gstring,'(i4.4)') NG
write(sstring,'(i4.4)') step
allocate(def(ng1,ng2,ng3))
allocate(den(ng1,ng2,ng3))
allocate(eigen(ng1,ng2,ng3,3))
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
      phi(1,1)=(def(ip,j,k)-2*def(i,j,k)+def(im,j,k))
      phi(2,2)=(def(i,jp,k)-2*def(i,j,k)+def(i,jm,k))
      phi(3,3)=(def(i,j,kp)-2*def(i,j,k)+def(i,j,km))
      phi(1,2)=(def(ip,jp,k)-def(im,jp,k)-def(ip,jm,k)+def(im,jm,k))/4
      phi(2,3)=(def(i,jp,kp)-def(i,jp,km)-def(i,jm,kp)+def(i,jm,km))/4
      phi(1,3)=(def(ip,j,kp)-def(im,j,kp)-def(ip,j,km)+def(im,j,km))/4
      phi(2,1)=phi(1,2)
      phi(3,2)=phi(2,3)
      phi(3,1)=phi(1,3)
      den(i,j,k)=-(phi(1,1)+phi(2,2)+phi(3,3))
      call jacobi(phi,3,3,lambda,eigenvector,nrot)
      call sort2_real_3real(3,lambda,eigenvector)
      eigen(i,j,k,1:3)=lambda(1:3)
    enddo
  enddo
enddo
! ascend
filename=trim(workdir)//trim(zstring)//'eigen1-'//sstring//'-'//gstring//'.bin'
write(*,*) 'writing: ',trim(filename)
open(31,file=filename,status='replace',access='stream')
write(31) eigen(:,:,:,1)
close(31)

filename=trim(workdir)//trim(zstring)//'eigen2-'//sstring//'-'//gstring//'.bin'
write(*,*) 'writing: ',trim(filename)
open(31,file=filename,status='replace',access='stream')
write(31) eigen(:,:,:,2)
close(31)

filename=trim(workdir)//trim(zstring)//'eigen3-'//sstring//'-'//gstring//'.bin'
write(*,*) 'writing: ',trim(filename)
open(31,file=filename,status='replace',access='stream')
write(31) eigen(:,:,:,3)
close(31)

filename=trim(workdir)//trim(zstring)//'eigen0-'//sstring//'-'//gstring//'.bin'
write(*,*) 'writing: ',trim(filename)
open(31,file=filename,status='replace',access='stream')
write(31) den(:,:,:)
close(31)

deallocate(def)
deallocate(den)
deallocate(eigen)

endprogram caleigen

include 'recipes08.f90'
include 'recipes11.f90'
