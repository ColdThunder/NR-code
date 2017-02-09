! use with 'recipes08.f90'
! 8 Feb 2017, Yu Yu

subroutine regrid(hpos,np,def,ng,newpos)
use omp_lib
implicit none
! io
integer(4)::np
real(4)::hpos(3,np)
integer(4)::ng
real(4)::def(ng,ng,ng)
real(4)::newpos(3,np)
! curlilinear coordinate
real(4)::cpos(3,ng**3)
! nearest particle
integer(4)::pid,sid,gid
integer(4)::p2grid(ng**3) 
integer(4)::sorted(ng**3)
integer(4)::gpro(2,ng**3)
integer(4)::pm
logical(4)::found
real(4)::distance,distancemin
! regrid 
integer(4)::hid
real(4)::dis(3),np_dis(3)
integer(4)::np_gid
real(4)::phixx,phiyy,phizz
! misc
integer(4)::i,j,k,ip,im,jp,jm,kp,km,i0,j0,k0,i1,j1,k1,i2,j2,k2
logical(4)::debug=.false.
real(4)::t1,t2
real(8)::ot1,ot2

write(*,*) ''
write(*,*) 'REGRID'
! obtain the coordinate of each curvilinear frame
call cpu_time(t1)
do k=1,ng
  kp=mod(k,ng)+1
  km=mod(k+ng-2,ng)+1
  do j=1,ng
    jp=mod(j,ng)+1
    jm=mod(j+ng-2,ng)+1
    do i=1,ng
      ip=mod(i,ng)+1
      im=mod(i+ng-2,ng)+1
      gid=(k-1)*ng*ng+(j-1)*ng+i
      cpos(1,gid)=(def(ip,j ,k )-def(im,j ,k ))/2+i-0.5
      cpos(2,gid)=(def(i ,jp,k )-def(i ,jm,k ))/2+j-0.5
      cpos(3,gid)=(def(i ,j ,kp)-def(i ,j ,km))/2+k-0.5
    enddo
  enddo
enddo
where (cpos.gt.ng) cpos=cpos-ng
where (cpos.le.0.) cpos=cpos+ng
call cpu_time(t2)
write(*,*) 'time consumed for curvilinear position:',t2-t1,'seconds'

if (debug) then
  open(31,file='newcor.bin',access='stream',status='replace')
  write(31) cpos
  close(31)
endif

if (debug) write(*,*) 'curvilinear coordinate OK'

call cpu_time(t1)
gpro(:,:)=0
do pid=1,ng**3
  i=ceiling(cpos(1,pid))
  j=ceiling(cpos(2,pid))
  k=ceiling(cpos(3,pid))
  gid=(k-1)*ng*ng+(j-1)*ng+i
  p2grid(pid)=gid
  sorted(pid)=pid
  gpro(1,gid)=gpro(1,gid)+1
enddo
call sort2_int(ng**3,p2grid,sorted)
sid=1
do gid=1,ng**3
  if (gpro(1,gid).eq.0) cycle
  gpro(2,gid)=sid
  sid=sid+gpro(1,gid)
enddo
! gpro(1) store the number of particles belong to this grid
! gpro(2) store the begin location of particles in "sorted"
call cpu_time(t2)
write(*,*) 'time consumed for curvilinear information:',t2-t1,'seconds'

if (debug) write(*,*) 'linklist OK'
if (debug) then
  do sid=1,5
    write(*,*) gpro(1:2,sid)
  enddo
  do sid=ng**3-5,ng**3
    write(*,*) gpro(1:2,sid)
  enddo
endif

call cpu_time(t1)
ot1=omp_get_wtime()
! loop over all particles
do hid=1,np
! find the nearest curvilinear frame for each particle 
  i=ceiling(hpos(1,hid))
  j=ceiling(hpos(2,hid))
  k=ceiling(hpos(3,hid))
  if (debug) write(*,*) 'halo',hid,'belong to',i,j,k
  if (debug) write(*,*) 'halo position',hpos(1:3,hid)
  ! search nearby grid for nearest curvilinear frame
  found=.false.
  distancemin=ng**2/4
  do pm=1,ng/2-1 ! begin with 1 to prevent boundary effect
    if (debug) write(*,*) 'serach range',pm
    do k1=k-pm,k+pm
      k2=mod(k1+ng-1,ng)+1
      do j1=j-pm,j+pm
        j2=mod(j1+ng-1,ng)+1
        do i1=i-pm,i+pm
          i2=mod(i1+ng-1,ng)+1
          gid=(k2-1)*ng*ng+(j2-1)*ng+i2
          if (debug) write(*,*) 'searching',i2,j2,k2
          if (gpro(1,gid).eq.0) cycle
          found=.true.
          if (debug) write(*,*) 'found particles',gpro(1,gid)
          ! check all the particles in grid gid
          do pid=gpro(2,gid),gpro(2,gid)+gpro(1,gid)-1
            if (debug) write(*,*) 'particle',pid,'in gpro has position',cpos(1:3,sorted(pid))
            dis(1:3)=hpos(1:3,hid)-cpos(1:3,sorted(pid))
            if (debug) write(*,*) 'distance 3d',dis
            where (dis.gt.ng/2) dis=dis-ng
            where (dis.lt.-ng/2) dis=dis+ng
            if (debug) write(*,*) 'distance 3d',dis
            distance=sum(dis**2)
            if (debug) write(*,*) 'with distace',sqrt(distance),'distancemin=',sqrt(distancemin)
            if (distance.lt.distancemin) then
              distancemin=distance
              np_dis=dis
              np_gid=sorted(pid)
            endif
          enddo ! end for search in one grid
        enddo
      enddo
    enddo ! end for search range pm
    if (found) exit
  enddo 
  k0=(np_gid-1)/ng**2+1
  j0=(np_gid-(k0-1)*ng*ng-1)/ng+1
  i0=np_gid-(k0-1)*ng*ng-(j0-1)*ng
  if (debug) write(*,*) 'halo',hid,'find',np_gid,'as nearest',i0,j0,k0
  if (debug) write(*,*) 'dis',np_dis
! loop nearby hexahedron to determin which one the particle belongs to
  kp=mod(k0,ng)+1
  km=mod(k0+ng-2,ng)+1
  jp=mod(j0,ng)+1
  jm=mod(j0+ng-2,ng)+1
  ip=mod(i0,ng)+1
  im=mod(i0+ng-2,ng)+1
  phixx=(def(ip,j0,k0)-2*def(i0,j0,k0)+def(im,j0,k0))
  phiyy=(def(i0,jp,k0)-2*def(i0,j0,k0)+def(i0,jm,k0))
  phizz=(def(i0,j0,kp)-2*def(i0,j0,k0)+def(i0,j0,km))
  if (debug) write(*,*) 'phi',phixx,phiyy,phizz
  np_dis(1:3)=np_dis(1:3)/(1.+(/phixx,phiyy,phizz/))
  if (debug) write(*,*) 'np_dis',np_dis
  newpos(1:3,hid)=(/i0,j0,k0/)-(/0.5,0.5,0.5/)+np_dis(1:3)
  if (debug) write(*,*) ''
  if (debug) write(*,*) ''
  if (debug) write(*,*) ''
enddo
where (newpos.gt.ng) newpos=newpos-ng
where (newpos.le.0.) newpos=newpos+ng
call cpu_time(t2)
ot2=omp_get_wtime()
write(*,*) 'cpu time consumed for calculating new position',t2-t1,'seconds'
write(*,*) 'real time consumed for calculating new position',real(ot2-ot1),'seconds'
endsubroutine regrid
