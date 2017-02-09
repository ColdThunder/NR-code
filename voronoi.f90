! use with 'recipes08.f90'
! 8 Feb 2017, Yu Yu

subroutine voronoi(pos,np,den,ng)
use omp_lib
implicit none
! io
integer(4)::np
real(4)::pos(3,np) ! (0,ng]
integer(4)::ng
real(4)::den(ng,ng,ng)
! varible kind
! kind of p2grid,pid should cover ng^3
! kind of gpro,sorted,gid,sid should be same as np 
! grid list
integer(4)::p2grid(np) ! 
integer(4)::sorted(np) ! store particle id in order of grid index
integer(4)::gpro(2,ng**3) ! information for array "sorted"
integer(4)::gid,sid,pid
! nearest particle
integer(4)::dentag(ng,ng,ng) ! record the nearest particle id for grid
integer(4)::ptag(np) ! store how many grid treat this particle as nearest
real(4)::distance,distancemin
integer(4)::partlist  ! should be same kind with np
real(4)::gpos(3),ppos(3)
integer(4)::pm,i1,j1,k1,i2,j2,k2
logical(4)::found
! assignmass
real(4)::dsum,mean,sigma,mmin,mmax
! misc
integer(4)::i,j,k
real(4)::rdummy
logical(4)::debug=.false.

real(4)::t1,t2
write(*,*) ''
write(*,*) 'Voronoi mass assignment'
! produce grid list
if (debug) write(*,*) 'producing grid particle list'
gpro(:,:)=0
do pid=1,np
  i1=ceiling(pos(1,pid))
  j1=ceiling(pos(2,pid))
  k1=ceiling(pos(3,pid))
!  i1=mod(i1+ng-1,ng)+1   ! for safety
!  j1=mod(j1+ng-1,ng)+1
!  k1=mod(k1+ng-1,ng)+1
  gid=(k1-1)*ng*ng+(j1-1)*ng+i1
  p2grid(pid)=gid ! "pid"-th particle belong to grid "gid"
  sorted(pid)=pid ! a trival integer array for particle id 1-np
  gpro(1,gid)=gpro(1,gid)+1 ! number count for grid "gid"
enddo
call sort2_int(np,p2grid,sorted)
sid=1
do gid=1,ng**3
  if (gpro(1,gid).eq.0) cycle
  gpro(2,gid)=sid ! begin location for grid "gid" in "sorted"
  sid=sid+gpro(1,gid)
enddo
! gpro(1) store the number of particles belong to this grid
! gpro(2) store the begin location of particles in "sorted"

! being nearest particle
if (debug) write(*,*) 'begin finding nearest particle'
ptag=0
dentag=0
!$omp parallel do default(private) shared(debug,ng,gpro,pos,sorted,dentag,ptag)
do k=1,ng
  if (debug) write(*,*) 'calculating k=',k
  call cpu_time(t1)
  do j=1,ng
    do i=1,ng
      gpos=((/i,j,k/)-0.5)
      distancemin=ng**2/4.
      found=.false.
      do pm=0,ng/2-1
        do k1=k-pm,k+pm
          k2=mod(k1+ng-1,ng)+1
          do j1=j-pm,j+pm
            j2=mod(j1+ng-1,ng)+1
            do i1=i-pm,i+pm
              i2=mod(i1+ng-1,ng)+1
              gid=(k2-1)*ng*ng+(j2-1)*ng+i2
              if (gpro(1,gid).eq.0) cycle
              found=.true.
              ! check all the particles in grid gid
              do pid=gpro(2,gid),gpro(2,gid)+gpro(1,gid)-1
                ppos(1:3)=pos(1:3,sorted(pid))-gpos(1:3)
                where(ppos.gt.ng/2) ppos=ppos-ng
                where(ppos.lt.-ng/2) ppos=ppos+ng
                distance=sum(ppos**2)
                if (distance.lt.distancemin) then
                  distancemin=distance
                  partlist=sorted(pid)
                endif
              enddo
              ! end check all the particles in grid id
            enddo
          enddo
        enddo
        if (found) then
          ! ptag is the number of grid theat this particle as nearest particle
          ! dentag record the nearest particle to this grid
          !$omp atomic
          ptag(partlist)=ptag(partlist)+1
          dentag(i,j,k)=partlist
           exit
        endif
        ! end search for gird (i,j,k)
      enddo
    enddo
  enddo
  call cpu_time(t2)
  if (debug) write(*,*) 'TIME:',t2-t1,'seconds'
enddo
!$omp end parallel do

! begin assignmass
if (debug) write(*,*) 'assigning mass'
den=0.
! assignment mass to gird for untagged particles 
do pid=1,np
  if (ptag(pid).ne.0) cycle
  i1=ceiling(pos(1,pid))
  j1=ceiling(pos(2,pid))
  k1=ceiling(pos(3,pid))
  den(i1,j1,k1)=den(i1,j1,k1)+1.
enddo
! find assigned mass for the grid from tagged particles
do k=1,ng
  do j=1,ng
    do i=1,ng
      den(i,j,k)=den(i,j,k)+1./ptag(dentag(i,j,k))
    enddo
  enddo
enddo
dsum=sum(real(den,8))
write(*,*) 'sum:',dsum
den=den/sum(den)*float(ng)**3
dsum=sum(real(den,8))
mean=dsum/float(ng)**3
write(*,*) 'mean:',mean
sigma=sqrt(sum((real(den,8)-mean)**2)/float(ng)**3)
write(*,*) 'sigma:',sigma
mmin=minval(den)
mmax=maxval(den)
write(*,*) 'min,max:',mmin,mmax
write(*,*) 'number of den=min:',count(real(den).eq.mmin)
endsubroutine voronoi
