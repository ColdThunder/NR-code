!-*- Fortran -*-	

subroutine mg64(u,rhs,defpot,ubig,dx,nu,iopt)
implicit none
integer iopt,nu
real u(64,64,64),rhs(64,64,64),defpot(64,64,64),ubig(nu,64,64,64),dx
!dir$ shared *u(:block,:block,:),*defpot(:block,:block,:)
!dir$ shared *rhs(:block,:block,:)      
! locals
integer isolver,imgalg,irelaxopt,nrelax,i
#define N2 64/2
real uhalf(N2,N2,N2),rhshalf(N2,N2,N2),defhalf(N2,N2,N2),usmall(1,N2,N2,N2)
!dir$ shared uhalf(:block,:block,:),defhalf(:block,:block,:)
!dir$ shared rhshalf(:block,:block,:),usmall(:,:block,:block,:)
real relaxerr,r1

isolver=mod(iopt,8)

#ifdef _ALPHA
! first touch placement of the subarrays
!cdec$ MIGRATE_NEXT_TOUCH_NOPRESERVE(uhalf,defhalf,rhshalf,usmall)
call szero(uhalf,N2,N2,N2)
call szero(defhalf,N2,N2,N2)
call szero(rhshalf,N2,N2,N2)
call szero(usmall,N2,N2,N2)
#endif
nrelax=4*NG/64
irelaxopt=isolver*8+3
! calculate residual on half grid      
call relax(u,rhshalf,rhs,defpot,ubig,r1,dx,64,64,64,nu,nrelax,irelaxopt) 

! goto 59
call restrict(defhalf,defpot,64,64,64)
call restrictu(usmall,ubig,64,64,64,nu)
call szero(uhalf,N2,N2,N2)

call mg32(uhalf,rhshalf,defhalf,usmall,dx*2,1,iopt)
      
call inject(u,uhalf,64,64,64)
59    continue
irelaxopt=isolver*8+2
call relax(u,rhshalf,rhs,defpot,ubig,relaxerr,dx,64,64,64,nu,nrelax,irelaxopt) 
if (64 .eq. NG) then
   write(*,*) 'relaxerr=',relaxerr
endif
if (64 .ne. NG .and. relaxerr .gt. r1 ) then
  write(*,73) 64,r1,relaxerr
73      format(' mg: relaxation diverged, level=',i3,' r1,2=',2e12.4)
! prolong resets u: forget about the current level relaxation!
!  write(*,*) 'call prolong 64'
!  call prolong(u,uhalf,64,64,64)
  call szero(u,64,64,64)
else
  call zerosum(u,64,64,64)
endif
return
endsubroutine mg64
