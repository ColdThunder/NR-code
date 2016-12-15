! recipes for variations of sorting
! SORT2_REAL: sort the real array, meanwhile changing another real array
! SORT2_INT: sort the integer array, meanwhile changing another integer array
! SORT2_REAL_INT: sort the real array, meanwhile changing the integer array
! SORT2_INT_NREAL: sort the integer array, meanwhile changing the n column real array
! SORT2_INT_REAL: call SORT2_INT_NREAL with n=1
! SORT2_INT_3REAL: call SORT2_INT_NREAL with n=3
! SORT2_REAL_NREAL: sort the real array, meanwhile changing the n column real array
! SORT2_REAL_REAL: call SORT2_REAL_NREAL with n=1
! SORT2_REAL_3REAL: call SORT2_REAL_NREAL with n=3
! SORT_INT: sort the integer array
! SORT_REAL: sort the real array
! SORT_REAL_TABLE: sort the 2 dimensional real array according to a specifid column
! last modified: Mar 2 2016

SUBROUTINE sort2_real(n,arr,brr) 
INTEGER n,M,NSTACK
REAL arr(n),brr(n) 
PARAMETER (M=7,NSTACK=50)
INTEGER i,ir,j,jstack,k,l,istack(NSTACK) 
REAL a,b,temp

jstack=0
l=1
ir=n
1  if(ir-l.lt.M)then
  do j=l+1,ir 
    a=arr(j) 
    b=brr(j)
    do i=j-1,l,-1
      if(arr(i).le.a)goto 2 
      arr(i+1)=arr(i) 
      brr(i+1)=brr(i)
    enddo
    i=l-1 
2   arr(i+1)=a 
    brr(i+1)=b
  enddo
  if(jstack.eq.0)return 
  ir=istack(jstack) 
  l=istack(jstack-1) 
  jstack=jstack-2
else 
  k=(l+ir)/2
  temp=arr(k) 
  arr(k)=arr(l+1) 
  arr(l+1)=temp
  temp=brr(k) 
  brr(k)=brr(l+1) 
  brr(l+1)=temp 
  if(arr(l).gt.arr(ir))then
    temp=arr(l) 
    arr(l)=arr(ir) 
    arr(ir)=temp 
    temp=brr(l) 
    brr(l)=brr(ir) 
    brr(ir)=temp
  endif 
  if(arr(l+1).gt.arr(ir))then
    temp=arr(l+1) 
    arr(l+1)=arr(ir) 
    arr(ir)=temp 
    temp=brr(l+1) 
    brr(l+1)=brr(ir) 
    brr(ir)=temp
  endif 
  if(arr(l).gt.arr(l+1))then
    temp=arr(l) 
    arr(l)=arr(l+1) 
    arr(l+1)=temp 
    temp=brr(l) 
    brr(l)=brr(l+1) 
    brr(l+1)=temp
  endif 
  i=l+1
  j=ir 
  a=arr(l+1) 
  b=brr(l+1)
3 continue 
  i=i+1
  if(arr(i).lt.a)goto 3
4 continue
  j=j-1 
  if(arr(j).gt.a)goto 4 
  if(j.lt.i)goto 5 
  temp=arr(i) 
  arr(i)=arr(j) 
  arr(j)=temp 
  temp=brr(i) 
  brr(i)=brr(j) 
  brr(j)=temp
  goto 3
5 arr(l+1)=arr(j) 
  arr(j)=a
  brr(l+1)=brr(j) 
  brr(j)=b 
  jstack=jstack+2
  if(jstack.gt.NSTACK)pause'NSTACK too small in sort2' 
  if(ir-i+1.ge.j-l)then
    istack(jstack)=ir 
    istack(jstack-1)=i 
    ir=j-1
  else 
    istack(jstack)=j-1 
    istack(jstack-1)=l 
    l=i
  endif 
endif
goto 1 
ENDSUBROUTINE sort2_real

SUBROUTINE sort2_int(n,arr,brr) 
INTEGER n,M,NSTACK
INTEGER arr(n),brr(n) 
PARAMETER (M=7,NSTACK=50)
INTEGER i,ir,j,jstack,k,l,istack(NSTACK) 
INTEGER a,b,temp

jstack=0
l=1
ir=n
1  if(ir-l.lt.M)then
  do j=l+1,ir 
    a=arr(j) 
    b=brr(j)
    do i=j-1,l,-1
      if(arr(i).le.a)goto 2 
      arr(i+1)=arr(i) 
      brr(i+1)=brr(i)
    enddo
    i=l-1 
2   arr(i+1)=a 
    brr(i+1)=b
  enddo
  if(jstack.eq.0)return 
  ir=istack(jstack) 
  l=istack(jstack-1) 
  jstack=jstack-2
else 
  k=(l+ir)/2
  temp=arr(k) 
  arr(k)=arr(l+1) 
  arr(l+1)=temp
  temp=brr(k) 
  brr(k)=brr(l+1) 
  brr(l+1)=temp 
  if(arr(l).gt.arr(ir))then
    temp=arr(l) 
    arr(l)=arr(ir) 
    arr(ir)=temp 
    temp=brr(l) 
    brr(l)=brr(ir) 
    brr(ir)=temp
  endif 
  if(arr(l+1).gt.arr(ir))then
    temp=arr(l+1) 
    arr(l+1)=arr(ir) 
    arr(ir)=temp 
    temp=brr(l+1) 
    brr(l+1)=brr(ir) 
    brr(ir)=temp
  endif 
  if(arr(l).gt.arr(l+1))then
    temp=arr(l) 
    arr(l)=arr(l+1) 
    arr(l+1)=temp 
    temp=brr(l) 
    brr(l)=brr(l+1) 
    brr(l+1)=temp
  endif 
  i=l+1
  j=ir 
  a=arr(l+1) 
  b=brr(l+1)
3 continue 
  i=i+1
  if(arr(i).lt.a)goto 3
4 continue
  j=j-1 
  if(arr(j).gt.a)goto 4 
  if(j.lt.i)goto 5 
  temp=arr(i) 
  arr(i)=arr(j) 
  arr(j)=temp 
  temp=brr(i) 
  brr(i)=brr(j) 
  brr(j)=temp
  goto 3
5 arr(l+1)=arr(j) 
  arr(j)=a
  brr(l+1)=brr(j) 
  brr(j)=b 
  jstack=jstack+2
  if(jstack.gt.NSTACK)pause'NSTACK too small in sort2' 
  if(ir-i+1.ge.j-l)then
    istack(jstack)=ir 
    istack(jstack-1)=i 
    ir=j-1
  else 
    istack(jstack)=j-1 
    istack(jstack-1)=l 
    l=i
  endif 
endif
goto 1 
ENDSUBROUTINE sort2_int

SUBROUTINE sort2_real_int(n,arr,brr) 
INTEGER n,M,NSTACK
REAL arr(n)
INTEGER brr(n) 
PARAMETER (M=7,NSTACK=50)
INTEGER i,ir,j,jstack,k,l,istack(NSTACK) 
REAL a,atemp
INTEGER b,btemp

jstack=0
l=1
ir=n
1  if(ir-l.lt.M)then
  do j=l+1,ir 
    a=arr(j) 
    b=brr(j)
    do i=j-1,l,-1
      if(arr(i).le.a)goto 2 
      arr(i+1)=arr(i) 
      brr(i+1)=brr(i)
    enddo
    i=l-1 
2   arr(i+1)=a 
    brr(i+1)=b
  enddo
  if(jstack.eq.0)return 
  ir=istack(jstack) 
  l=istack(jstack-1) 
  jstack=jstack-2
else 
  k=(l+ir)/2
  atemp=arr(k) 
  arr(k)=arr(l+1) 
  arr(l+1)=atemp
  btemp=brr(k) 
  brr(k)=brr(l+1) 
  brr(l+1)=btemp 
  if(arr(l).gt.arr(ir))then
    atemp=arr(l) 
    arr(l)=arr(ir) 
    arr(ir)=atemp 
    btemp=brr(l) 
    brr(l)=brr(ir) 
    brr(ir)=btemp
  endif 
  if(arr(l+1).gt.arr(ir))then
    atemp=arr(l+1) 
    arr(l+1)=arr(ir) 
    arr(ir)=atemp 
    btemp=brr(l+1) 
    brr(l+1)=brr(ir) 
    brr(ir)=btemp
  endif 
  if(arr(l).gt.arr(l+1))then
    atemp=arr(l) 
    arr(l)=arr(l+1) 
    arr(l+1)=atemp 
    btemp=brr(l) 
    brr(l)=brr(l+1) 
    brr(l+1)=btemp
  endif 
  i=l+1
  j=ir 
  a=arr(l+1) 
  b=brr(l+1)
3 continue 
  i=i+1
  if(arr(i).lt.a)goto 3
4 continue
  j=j-1 
  if(arr(j).gt.a)goto 4 
  if(j.lt.i)goto 5 
  atemp=arr(i) 
  arr(i)=arr(j) 
  arr(j)=atemp 
  btemp=brr(i) 
  brr(i)=brr(j) 
  brr(j)=btemp
  goto 3
5 arr(l+1)=arr(j) 
  arr(j)=a
  brr(l+1)=brr(j) 
  brr(j)=b 
  jstack=jstack+2
  if(jstack.gt.NSTACK)pause'NSTACK too small in sort2' 
  if(ir-i+1.ge.j-l)then
    istack(jstack)=ir 
    istack(jstack-1)=i 
    ir=j-1
  else 
    istack(jstack)=j-1 
    istack(jstack-1)=l 
    l=i
  endif 
endif
goto 1 
ENDSUBROUTINE sort2_real_int

SUBROUTINE sort2_int_nreal(n,arr,brr,nc) 
INTEGER n,M,NSTACK
INTEGER arr(n)
REAL brr(nc,n) 
PARAMETER (M=7,NSTACK=50)
INTEGER i,ir,j,jstack,k,l,istack(NSTACK) 
INTEGER a,atemp
REAL b(nc),btemp(nc)

jstack=0
l=1
ir=n
1  if(ir-l.lt.M)then
  do j=l+1,ir 
    a=arr(j) 
    b(1:nc)=brr(1:nc,j)
    do i=j-1,l,-1
      if(arr(i).le.a)goto 2 
      arr(i+1)=arr(i) 
      brr(1:nc,i+1)=brr(1:nc,i)
    enddo
    i=l-1 
2   arr(i+1)=a 
    brr(1:nc,i+1)=b(1:nc)
  enddo
  if(jstack.eq.0)return 
  ir=istack(jstack) 
  l=istack(jstack-1) 
  jstack=jstack-2
else 
  k=(l+ir)/2
  atemp=arr(k) 
  arr(k)=arr(l+1) 
  arr(l+1)=atemp
  btemp(1:nc)=brr(1:nc,k) 
  brr(1:nc,k)=brr(1:nc,l+1) 
  brr(1:nc,l+1)=btemp(1:nc)
  if(arr(l).gt.arr(ir))then
    atemp=arr(l) 
    arr(l)=arr(ir) 
    arr(ir)=atemp 
    btemp(1:nc)=brr(1:nc,l) 
    brr(1:nc,l)=brr(1:nc,ir) 
    brr(1:nc,ir)=btemp(1:nc)
  endif 
  if(arr(l+1).gt.arr(ir))then
    atemp=arr(l+1) 
    arr(l+1)=arr(ir) 
    arr(ir)=atemp 
    btemp(1:nc)=brr(1:nc,l+1) 
    brr(1:nc,l+1)=brr(1:nc,ir) 
    brr(1:nc,ir)=btemp(1:nc)
  endif 
  if(arr(l).gt.arr(l+1))then
    atemp=arr(l) 
    arr(l)=arr(l+1) 
    arr(l+1)=atemp 
    btemp(1:nc)=brr(1:nc,l) 
    brr(1:nc,l)=brr(1:nc,l+1) 
    brr(1:nc,l+1)=btemp(1:nc)
  endif 
  i=l+1
  j=ir 
  a=arr(l+1) 
  b(1:nc)=brr(1:nc,l+1)
3 continue 
  i=i+1
  if(arr(i).lt.a)goto 3
4 continue
  j=j-1 
  if(arr(j).gt.a)goto 4 
  if(j.lt.i)goto 5 
  atemp=arr(i) 
  arr(i)=arr(j) 
  arr(j)=atemp 
  btemp(1:nc)=brr(1:nc,i) 
  brr(1:nc,i)=brr(1:nc,j) 
  brr(1:nc,j)=btemp(1:nc)
  goto 3
5 arr(l+1)=arr(j) 
  arr(j)=a
  brr(1:nc,l+1)=brr(1:nc,j) 
  brr(1:nc,j)=b(1:nc)
  jstack=jstack+2
  if(jstack.gt.NSTACK)pause'NSTACK too small in sort2' 
  if(ir-i+1.ge.j-l)then
    istack(jstack)=ir 
    istack(jstack-1)=i 
    ir=j-1
  else 
    istack(jstack)=j-1 
    istack(jstack-1)=l 
    l=i
  endif 
endif
goto 1 
ENDSUBROUTINE sort2_int_nreal

SUBROUTINE sort2_int_real(n,arr,brr)
INTEGER n
INTEGER arr(n)
REAL brr(1,n)
call sort2_int_nreal(n,arr,brr,1)
ENDSUBROUTINE sort2_int_real

SUBROUTINE sort2_int_3real(n,arr,brr)
INTEGER n
INTEGER arr(n)
REAL brr(3,n)
call sort2_int_nreal(n,arr,brr,3)
ENDSUBROUTINE sort2_int_3real

SUBROUTINE sort2_real_nreal(n,arr,brr,nc) 
INTEGER n,M,NSTACK
REAL arr(n)
REAL brr(nc,n) 
PARAMETER (M=7,NSTACK=50)
INTEGER i,ir,j,jstack,k,l,istack(NSTACK) 
REAL a,atemp
REAL b(nc),btemp(nc)

jstack=0
l=1
ir=n
1  if(ir-l.lt.M)then
  do j=l+1,ir 
    a=arr(j) 
    b(1:nc)=brr(1:nc,j)
    do i=j-1,l,-1
      if(arr(i).le.a)goto 2 
      arr(i+1)=arr(i) 
      brr(1:nc,i+1)=brr(1:nc,i)
    enddo
    i=l-1 
2   arr(i+1)=a 
    brr(1:nc,i+1)=b(1:nc)
  enddo
  if(jstack.eq.0)return 
  ir=istack(jstack) 
  l=istack(jstack-1) 
  jstack=jstack-2
else 
  k=(l+ir)/2
  atemp=arr(k) 
  arr(k)=arr(l+1) 
  arr(l+1)=atemp
  btemp(1:nc)=brr(1:nc,k) 
  brr(1:nc,k)=brr(1:nc,l+1) 
  brr(1:nc,l+1)=btemp(1:nc)
  if(arr(l).gt.arr(ir))then
    atemp=arr(l) 
    arr(l)=arr(ir) 
    arr(ir)=atemp 
    btemp(1:nc)=brr(1:nc,l) 
    brr(1:nc,l)=brr(1:nc,ir) 
    brr(1:nc,ir)=btemp(1:nc)
  endif 
  if(arr(l+1).gt.arr(ir))then
    atemp=arr(l+1) 
    arr(l+1)=arr(ir) 
    arr(ir)=atemp 
    btemp(1:nc)=brr(1:nc,l+1) 
    brr(1:nc,l+1)=brr(1:nc,ir) 
    brr(1:nc,ir)=btemp(1:nc)
  endif 
  if(arr(l).gt.arr(l+1))then
    atemp=arr(l) 
    arr(l)=arr(l+1) 
    arr(l+1)=atemp 
    btemp(1:nc)=brr(1:nc,l) 
    brr(1:nc,l)=brr(1:nc,l+1) 
    brr(1:nc,l+1)=btemp(1:nc)
  endif 
  i=l+1
  j=ir 
  a=arr(l+1) 
  b(1:nc)=brr(1:nc,l+1)
3 continue 
  i=i+1
  if(arr(i).lt.a)goto 3
4 continue
  j=j-1 
  if(arr(j).gt.a)goto 4 
  if(j.lt.i)goto 5 
  atemp=arr(i) 
  arr(i)=arr(j) 
  arr(j)=atemp 
  btemp(1:nc)=brr(1:nc,i) 
  brr(1:nc,i)=brr(1:nc,j) 
  brr(1:nc,j)=btemp(1:nc)
  goto 3
5 arr(l+1)=arr(j) 
  arr(j)=a
  brr(1:nc,l+1)=brr(1:nc,j) 
  brr(1:nc,j)=b(1:nc)
  jstack=jstack+2
  if(jstack.gt.NSTACK)pause'NSTACK too small in sort2' 
  if(ir-i+1.ge.j-l)then
    istack(jstack)=ir 
    istack(jstack-1)=i 
    ir=j-1
  else 
    istack(jstack)=j-1 
    istack(jstack-1)=l 
    l=i
  endif 
endif
goto 1 
ENDSUBROUTINE sort2_real_nreal

SUBROUTINE sort2_real_real(n,arr,brr)
INTEGER n
REAL arr(n)
REAL brr(1,n)
call sort2_real_nreal(n,arr,brr,1)
ENDSUBROUTINE sort2_real_real

SUBROUTINE sort2_real_3real(n,arr,brr)
INTEGER n
REAL arr(n)
REAL brr(3,n)
call sort2_real_nreal(n,arr,brr,3)
ENDSUBROUTINE sort2_real_3real

SUBROUTINE sort_int(n,arr) 
INTEGER n,M,NSTACK
INTEGER arr(n) 
PARAMETER (M=7,NSTACK=50)
INTEGER i,ir,j,jstack,k,l,istack(NSTACK) 
INTEGER a,temp

jstack=0
l=1
ir=n
1  if(ir-l.lt.M)then
  do j=l+1,ir 
    a=arr(j) 
    do i=j-1,l,-1
      if(arr(i).le.a)goto 2 
      arr(i+1)=arr(i) 
    enddo
    i=l-1 
2   arr(i+1)=a 
  enddo
  if(jstack.eq.0)return 
  ir=istack(jstack) 
  l=istack(jstack-1) 
  jstack=jstack-2
else 
  k=(l+ir)/2
  temp=arr(k) 
  arr(k)=arr(l+1) 
  arr(l+1)=temp
  if(arr(l).gt.arr(ir))then
    temp=arr(l) 
    arr(l)=arr(ir) 
    arr(ir)=temp 
  endif 
  if(arr(l+1).gt.arr(ir))then
    temp=arr(l+1) 
    arr(l+1)=arr(ir) 
    arr(ir)=temp 
  endif 
  if(arr(l).gt.arr(l+1))then
    temp=arr(l) 
    arr(l)=arr(l+1) 
    arr(l+1)=temp 
  endif 
  i=l+1
  j=ir 
  a=arr(l+1) 
3 continue 
  i=i+1
  if(arr(i).lt.a)goto 3
4 continue
  j=j-1 
  if(arr(j).gt.a)goto 4 
  if(j.lt.i)goto 5 
  temp=arr(i) 
  arr(i)=arr(j) 
  arr(j)=temp 
  goto 3
5 arr(l+1)=arr(j) 
  arr(j)=a
  jstack=jstack+2
  if(jstack.gt.NSTACK)pause'NSTACK too small in sort2' 
  if(ir-i+1.ge.j-l)then
    istack(jstack)=ir 
    istack(jstack-1)=i 
    ir=j-1
  else 
    istack(jstack)=j-1 
    istack(jstack-1)=l 
    l=i
  endif 
endif
goto 1 
ENDSUBROUTINE sort_int

SUBROUTINE sort_real(n,arr) 
INTEGER n,M,NSTACK
REAL arr(n) 
PARAMETER (M=7,NSTACK=50)
INTEGER i,ir,j,jstack,k,l,istack(NSTACK) 
REAL a,temp

jstack=0
l=1
ir=n
1  if(ir-l.lt.M)then
  do j=l+1,ir 
    a=arr(j) 
    do i=j-1,l,-1
      if(arr(i).le.a)goto 2 
      arr(i+1)=arr(i) 
    enddo
    i=l-1 
2   arr(i+1)=a 
  enddo
  if(jstack.eq.0)return 
  ir=istack(jstack) 
  l=istack(jstack-1) 
  jstack=jstack-2
else 
  k=(l+ir)/2
  temp=arr(k) 
  arr(k)=arr(l+1) 
  arr(l+1)=temp
  if(arr(l).gt.arr(ir))then
    temp=arr(l) 
    arr(l)=arr(ir) 
    arr(ir)=temp 
  endif 
  if(arr(l+1).gt.arr(ir))then
    temp=arr(l+1) 
    arr(l+1)=arr(ir) 
    arr(ir)=temp 
  endif 
  if(arr(l).gt.arr(l+1))then
    temp=arr(l) 
    arr(l)=arr(l+1) 
    arr(l+1)=temp 
  endif 
  i=l+1
  j=ir 
  a=arr(l+1) 
3 continue 
  i=i+1
  if(arr(i).lt.a)goto 3
4 continue
  j=j-1 
  if(arr(j).gt.a)goto 4 
  if(j.lt.i)goto 5 
  temp=arr(i) 
  arr(i)=arr(j) 
  arr(j)=temp 
  goto 3
5 arr(l+1)=arr(j) 
  arr(j)=a
  jstack=jstack+2
  if(jstack.gt.NSTACK)pause'NSTACK too small in sort2' 
  if(ir-i+1.ge.j-l)then
    istack(jstack)=ir 
    istack(jstack-1)=i 
    ir=j-1
  else 
    istack(jstack)=j-1 
    istack(jstack-1)=l 
    l=i
  endif 
endif
goto 1 
ENDSUBROUTINE sort_real

SUBROUTINE sort_real_table(nc,n,brr,isort) 
INTEGER n,M,NSTACK
REAL brr(nc,n),arr(n) 
PARAMETER (M=7,NSTACK=50)
INTEGER i,ir,j,jstack,k,l,istack(NSTACK) 
REAL b(nc),tempc(nc),a,temp

arr(1:n)=brr(isort,1:n)

jstack=0
l=1
ir=n
1  if(ir-l.lt.M)then
  do j=l+1,ir 
    a=arr(j) 
    b(1:nc)=brr(1:nc,j) 
    do i=j-1,l,-1
      if(arr(i).le.a)goto 2 
      arr(i+1)=arr(i) 
      brr(1:nc,i+1)=brr(1:nc,i) 
    enddo
    i=l-1 
2   arr(i+1)=a 
    brr(1:nc,i+1)=b(1:nc)
  enddo
  if(jstack.eq.0)return 
  ir=istack(jstack) 
  l=istack(jstack-1) 
  jstack=jstack-2
else 
  k=(l+ir)/2
  temp=arr(k) 
  arr(k)=arr(l+1) 
  arr(l+1)=temp
  tempc(1:nc)=brr(1:nc,k) 
  brr(1:nc,k)=brr(1:nc,l+1) 
  brr(1:nc,l+1)=tempc(1:nc)
  if(arr(l).gt.arr(ir))then
    temp=arr(l) 
    arr(l)=arr(ir) 
    arr(ir)=temp 
    tempc(1:nc)=brr(1:nc,l) 
    brr(1:nc,l)=brr(1:nc,ir) 
    brr(1:nc,ir)=tempc(1:nc) 
  endif 
  if(arr(l+1).gt.arr(ir))then
    temp=arr(l+1) 
    arr(l+1)=arr(ir) 
    arr(ir)=temp 
    tempc(1:nc)=brr(1:nc,l+1) 
    brr(1:nc,l+1)=brr(1:nc,ir) 
    brr(1:nc,ir)=tempc(1:nc) 
  endif 
  if(arr(l).gt.arr(l+1))then
    temp=arr(l) 
    arr(l)=arr(l+1) 
    arr(l+1)=temp 
    tempc(1:nc)=brr(1:nc,l) 
    brr(1:nc,l)=brr(1:nc,l+1) 
    brr(1:nc,l+1)=tempc(1:nc) 
  endif 
  i=l+1
  j=ir 
  a=arr(l+1) 
  b(1:nc)=brr(1:nc,l+1) 
3 continue 
  i=i+1
  if(arr(i).lt.a)goto 3
4 continue
  j=j-1 
  if(arr(j).gt.a)goto 4 
  if(j.lt.i)goto 5 
  temp=arr(i) 
  arr(i)=arr(j) 
  arr(j)=temp 
  tempc(1:nc)=brr(1:nc,i) 
  brr(1:nc,i)=brr(1:nc,j) 
  brr(1:nc,j)=tempc(1:nc) 
  goto 3
5 arr(l+1)=arr(j) 
  arr(j)=a
  brr(1:nc,l+1)=brr(1:nc,j) 
  brr(1:nc,j)=b(1:nc)
  jstack=jstack+2
  if(jstack.gt.NSTACK)pause'NSTACK too small in sort2' 
  if(ir-i+1.ge.j-l)then
    istack(jstack)=ir 
    istack(jstack-1)=i 
    ir=j-1
  else 
    istack(jstack)=j-1 
    istack(jstack-1)=l 
    l=i
  endif 
endif
goto 1 
ENDSUBROUTINE sort_real_table

