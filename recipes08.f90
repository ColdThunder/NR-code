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
END

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
END

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
END

SUBROUTINE sort2_int_real(n,arr,brr) 
INTEGER n,M,NSTACK
INTEGER arr(n)
REAL brr(n) 
PARAMETER (M=7,NSTACK=50)
INTEGER i,ir,j,jstack,k,l,istack(NSTACK) 
INTEGER a,atemp
REAL b,btemp

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
END


SUBROUTINE reverse_real(n,arr)
INTEGER n
REAL arr(n),brr(n)
INTEGER i
do i=1,n
  brr(i)=arr(n+1-i)
enddo
arr=brr
END

SUBROUTINE reverse_int(n,arr)
INTEGER n
INTEGER arr(n),brr(n)
INTEGER i
do i=1,n
  brr(i)=arr(n+1-i)
enddo
arr=brr
END


SUBROUTINE sort2_real_3real(n,arr,brr) 
INTEGER n,M,NSTACK
REAL arr(n)
REAL brr(3,n) 
PARAMETER (M=7,NSTACK=50)
INTEGER i,ir,j,jstack,k,l,istack(NSTACK) 
REAL a,atemp
real b(3),btemp(3)

jstack=0
l=1
ir=n
1  if(ir-l.lt.M)then
  do j=l+1,ir 
    a=arr(j) 
    b(:)=brr(:,j)
    do i=j-1,l,-1
      if(arr(i).le.a)goto 2 
      arr(i+1)=arr(i) 
      brr(:,i+1)=brr(:,i)
    enddo
    i=l-1 
2   arr(i+1)=a 
    brr(:,i+1)=b(:)
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
  btemp(:)=brr(:,k) 
  brr(:,k)=brr(:,l+1) 
  brr(:,l+1)=btemp(:) 
  if(arr(l).gt.arr(ir))then
    atemp=arr(l) 
    arr(l)=arr(ir) 
    arr(ir)=atemp 
    btemp(:)=brr(:,l) 
    brr(:,l)=brr(:,ir) 
    brr(:,ir)=btemp(:)
  endif 
  if(arr(l+1).gt.arr(ir))then
    atemp=arr(l+1) 
    arr(l+1)=arr(ir) 
    arr(ir)=atemp 
    btemp(:)=brr(:,l+1) 
    brr(:,l+1)=brr(:,ir) 
    brr(:,ir)=btemp(:)
  endif 
  if(arr(l).gt.arr(l+1))then
    atemp=arr(l) 
    arr(l)=arr(l+1) 
    arr(l+1)=atemp 
    btemp(:)=brr(:,l) 
    brr(:,l)=brr(:,l+1) 
    brr(:,l+1)=btemp(:)
  endif 
  i=l+1
  j=ir 
  a=arr(l+1) 
  b(:)=brr(:,l+1)
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
  btemp(:)=brr(:,i) 
  brr(:,i)=brr(:,j) 
  brr(:,j)=btemp(:)
  goto 3
5 arr(l+1)=arr(j) 
  arr(j)=a
  brr(:,l+1)=brr(:,j) 
  brr(:,j)=b(:) 
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
END

SUBROUTINE sort2_int8_long(n,arr,brr) 
INTEGER(8) n,M,NSTACK
INTEGER(8) arr(n),brr(n) 
PARAMETER (M=7,NSTACK=50)
INTEGER(8) i,ir,j,jstack,k,l,istack(NSTACK) 
INTEGER(8) a,b,temp

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
END
