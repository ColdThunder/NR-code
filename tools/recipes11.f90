SUBROUTINE elmhes(a,n,np)
INTEGER n,np
REAL a(np,np)
INTEGER i,j,m
REAL x,y
do m=2,n-1
  x=0.
  i=m
  do j=m,n
    if(abs(a(j,m-1)).gt.abs(x))then
      x=a(j,m-1)
      i=j
    endif
  enddo
  if(i.ne.m)then
    do j=m-1,n
      y=a(i,j)
      a(i,j)=a(m,j)
      a(m,j)=y
    enddo
    do j=1,n
      y=a(j,i)
      a(j,i)=a(j,m)
      a(j,m)=y
    enddo
  endif
  if(x.ne.0.)then
    do i=m+1,n
      y=a(i,m-1)
      if(y.ne.0.)then
        y=y/x
        a(i,m-1)=y
        do j=m,n
          a(i,j)=a(i,j)-y*a(m,j)
        enddo
        do j=1,n
          a(j,m)=a(j,m)+y*a(j,i)
        enddo
      endif
    enddo
  endif
enddo
return
ENDSUBROUTINE elmhes

SUBROUTINE hqr(a,n,np,wr,wi)
INTEGER n,np
REAL a(np,np),wi(np),wr(np)
INTEGER i,its,j,k,l,m,nn
REAL anorm,p,q,r,s,t,u,v,w,x,y,z
anorm=abs(a(1,1))
do i=2,n
  do j=i-1,n
    anorm=anorm+abs(a(i,j))
  enddo
enddo
nn=n
t=0.
1 if(nn.ge.1)then
 its=0
2  do l=nn,2,-1
    s=abs(a(l-1,l-1))+abs(a(l,l))
    if(s.eq.0.)s=anorm
    if(abs(a(l,l-1))+s.eq.s)goto 3
  enddo
  l=1
3 x=a(nn,nn)
  if(l.eq.nn)then
    wr(nn)=x+t
    wi(nn)=0.
    nn=nn-1
  else
    y=a(nn-1,nn-1)
    w=a(nn,nn-1)*a(nn-1,nn)
    if(l.eq.nn-1)then
      p=0.5*(y-x)
      q=p**2+w
      z=sqrt(abs(q))
      x=x+t
      if(q.ge.0.)then
        z=p+sign(z,p)
        wr(nn)=x+z
        wr(nn-1)=wr(nn)
        if(z.ne.0.)wr(nn)=x-w/z
        wi(nn)=0.
        wi(nn-1)=0.
      else
        wr(nn)=x+p
        wr(nn-1)=wr(nn)
        wi(nn)=z
        wi(nn-1)=-z
      endif
      nn=nn-2
    else
      if(its.eq.30)pause 'too many iterations in hqr'
      if(its.eq.10.or.its.eq.20)then
        t=t+x
        do i=1,nn
          a(i,i)=a(i,i)-x
        enddo
        s=abs(a(nn,nn-1))+abs(a(nn-1,nn-2))
        x=0.75*s
        y=x
        w=-0.4375*s**2
      endif
      its=its+1
      do m=nn-2,l,-1
        z=a(m,m)
        r=x-z
        s=y-z
        p=(r*s-w)/a(m+1,m)+a(m,m+1)
        q=a(m+1,m+1)-z-r-s
        r=a(m+2,m+1)
        s=abs(p)+abs(q)+abs(r)
        p=p/s
        q=q/s
        r=r/s
        if(m.eq.l)goto 4
        u=abs(a(m,m-1))*(abs(q)+abs(r))
        v=abs(p)*(abs(a(m-1,m-1))+abs(z)+abs(a(m+1,m+1)))
        if(u+v.eq.v)goto 4
      enddo
4     do i=m+2,nn
        a(i,i-2)=0.
        if (i.ne.m+2) a(i,i-3)=0.
      enddo
      do k=m,nn-1
        if(k.ne.m)then
          p=a(k,k-1)
          q=a(k+1,k-1)
          r=0.
          if(k.ne.nn-1)r=a(k+2,k-1)
          x=abs(p)+abs(q)+abs(r)
          if(x.ne.0.)then
            p=p/x
            q=q/x
            r=r/x
          endif
        endif
        s=sign(sqrt(p**2+q**2+r**2),p)
        if(s.ne.0.)then
          if(k.eq.m)then
            if(l.ne.m)a(k,k-1)=-a(k,k-1)
          else
            a(k,k-1)=-s*x
          endif
          p=p+s
          x=p/s
          y=q/s
          z=r/s
          q=q/p
          r=r/p
          do j=k,nn
            p=a(k,j)+q*a(k+1,j)
            if(k.ne.nn-1)then
              p=p+r*a(k+2,j)
              a(k+2,j)=a(k+2,j)-p*z
            endif
            a(k+1,j)=a(k+1,j)-p*y
            a(k,j)=a(k,j)-p*x
          enddo
          do i=l,min(nn,k+3)
            p=x*a(i,k)+y*a(i,k+1)
            if(k.ne.nn-1)then
              p=p+z*a(i,k+2)
              a(i,k+2)=a(i,k+2)-p*r
            endif
            a(i,k+1)=a(i,k+1)-p*q
            a(i,k)=a(i,k)-p
          enddo
        endif
      enddo
      goto 2
    endif
  endif
goto 1
endif
return
ENDSUBROUTINE hqr

SUBROUTINE jacobi(a,n,np,d,v,nrot)
INTEGER n,np,nrot,NMAX
REAL a(np,np),d(np),v(np,np)
PARAMETER (NMAX=500)
INTEGER i,ip,iq,j
REAL c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
do ip=1,n
  do iq=1,n
    v(ip,iq)=0.
  enddo
  v(ip,ip)=1.
enddo
do ip=1,n
  b(ip)=a(ip,ip)
  d(ip)=b(ip)
  z(ip)=0.
enddo
nrot=0
do i=1,50
  sm=0.
  do ip=1,n-1
    do iq=ip+1,n
      sm=sm+abs(a(ip,iq))
    enddo
  enddo
  if(sm.eq.0.)return
  if(i.lt.4)then
    tresh=0.2*sm/n**2
  else
    tresh=0.
  endif
  do ip=1,n-1
    do iq=ip+1,n
      g=100.*abs(a(ip,iq))
      if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
        a(ip,iq)=0.
      else if(abs(a(ip,iq)).gt.tresh)then
        h=d(iq)-d(ip)
        if(abs(h)+g.eq.abs(h))then
          t=a(ip,iq)/h
        else
          theta=0.5*h/a(ip,iq)
          t=1./(abs(theta)+sqrt(1.+theta**2))
          if(theta.lt.0.)t=-t
        endif
        c=1./sqrt(1+t**2)
        s=t*c
        tau=s/(1.+c)
        h=t*a(ip,iq)
        z(ip)=z(ip)-h
        z(iq)=z(iq)+h
        d(ip)=d(ip)-h
        d(iq)=d(iq)+h
        a(ip,iq)=0.
        do j=1,ip-1
          g=a(j,ip)
          h=a(j,iq)
          a(j,ip)=g-s*(h+g*tau)
          a(j,iq)=h+s*(g-h*tau)
        enddo
        do j=ip+1,iq-1
          g=a(ip,j)
          h=a(j,iq)
          a(ip,j)=g-s*(h+g*tau)
          a(j,iq)=h+s*(g-h*tau)
        enddo
        do j=iq+1,n
          g=a(ip,j)
          h=a(iq,j)
          a(ip,j)=g-s*(h+g*tau)
          a(iq,j)=h+s*(g-h*tau)
        enddo
        do j=1,n
          g=v(j,ip)
          h=v(j,iq)
          v(j,ip)=g-s*(h+g*tau)
          v(j,iq)=h+s*(g-h*tau)
        enddo
        nrot=nrot+1
      endif
    enddo
  enddo
  do ip=1,n
    b(ip)=b(ip)+z(ip)
    d(ip)=b(ip)
    z(ip)=0.
  enddo
enddo
pause 'too many iterations in jacobi'
return
ENDSUBROUTINE jacobi

