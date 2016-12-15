subroutine defp2delta(def,delta,nn)
implicit none
integer(4)::nn(3)
real(4)::def(nn(1),nn(2),nn(3))
real(4)::delta(nn(1),nn(2),nn(3))
real(4)::phixx,phiyy,phizz
integer(4)::i,j,k,ip,im,jp,jm,kp,km

do k=1,nn(3)
  kp=mod(k,nn(3))+1
  km=mod(k+nn(3)-2,nn(3))+1
  do j=1,nn(2)
    jp=mod(j,nn(2))+1
    jm=mod(j+nn(2)-2,nn(2))+1
    do i=1,nn(1)
      ip=mod(i,nn(1))+1
      im=mod(i+nn(1)-2,nn(1))+1
      phixx=(def(ip,j,k)-2*def(i,j,k)+def(im,j,k))
      phiyy=(def(i,jp,k)-2*def(i,j,k)+def(i,jm,k))
      phizz=(def(i,j,kp)-2*def(i,j,k)+def(i,j,km))
      delta(i,j,k)=-(phixx+phiyy+phizz)
    enddo
  enddo
enddo

return
endsubroutine defp2delta
