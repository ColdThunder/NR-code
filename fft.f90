subroutine fft3(a,L,command)
implicit none
integer(4)::L(3)
real(4)::a(L(1)+2,L(2),L(3))
integer(4)::command
call rlft3(a(1:L(1),:,:),a(L(1)+1:L(1)+2,:,:),L,command)
if (command.eq.-1) then
  a(1:L(1),:,:)=a(1:L(1),:,:)*2./L(1)/L(2)/L(3)
endif
endsubroutine fft3
