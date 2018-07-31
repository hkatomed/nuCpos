subroutine HBA_3(inseq, logasc, &
                        &freqL1, tranL1, TtranL2, TtranL3, TtranL4, TfreqN4, TtranN4)

  implicit none
  character(147)   inseq
  integer     c(147), w(147), i, j, k, l, n, t, z
  real(8)    asc
  real(8)    logasc
  real(8)    TtranL2(16,4),TtranL3(64,4),TtranL4(256,4)
  real(8)    TtranN4((147-4)*256,4), TfreqN4(64,4)
  real(8)    freqL1(4),tranL1(4,4),tranL2(4,4,4),tranL3(4,4,4,4),tranL4(4,4,4,4,4)
  real(8)    freqN4(4,4,4,4)
  real(8)    tranN4(5:147,4,4,4,4,4),freqL4(4,4,4,4)


  do i=1,147
    if(inseq(i:i)=='A'.or.inseq(i:i)=='a') then
      w(i) = 1
    elseif(inseq(i:i)=='C'.or.inseq(i:i)=='c') then
      w(i) = 2
    elseif(inseq(i:i)=='G'.or.inseq(i:i)=='g') then
      w(i) = 3
    elseif(inseq(i:i)=='T'.or.inseq(i:i)=='t') then
      w(i) = 4
    else
      return
    endif
  enddo

  do i=1,147
    c(i) = 5_1 - w(147 - i + 1)
  end do

  do i=1,4; do j=1,4; tranL2(i,j,:)=TtranL2((i-1)*4+j,:)
  do k=1,4; tranL3(i,j,k,:)=TtranL3((i-1)*16+(j-1)*4+k,:)
  do l=1,4; tranL4(i,j,k,l,:)=TtranL4((i-1)*64+(j-1)*16+(k-1)*4+l,:)
  end do; end do; end do; end do

  do i=1,4; do j=1,4; do k=1,4
    freqN4(i,j,k,:)=TfreqN4((i-1)*16+(j-1)*4+k,:)
  end do; end do; end do

  do n=5,147; do i=1,4; do j=1,4; do k=1,4; do l=1,4
    tranN4(n,i,j,k,l,:)=TtranN4((n-5)*256+(i-1)*64+(j-1)*16+(k-1)*4+l,:)
  end do; end do; end do; end do; end do

  do i=1,4; do j=1,4; do k=1,4; do l=1,4
    freqL4(i,j,k,l)=freqL1(i)*tranL1(i,j)*tranL2(i,j,k)*tranL3(i,j,k,l)
  end do; end do; end do; end do

! 1:147 (asc)
  t = 74; z = 147
  asc=freqN4(w(t-73),w(t-72),w(t-71),w(t-70))&
      &/freqL4(w(t-73),w(t-72),w(t-71),w(t-70))&
      &*freqN4(c(z-t-72),c(z-t-71),c(z-t-70),c(z-t-69))&
      &/freqL4(c(z-t-72),c(z-t-71),c(z-t-70),c(z-t-69))
  do i=5,147
      asc=asc*tranN4(i,w(t-78+i),w(t-77+i),w(t-76+i),w(t-75+i),w(t-74+i))&
      &/tranL4(w(t-78+i),w(t-77+i),w(t-76+i),w(t-75+i),w(t-74+i))&
      &*tranN4(i,c(z-t-77+i),c(z-t-76+i),c(z-t-75+i),c(z-t-74+i),c(z-t-73+i))&
      &/tranL4(c(z-t-77+i),c(z-t-76+i),c(z-t-75+i),c(z-t-74+i),c(z-t-73+i))
  end do

  logasc = log(asc)

end subroutine HBA_3

