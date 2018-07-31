
subroutine nuCpos2_2(inseq,seqlen,freqL1,tranL1,TtranL2,TtranL3,TtranL4,TfreqN4,TtranN4,mlL,Pd&
  &,std,pstart,nucoccup,viterbi,affinity)

  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer   t,d,i,j,k,l,m,n,z,mlL,nvt,std,seqlen
  real*8    tempN,tempL,temp,tp,Pd(mlL)
  real*8    TtranL2(16,4),TtranL3(64,4),TtranL4(256,4),TfreqN4(64,4),TtranN4((147-4)*256,4)
  real*8    freqL1(4),tranL1(4,4),tranL2(4,4,4),tranL3(4,4,4,4),tranL4(4,4,4,4,4)
  real*8    freqN4(4,4,4,4),tranN4(5:147,4,4,4,4,4),freqL4(4,4,4,4)

  integer*1,allocatable:: w(:),c(:)
  integer,allocatable:: change(:),Nst(:); real*8,allocatable:: ppEndN(:),ppN(:),asc(:)
  real*8,allocatable:: hatFN(:),hatFL(:),r(:),hatAN(:),hatAL(:),ra(:),hatBN(:),hatBL(:),rb(:)
  character(len=seqlen) inseq
  ! real(8) pstart(seqlen),nucoccup(seqlen),viterbi(seqlen),affinity(seqlen)
  double precision,dimension(seqlen):: pstart,nucoccup,viterbi,affinity
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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


  z=seqlen
  allocate(w(z),c(z),change(z))
  allocate(r(z),hatFN(0:z),hatFL(z))
  allocate(ra(z),hatAN(0:z),hatAL(z),rb(0:(z-1)),hatBN(0:(z-1)),hatBL(0:z))

  ! print *, inseq(1:10)

  do i=1,z
    if(inseq(i:i)=='A'.or.inseq(i:i)=='a') then
      w(i) = 1
    elseif(inseq(i:i)=='C'.or.inseq(i:i)=='c') then
      w(i) = 2
    elseif(inseq(i:i)=='G'.or.inseq(i:i)=='g') then
      w(i) = 3
    elseif(inseq(i:i)=='T'.or.inseq(i:i)=='t') then
      w(i) = 4
    elseif(inseq(i:i)=='N'.or.inseq(i:i)=='n') then
      w(i) = 0
    else
      return
    endif
  enddo

  do i=1,z
    c(i) = 5_1 - w(z - i + 1)
  end do

  ! print *, w(1:10)
  ! print *, c(1:10)

  !!!!!!!!!!!!!!!!!!!!!!!  viterbi algorithm  !!!!!!!!!!!!!!!!!!!!!!!!

  change=0; hatFN(0)=1.0; hatAN(0)=1.0; hatBL(z)=1.0
  !!!!!!!!!!!!!! 1 -- 147 !!!!!!!!!!!!!!

  hatFN(1)=0.0; hatFL(1)=1.0
  r(1)=Pd(1)*freqL1(w(1))*freqL1(c(z))

  hatFN(2)=0.0; hatFL(2)=1.0
  r(2)=Pd(2)*freqL1(w(1))*freqL1(c(z-1))*tranL1(w(1),w(2))*tranL1(c(z-1),c(z))/r(1)

  hatFN(3)=0.0; hatFL(3)=1.0
  r(3)=Pd(3)*freqL1(w(1))*freqL1(c(z-2))*tranL1(w(1),w(2))*tranL1(c(z-2),c(z-1))/r(1)&
  &*tranL2(w(1),w(2),w(3))*tranL2(c(z-2),c(z-1),c(z))/r(2)

  hatFN(4)=0.0; hatFL(4)=1.0
  r(4)=Pd(4)*freqL4(w(1),w(2),w(3),w(4))*freqL4(c(z-3),c(z-2),c(z-1),c(z))/r(1)/r(2)/r(3)

  do t=5,min(mlL,147)
    hatFN(t)=0.0; hatFL(t)=1.0

    r(t)=Pd(t)*freqL4(w(1),w(2),w(3),w(4))*freqL4(c(z-t+1),c(z-t+2),c(z-t+3),c(z-t+4))/r(1)/r(2)/r(3)
    do i=5,t
      r(t)=r(t)*tranL4(w(i-4),w(i-3),w(i-2),w(i-1),w(i))*tranL4(c(z-t+i-4),c(z-t+i-3),c(z-t+i-2),c(z-t+i-1),c(z-t+i))/r(i-1)
    end do
  end do

  do t=mlL+1,147; hatFN(t)=0.0; hatFL(t)=0.0; r(t)=1.0/16.0; end do



  hatAN(1:147)=hatFN(1:147); hatAL(1:147)=hatFL(1:147); ra(1:147)=r(1:147)
  !!!!!!!!!!!!!!!! 148 -- z !!!!!!!!!!!!!
      
  do t=148,z

    tempN=0.0
    if(hatFL(t-147)/=0.0.and.t<z) then
      tempN=hatFL(t-147)*freqN4(w(t-146),w(t-145),w(t-144),w(t-143))*freqN4(c(z-t+1),c(z-t+2),c(z-t+3),c(z-t+4))&
      &/r(t-146)/r(t-145)/r(t-144)
      do i=5,147
        tempN=tempN*tranN4(i,w(t-151+i),w(t-150+i),w(t-149+i),w(t-148+i),w(t-147+i))/r(t-148+i)&
        &*tranN4(i,c(z-t+i-4),c(z-t+i-3),c(z-t+i-2),c(z-t+i-1),c(z-t+i))
      end do
    end if

    tempL=0.0
    do d=1,min(mlL,t)
      temp=0.0

      if(d<=3.and.hatFN(t-d)/=0.0) then
        temp=hatFN(t-d)*Pd(d)*freqL1(w(t-d+1))*freqL1(c(z-t+1))
        if(d>=2) temp=temp*tranL1(w(t-d+1),w(t-d+2))*tranL1(c(z-t+1),c(z-t+2))/r(t-d+1)
        if(d>=3) temp=temp*tranL2(w(t-d+1),w(t-d+2),w(t-d+3))*tranL2(c(z-t+1),c(z-t+2),c(z-t+3))/r(t-d+2)  
      end if

      if(d>=4) then
        if(d==4) then
          tp=freqL4(c(z-t+1),c(z-t+2),c(z-t+3),c(z-t+4))
        else if(d>4) then
          tp=tranL4(w(t-d+1),w(t-d+2),w(t-d+3),w(t-d+4),w(t-d+5))*tranL4(c(z-t+d-4),c(z-t+d-3),c(z-t+d-2),c(z-t+d-1),c(z-t+d))&
          &*tp/r(t-d+4)
        end if
        if(hatFN(t-d)/=0.0) temp=hatFN(t-d)*Pd(d)*freqL4(w(t-d+1),w(t-d+2),w(t-d+3),w(t-d+4))*tp/r(t-d+1)/r(t-d+2)/r(t-d+3)
      end if

      if(temp>tempL) then
        tempL=temp; change(t)=t-d
      end if
    end do

    r(t)=tempN+tempL
    if(r(t)==0.0) then
      r(t)=1.0/16.0; hatFN(t)=0.0; hatFL(t)=0.0
    else
      hatFN(t)=tempN/r(t); hatFL(t)=tempL/r(t)
    end if


    !FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    tempN=0.0
    if(hatAL(t-147)/=0.0.and.t<z) then
      tempN=hatAL(t-147)*freqN4(w(t-146),w(t-145),w(t-144),w(t-143))*freqN4(c(z-t+1),c(z-t+2),c(z-t+3),c(z-t+4))&
      &/ra(t-146)/ra(t-145)/ra(t-144)
      do i=5,147
        tempN=tempN*tranN4(i,w(t-151+i),w(t-150+i),w(t-149+i),w(t-148+i),w(t-147+i))/ra(t-148+i)&
        &*tranN4(i,c(z-t+i-4),c(z-t+i-3),c(z-t+i-2),c(z-t+i-1),c(z-t+i))
      end do
    end if

    tempL=0.0
    do d=1,min(mlL,t)
      temp=0.0

      if(d<=3.and.hatAN(t-d)/=0.0) then
        temp=hatAN(t-d)*Pd(d)*freqL1(w(t-d+1))*freqL1(c(z-t+1))
        if(d>=2) temp=temp*tranL1(w(t-d+1),w(t-d+2))*tranL1(c(z-t+1),c(z-t+2))/ra(t-d+1)
        if(d>=3) temp=temp*tranL2(w(t-d+1),w(t-d+2),w(t-d+3))*tranL2(c(z-t+1),c(z-t+2),c(z-t+3))/ra(t-d+2)  
      end if

      if(d>=4) then
        if(d==4) then
          tp=freqL4(c(z-t+1),c(z-t+2),c(z-t+3),c(z-t+4))
        else if(d>4) then
          tp=tranL4(w(t-d+1),w(t-d+2),w(t-d+3),w(t-d+4),w(t-d+5))*tranL4(c(z-t+d-4),c(z-t+d-3),c(z-t+d-2),c(z-t+d-1),c(z-t+d))*tp&
          &/ra(t-d+4)
        end if
        if(hatAN(t-d)/=0.0) temp=hatAN(t-d)*Pd(d)*freqL4(w(t-d+1),w(t-d+2),w(t-d+3),w(t-d+4))*tp/ra(t-d+1)/ra(t-d+2)/ra(t-d+3)
      end if

      tempL=tempL+temp
    end do

    ra(t)=tempN+tempL
    if(ra(t)==0.0) then
      ra(t)=1.0/16.0; hatAN(t)=0.0; hatAL(t)=0.0
    else
      hatAN(t)=tempN/ra(t); hatAL(t)=tempL/ra(t)
    end if
    !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

  end do
  !!!!!!!!!!!!!!!!!!!!!!! backward !!!!!!!!!!!!!!!!!!!!!!!!!

  !FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
  do t=z-1,max(z-mlL,z-147),-1
    hatBN(t)=1.0; hatBL(t)=0.0
    rb(t)=Pd(z-t)*freqL1(w(t+1))*freqL1(c(1))
    if(z-t>=2) rb(t)=rb(t)*tranL1(w(t+1),w(t+2))*tranL1(c(1),c(2))/rb(t+1)
    if(z-t>=3) rb(t)=rb(t)*tranL2(w(t+1),w(t+2),w(t+3))*tranL2(c(1),c(2),c(3))/rb(t+2)
    if(z-t>=4) rb(t)=rb(t)*tranL3(w(t+1),w(t+2),w(t+3),w(t+4))*tranL3(c(1),c(2),c(3),c(4))/rb(t+3)
    do i=5,z-t
      rb(t)=rb(t)*tranL4(w(t+i-4),w(t+i-3),w(t+i-2),w(t+i-1),w(t+i))&
      &*tranL4(c(i-4),c(i-3),c(i-2),c(i-1),c(i))/rb(t+i-1)
    end do
  end do

  do t=z-mlL-1,z-147,-1; hatBN(t)=0.0; hatBL(t)=0.0; rb(t)=1.0/16.0; end do

  do t=z-148,0,-1

    tempL=0.0
    if(hatBN(t+147)>0.0.and.t>0) then
      tempL=hatBN(t+147)*freqN4(w(t+1),w(t+2),w(t+3),w(t+4))*freqN4(c(z-t-146),c(z-t-145),c(z-t-144),c(z-t-143))&
      &/rb(t+1)/rb(t+2)/rb(t+3)
      do i=5,147
        tempL=tempL*tranN4(i,w(t+i-4),w(t+i-3),w(t+i-2),w(t+i-1),w(t+i))/rb(t+i-1)&
        &*tranN4(i,c(z-t-151+i),c(z-t-150+i),c(z-t-149+i),c(z-t-148+i),c(z-t-147+i))
      end do
    end if

    tempN=0.0
    do d=1,min(mlL,z-t)
      temp=0.0

      if(d<=3.and.hatBL(t+d)>0.0) then
        temp=Pd(d)*freqL1(w(t+1))*freqL1(c(z-t-d+1))*hatBL(t+d)
        if(d>=2) temp=temp*tranL1(w(t+1),w(t+2))*tranL1(c(z-t-d+1),c(z-t-d+2))/rb(t+1)
        if(d>=3) temp=temp*tranL2(w(t+1),w(t+2),w(t+3))*tranL2(c(z-t-d+1),c(z-t-d+2),c(z-t-d+3))/rb(t+2)  
      end if

      if(d>=4) then
        if(d==4) then
          tp=freqL4(w(t+1),w(t+2),w(t+3),w(t+4))/rb(t+1)/rb(t+2)/rb(t+3)
        else if(d>4) then
          tp=tranL4(w(t+d-4),w(t+d-3),w(t+d-2),w(t+d-1),w(t+d))*tranL4(c(z-t-d+1),c(z-t-d+2),c(z-t-d+3),c(z-t-d+4),c(z-t-d+5))&
          &*tp/rb(t+d-1)
        end if
        if(hatBL(t+d)>0.0) temp=tp*freqL4(c(z-t-d+1),c(z-t-d+2),c(z-t-d+3),c(z-t-d+4))*Pd(d)*hatBL(t+d)
      end if

      tempN=tempN+temp
    end do

    rb(t)=tempN+tempL
    if(rb(t)==0.0) then
      rb(t)=1.0/16.0; hatBN(t)=0.0; hatBL(t)=0.0
    else
      hatBN(t)=tempN/rb(t); hatBL(t)=tempL/rb(t)
    end if

  end do
  !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
  !^^^^^^^^^^^^^^^^^^^ end of backward ^^^^^^^^^^^^^^^^^^^

  deallocate(r,hatFN,hatFL,hatAL,hatBL); allocate(ppEndN(z+146),ppN(z),asc(74:z-73),Nst(0:z/147))

    ppEndN=0.0; ppN=0.0; temp=1.0
    do t=z-1,1,-1; temp=temp*rb(t)/ra(t+1); ppEndN(t)=hatAN(t)*hatBN(t)*temp; end do
    do t=1,z; do i=t,t+146; ppN(t)=ppN(t)+ppEndN(i); end do; end do

    do t=74,z-73
      asc(t)=freqN4(w(t-73),w(t-72),w(t-71),w(t-70))/freqL4(w(t-73),w(t-72),w(t-71),w(t-70))&
      &*freqN4(c(z-t-72),c(z-t-71),c(z-t-70),c(z-t-69))/freqL4(c(z-t-72),c(z-t-71),c(z-t-70),c(z-t-69))
      do i=5,147
        asc(t)=asc(t)*tranN4(i,w(t-78+i),w(t-77+i),w(t-76+i),w(t-75+i),w(t-74+i))&
                    &/tranL4(w(t-78+i),w(t-77+i),w(t-76+i),w(t-75+i),w(t-74+i))&
                    &*tranN4(i,c(z-t-77+i),c(z-t-76+i),c(z-t-75+i),c(z-t-74+i),c(z-t-73+i))&
                    &/tranL4(c(z-t-77+i),c(z-t-76+i),c(z-t-75+i),c(z-t-74+i),c(z-t-73+i))
      end do
      asc(t)=log(asc(t))
    end do

    if(std ==1) then
    temp=sum(asc(:))/real(z-146); tp=sqrt(sum((asc(:)-temp)**2)/(z-147.0))
    asc=(asc-temp)/tp
    end if

    n=z; m=change(z); k=0; w=0
    do while(m>0)
      k=k+1; Nst(k)=m-146; w(m-146:m)=1
      n=m-147; m=change(n)
    end do; nvt=k; Nst(0)=z+1

  ! print *, pstart(1:3)
  ! print *, pstart(1000:1003)

  pstart(1:z)=ppEndN(147:z+146)
  nucoccup(1:z)=ppN(1:z)
  viterbi(1:z)=w(1:z)
  affinity(1:73)=0
  affinity(z-72:z)=0
  affinity(74:z-73)=asc(74:z-73)


  ! open(1,file='testXYZ.txt',status='replace')
  ! do t=1,z
  ! if(t<74.or.z-73<t) then
  ! write(1,'(i9,2x,f5.3,3x,f5.3,3x,i1,4x,a)') t,ppEndN(t+146),ppN(t),w(t),'NA'
  ! else
  ! write(1,'(i9,2x,f5.3,3x,f5.3,3x,i1,1x,f8.3)') t,ppEndN(t+146),ppN(t),w(t),asc(t)
  ! end if
  ! end do
  ! close(1)



  deallocate(w,c,change,ra,hatAN,rb,hatBN,ppEndN,ppN,asc,Nst)

  !!!!!!!!
  !!!!!!!!

end subroutine nuCpos2_2

