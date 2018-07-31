subroutine localHBA_3(inseq, logascSA, logascSB, logascSC, &
			&logascSD, logascSE, logascSF, logascSG, &
			&logascSH, logascSI, logascSJ, &
			&logascSK, logascSL, logascSM, &
                        &freqL1, tranL1, TtranL2, TtranL3, TtranL4, &
                        &TfreqN4SA, TfreqN4SB, TfreqN4SC, TfreqN4SD, TfreqN4SE, &
			&TfreqN4SF, TfreqN4SG, &
                        &TfreqN4SH, TfreqN4SI, TfreqN4SJ, TfreqN4SK, TfreqN4SL, &
			&TfreqN4SM, TtranN4)

  implicit none
  character(147)   inseq
  integer     c(147), w(147), i, j, k, l, n
  real(8)    ascSA, ascSB, ascSC, ascSD, ascSE, ascSF, ascSG
  real(8)    ascSH, ascSI, ascSJ, ascSK, ascSL, ascSM
  real(8)    logascSA, logascSB, logascSC, logascSD
  real(8)    logascSE, logascSF, logascSG
  real(8)    logascSH, logascSI, logascSJ, logascSK, logascSL, logascSM
  real(8)    TtranL2(16,4),TtranL3(64,4),TtranL4(256,4)
  real(8)    TfreqN4SA(64,4),TfreqN4SB(64,4),TfreqN4SC(64,4)
  real(8)    TfreqN4SD(64,4),TfreqN4SE(64,4),TfreqN4SF(64,4),TfreqN4SG(64,4)
  real(8)    TfreqN4SH(64,4),TfreqN4SI(64,4),TfreqN4SJ(64,4)
  real(8)    TfreqN4SK(64,4),TfreqN4SL(64,4),TfreqN4SM(64,4)
  real(8)    TtranN4((147-4)*256,4)
  real(8)    freqL1(4),tranL1(4,4),tranL2(4,4,4),tranL3(4,4,4,4),tranL4(4,4,4,4,4)
  real(8)    freqN4SA(4,4,4,4),freqN4SB(4,4,4,4),freqN4SC(4,4,4,4)
  real(8)    freqN4SD(4,4,4,4),freqN4SE(4,4,4,4),freqN4SF(4,4,4,4),freqN4SG(4,4,4,4)
  real(8)    freqN4SH(4,4,4,4),freqN4SI(4,4,4,4),freqN4SJ(4,4,4,4)
  real(8)    freqN4SK(4,4,4,4),freqN4SL(4,4,4,4),freqN4SM(4,4,4,4)
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
    freqN4SA(i,j,k,:)=TfreqN4SA((i-1)*16+(j-1)*4+k,:)
  end do; end do; end do
  do i=1,4; do j=1,4; do k=1,4
    freqN4SB(i,j,k,:)=TfreqN4SB((i-1)*16+(j-1)*4+k,:)
  end do; end do; end do
  do i=1,4; do j=1,4; do k=1,4
    freqN4SC(i,j,k,:)=TfreqN4SC((i-1)*16+(j-1)*4+k,:)
  end do; end do; end do
  do i=1,4; do j=1,4; do k=1,4
    freqN4SD(i,j,k,:)=TfreqN4SD((i-1)*16+(j-1)*4+k,:)
  end do; end do; end do
  do i=1,4; do j=1,4; do k=1,4
    freqN4SE(i,j,k,:)=TfreqN4SE((i-1)*16+(j-1)*4+k,:)
  end do; end do; end do
  do i=1,4; do j=1,4; do k=1,4
    freqN4SF(i,j,k,:)=TfreqN4SF((i-1)*16+(j-1)*4+k,:)
  end do; end do; end do
  do i=1,4; do j=1,4; do k=1,4
    freqN4SG(i,j,k,:)=TfreqN4SG((i-1)*16+(j-1)*4+k,:)
  end do; end do; end do

  do i=1,4; do j=1,4; do k=1,4
    freqN4SH(i,j,k,:)=TfreqN4SH((i-1)*16+(j-1)*4+k,:)
  end do; end do; end do
  do i=1,4; do j=1,4; do k=1,4
    freqN4SI(i,j,k,:)=TfreqN4SI((i-1)*16+(j-1)*4+k,:)
  end do; end do; end do
  do i=1,4; do j=1,4; do k=1,4
    freqN4SJ(i,j,k,:)=TfreqN4SJ((i-1)*16+(j-1)*4+k,:)
  end do; end do; end do
  do i=1,4; do j=1,4; do k=1,4
    freqN4SK(i,j,k,:)=TfreqN4SK((i-1)*16+(j-1)*4+k,:)
  end do; end do; end do
  do i=1,4; do j=1,4; do k=1,4
    freqN4SL(i,j,k,:)=TfreqN4SL((i-1)*16+(j-1)*4+k,:)
  end do; end do; end do
  do i=1,4; do j=1,4; do k=1,4
    freqN4SM(i,j,k,:)=TfreqN4SM((i-1)*16+(j-1)*4+k,:)
  end do; end do; end do

  do n=5,147; do i=1,4; do j=1,4; do k=1,4; do l=1,4
    tranN4(n,i,j,k,l,:)=TtranN4((n-5)*256+(i-1)*64+(j-1)*16+(k-1)*4+l,:)
  end do; end do; end do; end do; end do

  do i=1,4; do j=1,4; do k=1,4; do l=1,4
    freqL4(i,j,k,l)=freqL1(i)*tranL1(i,j)*tranL2(i,j,k)*tranL3(i,j,k,l)
  end do; end do; end do; end do

! 1:21 (ascSA)
  ascSA=freqN4SA(w(1),w(2),w(3),w(4))&
      &/freqL4(w(1),w(2),w(3),w(4))&
      &*freqN4SG(c(127),c(128),c(129),c(130))&
      &/freqL4(c(127),c(128),c(129),c(130))
  do i=5,21
    ascSA=ascSA*tranN4(i,w(i-4),w(i-3),w(i-2),w(i-1),w(i))&
      &/tranL4(w(i-4),w(i-3),w(i-2),w(i-1),w(i))
  end do
  do i=131,147
    ascSA=ascSA*tranN4(i,c(i-4),c(i-3),c(i-2),c(i-1),c(i))&
      &/tranL4(c(i-4),c(i-3),c(i-2),c(i-1),c(i))
  end do

! 22:42 (ascSB)
  ascSB=freqN4SB(w(22),w(23),w(24),w(25))&
      &/freqL4(w(22),w(23),w(24),w(25))&
      &*freqN4SF(c(106),c(107),c(108),c(109))&
      &/freqL4(c(106),c(107),c(108),c(109))
  do i=26,42
    ascSB=ascSB*tranN4(i,w(i-4),w(i-3),w(i-2),w(i-1),w(i))&
      &/tranL4(w(i-4),w(i-3),w(i-2),w(i-1),w(i))
  end do
  do i=110,126
    ascSB=ascSB*tranN4(i,c(i-4),c(i-3),c(i-2),c(i-1),c(i))&
      &/tranL4(c(i-4),c(i-3),c(i-2),c(i-1),c(i))
  end do

! 43:63 (ascSC)
  ascSC=freqN4SC(w(43),w(44),w(45),w(46))&
      &/freqL4(w(43),w(44),w(45),w(46))&
      &*freqN4SE(c(85),c(86),c(87),c(88))&
      &/freqL4(c(85),c(86),c(87),c(88))
  do i=47,63
    ascSC=ascSC*tranN4(i,w(i-4),w(i-3),w(i-2),w(i-1),w(i))&
      &/tranL4(w(i-4),w(i-3),w(i-2),w(i-1),w(i))
  end do
  do i=89,105
    ascSC=ascSC*tranN4(i,c(i-4),c(i-3),c(i-2),c(i-1),c(i))&
      &/tranL4(c(i-4),c(i-3),c(i-2),c(i-1),c(i))
  end do

! 64:84 (ascSD)
  ascSD=freqN4SD(w(64),w(65),w(66),w(67))&
      &/freqL4(w(64),w(65),w(66),w(67))&
      &*freqN4SD(c(64),c(65),c(66),c(67))&
      &/freqL4(c(64),c(65),c(66),c(67))
  do i=68,84
    ascSD=ascSD*tranN4(i,w(i-4),w(i-3),w(i-2),w(i-1),w(i))&
      &/tranL4(w(i-4),w(i-3),w(i-2),w(i-1),w(i))&
      &*tranN4(i,c(i-4),c(i-3),c(i-2),c(i-1),c(i))&
      &/tranL4(c(i-4),c(i-3),c(i-2),c(i-1),c(i))
  end do

! 85:105 (ascSE)
  ascSE=freqN4SE(w(85),w(86),w(87),w(88))&
      &/freqL4(w(85),w(86),w(87),w(88))&
      &*freqN4SC(c(43),c(44),c(45),c(46))&
      &/freqL4(c(43),c(44),c(45),c(46))
  do i=89,105
    ascSE=ascSE*tranN4(i,w(i-4),w(i-3),w(i-2),w(i-1),w(i))&
      &/tranL4(w(i-4),w(i-3),w(i-2),w(i-1),w(i))
  end do
  do i=47,63
    ascSE=ascSE*tranN4(i,c(i-4),c(i-3),c(i-2),c(i-1),c(i))&
      &/tranL4(c(i-4),c(i-3),c(i-2),c(i-1),c(i))
  end do

! 106:126 (ascSF)
  ascSF=freqN4SF(w(106),w(107),w(108),w(109))&
      &/freqL4(w(106),w(107),w(108),w(109))&
      &*freqN4SB(c(22),c(23),c(24),c(25))&
      &/freqL4(c(22),c(23),c(24),c(25))
  do i=110,126
    ascSF=ascSF*tranN4(i,w(i-4),w(i-3),w(i-2),w(i-1),w(i))&
      &/tranL4(w(i-4),w(i-3),w(i-2),w(i-1),w(i))
  end do
  do i=26,42
    ascSF=ascSF*tranN4(i,c(i-4),c(i-3),c(i-2),c(i-1),c(i))&
      &/tranL4(c(i-4),c(i-3),c(i-2),c(i-1),c(i))
  end do

! 127:147 (ascSG)
  ascSG=freqN4SG(w(127),w(128),w(129),w(130))&
      &/freqL4(w(127),w(128),w(129),w(130))&
      &*freqN4SA(c(1),c(2),c(3),c(4))&
      &/freqL4(c(1),c(2),c(3),c(4))
  do i=131,147
    ascSG=ascSG*tranN4(i,w(i-4),w(i-3),w(i-2),w(i-1),w(i))&
      &/tranL4(w(i-4),w(i-3),w(i-2),w(i-1),w(i))
  end do
  do i=5,21
    ascSG=ascSG*tranN4(i,c(i-4),c(i-3),c(i-2),c(i-1),c(i))&
      &/tranL4(c(i-4),c(i-3),c(i-2),c(i-1),c(i))
  end do


! 12:31 (ascSH)
  ascSH=freqN4SH(w(12),w(13),w(14),w(15))&
      &/freqL4(w(12),w(13),w(14),w(15))&
      &*freqN4SM(c(117),c(118),c(119),c(120))&
      &/freqL4(c(117),c(118),c(119),c(120))
  do i=16,31
    ascSH=ascSH*tranN4(i,w(i-4),w(i-3),w(i-2),w(i-1),w(i))&
      &/tranL4(w(i-4),w(i-3),w(i-2),w(i-1),w(i))
  end do
  do i=121,136
    ascSH=ascSH*tranN4(i,c(i-4),c(i-3),c(i-2),c(i-1),c(i))&
      &/tranL4(c(i-4),c(i-3),c(i-2),c(i-1),c(i))
  end do

! 33:52 (ascSI)
  ascSI=freqN4SI(w(33),w(34),w(35),w(36))&
      &/freqL4(w(33),w(34),w(35),w(36))&
      &*freqN4SL(c(96),c(97),c(98),c(99))&
      &/freqL4(c(96),c(97),c(98),c(99))
  do i=37,52
    ascSI=ascSI*tranN4(i,w(i-4),w(i-3),w(i-2),w(i-1),w(i))&
      &/tranL4(w(i-4),w(i-3),w(i-2),w(i-1),w(i))
  end do
  do i=100,115
    ascSI=ascSI*tranN4(i,c(i-4),c(i-3),c(i-2),c(i-1),c(i))&
      &/tranL4(c(i-4),c(i-3),c(i-2),c(i-1),c(i))
  end do

! 54:73 (ascSJ)
  ascSJ=freqN4SJ(w(54),w(55),w(56),w(57))&
      &/freqL4(w(54),w(55),w(56),w(57))&
      &*freqN4SK(c(75),c(76),c(77),c(78))&
      &/freqL4(c(75),c(76),c(77),c(78))
  do i=58,73
    ascSJ=ascSJ*tranN4(i,w(i-4),w(i-3),w(i-2),w(i-1),w(i))&
      &/tranL4(w(i-4),w(i-3),w(i-2),w(i-1),w(i))
  end do
  do i=79,94
    ascSJ=ascSJ*tranN4(i,c(i-4),c(i-3),c(i-2),c(i-1),c(i))&
      &/tranL4(c(i-4),c(i-3),c(i-2),c(i-1),c(i))
  end do

! 75:94 (ascSK)
  ascSK=freqN4SK(w(75),w(76),w(77),w(78))&
      &/freqL4(w(75),w(76),w(77),w(78))&
      &*freqN4SJ(c(54),c(55),c(56),c(57))&
      &/freqL4(c(54),c(55),c(56),c(57))
  do i=79,94
    ascSK=ascSK*tranN4(i,w(i-4),w(i-3),w(i-2),w(i-1),w(i))&
      &/tranL4(w(i-4),w(i-3),w(i-2),w(i-1),w(i))
  end do
  do i=58,73
    ascSK=ascSK*tranN4(i,c(i-4),c(i-3),c(i-2),c(i-1),c(i))&
      &/tranL4(c(i-4),c(i-3),c(i-2),c(i-1),c(i))
  end do

! 96:115 (ascSL)
  ascSL=freqN4SL(w(96),w(97),w(98),w(99))&
      &/freqL4(w(96),w(97),w(98),w(99))&
      &*freqN4SI(c(33),c(34),c(35),c(36))&
      &/freqL4(c(33),c(34),c(35),c(36))
  do i=100,115
    ascSL=ascSL*tranN4(i,w(i-4),w(i-3),w(i-2),w(i-1),w(i))&
      &/tranL4(w(i-4),w(i-3),w(i-2),w(i-1),w(i))
  end do
  do i=37,52
    ascSL=ascSL*tranN4(i,c(i-4),c(i-3),c(i-2),c(i-1),c(i))&
      &/tranL4(c(i-4),c(i-3),c(i-2),c(i-1),c(i))
  end do

! 117:136 (ascSM)
  ascSM=freqN4SM(w(117),w(118),w(119),w(120))&
      &/freqL4(w(117),w(118),w(119),w(120))&
      &*freqN4SH(c(12),c(13),c(14),c(15))&
      &/freqL4(c(12),c(13),c(14),c(15))
  do i=121,136
    ascSM=ascSM*tranN4(i,w(i-4),w(i-3),w(i-2),w(i-1),w(i))&
      &/tranL4(w(i-4),w(i-3),w(i-2),w(i-1),w(i))
  end do
  do i=16,31
    ascSM=ascSM*tranN4(i,c(i-4),c(i-3),c(i-2),c(i-1),c(i))&
      &/tranL4(c(i-4),c(i-3),c(i-2),c(i-1),c(i))
  end do

  logascSA = log(ascSA)
  logascSB = log(ascSB)
  logascSC = log(ascSC)
  logascSD = log(ascSD)
  logascSE = log(ascSE)
  logascSF = log(ascSF)
  logascSG = log(ascSG)
  logascSH = log(ascSH)
  logascSI = log(ascSI)
  logascSJ = log(ascSJ)
  logascSK = log(ascSK)
  logascSL = log(ascSL)
  logascSM = log(ascSM)

end subroutine localHBA_3

