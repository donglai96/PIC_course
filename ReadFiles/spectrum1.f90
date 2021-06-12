!-----------------------------------------------------------------------
! this program finds spectrum and autocorrelation for 1d periodic data
! interactively modifies input parameters
! written for macintosh - viktor k. decyk, ucla
! copyright 1990, regents of the university of california
! character*8 code(nc), where nc is the number of input parameters
! character*8 cp(ncc), where ncc is the number of character parameters
! dimension ip(nci), ap(ncr)
! where nci, ncr are the number of integer and real parameters
! update: may 23, 2021
      program spectrum1
      use in1, only: idrun, indx, tend, dt, ntp, modesxp, nprec, ntdi,  &
     &modesxdi, ndirec, t0, ceng, fpname, pot1d, fdiname, deni1d
      use cmfield1
      use graf1
      use graf2
      implicit none
      integer, parameter :: psolve = 1
      integer, parameter :: ncc=0, nci=9, ncr=4, nc=ncc+nci+ncr, ns=2
      integer :: batch
      integer :: ndd, ndp, nrec, modesx, lts, its, nts
      integer :: kxmax, kxmin, kxmx, kxmn
      integer :: ntd, ntq, ntc, nplot, nvf, istyle, itwo, iwmax
      integer :: nx, nxh, nxe, nxeh, nx1
      integer :: it1, mtp, nt2d, it, nt, nt2, inft, nft, nfth, isign
      integer :: i, j, k, ii, i1, j1, irc, nts1, nts2, ntc1, ntc2
      integer :: jk, iw, iw1, iw2, ncd, ntr
      integer :: inx, nxd, nxdh, npot, dmap
      integer :: nf, nodesx, mt, nl
      real :: wmin, wmax, dw, pmin, qmin, fmin, anorm
      real :: at1, at2, dtp, dnx, dts, tt0, dkx, ak
      double precision :: sum1, sum2
      complex :: zt1, zsc, zst
      integer, dimension(1) :: itm
      integer, dimension(nci) :: ip
      real, dimension(ncr) :: ap
      real, dimension(:), pointer :: potb, pots, potr, pk
      real, dimension(:), pointer :: potd
      complex, dimension(:), pointer :: zpotx, zpoty
      real, dimension(:,:), pointer :: pwk, pck, pak, pkw, pckw
      real, dimension(:), pointer :: potc, time
      complex, dimension(:,:), pointer ::  zpott
      real, dimension(:), pointer :: g
      integer, dimension(:), pointer :: mixup, mixupt
      complex, dimension(:), pointer :: sct, sctt
      real, dimension(:), pointer :: wm, p, pc, ps, pcs
!
      character(len=64) :: dmetaf
      character(len=1) :: dtype
      character(len=60) :: prompt
      character(len=32) :: fname, ftname
      character(len=10) :: cdrun0
      character(len=20) :: cdrun, runid
      character(len=64) :: label, chr
      character(len=10), dimension(4) :: chrs
      character(len=4), dimension(2) :: chws
      character(len=8), dimension(nc) :: code
      character(len=8), dimension(1) :: cp
!
      equivalence (lts, ip(1)), (its, ip(2)), (nts, ip(3))
      equivalence (kxmin, ip(4)), (kxmax, ip(5))
      equivalence (ntq, ip(6)), (ntc, ip(7)), (nplot, ip(8))
      equivalence (dmap, ip(9))
      equivalence (wmin, ap(1)), (wmax, ap(2)), (dw, ap(3))
      equivalence (fmin, ap(4))
!
! define namelist
      namelist /inspect1/ batch, dmetaf, dtype, lts, its, nts, kxmin,   &
     &kxmax, ntd, ntc, nplot, nvf, ntr, wmin, wmax, dw, dmap, fmin
!
  981 format (' ACTUAL LENGTH OF DATA = ',i6,' DATA EXPECTED = ',i6)
  982 format (' ERROR IN PARAMETERS, LTS,ITS,NTS,NT = ',4i7)
  983 format (' ERROR IN MODE PARAMETERS, KXMIN,KXMAX,MODESX = ',3i7)
  985 format (' ERROR IN DISPLAY PARAMETERS, NTC,NTD,NTS = ',3i7)
  987 format (' ERROR IN FREQUENCY PARAMETERS, IW,IWMAX = ',2i7)
!
  991 format (' <ENERGY> = ',e14.7,',')
  992 format (', SUM = ',e14.7)
  996 format (' KX ',i4,', AKX=',f7.4)
! istyle = (0,1) = use (brief,full) menu style
      data istyle /1/
      data chrs /'   REAL   ','IMAGINARY ',' POSITIVE ',' NEGATIVE '/
      data chws /' W>0','W<0' /
      data code /'LTS   ','ITS   ','NTS   ','KMIN  ','KMAX  ','NTD  ',  &
     &'NTC   ','NPLOT ','DMAP  ','WMIN  ','WMAX  ','DW    ','FMIN  '/
      data itwo /2/
      data pmin,qmin /1.0e-14,1.0e-6/
      data nl /128/
      data iwmax /5001/
! default namelist values
! batch = (0,1) = (interactive,batch)
      data batch /0/
! lts,its,nts = initial data point, increment, and number of data points
      data lts,its,nts /1,1,1/
! kxmin,kxmax = initial mode number and number of modes
      data kxmin,kxmax /0,1/
! ntq = number of points in time display
! ntc =  number of lag times in correlation
      data ntq, ntc /1,0/
! wmin,wmax,dw = initial, final frequency and increment
      data wmin,wmax,dw /0.,2.0,.01/
! nplot = number of plots per page
      data nplot /4/
! ntr = display raw data every ntr steps 
      data ntr /0/
! dmap = (0,1) = (no,yes) display spectrum map
      data dmap /0/
      fmin = qmin
! open graphics device
      call GROPEN0(1)
      call STPALIT(3)
! check if running in batch mode
      open(unit=8,file='inspect1',form='formatted',status='old',iostat= &
     &irc)
      if (irc==0) then
         read (8,inspect1,iostat=irc)
         if ((irc==0).and.(batch/=0)) chr = dmetaf
         ntq = ntd
         rewind 8
      endif
! get file name
      label = ' '
    5 prompt = 'enter runid, or q to quit:'
      if (batch==0) call GTINPUT(label,prompt,chr,irc)
      if ((chr=='q').or.(chr=='Q')) go to 190
      nf = 1
      i = index(chr,'/')
      if (i > 1) nf = 2
      fname = 'diag1.'//chr(i+1:)
! attempt to open diagnostic metafile for given runid
      do j = 1, nf
         if (j==2) then
            write (cdrun0,'(i10)') idrun
            cdrun0 = adjustl(cdrun0)
            if (ndp==0) then
               nodesx = modesxp
               ftname = fpname
            else if (ndd==0) then
               nodesx = modesxdi
               ftname = fdiname
               cdrun0 = trim(cdrun0)//'d'
            endif
            nxh = 2**(indx-1)
            if (nodesx > nxh) nodesx = nxh
            allocate(zpoty(nodesx))
            fname = 'diag1.'//chr(1:i-1)
            close(unit=19)
            chr = ftname
         endif
         open(unit=19,file=trim(fname),form='formatted',status='old',   &
     &iostat=irc)
         if (irc /= 0) then
            label = 'open error for '//trim(fname)
            if (batch==0) then
               go to 5
            else
               write (*,*) label
               go to 190
            endif
         endif
! attempt to read diagnostic data from diagnostic namelists
! ndd = 0 means ion density data not found
         read (19,deni1d,iostat=ndd)
         rewind 19
! ndp = 0 means potential data not found
         read (19,pot1d,iostat=ndp)
         rewind 19
! if both diagnostics found, ask user to choose
         if ((ndp==0).and.(ndd==0)) then
            label = ' '
            if (nf==2) then
               label = 'for second file'
               if (j==2) label = 'for first file'
            endif
            prompt = 'enter p or d for potential or ion density diagnost&
     &ic:'
    6       if (batch==0) call GTINPUT(label,prompt,ftname,irc)
            if ((ftname=='q').or.(ftname=='Q')) then
               close(unit=19)
               label = ' '
               go to 5
            endif
! if potential is chosen set ion density flag to ndd = -1
            if ((ftname=='p').or.(ftname=='P')) then
               ndd = -1
            else if ((ftname=='d').or.(ftname=='D')) then
! if ion density is chosen set potential flag ndp = -1
               ndp = -1
            else
               label = 'Invalid string'
               go to 6
            endif
         else if ((ndp /= 0).and.(ndd /=0)) then
            label = 'No valid Namelists found'
            close(unit=19)
            if (batch==0) then
               go to 5
            else
               write (*,*) label
               go to 190
            endif
         endif
      enddo
! create string from idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
! reread chosen diagnostic data from diagnostic namelist
! ntp becomes generic time step for both potential and ion density
      if (ndp==0) then
         read (19,pot1d,iostat=ndp)
         modesx = modesxp
         nrec = nprec
         fname = fpname
         anorm = real((2**indx))/ceng
      else if (ndd==0) then
         read (19,deni1d,iostat=ndd)
         modesx = modesxdi
         nrec = ndirec; ntp = ntdi
         cdrun = trim(cdrun)//'d'
         fname = fdiname
         anorm = 1.0
      endif
      if (nf==2) cdrun = trim(cdrun)//'-'//trim(cdrun0)
      runid = 'spectrum.'//trim(cdrun)
! calculate initial parameters
      nx = 2**indx
      nxh = nx/2
      nxe = nx + 2
      nxeh = nxe/2
      nx1 = nx + 1
      if (modesx > nxh) modesx = nxh
      kxmin = 0
      kxmax = modesx - 1
      dnx = 6.28318530717959/real(nx)
      dtp = dt*real(ntp)
! time parameters
      if (nrec > 0) then
         mtp = nrec
      else
         it1 = (tend - t0)/dt + .0001
         mtp = (it1 - 1)/ntp + 1
      endif
      nt2d = mtp + mtp
! allocate space data
      allocate(zpotx(modesx))
      allocate(potb(modesx),pots(modesx))
      allocate(potr(modesx))
      allocate(pk(modesx))
      allocate(pwk(modesx,ns),pck(modesx,ns),pak(modesx,ns))
! allocate time data
      allocate(potc(nt2d),time(mtp))
      allocate(zpott(mtp,modesx),stat=npot)
! allocate data for periodic boundary conditions
      inx = indx
      nxd = nx
      nxdh = nxd/2
! allocate spatial data
      allocate(potd(nxd))
! allocate data for ffts
      allocate(mixup(nxh),sct(nxh))
! prepare fft tables
      call mfft1_init(mixup,sct,inx)
! enter raw data display parameter
      label = ' '
      prompt = 'enter N to display raw data every N steps (0 for none)'
      if (batch==0) then
         call GTINPUT(label,prompt,ftname,irc)
         if (irc==1) go to 5
         read (ftname,*,iostat=irc) ntr
         if (irc /= 0) ntr = 1
      endif
! clear screen
      call CLRSCRN
      if (ntr > 0) then
! set number of plots per page
         call SETNPLT(4,irc)
      endif
! get wave numbers
      call SAK1(pk,nx,modesx,modesx)
      zsc = cmplx(999.0,999.0)
      zst = cmplx(2.0,2.0)
      if (ndp==0) then
         label = ' POTENTIAL ENERGY VERSUS |K|'
      else if (ndd==0) then
         label = ' ION DENSITY**2 VERSUS |K|'
      endif
! read diagnostic input
      call fsopen1(11,trim(fname))
      nt = 1
      if (nt /= 1) go to 20
      if (nf==2) then
         call fsopen1(12,trim(chr))
         mt = 1
         if (mt /= 1) go to 20
      endif
      zpotx = cmplx(999.0,-999.0)
   10 call freadc1(11,zpotx,modesx)
      nt = nt + 1
      if (ircfio /= 0) go to 20
      if (nf==2) then
         call freadc1(12,zpoty,nodesx)
         if (ircfio /= 0) go to 20
         zpotx = zpotx - zpoty
      endif
      it = nt - 1
      if (npot==0) zpott(it,:) = zpotx
! create time array
      time(it) = t0 + dtp*real(it - 1)
! calculate potential energy
      if (ndp==0) then
         call WEL1(zpotx,potb,ceng,potc(it),nx,modesx,modesx)
! calculate average of potential energy
         if (it==1) then
            sum1 = potc(it)
            pots = potb
         else
            sum1 = sum1 + potc(it)
            pots = pots + potb
         endif
      else if (ndd==0) then
         call SQ2POT1(zpotx,potr,nx,modesx,modesx)
      endif
! display raw data
      if (ntr > 0) then
      it1 = (it-1)/ntr
      if ((it-1)==ntr*it1) then
         write (chr,*) 'NT = ', it
         call mwrmodes1(potd,zpotx,nx,modesx)
         isign = 1
         call mfft1r(potd,isign,mixup,sct,indx)
! display potential
         if (ndp==0) then
            call dscaler1(potd,' POTENTIAL VS X',it,999,0,nx,irc)
! display ion density
         else if (ndd==0) then
            call dscaler1(potd,' ION DENSITY VS X',it,999,2,nx,irc)
         endif
! display potential energy
         if (ndp==0) then
            call dscaler1(potb,' POTENTIAL ENERGY VS MODE NUMBER',it,999&
     &,1,modesx,irc)
! display ion density
         else if (ndd==0) then
            call dscaler1(potr,' ION DENSITY**2 VS MODE NUMBER',it,999, &
     &1,modesx,irc)
         endif
         if (irc > 127) then
            ntr = irc - 128
            if (ntr==0) call CLRSCRN
            irc = 0
         endif
         if (irc==1) go to 190
      endif
      endif
      if (it < mtp) go to 10
   20 call RSTSCRN
      nt = it
      if (nt /= mtp) then
         write (*,981) nt, mtp
      endif
      nt = min(nt,mtp)
      call SETNPLT(1,irc)
! display potential energy
      if (ndp==0) then
         label = ' POTENTIAL ENERGY VERSUS TIME'
         chr = ' RUNID='
         chr = trim(chr)//trim(cdrun)
         call DISPS(potc(1),label,time(1),time(nt),999,2,nt,chr,irc)
         if (irc==1) go to 180
! display average potential energy versus |k|
         pots = pots/real(nt)
         sum1 = sum1/real(nt)
         write (ftname,991) sum1
         chr = ' RUNID='//trim(cdrun)
         chr = trim(ftname)//chr
         label = ' TIME-AVERAGED POTENTIAL ENERGY VERSUS |K|'
         call DISPC(pots,pk(1),label,zsc,zst,2,modesx,modesx,1,chr,chrs,&
     &irc)
         if (irc==1) go to 180
      endif
      if (irc==1) go to 180
! quit if no memory available
      if (npot/=0) then
         write (*,*) 'fatal zpott allocation error=', npot
         go to 190
      endif
! allocate correlation time data
      nt2 = nt + nt
      inft = 0
      nft = 1
   30 inft = inft + 1
      nft = 2*nft
      if (nt2 > nft) go to 30
      nfth = nft/2
      allocate(g(2*nft))  
      allocate(mixupt(nft),sctt(nfth))
      nts = nt
      ntq = nts
      ntc = nt/3
      if (ntq > nt) ntq = nt
! create fft table for correlations
      call WFFT1CINIT(mixupt,sctt,inft,nft,nfth)
! allocate frequency data 
      allocate(wm(iwmax),p(2*iwmax),pc(2*iwmax))
      allocate(ps(2*iwmax),pcs(2*iwmax))
      allocate(pkw(modesx,2*iwmax),pckw(modesx,2*iwmax))
! enter analysis parameters
! nts = number of data points used
! ntd = number of points used in time display
   50 if (batch==0) then
         call MENUCR1(code,cp,ip,ap,nc,ncc,nci,ncr,istyle,irc)
         ntd = ntq
      else
         read (8,inspect1,iostat=irc)
         close(unit=8)
      endif
      if (irc==1) go to 190
      if (irc==3) go to 50
      prompt = ' Hit carrriage return or enter key to continue '
! make sure time display parameters make sense
      if ((lts < 1).or.(((nts - 1)*its + lts) > nt)) then
         write (label,982) lts, its, nts, nt
         if (batch==0) then
            call GTINPUT(label,prompt,ftname,irc)
            go to 50
         else
            write (*,*) label
            go to 190
         endif
      endif
      nts1 = nts + 1
      nts2 = nts + nts
      dts = dtp*real(its)
      tt0 = t0 + dtp*real(lts - 1)
! make sure mode number parameters make sense
      if ((kxmin < 0) .or. (kxmax >= modesx)) then
         write (label,983) kxmin, kxmax, modesx
         if (batch==0) then
            call GTINPUT(label,prompt,ftname,irc)
            go to 50
         else
            write (*,*) label
            go to 190
         endif
!        stop
      endif
! make sure lag time and time display parameters make sense
      if ((ntd < 0).or.(ntd > nts).or.(ntc < 0).or.(ntc > nts)) then
         write (label,985) ntc, ntd, nts
         if (batch==0) then
            call GTINPUT(label,prompt,ftname,irc)
            go to 50
         else
            write (*,*) label
            go to 190
         endif
      endif
      ntc1 = ntc + 1
      ntc2 = 2*ntc
! make sure frequency parameters make sense
      if (dw==0.0) then
         iw = 0
      else
         iw = (wmax - wmin)/dw + 1.001
      endif
      if (iw < 0) then
         write (label,987) iw, iwmax
         if (batch==0) then
            call GTINPUT(label,prompt,ftname,irc)
            go to 50
         else
            write (*,*) label
            go to 190
         endif
      else if (iw > iwmax) then
         deallocate(wm,p,pc,ps,pcs)
         iwmax = iw
         write (*,*) 'reallocating frequency data, iwmax=', iwmax
         allocate(wm(iwmax),p(2*iwmax),pc(2*iwmax))
         allocate(ps(2*iwmax),pcs(2*iwmax))
         allocate(pkw(modesx,2*iwmax),pckw(modesx,2*iwmax))
      endif
      iw2 = 2*iw
      iw1 = iw + 1
      do 60 k = 1, iw
      wm(k) = wmin + dw*real(k - 1)
   60 continue
! set number of plots per page
      call SETNPLT(nplot,irc)
! display analysis parameters
      call WPCORR1(runid,indx,ntp,modesx,psolve,t0,tend,dt,lts,its,nts, &
     &kxmin,kxmax,ntd,ntc,wmin,wmax,dw,irc)
      if ((irc==1).or.(irc==3)) go to 50
      kxmn = kxmin + 1
      kxmx = kxmax + 1
! clear screen
      call CLRSCRN
      if (iw==0) go to 120
!
! begin main loop to find dispersion relation
!
      ps = 0.0; pcs = 0.0
      do 100 j = kxmn, kxmx
      j1 = j - 1
      jk = j1 - kxmn + 2
      dkx = dnx*real(j1)
      ak = dkx
! load data in proper format
      do 70 i = 1, nts
      ii = lts + (i - 1)*its
! potential
      zt1 = zpott(ii,j)
      potc(i) = real(zt1)
      potc(i+nts) = aimag(zt1)
   70 continue
! perform frequency analysis of raw data
      call SPECT(potc,wm,p,tt0,dts,nts,iw,nts2,iw2)
! add contributions from different modes
      ps = ps + p
      pkw(jk,1:iw2) = p(1:iw2)
! extract maximum frequency value
      itm = maxloc(p(1:iw))
      at1 = wm(itm(1))
      itm = maxloc(p(iw1:iw2))
      at2 = wm(itm(1))
      pwk(jk,1) = at1
      pwk(jk,2) = at2
      pak(jk,:) = ak
      if (ntc==0) then
         go to 100
      endif
! auto-correlation function
      call CORRE(potc,g,g,potc,mixupt,sctt,inft,nts,ntc,nts2,nft,nfth,  &
     &ntc2)
! perform frequency analysis of correlated data
      call SPECT(potc,wm,pc,tt0,dts,ntc,iw,ntc2,iw2)
! add contributions from different modes
      pcs = pcs + pc
      pckw(jk,1:iw2) = pc(1:iw2)
! extract maximum frequency value
      itm = maxloc(pc(1:iw))
      at1 = wm(itm(1))
      itm = maxloc(pc(iw1:iw2))
      at2 = wm(itm(1))
      pck(jk,1) = at1
      pck(jk,2) = at2
  100 continue
!
! end main loop to find dispersion relation
!
! display frequency information
      call SETNPLT(1,irc)
!
      if (dmap > 0) then
         label = ' SPECTRUM OMEGA VERSUS K'
         if (ndp==0) then
            chr = ' POTENTIAL, W > 0, RUNID='//trim(cdrun)
         else if (ndd==0) then
            chr = ' ION DENSITY, W > 0, RUNID='//trim(cdrun)
         endif
         jk = kxmx - kxmn + 1
         at1 = dnx*real(kxmn-1)
         at2 = dnx*real(kxmx-1)
         call CARPETL(pkw(1,1),label,at1,at2,wm(1),wm(iw),999,1,jk,iw,  &
     &modesx,chr,nl,irc)
         if ((irc==1).or.(irc==3)) go to 50
         if (ndp==0) then
            chr = ' POTENTIAL, W < 0, RUNID='//trim(cdrun)
         else if (ndd==0) then
            chr = ' ION DENSITY, W < 0, RUNID='//trim(cdrun)
         endif
         call CARPETL(pkw(1,iw1),label,at1,at2,wm(1),wm(iw),999,1,jk,iw,&
     &modesx,chr,nl,irc)
         if ((irc==1).or.(irc==3)) go to 50
! display logarithm
         do i = 1, iw2
         do j = 1, jk
         at1 = pkw(j,i)
         if (at1 < fmin) at1 = fmin
         pkw(j,i) = alog(at1)
         enddo
         enddo
         label = ' LOG SPECTRUM OMEGA VERSUS K'
         if (ndp==0) then
            chr = ' POTENTIAL, W > 0, RUNID='//trim(cdrun)
         else if (ndd==0) then
            chr = ' ION DENSITY, W > 0, RUNID='//trim(cdrun)
         endif
         call CARPETL(pkw(1,1),label,at1,at2,wm(1),wm(iw),999,2,jk,iw,  &
     &modesx,chr,nl,irc)
         if ((irc==1).or.(irc==3)) go to 50
         if (ndp==0) then
            chr = ' POTENTIAL, W < 0, RUNID='//trim(cdrun)
         else if (ndd==0) then
            chr = ' ION DENSITY, W < 0, RUNID='//trim(cdrun)
         endif
         call CARPETL(pkw(1,iw1),label,at1,at2,wm(1),wm(iw),999,2,jk,iw,&
     &modesx,chr,nl,irc)
         if ((irc==1).or.(irc==3)) go to 50
      endif
!
! normalize and display integrated spectrum versus frequency
      sum1 = 0.0d0
      sum2 = 0.0d0
      sum1 = sum1 + sum2 - ps(iw+1)
      write (ftname,992) sum1
      label = ' INTEGRATED SPECTRUM VERSUS OMEGA'//trim(ftname)
      if (ndp==0) then
         chr = ' POTENTIAL, RUNID='//trim(cdrun)
      else if (ndd==0) then
         chr = ' ION DENSITY, RUNID='//trim(cdrun)
      endif
      ncd = 58
      call DISPR(ps(1),label,wm(1),wm(iw),999,1,0,iw,iw,itwo,chr(1:ncd),&
     &chrs(3),irc)
      if ((irc==1).or.(irc==3)) go to 50
! display maximum frequency
      zsc = cmplx(999.0,999.0)
      zst = cmplx(2.0,2.0)
      label = ' FREQUENCY W VERSUS |K|'
      it1 = (kxmx - kxmn + 1)
      call DISPC(pwk,pak(1,1),label,zsc,zst,2,it1,modesx,ns,chr,chws,   &
     &irc)
      if ((irc==1).or.(irc==3)) go to 50
      if (ntc > 0) then
!
         if (dmap > 0) then
            label = ' CORRELATION SPECTRUM OMEGA VERSUS K'
            if (ndp==0) then
               chr = ' POTENTIAL, W > 0, RUNID='//trim(cdrun)
            else if (ndd==0) then
               chr = ' ION DENSITY, W > 0, RUNID='//trim(cdrun)
            endif
            jk = kxmx - kxmn + 1
            at1 = dnx*real(kxmn-1)
            at2 = dnx*real(kxmx-1)
            call CARPETL(pckw(1,1),label,at1,at2,wm(1),wm(iw),999,1,jk, &
     &iw,modesx,chr,nl,irc)
            if ((irc==1).or.(irc==3)) go to 50
            if (ndp==0) then
               chr = ' POTENTIAL, W < 0, RUNID='//trim(cdrun)
            else if (ndd==0) then
               chr = ' ION DENSITY, W < 0, RUNID='//trim(cdrun)
            endif
            call CARPETL(pckw(1,iw1),label,at1,at2,wm(1),wm(iw),999,1,jk&
     &,iw,modesx,chr,nl,irc)
            if ((irc==1).or.(irc==3)) go to 50
! display logarithm
            do i = 1, iw2
            do j = 1, jk
            at1 = sqrt(pckw(j,i))
            if (at1 < fmin) at1 = fmin
            pckw(j,i) = alog(at1)
            enddo
            enddo
            label = ' LOG CORRELATION SPECTRUM OMEGA VERSUS K'
            if (ndp==0) then
               chr = ' POTENTIAL, W > 0, RUNID='//trim(cdrun)
            else if (ndd==0) then
               chr = ' ION DENSITY, W > 0, RUNID='//trim(cdrun)
            endif
            call CARPETL(pckw(1,1),label,at1,at2,wm(1),wm(iw),999,2,jk, &
     &iw,modesx,chr,nl,irc)
            if ((irc==1).or.(irc==3)) go to 50
            if (ndp==0) then
               chr = ' POTENTIAL, W < 0, RUNID='//trim(cdrun)
            else if (ndd==0) then
               chr = ' ION DENSITY, W < 0, RUNID='//trim(cdrun)
            endif
            call CARPETL(pckw(1,iw1),label,at1,at2,wm(1),wm(iw),999,2,jk&
     &,iw,modesx,chr,nl,irc)
            if ((irc==1).or.(irc==3)) go to 50
         endif
!
! normalize and display integrated correlation spectrum versus frequency
         sum1 = 0.0d0
         sum2 = 0.0d0
         sum1 = sum1 + sum2 - pcs(iw+1)
         write (ftname,992) sum1
         label = ' INTEGRATED CORRELATION SPECTRUM VS OMEGA'//trim(ftnam&
     &e)
         if (ndp==0) then
            chr = ' POTENTIAL CORRELATION, RUNID='//trim(cdrun)
         else if (ndd==0) then
            chr = ' ION DENSITY CORRELATION, RUNID='//trim(cdrun)
         endif
         ncd = 58
         call DISPR(pcs(1),label,wm(1),wm(iw),999,1,0,iw,iw,itwo,       &
     &chr(1:ncd),chrs(3),irc)
         if ((irc==1).or.(irc==3)) go to 50
! display maximum correlation frequency
         label = ' FREQUENCY W VERSUS |K|'
         call DISPC(pck,pak(1,1),label,zsc,zst,2,it1,modesx,ns,chr,chws,&
     &irc)
         if ((irc==1).or.(irc==3)) go to 50
      endif
      call SETNPLT(nplot,irc)
!
! begin main loop to display individual modes
!
      irc = 0
      j = kxmn - 1
  120 if (irc < 128) then
         j = j + 1
      else
         j = irc - 127
         if (j < kxmn) then
            kxmn = 1
         else if (j > kxmx) then
            kxmx = modesx
            j = modesx
         endif
         irc = 0
      endif
      if ((j < kxmn) .or. (j > kxmx)) go to 160
      j1 = j - 1
      dkx = dnx*real(j1)
      ak = dkx
! load data in proper format
      do 130 i = 1, nts
      ii = lts + (i - 1)*its
! potential
      zt1 = zpott(ii,j)
      potc(i) = real(zt1)
      potc(i+nts) = aimag(zt1)
      time(i) = tt0 + dts*real(i - 1)
  130 continue
! time display of raw data
      if (ntd > 0) then
         if (ndp==0) then
            label = ' POTENTIAL VERSUS TIME, RUNID='//trim(cdrun)
         else if (ndd==0) then
            label = ' ION DENSITY VERSUS TIME, RUNID='//trim(cdrun)
         endif
         write (chr,996) j1, dkx
         ncd = 58
         call DISPR(potc(1),label,time(1),time(ntd),999,0,0,ntd,nts,    &
     &itwo,chr(1:ncd),chrs,irc)
         if (irc==1) go to 180
         if (irc >= 128) go to 120
! plot complex phase
         if (ndp==0) then
            label = ' REAL VERSUS IMAGINARY PART OF POTENTIAL'
         else if (ndd==0) then
            label = ' REAL VERSUS IMAGINARY PART OF ION DENSITY'
         endif
         write (chr,996) j1, dkx
         ncd = 20
         zsc = cmplx(999.,999.)
         zst = cmplx(0.,0.)
         call DISPC(potc,potc(nts1),label,zsc,zst,1,ntd,ntd,1,          &
     &chr(1:ncd),' ',irc)
         if (irc==1) go to 180
         if (irc >= 128) go to 120
! log display of raw data
         do 140 i = 1, nts
         i1 = nts + i
         at1 = sqrt(potc(i)**2 + potc(i1)**2)
         if (at1 < pmin) at1 = pmin
         g(i) = alog(at1)
  140    continue
         if (ndp==0) then
            label = ' LN |POTENTIAL| VERSUS TIME, RUNID='//trim(cdrun)
         else if (ndd==0) then
            label = ' LN |ION DENSITY| VERSUS TIME, RUNID='//trim(cdrun)
         endif
         write (chr,996) j1, dkx
         ncd = 20
         call DISPS(g(1),label,time(1),time(ntd),999,2,ntd,chr(1:ncd),  &
     &irc)
         if (irc==1) go to 180
         if (irc >= 128) go to 120
      endif
      if (iw > 0) then
! perform frequency analysis of raw data
         call SPECT(potc,wm,p,tt0,dts,nts,iw,nts2,iw2)
! power spectrum of raw data
         if (ndp==0) then
            label = ' POTENTIAL VERSUS OMEGA, RUNID='//trim(cdrun)
         else if (ndd==0) then
            label = ' ION DENSITY VERSUS OMEGA, RUNID='//trim(cdrun)
         endif
         write (chr,996) j1, dkx
         ncd = 20
         call DISPR(p(1),label,wm(1),wm(iw),999,1,0,iw,iw,itwo,         &
     &chr(1:ncd),chrs(3),irc)
         if (irc==1) go to 180
         if (irc >= 128) go to 120
      endif
      if (ntc==0) go to 120
! auto-correlation function
      call CORRE(potc,g,g,potc,mixupt,sctt,inft,nts,ntc,nts2,nft,nfth,  &
     &ntc2)
      if (ntd > 0) then
! time display of correlated data
         it1 = ntd
         if (it1 > ntc) it1 = ntc
         if (ndp==0) then
            label = ' POT CORRELATION VERSUS TIME, RUNID='//trim(cdrun)
         else if (ndd==0) then
            label = ' ION CORRELATION VERSUS TIME, RUNID='//trim(cdrun)
         endif
         write (chr,996) j1, dkx
         ncd = 20
         call DISPR(potc(1),label,time(1),time(it1),999,0,0,it1,ntc,itwo&
     &,chr(1:ncd),chrs,irc)
         if (irc==1) go to 180
         if (irc >= 128) go to 120
! plot complex phase
         if (ndp==0) then
            label = ' REAL VERSUS IMAGINARY PART OF POT CORRELATION'
         else if (ndd==0) then
            label = ' REAL VERSUS IMAGINARY PART OF ION CORRELATION'
         endif
         write (chr,996) j1, dkx
         ncd = 20
         zsc = cmplx(999.,999.)
         zst = cmplx(0.,0.)
         call DISPC(potc,potc(ntc1),label,zsc,zst,1,it1,it1,1,chr(1:ncd)&
     &,' ',irc)
         if (irc==1) go to 180
         if (irc >= 128) go to 120
! log display of correlated data
         do 150 i = 1, ntc
         i1 = ntc + i
         at1 = sqrt(potc(i)**2 + potc(i1)**2)
         if (at1 < pmin) at1 = pmin
         g(i) = alog(at1)
  150    continue
         if (ndp==0) then
            label = ' LN |POT CORRELATION| VS TIME, RUNID='//trim(cdrun)
         else if (ndd==0) then
            label = ' LN |ION CORRELATION| VS TIME, RUNID='//trim(cdrun)
         endif
         write (chr,996) j1, dkx
         ncd = 20
         call DISPS(g(1),label,time(1),time(it1),999,2,it1,chr(1:ncd),  &
     &irc)
         if (irc==1) go to 180
         if (irc >= 128) go to 120
      endif
      if (iw > 0) then
! perform frequency analysis of correlated data
         call SPECT(potc,wm,pc,tt0,dts,ntc,iw,ntc2,iw2)
! power spectrum of correlated data
         if (ndp==0) then
            label = ' POT CORRELATION VS OMEGA, RUNID='//trim(cdrun)
         else if (ndd==0) then
            label = ' ION CORRELATION VS OMEGA, RUNID='//trim(cdrun)
         endif
         write (chr,996) j1, dkx
         ncd = 20
         call DISPR(pc(1),label,wm(1),wm(iw),999,1,0,iw,iw,itwo,        &
     &chr(1:ncd),chrs(3),irc)
         if (irc==1) go to 180
         if (irc >= 128) go to 120
      endif
      go to 120
  160 continue
!
! end main loop to display individual modes
!
  180 continue
      if (batch==0) go to 50
! close graphics device
  190 call GRCLOSE
      end program
