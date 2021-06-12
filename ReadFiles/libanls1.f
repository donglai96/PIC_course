c library for simulation post-processors
c copyright 1990, regents of the university of california
c update: july 18, 2011
      subroutine MENUCR1(code,cp,ip,ap,nc,ncc,nci,ncr,istyle,irc)
c this program sets up menu to enter parameters for correlation program
c the parameter names must be in table code (nc of them), and created
c in order of type character*8 (ncc of them), integer (nci of them), and
c real (ncr of them).  the values go into the array cp, ip, and ap,
c for character, integer, and real data types, respectively.
c input arguments: code, istyle
c code = variable symbol character array
c cp = array containing character variables
c ip = array containing integer variables
c ap = array containing real variables
c nc = total number of variables
c ncc = number of character variables
c nci = number of integer variables
c ncr = number of real variables
c istyle = (0,1) = use (brief,full) menu style
c irc = return code (0 = normal return)
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c idwk = workstation identifier
c rx, ry = ndc coordinates of upper-right corner of workstation window
c idstr = string device number, 0 if no string device available
      character*8 code(nc)
      character*8 cp(*)
      dimension ip(nci), ap(ncr)
      character*40 prompt, chr
      character*24 input
      character*20 cval
      character*8 codev, crun, cquit
c     equivalence (lts, ip(1)), (its, ip(2)), (nts, ip(3))
c     equivalence (kmin, ip(4)), (kmax, ip(5))
c     equivalence (ntd, ip(6)), (ntc, ip(7)), (nplot, ip(8))
c     equivalence (dmap, ip(9))
c     equivalence (wmin, ap(1)), (wmax, ap(2)), (dw, ap(3))
c     equivalence (fmin,ap(4))
   91 format (a24)
   92 format (' LTS=',i6,' ITS=',i6,' NTS=',i6)
   93 format (' KMIN=',i6,' KMAX=',i6)
   94 format (' NTD=',i6,' NTC=',i6)
   95 format (' WMIN=',f8.4,' WMAX=',f8.4,' DW=',f8.4)
   96 format (' NPLOT=',i6)
   97 format (' DMAP=',i6,' FMIN=',e10.3)
      data crun,cquit /'RUN     ','QUIT    '/
      data prompt /' enter parameter, or type quit or run   '/
      data csize /0.024/
      irc = 0
c exit if not input device
      if (idstr.eq.0) return
      if (ry.lt..3) csize = .032
      chh = ry*csize
c space = vertical spacing between characters
      space = 1.5*chh
      ay = ry
c exit if character height too large
      if (ay.lt.space) return
      ax = 0.
c define and select identity transformation
      call idtran
c set default plotting parameters
c 1 = default font, 0 = string precision
      call dfplps (chh,1,0)
c display current parameters
c clear workstation, always
   20 call gclrwk(idwk,1)
      ay = ry - space
      if ((istyle.eq.1).and.(ay.ge.0.)) call gtx(ax,ay,'    TIME INPUT P
     1ARAMETERS')
      write (chr,92) ip(1), ip(2), ip(3)
      ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,chr)
      ay = ay - space
      if ((istyle.eq.1).and.(ay.ge.0.)) call gtx(ax,ay,'    MODE PARAMET
     1ERS')
      write (chr,93) ip(4), ip(5)
      ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,chr)
      ay = ay - space
      if ((istyle.eq.1).and.(ay.ge.0.)) call gtx(ax,ay,'    TIME DISPLAY
     1 AND CORRELATION PARAMETERS')
      write (chr,94) ip(6), ip(7)
      ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,chr)
      ay = ay - space
      if ((istyle.eq.1).and.(ay.ge.0.)) call gtx(ax,ay,'    FREQUENCY PA
     1RAMETERS')
      write (chr,95) ap(1), ap(2), ap(3)
      ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,chr)
      ay = ay - space
      if ((istyle.eq.1).and.(ay.ge.0.)) call gtx(ax,ay,'    DISPLAY PARA
     1METER')
      write (chr,96) ip(8)
      ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,chr)
      ay = ay - space
      if ((istyle.eq.1).and.(ay.ge.0.)) call gtx(ax,ay,'    SPECTRUM MAP
     1 PARAMETER')
      write (chr,97) ip(9), ap(4)
      ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,chr)
c display prompt
   30 ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,prompt)
c read input
c update workstation, perform
      call guwk(idwk,1)
      input = ' '
c request string
      call grqst(idwk,idstr,istat,lens,input)
c     read (5,91,end=35) input
c abort
      if ((input(1:1).eq.'q').or.(input(1:1).eq.'Q')) then
         if (input(2:2).eq.' ') then
            irc = 1
            return
         endif
      endif
c parse input string
      call vpars(input,codev,ic,val,id,cval,ird)
c re-display current page
      if (ird.eq.1) go to 20
c check if known command
c normal exit
      if (codev.eq.crun) then
         return
c abort
      elseif (codev.eq.cquit) then
         irc = 1
         return
      endif
c check if known variable
      call findvn(code,codev,nc,nvar)
c re-display current page
      if (nvar.eq.0) go to 20
c help request
      if (ird.gt.1) then
         call helpcr1(nvar,ax,ay,space)
         go to 30
      endif
      nct = ncc + nci
      if (nvar.gt.nct) then
         ap(nvar-nct) = val
      elseif (nvar.le.ncc) then
         cp(nvar) = cval
      else
         ip(nvar-ncc) = val
      endif
      go to 20
      end
      subroutine MENUCRV1(code,cp,ip,ap,nc,ncc,nci,ncr,istyle,irc)
c this program sets up menu to enter parameters for correlation program
c the parameter names must be in table code (nc of them), and created
c in order of type character*8 (ncc of them), integer (nci of them), and
c real (ncr of them).  the values go into the array cp, ip, and ap,
c for character, integer, and real data types, respectively.
c for vector data
c input arguments: code, istyle
c code = variable symbol character array
c cp = array containing character variables
c ip = array containing integer variables
c ap = array containing real variables
c nc = total number of variables
c ncc = number of character variables
c nci = number of integer variables
c ncr = number of real variables
c istyle = (0,1) = use (brief,full) menu style
c irc = return code (0 = normal return)
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c idwk = workstation identifier
c rx, ry = ndc coordinates of upper-right corner of workstation window
c idstr = string device number, 0 if no string device available
      character*8 code(nc)
      character*8 cp(*)
      dimension ip(nci), ap(ncr)
      character*40 prompt, chr
      character*24 input
      character*20 cval
      character*8 codev, crun, cquit
c     equivalence (lts, ip(1)), (its, ip(2)), (nts, ip(3))
c     equivalence (kmin, ip(4)), (kmax, ip(5))
c     equivalence (ntd, ip(6)), (ntc, ip(7)), (nplot, ip(8))
c     equivalence (nvf, ip(9)), (dmap, ip(10))
c     equivalence (wmin, ap(1)), (wmax, ap(2)), (dw, ap(3))
c     equivalence (fmin,ap(4))
   91 format (a24)
   92 format (' LTS=',i6,' ITS=',i6,' NTS=',i6)
   93 format (' KMIN=',i6,' KMAX=',i6)
   94 format (' NTD=',i6,' NTC=',i6)
   95 format (' WMIN=',f8.4,' WMAX=',f8.4,' DW=',f8.4)
   96 format (' NPLOT=',i6,' NVF=',i6)
   97 format (' DMAP=',i6,' FMIN=',e10.3)
      data crun,cquit /'RUN     ','QUIT    '/
      data prompt /' enter parameter, or type quit or run   '/
      data csize /0.024/
      irc = 0
c exit if not input device
      if (idstr.eq.0) return
      if (ry.lt..3) csize = .032
      chh = ry*csize
c space = vertical spacing between characters
      space = 1.5*chh
      ay = ry
c exit if character height too large
      if (ay.lt.space) return
      ax = 0.
c define and select identity transformation
      call idtran
c set default plotting parameters
c 1 = default font, 0 = string precision
      call dfplps (chh,1,0)
c display current parameters
c clear workstation, always
   20 call gclrwk(idwk,1)
      ay = ry - space
      if ((istyle.eq.1).and.(ay.ge.0.)) call gtx(ax,ay,'    TIME INPUT P
     1ARAMETERS')
      write (chr,92) ip(1), ip(2), ip(3)
      ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,chr)
      ay = ay - space
      if ((istyle.eq.1).and.(ay.ge.0.)) call gtx(ax,ay,'    MODE PARAMET
     1ERS')
      write (chr,93) ip(4), ip(5)
      ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,chr)
      ay = ay - space
      if ((istyle.eq.1).and.(ay.ge.0.)) call gtx(ax,ay,'    TIME DISPLAY
     1 AND CORRELATION PARAMETERS')
      write (chr,94) ip(6), ip(7)
      ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,chr)
      ay = ay - space
      if ((istyle.eq.1).and.(ay.ge.0.)) call gtx(ax,ay,'    FREQUENCY PA
     1RAMETERS')
      write (chr,95) ap(1), ap(2), ap(3)
      ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,chr)
      ay = ay - space
      if ((istyle.eq.1).and.(ay.ge.0.)) call gtx(ax,ay,'    DISPLAY PARA
     1METER')
      write (chr,96) ip(8), ip(9)
      ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,chr)
      ay = ay - space
      if ((istyle.eq.1).and.(ay.ge.0.)) call gtx(ax,ay,'    SPECTRUM MAP
     1 PARAMETER')
      write (chr,97) ip(10), ap(4)
      ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,chr)
c display prompt
   30 ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,prompt)
c read input
c update workstation, perform
      call guwk(idwk,1)
      input = ' '
c request string
      call grqst(idwk,idstr,istat,lens,input)
c     read (5,91,end=35) input
c abort
      if ((input(1:1).eq.'q').or.(input(1:1).eq.'Q')) then
         if (input(2:2).eq.' ') then
            irc = 1
            return
         endif
      endif
c parse input string
      call vpars(input,codev,ic,val,id,cval,ird)
c re-display current page
      if (ird.eq.1) go to 20
c check if known command
c normal exit
      if (codev.eq.crun) then
         return
c abort
      elseif (codev.eq.cquit) then
         irc = 1
         return
      endif
c check if known variable
      call findvn(code,codev,nc,nvar)
c re-display current page
      if (nvar.eq.0) go to 20
c help request
      if (ird.gt.1) then
         call helpcrv1(nvar,ax,ay,space)
         go to 30
      endif
      nct = ncc + nci
      if (nvar.gt.nct) then
         ap(nvar-nct) = val
      elseif (nvar.le.ncc) then
         cp(nvar) = cval
      else
         ip(nvar-ncc) = val
      endif
      go to 20
      end
      subroutine helpcr1(nvar,ax,ay,space)
c this subroutine displays help for correlation parameters
c input arguments: all
c nvar = variable number for requested help
c ax/ay = location of last character write
c space = vertical spacing between characters
      parameter(nhelp=13)
      character*48 chr(nhelp)
      data chr /' LTS = INITIAL DATA POINT USED                  ',
     1' ITS = INCREMENT BETWEEN DATA POINTS USED       ',
     2' NTS = NUMBER OF DATA POINTS USED               ',
     3' KMIN = MINIMUM MODE NUMBER ANALYZED            ',
     4' KMAX = MAXIMUM MODE NUMBER ANALYZED            ',
     7' NTD = NUMBER OF POINTS USED IN TIME DISPLAY    ',
     8' NTC = NUMBER OF LAG TIMES IN CORRELATION       ',
     9' NPLOT = NUMBER OF PLOTS PER PAGE               ',
     a' DMAP = (0,1) = (NO,YES) DISPLAY SPECTRUM MAP   ',
     b' WMIN = MINIMUM FREQUENCY USED IN POWER SPECTRUM',
     c' WMAX = MAXIMUM FREQUENCY USED IN POWER SPECTRUM',
     d' DW = FREQUENCY INCREMENT USED IN POWER SPECTRUM',
     e' FMIN = MINIMUM CUTOFF FOR LOG                  '/
      if ((nvar.lt.1).or.(nvar.gt.nhelp)) return
      ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,chr(nvar))
      end
      subroutine helpcrv1(nvar,ax,ay,space)
c this subroutine displays help for correlation parameters
c input arguments: all
c nvar = variable number for requested help
c ax/ay = location of last character write
c space = vertical spacing between characters
      parameter(nhelp=14)
      character*48 chr(nhelp)
      data chr /' LTS = INITIAL DATA POINT USED                  ',
     1' ITS = INCREMENT BETWEEN DATA POINTS USED       ',
     2' NTS = NUMBER OF DATA POINTS USED               ',
     3' KMIN = MINIMUM MODE NUMBER ANALYZED            ',
     4' KMAX = MAXIMUM MODE NUMBER ANALYZED            ',
     7' NTD = NUMBER OF POINTS USED IN TIME DISPLAY    ',
     8' NTC = NUMBER OF LAG TIMES IN CORRELATION       ',
     9' NPLOT = NUMBER OF PLOTS PER PAGE               ',
     a' NVF = VECTOR FIELD (A=1,E=2,B=3)               ',
     b' DMAP = (0,1) = (NO,YES) DISPLAY SPECTRUM MAP   ',
     c' WMIN = MINIMUM FREQUENCY USED IN POWER SPECTRUM',
     d' WMAX = MAXIMUM FREQUENCY USED IN POWER SPECTRUM',
     e' DW = FREQUENCY INCREMENT USED IN POWER SPECTRUM',
     f' FMIN = MINIMUM CUTOFF FOR LOG                  '/
      if ((nvar.lt.1).or.(nvar.gt.nhelp)) return
      ay = ay - space
      if (ay.ge.0.) call gtx(ax,ay,chr(nvar))
      end
      subroutine WPCORR1(runid,indx,ntp,modes,psolve,t0,tend,dt,lts,its,
     1nts,kmin,kmax,ntd,ntc,wmin,wmax,dw,irc)
c this subroutine displays parameters for correlation program
c input arguments: all except irc
c irc = return code (0 = normal return)
      character*20 runid
      integer psolve
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c idwk = workstation identifier
c rx, ry = ndc coordinates of upper-right corner of workstation window
      character*48 chr
   91 format (' SPECTRUM ANALYSIS FOR 1D DATA')
   92 format (' RUNID= ',a20,' INDX=',i3)
   93 format (' NTP=',i6,' MODES=',i5,' PSOLVE=',i2)
   94 format (' T0=',f8.1,' TEND=',f8.1,' DT=',f8.5)
   95 format (' LTS=',i6,' ITS=',i6,' NTS=',i6)
   96 format (' KMIN=',i6,' KMAX=',i6)
   97 format (' NTD=',i6,' NTC=',i6)
   98 format (' WMIN=',f8.4,' WMAX=',f8.4,' DW=',f8.4)
      data csize /0.024/
      irc = 0
      if (ry.lt..3) csize = .032
      chh = ry*csize
c space = vertical spacing between characters
      space = 1.5*chh
c initiate plot
c clear workstation, always
      call gclrwk(idwk,1)
c define and select identity transformation
      call idtran
c set default plotting parameters
c 1 = default font, 0 = string precision
      call dfplps (chh,1,0)
      ay = ry
      if (ay.lt.chh) return
      ax = csize*rx
      write (chr,91)
      ay = ay - space
c draw text
      call gtx(ax,ay,chr(1:30))
      write (chr,92) runid, indx
      ay = ay - space
c draw text
      call gtx(ax,ay,chr(1:37))
      write (chr,93) ntp, modes, psolve
      ay = ay - space
c draw text
      call gtx(ax,ay,chr(1:33))
      write (chr,94) t0, tend, dt
      ay = ay - space
c draw text
      call gtx(ax,ay,chr(1:38))
      write (chr,95) lts, its, nts
      ay = ay - space
c draw text
      call gtx(ax,ay,chr(1:33))
      write (chr,96) kmin, kmax
      ay = ay - space
c draw text
      call gtx(ax,ay,chr(1:24))
      write (chr,97) ntd, ntc
      ay = ay - space
c draw text
      call gtx(ax,ay,chr(1:22))
      write (chr,98) wmin, wmax, dw
      ay = ay - space
c draw text
      call gtx(ax,ay,chr(1:40))
c update workstation, perform
   40 call guwk(idwk,1)
c read code from input device, if present
      call readrc(irc)
      return
      end
      subroutine WPCORRV1(runid,indx,ntp,modes,psolve,t0,tend,dt,lts,its
     1,nts,kmin,kmax,ntd,ntc,nvf,wmin,wmax,dw,irc)
c this subroutine displays parameters for correlation program
c input arguments: all except irc
c irc = return code (0 = normal return)
      character*20 runid
      integer psolve
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c idwk = workstation identifier
c rx, ry = ndc coordinates of upper-right corner of workstation window
      character*48 chr
   91 format (' SPECTRUM ANALYSIS FOR 1D DATA')
   92 format (' RUNID= ',a20,' INDX=',i3)
   93 format (' NTP=',i6,' MODES=',i5,' PSOLVE=',i2)
   94 format (' T0=',f8.1,' TEND=',f8.1,' DT=',f8.5)
   95 format (' LTS=',i6,' ITS=',i6,' NTS=',i6)
   96 format (' KMIN=',i6,' KMAX=',i6)
   97 format (' NTD=',i6,' NTC=',i6,' NVF=',i6)
   98 format (' WMIN=',f8.4,' WMAX=',f8.4,' DW=',f8.4)
      data csize /0.024/
      irc = 0
      if (ry.lt..3) csize = .032
      chh = ry*csize
c space = vertical spacing between characters
      space = 1.5*chh
c initiate plot
c clear workstation, always
      call gclrwk(idwk,1)
c define and select identity transformation
      call idtran
c set default plotting parameters
c 1 = default font, 0 = string precision
      call dfplps (chh,1,0)
      ay = ry
      if (ay.lt.chh) return
      ax = csize*rx
      write (chr,91)
      ay = ay - space
c draw text
      call gtx(ax,ay,chr(1:30))
      write (chr,92) runid, indx
      ay = ay - space
c draw text
      call gtx(ax,ay,chr(1:37))
      write (chr,93) ntp, modes, psolve
      ay = ay - space
c draw text
      call gtx(ax,ay,chr(1:33))
      write (chr,94) t0, tend, dt
      ay = ay - space
c draw text
      call gtx(ax,ay,chr(1:38))
      write (chr,95) lts, its, nts
      ay = ay - space
c draw text
      call gtx(ax,ay,chr(1:33))
      write (chr,96) kmin, kmax
      ay = ay - space
c draw text
      call gtx(ax,ay,chr(1:24))
      write (chr,97) ntd, ntc, nvf
      ay = ay - space
c draw text
      call gtx(ax,ay,chr(1:33))
      write (chr,98) wmin, wmax, dw
      ay = ay - space
c draw text
      call gtx(ax,ay,chr(1:40))
c update workstation, perform
   40 call guwk(idwk,1)
c read code from input device, if present
      call readrc(irc)
      return
      end
      subroutine CORRE(pot,f,g,c,mixup,sct,inft,nt,ntc,nt2,nft,nfth,ntc2
     1)
c this program calculates the correlation function
c c(tau) = (1/(nt - tau))*sum{conjg(f(t))*g(t+tau)},
c where f(t) = g(t) = pot(t) - pot0, pot0 = sum{pot(t)}/nt, and
c the sums are over the nt values of t, and ntc values of tau are found
c all the real data in pot is stored before all the imaginary data
c complex fft tables must have been prepared for length nft,
c where nft must be greater than or equal to nt + ntc
c nfth = dimension of sine-cosine table, nfth >= nft/2
c nt2 = dimension of input array pot, nt2 >= 2*nt
c ntc2 = dimension of output array c, ntc2 >= 2*ntc
      complex f, g, sct, zero, zum1
      dimension pot(nt2), f(nft), g(nft), c(ntc2)
      dimension mixup(nft), sct(nfth)
      it2 = nt + 1
      zero = cmplx(0.,0.)
      zum1 = zero
      do 10 i = 1, nt
      g(i) = cmplx(pot(i),pot(i+nt))
      zum1 = zum1 + g(i)
   10 continue
      do 20 i = it2, nft
      g(i) = zero
   20 continue
      zum1 = zum1/real(nt)
      do 30 i = 1, nt
      g(i) = g(i) - zum1
   30 continue
      isign = 1
      call FFT1C(g,isign,mixup,sct,inft,nft,nfth)
      do 40 i = 1, nft
      g(i) = g(i)*conjg(f(i))
   40 continue
      isign = -1
      call FFT1C(g,isign,mixup,sct,inft,nft,nfth)
c     at2 = real(nt)/real(g(1))
      at2 = 1.
      do 50 i = 1, ntc
      at1 = at2/real(nt - i + 1)
      c(i) = at1*real(g(i))
   50 continue
      do 60 i = 1, ntc
      at1 = at2/real(nt - i + 1)
      c(i+ntc) = at1*aimag(g(i))
   60 continue
      return
      end
      subroutine SPECT(f,wm,p,t0,dt,nt,iw,nt2,iw2)
c this subroutine performs frequency analysis of complex time series,
c p(w) = |(1/nt)*sum {f(t)*exp(sqrt(-1)*wm*(t-t0))}|**2
c it is an sft (slow fourier transform), but you can pick your frequency
c on input, f contains the data to be analyzed, real and imaginary
c parts stored adjacent, and wm(w) contains the (positive) frequencies.
c on output, the first have of the array p contains result for positive
c frequencies, the second half the negative frequencies.
c t0 = starting time value
c dt = time step
c nt = number of input data points, 
c iw = number of (positive) frequencies
c nt2 = dimension of input array, nts >= 2*nt
c iw2 = dimension of power spectrum array iw2 >= 2*iw
      double precision at1,at2,at3,sum1,sum2,sum3,cwdt,swdt
      dimension f(nt2), wm(iw), p(iw2)
      anl = 1./real(nt)
      do 20 k = 1, iw
      k1 = k + iw
      at3 = wm(k)*t0
      at1 = dcos(at3)
      at2 = dsin(at3)
      at3 = wm(k)*dt
      cwdt = dcos(at3)
      swdt = dsin(at3)
      sum1 = 0.
      sum2 = 0.
      sum3 = 0.
      sum4 = 0.
      do 10 i = 1, nt
      i1 = i + nt
      sum1 = sum1 + f(i)*at1
      sum2 = sum2 + f(i1)*at1
      sum3 = sum3 + f(i1)*at2
      sum4 = sum4 + f(i)*at2
      at3 = at1*cwdt - at2*swdt
      at2 = at2*cwdt + at1*swdt
      at1 = at3
   10 continue
      at1 = anl*(sum1 - sum3)
      at2 = anl*(sum2 + sum4)
      p(k) = at1*at1 + at2*at2
      at1 = anl*(sum1 + sum3)
      at2 = anl*(sum2 - sum4)
      p(k1) = at1*at1 + at2*at2
   20 continue
      return
      end
      subroutine xspect(f,wm,p,t0,dt,nt2,iw,iw2)
c this subroutine performs frequency analysis of time series data,
c p(w) = |(1/nt)*sum {f(t)*exp(sqrt(-1)*wm*(t-t0))}|**2
c it is an sft (slow fourier transform), but you can pick your frequency
      double precision at1,at2,at3,sum1,sum2,sum3,wt0,wdt
      dimension f(nt2), wm(iw), p(iw2)
      nt = nt2/2
      anl = 1./real(nt)
      do 20 k = 1, iw
      k1 = k + iw
      wt0 = wm(k)*t0
      wdt = wm(k)*dt
      sum1 = 0.
      sum2 = 0.
      sum3 = 0.
      sum4 = 0.
      do 10 i = 1, nt
      i1 = i + nt
      at3 = wt0 + wdt*real(i - 1)
      at1 = dcos(at3)
      at2 = dsin(at3)
      sum1 = sum1 + f(i)*at1
      sum2 = sum2 + f(i1)*at1
      sum3 = sum3 + f(i1)*at2
      sum4 = sum4 + f(i)*at2
   10 continue
      at1 = anl*(sum1 - sum3)
      at2 = anl*(sum2 + sum4)
      p(k) = at1*at1 + at2*at2
      at1 = anl*(sum1 + sum3)
      at2 = anl*(sum2 - sum4)
      p(k1) = at1*at1 + at2*at2
   20 continue
      return
      end
      subroutine SAK1(vpk,nx,modesx,modesxd)
c this subroutine calculates amplitude of the wave number for each
c fourier mode, using the equations:
c vpk(kx) = kx*
c where kx = 2pi*j/nx, and j = fourier mode number
c nx = system length in x direction
c modesx = number of modes to store in x direction,
c where modesx <= nx/2+1
c modesxd = second dimension of output array vpott, modesxd >= modesx
      implicit none
      integer nx, modesx, modesxd
      real vpk
      dimension vpk(modesxd)
c local data
      integer nxh, jmx, j
      real dnx, dkx
      nxh = nx/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      jmx = min0(modesx,nxh)
      dnx = 6.28318530717959/real(nx)
c mode numbers 0 < kx < nx/2
      do 10 j = 2, jmx
      dkx = dnx*real(j - 1)
      vpk(j) = dkx
   10 continue
      vpk(1) = 0.0
      return
      end
      subroutine SQ2POT1(pott,potr,nx,modesx,modesxd)
c this subroutine calculates the square of the amplitude of
c complex scalar array pott and stores it in real array potr
c nx = system length in x direction
c modesx = number of modes to store in x direction,
c where modesx <= nx/2+1
c modesxd = second dimension of output array pott, modesxd >= modesx
      implicit none
      integer nx, modesx,  modesxd
      real potr
      complex pott
      dimension pott(modesxd), potr(modesxd)
c local data
      integer nxh, jmx, j
      complex zt1
      nxh = nx/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      jmx = min0(modesx,nxh)
c calculate the square of the amplitude
c mode numbers 0 < kx < nx/2
      do 10 j = 2, jmx
      zt1 = pott(j)
      potr(j) = zt1*conjg(zt1)
   10 continue
      potr(1) = 0.0
      return
      end
      subroutine WEL1(pott,potb,affp,we,nx,modesx,modesxd)
c this subroutine calculates longitudinal electric field energy from
c potential for each fourier mode
c the energy is calculated using the equations:
c potb(kx) = fx(kx)*conjg(fx(kx)/affp
c fx(kx,ky) = -sqrt(-1)*kx*pott(kx)
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c except for fx(kx=pi) = fx(kx=0) = 0.
c affp = normalization constant = nx/np, where np=number of particles
c we = electric field energy
c nx = system length in x direction
c modesx = number of modes to store in x direction,
c where modesx <= nx/2+1
c modesxd = second dimension of output array pott, modesxd >= modesx
      implicit none
      integer nx, modesx, modesxd
      real potb, affp, we
      complex pott
      dimension pott(modesxd), potb(modesxd)
      integer nxh, jmx, j
      real dnx, dkx, anorm
      complex zt1
      double precision wp
      nxh = nx/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      jmx = min0(modesx,nxh)
      dnx = 6.28318530717959/real(nx)
      anorm = real(nx)/affp
      wp = 0.0d0
c calculate the energy per mode
c mode numbers 0 < kx < nx/2
      do 10 j = 2, jmx
      dkx = dnx*real(j - 1)
      zt1 = dkx*cmplx(aimag(pott(j)),-real(pott(j)))
      we = anorm*(zt1*conjg(zt1))
      potb(j) = we
      wp = wp + we
   10 continue
      potb(1) = 0.0
      we = wp
      return
      end
      subroutine SQ2VPOT1(vpott,vpotr,nx,modesx,modesxd)
c this subroutine calculates the square of the amplitude of
c complex vector array vpott and stores it in real array vpotr
c nx = system length in x direction
c modesx = number of modes to store in x direction,
c where modesx <= nx/2+1
c modesxd = second dimension of output array vpott, modesxd >= modesx
      implicit none
      integer nx, modesx, modesxd
      real vpotr
      complex vpott
      dimension vpott(2,modesxd), vpotr(modesxd)
      integer nxh, jmx, j
      complex zt1, zt2
      nxh = nx/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      jmx = min0(modesx,nxh)
c calculate the square of the amplitude
c mode numbers 0 < kx < nx/2
      do 10 j = 2, jmx
      zt1 = vpott(1,j)
      zt2 = vpott(2,j)
      vpotr(j) = zt1*conjg(zt1) + zt2*conjg(zt2) 
   10 continue
      vpotr(1) = 0.0
      return
      end
      subroutine WBR1(vpott,vpotb,affp,ci,wm,nx,modesx,modesxd)
c this subroutine calculates magnetic field energy from vector potential
c for each fourier mode
c the energy is calculated using the equations:
c vpotb(kx) = (gy(kx)*conjg(gy(kx)+gz(kx)*conjg(gz(kx)))/(affp*ci**2)
c gy(kx) = -sqrt(-1)*kx*vpott(2,kx)
c gz(kxy) = sqrt(-1)*(kx*vpott(1,kx))
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c except for gy(kx=pi) = gz(kx=pi) = 0, and gy(kx=0) = gz(kx=0) = 0.
c affp = normalization constant = nx/np, where np=number of particles
c ci = reciprical of velocity of light
c wm = total magnetic field energy
c nx = system length in x direction
c modesx = number of modes to store in x direction,
c where modesx <= nx/2+1
c modesxd = second dimension of output array vpott, modesxd >= modesx
      implicit none
      integer nx, modesx, modesxd
      real vpotb, affp, ci, wm
      complex vpott
      dimension vpott(2,modesxd), vpotb(modesxd)
      integer nxh, jmx, j
      real dnx, dkx, anorm
      complex zt1, zt2, zt3
      double precision wp
      nxh = nx/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      jmx = min0(modesx,nxh)
      dnx = 6.28318530717959/real(nx)
      anorm = real(nx)/(affp*ci*ci)
      wp = 0.0d0
c calculate the energy per mode
c mode numbers 0 < kx < nx/2
      do 10 j = 2, jmx
      dkx = dnx*real(j - 1)
      zt1 = cmplx(-aimag(vpott(2,j)),real(vpott(2,j)))
      zt2 = cmplx(-aimag(vpott(1,j)),real(vpott(1,j)))
      zt3 = dkx*zt2
      zt2 = -dkx*zt1
      wm = anorm*(zt2*conjg(zt2) + zt3*conjg(zt3))
      vpotb(j) = wm
      wp = wp + wm
   10 continue
      vpotb(1) = 0.0
      wm = wp
      return
      end
      subroutine WER1(vpotp,vpote,affp,wf,nx,modesx,modesxd)
c this subroutine calculates transverse electric field energy from
c difference of vector potential for each fourier mode
c the energy is calculated using the equations:
c vpote(kx) = gy(kx)*conjg(gy(kx)+gz(kx)*conjg(gz(kx))), where
c gy(kx) = vpotp(1,kx), and gz(kx) = vpotp(2,kx)
c where kx = 2pi*j/nx, and j = fourier mode numbers,
c except for gy(kx=pi) = gz(kx=pi) = 0, and gy(kx=0) = gz(kx=0) = 0.
c affp = normalization constant = nx/np, where np=number of particles
c wf = total transverse electric field energy
c nx = system length in x direction
c modesx = number of modes to store in x direction,
c where modesx <= nx/2+1
c modesxd = second dimension of output array vpotp, modesxd >= modesx
      implicit none
      integer nx, modesx, modesxd
      real vpote, affp, wf
      complex vpotp
      dimension vpotp(2,modesxd), vpote(modesxd)
      integer nxh, jmx, j
      real anorm
      complex zt2, zt3
      double precision wp
      nxh = nx/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      jmx = min0(modesx,nxh)
      anorm = real(nx)/affp
      wp = 0.0d0
c calculate the energy per mode
c mode numbers 0 < kx < nx/2
      do 10 j = 2, jmx
      zt2 = vpotp(1,j)
      zt3 = vpotp(2,j)
      wf = anorm*(zt2*conjg(zt2) + zt3*conjg(zt3))
      vpote(j) = wf
      wp = wp + wf
   10 continue
      vpote(1) = 0.0
      wf = wp
      return
      end
