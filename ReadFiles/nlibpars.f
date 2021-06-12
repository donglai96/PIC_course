c parsing library
c copyright 1991, regents of the university of california
c update: july 12, 2011
      subroutine vinput(iunit,junit,kunit,icmpl,ircopy)
c this subroutine interactively modifies and updates input variables.
c current values of the input variables, in namelist format, can be read
c from a file connected to unit=junit.  a variable description file
c (read from a file connected to unit=iunit) is required.  If junit=0,
c then default values for the input variables are read from the
c variable description file rather than from the namelist input file.
c after modification, updated values in namelist format are written back
c to unit=junit, if junit > 0.
c an include file consisting of parameter statements (for variables
c flagged in the variable description file) and a namelist declaration
c statement (for variables not flagged) is written to a file connected
c to unit=kunit, if kunit > 0.  this file is also not written if
c iwrt = 0 and icmpl = 0.
c input arguments: iunit, junit, kunit
c iunit/junit/kunit = file unit numbers for input data, variable
c description file and include file.  the files corresponding to these
c unit numbers are assumed to be already open.
c two additional flags are also returned:
c icmpl = (0,1) = (no,yes) a variable needed for dimensions has changed.
c ircopy = (0,1) = (no,yes) request to save current input in new file
c written by viktor k. decyk, ucla
c lpmax/lgmax/nvmax = maximum number of pages/groups/variables
c ncmax = maximum number of character variables
      parameter (lpmax=20,lgmax=100,nvmax=200,ncmax=10,itwo=2)
      character*8 code(nvmax)
      character*20 cvalue(ncmax)
      character*72 helpv(nvmax), group(lgmax), page(lpmax)
      dimension icode(nvmax), ig(lgmax), ip(lpmax)
      dimension value(nvmax), range(itwo,nvmax)
      save code,cvalue,helpv,group,page,icode,ig,ip,value,range,iwrt
c iwrt = (0,1) = write output parameter file (only if icmpl=1, always)
      data iwrt /1/
c create variable list and headings from variable list file
      call PARSVF(page,ip,group,ig,code,helpv,icode,value,range,cvalue,l
     1pmax,lgmax,nvmax,ncmax,itwo,npage,ngroup,nvars,iunit)
      rewind iunit
      if (nvars.eq.0) then
         write (6,'(1x,25a)') 'error: no variables found'
         return
      endif
c read namelist input file
      call RDNLST(code,icode,value,range,cvalue,nvars,ncmax,itwo,junit)
      rewind junit
c modify variables
      call MENUVFD(page,ip,group,ig,code,helpv,icode,value,range,cvalue,
     1npage,ngroup,nvars,ncmax,itwo,icmpl,ircopy,iunit,irc)
c quit
      if (irc.eq.1) return
c update namelist input file
      call WRNLST(code,icode,value,cvalue,nvars,ncmax,junit)
c write the include file
      if ((iwrt.eq.1).or.(icmpl.eq.1)) call WRITVF(code,icode,value,cval
     1ue,nvars,ncmax,kunit)
      return
      end
      subroutine PARSVF(page,ip,group,ig,code,helpv,icode,value,range,cv
     1alue,lpmax,lgmax,nvmax,ncmax,itwo,npage,ngroup,nvars,iunit)
c this subroutine reads a variable list file and creates variable
c lists and headings for a menu driven pre-processor for fortran
c each line of variable list file must be preceeded by a special symbol
c @ = new page header, ! = new group header. % = integer variable,
c & = floating point variable, $ = character variable,
c * = special variable, / = range specification, ? = help for variable
c # = item needed for dimensions,
c lines without special symbols are considered comments and are ignored
c input arguments: lpmax, lgmax, nvmax, ncmax, itwo, iunit
c lpmax/lgmax/nvmax = maximum number of pages/groups/variables
c ncmax = maximum number of character variables
c itwo = first dimension of range array
c iunit = file unit number for variable description file, assumed to
c be already open
c on exit the following are defined:
c page = page heading character array
c ip = array containing the number of groups in each page
c group = group heading character array
c ig = array containing the number of variables in each group
c code = variable symbol character array
c helpv = character array containing one line definition of variable
c icode = code array for variable types,1=integer,2=floating point,
c    3=character,4=special, plus 8 if variable needed for dimensions
c value = array containing default value of numeric variable
c range = array containing minimum/maximum allowed range of value
c cvalue = array containing value of character or special variable
c npage/ngroup/nvars = number of pages/groups/variables
      character*8 code(nvmax)
      character*20 cvalue(ncmax)
      character*72 helpv(nvmax), group(lgmax), page(lpmax)
      dimension icode(nvmax), ig(lgmax), ip(lpmax)
      dimension value(nvmax), range(itwo,nvmax)
      character*80 c
      character*20 cval
      character*1 c1,c2,c3,c4,c5,c6,c7,c8,c9,v
      save c1,c2,c3,c4,c5,c6,c7,c8,c9
   91 format (a80)
c special symbols
      data c1,c2,c3,c4,c5,c6,c7,c8 /'@','!','%','&','$','*','/','?'/
      data c9 /'#'/
      data i,j,k,nc,newp,nhelp /0,0,0,0,1,-1/
c initialize
      lenc = len(c)
      page(1) = 'menu 1'
      ip1 = 1
      ig1 = 1
      if (iunit.le.0) then
         write (6,'(1x,a21,i8)') 'invalid input unit = ',iunit
         go to 20
      endif
c main loop
   10 read (iunit,91,end=20) c
      v = c(1:1)
c check if new page
      if (v.eq.c1) then
         if (i.gt.0) ip(i) = j - ip1 + 1
         i = i + 1
         if (i.gt.lpmax) then
            write (6,'(1x,a20)') 'page memory exceeded'
            i = i - 1
            go to 20
         endif
         ip1 = j + 1
         page(i) = c(2:lenc)
         newp = 1
         nhelp = -1
         go to 10
c check if new group
      elseif (v.eq.c2) then
         if (j.gt.0) ig(j) = k - ig1 + 1
         j = j + 1
         if (j.gt.lgmax) then
            write (6,'(1x,a21)') 'group memory exceeded'
            j = j - 1
            go to 20
         endif
         ig1 = k + 1
         group(j) = c(2:lenc)
         newp = 0
         nhelp = -1
         go to 10
c check if integer variable
      elseif (v.eq.c3) then
         it = 1
c check if floating point variable
      elseif (v.eq.c4) then
         it = 2
c check if character variable
      elseif (v.eq.c5) then
         it = 3
c check if special variable
      elseif (v.eq.c6) then
         it = 4
c check if range specified
      elseif (v.eq.c7) then
         call rpars (c(2:lenc),icode(k),range(1,k),itwo)
         if ((icode(k) - (icode(k)/8)*8).lt.3) then
            if (value(k).lt.range(1,k)) value(k) = range(1,k)
            if (value(k).gt.range(2,k)) value(k) = range(2,k)
         endif
         go to 10
c check if variable help
      elseif (v.eq.c8) then
         nhelp = nhelp + 1
         if (nhelp.eq.1)  helpv(k) = c(2:lenc)
         go to 10
      else
         go to 10
      endif
c calculate codes for variables
      if (i.eq.0) i = 1
c force new group
      if (newp.ne.0) then
         if (j.gt.0) ig(j) = k - ig1 + 1
         j = j + 1
         if (j.gt.lgmax) then
            write (6,'(1x,a21)') 'group memory exceeded'
            j = j - 1
            go to 20
         endif
         ig1 = k + 1
         group(j) = ' '
      endif
      k = k + 1
      if (k.gt.nvmax) then
         write (6,'(1x,a24)') 'variable memory exceeded'
         k = k - 1
         go to 20
      endif
      helpv(k) = 'definition unavailable'
c check if item needed for dimensions
      is = 0
      if (c(2:2).eq.c9) is = 1
      icode(k) = 8*is + it
      is = is + 2
      call vpars(c(is:lenc),code(k),ic,value(k),id,cval,irc)
      call rpars (c7,icode(k),range(1,k),itwo)
c integer or floating point variable
      if (it.lt.3) then
         if (value(k).lt.range(1,k)) value(k) = range(1,k)
         if (value(k).gt.range(2,k)) value(k) = range(2,k)
c character or special character
      else
         nc = nc + 1
         if (nc.gt.ncmax) then
            write (6,'(1x,a25)') 'character memory exceeded'
            nc = nc - 1
            go to 20
         endif
         value(k) = nc
         cvalue(nc) = cval
      endif
      newp = 0
      nhelp = 0
      if (irc.ne.1) go to 10
c invalid variable name
      k = k - 1
      if (it.ge.3) nc = nc - 1
      go to 10
c finish
   20 npage = i
      ngroup = j
      nvars = k
      ip(i) = j - ip1 + 1
      ig(j) = k - ig1 + 1
      return
      end
      subroutine RDNLST(code,icode,value,range,cvalue,nvars,ncmax,itwo,i
     1unit)
c this subroutine reads input data for variables in input data file.
c the value read replaces the default value if the symbol is found and
c its value is a valid number within range.
c input arguments: code, range, nvars, ncmax, itwo, iunit
c code = variable symbol character array
c icode = code array for variable types,1=integer,2=floating point,
c    3=character,4=special, plus 8 if variable needed for dimensions
c range = array containing minimum/maximum allowed range of value
c nvars = number of variables
c ncmax = maximum number of character variables
c itwo = first dimension of range array
c iunit = file unit number for input data, assumed to be open
c output variables written:
c value = array containing default value of numeric variable
c cvalue = array containing value of character or special variable
      character*8 code(nvars)
      character*20 cvalue(ncmax)
      character*80 c
      character*22 cval
      character*8 codev
      dimension icode(nvars), value(nvars), range(itwo,nvars)
   91 format (a80)
      if (iunit.le.0) return
      im = len(c)
c read next input line
   10 read (iunit,91,end=40) c
c find comma
      i = 0
      ns = 1
      iflg = 0
   20 i = i + 1
c end of line found
      if (i.gt.im) then
         iflg = 1
c try again
      elseif (c(i:i).ne.',') then
         go to 20
      endif
c comma or end of line found
      nv = i - 1
      if ((nv-ns).lt.0) go to 30
      call vpars(c(ns:nv),codev,ic,val,id,cval,irc)
      if (irc.eq.1) go to 30
c find symbol
      call findvn(code,codev,nvars,nvar)
      if (nvar.eq.0) then
c error message
         write (6,'(1x,a24,a8)') 'unknown input variable: ', codev
         go to 30
      endif
c check variable type
      itype = icode(nvar) - (icode(nvar)/8)*8
c integer or floating point variable
      if (itype.lt.3) then
c valid number, check if in range
         if ((val.ge.range(1,nvar)).and.(val.le.range(2,nvar))) then
            value(nvar) = val
         endif
c character or special variable
      else
         if (irc.ne.2) then
            nc = value(nvar)
            is = 1
            if (cval(1:1).eq.'''') is = 2
            cvalue(nc) = cval(is:is+20)
         endif
      endif
c find next variable in current line
   30 if (iflg.eq.0) then
         ns = nv + 2
         go to 20
c read next line
      else
         go to 10
      endif
c end of file
   40 return
      end
      subroutine WRNLST(code,icode,value,cvalue,nvars,ncmax,iunit)
c this subroutine creates a file with two namelists. the first namelist,
c called input, is for variables which are not flagged with icdoe > 8.
c the second namelist, called parms, is for those which are.
c input arguments: all
c code = variable symbol character array
c icode = code array for variable types,1=integer,2=floating point,
c    3=character,4=special, plus 8 if variable needed for dimensions
c value = array containing default value of numeric variable
c cvalue = array containing value of character or special variable
c nvars = number of variables
c ncmax = maximum number of character variables
c iunit = file unit number for output data, assumed to be open
      character*8 code(nvars)
      character*20 cvalue(ncmax)
      dimension icode(nvars), value(nvars)
      save kstrt,iquote,lvm
c iquote = (0,1) = (no,yes) include quote marks in character strings
c lvm = maximum length of variable display line
      data kstrt,iquote,lvm /1,1,59/
      if (iunit.le.0) return
c write header for first namelist
      write (iunit,'(1x,a6)') '&input'
      isel = 1
c write body of first namelist
      call wrvarf(code,icode,value,cvalue,nvars,ncmax,iunit,kstrt,nvars,
     1isel,iquote,lvm)
c write footer for first namelist
      write (iunit,'(1x,a4)') '&end'
c     write (iunit,'(1x,a1)') '/'
c write header for second namelist
      write (iunit,'(1x,a6)') '&parms'
      isel = 0
c write body of second namelist
      call wrvarf(code,icode,value,cvalue,nvars,ncmax,iunit,kstrt,nvars,
     1isel,iquote,lvm)
c write footer for second namelist
      write (iunit,'(1x,a4)') '&end'
      return
      end
      subroutine wrvarf(code,icode,value,cvalue,nvars,ncmax,iunit,kstrt,
     1knum,isel,iquote,lvm)
c this subroutine writes selected variables in a namelist format, that
c is, in the form variable = value.
c input arguments: all
c code = variable symbol character array
c icode = code array for variable types,1=integer,2=floating point,
c    3=character,4=special, plus 8 if variable needed for dimensions
c value = array containing default value of numeric variable
c cvalue = array containing value of character or special variable
c nvars = number of variables
c ncmax = maximum number of character variables
c iunit = file unit number for output data, assumed to be open
c kstrt = index of first variable to be searched
c knum = number of variables to be searched
c isel = (-1,0,1) = select (all,icode>8,icode<8) variables,
c iquote = (0,1) = (no,yes) include quote marks in character strings
c lvm = maximum length of variable display line
      character*8 code(nvars)
      character*20 cvalue(ncmax)
      character*80 chr
      character*32 chrv
      character*1 c1,c2
      dimension icode(nvars), value(nvars)
   91 format (a)
      save c1,c2,komma
      data c1,c2 /',',' '/
c komma = (0,1) = delimit variables with (blank,comma)
      data komma /0/
      if (iunit.le.0) return
      nc = 1
      lvs = 0
      chr = ' '
      ks = kstrt - 1
c main loop
      do 10 k = 1, knum
      k1 = k + ks
      itc = icode(k1)/8
c skip unselected items
      if (itc.eq.isel) go to 10
c if character find location in cvalue array
      if ((icode(k1) - 8*itc).ge.3) nc = value(k1)
c get variable name and value string
      call frmtv(code(k1),icode(k1),value(k1),cvalue(nc),iquote,chrv,lv)
      if (lv.gt.lvm) lv = lvm
      lvt = lvs + lv
c write out line to avoid overflow
      if (lvt.gt.lvm) then
c replace comma with blank
         if (komma.eq.0) chr(lvs:lvs) = c2
         write (iunit,91) chr(1:lvs)
         lvs = 0
         chr = ' '
         lvt = lv
      endif
      ls1 = lvs + 1
c add variable name and value string to output string
      chr(ls1:ls1+lv) = chrv(1:lv)
      lvs = lvt + 1
c add comma
      chr(lvs:lvs) = c1
   10 continue
c flush last line
      if (lvs.gt.0) then
c replace comma with blank
         if (komma.eq.0) chr(lvs:lvs) = c2
         write (iunit,91) chr(1:lvs)
      endif
      return
      end
      subroutine WRITVF(code,icode,value,cvalue,nvars,ncmax,iunit)
c this subroutine creates a file of parameter statements for variables
c in a variable list file which are flagged with icode > 8, and a
c namelist declaration for the remaining variables.
c input arguments: all
c code = variable symbol character array
c icode = code array for variable types,1=integer,2=floating point,
c    3=character,4=special, plus 8 if variable needed for dimensions
c value = array containing default value of numeric variable
c cvalue = array containing value of character or special variable
c nvars = number of variables
c ncmax = maximum number of character variables
c iunit = file unit number for output data, assumed to be open
      character*8 code(nvars)
      character*20 cvalue(ncmax)
      character*80 chr
      character*32 chrv
      character*10 parml
      character*1 c1,c2,c3
      dimension icode(nvars), value(nvars)
   91 format (6x,'namelist /input/',a50)
   92 format (5x,i1,a66)
   93 format (6x,'namelist /parms/',a50)
      save iquote,iparm,c1,c2,c3,parml
c iquote = (0,1) = (no,yes) include quote marks in character strings
c iparm = (0,1) = write parameter variables as (namelist,parameters)
      data iquote,iparm /1,0/
      data c1,c2,c3 /')',',',' '/
      data parml /'parameter('/
      if (iunit.le.0) return
c check whether to write parameter variables as namelist or parameters
      if (iparm.eq.0) go to 20
c write parameter variables as parameters
      nc = 1
      lvs = 0
      lvm = 55
      chr = ' '
c main loop for parameter variables
      do 10 k = 1, nvars
      itc = icode(k)/8
c skip variables which are not flagged as parameters
      if (itc.eq.0) go to 10
c if character find location in cvalue array
      if ((icode(k) - 8*itc).ge.3) nc = value(k)
c get variable name and value string
      call frmtv(code(k),icode(k),value(k),cvalue(nc),iquote,chrv,lv)
      if (lv.gt.lvm) lv = lvm
      lvt = lvs + lv
c write out line to avoid overflow
      if (lvt.gt.lvm) then
c replace comma with parenthesis
         chr(lvs:lvs) = c1
         write (iunit,'(6x,a10,a56)') parml, chr
         lvs = 0
         chr = ' '
         lvt = lv
      endif
      ls1 = lvs + 1
c add variable name and value string to output string
      chr(ls1:ls1+lv) = chrv(1:lv)
      lvs = lvt + 1
c add comma
      chr(lvs:lvs) = c2
   10 continue
c flush last line
      if (lvs.gt.0) then
c replace comma with parenthesis
         chr(lvs:lvs) = c1
         write (iunit,'(6x,a10,a56)') parml, chr
      endif
c write namelist declaration
   20 nc = 1
      lvs = 0
      lvm = 49
      nvd = 0
      chr = ' '
c main loop for namelist variables
      do 30 k = 1, nvars
      itc = icode(k)/8
c skip variables which are not namelist variables
      if (itc.ge.1) go to 30
c if character find location in cvalue array
      if ((icode(k) - 8*itc).ge.3) nc = value(k)
c get variable name string
      call frmtn(code(k),chrv,lv)
      if (lv.gt.lvm) lv = lvm
      lvt = lvs + lv
c write out line to avoid overflow
      if (lvt.gt.lvm) then
         if (nvd.eq.0) write (iunit,91) chr
         if (nvd.gt.0) write (iunit,92) nvd, chr
         lvs = 0
         lvt = lv
         lvm = 65
         nvd = nvd + 1
         chr = ' '
      endif
      ls1 = lvs + 1
c add variable name string to output string
      chr(ls1:ls1+lv) = chrv(1:lv)
      lvs = lvt + 1
c add comma
      chr(lvs:lvs) = c2
   30 continue
c flush last line
      if (lvs.gt.0) then
c replace comma with blank
         chr(lvs:lvs) = c3
         if (nvd.eq.0) write (iunit,91) chr
         if (nvd.gt.0) write (iunit,92) nvd, chr
      endif
c check whether to write parameter variables as namelist or parameters
      if (iparm.eq.1) return
c write parameter variables as namelist
      nc = 1
      lvs = 0
      lvm = 49
      nvd = 0
      chr = ' '
c main loop for namelist variables
      do 40 k = 1, nvars
      itc = icode(k)/8
c skip variables which are not namelist variables
      if (itc.eq.0) go to 40
c if character find location in cvalue array
      if ((icode(k) - 8*itc).ge.3) nc = value(k)
c get variable name string
      call frmtn(code(k),chrv,lv)
      if (lv.gt.lvm) lv = lvm
      lvt = lvs + lv
c write out line to avoid overflow
      if (lvt.gt.lvm) then
         if (nvd.eq.0) write (iunit,93) chr
         if (nvd.gt.0) write (iunit,92) nvd, chr
         lvs = 0
         lvt = lv
         lvm = 65
         nvd = nvd + 1
         chr = ' '
      endif
      ls1 = lvs + 1
c add variable name string to output string
      chr(ls1:ls1+lv) = chrv(1:lv)
      lvs = lvt + 1
c add comma
      chr(lvs:lvs) = c2
   40 continue
c flush last line
      if (lvs.gt.0) then
c replace comma with blank
         chr(lvs:lvs) = c3
         if (nvd.eq.0) write (iunit,93) chr
         if (nvd.gt.0) write (iunit,92) nvd, chr
      endif
      return
      end
      subroutine MENUVFD(page,ip,group,ig,code,helpv,icode,value,range,c
     1value,npage,ngroup,nvars,ncmax,itwo,icmpl,ircopy,iunit,irc)
c this subroutine displays variables and their values read from a
c variable list file by the subroutine parsvf, and allows users to
c modify them.  the subroutine also checks if the the values are in
c the allowed range.
c input arguments: all except icmpl, ircopy
c page = page heading character array
c ip = array containing the number of groups in each page
c group = group heading character array
c ig = array containing the number of variables in each group
c code = variable symbol character array
c helpv = character array containing one line definition of variable
c icode = code array for variable types,1=integer,2=floating point,
c    3=character,4=special, plus 8 if variable needed for dimensions
c value = array containing default value of numeric variable
c range = array containing minimum/maximum allowed range of value
c cvalue = array containing value of character or special variable
c npage/ngroup/nvars = number of pages/groups/variables
c ncmax = maximum number of character variables
c itwo = first dimension of range array
c icmpl = (0,1) = (no,yes) variable needed for dimensions has changed.
c ircopy = (0,1) = (no,yes) save current input in new file
c iunit = file unit number for variable description file, assumed open
c irc = return code (0 = normal return)
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c idwk = workstation identifier
c rx, ry = ndc coordinates of upper-right corner of workstation window
c idstr = string device number, 0 if no string device available
      character*8 code(nvars)
      character*20 cvalue(ncmax)
      character*72 helpv(nvars), group(ngroup), page(npage)
      character*60 chr, prompt1
      character*32 input
      character*20 cval
      character*8 codev, cmenu, crun, csave, cquit, chelp, cstyle
      character*1 v
      dimension icode(nvars), ig(ngroup), ip(npage)
      dimension value(nvars), range(itwo,nvars)
      save iquote,isel,istyle,lvm
      save cmenu,crun,csave,cquit,chelp,cstyle
   91 format (a32)
   92 format (i9,' out of range:',i9,'<',a8,'<',i9)
   93 format (e12.5,' out of range:',e12.5,'<',a8,'<',e12.5)
   94 format (' value= ',i9,', range= ',i9,'<value<',i9)
   95 format (' value= ',e12.5,', range= ',e12.5,'<value<',e12.5)
   96 format (' current value = ',a20)
c iquote = (0,1) = (no,yes) include quote marks in character strings
c lvm = maximum length of variable display line
      data iquote,isel,istyle,lvm /0,-1,2,59/
c the following commands should be in upper case:
      data cmenu,crun,csave /'MENU    ','RUN     ','SAVE    '/
      data cquit,chelp,cstyle /'QUIT    ','HELP    ','STYLE   '/
      data prompt1 /' enter menu number, or type run, save, help, style,
     1 or quit '/
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
      nc = 1
      icmpl = 0
      ircopy = 0
c display main menu
      itc = 1
c display new page
      if (npage.eq.1) go to 40
      is = 0
      in = 0
      js0 = 0
      jn0 = 0
c display page headers of main menu
c clear workstation, always
   10 call gclrwk(idwk,1)
      ay = ry
      irt = 1
      iz = ichar('0')
      do 20 i = 1, npage
      ay = ay - space
c draw text
      if (ay.ge.0.) call gtx(ax,ay,' '//char(i+iz)//': '//page(i))
   20 continue
c display prompt
      ay = ay - space
c draw text
      if (ay.ge.0.) call gtx(ax,ay,prompt1)
c read input
c update workstation, perform
   30 call guwk(idwk,1)
      input = ' '
c request string
      call grqst(idwk,idstr,istat,lens,input)
c     read (5,91,end=35) input
      v = input(1:1)
c check if new page is requested
      itc = ichar(v) - ichar('0')
c re-display page headers of main menu
      if (itc.eq.0) go to 10
c if valid page number entered, display new page
      if ((itc.ge.1).and.(itc.le.npage)) go to 40
c check if system command
      im = len(input)
      if (v.eq.'!') then
         call opsysc(input(2:im))
c abort
      elseif ((v.eq.'q').or.(v.eq.'Q')) then
         if (input(2:2).eq.' ') then
           irc = 1
           return
         endif
c check if special help command
      elseif (v.eq.'?') then
c clear workstation, always
         call gclrwk(idwk,1)
         ay = ry
c general help
         call helper(ax,ay,space)
c read more input
         go to 30
      endif
c parse input string
      call vpars(input,codev,ic,val,id,cval,ird)
c re-display current page
      if (ird.eq.1) go to 110
c check if known command
c main menu request
      if (codev.eq.cmenu) then
         go to 10
c normal exit
      elseif (codev.eq.crun) then
         return
c set save request flag
      elseif (codev.eq.csave) then
         ircopy = 1
c abort
      elseif (codev.eq.cquit) then
         irc = 1
         return
c general help
      elseif (codev.eq.chelp) then
c clear workstation, always
         call gclrwk(idwk,1)
         ay = ry
c general help
         call helper(ax,ay,space)
c read more input
         go to 30
c set style
      elseif (codev.eq.cstyle) then
         ival = val
c valid response
         if ((ird.eq.0).and.(ival.ge.0).and.(ival.le.2)) then
c do nothing if style does not change
            if (ival.ne.istyle) then
c reset style variable
               istyle = ival
c re-display current page
               go to 110
            endif
c invalid response
         else
            ay = ay - space
            if (ay.ge.0.) then
c draw text
               call gtx(ax,ay,' usage: style = (0,1,2) for (terse,normal
     1,verbose) displays')
            endif
         endif
c write prompt
         go to 100
      endif
c check if known variable
      call findvn(code,codev,nvars,nvar)
c re-display current page
      if (nvar.eq.0) go to 110
c examine type
      itype = icode(nvar) - (icode(nvar)/8)*8
c integer or floating point variable
      if (itype.lt.3) then
c help request
         if (ird.gt.1) go to 90
c valid number, check if in range
         if ((val.ge.range(1,nvar)).and.(val.le.range(2,nvar))) then
c value in range
            value(nvar) = val
c set flag to indicate parameter has been modified
            if ((icode(nvar)/8).ge.1) icmpl = 1
c re-display current page
            go to 110
         else
c value out of range
c floating point type
            if (itype.eq.2) then
            write (chr,93) val, range(1,nvar), code(nvar), range(2,nvar)
c integer type
            else
               ival = val
               ir1 = range(1,nvar)
               ir2 = range(2,nvar)
               write (chr,92) ival, ir1, code(nvar), ir2
            endif
            ay = ay - space
c draw text
            if (ay.ge.0.) call gtx(ax,ay,chr)
         endif
c character or special variable
      else
c help request
         if (ird.eq.2) go to 90
         nc = value(nvar)
         cvalue(nc) = cval
c set flag to indicate parameter has been modified
         if ((icode(nvar)/8).ge.1) icmpl = 1
c re-display current page
         go to 110
      endif
c write prompt
      go to 100
c new page
   40 ipage = itc
      is = 0
      in = 0
      js0 = 0
      jn0 = 0
c find offsets
      do 60 i = 1, ipage
      is = is + in
      in = ip(i)
      if (i.lt.ipage) then
         do 50 j = 1, in
         js0 = js0 + jn0
         jn0 = ig(j+is)
   50    continue
      endif
   60 continue
c display page
c clear workstation, always
   70 call gclrwk(idwk,1)
      ay = ry
      irt = 2
      js = js0
      jn = jn0
c loop over groups
      do 80 j = 1, in
      j1 = j + is
      if (istyle.ge.1) then
         ay = ay - space
c draw text
         if (ay.ge.0.) call gtx(ax,ay,group(j1))
      endif
      js = js + jn
      jn = ig(j1)
c display variables and values
      ks = js + 1
      call dsvarf(code,icode,value,cvalue,nvars,ncmax,ks,jn,isel,iquote,
     1lvm,ax,ay,space)
   80 continue
c write prompt
      go to 100
c variable help
   90 if (istyle.lt.2) then
         ay = ay - space
c draw text
         if (ay.ge.0.) call gtx(ax,ay,helpv(nvar))
c write prompt
         go to 100
      elseif (istyle.ne.2) then
c write prompt
         go to 100
c extended variable help
      else
c clear workstation, always
         call gclrwk(idwk,1)
         ay = ry
         call helpvf(codev,ax,ay,space,iunit)
         if (iunit.gt.0) rewind iunit
      endif
c integer type
      if (itype.eq.1) then
         ival = value(nvar)
         ir1 = range(1,nvar)
         ir2 = range(2,nvar)
         write (chr,94) ival, ir1, ir2
c floating point type
      elseif (itype.eq.2) then
         write (chr,95) value(nvar), range(1,nvar), range(2,nvar)
c character or special type
      elseif (itype.ge.3) then
         nc = value(nvar)
         write (chr,96) cvalue(nc)
      endif
      ay = ay - space
c draw text
      if (ay.ge.0.) call gtx(ax,ay,chr)
c write prompt
  100 ay = ay - space
      if (ay.ge.0.) then
c draw text
         if (irt.eq.1) then
            call gtx(ax,ay,prompt1)
         elseif (irt.eq.2) then
            call gtx(ax,ay,' enter variable=value, or type menu, run, sa
     1ve, or quit')
         endif
      endif
c read more input
      go to 30
c re-display current page
  110 go to (10,70), irt
      return
      end
      subroutine vpars(input,code,ic,value,id,cval,irc)
c this subroutine parses input and looks for the form: variable = value.
c the variable name found is returned in code, and its length is in ic.
c blanks are ignored, and lower case is translated to upper case.
c if no variable name is found or if the name contains characters other
c than letters or numbers, then an error code of 1 is returned.
c if a valid variable is found, a value is searched for.  if found, the
c characters are returned in cval, the numerical value is returned in
c value, and the number of numerical digits is returned in id.
c if no value is found, an error code of 2 is returned.
c if the value has no digits, then an error code of 3 is returned.
c input argument: input
c input = input character string to be parsed
c code = character string containing variable name
c ic = length of variable name
c value = numeric value of variable
c id = number of digits of numeric value of variable
c cval = string containing numeric value of variable
c irc = (0,1,2,3) = (normal,illegal,valueless,empty) variable
      character*(*) input
      character*(*) code
      character*(*) cval
      save ltou
c ltou = (0,1) = (no,yes) convert lower case to upper
      data ltou /1/
      code = ' '
      ic = 0
      value = 0.
      id = 0
      cval = ' '
      irc = 0
      im = len(input)
c find if there is an equal sign
      i = 0
   10 i = i + 1
      if (i.gt.im) go to 20
      if (input(i:i).eq.'=') go to 20
      go to 10
   20 nv = i - 1
c no variable found
      if (nv.lt.1) then
         irc = 1
         return
      endif
      call findv(input(1:nv),code,ltou,ic,irc)
      if (irc.eq.1) return
      ns = nv + 2
c no value found
      if (ns.gt.im) then
         irc = 2
c evaluate value
      else
         cval = input(ns:im)
         call evalc(input(ns:im),ival,val,id)
         value = val
c no digits found
         if (id.eq.0) irc = 3
      endif
      return
      end
      subroutine rpars(input,icode,range,itwo)
c this subroutine parses input and looks for the form: lower < upper.
c the numerical value of lower and upper are returned in range(1) and
c range(2), respectively, ignoring blanks. if the '<' symbol is missing,
c then the value is assumed to be lower.
c if no numerical value is found for either lower or upper, or if the
c value contains characters other than numbers, or if upper < lower,
c then default values for lower and upper are used.
c input arguments: input, icode, itwo
c input = input character string to be parsed
c icode = code for variable types,1=integer,2=floating point,
c    3=character,4=special, plus 8 if variable needed for dimensions
c range = minimum/maximum allowed range of value
c itwo = dimension of range array
      character*(*) input
      dimension range(itwo)
      i = 0
      im = len(input)
      itype = icode - (icode/8)*8
c integer type
      if (itype.eq.1) then
         range(1) = -16777215
         range(2) = 16777215
c defaults
      else
         range(1) = -1.0e+35
         range(2) = 1.0e+35
      endif
c find if there is a less than sign
   10 i = i + 1
      if (i.gt.im) go to 20
      if (input(i:i).eq.'<') go to 20
      go to 10
   20 nv = i - 1
      if (nv.ge.1) then
c evaluate lower
         call evalc(input(1:nv),ival,val,id)
         if ((id.gt.0).and.(val.le.range(2))) range(1) = val
      endif
      ns = nv + 2
      if (ns.gt.im) return
c evaluate upper
      call evalc(input(ns:im),ival,val,id)
      if ((id.gt.0).and.(val.ge.range(1))) range(2) = val
      return
      end
      subroutine findvn(code,codev,nvars,nvar)
c this subroutine returns variable index into array of variable names.
c code is input array which has nvars valid entries.
c nvar is output index, zero if no entry is found.
c performs a circular search from last reference
c input arguments: code, codev, nvars
c code = variable symbol character array
c codev = variable symbol to be located
c nvars = number of variables
c nvar = index into code array
      character*8 code(nvars)
      character*(*) codev
      save iv
      data iv /0/
      i = 0
c main loop
   10 i = i + 1
      iv = iv + 1
c found it
      if (code(iv).eq.codev) then
         nvar = iv
         if (iv.ge.nvars) iv = 0
         return
      endif
      if (iv.ge.nvars) iv = 0
      if (i.lt.nvars) go to 10
c not found
      nvar = 0
      return
      end
      subroutine findv(input,code,ltou,ic,irc)
c this subroutine parses input and looks for valid variable name
c the variable found is returned in code, and its length is in ic.
c blanks are ignored, lower case is translated to upper case if ltou
c switch is set.  characters other than letters and numbers stop the
c parsing and return an error code of 1, as does null input.
c input arguments: input, ltou
c input = input string to be parsed
c code = character string containing legal variable name
c ltou = (0,1) = (no,yes) convert lower case to upper
c ic = length of variable name
c irc = (0,1) = (normal,illegal) variable
      character*(*) input
      character*(*) code
      character*1 blank,c0,c9,ca,cz,cl
      save blank,c0,c9,ca,cz,cl
c ca and cz must be in upper case, cl must be in lower case
      data blank,c0,c9,ca,cz,cl /' ','0','9','A','Z','a'/
      irc = 0
      ic = 0
      code = blank
c set limit values for valid codes
      ib = ichar(blank)
      i0 = ichar(c0)
      i9 = ichar(c9)
      ia = ichar(ca)
      iz = ichar(cz)
      il = ichar(cl) - ia
      nv = len(input)
      nc = len(code)
      i = 0
c scan through input
   10 i = i + 1
      if (i.gt.nv) go to 30
      iv = ichar(input(i:i))
c ignore blanks
      if (iv.eq.ib) go to 10
c numbers are ok
      if ((iv.ge.i0).and.(iv.le.i9)) go to 20
c upper case is ok
      if ((iv.ge.ia).and.(iv.le.iz)) go to 20
c convert lower case to upper case and check if ok
      ivt = iv - il
      if (ltou.eq.1) iv = ivt
      if ((ivt.ge.ia).and.(ivt.le.iz)) go to 20
c illegal variable name
      irc = 1
      return
c character is ok
   20 ic = ic + 1
      if (ic.le.nc) then
         code(ic:ic) = char(iv)
         go to 10
      endif
   30 if (ic.eq.0) irc = 1
      return
      end
      subroutine evalc(cval,ival,val,id)
c this subroutine evaluates numerical value of character string.
c both integer, floating point, and e-format numbers are accepted.
c for integer, returns both integer and real values.
c for floating point, returns both real result and its integer part.
c for e-format, returns real results and integer exponent.
c id is the number of digits found.  illegal character terminates
c evaluation.
c input argument: cval
c cval = input string to be evaluated
c ival/val = integer/real numeric value of string
c id = number of digits of numeric value of string
      character*1 v
      character*(*) cval
      data ib /10/
      is = ichar('0')
      im = len(cval)
      num = 0
      norm = 1
      mnum = num
      mnorm = norm
      if = 0
      ie = 0
      id = 0
      i = 0
c main loop
   10 i = i + 1
      if (i.gt.im) go to 20
      v = cval(i:i)
      iv = ichar(v) - is
      if ((iv.ge.0).and.(iv.le.9)) then
         num = iv + ib*num
         id = id + 1
         if (if.eq.1) norm=ib*norm
         go to 10
      endif
      if ((v.eq.' ').or.(v.eq.'+')) go to 10
      if (v.eq.'-') then
         norm = -norm
         go to 10
      endif
      if (ie.eq.1) go to 20
      if ((v.eq.'.').and.(if.eq.0)) then
         if = 1
         go to 10
      endif
      if ((v.eq.'E').or.(v.eq.'e')) then
         mnum = num
         mnorm = norm
         num = 0
         norm = 1
         ie = 1
         if = 0
         go to 10
      endif
   20 ival = num/norm
      if (ie.eq.0) val = real(num)/real(norm)
      if (ie.eq.1) val = (real(mnum)*10.**ival)/real(mnorm)
      return
      end
      subroutine frmtv(code,icode,value,cvalue,iquote,chrv,lv)
c this subroutine writes a variable name and its value into a character.
c variable symbol is in code, its type in icode, and its value in value,
c for numerical variables and in cvalue for character or special types.
c output is in character chrv, and the length of output is in lv.
c input arguments: code, icode, value, cvalue, iquote
c code = variable symbol character
c icode = code for variable types,1=integer,2=floating point,
c    3=character,4=special, plus 8 if variable needed for dimensions
c value = default value of numeric variable
c cvalue = value of character or special variable
c iquote = (0,1) = (no,yes) include quote marks in character strings
c chrv = string output, in form the variable = value
c lv = length of output string
      character*(*) code
      character*(*) cvalue
      character*(*) chrv
      character*20 cfmt
      character*1 blank,ca,cz,cl
      save blank,ca,cz,cl,kutol,big
c ca and cz must be in upper case, cl must be in lower case
      data blank,ca,cz,cl /' ','A','Z','a'/
c kutol = (0,1) = (no,yes) convert upper case to lower
      data kutol /0/
      data big /1000./
c set limit values for valid codes
      ia = ichar(ca)
      iz = ichar(cz)
      il = ichar(cl) - ia
      cfmt = '(1x,a ,2h =,        '
      chrv = blank
      lv = 0
      im = len(code)
      nv = 0
      itype = icode - (icode/8)*8
c illegal type
      if ((itype.lt.1).or.(itype.gt.4)) return
c find end of variable
   10 nv = nv + 1
      if (nv.gt.im) go to 20
      if (code(nv:nv).eq.blank) go to 20
      if (kutol.eq.0) go to 10
c convert to lower case, if requested
      iv = ichar(code(nv:nv))
      if ((iv.ge.ia).and.(iv.le.iz)) code(nv:nv) = char(iv+il)
      go to 10
   20 if (nv.gt.1) nv = nv - 1
      if (nv.gt.9) nv = 9
      cfmt(6:6) = char(nv + ichar('0'))
      absv = abs(value)
c integer type
      if (itype.eq.1) then
c small integer
         if (absv.lt.big) then
            cfmt(13:15) = 'i4)'
            lv = nv + 7
c large integer
         else
            cfmt(13:15) = 'i8)'
            lv = nv + 11
         endif
         ivalue = value
         write (chrv,cfmt) code(1:nv), ivalue
c floating point type
      elseif (itype.eq.2) then
c few digits
        if ((absv.lt.big).and.(absv.eq.(real(int(big*absv))/big))) then
            cfmt(13:17) = 'f8.3)'
            lv = nv + 11
c many digits
         else
            cfmt(13:18) = 'e14.7)'
            lv = nv + 17
         endif
         write (chrv,cfmt) code(1:nv), value
c character or special type
      else
c no quote marks
         if (iquote.eq.0) then
            cfmt(13:16) = 'a20)'
            lv = nv + 23
c put quotes around character strings
         else
            cfmt(10:11) = '='''
            cfmt(13:20) = 'a20,1h'')'
            lv = nv + 24
         endif
         write (chrv,cfmt) code(1:nv), cvalue
      endif
      return
      end
      subroutine frmtn(code,chrv,lv)
c this subroutine writes a variable name into a character.
c variable symbol is in code.
c output is in character chrv, and the length of output is in lv.
c input argument: code
c code = variable symbol character
c chrv = string output, in form the variable = value
c lv = length of output string
      character*(*) code
      character*(*) chrv
      character*8 cfmt
      character*1 blank,ca,cz,cl
      save blank,ca,cz,cl,kutol
c ca and cz must be in upper case, cl must be in lower case
      data blank,ca,cz,cl /' ','A','Z','a'/
c kutol = (0,1) = (no,yes) convert upper case to lower
      data kutol /1/
c set limit values for valid codes
      ia = ichar(ca)
      iz = ichar(cz)
      il = ichar(cl) - ia
      cfmt = '(1x,a ) '
      chrv = blank
      lv = 0
      im = len(code)
      nv = 0
c find end of variable
   10 nv = nv + 1
      if (nv.gt.im) go to 20
      if (code(nv:nv).eq.blank) go to 20
      if (kutol.eq.0) go to 10
c convert to lower case, if requested
      iv = ichar(code(nv:nv))
      if ((iv.ge.ia).and.(iv.le.iz)) code(nv:nv) = char(iv+il)
      go to 10
   20 if (nv.gt.1) nv = nv - 1
      if (nv.gt.9) nv = 9
      cfmt(6:6) = char(nv + ichar('0'))
      lv = nv + 1
      write (chrv,cfmt) code(1:nv)
      return
      end
      subroutine dsvarf(code,icode,value,cvalue,nvars,ncmax,kstrt,knum,i
     1sel,iquote,lvm,ax,ay,space)
c this subroutine writes selected variables in a namelist format, that
c is, in the form variable = value.
c input arguments: all
c code = variable symbol character array
c icode = code array for variable types,1=integer,2=floating point,
c    3=character,4=special, plus 8 if variable needed for dimensions
c value = array containing default value of numeric variable
c cvalue = array containing value of character or special variable
c nvars = number of variables
c ncmax = maximum number of character variables
c kstrt = index of first variable to be searched
c knum = number of variables to be searched
c isel = (-1,0,1) = select (all,icode>8,icode<8) variables,
c iquote = (0,1) = (no,yes) include quote marks in character strings
c lvm = maximum length of variable display line
c ax/ay = location of last character write
c space = vertical spacing between characters
      character*8 code(nvars)
      character*20 cvalue(ncmax)
      character*80 chr
      character*32 chrv
      character*1 c1,c2
      dimension icode(nvars), value(nvars)
      save c1,c2,komma
      data c1,c2 /',',' '/
c komma = (0,1) = delimit variables with (blank,comma)
      data komma /0/
      irc = 0
      nc = 1
      lvs = 0
      chr = ' '
      ks = kstrt - 1
c main loop
      do 10 k = 1, knum
      k1 = k + ks
      itc = icode(k1)/8
c skip unselected items
      if (itc.eq.isel) go to 10
c if character find location in cvalue array
      if ((icode(k1) - 8*itc).ge.3) nc = value(k1)
c get variable name and value string
      call frmtv(code(k1),icode(k1),value(k1),cvalue(nc),iquote,chrv,lv)
      if (lv.gt.lvm) lv = lvm
      lvt = lvs + lv
c write out line to avoid overflow
      if (lvt.gt.lvm) then
c replace comma with blank
         if (komma.eq.0) chr(lvs:lvs) = c2
         ay = ay - space
         if (ay.lt.0.) return
c draw text
         call gtx(ax,ay,chr(1:lvs))
         lvs = 0
         chr = ' '
         lvt = lv
      endif
      ls1 = lvs + 1
c add variable name and value string to output string
      chr(ls1:ls1+lv) = chrv(1:lv)
      lvs = lvt + 1
c add comma
      chr(lvs:lvs) = c1
   10 continue
c flush last line
      if (lvs.gt.0) then
c replace comma with blank
         if (komma.eq.0) chr(lvs:lvs) = c2
         ay = ay - space
         if (ay.lt.0.) return
c draw text
         call gtx(ax,ay,chr(1:lvs))
      endif
      return
      end
      subroutine opsysc(input)
c this subroutine executes system command
      character*(*) input
c     character*1 c
   91 format (a1)
      lenc = len(input)
c     call tsolnk(5,input,lenc,irc,irn,iabend)
c     ierr = lib$spawn(input)
c     ierr = system(input)
c     read (5,91,end=10) c
   10 return
      end
      subroutine helper(ax,ay,space)
c this subroutine displays help for simulation program
c input arguments: all
c ax/ay = location of last character write
c space = vertical spacing between characters
      parameter(nhelp=14)
      character*60 chr(nhelp)
      data chr /'     this program modifies and updates input parameters
     1  ',' type a menu number to select the corresponding menu page   '
     2,' to change parameters, type name=number, where name is the  '
     3,' parameter to be changed and number is the new value, e.g.,'
     4,' indx=7                                                    '
     5,' to obtain the meaning of a parameter, type its name alone '
     6,' to display all the current parameters, hit carriage return'
     7,' to omit headings in display of parameters, type style=0   '
     8,' to include headings in display of parameters, type style=1'
     9,' to obtain extended help for variables, type style=2       '
     a,' to return to the main input parameter menu, type menu     '
     b,' to abort the modification process, type quit              '
     c,' to run program when parameters are satisfactory, type run '
     d,' to save a copy current input to new file, type save       '/
      do 10 i = 1, nhelp
      ay = ay - space
c draw text
      if (ay.ge.0.) call gtx(ax,ay,chr(i))
   10 continue
      return
      end
      subroutine helpvf(code,ax,ay,space,iunit)
c this subroutine finds extended help for variable code in variable
c definition file
c input arguments: all
c code = variable symbol character array
c ax/ay = location of last character write
c space = vertical spacing between characters
c iunit = file unit number for variable description file, assumed open
      character*(*) code
      character*80 c
      character*8 codev,cval
      character*1 c1,c2,c3,c4,c5,c6,c8,c9,v
      save c1,c2,c3,c4,c5,c6,c8,c9
   91 format (a80)
      data c1,c2,c3,c4,c5,c6,c8,c9 /'@','!','%','&','$','*','?','#'/
      lenc = len(c)
      iflg = 0
      if (iunit.le.0) then
         iflg = 1
         go to 40
      endif
c find variable in file
   10 read (iunit,91,end=40) c
      v = c(1:1)
      if ((v.eq.c3).or.(v.eq.c4).or.(v.eq.c5).or.(v.eq.c6)) go to 20
      go to 10
   20 is = 2
      if (c(2:2).eq.c9) is = 3
      call vpars(c(is:lenc),codev,ic,value,id,cval,irc)
      if (irc.eq.1) go to 10
      if (code.ne.codev) go to 10
c variable found
      iflg = 1
   30 read (iunit,91,end=40) c
      v = c(1:1)
      if ((v.eq.c3).or.(v.eq.c4).or.(v.eq.c5).or.(v.eq.c6)) go to 40
      if ((v.eq.c1).or.(v.eq.c2)) go to 40
      if (v.ne.c8) go to 30
c found help
      iflg = 2
      is = 2
      if (c(2:2).eq.c9) is = 3
      ay = ay - space
c draw text
      if (ay.ge.0.) call gtx(ax,ay,' '//c(is:lenc))
      go to 30
c errors
   40 if (iflg.eq.0) then
         ay = ay - space
c draw text
         if (ay.ge.0.) call gtx(ax,ay,' variable not found')
      elseif (iflg.eq.1) then
         ay = ay - space
c draw text
         if (ay.ge.0.) call gtx(ax,ay,' help not found')
      endif
      return
      end
      subroutine WPARMS(iunit,irc)
c this subroutine displays variables in namelist file with unit=iunit
c input argument: iunit
c iunit = file unit number for output data, assumed to be open
c irc = return code (0 = normal return)
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c idwk = workstation identifier
c rx, ry = ndc coordinates of upper-right corner of workstation window
      character*60 chr
      save lvm,csize
   91 format (a60)
c lvm = maximum length of variable display line
      data lvm,csize /59,0.024/
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
      if (ay.lt.chh) go to 40
      ax = csize*rx
c write namelist variables
   10 read (iunit,91,end=40) chr
      if ((chr(2:2).eq.'&').or.(chr(2:2).eq.'/')) go to 10
   20 ay = ay - 1.5*chh
      if (ay.lt.0.) go to 30
c draw text
      call gtx(ax,ay,chr(2:lvm))
      go to 10
c page full
c update workstation, perform
   30 call guwk(idwk,1)
c read code from input device, if present
      call readrc(irc)
c clear workstation, always
      call gclrwk(idwk,1)
      ay = ry
      go to 20
c end of data
c update workstation, perform
   40 call guwk(idwk,1)
c read code from input device, if present
      call readrc(irc)
      return
      end
      subroutine MENUPS(iunit,bvl,bvr,ndp,ndv,ndd,jnmv,jpro,jps,nplot,am
     1odex,freq,trmp,toff,ntr,vtest,ibcs,anx,edge1,edge2)
c this subroutine interactively modifies and updates input variables.
c a variable description file (read from unit=iunit) is required.
c written for the sun - viktor k. decyk, ucla
      parameter (lpmax=3,lgmax=10,nvmax=50,ncmax=5,itwo=2)
      parameter (nvarg=16)
      character*8 code(nvmax)
      character*20 cvalue(ncmax)
      character*72 helpv(nvmax), group(lgmax), page(lpmax)
      dimension icode(nvmax), ig(lgmax), ip(lpmax)
      dimension value(nvmax), range(itwo,nvmax)
      dimension nvp(nvarg)
      save code,cvalue,helpv,group,page,icode,ig,ip,value,range,nvp
      save istart
   91 format (' error: no variables found')
      data istart /0/
      if (istart.eq.0) then
         open(unit=iunit,file='varinr',form='formatted',status='old')
c create variable list and headings from variable list file
         call PARSVF(page,ip,group,ig,code,helpv,icode,value,range,cvalu
     1e,lpmax,lgmax,nvmax,ncmax,itwo,npage,ngroup,nvars,iunit)
         rewind iunit
         if (nvars.le.0) then
            write (6,91)
            return
         endif
c find indices for variables passed in argument
         call findvn(code,'BVL     ',nvars,nvp(1))
         call findvn(code,'BVR     ',nvars,nvp(2))
         call findvn(code,'NDP     ',nvars,nvp(3))
         call findvn(code,'NDV     ',nvars,nvp(4))
         call findvn(code,'NDD     ',nvars,nvp(5))
         call findvn(code,'JNMV    ',nvars,nvp(6))
         call findvn(code,'JPRO    ',nvars,nvp(7))
         call findvn(code,'JPS     ',nvars,nvp(8))
         call findvn(code,'NPLOT   ',nvars,nvp(9))
         call findvn(code,'AMODEX  ',nvars,nvp(10))
         call findvn(code,'FREQ    ',nvars,nvp(11))
         call findvn(code,'TRMP    ',nvars,nvp(12))
         call findvn(code,'TOFF    ',nvars,nvp(13))
         call findvn(code,'EDGE    ',nvars,nvp(14))
         call findvn(code,'NTR     ',nvars,nvp(15))
         call findvn(code,'VTEST   ',nvars,nvp(16))
         istart = 1
      endif
c update variable list with values passed in argument
      edge = edge1
      if (nvp(1).gt.0) value(nvp(1)) = bvl
      if (nvp(2).gt.0) value(nvp(2)) = bvr
      if (nvp(3).gt.0) value(nvp(3)) = ndp
      if (nvp(4).gt.0) value(nvp(4)) = ndv
      if (nvp(5).gt.0) value(nvp(5)) = ndd
      if (nvp(6).gt.0) value(nvp(6)) = jnmv
      if (nvp(7).gt.0) value(nvp(7)) = jpro
      if (nvp(8).gt.0) value(nvp(8)) = jps
      if (nvp(9).gt.0) value(nvp(9)) = nplot
      if (nvp(10).gt.0) value(nvp(10)) = amodex
      if (nvp(11).gt.0) value(nvp(11)) = freq
      if (nvp(12).gt.0) value(nvp(12)) = trmp
      if (nvp(13).gt.0) value(nvp(13)) = toff
      if (nvp(14).gt.0) value(nvp(14)) = edge
      if (nvp(15).gt.0) value(nvp(15)) = ntr
      if (nvp(16).gt.0) value(nvp(16)) = vtest
c modify variables
   10 call MENUVFD(page,ip,group,ig,code,helpv,icode,value,range,cvalue,
     1npage,ngroup,nvars,ncmax,itwo,icmpl,ircopy,iunit,irc)
c update variables passed in argument with values from variable list
      if (nvp(1).gt.0) bvl = value(nvp(1))
      if (nvp(2).gt.0) bvr = value(nvp(2))
      if (nvp(3).gt.0) ndp = value(nvp(3))
      if (nvp(4).gt.0) ndv = value(nvp(4))
      if (nvp(5).gt.0) ndd = value(nvp(5))
      if (nvp(6).gt.0) jnmv = value(nvp(6))
      if (nvp(6).gt.0) jpro = value(nvp(7))
      if (nvp(8).gt.0) jps = value(nvp(8))
      if (nvp(9).gt.0) nplot = value(nvp(9))
      if (nvp(10).gt.0) amodex = value(nvp(10))
      if (nvp(11).gt.0) freq = value(nvp(11))
      if (nvp(12).gt.0) trmp = value(nvp(12))
      if (nvp(13).gt.0) toff = value(nvp(13))
      if (nvp(14).gt.0) edge = value(nvp(14))
      if (nvp(15).gt.0) ntr = value(nvp(15))
      if (nvp(16).gt.0) vtest = value(nvp(16))
c edge cannot be increased
      if (edge.gt.edge1) then
         edge = edge1
         if (nvp(14).gt.0) value(nvp(14)) = edge
c        go to 10
      endif
      edge1 = edge
      if (ibcs.eq.7) edge1 = 0.
      edge2 = anx - edge1
      return
      end
      subroutine MENUP2(iunit,ndp,ndv,ndd,jnmv,jpro,jps,nplot,ntr,ibcs,a
     1nx,edge1,edge2)
c this subroutine interactively modifies and updates input variables.
c a variable description file (read from unit=iunit) is required.
c written for the sun - viktor k. decyk, ucla
      parameter (lpmax=3,lgmax=10,nvmax=50,ncmax=5,itwo=2)
      parameter (nvarg=9)
      character*8 code(nvmax)
      character*20 cvalue(ncmax)
      character*72 helpv(nvmax), group(lgmax), page(lpmax)
      dimension icode(nvmax), ig(lgmax), ip(lpmax)
      dimension value(nvmax), range(itwo,nvmax)
      dimension nvp(nvarg)
      save code,cvalue,helpv,group,page,icode,ig,ip,value,range,nvp
      save istart
   91 format (' error: no variables found')
      data istart /0/
      if (istart.eq.0) then
         open(unit=iunit,file='varinr',form='formatted',status='old')
c create variable list and headings from variable list file
         call PARSVF(page,ip,group,ig,code,helpv,icode,value,range,cvalu
     1e,lpmax,lgmax,nvmax,ncmax,itwo,npage,ngroup,nvars,iunit)
         rewind iunit
         if (nvars.le.0) then
         write (6,91)
         return
         endif
c find indices for variables passed in argument
         call findvn(code,'NDP     ',nvars,nvp(1))
         call findvn(code,'NDV     ',nvars,nvp(2))
         call findvn(code,'NDD     ',nvars,nvp(3))
         call findvn(code,'JNMV    ',nvars,nvp(4))
         call findvn(code,'JPRO    ',nvars,nvp(5))
         call findvn(code,'JPS     ',nvars,nvp(6))
         call findvn(code,'NPLOT   ',nvars,nvp(7))
         call findvn(code,'EDGE    ',nvars,nvp(8))
         call findvn(code,'NTR     ',nvars,nvp(9))
         istart = 1
      endif
c update variable list with values passed in argument
      edge = edge1
      if (nvp(1).gt.0) value(nvp(1)) = ndp
      if (nvp(2).gt.0) value(nvp(2)) = ndv
      if (nvp(3).gt.0) value(nvp(3)) = ndd
      if (nvp(4).gt.0) value(nvp(4)) = jnmv
      if (nvp(5).gt.0) value(nvp(5)) = jpro
      if (nvp(6).gt.0) value(nvp(6)) = jps
      if (nvp(7).gt.0) value(nvp(7)) = nplot
      if (nvp(8).gt.0) value(nvp(8)) = edge
      if (nvp(9).gt.0) value(nvp(9)) = ntr
c modify variables
   10 call MENUVFD(page,ip,group,ig,code,helpv,icode,value,range,cvalue,
     1npage,ngroup,nvars,ncmax,itwo,icmpl,ircopy,iunit,irc)
c update variables passed in argument with values from variable list
      if (nvp(1).gt.0) ndp = value(nvp(1))
      if (nvp(2).gt.0) ndv = value(nvp(2))
      if (nvp(3).gt.0) ndd = value(nvp(3))
      if (nvp(4).gt.0) jnmv = value(nvp(4))
      if (nvp(5).gt.0) jpro = value(nvp(5))
      if (nvp(6).gt.0) jps = value(nvp(6))
      if (nvp(7).gt.0) nplot = value(nvp(7))
      if (nvp(8).gt.0) edge = value(nvp(8))
      if (nvp(9).gt.0) ntr = value(nvp(9))
c edge cannot be increated
      if (edge.gt.edge1) then
         edge = edge1
         if (nvp(8).gt.0) value(nvp(8)) = edge
c        go to 10
      endif
      edge1 = edge
      if (ibcs.eq.7) edge1 = 0.
      edge2 = anx - edge1
      return
      end
