C
C *******************************************************************   
C * This program reads and aanalyzes 21x13 pixel summary files.     *   
C * Create cluster length templates for x-/y-fitting procedure      *
C * Implement 2-pass appproach (3/6/06)                             *
C * Store templates at Nx x Ny                                       *
C * Assume that data are corrected for ROC response.                *
C * Use input thickness (zsize).                                    *
C * Add (0,0) point to each fit (8 Feb 2007).                       *
C * Fake z error parameters for short clusters if possible (8 Feb)  *
C * Reverse run processing direction to make faking work with cotb  *
C * reflection scheme (cotbeta > 0 only is stored temp) [1 Mar 07]  *
C * Add pixmax information to both output files.       [21 Oct 08]  *
C * Make adjustable array and template sizes.          [21 Oct 08]  *
C * Add digits to the error parameterization.          [22 Jul 09]  *
C * Modify small cluster error parameterization.       [22 Jul 09]  *
C * Add single pixel processing to estimate Lor widths [31 Oct 09]  *
C * Update for new response modeling input files       [28 Jan 10]  *
C *******************************************************************
C                                                                       
c
      program gen_xy_template
      implicit double precision (a-h,o-z)
      external fqfit
      parameter (nev = 32000, npt = 300, nprm = 5, maxarg = 5)
      parameter (nqpt = 150000)
C  Nx and Nx are the dimensions of the x and y arrays
      parameter (Nx = 21, Ny = 13, NOPT = 2) 
c
      dimension qsum(nev),
     > pixev(Nx,Ny,nev),xh(nev),yh(nev),i2d(nev), j2d(nev)
      dimension ier(nprm), xpar(nprm,4)
      dimension bixin(Nx,Ny)
      dimension pvec(6), dgauss(21)
      dimension xytemp(Nx,Ny,7,7),xytmp2(Nx,Ny,7,7),nxytry(7,7)
      character*30 fname
      character*80 header
      real*4 mpv, sigma, arg, prob, dummy
      real*8 noise, arglis(maxarg)
      common /qfit/ qsig(nqpt), qexp(nqpt), nqfit
      COMMON / UNIT / ISYSRD,ISYSWR,ISYSSA,ISUMMR
      common /fparams/ par(nprm),parerr(nprm), 
     >     emat(nprm,nprm),f0,ipar(nprm),npr,istat,chnam(nprm)
      character*10 chnam
      logical storep, tstsig
      integer NHx, NHy
C
C  Define the storage for the HBOOK histograms
C
      PARAMETER (NPAWC=20000)
      COMMON /PAWC/ HMEMOR(NPAWC)
      CHARACTER*4 CHOPT(NOPT)
      DATA CHOPT/'STA ','FIT '/
C
C  Initialize HBOOK histograms
C
      CALL HLIMIT(NPAWC)
      CALL HBOOK1(101,'Landau prob (i lt NHx)',100,0.,1.,0.)
      CALL HBOOK1(102,'Landau prob (i gt NHx)',100,0.,1.,0.)
      do i=101,102
        CALL HIDOPT(i,'STAT')
      enddo
C
C  Define middle of dimension ranges
C
      NHx = Nx/2+1
      NHy = Ny/2+1 
C
C  Define poldata to process, and read the max error on a run
C
      open(16,file='pix_2t.proc',status='old',err=11111)
      read(16,*) nfile,numrun,noise,q100,q101,fq100,fcls,
     > fgain,rnoise,linear
      close(16)
      print 10, numrun,nfile,q100,q101,fq100,noise,fcls,
     > fgain,rnoise,linear
10    format(1x,'processing',i3,' files starting at file ',i5,/,
     >1x,'thr0 = ',f6.1,', thr1 = ',f6.1,' rms thresh = ',f6.3,
     >' preamp noise ',f6.1,/,1x,'cluster noise frac = ',f6.3,
     >', gain noise frac = ',f6.3,', readout noise = ',f6.1,/,1x,
     >'linear resp = ',i1)
      tstsig = .true.
      if(linear.eq.2) tstsig = .false.
C
C  Define rms of noise
C
      rten=10.d0
      thr= q100
      thr10=0.1*thr
C
C  Initialize random number generator
C
      lux = 3
      ntotin = 0
      nto2in = 0
      call rluxgo(lux,1234567,ntotin,nto2in)
C
C  Flag if any parameters are set aside for short template usage
C
      storep = .false.
C
C  Open the pixel summary files 
C
      lfile = nfile+numrun-1
      do ifile=lfile,nfile,-1
        if(ifile.lt.10000) then
          write(fname,90) ifile
90        format('template_events_d',I4.4,'.out')
        else
          write(fname,91) ifile
91        format('template_events_d',I5.5,'.out')
        endif
        open(15,file=fname,status='old',err=11111)
C
C  Read in pixel summary data
C
        read(15,100) header
100     format(A80)
C
C  Read in the x and y sizes of the pixels
C
        read(15,*) xsize, ysize, zsize
        print 110, xsize, ysize, zsize
110     format(1x,'xsize = ',f8.2,' um, ysize = ',f8.2,
     >  ' um, zsize = ',f8.2)
        hxsize=xsize/2.d0
        hysize=ysize/2.d0
C
C  Define center z position of the silicon
C
        zcen = zsize/2.d0
C
C  Define center z position of the silicon
C
        zcen = zsize/2.d0
C
C  Begin data loop
C
        n=0
	qavg=0.
200     n=n+1
        if(n.gt.nev) go to 500
201     read(15,*,end=500,err=33333) (pvec(i), i=1,6), nelec
        if(n.eq.1) then
	  cosx=pvec(4)
	  cosy=pvec(5)
	  cosz=pvec(6)
	  clslnx = dabs(zsize*cosx/cosz)
	  clslny = dabs(zsize*cosy/cosz)
	endif
	do j=1,Ny
	  read(15,*) (bixin(i,j),i=1,Nx)
	enddo
C
C  Calculate the hit position at the pixel center plane
C
	xhit = pvec(1) + (zcen-pvec(3))*pvec(4)/pvec(6)
	yhit = pvec(2) + (zcen-pvec(3))*pvec(5)/pvec(6)
        xh(n)=xhit
        yh(n)=yhit
	i2d(n)=int(xhit/xsize*6.+4.5)
        if(i2d(n).lt.1) i2d(n) = 1
        if(i2d(n).gt.7) i2d(n) = 7
	j2d(n)=int(yhit/ysize*6.+4.5)
        if(j2d(n).lt.1) j2d(n) = 1
        if(j2d(n).gt.7) j2d(n) = 7
	qsum(n)=0.
	do j=1,Ny
	  do i=1,Nx
	    qin=bixin(i,j)*rten
	    if(qin.lt.0.d0) qin=0.d0
            pixev(i,j,n)=qin
            qsum(n)=qsum(n)+qin
 	  enddo
        enddo
	qavg=qavg+qsum(n)
        go to 200
500     n=n-1
        nmc=n
        close(15)
	qavg=qavg/float(nmc)
        print 501, ifile, nmc, qavg
501     format(1x,'finish run ',i5,' read ',i5,' events, qavg = ',f9.1)
      do k=1,7
        do l=1,7
          nxytry(k,l) = 0
          do i=1,Nx
            do j=1,Ny
	      xytemp(i,j,k,l)=0.
              xytmp2(i,j,k,l)=0.
            enddo
          enddo
        enddo
      enddo
C
      do n=1,nmc
	if(qsum(n).lt.qavg) then
C
C  count the pixel signals for events with less than average charge
C
C  Make 2-d template shapes
C
           k=i2d(n)
           l=j2d(n)
           nxytry(k,l)=nxytry(k,l)+1
           do i=1,Nx
	     do j=1,Ny
	       xytemp(i,j,k,l)=xytemp(i,j,k,l)+pixev(i,j,n)
	       xytmp2(i,j,k,l)=xytmp2(i,j,k,l)+pixev(i,j,n)**2
       	     enddo
           enddo
           if(k.eq.1) then
             nxytry(7,l)=nxytry(7,l)+1
             do i=2,Nx
	       do j=1,Ny
		 xytemp(i,j,7,l)=xytemp(i,j,7,l)+pixev(i-1,j,n)
		 xytmp2(i,j,7,l)=xytmp2(i,j,7,l)+pixev(i-1,j,n)**2
       	       enddo
             enddo
           endif
           if(k.eq.7) then
             nxytry(1,l)=nxytry(1,l)+1
	     do i=2,Nx
	       do j=1,Ny
		 xytemp(i-1,j,1,l)=xytemp(i-1,j,1,l)+pixev(i,j,n)
		 xytmp2(i-1,j,1,l)=xytmp2(i-1,j,1,l)+pixev(i,j,n)**2
       	       enddo
             enddo
           endif
           if(l.eq.1) then
             nxytry(k,7)=nxytry(k,7)+1
             do i=1,Nx
	       do j=2,Ny
		 xytemp(i,j,k,7)=xytemp(i,j,k,7)+pixev(i,j-1,n)
		 xytmp2(i,j,k,7)=xytmp2(i,j,k,7)+pixev(i,j-1,n)**2
       	       enddo
             enddo
           endif
           if(l.eq.7) then
             nxytry(k,1)=nxytry(k,1)+1
	     do i=1,Nx
	       do j=2,Ny
		 xytemp(i,j-1,k,1)=xytemp(i,j-1,k,1)+pixev(i,j,n)
 		 xytmp2(i,j-1,k,1)=xytmp2(i,j-1,k,1)+pixev(i,j,n)**2
      	       enddo
             enddo
           endif
           if(k.eq.1.and.l.eq.1) then
             nxytry(7,7)=nxytry(7,7)+1
	     do i=2,Nx
	       do j=2,Ny
		 xytemp(i,j,7,7)=xytemp(i,j,7,7)+pixev(i-1,j-1,n)
		 xytmp2(i,j,7,7)=xytmp2(i,j,7,7)+pixev(i-1,j-1,n)**2
       	       enddo
             enddo
           endif
           if(k.eq.7.and.l.eq.7) then
             nxytry(1,1)=nxytry(1,1)+1
             do i=2,Nx
	       do j=2,Ny
		 xytemp(i-1,j-1,1,1)=xytemp(i-1,j-1,1,1)+pixev(i,j,n)
		 xytmp2(i-1,j-1,1,1)=xytmp2(i-1,j-1,1,1)+pixev(i,j,n)**2
       	       enddo
             enddo
           endif
           if(k.eq.1.and.l.eq.7) then
             nxytry(7,1)=nxytry(7,1)+1
             do i=2,Nx
	       do j=2,Ny
		 xytemp(i,j-1,7,1)=xytemp(i,j-1,7,1)+pixev(i-1,j,n)
		 xytmp2(i,j-1,7,1)=xytmp2(i,j-1,7,1)+pixev(i-1,j,n)**2
       	       enddo
             enddo
           endif
           if(k.eq.7.and.l.eq.1) then
             nxytry(1,7)=nxytry(1,7)+1
             do i=2,Nx
	       do j=2,Ny
		 xytemp(i-1,j,1,7)=xytemp(i-1,j,1,7)+pixev(i,j-1,n)
		 xytmp2(i-1,j,1,7)=xytmp2(i-1,j,1,7)+pixev(i,j-1,n)**2
       	       enddo
             enddo
           endif
	endif
      enddo
C
C  Renorm the distributions, find the maximum average pixel signal
C
      pixmax = 0.
      imin = Nx
      imax = 1
      jmin = Ny
      jmax = 1
      do k=1,7
        do l=1,7
          do i=1,Nx
            do j=1,Ny
              sigxy=xytemp(i,j,k,l)/float(nxytry(k,l))
              xytemp(i,j,k,l) = sigxy
	      if(sigxy.gt.pixmax) pixmax=sigxy
              xytmp2(i,j,k,l) = xytmp2(i,j,k,l)/float(nxytry(k,l))
     >        - sigxy**2
               if(sigxy.gt.thr10) then
                if(i.lt.imin) imin = i
                if(i.gt.imax) imax = i
                if(j.lt.jmin) jmin = j
                if(j.gt.jmax) jmax = j
              endif
            enddo
          enddo
        enddo
      enddo
C
C  Analyze error vs signals for entry and exit sides of cluster
C
      ISYSRD=5
      ISYSWR=6
      ISYSSA=7
      call mninit(isysrd,isyswr,isyssa)
C
C  Do the Landau fits
C
      nqfit=0
      do n=1,nmc
        k=i2d(n)
        l=j2d(n)
        if((k.gt.1.and.k.lt.7).and.(l.gt.1.and.l.lt.7)) then
          do i=1,NHx
            do j=1,Ny
	      if(xytemp(i,j,k,l).gt.thr
     >        .and.pixev(i,j,n).gt.0.) then
		if(nqfit.lt.nqpt) nqfit=nqfit+1
		qexp(nqfit)=xytemp(i,j,k,l)
                call doublg(dgauss)
                qsig(nqfit)=pixev(i,j,n)+dgauss(1)*noise
	      endif
	    enddo
	  enddo
        endif
      enddo
      print 900, nqfit
900   format(1x,'Fitting ',i6,' pixels')
c
c     Clear old parameter settings
c
        arglis(1)=0.
        call mnexcm(fqfit,'CLEAR',arglis,0,ierflg,dummy)
        call mnseti('Pixel Landau Fit')
        call mnparm(1,'off',0.d0,0.d0,0.0d0,0.0d0,ier(1))
        call mnparm(2,'gain',0.8d0,0.02d0,0.0d0,0.0d0,ier(2))
	call mnparm(3,'sigf',0.0d0,0.00d0,0.0d0,0.0d0,ier(3))
        call mnparm(4,'noise',1000.d0,10.0d0,0.0d0,0.0d0,ier(4))
	call mnparm(5,'xb',0.0d0,0.00d0,0.0d0,0.0d0,ier(5))
c
c   Do the fit
c
        arglis(1)=2.
        call mnexcm(fqfit,'SET STR',arglis,1,ierflg,dummy)
        arglis(1)=0.
        call mnexcm(fqfit,'SET NOGRADIENT',arglis,0,ierflg,dummy)
        arglis(1)=10000.
        call mnexcm(fqfit,'MIGRAD',arglis,1,ierflg,dummy)
        call mnexcm(fqfit,'HESSE',arglis,1,ierflg,dummy)
C
C  Get the fit status
C
        CALL MNSTAT(FMIN,FEDM,ERRDEF,NPARI,NPARX,istat)
        arglis(1)=0.
        call mnexcm(fqfit,'RETURN',arglis,0,ierflg,dummy)
        do i=1,nprm
          xpar(i,3)=par(i)
        enddo
      do n=1,nqfit
	mpv = par(1) + par(2)*qexp(n)
	sigma = sqrt((mpv*par(3))**2+par(4)**2)
        if(sigma.lt.1.4*thr10) sigma = 1.4*thr10
	arg = (qsig(n) - mpv)/sigma - 0.22278
	prob = dislan(arg)
        call hfill(101,prob,dummy,1.)
      ENDDO
      nqfit=0
      do n=1,nmc
        k=i2d(n)
        l=j2d(n)
        if((k.gt.1.and.k.lt.7).and.(l.gt.1.and.l.lt.7)) then
          do i=NHx,Nx
            do j=1,Ny
	      if(xytemp(i,j,k,l).gt.thr
     >        .and.pixev(i,j,n).gt.0.) then
		if(nqfit.lt.nqpt) nqfit=nqfit+1
		qexp(nqfit)=xytemp(i,j,k,l)
                call doublg(dgauss)
                qsig(nqfit)=pixev(i,j,n)+dgauss(1)*noise
	      endif
	    enddo
	  enddo
        endif
      enddo
      print 900, nqfit
c
c     Clear old parameter settings
c
        arglis(1)=0.
        call mnexcm(fqfit,'CLEAR',arglis,0,ierflg,dummy)
        call mnseti('Pixel Landau Fit')
        call mnparm(1,'off',0.d0,0.d0,0.0d0,0.0d0,ier(1))
        call mnparm(2,'gain',0.8d0,0.02d0,0.0d0,0.0d0,ier(2))
	call mnparm(3,'sigf',0.0d0,0.00d0,0.0d0,0.0d0,ier(3))
        call mnparm(4,'noise',1000.d0,10.d0,0.0d0,0.0d0,ier(4))
	call mnparm(5,'xb',0.0d0,0.00d0,0.0d0,0.0d0,ier(5))
c
c   Do the fit
c
        arglis(1)=2.
        call mnexcm(fqfit,'SET STR',arglis,1,ierflg,dummy)
        arglis(1)=0.
        call mnexcm(fqfit,'SET NOGRADIENT',arglis,0,ierflg,dummy)
        arglis(1)=10000.
        call mnexcm(fqfit,'MIGRAD',arglis,1,ierflg,dummy)
        call mnexcm(fqfit,'HESSE',arglis,1,ierflg,dummy)
C
C  Get the fit status
C
        CALL MNSTAT(FMIN,FEDM,ERRDEF,NPARI,NPARX,istat)
        arglis(1)=0.
        call mnexcm(fqfit,'RETURN',arglis,0,ierflg,dummy)
        do i=1,nprm
          xpar(i,4)=par(i)
        enddo
      do n=1,nqfit
	mpv = par(1) + par(2)*qexp(n)
	sigma = sqrt((mpv*par(3))**2+par(4)**2)
        if(sigma.lt.1.4*thr10) sigma = 1.4*thr10
	arg = (qsig(n) - mpv)/sigma - 0.22278
	prob = dislan(arg)
        call hfill(102,prob,dummy,1.)
      ENDDO
C
C  Open Summary files
C
      write(fname,1000) ifile
1000  format('zptemp_',I5.5,'.txt')
      open(16,file=fname,status='unknown',err=11111)
      write(16,1005) cosx, cosy, cosz
1005  format(3(1x,f9.6))
      write(16,1010) qavg,pixmax,imin,imax,jmin,jmax
1010  format(2(1x,f8.1),4(1x,i2))
        do k=3,4
          write(16,1020) (xpar(i,k), i=1,nprm)
1020      format(5(1x,D15.8))
        enddo
        do l=1,7
          do k=1,7
	    write(16,1030) l, (l*0.166667-0.666667)*ysize, 
     >      k, (k*0.166667-0.666667)*xsize
1030        format(1x,'biny ',i2,', ycenter = ',f8.2,' um,'
     >      ' binx ',i2,', xcenter = ',f8.2,' um')
            do j=1,Ny
	      write(16,1040) (xytemp(i,j,k,l),i=1,Nx)
1040          format(21(1x,f8.1))
            enddo
          enddo
        enddo
        close(16)
C
C  Open a PS meta-file 
C
        write(fname,1100) ifile
1100    format('landau_probs_',I5.5,'.ps')
        OPEN(UNIT=10,file=fname,form='formatted',status='unknown')
C
C  Initialize HPLOT output
C
        CALL HPLINT(0)
C  
C  Increase lines widths and use bold font
C
        CALL HPLSET('*WID',8.0)
        CALL HPLSET('*FON',-60.)
	     CALL HPLSET('HCOL',1203.) 
	     CALL HPLSET('FCOL',0004.) 
C
C  Add statistics options
C
        CALL HPLOPT(CHOPT,NOPT)
C
C  Define and open the meta-file type
C
        CALL IGMETA(-10,-100111)
C
C  Plot all histograms
C
        CALL HPLOT(0,'H',' ',0)
C
C  Close meta-file output
C
        CALL IGMETA(999,0)
        CLOSE(10)
        CALL HPLEND
      enddo
      stop
11111 print 22222
22222 format(1x,'Error opening input file')
      stop
33333 print 44444, runnum
44444 format(1x,'error reading run ',i6)
      stop
      end
C
C
      SUBROUTINE FQFIT(NPAR,GPAR,F,XPAR,IFLAG)
C
C
C  ************************************************************
C  * This routine returns -2ln(L) for the landau fit          *
C  ************************************************************
C
C
      parameter (nqpt = 150000, nprm=5, maxarg=5)
      implicit double precision (a-h,o-z)
      common /qfit/ qsig(nqpt), qexp(nqpt), nqfit
      COMMON / UNIT / ISYSRD,ISYSWR,ISYSSA,ISUMMR
      common /fparams/ par(nprm),parerr(nprm), 
     >     emat(nprm,nprm),f0,ipar(nprm),npr,istat,chnam(nprm)
      character*10 chnam
C
C  Dimension the parameter arrays
C
      DIMENSION XPAR(*), GPAR(*)
      REAL*4 MPV, SIGMA, ARG
C
C  Beginning of executable code
C
      GO TO (1000,2000,4000,4000,9999), IFLAG
C
1000  CONTINUE
C
C  Initialization Section
C
      GO TO 4000
2000  CONTINUE
C
C  Calculate parameter derivatives if necessary
C
 
4000  CONTINUE
C
C  Evaluate the chi-squared
C
      CHISQD=0.
C
C  Loop over all data points
C
      DO i=1,nqfit
	mpv = xpar(1) + xpar(2)*qexp(i)
	sigma = sqrt((mpv*xpar(3))**2+xpar(4)**2)
	arg = (qsig(i) - mpv)/sigma - 0.22278
	pdf = denlan(arg)/sigma
        if(pdf.lt.1.d-20) pdf = 1.d-20
        CHISQD=CHISQD+log(pdf) 
      ENDDO
C
C  Load the answer into the parameter list
C
      F=-2.d0*CHISQD
      IF(IFLAG.NE.3) GO TO 9999
C
C  Termination Section...print and save any results
C
      F0=CHISQD
      DO I=1,nprm
        CALL MNPOUT(I,CHNAM(I),PAR(I),PARERR(I),Z1,Z2,IV)
      ENDDO
C
C  Get the error matrix
C
      CALL MNEMAT(EMAT,nprm)
C
C  Calculate the number of degrees of freedom
C
C      NDOF=nptot-nprm
C
C  Print all points and the fit
C
9999  RETURN
      END
C
      SUBROUTINE DOUBLG(X)
C
C ****************************************************************
C * This routine calculates 2 gaussianly-distributed random nums *
C * Parameters:  X(2) - 2 random numbers (rms = 1.)              *
C ****************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(21), RBUFF(210)
      REAL*4 RVEC(210)
      LOGICAL FCALL
      DATA FCALL/.TRUE./
C
C  Initalize the parameters
C
      IF(FCALL) THEN
        TWOPI=2.d0*ACOS(-1.d0)
        IBASE=210
        FCALL=.FALSE.
      ENDIF
C
C  If all random numbers used up, generate 210 more
C
      IF(IBASE.EQ.210) THEN
        CALL RANLUX(RVEC,210)
        DO I=1,209,2
          ARG=1.d0-RVEC(I)
          IF(ARG.LT.1.d-30) ARG=1.d-30
          R=SQRT(-2.d0*LOG(ARG))
          PHI=TWOPI*RVEC(I+1)
          RBUFF(I)=R*COS(PHI)
          RBUFF(I+1)=R*SIN(PHI)
        ENDDO
        IBASE=0
      ENDIF
      DO I=1,2
        X(I)=RBUFF(IBASE+I)
      ENDDO
      IBASE=IBASE+2
      RETURN
      END
C
C
