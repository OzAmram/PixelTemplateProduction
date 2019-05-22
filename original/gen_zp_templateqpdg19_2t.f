C
C **************************************************************   
C * This program reads and analyzes Nx x Ny pixel summary and  *   
C * template files.  It performs the following tasks:          *
C * 1) generates error information versus cluster charge       *   
C * 2) converts between pixelav and flipped barrel coordinates *   
C * 3) merges all templates into single phi and z files        *   
C * This stage depends upon the readout threshold and noise    *   
C * Make a combined z,phi template.                            *   
C * Analyze errors by topology, use 41 bin (5 pix wdth) techn  *   
C * Add automatic cluster centering and double-size single pix *   
C * Add PSI46v2 response function (10/17/06)                   *
C * Fix possible problem with zero-suppressed ROC signals      *
C * Simulated full charge correction scheme.                   *
C * Add chi**2 to anaysis, plots, and template.                *
C * Reduce outer pseudopixel error by factor of 5              *
C * Try natural normalization for Chi2y(x)                     *
C * Fix pseudopixel norm and remove size 1 clusters from calc  *
C * Fix normalization to remove renorm of signal size          *
C * Determine and store offset information                     *
C * Determine chi**2 for one pixel+pseudopixel clusters        *
C * Try 2-d charge smoothing algorithm before fitting          *
C * Make Q_LF profile histos and fit for the functional form.  *
C * Try 5th order polynomial for Q_LF characterizarion.        *
C * Add minimum charge information to the template             *
C * Add charge scale matching to accommodate cmssw simulation. *
C * Change prompts to reflect larger pt coverage (>0.2 GeV/c). *
C * Add cluster length at fmax*Smax info for both projections. *
C * Allow for high thresholds and broken clusters and add Fpix *
C * bias information.                                          *
C * Use new ptemp and ztemp formats that contain single pixel  *
C * truncation information.                                    *
C * Add CPEGeneric error and bias information to the template  *
C * object.                                                    *
C * Use with version 5.2 of the template code.                 *
C * Change sing-pixel output to use defaults when zero entries *
C * Add single pixel chisquare info to the output              *
C * Add add second (tighter) qmin                              *
C * Add add second threshold and qscale correction for FPix    *
C * Add add charge distribution info.                          *
C * Bigger templates for cosmic y-reco.                        *
C * Add digits to the error parameterization.     [22 Jul 09]  *
C * Fix new Lorentz width problem (expanded temps).[23 Jul 09] *
C * Add Lorentz widths and pixel sizes to spares.  [23 Jul 09] *
C * Add qBin distribution information              [21 Sep 09] *
C * Add single pixel distribution info             [24 Sep 09] *
C * Fix definition of qbin scale for consistency   [26 Sep 09] *
C * Change to new version 16 format to split BPix and FPix     *
C * into different IDs.                            [29 Sep 09] *
C * Get single pixel info from external file       [01 Oct 09] *
C * (use with gen_xy_templateqd4 or higher)        [01 Oct 09] *
C * Use older threshold-sensitive bin definition   [26 Dec 09] *
C * Use new data-based response function           [22 Jan 10] *
C * Clustering algorithm to improve high eta resp  [26 Jan 10] *
C * Change response function to include pixel by pixel and     *
C * cluster by cluster smearing                    [17 Apr 10] *
C * Add Vavilov parameters for merged clusters     [28 Jan 11] *
C * Add chisquare calibration for merged clusters  [31 Jan 11] *
C * Add second threshold for single dcol hits      [14 Dec 11] *
C * Add generic template generation also.          [11 Mar 14] *
C * Expanded Header w/ 2 thresholds, lorentz bias, Qbin defs   *
C * Use old defs of Qbin [1.5,1.0,0.85]            [21 Mar 14] *
C * Update for Phase 1 electronic response         [11 Jan 17] *
C * Update for Phase VCAL/offset definitions       [06 Jan 17] *
C **************************************************************
C                                                                       
c
      program barrel_analz
c
      parameter (nev = 30000, nprm=5, maxarg=5, ns=12, nhprm=6)
C Nx x Ny are the pixelav cluster sizes, Lx + Ly are the template sizes
      parameter (Nx=21, Ny=13, Lx=25, LHx=12, Ly=17, LHy=8)
      implicit double precision (a-h,o-z)
      common /pixdat/ pvec(6,nev), pixel(Lx,Ly,nev), xsum(Lx,nev),
     >                ysum(Ly,nev),xhit(nev),yhit(nev),nelec(nev),
     >                ixc(nev),ixw(nev),iyc(nev),iyw(nev),
     >                ysum1(LHy,nev),ysum2(LHy,nev),xsum1(LHx,nev),
     >                xsum2(LHx,nev),dyx(2,nev),dyxc2(2,nev),
     >                dyxc1(2,nev),xsunc(Lx,nev),ysunc(Ly,nev),
     >                ixc1(nev),ixw1(nev),ixc2(nev),ixw2(nev),
     >                iyc1(nev),iyw1(nev),iyc2(nev),iyw2(nev),
     >                xm(nev),ym(nev),xcen(nev),xwdth(nev),
     >                xfrac(3,nev),yfrac(3,nev),qtotal(nev),
     >                qsmear(nev),qmerge(nev),npix(nev),
     >                nxw(3),nyw(3),xsize,ysize,rnelec,nmc
      real*8  noise, lorwdy, lorwdx, lorbsx, lorbsy
c
	   dimension bixin(Nx,Ny), sigraw(Lx,Ly), xgraw(Lx,Ly)
	   dimension ygraw(Lx,Ly), zgraw(Lx,Ly), pixraw(Lx,Ly)
	   dimension wgraw(Lx,Ly), vgraw(Lx,Ly)
	   dimension xtemp(Nx,9), ztemp(Lx,41)
	   dimension wgauss(21),vgauss(21),xgauss(21), 
     > ygauss(21),zgauss(21)
	   dimension iclust(2,200), ndcol(LHx)
	   dimension ytemp(Ny,9), ptemp(Ly,41), xpar(nprm,4)
	   dimension idhist(50), nqbin(4)
          dimension fqbin(4)
	   real*4 xshist(9,55), dqlfpr(6,10), qmsort(60)
	   real*4 qsort(Lx), qlfy, qlfx, dfy, dfx, chi2h
      dimension rnum(ns), rlimit(ns)
	   integer*4 isort(Lx), imsort(60)
      character*30 fname
	   character*1 ans
      character*80 header	  
	   logical tstsig, lclust(Lx,Ly)
C
C  Define the storage for the HBOOK histograms
C
      PARAMETER (NPAWC=30000)
      COMMON /PAWC/ HMEMOR(NPAWC)
C
C  Keep track of minimum chi**2 values
C
      common /minx/ xmin(600)
	   double precision xmin
C
C  Fitting stuff for the profile histo fits
C
      real*4 par(nhprm),step(nhprm),pmin(nhprm),
     >       pmax(nhprm),sigpar(nhprm)
      CHARACTER*4 CHPOL,CHOPOL
	   DATA CHPOL/'P5  '/, CHOPOL/'FQW '/
C	  DATA CHPOL/'P3  '/, CHOPOL/'FQ '/
      DATA tstsig/.false./
C
C  Store version number in template
C
      data nvers/21/, fmax/0.5d0/, ngver/1/
      ISYSRD=5
      ISYSWR=6
      ISYSSA=7
C
C  Define Useful constants
C
      Lxm1=Lx-1
      Lxm2=Lx-2
      Lxm3=Lx-3
      LHxp1=LHx+1
      Lym1=Ly-1
      Lym2=Ly-2
      Lym3=Ly-3
      LHyp1=LHy+1
      root12 = sqrt(12.)
C
C  Define poldata to process, and read the max error on a run
C
      open(15,file='pix_2t.proc',status='old',err=11111)
      read(15,*) nfile,numrun,noise,q100,q101,fq100,fclus,
     > fgain,rnoise,linear
      close(15)
      print 10, numrun,nfile,q100,q101,fq100,noise,fclus,
     > fgain,rnoise,linear
10    format(1x,'processing ',i3,' files starting at file ',i5,/,
     >1x,'thr0 = ',f6.1,', thr1 = ',f6.1,' rms thresh = ',f6.3,
     >' preamp noise ',f6.1,/,1x,'cluster noise frac = ',f6.3,
     >', gain noise frac = ',f6.3,', readout noise = ',f6.1,/,1x,
     >'linear resp = ',i1)
C
C  Define sampling correction for charge conversion
C
	   rten=10.d0
C
C  Initialize random number generator
C
      lux = 3
      ntotin = 0
      nto2in = 0
      call rluxgo(lux,1234567,ntotin,nto2in)
C	  
C  Open output files for the template summaries
C
	   write(fname,50) nfile
50    format('template_summary_zp',I5.5,'.out')
	   open(21,file=fname,status='unknown',err=11111)
C
C  Open output files for the generic summaries
C
      write(fname,55) nfile
55    format('generror_summary_zp',I5.5,'.out')
      open(22,file=fname,status='unknown',err=11111)
C
C Calculate 50% of threshold in q units and enc noise in adc units
C
      q50=0.5*q100
      q51=0.5*q101
C	  
C  Open input file for the Lorentz Widths
C
C      write(fname,86) nfile
C86    format('lorentz_widths_xy',I5.5,'.out')
C      open(15,file=fname,status='unknown',err=11111)
C      read(15,*) lorwdx, slorwx, nlorwx, cotamn,
C     > lorwdy, slorwy, nlorwy, cotbmn
C      close(15)
C	  
C  Open new style input file for the Lorentz Widths and Biases
C
      write(fname,87) nfile
87    format('lorentz_widths_new',I5.5,'.out')
      open(15,file=fname,status='unknown',err=11111)
      read(15,*) lorwdx, lorbsx, lorwdy, lorbsy
      close(15)
C		  
C	define bpix L1 or L2-4/FPix response parameters  
C		  
      IDtype = 0
      gain = 3.19
      ped = 16.46
      p0 = 0.01218
      p1 = 0.711
      p2 = 203.
      p3 = 148.		
      vcal = 47.
      voffst = 60.
      print 88
88    format(1x,'Use L1 vcal/offset? (y=yes)')
      read 125, ans
      if(ans.eq.'y') then
        vcal = 50.
        voffst = 670.
        print 89
89      format(1x,'Using L1 vcal/offset')
      endif
CC
C  Loop over all runs
C
      lfile = nfile+numrun-1
      do ifile=nfile,lfile
C
C  Open the pixel summary files 
C
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
110     format(1x,'xsize = ',f7.2,' um, ysize = ',f7.2,
     >         ' um, zsize = ',f7.2)
        hxsize=xsize/2.d0
        hysize=ysize/2.d0
C
C  Define center z position of the silicon
C
        zcen = zsize/2.d0
        q10=0.2*q50
        q02=0.2*q10
        if(ifile.eq.nfile) then
C
C  Initialize HBOOK histograms
C
          CALL HLIMIT(NPAWC)
          CALL HBOOKM(header)
		    print 120
120       format(1x,'Write template header info? (y=yes)')
          read 125, ans
125       format(A1)
          if(ans.eq.'y') then
            write(21,100) header
            write(22,100) header
	    print 130
130         format(1x,'Enter ID, NTy(60/28), NTyx(5/3), NTxx(29)',
     >      ' Dtype (0/1), Bfield, Bias Voltage, temperature,',
     >      ' fluence, q-scale')
            read *, id,NTy,NTyx,NTxx,IDt,Bfield,Vbias,temp,
     >              fluenc,qscale
            if(IDt.gt.1) then
              IDtype = IDt
            endif
       write(21,135)id,nvers,Bfield,NTy,NTyx,NTxx,IDtype,Vbias,
     >      temp,fluenc,qscale,q50,-lorwdx,-lorwdy,xsize,ysize,zsize,
     >      q51,-lorbsx,-lorbsy,1.5,1.0,0.85
135         format(1x,I5,1x,I3,1x,f5.2,4(1x,I3),3(1x,f7.2)
     >             ,1x,f7.5,1x,f7.1,2(1x,f7.2),3(1x,f8.3)
     >             ,1x,f7.1,2(1x,f7.2),3(1x,f7.2),)
       write(22,135)id,ngver,Bfield,NTy,NTyx,NTxx,IDtype,Vbias,
     >      temp,fluenc,qscale,q50,-lorwdx,-lorwdy,xsize,ysize,zsize,
     >      q51,-lorbsx,-lorbsy,1.5,1.0,0.85
          endif		
	    else
	      call hreset(0,' ')
	    endif
C
C  Initialize minimum counters
C
		do k=1,600
		  xmin(k)=1.d10
		enddo
		do k=1,60
		   qmsort(k) = (1.-0.01*k)*1.e11
		   imsort(k) = k
		enddo
C
C  Open template files
C
	    write(fname,160) ifile
160     format('ztemp_',I5.5,'.txt')
	    open(16,file=fname,status='unknown',err=11111)
	    read(16,*) cosx, cosy, cosz
161     format(3(1x,f8.5))
	    read(16,*) qavg, sxmax, pixmax
162     format(3(1x,f8.1))
	    print 162, qavg, sxmax, pixmax
C
C  Define sxmaxx to be fmax*sxmax to determine widths of the cluster
C
            sxmaxx = fmax*sxmax
	    do k=1,2
		  read(16,*) (xpar(i,k), i=1,nprm)
164       format(6(1x,D15.8))
	      print 164, (xpar(i,k), i=1,nprm)
	    enddo
	    do k=1,9
		  read(16,165) fname
165       format(a18)
		  read(16,*) (xtemp(i,k),i=1,Nx)
170       format(25(1x,f8.1))
	    enddo
	    close(16)
	    do k=1,9
		  print 170, (xtemp(i,k),i=1,Nx)
	    enddo		
	    zmax=0.
	    do k=1,9
	      do i=1,Nx
	        if(zmax.lt.xtemp(i,k)) zmax=xtemp(i,k)
	      enddo
	    enddo
C
C  Find cluster ends in x at sxmaxx in pixel units
C
            pxfrst=0.
            pxlast=0.
	    do i=1,Nx
              do k=8,1,-1
                if(xtemp(i,k).gt.sxmaxx) then
                  dxsig=xtemp(i,k) - xtemp(i,k+1)
                  if(dxsig.gt.0.) then
                    frac=(xtemp(i,k)-sxmaxx)/dxsig
                    if(frac.gt.1.) frac=1.
                    if(frac.lt.0.) frac=0.
                  else
                    frac=0.
                  endif
                  pxfrst=i-(k+frac-5.)/8.
                  go to 172
                endif
	      enddo
	    enddo
172         continue
	    do i=Nx,1,-1
              do k=2,9
                if(xtemp(i,k).gt.sxmaxx) then
                  dxsig=xtemp(i,k) - xtemp(i,k-1)
                  if(dxsig.gt.0.) then
                    frac=(xtemp(i,k)-sxmaxx)/dxsig
                    if(frac.gt.1.) frac=1.
                    if(frac.lt.0.) frac=0.
                  else
                    frac=0.
                  endif
                  pxlast=i-(k-frac-5.)/8.
                  go to 174
                endif
	      enddo
	    enddo
174         clslnx = pxlast-pxfrst
            if(clslnx.lt.0.) clslnx=0.
C
C  layout 41 entry template by shifting the basic one
C
	    do k=1,9
	      ztemp(1,k+16)=0.d0
	      ztemp(2,k+16)=0.d0
	      ztemp(Lxm1,k+16)=0.d0
	      ztemp(Lx,k+16)=0.d0
	      do i=1,Nx
C
C  This is a shift of template by two pixels
C
		    ztemp(i+2,k+16) = xtemp(i,k)
		  enddo
		enddo
		do k=1,8
		  ztemp(1,k+8)=0.d0
		  ztemp(Lxm2,k+8)=0.d0
		  ztemp(Lxm1,k+8)=0.d0
		  ztemp(Lx,k+8)=0.d0
	      do i=1,Nx
C
C  one pixel shift here
C
		    ztemp(i+1,k+8) = xtemp(i,k)
		  enddo
		enddo
		do k=1,8
		  ztemp(Lxm3,k)=0.d0
		  ztemp(Lxm2,k)=0.d0
		  ztemp(Lxm1,k)=0.d0
		  ztemp(Lx,k)=0.d0
	      do i=1,Nx
C
C  no shift here
C
		    ztemp(i,k) = xtemp(i,k)
		  enddo
		enddo
		do k=2,9
		  ztemp(1,k+24)=0.d0
		  ztemp(2,k+24)=0.d0
		  ztemp(3,k+24)=0.d0
		  ztemp(Lx,k+24)=0.d0
	      do i=1,Nx
C
C  Three pixel shift here
C
		    ztemp(i+3,k+24) = xtemp(i,k)
		  enddo
		enddo
		do k=2,9
		  ztemp(1,k+32)=0.d0
		  ztemp(2,k+32)=0.d0
		  ztemp(3,k+32)=0.d0
		  ztemp(Lx,k+32)=0.d0
	      do i=1,Nx
C
C  Four pixel shift here
C
		    ztemp(i+4,k+32) = xtemp(i,k)
		  enddo
		enddo
	    print 175
175     format(1x)
	    do k=1,41
		  print 170, (ztemp(i,k),i=1,Lx)
	    enddo		
C
C y-template
C
        write(fname,180) ifile
180     format('ptemp_',I5.5,'.txt')
        open(16,file=fname,status='unknown',err=11111)
	    read(16,*) cosx, cosy, cosz
	    read(16,*) qavg, symax, pixmax
	    print 162, qavg, symax, pixmax
C
C  Define symaxx to be fmax*symax to determine widths of the cluster
C
            symaxx = fmax*symax
	    do k=3,4
	      read(16,*) (xpar(i,k), i=1,nprm)
	      print 164, (xpar(i,k), i=1,nprm)
	    enddo
	    do k=1,9
	      read(16,165) fname
	      read(16,*) (ytemp(i,k),i=1,Ny)
	    enddo
	    close(16)
		 print 175
	    do k=1,9
	      print 170, (ytemp(i,k),i=1,Ny)
       enddo		
C
C  Find cluster ends in y at symax/2 in pixel units
C
       pyfrst=0.
       pylast=0.
	    do i=1,Ny
              do k=8,1,-1
                if(ytemp(i,k).gt.symaxx) then
                  dysig=ytemp(i,k) - ytemp(i,k+1)
                  if(dysig.gt.0.) then
                    frac=(ytemp(i,k)-symaxx)/dysig
                    if(frac.gt.1.) frac=1.
                    if(frac.lt.0.) frac=0.
                  else
                    frac=0.
                  endif
                  pyfrst=i-(k+frac-5.)/8.
                  go to 182
                endif
	      enddo
	    enddo
182         continue
	    do i=Ny,1,-1
              do k=2,9
                if(ytemp(i,k).gt.symaxx) then
                  dysig=ytemp(i,k) - ytemp(i,k-1)
                  if(dysig.gt.0.) then
                    frac=(ytemp(i,k)-symaxx)/dysig
                    if(frac.gt.1.) frac=1.
                    if(frac.lt.0.) frac=0.
                  else
                    frac=0.
                  endif
                  pylast=i-(k-frac-5.)/8.
                  go to 184
                endif
	      enddo
	    enddo
184         clslny = pylast-pyfrst
            if(clslny.lt.0.) clslny=0.
C
C  layout 41 entry template by shifting the basic one
C
		do k=1,9
		  ptemp(1,k+16)=0.d0
		  ptemp(2,k+16)=0.d0
		  ptemp(Lym1,k+16)=0.d0
		  ptemp(Ly,k+16)=0.d0
		  do i=1,Ny
		    ptemp(i+2,k+16) = ytemp(i,k)
		  enddo
		enddo
		do k=1,8
		  ptemp(1,k+8)=0.d0
		  ptemp(Lym2,k+8)=0.d0
		  ptemp(Lym1,k+8)=0.d0
		  ptemp(Ly,k+8)=0.d0
		  do i=1,Ny
		    ptemp(i+1,k+8) = ytemp(i,k)
		  enddo
		enddo
		do k=1,8
		  ptemp(Lym3,k)=0.d0
		  ptemp(Lym2,k)=0.d0
		  ptemp(Lym1,k)=0.d0
		  ptemp(Ly,k)=0.d0
		  do i=1,Ny
		    ptemp(i,k) = ytemp(i,k)
		  enddo
		enddo
		do k=2,9
		  ptemp(1,k+24)=0.d0
		  ptemp(2,k+24)=0.d0
		  ptemp(3,k+24)=0.d0
		  ptemp(Ly,k+24)=0.d0
		  do i=1,Ny
		    ptemp(i+3,k+24) = ytemp(i,k)
		  enddo
		enddo
		do k=2,9
		  ptemp(1,k+32)=0.d0
		  ptemp(2,k+32)=0.d0
		  ptemp(3,k+32)=0.d0
		  ptemp(4,k+32)=0.d0
		  do i=1,Ny
		    ptemp(i+4,k+32) = ytemp(i,k)
		  enddo
		enddo
		print 175
		do k=1,41
		  print 185, k, (ptemp(i,k),i=1,Ly)
185     format(1x,i2,25(1x,f7.0))
		enddo		
C
C  Begin data loop
C
        do i=1,3
          nxw(i)=0
          nyw(i)=0
        enddo
        xpsum=0.d0
        npsum=0
        rnelec=0.d0
        n=0
200     n=n+1
        if(n.gt.nev) go to 500
201     read(15,*,end=500,err=33333) (pvec(i,n), i=1,6), nelec(n)
		do j=1,Ny
		  read(15,*) (bixin(i,j),i=1,Nx)
        enddo
C
C  Get gaussianly-distributed random numbers
C
        call triplg(vgauss)
        qsmear(n)=1.+vgauss(1)*fclus
        icol = 0
        if(vgauss(2).lt.0.) icol = 1
        do i=1,LHx
          ndcol(i) = 0
        enddo
        do j=1,Ly
	       do i=1,Lx
	         pixel(i,j,n) = 0.d0
				lclust(i,j) = .true.
	         sigraw(i,j)=0.d0
	         vgraw(i,j)=0.d0
	         wgraw(i,j)=0.d0
	         xgraw(i,j)=0.d0
	         ygraw(i,j)=0.d0
	         zgraw(i,j)=0.d0
	       enddo
	       if(j.gt.2.and.j.lt.Lym1) then
	         call triplg(wgauss)
	         call triplg(xgauss)
			   call triplg(ygauss)
			   call triplg(zgauss)
	         do i=3,Lxm2
		        sigraw(i,j) = rten*bixin(i-2,j-2)
		        wgraw(i,j) = wgauss(i-2)
		        vgraw(i,j) = vgauss(i-2)
		        xgraw(i,j) = xgauss(i-2)
		        ygraw(i,j) = ygauss(i-2)
		        zgraw(i,j) = zgauss(i-2)
		        if(sigraw(i,j).gt.200.) then
				    sig=(sigraw(i,j)+xgraw(i,j)*noise)
				  else
			       sig = 0.
		       endif
               pixraw(i,j) = sig
		       if(sig.lt.(1.+wgraw(i,j)*fq100)*q100) then
					pixel(i,j,n) = 0.d0
               else
                 idcol = (i+icol)/2
                 ndcol(idcol) = ndcol(idcol)+1
		         if(linear.eq.0) then
		           sig=(1.+fgain*ygauss(i-2))*sig+zgauss(i-2)*rnoise
			      else
		           adc=int(p3+p2*tanh(p0*(sig+voffst)/(7.0*vcal)-p1))
		           sig = (1.+fgain*ygauss(i-2))*vcal*gain*(adc-ped)
     >               -voffst+zgauss(i-2)*rnoise
			      endif
                 pixel(i,j,n) = (1.+vgauss(1)*fclus)*sig
		       endif
			  enddo
			endif
		 enddo
C
C  Apply a second threshold to single double column hits
C	
        do i=1,Lx
          idcol = (i+icol)/2
          if(ndcol(idcol).eq.1) then
            do j=1,Ly
              if(pixel(i,j,n) > 0.) then
                if(pixraw(i,j).lt.(1.+wgraw(i,j)*fq100)*q101) 
     >           pixel(i,j,n)=0.
              endif
            enddo
          endif
        enddo
C
C  Look for clustering seed
C
        sigmax = 0.
        do k=1,200
          iclust(1,k) = 0
          iclust(2,k) = 0
        enddo
        do i=1,Lx
          do j=1,Ly
            sig=pixel(i,j,n)/qsmear(n)
            if(sig.gt.sigmax) then
              sigmax = sig
              iclust(1,1) = i
              iclust(2,1) = j
            endif
          enddo
        enddo
C
C  Interpose a clustering algorithm
C
       if(sigmax.lt.5000.) go to 201
       lclust(iclust(1,1),iclust(2,1)) = .false.
		 nclust = 1
		 kmin = 1
300    kmax = nclust
		 do k=kmin,kmax
		   imin = iclust(1,k) - 1
			imax = iclust(1,k) + 1
			jmin = iclust(2,k) - 1
			jmax = iclust(2,k) + 1
		   do j=jmin,jmax
		     do i=imin,imax
			    if((pixel(i,j,n).gt.0.).and.lclust(i,j)) then
				   lclust(i,j) = .false.
					nclust=nclust+1
					iclust(1,nclust) = i
					iclust(2,nclust) = j
					if(nclust.ge.200) go to 350
				 endif
			  enddo
			enddo
		 enddo
		 if(nclust.gt.kmax) then
		   kmin = kmax + 1
		   go to 300
		 endif
350    continue
		 do j=3,Lym2
		   do i=3,Lxm2
			  if(lclust(i,j)) pixel(i,j,n) = 0.
		   enddo
		 enddo
C		 
C  Done clustering		 
C		 
		 qmeas=0.d0
		 npix(n)=0
		 do j=1,Ly
	       if(j.gt.2.and.j.lt.Lym1) then
	         do i=3,Lxm2
				 if(pixel(i,j,n).gt.0.) npix(n)=npix(n)+1
		       qmeas=qmeas+pixel(i,j,n)
	        enddo
			endif
		 enddo
         if(qmeas.lt.q100) go to 201
		 qtotal(n)=qmeas
		 if(qmeas.lt.qmsort(imsort(1))) qmsort(imsort(1)) = qmeas
		 call sortzv(qmsort,imsort,60,1,1,0)
C
C  Histogram the number of generated electrons
C
        relec=float(nelec(n))
        call hfdp1(100,relec,1.d0)
        rnelec=rnelec+qtotal(n)
C
C  Calculate the hit position at the pixel center plane
C
        xhit(n) = pvec(1,n) + (zcen-pvec(3,n))*pvec(4,n)/pvec(6,n)
        yhit(n) = pvec(2,n) + (zcen-pvec(3,n))*pvec(5,n)/pvec(6,n)
C
C  Do the special uncorrected projections to simulate PixelCPEGeneric
C
		do i=1,Lx
		  xsunc(i,n) = 0.d0
		  do j=1,Ly
C
C  Try smoothing pixels after getting correct charge
C
			if(pixel(i,j,n).gt.0.d0) then
			 xsunc(i,n)=xsunc(i,n)+pixel(i,j,n)
			endif
		  enddo
		enddo
		do j=1,Ly
		  ysunc(j,n) = 0.d0
		  do i=1,Lx
			if(pixel(i,j,n).gt.0.d0) then
			  ysunc(j,n)=ysunc(j,n)+pixel(i,j,n)
			endif
		  enddo
		enddo
C
C  Now truncate the signals and do the correct projections
C
		do i=1,Lx
		  xsum(i,n) = 0.d0
		  do j=1,Ly
C
C  Try smoothing pixels after getting correct charge
C
			if(pixel(i,j,n).gt.0.d0) then
			 if(pixel(i,j,n).gt.pixmax) pixel(i,j,n)=pixmax
			 xsum(i,n)=xsum(i,n)+pixel(i,j,n)
			endif
		  enddo
		enddo
		do j=1,Ly
		  ysum(j,n) = 0.d0
		  do i=1,Lx
			if(pixel(i,j,n).gt.0.d0) then
			  ysum(j,n)=ysum(j,n)+pixel(i,j,n)
			endif
		  enddo
		enddo
C
C  calculate the double pixel column and row sums
C
		do i=1,LHx
		  xsum1(i,n) = 0.d0
		  xsum2(i,n) = 0.d0
		  i1=2*i-1
		  i2=2*i
		  do j=1,Ly
		    qin=sigraw(i1,j)+sigraw(i1+1,j)
			 sig=(qin+xgraw(i1,j)*noise)
			 if(sig.ge.(1.+wgraw(i1,j)*fq100)*q100) then
				if(linear.eq.0) then
				  sig=(1.+fgain*ygraw(i1,j))*sig+zgraw(i1,j)*rnoise
				else
				  adc=int(p3+p2*tanh(p0*(sig+voffst)/(7.0*vcal)-p1))
					sig = (1.+fgain*ygraw(i1,j))*vcal*gain*(adc-ped)
     >             -voffst+zgraw(i1,j)*rnoise
				endif
 			  xsum1(i,n)=xsum1(i,n)+(1.+vgraw(3,j)*fclus)*sig
			endif
			qin=sigraw(i2,j)+sigraw(i2+1,j)
			sig=(qin+xgraw(i2,j)*noise)
			 if(sig.ge.(1.+wgraw(i2,j)*fq100)*q100) then
				if(linear.eq.0) then
				  sig=(1.+fgain*ygraw(i2,j))*sig+zgraw(i2,j)*rnoise
				else
				  adc=int(p3+p2*tanh(p0*(sig+voffst)/(7.0*vcal)-p1))
					sig = (1.+fgain*ygraw(i2,j))*vcal*gain*(adc-ped)
     >             -voffst+zgraw(i2,j)*rnoise
				endif
 			  xsum2(i,n)=xsum2(i,n)+(1.+vgraw(3,j)*fclus)*sig
			endif
		  enddo
		enddo
		do j=1,LHy
		  ysum1(j,n) = 0.d0
		  ysum2(j,n) = 0.d0
		  j1=2*j-1
		  j2=2*j
		  do i=1,Lx
		    qin=sigraw(i,j1)+sigraw(i,j1+1)
			 sig=(qin+xgraw(i,j1)*noise)
			 if(sig.ge.(1.+wgraw(i,j1)*fq100)*q100) then
				if(linear.eq.0) then
				  sig=(1.+fgain*ygraw(i,j1))*sig+zgraw(i,j1)*rnoise
				else
				  adc=int(p3+p2*tanh(p0*(sig+voffst)/(7.0*vcal)-p1))
					sig = (1.+fgain*ygraw(i,j1))*vcal*gain*(adc-ped)
     >             -voffst+zgraw(i,j1)*rnoise
				endif
			   ysum1(j,n)=ysum1(j,n)+(1.+vgraw(3,j)*fclus)*sig
			 endif
		    qin=sigraw(i,j2)+sigraw(i,j2+1)
			 sig=(qin+xgraw(i,j2)*noise)
			 if(sig.ge.(1.+wgraw(i,j2)*fq100)*q100) then
				if(linear.eq.0) then
				  sig=(1.+fgain*ygraw(i,j2))*sig+zgraw(i,j2)*rnoise
				else
				  adc=int(p3+p2*tanh(p0*(sig+voffst)/(7.0*vcal)-p1))
					sig = (1.+fgain*ygraw(i,j2))*vcal*gain*(adc-ped)
     >             -voffst+zgraw(i,j2)*rnoise
				endif
			   ysum2(j,n)=ysum2(j,n)+(1.+vgraw(3,j)*fclus)*sig
			 endif
		  enddo
		enddo
C
C  Find cluster beginning and width
C
		ixc(n) = 0
		ixw(n) = 0
		ixl = 0
		do i=1,Lx
		  if(xsum(i,n).gt.0.d0) then
			if(ixw(n).eq.0) ixc(n) = i
		    ixw(n)=ixw(n)+1
			 ixl = i
		  endif
		enddo
C
C  Skip over (and count) any zero pixels
C
		if(ixw(n).ne.(ixl-ixc(n)+1)) go to 201
		iyc(n) = 0
		iyw(n) = 0
		iyl = 0
		do j=1,Ly
		  if(ysum(j,n).gt.0.d0) then
			if(iyw(n).eq.0) iyc(n) = j
			iyw(n)=iyw(n)+1
			iyl = j
		  endif
		enddo
C
C  Skip over (and count) any zero pixels
C
		if(iyw(n).ne.(iyl-iyc(n)+1)) go to 201
C
C  Find cluster beginning and width
C
		ixc1(n) = 0
		ixw1(n) = 0
		ixc2(n) = 0
		ixw2(n) = 0
		do i=1,LHx
		  if(xsum1(i,n).gt.0.d0) then
			if(ixw1(n).eq.0) ixc1(n) = i
		    ixw1(n)=ixw1(n)+1
		  endif
		  if(xsum2(i,n).gt.0.d0) then
			if(ixw2(n).eq.0) ixc2(n) = i
		    ixw2(n)=ixw2(n)+1
		  endif
		enddo
		iyc1(n) = 0
		iyw1(n) = 0
		iyc2(n) = 0
		iyw2(n) = 0
		do j=1,LHy
		  if(ysum1(j,n).gt.0.d0) then
			if(iyw1(n).eq.0) iyc1(n) = j
			iyw1(n)=iyw1(n)+1
		  endif
		  if(ysum2(j,n).gt.0.d0) then
			if(iyw2(n).eq.0) iyc2(n) = j
			iyw2(n)=iyw2(n)+1
		  endif
		enddo
C
C  Calculate the front and back signal fractions
C
		if(ixw(n).eq.1) then
		  xfrac(1,n)=0.d0
		  xfrac(2,n)=0.d0
		  xfrac(3,n)=0.d0
		elseif(ixw(n).eq.2) then
		  xfrac(1,n)=xsum(ixc(n)+1,n)
     >              /(xsum(ixc(n)+1,n)+xsum(ixc(n),n))
		  xfrac(2,n)=xfrac(1,n)
		  xfrac(3,n)=xfrac(1,n)
		else
		  xfrac(1,n)=xsum(ixc(n)+1,n)
     >              /(xsum(ixc(n)+1,n)+xsum(ixc(n),n))
		  il=ixc(n)+ixw(n)-1
		  xfrac(2,n)=xsum(il,n)
     >              /(xsum(il,n)+xsum(il-1,n))
		  xfrac(3,n)=xsum(il,n)
     >              /(xsum(il,n)+xsum(ixc(n),n))
        endif
		if(iyw(n).eq.1) then
		  yfrac(1,n)=0.d0
		  yfrac(2,n)=0.d0
		  yfrac(3,n)=0.d0
		elseif(iyw(n).eq.2) then
		  yfrac(1,n)=ysum(iyc(n)+1,n)
     >              /(ysum(iyc(n)+1,n)+ysum(iyc(n),n))
		  yfrac(2,n)=yfrac(1,n)
		  yfrac(3,n)=yfrac(1,n)
		else
		  yfrac(1,n)=ysum(iyc(n)+1,n)
     >              /(ysum(iyc(n)+1,n)+ysum(iyc(n),n))
		  il=iyc(n)+iyw(n)-1
		  yfrac(2,n)=ysum(il,n)
     >              /(ysum(il,n)+ysum(il-1,n))
		  yfrac(3,n)=ysum(il,n)
     >              /(ysum(il,n)+ysum(iyc(n),n))
		endif
C
C  Catalog cluster sizes
C
		iw=ixw(n)
		if(iw.lt.1) go to 201
		if(iw.gt.3) iw=3
		nxw(iw)=nxw(iw)+1
		iw=iyw(n)
		if(iw.lt.1) go to 201
		if(iw.gt.3) iw=3
		nyw(iw)=nyw(iw)+1

C
C  Calculate the weighted means
C
		xm(n) = 0.d0
		tot=0.d0
		np=ixw(n)
		do i=1,np
		  x=float(ixc(n)+i-13)*xsize
		  xm(n)=xm(n)+xsum(ixc(n)+i-1,n)*x
		  tot=tot+xsum(ixc(n)+i-1,n)
		enddo
		xm(n)=xm(n)/tot
C
C  do y
C
		ym(n) = 0.d0
		tot=0.d0
		np=iyw(n)
		do i=1,np
		  y=float(iyc(n)+i-6)*ysize
		  ym(n)=ym(n)+ysum(iyc(n)+i-1,n)*y
		  tot=tot+ysum(iyc(n)+i-1,n)
		enddo
		ym(n)=ym(n)/tot
C
C  Calculate x position from the cluster center
C
		if(xfrac(1,n).gt.0.d0) then
		  R12=1.d0/xfrac(1,n)-1.d0
		  if(R12.lt.0) R12=0.d0
		  if(R12.gt.1.d0) R12=1.d0
		else
		  R12=1.d0
		endif
		xf=(float(ixc(n))-11.5-R12)*xsize
		R21=xfrac(2,n)/(1.d0-xfrac(2,n))
		if(R21.lt.0) R21=0.d0
		if(R21.gt.1.d0) R21=1.d0
		xb=(float(ixc(n)+ixw(n))-13.5+R21)*xsize
		xcen(n)=(xf+xb)/2.d0
		xwdth(n)=xb-xf
        xpos=xhit(n)
        go to 200
500     n=n-1
        nmc=n
        close(15)
        print 501, ifile, nmc
501     format(1x,'run ',i5,', read ',i5,' events')
        rnelec=rnelec/float(nmc)
C Restore old def to check for consistency
        rlimit(1) = 1.5d0*rnelec
        rlimit(2) = 1.4d0*rnelec
        rlimit(3) = 1.3d0*rnelec
        rlimit(4) = 1.2d0*rnelec
        rlimit(5) = 1.1d0*rnelec
        rlimit(6) = 1.0d0*rnelec
        rlimit(7) = 0.95d0*rnelec
        rlimit(8) = 0.90d0*rnelec
        rlimit(9) = 0.85d0*rnelec
        rlimit(10) = 0.80d0*rnelec
        rlimit(11) = 0.75d0*rnelec
        rlimit(12) = 0.60d0*rnelec
C        rlimit(1) = 1.5d0*qavg
C        rlimit(2) = 1.4d0*qavg
C        rlimit(3) = 1.3d0*qavg
C        rlimit(4) = 1.2d0*qavg
C        rlimit(5) = 1.1d0*qavg
C        rlimit(6) = 1.0d0*qavg
C        rlimit(7) = 0.95d0*qavg
C        rlimit(8) = 0.90d0*qavg
C        rlimit(9) = 0.85d0*qavg
C        rlimit(10) = 0.80d0*qavg
C        rlimit(11) = 0.75d0*qavg
C        rlimit(12) = 0.70d0*qavg
C
C  Determine the average offset of single-pixel clusters
C
        nyone=0
	    dyone=0.
		syone=0.
        nytwo=0
	    dytwo=0.
		sytwo=0.
	    do l=1,ns
	      rnum(l)=0.d0
	    enddo
	    do l=1,4
	      nqbin(l)=0
	    enddo
        do k=1,nmc
	      call hfdp1(500,qtotal(k),1.d0)
          rcorr = qtotal(k)/float(nelec(k))
          call hfdp1(606,rcorr,1.d0)	      
         if(k.gt.1) then
			  qmerge(k) = qtotal(k)+qsmear(k)*qtotal(k-1)/qsmear(k-1)
			  call hfdp1(505,qmerge(k),1.d0)
			endif
	      do l=1,ns
			  if(k.gt.1) then
		       if(qmerge(k).lt.2.*rlimit(l)) rnum(l)=rnum(l)+1.d0
			  endif
	      enddo
	      if(qtotal(k).gt.rlimit(1)) then
                nqbin(1) = nqbin(1)+1
	      endif
	      if(qtotal(k).gt.rlimit(6).and.qtotal(k).le.rlimit(1)) then
                nqbin(2) = nqbin(2)+1
	      endif
	      if(qtotal(k).gt.rlimit(9).and.qtotal(k).le.rlimit(6)) then
                 nqbin(3) = nqbin(3)+1
	       endif
               if(qtotal(k).le.rlimit(9)) then
                 nqbin(4) = nqbin(4)+1
	       endif
C
C  calculate offset from center of pixel
C
          if(iyw(k).eq.1) then	
		    nyone=nyone+1    
	        y0=(iyc(k)-LHyp1)*ysize
			deltay=y0-yhit(k)
	        dyone=dyone+deltay
			syone=syone+deltay**2
		  endif
          if(iyw1(k).eq.1) then	
		    nytwo=nytwo+1    
	        y0=(iyc1(k)*2-(float(Ly)/2.+1.))*ysize
			deltay=y0-yhit(k)
	        dytwo=dytwo+deltay
			sytwo=sytwo+deltay**2
		  endif
          if(iyw2(k).eq.1) then	
		    nytwo=nytwo+1    
	        y0=(iyc2(k)*2-(float(Ly)/2.))*ysize
			deltay=y0-yhit(k)
	        dytwo=dytwo+deltay
			sytwo=sytwo+deltay**2
		  endif
	    enddo
	    if(nyone.le.10) then
	      dyone=lorbsy
	      syone=ysize/root12
	    else
	      dyone=dyone/float(nyone)
	      syone=syone/float(nyone)-dyone**2
		  if(syone.lt.0.) syone=0.
		  syone=sqrt(syone)
	    endif
	    if(nytwo.le.10) then
	      dytwo=lorbsy
	      sytwo=2.*ysize/root12
	    else
	      dytwo=dytwo/float(nytwo)
	      sytwo=sytwo/float(nytwo)-dytwo**2
		  if(sytwo.lt.0.) sytwo=0.
		  sytwo=sqrt(sytwo)
	    endif
	    print 600, nyone, dyone, syone, q50
600     format(1x,i5,' single pixel y-clusters, average offset = ',f7.1,
     >         ', rms = ',f7.1,', 0.5*threshold = ',f7.1)
	    print 605, nytwo, dytwo, sytwo
605     format(1x,i5,' single big pixel y-clusters, average offset = ',
     >         f7.1,', rms = ',f7.1,/)
	   do l=1,ns
		  print 606, rnum(l), l
606     format(1x,f7.0,' merged events less than limit ',i2)
      enddo
C
C  Loop over events and determine y position from linear fit
C
        do k=1,nmc
          relec=qtotal(k)
C
C  Get first position estimate from ORCA estimator
C	    
	      y0=(iyc(k)-LHyp1)*ysize
	      yc = y0+float(iyw(k)-1)*(ysize/2.)
		  ywidth=zsize*pvec(5,k)/pvec(6,k)-lorwdy
		  if(iyw(k).lt.3) then
		      wzin=0.
	       else
		      wzin=(iyw(k)-2)*ysize
	       endif
	       if(iyw(k).eq.1) then
	          shift=0.
	       elseif(iyw(k).ge.2) then
    	      shift=(ysunc(iyc(k)+iyw(k)-1,k)-ysunc(iyc(k),k))
     >        /(ysunc(iyc(k)+iyw(k)-1,k)+ysunc(iyc(k),k))
     >        /2.*(abs(ywidth)-wzin)
	       endif
		   yrec1=yc+shift-lorbsy
C
C  Shift the y-clusters to put the middle in pixels 5 (or 4 if even)
C
           midpix = iyc(k)+iyw(k)/2
		   ishifty = LHyp1 - midpix
		   lstpix = iyc(k)+iyw(k)-1
		   if(ishifty.gt.0) then
		     do i=lstpix,iyc(k),-1
			   ysum(i+ishifty,k) = ysum(i,k)
			   ysum(i,k) = 0.
			  enddo
		   elseif(ishifty.lt.0) then
		     do i=iyc(k), lstpix
			   ysum(i+ishifty,k) = ysum(i,k)
			   ysum(i,k)=0.
			  enddo
		   endif
		   iyc(k)=iyc(k)+ishifty
C
C  Calculate the min first-pass chi**2 in all cases (not just cluster size > 1)
C
C  Search for template that minimizes chi**2
C
			chimny=1.d15
		    lmin=0
			sythr=1.1*symax
			nsort=0
			do i=1,Ly
			  if(ysum(i,k).gt.0.) then
			    nsort=nsort+1
			    qsort(nsort)=ysum(i,k)
			  endif
			enddo
			call sortzv(qsort,isort,nsort,1,0,0)
			if(iyw(k).gt.1) then
		       if(qsort(isort(2)).gt.sythr) sythr=1.01*qsort(isort(2))
			else
	           sythr=1.01*qsort(isort(1))
			endif
C
C  add 1/2 threshold signals to end pixels
C   
			ysum(iyc(k)-2,k)=q10
			ysum(iyc(k)-1,k)=q10
			ysum(iyc(k)+iyw(k),k)=q10
			ysum(iyc(k)+iyw(k)+1,k)=q10
			do l=3,39
C
C  Zero linear fit counters
C	    
			   ss2=0.
			   ssa=0.
			   sa2=0.
			   do i=1,Ly
			     if(ysum(i,k).gt.0.) then
			       ysig=min(ysum(i,k), symax)
			       ysig=max(ysig, q50)
			       if(i.le.LHyp1) then
				  sigi2=xpar(1,3)+xpar(2,3)*ysig+xpar(3,3)*ysig**2
     >                     +xpar(4,3)*ysig**3+xpar(5,3)*ysig**4
				    else
			          sigi2=xpar(1,4)+xpar(2,4)*ysig+xpar(3,4)*ysig**2
     >                     +xpar(4,4)*ysig**3+xpar(5,4)*ysig**4
				    endif
			       if(sigi2.le.0.) sigi2 = q50**2
			       if(sigi2.lt.0..and.ysig.gt.q10) then
	              print 8887, i, sigi2, ysig
8887       format(1x,'i = ',i2,' sigy2 = ',e11.4,' ysignal = ',e11.4)
		     print 8886, 3,xpar(1,3),xpar(2,3),xpar(3,3),
     >               xpar(4,3),xpar(5,3)
		     print 8886, 4,xpar(1,4),xpar(2,4),xpar(3,4),
     >               xpar(4,4),xpar(5,4)
8886       format(1x,'xpar(1-5,',i1,') = ',5(e11.4,1x))
C
C  If tstsig is set true, quit on this error
C
                     if(tstsig) then
                       print 8889, ifile
                       go to 55555
					 else
					   sigi2 = q50**2
					 endif
				    endif
				    if(i.le.(iyc(k)-1).or.i.ge.(iyc(k)+iyw(k))) 
     >                 sigi2=q50**2
			        if(ysum(i,k).gt.sythr) sigi2=(10000.)**2
		            ss2=ss2+ysum(i,k)**2/sigi2
			        ssa=ssa+ysum(i,k)*ptemp(i,l)/sigi2
			        sa2=sa2+ptemp(i,l)**2/sigi2
		          endif
			   enddo
C
C  rnorm1 is the normalization factor
C
			  rnorm1=ssa/ss2
			  chi2=ss2-2.*ssa/rnorm1+sa2/rnorm1**2
			  if(chi2.lt.0.) then
				print 7777, chi2, rnorm1, ss2, ssa, sa2
7777            format(1x,'negative chi2y:',5(2x,e11.4))
			  endif
			  if(chi2.lt.chimny) then
				 chimny=chi2
				 lmin=l
				 ss2min=ss2
				 ssamin=ssa
				 sa2min=sa2
			  endif
			enddo
			if(chimny.lt.0.) chimny=0.
			if(k.eq.1) then
			  ss20=ss2min
			  ssa0=ssamin
			  sa20=sa2min
			else
			  ss2=ss20+ss2min
			  ssa=ssa0+ssamin
			  sa2=sa20+sa2min
			  rnorm1=ssa/ss2
			  chi2=ss2-2.*ssa/rnorm1+sa2/rnorm1**2
	        if(qmerge(k).ge.2.*rlimit(1)) then
		       call hfminx(506,chi2,1.d0) 
		     endif
		     if(qmerge(k).gt.2.*rlimit(6)
     >        .and.qmerge(k).le.2.*rlimit(1)) then
		       call hfminx(507,chi2,1.d0) 
		     endif
		     if(qmerge(k).gt.2.*rlimit(9)
     >        .and.qmerge(k).le.2.*rlimit(6)) then
		       call hfminx(508,chi2,1.d0) 
		     endif
		     if(qmerge(k).le.2.*rlimit(9)) then
		       call hfminx(509,chi2,1.d0) 
		     endif
			  ss20=ss2min
			  ssa0=ssamin
			  sa20=sa2min
			endif
C
C  Use the second stage fit for larger > clusters
C
		   if(iyw(k).gt.1) then
		      yrec2=(0.125*lmin-2.625-ishifty)*ysize
C
C  Make an interpolating pass
C
		     ibinl=lmin-1
		     ibinh=ibinl+2
		     if(ibinl.lt.1) ibinl=1
		     if(ibinh.gt.41) ibinh=41
C
C  Zero linear fit counters
C	    
		     ss2=0.
		     ssa=0.
			 sa2=0.
                 ssba=0.
	         saba=0.
	         sba2=0. 
C
C  Take disributions from first estimate
C	
		     do i=1,Ly
		       if(ysum(i,k).gt.0.) then
			     ysig=min(ysum(i,k), symax)
			     ysig=max(ysig, q50)
			     if(i.le.LHyp1) then
			       sigi2=xpar(1,3)+xpar(2,3)*ysig+xpar(3,3)*ysig**2
     >                  +xpar(4,3)*ysig**3+xpar(5,3)*ysig**4
			     else
			       sigi2=xpar(1,4)+xpar(2,4)*ysig+xpar(3,4)*ysig**2
     >                  +xpar(4,4)*ysig**3+xpar(5,4)*ysig**4
			     endif
				 if(i.le.(iyc(k)-1).or.i.ge.(iyc(k)+iyw(k))) 
     >              sigi2=q50**2
			     if(ysum(i,k).gt.sythr) sigi2=(10000.)**2
			     if(sigi2.le.0.) sigi2 = q50**2
		         ss2=ss2+ysum(i,k)**2/sigi2
			     ssa=ssa+ysum(i,k)*ptemp(i,ibinl)/sigi2
                 sa2=sa2+ptemp(i,ibinl)**2/sigi2
			     ssba=ssba+ysum(i,k)*(ptemp(i,ibinh)-ptemp(i,ibinl))
     >                /sigi2
			     saba=saba+ptemp(i,ibinl)*(ptemp(i,ibinh)-
     >                ptemp(i,ibinl))/sigi2
			     sba2=sba2+(ptemp(i,ibinh)-ptemp(i,ibinl))**2/sigi2
		       endif
	         enddo
		     rat3=(ssba*ssa-ss2*saba)/(ss2*sba2-ssba**2)
		     if(rat3.lt.0.) rat3=0.
		     if(rat3.gt.1.) rat3=1.0
			 rnorm3=(ssa+rat3*ssba)/ss2
			 yrec3=(0.125*ibinl-2.625+rat3*(ibinh-ibinl)*0.125
     >             -ishifty)*ysize
			 chi2y=ss2-2./rnorm3*ssa-2.*rat3/rnorm3*ssba
     >            +(sa2+2.*rat3*saba+rat3**2*sba2)/rnorm3**2
             if(chi2y.lt.0.) chi2y=0.
		   else
		     yrec2=(iyc(k)-LHyp1-ishifty)*ysize-dyone
	             yrec3=yrec2
		     chi2y=chimny
		   endif
		   dy1=yrec1-yhit(k)
		   dy2=yrec2-yhit(k)
		   dy3=yrec3-yhit(k)
		  if(iyw(k).eq.1) then
			call hfminx(106,chi2y,1.d0) 
		  else
	        if(relec.ge.rlimit(1)) then
		       call hfminx(107,chi2y,1.d0) 
		    endif
		    if(relec.gt.rlimit(6).and.relec.le.rlimit(1)) then
		       call hfminx(108,chi2y,1.d0) 
		    endif
		    if(relec.gt.rlimit(9).and.relec.le.rlimit(6)) then
		       call hfminx(109,chi2y,1.d0) 
		    endif
		    if(relec.le.rlimit(9)) then
		       call hfminx(110,chi2y,1.d0) 
		    endif
		  endif
		  dyx(1,k) = dy3
                  dyxc2(1,k) = dy2
                  dyxc1(1,k) = dy1
		  qlfy = 2.d0*yfrac(3,k) - 1.d0
		  dfy = dy3
		  if(iyw(k).gt.1) then 
			call hfill(111,qlfy,dfy,1.)
		    if(relec.gt.rlimit(1)) then
			  call hfill(112,qlfy,dfy,1.)
		    endif
		    if(relec.gt.rlimit(6).and.relec.le.rlimit(1)) then
			  call hfill(113,qlfy,dfy,1.)
		    endif
		    if(relec.gt.rlimit(9).and.relec.le.rlimit(6)) then
			  call hfill(114,qlfy,dfy,1.)
		    endif
		    if(relec.le.rlimit(9)) then
			  call hfill(115,qlfy,dfy,1.)
			endif
		  endif
		enddo  
		ihist=0
		DO I=11,15
		  ID=100+I
		  PAR(1)=0.d0
		  STEP(1)=0.d0
		  PAR(2)=HMAX(ID)/2.d0
		  STEP(2)=0.1*PAR(2)
		  PAR(3)=0.d0
		  STEP(3)=0.d0
		  PAR(4)=0.d0
		  STEP(4)=0.01d0
		  PAR(5)=0.d0
		  STEP(5)=0.d0
		  PAR(6)=0.d0
		  STEP(6)=0.01d0
 		  ihist=ihist+1
		  idhist(ihist)=ID
                  CALL HNOENT(ID, NOENT)
                  IF(NOENT.gt.50) THEN
		    CALL HFITHN(ID,CHPOL,CHOPOL,6,PAR,STEP,
     >                   PMIN,PMAX,SIGPAR,CHI2H)
                    if(CHI2H.lt.1.d6) then
		      dqlfpr(1,ihist) = par(1)
		      dqlfpr(2,ihist) = par(2)
		      dqlfpr(3,ihist) = par(3)
		      dqlfpr(4,ihist) = par(4)
		      dqlfpr(5,ihist) = par(5)
		      dqlfpr(6,ihist) = par(6)
                    else
		      dqlfpr(1,ihist) = 0.
		      dqlfpr(2,ihist) = 0.
		      dqlfpr(3,ihist) = 0.
		      dqlfpr(4,ihist) = 0.
		      dqlfpr(5,ihist) = 0.
		      dqlfpr(6,ihist) = 0.
                    endif
                  else
		    dqlfpr(1,ihist) = 0.
		    dqlfpr(2,ihist) = 0.
		    dqlfpr(3,ihist) = 0.
		    dqlfpr(4,ihist) = 0.
		    dqlfpr(5,ihist) = 0.
		    dqlfpr(6,ihist) = 0.
                  endif
	        ENDDO
		do k=1,nmc
		 if(iyw(k).gt.1) then
                   relec=qtotal(k)
		   dy2 = dyxc2(1,k)
		   dy1 = dyxc1(1,k)
		  qlfy = 2.d0*yfrac(3,k) - 1.d0
		  if(qlfy.gt.0.9) qlfy = 0.9
		  if(qlfy.lt.-0.9) qlfy = -0.9
	      if(relec.ge.rlimit(1)) then
			 dyc=dqlfpr(1,2)+dqlfpr(2,2)*qlfy+dqlfpr(3,2)*qlfy**2
     >          +dqlfpr(4,2)*qlfy**3+dqlfpr(5,2)*qlfy**4
     >          +dqlfpr(6,2)*qlfy**5
			 dy3 = dyx(1,k) - dyc
			 call hfdp1(102,dy3,1.d0)
			 call hfdp1(302,dy2,1.d0)
			 call hfdp1(402,dy1,1.d0)
		  endif
		  if(relec.gt.rlimit(6).and.relec.le.rlimit(1)) then
			 dyc=dqlfpr(1,3)+dqlfpr(2,3)*qlfy+dqlfpr(3,3)*qlfy**2
     >          +dqlfpr(4,3)*qlfy**3+dqlfpr(5,3)*qlfy**4
     >          +dqlfpr(6,3)*qlfy**5
			 dy3 = dyx(1,k) - dyc
			 call hfdp1(103,dy3,1.d0)
			 call hfdp1(303,dy2,1.d0)
			 call hfdp1(403,dy1,1.d0)
		  endif
		  if(relec.gt.rlimit(9).and.relec.le.rlimit(6)) then
			 dyc=dqlfpr(1,4)+dqlfpr(2,4)*qlfy+dqlfpr(3,4)*qlfy**2
     >          +dqlfpr(4,4)*qlfy**3+dqlfpr(5,4)*qlfy**4
     >          +dqlfpr(6,4)*qlfy**5
			 dy3 = dyx(1,k) - dyc
			 call hfdp1(104,dy3,1.d0)
			 call hfdp1(304,dy2,1.d0)
			 call hfdp1(404,dy1,1.d0)
		  endif
		  if(relec.le.rlimit(9)) then
			 dyc=dqlfpr(1,5)+dqlfpr(2,5)*qlfy+dqlfpr(3,5)*qlfy**2
     >          +dqlfpr(4,5)*qlfy**3+dqlfpr(5,5)*qlfy**4
     >          +dqlfpr(6,5)*qlfy**5
			 dy3 = dyx(1,k) - dyc
			 call hfdp1(105,dy3,1.d0)
			 call hfdp1(305,dy2,1.d0)
			 call hfdp1(405,dy1,1.d0)
		  endif
		  call hfdp1(101,dy3,1.d0) 
		  call hfdp1(301,dy2,1.d0)
		  call hfdp1(401,dy1,1.d0)
		 endif
		enddo
C
C  Determine the average offset of single-pixel clusters
C
        nxone=0
	    dxone=0.
		sxone=0.
        nxtwo=0
	    dxtwo=0.
		sxtwo=0.
        do k=1,nmc
C
C  calculate offset from center of pixel
C
          if(ixw(k).eq.1) then	
		    nxone=nxone+1    
	        x0=(ixc(k)-LHxp1)*xsize
			deltax=x0-xhit(k)
	        dxone=dxone+deltax
			sxone=sxone+deltax**2
		  endif
          if(ixw1(k).eq.1) then	
		    nxtwo=nxtwo+1    
	        x0=(ixc1(k)*2-(float(Lx)/2.+1.))*xsize
			deltax=x0-xhit(k)
	        dxtwo=dxtwo+deltax
			sxtwo=sxtwo+deltax**2
		  endif
          if(ixw2(k).eq.1) then	
		    nxtwo=nxtwo+1    
	        x0=(ixc2(k)*2-(float(Lx)/2.))*xsize
			deltax=x0-xhit(k)
	        dxtwo=dxtwo+deltax
			sxtwo=sxtwo+deltax**2
		  endif
	    enddo
	    if(nxone.le.10) then
	      dxone=lorbsx
	      sxone=xsize/root12
	    else
	      dxone=dxone/float(nxone)
	      sxone=sxone/float(nxone)-dxone**2
		  if(sxone.lt.0.) sxone=0.
		  sxone=sqrt(sxone)
	    endif
	    if(nxtwo.le.10) then
	      dxtwo=lorbsx
	      sxtwo=2.*xsize/root12
	    else
	      dxtwo=dxtwo/float(nxtwo)
	      sxtwo=sxtwo/float(nxtwo)-dxtwo**2
		  if(sxtwo.lt.0.) sxtwo=0.
		  sxtwo=sqrt(sxtwo)
	    endif
	    print 610, nxone, dxone, sxone
610     format(1x,i5,' single pixel x-clusters, average offset = ',
     >         f7.1,', rms = ',f7.1)
	    print 615, nxtwo, dxtwo, sxtwo
615     format(1x,i5,' single big pix x-clusters, average offset = ',
     >         f7.1,', rms = ',f7.1)
C
C  Loop over events and determine x position from linear fit
C
        do k=1,nmc
          relec=qtotal(k)
C
C  Keep track of total cluster size in pixels vs charge band
C
          dnpix=npix(k)
	  if(relec.ge.rlimit(1)) then
            call hfdp1(501,dnpix,1.d0) 
          endif
	  if(relec.gt.rlimit(6).and.relec.le.rlimit(1)) then
	    call hfdp1(502,dnpix,1.d0) 
	  endif
          if(relec.gt.rlimit(9).and.relec.le.rlimit(6)) then
	    call hfdp1(503,dnpix,1.d0)
	  endif
	  if(relec.le.rlimit(9)) then
	    call hfdp1(504,dnpix,1.d0)
	  endif
C
C  Get first position estimate from ORCA estimator
C	    
	      x0=(ixc(k)-LHxp1)*xsize
	      xc = x0+float(ixw(k)-1)*(xsize/2.)
		  xwidth=zsize*pvec(4,k)/pvec(6,k)-lorwdx
	  	  if(ixw(k).lt.3) then
		      wzin=0.
		  else
		      wzin=(ixw(k)-2)*xsize
		  endif
	      if(ixw(k).eq.1) then
	         shift=0.
	      elseif(ixw(k).eq.2) then
    	        shift=(xsunc(ixc(k)+ixw(k)-1,k)-xsunc(ixc(k),k))
     >          /(xsunc(ixc(k)+ixw(k)-1,k)+xsunc(ixc(k),k))
     >          /2.*(abs(xwidth)-wzin)
	      else
    	        shift=(xsunc(ixc(k)+ixw(k)-1,k)-xsunc(ixc(k),k))
     >          /(xsunc(ixc(k)+ixw(k)-1,k)+xsunc(ixc(k),k))/2.*xsize
	      endif
	      xrec1=xc+shift-lorbsx
C
C  Shift the x-clusters to put the middle in pixels 12 (or 11 if even)
C
           midpix = ixc(k)+ixw(k)/2
		   ishiftx = LHxp1 - midpix
		   lstpix = ixc(k)+ixw(k)-1
		   if(ishiftx.gt.0) then
		     do i=lstpix,ixc(k),-1
			   xsum(i+ishiftx,k) = xsum(i,k)
			   xsum(i,k) = 0.
			  enddo
		   elseif(ishiftx.lt.0) then
		     do i=ixc(k), lstpix
			   xsum(i+ishiftx,k) = xsum(i,k)
			   xsum(i,k)=0.
			  enddo
		   endif
		   ixc(k)=ixc(k)+ishiftx
C
C  Calculate the min first-pass chi**2 in all cases (not just cluster size > 1)
C
C  Take disributions from first estimate
C
C  Zero linear fit counters
C	     
		    chimnx=1.d15
		    lmin=0
		    sxthr=1.1*sxmax
		    nsort=0
		    do i=1,Lx
		      if(xsum(i,k).gt.0.) then
		        nsort=nsort+1
		        qsort(nsort)=xsum(i,k)
		      endif
		    enddo
		    call sortzv(qsort,isort,nsort,1,0,0)
			if(ixw(k).gt.1) then
		      if(qsort(isort(2)).gt.sxthr) sxthr=1.01*qsort(isort(2))
			else
			  sxthr=1.01*qsort(isort(1))
			endif
C	
C  add 1/2 threshold signals to end pixels
C   
		    xsum(ixc(k)-2,k)=q10
		    xsum(ixc(k)-1,k)=q10
		    xsum(ixc(k)+ixw(k),k)=q10
		    xsum(ixc(k)+ixw(k)+1,k)=q10
            do l=3,39
	          ss2=0.
	          ssa=0.
		      sa2=0.
		      do i=1,Lx
		        if(xsum(i,k).gt.0.) then
			      xsig=min(xsum(i,k), sxmax)
			      xsig=max(xsig, q50)
			      if(i.le.LHxp1) then
			        sigi2=xpar(1,1)+xpar(2,1)*xsig+xpar(3,1)*xsig**2
     >                   +xpar(4,1)*xsig**3+xpar(5,1)*xsig**4
			      else
			        sigi2=xpar(1,2)+xpar(2,2)*xsig+xpar(3,2)*xsig**2
     >                   +xpar(4,2)*xsig**3+xpar(5,2)*xsig**4
			      endif
				  if(sigi2.lt.0..and.xsig.gt.q10) then
				     print 8888, i, sigi2, xsig
8888       format(1x,'i = ',i2,' sigx2 = ',e11.4,' xsignal = ',e11.4)
		     print 8886, 1,xpar(1,1),xpar(2,1),xpar(3,1),
     >               xpar(4,1),xpar(5,1)
		     print 8886, 2,xpar(1,2),xpar(2,2),xpar(3,2),
     >               xpar(4,2),xpar(5,2)
C
C  If tstsig is set true, quit on this error
C
                     if(tstsig) then
                       print 8889, ifile
8889                   format(1x,'Critical error in run ',i5)
                       go to 55555
					 else
					   sigi2 = q50**2
					 endif
				  endif
			      if(i.le.(ixc(k)-1).or.i.ge.(ixc(k)+ixw(k))) 
     >               sigi2=q50**2
			      if(xsum(i,k).gt.sxthr) sigi2=(10000.)**2
			      if(sigi2.le.0.) sigi2 = q50**2
		          ss2=ss2+xsum(i,k)**2/sigi2
			      ssa=ssa+xsum(i,k)*ztemp(i,l)/sigi2
			      sa2=sa2+ztemp(i,l)**2/sigi2
		        endif
	          enddo
		      rat=ssa/ss2
		      chi2=ss2-2.*ssa/rat+sa2/rat**2
			  if(chi2.lt.0.) then
				print 7777, chi2, rat, ss2, ssa, sa2
7779            format(1x,'negative chi2x:',5(2x,e11.4))
			  endif
		      if(chi2.lt.chimnx) then
		        chimnx=chi2
			    lmin=l
				 ss2min=ss2
				 ssamin=ssa
				 sa2min=sa2
		      endif
		    enddo
			if(chimnx.lt.0.) chimnx = 0.
			if(k.eq.1) then
			  ss20=ss2min
			  ssa0=ssamin
			  sa20=sa2min
			else
			  ss2=ss20+ss2min
			  ssa=ssa0+ssamin
			  sa2=sa20+sa2min
			  rnorm1=ssa/ss2
			  chi2=ss2-2.*ssa/rnorm1+sa2/rnorm1**2
	        if(qmerge(k).ge.2.*rlimit(1)) then
		       call hfminx(510,chi2,1.d0) 
		     endif
		     if(qmerge(k).gt.2.*rlimit(6)
     >        .and.qmerge(k).le.2.*rlimit(1)) then
		       call hfminx(511,chi2,1.d0) 
		     endif
		     if(qmerge(k).gt.2.*rlimit(9)
     >        .and.qmerge(k).le.2.*rlimit(6)) then
		       call hfminx(512,chi2,1.d0) 
		     endif
		     if(qmerge(k).le.2.*rlimit(9)) then
		       call hfminx(513,chi2,1.d0) 
		     endif
			  ss20=ss2min
			  ssa0=ssamin
			  sa20=sa2min
			endif
		  if(ixw(k).gt.1) then
		    xrec2=(0.125*lmin-2.625-ishiftx)*xsize
C
C  Make a second. refined pass
C
C
C  Make a second. refined pass
C
		    ibinl=lmin-1
		    ibinh=ibinl+2
		    if(ibinl.lt.1) ibinl=1
		    if(ibinh.gt.41) ibinh=41
C
C  Zero linear fit counters
C	    
	        ss2=0.
	        ssa=0.
            sa2=0.
            ssba=0.
	        saba=0.
	        sba2=0. 
C
C  Take disributions from first estimate
C	
		    do i=1,Lx
		      if(xsum(i,k).gt.0.) then
			    xsig=min(xsum(i,k), sxmax)
			    xsig=max(xsig, q50)
			    if(i.le.LHxp1) then
			      sigi2=xpar(1,1)+xpar(2,1)*xsig+xpar(3,1)*xsig**2
     >                 +xpar(4,1)*xsig**3+xpar(5,1)*xsig**4
			    else
			      sigi2=xpar(1,2)+xpar(2,2)*xsig+xpar(3,2)*xsig**2
     >                 +xpar(4,2)*xsig**3+xpar(5,2)*xsig**4
			    endif
				if(i.le.(ixc(k)-1).or.i.ge.(ixc(k)+ixw(k))) 
     >            sigi2=q50**2
			    if(xsum(i,k).gt.sxthr) sigi2=(10000.)**2
			    if(sigi2.le.0.) sigi2 = q50**2
		        ss2=ss2+xsum(i,k)**2/sigi2
			    ssa=ssa+xsum(i,k)*ztemp(i,ibinl)/sigi2
			    sa2=sa2+ztemp(i,ibinl)**2/sigi2
			    ssba=ssba+xsum(i,k)*(ztemp(i,ibinh)-ztemp(i,ibinl))
     >               /sigi2
			    saba=saba+ztemp(i,ibinl)*(ztemp(i,ibinh)-ztemp(i,ibinl)
     >                    )/sigi2
			    sba2=sba2+(ztemp(i,ibinh)-ztemp(i,ibinl))**2/sigi2
		      endif
	        enddo
		    rat3=(ssba*ssa-ss2*saba)/(ss2*sba2-ssba**2)
		    if(rat3.lt.0.) rat3=0.
		    if(rat3.gt.1.) rat3=1.0
		    rnorm3=(ssa+rat3*ssba)/ss2
		    xrec3=(0.125*ibinl-2.625+rat3*(ibinh-ibinl)*0.125-ishiftx)
     >             *xsize
			chi2x=ss2-2./rnorm3*ssa-2.*rat3/rnorm3*ssba
     >            +(sa2+2.*rat3*saba+rat3**2*sba2)/rnorm3**2
	        if(chi2x.lt.0.) chi2x=0.
		  else
		    xrec2=(ixc(k)-LHxp1-ishiftx)*xsize-dxone
	            xrec3=xrec2
		    chi2x=chimnx
		  endif
		  dx1=xrec1-xhit(k)
		  dx2=xrec2-xhit(k)
		  dx3=xrec3-xhit(k)
		  if(ixw(k).eq.1) then
			call hfminx(206,chi2x,1.d0) 
		  else
	            if(relec.ge.rlimit(1)) then
		       call hfminx(207,chi2x,1.d0) 
		    endif
		    if(relec.gt.rlimit(6).and.relec.le.rlimit(1)) then
		       call hfminx(208,chi2x,1.d0) 
		    endif
		    if(relec.gt.rlimit(9).and.relec.le.rlimit(6)) then
		       call hfminx(209,chi2x,1.d0) 
		    endif
		    if(relec.le.rlimit(9)) then
		       call hfminx(210,chi2x,1.d0) 
		    endif
		  endif
		  dyx(2,k) = dx3
                  dyxc2(2,k) = dx2
                  dyxc1(2,k) = dx1
		  qlfx = 2.d0*xfrac(3,k) - 1.d0
		  dfx = dx3
		  if(ixw(k).gt.1) then 
			call hfill(211,qlfx,dfx,1.)
		    if(relec.gt.rlimit(1)) then
			  call hfill(212,qlfx,dfx,1.)
		    endif
		    if(relec.gt.rlimit(6).and.relec.le.rlimit(1)) then
			  call hfill(213,qlfx,dfx,1.)
		    endif
		    if(relec.gt.rlimit(9).and.relec.le.rlimit(6)) then
			  call hfill(214,qlfx,dfx,1.)
		    endif
		    if(relec.le.rlimit(9)) then
			  call hfill(215,qlfx,dfx,1.)
			endif
		  endif
		enddo  
		DO I=11,15
		  ID=200+I
		  PAR(1)=0.d0
		  STEP(1)=0.d0
		  PAR(2)=HMAX(ID)/2.d0
		  STEP(2)=0.1*PAR(2)
		  PAR(3)=0.d0
		  STEP(3)=0.d0
		  PAR(4)=0.d0
		  STEP(4)=0.01d0
		  PAR(5)=0.d0
		  STEP(5)=0.d0
		  PAR(6)=0.d0
		  STEP(6)=0.01d0
 		  ihist=ihist+1
		  idhist(ihist)=ID
                  CALL HNOENT(ID, NOENT)
                  IF(NOENT.gt.50) THEN
		    CALL HFITHN(ID,CHPOL,CHOPOL,6,PAR,STEP,
     >                   PMIN,PMAX,SIGPAR,CHI2H)
                    if(CHI2H.lt.1.d6) then
		      dqlfpr(1,ihist) = par(1)
		      dqlfpr(2,ihist) = par(2)
		      dqlfpr(3,ihist) = par(3)
		      dqlfpr(4,ihist) = par(4)
		      dqlfpr(5,ihist) = par(5)
		      dqlfpr(6,ihist) = par(6)
                    else
		      dqlfpr(1,ihist) = 0.
		      dqlfpr(2,ihist) = 0.
		      dqlfpr(3,ihist) = 0.
		      dqlfpr(4,ihist) = 0.
		      dqlfpr(5,ihist) = 0.
		      dqlfpr(6,ihist) = 0.
                    endif
                  else
		    dqlfpr(1,ihist) = 0.
		    dqlfpr(2,ihist) = 0.
		    dqlfpr(3,ihist) = 0.
		    dqlfpr(4,ihist) = 0.
		    dqlfpr(5,ihist) = 0.
		    dqlfpr(6,ihist) = 0.
                  endif
	        ENDDO
		do k=1,nmc
		 if(ixw(k).gt.1) then
                  relec=qtotal(k)
                  dx2 = dyxc2(2,k)
                  dx1 = dyxc1(2,k)
		  qlfx = 2.d0*xfrac(3,k) - 1.d0
		  if(qlfx.gt.0.9) qlfx = 0.9
		  if(qlfx.lt.-0.9) qlfx = -0.9
	      if(relec.ge.rlimit(1)) then
			 dxc=dqlfpr(1,7)+dqlfpr(2,7)*qlfx+dqlfpr(3,7)*qlfx**2
     >          +dqlfpr(4,7)*qlfx**3+dqlfpr(5,7)*qlfx**4
     >          +dqlfpr(6,7)*qlfx**5
			 dx3 = dyx(2,k) - dxc
			 call hfdp1(202,dx3,1.d0)
			 call hfdp1(307,dx2,1.d0)
			 call hfdp1(407,dx1,1.d0)
		  endif
		  if(relec.gt.rlimit(6).and.relec.le.rlimit(1)) then
			 dxc=dqlfpr(1,8)+dqlfpr(2,8)*qlfx+dqlfpr(3,8)*qlfx**2
     >          +dqlfpr(4,8)*qlfx**3+dqlfpr(5,8)*qlfx**4
     >          +dqlfpr(6,8)*qlfx**5
			 dx3 = dyx(2,k) - dxc
			 call hfdp1(203,dx3,1.d0)
			 call hfdp1(308,dx2,1.d0)
			 call hfdp1(408,dx1,1.d0)
		  endif
		  if(relec.gt.rlimit(9).and.relec.le.rlimit(6)) then
			 dxc=dqlfpr(1,9)+dqlfpr(2,9)*qlfx+dqlfpr(3,9)*qlfx**2
     >          +dqlfpr(4,9)*qlfx**3+dqlfpr(5,9)*qlfx**4
     >          +dqlfpr(6,9)*qlfx**5
			 dx3 = dyx(2,k) - dxc
			 call hfdp1(204,dx3,1.d0)
			 call hfdp1(309,dx2,1.d0)
			 call hfdp1(409,dx1,1.d0)
		  endif
		  if(relec.le.rlimit(9)) then
			 dxc=dqlfpr(1,10)+dqlfpr(2,10)*qlfx+dqlfpr(3,10)*qlfx**2
     >          +dqlfpr(4,10)*qlfx**3+dqlfpr(5,10)*qlfx**4
     >          +dqlfpr(6,10)*qlfx**5
			 dx3 = dyx(2,k) - dxc
			 call hfdp1(205,dx3,1.d0)
			 call hfdp1(310,dx2,1.d0)
			 call hfdp1(410,dx1,1.d0)
		  endif
		  call hfdp1(201,dx3,1.d0) 
		  call hfdp1(306,dx2,1.d0)
		  call hfdp1(406,dx1,1.d0)
		 endif
		enddo
C
C  Output everything
C
        print 1700, header
1700    format(1x,A80)
        print 1800, n
1800    format(1x,'number of events processed = ',i9)
        print 1805, lorwdx, lorwdy, lorbsx, lorbsy
1805    format(1x,'x/y lorentz widths = ',f6.1,'/',f6.1,
     >            ', x/y lorentz biases = ',f6.1,'/',f6.1)
        call hprintm(ifile,idhist,xshist)
C
C  Reformat and print template information
C
C  First the global z (local y) template
C
	    write(21,1890) ifile, -cosy, -cosx, -cosz
1890    format(1x,i5,3(1x,f9.6))
	    write(21,1900) rnelec,pixmax,sxmax,-dxone,sxone,symax,
     >                     -dyone,syone
1900    format(8(1x,f8.1))
	    write(21,1910) -dxtwo,sxtwo,-dytwo,sytwo,qmsort(imsort(30)), 
     >        clslnx,clslny
C     >        clslnx,clslny,xshist(2,41),xshist(3,41),xshist(4,41)
1910    format(5(1x,f8.1),2(1x,f8.2),2(1x,f8.1),1x,f8.6)
	    do k=2,1,-1
		  write(21,164) (xpar(i,k), i=1,nprm)
	    enddo
	    do k=9,1,-1
		  write(21,170) (xtemp(i,k),i=Nx,1,-1)
	    enddo
	    do k=4,3,-1
		  write(21,164) (xpar(i,k), i=1,nprm)
	    enddo
	    do k=9,1,-1
		  write(21,170) (ytemp(i,k),i=Ny,1,-1)
	    enddo
      do k=7,10
        if(xshist(2,k).lt.5.) xshist(2,k)=5.
        if(xshist(4,k).lt.5.) xshist(4,k)=5.
        write(21,2000) -xshist(1,k),xshist(2,k),-xshist(3,k),xshist(4,k)
2000      format(4(1x,f9.1))
      enddo
      do k=7,10
        write(21,164) -dqlfpr(1,k),-dqlfpr(2,k),-dqlfpr(3,k),
     >                -dqlfpr(4,k),-dqlfpr(5,k),-dqlfpr(6,k)
      enddo
      do k=2,5
        if(xshist(2,k).lt.3.) xshist(2,k)=3.
        if(xshist(4,k).lt.3.) xshist(4,k)=3.
        write(21,2000) -xshist(1,k),xshist(2,k),-xshist(3,k),xshist(4,k)
      enddo
      do k=2,5
        write(21,164) -dqlfpr(1,k),-dqlfpr(2,k),-dqlfpr(3,k),
     >                -dqlfpr(4,k),-dqlfpr(5,k),-dqlfpr(6,k)
       enddo
C
C  Add the average Chi**2 and rms' as functions of q to the template
C
      do k=12,15
	write(21,2100) xshist(1,k+5),xshist(2,k+5),xshist(1,k),
     >                 xshist(2,k)
2100      format(4(1x,f9.3))
      enddo
C
C  Replace the gaussian fit info for the first pass errors with merged cluster chi2 info
C
      do k=27,30
        write(21,2150) -xshist(1,k),xshist(2,k),
     >                 xshist(1,k+24),xshist(2,k+24)
2150      format(2(1x,f9.1),2(1x,f9.3))
      enddo
      do k=22,25
        write(21,2150) -xshist(1,k),xshist(2,k),
     >                 xshist(1,k+25),xshist(2,k+25)
      enddo
      do k=37,40
        write(21,2000) -xshist(1,k),xshist(2,k),-xshist(3,k),
     >                 xshist(4,k)
      enddo
      do k=32,35
        write(21,2000) -xshist(1,k),xshist(2,k),-xshist(3,k),
     >                 xshist(4,k)
      enddo
C
C  Add 20 spare words to each template entry
C	
C
C  Calculate and store the qbin fractions
C
      tqbin = float(nqbin(1)+nqbin(2)+nqbin(3)+nqbin(4))
      do k=1,4
        fqbin(k) = float(nqbin(k))/tqbin
      enddo	
      fyone = float(nyone)/float(nmc)
      fxone = float(nxone)/float(nmc)
      fytwo = float(nytwo)/float(2*nmc)
      fxtwo = float(nxtwo)/float(2*nmc)
      write(21,2200) xshist(1,16),xshist(2,16),
     >    xshist(1,11),xshist(2,11),qmsort(imsort(1)),
     >    xshist(2,41),xshist(3,41),xshist(4,41),
     >    xshist(1,55),1.
2200  format(4(1x,f9.3),3(1x,f8.1),1x,f8.5,1x,f8.5,1x,f8.1)
      write(21,2201) xshist(2,46),xshist(3,46),xshist(4,46),
     >  (fqbin(k), k=1,3),fxone,fyone,fxtwo,fytwo
2201  format(2(1x,f9.1),8(1x,f8.5))
C
C  Reformat and print generic information
C
C  First the global z (local y) template
C
      write(22,1890) ifile, -cosy, -cosx, -cosz
      write(22,1900) rnelec,pixmax,-dxone,sxone,
     >                     -dyone,syone
      write(22,2210) -dxtwo,sxtwo,-dytwo,sytwo,qmsort(imsort(30)),
     > qmsort(imsort(1))
2210  format(6(1x,f8.1))
      do k=37,40
        if(xshist(2,k).lt.5.) xshist(2,k)=5.
        if(xshist(2,k-5).lt.3.) xshist(2,k-5)=3.
                write(22,2000) -xshist(1,k),xshist(2,k),-xshist(1,k-5),
     >                 xshist(2,k-5)
      enddo
C
C  End the loop over runs
C
	  enddo
	  close(21)
      close(22)
      stop
11111 print 22222
22222 format(1x,'Error opening input file')
      stop
33333 print 44444, ifile, n
44444 format(1x,'error reading run ',i6, 'event ',i6)
      stop
55555 print 66666
66666 format(1x,'negative error squared, stop now')
      stop
      end
      SUBROUTINE HBOOKM(TITLE)
C
C *******************************
C  This routine books histograms
C *******************************
C
      implicit double precision (a-h,o-z)
      real*4 halfxs, halfys, chixmx, chiymx
      character*80 TITLE
      CALL HTITLE(TITLE)
C      halfxs=xsize/2.d0
C      nx=halfxs
      halfxs=300.
      nx=120
C      halfys=ysize/2.d0
C      ny=halfys
      halfys=300.
      ny=120
      ncx=300
	  chiymx=48.
	  chixmx=150.
      CALL HBOOK1(100,'Number of generated e',150,0.,500000.,0.)
      CALL HBOOK1(101,'dy_temp (all sig)',nx,-halfxs,halfxs,0.)
      CALL HBOOK1(102,'dy_temp (signal @> 1.5mn)',nx,-halfxs,halfxs,0.)      
      CALL HBOOK1(103,'dy_temp (1.5mn @> signal @> 1.0mn)',nx,-halfxs,
     > halfxs,0.)      
      CALL HBOOK1(104,'dy_temp (1.0mn @> signal @> 0.85mn)',nx,-halfxs,     
     > halfxs,0.)      
      CALL HBOOK1(105,'dy_temp (0.85mn @> signal)',nx,-halfxs,halfxs,0.)      
      CALL HBOOK1(106,'chi2y_temp (single pix)',nx,0.,chiymx,0.)
      CALL HBOOK1(107,'chi2y_temp (signal @> 1.5mn)',nx,0.,chiymx,0.)      
      CALL HBOOK1(108,'chi2y_temp (1.5mn @> signal @> 1.0mn)',nx,0.,
     > chiymx,0.)      
      CALL HBOOK1(109,'chi2y_temp (1.0mn @> signal @> 0.85mn)',nx,0.,     
     > chiymx,0.)      
      CALL HBOOK1(110,'chi2y_temp (0.85mn @> signal)',nx,0.,chiymx,0.)      
      CALL HBPROF(111,'dy_vs_qlf (all sig)',10,-1.,1.,-50.,50.,' ')
      CALL HBPROF(112,'dy_vs_qlf (signal @> 1.5mn)',10,-1.,1.,
     > -50.,50.,' ')
      CALL HBPROF(113,'dy_vs_qlf (1.5mn @> signal @> 1.0mn)',10,-1.,1.,
     > -50.,50.,' ')      
      CALL HBPROF(114,'dy_vs_qlf (1.0mn @> signal @> 0.85mn)',10,-1.,1.,
     > -50.,50.,' ')      
      CALL HBPROF(115,'dy_vs_qlf (0.85mn @> signal)',10,-1.,1.,
     > -50.,50.,' ')      
      CALL HBOOK1(201,'dx_temp (all sig)',nx,-halfxs,halfxs,0.)
      CALL HBOOK1(202,'dx_temp (signal @> 1.5mn)',nx,-halfxs,halfxs,0.)      
      CALL HBOOK1(203,'dx_temp (1.5mn @> signal @> 1.0mn)',nx,-halfxs,
     > halfxs,0.)      
      CALL HBOOK1(204,'dx_temp (1.0mn @> signal @> 0.85mn)',nx,-halfxs,     
     > halfxs,0.)      
      CALL HBOOK1(205,'dx_temp (0.85mn @> signal)',nx,-halfxs,halfxs,0.)      
      CALL HBOOK1(206,'chi2x_temp (single pix)',nx,0.,chiymx,0.)
      CALL HBOOK1(207,'chi2x_temp (signal @> 1.5mn)',nx,0.,chiymx,0.)      
      CALL HBOOK1(208,'chi2x_temp (1.5mn @> signal @> 1.0mn)',nx,0.,
     > chiymx,0.)
      CALL HBOOK1(209,'chi2x_temp (1.0mn @> signal @> 0.85mn)',nx,0.,     
     > chiymx,0.)      
      CALL HBOOK1(210,'chi2x_temp (0.85mn @> signal)',nx,0.,chiymx,0.)      
      CALL HBPROF(211,'dx_vs_qlf (all sig)',10,-1.,1.,-50.,50.,' ')
      CALL HBPROF(212,'dx_vs_qlf (signal @> 1.5mn)',10,-1.,1.,
     > -50.,50.,' ')
      CALL HBPROF(213,'dx_vs_qlf (1.5mn @> signal @> 1.0mn)',10,-1.,1.,
     > -50.,50.,' ')      
      CALL HBPROF(214,'dx_vs_qlf (1.0mn @> signal @> 0.85mn)',10,-1.,1.,
     > -50.,50.,' ')      
      CALL HBPROF(215,'dx_vs_qlf (0.85mn @> signal)',10,-1.,1.,
     > -50.,50.,' ')      
      CALL HBOOK1(301,'dyc2_temp (all sig)',nx,-halfxs,halfxs,0.)
      CALL HBOOK1(302,'dyc2_temp (signal @> 1.5mn)',nx,-halfxs,halfxs,
     > 0.)      
      CALL HBOOK1(303,'dyc2_temp (1.5mn @> signal @> 1.0mn)',nx,-halfxs,
     > halfxs,0.)      
      CALL HBOOK1(304,'dyc2_temp (1.0mn @> signal @> 0.85mn)',nx,
     > -halfxs,halfxs,0.)      
      CALL HBOOK1(305,'dyc2_temp (0.85mn @> signal)',nx,-halfxs,halfxs,
     > 0.)      
      CALL HBOOK1(306,'dxc2_temp (all sig)',nx,-halfxs,halfxs,0.)
      CALL HBOOK1(307,'dxc2_temp (signal @> 1.5mn)',nx,-halfxs,halfxs,
     > 0.)      
      CALL HBOOK1(308,'dxc2_temp (1.5mn @> signal @> 1.0mn)',nx,-halfxs,
     > halfxs,0.)      
      CALL HBOOK1(309,'dxc2_temp (1.0mn @> signal @> 0.85mn)',nx,
     > -halfxs,halfxs,0.)      
      CALL HBOOK1(310,'dxc2_temp (0.85mn @> signal)',nx,-halfxs,halfxs,
     > 0.)      
      CALL HBOOK1(401,'dyc1_std (all sig)',nx,-halfxs,halfxs,0.)
      CALL HBOOK1(402,'dyc1_std (signal @> 1.5mn)',nx,-halfxs,halfxs,
     > 0.)      
      CALL HBOOK1(403,'dyc1_std (1.5mn @> signal @> 1.0mn)',nx,-halfxs,
     > halfxs,0.)      
      CALL HBOOK1(404,'dyc1_std (1.0mn @> signal @> 0.85mn)',nx,
     > -halfxs,halfxs,0.)      
      CALL HBOOK1(405,'dyc1_std (0.85mn @> signal)',nx,-halfxs,halfxs,
     > 0.)      
      CALL HBOOK1(406,'dxc1_std (all sig)',nx,-halfxs,halfxs,0.)
      CALL HBOOK1(407,'dxc1_std (signal @> 1.5mn)',nx,-halfxs,halfxs,
     > 0.)      
      CALL HBOOK1(408,'dxc1_std (1.5mn @> signal @> 1.0mn)',nx,-halfxs,
     > halfxs,0.)      
      CALL HBOOK1(409,'dxc1_std (1.0mn @> signal @> 0.85mn)',nx,
     > -halfxs,halfxs,0.)      
      CALL HBOOK1(410,'dxc1_std (0.85mn @> signal)',nx,-halfxs,halfxs,
     > 0.)      
      CALL HBOOK1(500,'Cluster Charge',250,0.,500000.,0.)
      CALL HBOOK1(501,'npix(signal @> 1.5mn)',40,0.5,40.5,0.)
      CALL HBOOK1(502,'npix(1.5mn @> signal @> 1.0mn)',40,0.5,40.5,0.)
      CALL HBOOK1(503,'npix(1.0mn @> signal @> 0.85mn)',40,0.5,40.5,0.)
      CALL HBOOK1(504,'npix(0.85mn @> signal)',40,0.5,40.5,0.)
      CALL HBOOK1(505,'2 Cluster Merged Charge',500,0.,1000000.,0.)
      CALL HBOOK1(506,'chi2y_2cls (signal @> 1.5mn)',nx,0.,chiymx,0.)      
      CALL HBOOK1(507,'chi2y_2cls (1.5mn @> signal @> 1.0mn)',nx,0.,
     > chiymx,0.)      
      CALL HBOOK1(508,'chi2y_2cls (1.0mn @> signal @> 0.85mn)',nx,0.,     
     > chiymx,0.)      
      CALL HBOOK1(509,'chi2y_2cls (0.85mn @> signal)',nx,0.,chiymx,0.)      
      CALL HBOOK1(510,'chi2x_2cls (signal @> 1.5mn)',nx,0.,chiymx,0.)      
      CALL HBOOK1(511,'chi2x_2cls (1.5mn @> signal @> 1.0mn)',nx,0.,
     > chiymx,0.)
      CALL HBOOK1(512,'chi2x_2cls (1.0mn @> signal @> 0.85mn)',nx,0.,     
     > chiymx,0.)      
      CALL HBOOK1(513,'chi2x_2cls (0.85mn @> signal)',nx,0.,chiymx,0.)      
      CALL HBOOK1(606,'measured Q/generated Q',300,0.,1.5,0.)    
C
      do I=100,110
        CALL HIDOPT(I,'STAT')
		if(i.gt.105) call hidopt(i,'PFUN')
      enddo
      do I=111,115
        CALL HIDOPT(I,'PROE')
      enddo
      do I=201,210
        CALL HIDOPT(I,'STAT')
		if(i.gt.205) call hidopt(i,'PFUN')
      enddo
      do I=211,215
        CALL HIDOPT(I,'PROE')
      enddo
      do I=301,310
        CALL HIDOPT(I,'STAT')
      enddo
      do I=401,410
        CALL HIDOPT(I,'STAT')
      enddo
      do I=500,505
        CALL HIDOPT(I,'STAT')
      enddo
      do I=506,513
        CALL HIDOPT(I,'STAT')
		  call hidopt(i,'PFUN')
      enddo
      do I=606,606
        CALL HIDOPT(I,'STAT')
      enddo      
      RETURN
      END
      SUBROUTINE HPRINTM(nfile,idhist,xshist)
C
C **************************************
C * This routine prints the histograms *
C **************************************
C
      EXTERNAL ch2fcn, ufit
      real*4 ufit
      PARAMETER (NOPT=2, NP=5)
	  common /ch2par/ dmean, dtot, dshift
	  double precision dmean, dtot, dshift, xmean
      common /HCFITD/dpar(24),fitfcn
      real*8 dpar, fitfcn
      integer nfile, num, idhist(54), nx, ny, nwt, loc
	  real*4 xshist(9,54), x, xdummy, xmi, xma, ymi, yma
C
C  Keep track of minimum chi**2 values
C
      common /minx/ xmin(600)
	  double precision xmin
      character*23 fname
	  character*80 htitle
      dimension par(NP), step(NP), pmin(NP), pmax(NP), sigpar(NP)
      CHARACTER*4 CHOPT(NOPT), FCHOPT, CHFUN, choice
      CHARACTER*4 GCHOPT
      DATA CHOPT/'STA ','FIT '/, FCHOPT/'FQW '/
      DATA GCHOPT/'FUBV'/
	  DATA CHFUN/'G   '/
	  DATA choice/'HIST'/
C
C  Open a PS meta-file 
C
        write(fname,100) nfile
100     format('template_histos',I5.5,'.ps')
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
C  Fit double gaussians to each dx histogram (ihist = 1,10)
C
      ihist=0
      DO J=1,2
        DO I=1,5
          ID=100*J+I
          PAR(1)=HMAX(ID)
          STEP(1)=0.1*PAR(1)
          PAR(2) = hstati(id,1,choice,num)
          PAR(3) = hstati(id,2,choice,num)
          STEP(2)=0.05d0*PAR(3)
          STEP(3)=0.2*PAR(3)
C          PAR(4)=0.05d0*PAR(1)
C          STEP(4)=0.1*PAR(4)
C          PAR(5)=0.d0
C          STEP(5)=1.d0
C          PAR(6)=100.d0
C          STEP(6)=0.1*PAR(6)
 	        ihist=ihist+1
          idhist(ihist)=ID
          CALL HNOENT(ID, NOENT)
          IF(NOENT.gt.50) THEN
            CALL HFITHN(ID,CHFUN,FCHOPT,3,PAR,STEP,
     >                PMIN,PMAX,SIGPAR,CHI2)
	         xshist(1,ihist) = hstati(id,1,choice,num)
	         xshist(2,ihist) = hstati(id,2,choice,num)
            if(abs(par(2)).lt.300..and.abs(par(3)).lt.175.) then
	            xshist(3,ihist) = par(2)
               xshist(4,ihist) = abs(par(3))
            else
	           xshist(3,ihist) = xshist(1,ihist)
              xshist(4,ihist) = xshist(2,ihist)
            endif
          else
	         xshist(1,ihist) = hstati(id,1,choice,num)
				xshist(2,ihist) = hstati(id,2,choice,num)
	         xshist(3,ihist) = xshist(1,ihist)
            xshist(4,ihist) = xshist(2,ihist)
          endif
        ENDDO
	  ENDDO
C
C  Call ch2fcn to establish it's type for f2c use as an external
C
	  dtot = 1.0d0
	  dmean = 0.3d0
	  dshift = 0.d0
	  x = 0.1
	  xdummy = ch2fcn(x)
C
C  Plot the theoretical chi**2 function for the average value of each plot (ihist = 11,20)
C
      DO J=1,2
        DO I=6,10
          ID=100*J+I
		  ihist=ihist+1
		  idhist(ihist)=ID
		  total=hsum(id)
		  if(total.gt.20.) then
		    xmean = hstati(id,1,choice,num)
		    dshift = xmin(id)
		    if(dshift.lt.0.d0) dshift=0.d0
		    if(dshift.gt.xmean) dshift=0.d0
		    dmean = xmean-dshift
		    if(dmean.lt.0.1d0) dmean = xmean
		    xshist(1,ihist) = dmean
		    xshist(2,ihist) = dshift
		    call hgive(id,htitle,nx,xmi,xma,ny,ymi,yma,nwt,loc)
		    dtot=(xma-xmi)/float(nx)*total		  
		    call hfunc(id, ch2fcn)
		  else
		    xshist(1,ihist) = 0.10
		    xshist(2,ihist) = 0.16
		  endif    
        ENDDO
	  ENDDO
C  (ihist = 21,30)
        DO I=1,10
          ID=300+I
          PAR(1)=HMAX(ID)
          STEP(1)=0.1*PAR(1)
          PAR(2) = hstati(id,1,choice,num)
          PAR(3) = hstati(id,2,choice,num)
          STEP(2)=0.05d0*PAR(3)
          STEP(3)=0.2d0*PAR(3)
C          PAR(4)=0.05d0*PAR(1)
C          STEP(4)=0.1*PAR(4)
C          PAR(5)=0.d0
C          STEP(5)=1.d0
C          PAR(6)=100.d0
C          STEP(6)=0.1*PAR(6)
 	  ihist=ihist+1
          idhist(ihist)=ID
          CALL HNOENT(ID, NOENT)
          IF(NOENT.gt.50) THEN
            CALL HFITHN(ID,CHFUN,FCHOPT,3,PAR,STEP,
     >                PMIN,PMAX,SIGPAR,CHI2)
	    xshist(1,ihist) = hstati(id,1,choice,num)
	    xshist(2,ihist) = hstati(id,2,choice,num)
            if(abs(par(2)).lt.300..and.abs(par(3)).lt.175.) then
	       xshist(3,ihist) = par(2)
               xshist(4,ihist) = abs(par(3))
            else
	      xshist(3,ihist) = xshist(1,ihist)
              xshist(4,ihist) = xshist(2,ihist)
            endif
          else
	    xshist(1,ihist) = hstati(id,1,choice,num)
	    xshist(2,ihist) = hstati(id,2,choice,num)
	    xshist(3,ihist) = xshist(1,ihist)
            xshist(4,ihist) = xshist(2,ihist)
          endif
        ENDDO
C  (ihist = 31,40)
        DO I=1,10
          ID=400+I
          PAR(1)=HMAX(ID)
          STEP(1)=0.1*PAR(1)
          PAR(2) = hstati(id,1,choice,num)
          PAR(3) = hstati(id,2,choice,num)
          STEP(2)=0.05d0*PAR(3)
          STEP(3)=0.2d0*PAR(3)
C          PAR(4)=0.05d0*PAR(1)
C          STEP(4)=0.1*PAR(4)
C          PAR(5)=0.d0
C          STEP(5)=1.d0
C          PAR(6)=100.d0
C          STEP(6)=0.1*PAR(6)
 	       ihist=ihist+1
          idhist(ihist)=ID
          CALL HNOENT(ID, NOENT)
          IF(NOENT.gt.50) THEN
            CALL HFITHN(ID,CHFUN,FCHOPT,3,PAR,STEP,
     >                PMIN,PMAX,SIGPAR,CHI2)
	         xshist(1,ihist) = hstati(id,1,choice,num)
	         xshist(2,ihist) = hstati(id,2,choice,num)
            if(abs(par(2)).lt.300..and.abs(par(3)).lt.175.) then
	            xshist(3,ihist) = par(2)
               xshist(4,ihist) = abs(par(3))
            else
	           xshist(3,ihist) = xshist(1,ihist)
              xshist(4,ihist) = xshist(2,ihist)
            endif
          else
	         xshist(1,ihist) = hstati(id,1,choice,num)
	         xshist(2,ihist) = hstati(id,2,choice,num)
	         xshist(3,ihist) = xshist(1,ihist)
            xshist(4,ihist) = xshist(2,ihist)
          endif
        ENDDO
C
C  Call ufit to establish it's type for f2c use as an external
C
	     dpar(1) = 1.0d0
	     dpar(2) = 0.0d0
	     dpar(3) = 1.d0
		  dpar(4) = 0.1d0
	     x = 1.0
	     xdummy = ufit(x)
        print 7777, xdummy
7777    format(1x,'xdummy = ',e11.4)
C ihist = 41
 	     ihist=ihist+1
        idhist(ihist)=500
        CALL HNOENT(500, NOENT)
        IF(NOENT.gt.50) THEN
          par(2)=hstati(500,1,choice,num)
          par(1)=NOENT*20000./par(2)
          par(3)=0.1*par(2)
          par(4)=0.02*par(2)/20000.
          step(1) = 0.2*par(1)
          step(2) = 0.1*par(2)
          step(3) = par(3)*0.2
          step(4) = 0.2*par(4)
          pmin(1) = 0.
          pmin(2) = 0.
          pmin(3) = 0.
          pmin(4) = 0.01
          pmax(1) = 0.
          pmax(2) = 0.
          pmax(3) = 0.
          pmax(4) = 10.0
          CALL HFITH(500,ufit,GCHOPT,4,PAR,STEP,
     >                PMIN,PMAX,SIGPAR,CHI2)
          if(par(2).lt.0..or.par(3).lt.0.) then
            par(2)=0.5*hstati(500,1,choice,num)
            par(1)=NOENT*12500./par(2)
            par(3)=0.20*par(2)
            par(4)=0.01*par(2)/20000.
            step(1) = 0.2*par(1)
            step(2) = 0.1*par(2)
            step(3) = par(3)*0.2
            step(4) = 0.2*par(4)
            pmin(1) = 0.
            pmin(2) = 0.
            pmin(3) = 0.
            pmin(4) = 0.01
            pmax(1) = 0.
            pmax(2) = 0.
            pmax(3) = 0.
            pmax(4) = 10.
            CALL HFITH(500,ufit,GCHOPT,4,PAR,STEP,
     >           PMIN,PMAX,SIGPAR,CHI2)	       
           endif
	       xshist(1,ihist) = hstati(500,1,choice,num)
	       xshist(2,ihist) = par(2)
	       xshist(3,ihist) = par(3)
          xshist(4,ihist) = par(4)
        else
	       xshist(1,ihist) = hstati(500,1,choice,num)
	       xshist(2,ihist) = 0.
	       xshist(3,ihist) = 0.
          xshist(4,ihist) = 0.
        endif
C ihist = 42-45
        DO I=1,4
          ID=500+I
 	       ihist=ihist+1
          idhist(ihist)=ID
          CALL HNOENT(ID, NOENT)
          IF(NOENT.gt.50) THEN
            xshist(1,ihist) = hstati(id,1,choice,num)
	         xshist(2,ihist) = hstati(id,2,choice,num)
	         xshist(3,ihist) = 0.
	         xshist(4,ihist) = 0.
          else
	         xshist(1,ihist) = 0.
				xshist(2,ihist) = 0.
	         xshist(3,ihist) = 0.
	         xshist(4,ihist) = 0.
          endif
        ENDDO
C ihist = 46
 	     ihist=ihist+1
        idhist(ihist)=505
        CALL HNOENT(505, NOENT)
        IF(NOENT.gt.50) THEN
          par(2)=hstati(505,1,choice,num)
          par(1)=NOENT*20000./par(2)
          par(3)=0.1*par(2)
          par(4)=0.02*par(2)/20000.
          step(1) = 0.2*par(1)
          step(2) = 0.1*par(2)
          step(3) = par(3)*0.2
          step(4) = 0.2*par(4)
          pmin(1) = 0.
          pmin(2) = 0.
          pmin(3) = 0.
          pmin(4) = 0.01
          pmax(1) = 0.
          pmax(2) = 0.
          pmax(3) = 0.
          pmax(4) = 10.0
          CALL HFITH(505,ufit,GCHOPT,4,PAR,STEP,
     >                PMIN,PMAX,SIGPAR,CHI2)
          if(par(2).lt.0..or.par(3).lt.0.) then
            par(2)=0.4*hstati(500,1,choice,num)
            par(1)=NOENT*12500./par(2)
            par(3)=0.20*par(2)
            par(4)=0.01*par(2)/20000.
            step(1) = 0.2*par(1)
            step(2) = 0.1*par(2)
            step(3) = par(3)*0.2
            step(4) = 0.2*par(4)
            pmin(1) = 0.
            pmin(2) = 0.
            pmin(3) = 0.
            pmin(4) = 0.01
            pmax(1) = 0.
            pmax(2) = 0.
            pmax(3) = 0.
            pmax(4) = 10.
            CALL HFITH(505,ufit,GCHOPT,4,PAR,STEP,
     >                PMIN,PMAX,SIGPAR,CHI2)	       
           endif
	       xshist(1,ihist) = hstati(505,1,choice,num)
	       xshist(2,ihist) = par(2)
	       xshist(3,ihist) = par(3)
          xshist(4,ihist) = par(4)
        else
	       xshist(1,ihist) = hstati(505,1,choice,num)
	       xshist(2,ihist) = 0.
	       xshist(3,ihist) = 0.
          xshist(4,ihist) = 0.
        endif
C
C  Plot the theoretical chi**2 function for the average value of each plot
C
C ihist = 47,54
        DO id=506,513
 		    ihist=ihist+1
		    idhist(ihist)=ID
		    total=hsum(id)
		    if(total.gt.20.) then
		      xmean = hstati(id,1,choice,num)
		      dshift = xmin(id)
		      if(dshift.lt.0.d0) dshift=0.d0
		      if(dshift.gt.xmean) dshift=0.d0
		      dmean = xmean-dshift
		      if(dmean.lt.0.1d0) dmean = xmean
		      xshist(1,ihist) = dmean
		      xshist(2,ihist) = dshift
		      call hgive(id,htitle,nx,xmi,xma,ny,ymi,yma,nwt,loc)
		      dtot=(xma-xmi)/float(nx)*total		  
		      call hfunc(id, ch2fcn)
		    else
		      xshist(1,ihist) = 0.10
		      xshist(2,ihist) = 0.16
		    endif    
        ENDDO
C
C  ihist=55
C
      DO I=6,6
        ID=600+I
          ihist=ihist+1
          idhist(ihist)=ID
          CALL HNOENT(ID, NOENT)
          IF(NOENT.gt.50) THEN
            xshist(1,ihist) = hstati(id,1,choice,num)
            xshist(2,ihist) = hstati(id,2,choice,num)
            xshist(3,ihist) = 0.
            xshist(4,ihist) = 0.
          else
            xshist(1,ihist) = 0.
            xshist(2,ihist) = 0.
            xshist(3,ihist) = 0.
            xshist(4,ihist) = 0.
          endif
      ENDDO
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
      RETURN
      END
      SUBROUTINE HFMINX(ID,XVAL,WVAL)
C
C *****************************************************************
C * This routine is a double to single precision bridge for HBOOK *
C * Keep track of minimum value of x                              *
C *****************************************************************
C
      common /minx/ xmin(600)
      REAL*8 XVAL,WVAL, xmin
      REAL*4 X,Y,W
      DATA Y/0./
	  if(id.gt.0.and.id.le.600) then
	     if(xval.lt.xmin(id)) xmin(id)=xval
	  endif
      X=XVAL
      W=WVAL
      CALL HFILL(ID,X,Y,W)
      RETURN
      END
      SUBROUTINE HFDP1(ID,XVAL,WVAL)
C
C *****************************************************************
C * This routine is a double to single precision bridge for HBOOK *
C *****************************************************************
C
      REAL*8 XVAL,WVAL
      REAL*4 X,Y,W
      DATA Y/0./
      X=XVAL
      W=WVAL
      CALL HFILL(ID,X,Y,W)
      RETURN
      END
      SUBROUTINE HFDP2(ID,XVAL,YVAL,WVAL)
C
C *****************************************************************
C * THIS ROUTINE IS A DOUBLE TO SINGLE PRECISION BRIDGE FOR HBOOK *
C *****************************************************************
C
      REAL*8 XVAL,YVAL,WVAL
      REAL*4 X,Y,W
      X=XVAL
      Y=YVAL
      W=WVAL
      CALL HFILL(ID,X,Y,W)
      RETURN
      END
C
      SUBROUTINE TRIPLG(X)
C
C *****************************************************************
C * This routine calculates 21 gaussianly-distributed random nums *
C * Parameters:  X(13) - 21 random numbers (rms = 1.)             *
C *****************************************************************
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
      DO I=1,21
        X(I)=RBUFF(IBASE+I)
      ENDDO
      IBASE=IBASE+21
      RETURN
      END
C
C
C
	  real function ch2fcn(x)
	  common /ch2par/ dmean, dtot, dshift
	  real*4 x, ch2
	  double precision dmean, dtot, hmean, dshift, xs
	  hmean = dmean/2.d0
	  if(hmean.le.0.) hmean=0.01
	  if(hmean.gt.50.) hmean=50.
	  xs=x-dshift
	  if(xs.le.0.001d0) xs=0.001d0
	  ch2 = dtot*(0.5d0**hmean)*(xs**(hmean-1.d0))
     >    *exp(-xs/2.d0)/dgamma(hmean)
	  ch2fcn=ch2
	  return
	  end

      real function ufit(x)
      real*4 x, kappa, beta2, xl, xu
C
C  Landau function to fit to the cluster charge distribution
C
      common /HCFITD/dpar(24),fitfcn
      real*8 dpar, fitfcn
      real*4 arg
      data beta2/1.0/
      kappa = dpar(4)
      arg = (x-dpar(2))/dpar(3)
      xl = arg - 0.001
      xu = arg + 0.001
      call vviset(kappa, beta2, 0, xl, xu)
      fitfcn = dpar(1)*vviden(arg)
      ufit=fitfcn
      return
      end
