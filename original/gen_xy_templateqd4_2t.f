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
	  external fcn
	  parameter (nev = 30000, npt = 200, nprm = 5, maxarg = 5)
C  Nx and Nx are the dimensions of the x and y arrays
          parameter (Nx = 21, Ny = 13)
c
      dimension xproj(Nx,nev),yproj(Ny,nev),qsum(nev),
     > pixev(Nx,Ny,nev),iproj(nev),jproj(nev),xh(nev),
     > yh(nev),i2d(nev), j2d(nev)
      dimension xprojt(Nx,nev),yprojt(Ny,nev),pixelt(Nx,Ny)
      dimension ier(nprm),xpar(nprm,4), spar(nprm,4)
      dimension bixin(Nx,Ny), pixel(Nx,Ny)
      dimension pvec(6), xsum(Nx), ysum(Ny), jy(Ny)
      dimension xtemp(Nx,9), xtemp2(Nx,9), nxntry(9)  
      dimension ytemp(Ny,9), ytemp2(Ny,9), nyntry(9)  
      dimension xytemp(Nx,Ny,4,4), nxytry(4,4)
      character*30 fname
      character*80 header
      real*8 noise, arglis(maxarg), lorxw1
      real*8 loryw1
      common /fitdat/ signal(npt), ssignl(npt), nfit
      COMMON / UNIT / ISYSRD,ISYSWR,ISYSSA,ISUMMR
      common /fparams/ par(nprm),parerr(nprm), 
     >     emat(nprm,nprm),f0,ipar(nprm),npr,istat,chnam(nprm)
      character*10 chnam
      logical storep, tstsig
      integer NHx, NHy
C
C  Define middle of dimension ranges
C
      NHx = Nx/2+1
      NHy = Ny/2+1 
C
C  Define poldata to process, and read the max error on a run
C
      open(16,file='pix_2t.proc',status='old',err=11111)
      read(16,*) nfile,numrun,noise,q100,q101,fq100,fcal,
     > fgain,rnoise,linear
      close(16)
      print 10, numrun,nfile,q100,q101,fq100,noise,fcal,
     > fgain,rnoise,linear
10    format(1x,'processing ',i3,' files starting at file ',i5,/,
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
C
C  Search for best measured one-pixel offsets
C
      lorxw1 = 0.
      loryw1 = 0.
      nx1max=0
      ny1max=0
      dx1sig=1.d10
      dy1sig=1.d10
      cotamn = 1.d10
      cotbmn = 1.d10
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
                  cota = cosy/cosz
                  cotb = cosx/cosz
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
	    iproj(n)=int(xhit/xsize*8.+5.5)
	    jproj(n)=int(yhit/ysize*8.+5.5)
	    i2d(n)=int(xhit/xsize*4.+3.0)
            if(i2d(n).lt.1) i2d(n) = 1
            if(i2d(n).gt.4) i2d(n) = 4
	    j2d(n)=int(yhit/ysize*4.+3.0)
            if(j2d(n).lt.1) j2d(n) = 1
            if(j2d(n).gt.4) j2d(n) = 4
	    qsum(n)=0.
		do i=1,Nx
		  xsum(i)=0.
                  xprojt(i,n)=0.
		  do j=1,Ny
		    pixel(i,j)=0.
                    pixelt(i,j)=0.
       		  enddo
		    do j=1,Ny
		      qin=bixin(i,j)*rten
		      if(qin.lt.0.d0) qin=0.d0
		      pixel(i,j)=qin
                      pixev(i,j,n)=qin
                      if(qin.gt.q100) pixelt(i,j)=qin
	              xsum(i)=xsum(i)+pixel(i,j)
		      xprojt(i,n)=xprojt(i,n)+pixelt(i,j)
 		    enddo
		  xproj(i,n)=xsum(i)
		  qsum(n)=qsum(n)+xsum(i)			
		enddo
		do j=1,Ny
		  ysum(j) = 0.d0
                  yprojt(j,n) = 0.d0
			do i=1,Nx
			  ysum(j)=ysum(j)+pixel(i,j)
                          yprojt(j,n)=yprojt(j,n)+pixelt(i,j)
			enddo
		  yproj(j,n)=ysum(j)
		enddo
	    qavg=qavg+qsum(n)
        go to 200
500     n=n-1
        nmc=n
        close(15)
	qavg=qavg/float(nmc)
        print 501, ifile, nmc, qavg
501     format(1x,'finish run ',i4,' read ',i5,' events, qavg = ',f9.1)
	    do j=1,9
		  nxntry(j)=0
	 	  nyntry(j)=0
		  do i=1,Nx
		    xtemp(i,j)=0.
		    xtemp2(i,j)=0.
		  enddo
		  do i=1,Ny
		    ytemp(i,j)=0.
		    ytemp2(i,j)=0.
		  enddo
	    enddo
      do k=1,4
        do l=1,4
          nxytry(k,l) = 0
          do i=1,Nx
            do j=1,Ny
	      xytemp(i,j,k,l)=0.
            enddo
          enddo
        enddo
      enddo
      mx1=0
      sx1=0.d0
      sx12=0.d0
      my1=0
      sy1=0.d0
      sy12=0.d0
      do n=1,nmc
C
C  Do single pixel analysis here before charge cut
C
            mx=0
            ix=0
            do i=1,Nx
              if(xprojt(i,n).gt.0.) then
                mx=mx+1
                ix=i
              endif
            enddo
            if(mx.eq.1) then
              xrec=xsize*(ix-NHx)
              dx=xrec-xh(n)
              mx1=mx1+1
              sx1=sx1+dx
              sx12=sx12+dx**2
            endif
            my=0
            do j=1,Ny
              if(yprojt(j,n).gt.0.) then
                my=my+1
                jy(my)=j
              endif
            enddo
            if(my.eq.1) then
              yrec1=ysize*(jy(1)-NHy+0.5)
              yrec2=ysize*(jy(1)-NHy-0.5)
              dy1=yrec1-yh(n)
              dy2=yrec2-yh(n)
              my1=my1+2
              sy1=sy1+dy1+dy2
              sy12=sy12+dy1**2+dy2**2
            endif
            if(my.eq.2) then
              if((jy(2)-jy(1)).eq.1) then
                yrec1=ysize*(jy(1)-NHy)
                yrec2=ysize*(jy(2)-NHy)
                yrec=0.5*(yrec1+yrec2)
                dy=yrec-yh(n)
                my1=my1+1
                sy1=sy1+dy
                sy12=sy12+dy**2
              endif
            endif
	    if(qsum(n).lt.qavg) then
C
C  count the pixel signals for events with less than average charge
C
                k=i2d(n)
                l=j2d(n)
                nxytry(k,l)=nxytry(k,l)+1
		do i=1,Nx
		  do j=1,Ny
		    xytemp(i,j,k,l)=xytemp(i,j,k,l)+pixev(i,j,n)
       		  enddo
                enddo
C
C  Make template shapes
C
		  k=iproj(n)
C
C  Renorm distributions and print them
C
		  if(k.lt.1.or.k.gt.9) then
		    print 505, k
505         format(1x,'problem w/ event = ',i5)
                  endif
		  if(k.eq.1.or.k.eq.9) then
			nxntry(1)=nxntry(1)+1
			nxntry(9)=nxntry(9)+1
		  else
		    nxntry(k)=nxntry(k)+1
		  endif
          do i=1,Nx
			if(k.eq.1) then
			  xtemp(i,1)=xtemp(i,1)+xproj(i,n)
			  xtemp2(i,1)=xtemp2(i,1)+xproj(i,n)**2
			  if(i.gt.1) then
			    xtemp(i,9)=xtemp(i,9)+xproj(i-1,n)
			    xtemp2(i,9)=xtemp2(i,9)+xproj(i-1,n)**2				
			  endif
			elseif(k.eq.9) then
			  xtemp(i,9)=xtemp(i,9)+xproj(i,n)
			  xtemp2(i,9)=xtemp2(i,9)+xproj(i,n)**2
			  if(i.gt.1) then
			    xtemp(i-1,1)=xtemp(i-1,1)+xproj(i,n)
			    xtemp2(i-1,1)=xtemp2(i-1,1)+xproj(i,n)**2
			  endif
			else	
              xtemp(i,k)=xtemp(i,k)+xproj(i,n)
              xtemp2(i,k)=xtemp2(i,k)+xproj(i,n)**2
            endif
          enddo
		  l=jproj(n)
		  if(l.lt.1.or.l.gt.9) then
		    print 505, l
          endif
		  if(l.eq.1.or.l.eq.9) then
			nyntry(1)=nyntry(1)+1
			nyntry(9)=nyntry(9)+1
		  else
		    nyntry(l)=nyntry(l)+1
		  endif
          do j=1,Ny
			if(l.eq.1) then
			  ytemp(j,1)=ytemp(j,1)+yproj(j,n)
			  ytemp2(j,1)=ytemp2(j,1)+yproj(j,n)**2
			  if(j.gt.1) then
			    ytemp(j,9)=ytemp(j,9)+yproj(j-1,n)
			    ytemp2(j,9)=ytemp2(j,9)+yproj(j-1,n)**2				
			  endif
			elseif(l.eq.9) then
			  ytemp(j,9)=ytemp(j,9)+yproj(j,n)
			  ytemp2(j,9)=ytemp2(j,9)+yproj(j,n)**2
			  if(j.gt.1) then
			    ytemp(j-1,1)=ytemp(j-1,1)+yproj(j,n)
			    ytemp2(j-1,1)=ytemp2(j-1,1)+yproj(j,n)**2
			  endif
			else	
              ytemp(j,l)=ytemp(j,l)+yproj(j,n)
              ytemp2(j,l)=ytemp2(j,l)+yproj(j,n)**2
            endif
          enddo
	endif
      enddo
C
C  Find the maximum average pixel signal
C
      pixmax=0.
      do k=1,4
        do l=1,4
          do i=1,Nx
            do j=1,Ny
              sigxy=xytemp(i,j,k,l)/float(nxytry(k,l))
	      if(sigxy.gt.pixmax) pixmax=sigxy
            enddo
          enddo
        enddo
      enddo
      if(mx1.gt.10) then
        sx1=sx1/float(mx1)
        sx12=sx12/float(mx1)
        sx12=sqrt((sx12-sx1**2)/float(mx1))
      endif
      print 600, mx1, sx1, sx12
600   format(1x,'number of 1 x-clusters = ',i5,', dx = ',f7.2,
     >'+-',f7.2)
      if(mx1.gt.10.and.abs(cotb).lt.cotbmn) then
        lorxw1 = sx1
        dx1sig = sx12
        nx1max=mx1
        cotbmn = abs(cotb)
      endif
      if(my1.gt.10) then
        sy1=sy1/float(my1)
        sy12=sy12/float(my1)
        sy12=sqrt((sy12-sy1**2)/float(my1))
      endif
      print 601, my1, sy1, sy12
601   format(1x,'number of 1 y-clusters = ',i5,', dy = ',f7.2,
     >'+-',f7.2)
      if(my1.gt.10.and.abs(cota).le.cotamn) then
        if(my1.gt.ny1max) then
          loryw1 = sy1
          dy1sig = sy12
          ny1max=my1
          cotamn = abs(cota)
        endif
      endif
C
C  Renorm distributions and print them
C
      do k=1,9
	    do i=1,Nx
		  xtemp(i,k)=xtemp(i,k)/float(nxntry(k))
		  xtemp2(i,k)=xtemp2(i,k)/float(nxntry(k))-xtemp(i,k)**2
		  if(xtemp(i,k).lt.0.d0) xtemp(i,k)=0.d0
		enddo
	  enddo
      do l=1,9
	    do j=1,Ny
		  ytemp(j,l)=ytemp(j,l)/float(nyntry(l))
		  ytemp2(j,l)=ytemp2(j,l)/float(nyntry(l))-ytemp(j,l)**2
		  if(ytemp(j,l).lt.0.d0) ytemp(j,l)=0.d0
		enddo
	  enddo
C
C  Analyze error vs signals for entry and exit sides of cluster
C
      ISYSRD=5
      ISYSWR=6
      ISYSSA=7
      call mninit(isysrd,isyswr,isyssa)
	  signal(1)=0.
	  ssignl(1)=0.
      nfit=1
	  sxmax=0.
	  ssxmax=0.
	  do k=1,9
        do i=1,NHx
		  if(xtemp(i,k).gt.20.d0) then
		    if(nfit.lt.npt) nfit=nfit+1
			signal(nfit)=xtemp(i,k)
			ssignl(nfit)=sqrt(xtemp2(i,k))
			if(signal(nfit).gt.sxmax) sxmax=signal(nfit)
			if(ssignl(nfit).gt.ssxmax) ssxmax=ssignl(nfit)
		  endif
		enddo
	  enddo
	  xmax=(int(sxmax/100.)+1)*100.
	  ymax=(int(ssxmax/10.)+1)*10.
	  bguess=ssxmax**2/sxmax
C
C  Open Topdraw file
C
        write(fname,700) ifile
700     format('sigmax_',I5.5,'.top')
        open(18,file=fname,status='unknown',err=11111)
        write(18,710) xmax, ymax
710     format(
     >    1x,'SET FONT DUPLEX',/,
     >    1x,'SET TICKS RIGHT  OFF',/,
     >    1x,'SET TICKS TOP    OFF',/,
     >    1x,'SET LIMITS X FROM 0. TO ',f8.0,
     >       ' Y FROM 0.00 TO ',f8.0,/,
     >    1x,'SET TICKS SIZE 0.05',/,1x,'SET SYMBOL 2O',/,
     >    1x,'SET TITLE SIZE 3.0',/,
     >    1x,'SET ORDER X Y',/,
     >    1x,'PLOT AXES',/,1x,'SET COLOR BLUE') 
	   do j=1,nfit
		  write(18,720) signal(j), ssignl(j)
720       format(1x,f8.1,2x,f8.1)
       enddo
	   write(18,730)
730    format(1x,'PLOT')
c
c     Clear old parameter settings
c
        arglis(1)=0.
        call mnexcm(fcn,'CLEAR',arglis,0,ierflg,dummy)
        call mnseti('Pixel uncertainty fit, x-entrance')
        call mnparm(1,'a',0.d0,1.0d0,0.0d0,0.0d0,ier(1))
        call mnparm(2,'b',bguess,1.0d0,0.0d0,0.0d0,ier(2))
		call mnparm(3,'c',0.d0,1.0d0,0.0d0,0.0d0,ier(3))
		call mnparm(4,'d',0.d0,1.0d0,0.0d0,0.0d0,ier(4))
		call mnparm(5,'e',0.d0,0.0d0,0.0d0,0.0d0,ier(5))
c
c   Do the fit
c
        arglis(1)=2.
        call mnexcm(fcn,'SET STR',arglis,1,ierflg,dummy)
        arglis(1)=0.
        call mnexcm(fcn,'SET NOGRADIENT',arglis,0,ierflg,dummy)
        arglis(1)=10000.
        call mnexcm(fcn,'MIGRAD',arglis,1,ierflg,dummy)
        call mnexcm(fcn,'HESSE',arglis,1,ierflg,dummy)
C
C  Get the fit status
C
        CALL MNSTAT(FMIN,FEDM,ERRDEF,NPARI,NPARX,istat)
        arglis(1)=0.
        call mnexcm(fcn,'RETURN',arglis,0,ierflg,dummy)
        do i=1,nprm
          xpar(i,1)=par(i)
        enddo
C
C  Plot the fit
C
      xstep=xmax/30.
      do i=1,31
         x=(i-1)*xstep
	    write(18,720) x, fit(xpar(1,1),x)
      enddo
	  write(18,750)
750   format(1x,'JOIN SOLID',/,1x,'SET COLOR RED')
	  signal(1)=0.
	  ssignl(1)=0.
      nfit=1
	  do k=1,9
        do i=NHx,Nx
		  if(xtemp(i,k).gt.20.d0) then
		    if(nfit.lt.npt) nfit=nfit+1
		    signal(nfit)=xtemp(i,k)
		    ssignl(nfit)=sqrt(xtemp2(i,k))
		  endif
		enddo
	  enddo
	   do j=1,nfit
		  write(18,720) signal(j), ssignl(j)
       enddo
	   write(18,730)
c
c     Clear old parameter settings
c
        arglis(1)=0.
        call mnexcm(fcn,'CLEAR',arglis,0,ierflg,dummy)
        call mnseti('Pixel uncertainty fit, x-exit')
        call mnparm(1,'a',0.d0,1.0d0,0.0d0,0.0d0,ier(1))
        call mnparm(2,'b',bguess,1.0d0,0.0d0,0.0d0,ier(2))
		call mnparm(3,'c',0.d0,1.0d0,0.0d0,0.0d0,ier(3))
		call mnparm(4,'d',0.d0,1.0d0,0.0d0,0.0d0,ier(4))
		call mnparm(5,'e',0.d0,0.0d0,0.0d0,0.0d0,ier(5))
c
c   Do the fit
c
        arglis(1)=2.
        call mnexcm(fcn,'SET STR',arglis,1,ierflg,dummy)
        arglis(1)=0.
        call mnexcm(fcn,'SET NOGRADIENT',arglis,0,ierflg,dummy)
        arglis(1)=10000.
        call mnexcm(fcn,'MIGRAD',arglis,1,ierflg,dummy)
        call mnexcm(fcn,'HESSE',arglis,1,ierflg,dummy)
C
C  Get the fit status
C
        CALL MNSTAT(FMIN,FEDM,ERRDEF,NPARI,NPARX,istat)
        arglis(1)=0.
        call mnexcm(fcn,'RETURN',arglis,0,ierflg,dummy)
        do i=1,nprm
          xpar(i,2)=par(i)
        enddo
C
C  Plot the fit
C
      xstep=xmax/30.
      do i=1,31
         x=(i-1)*xstep
	    write(18,720) x, fit(xpar(1,2),x)
	  enddo
	  write(18,850) -xmax/10.,ymax/2., header
850   format(1x,'JOIN SOLID',/,1x,'SET COLOR BLUE',/,
     >  1x,'TITLE BOTTOM ''signal''',/,
     >  1x,'TITLE',2(1x,f8.1),' DATA CENTER ANGLE 90 ''sigma''',/,
     >  1x,'TITLE TOP ''',A80,'''')
	  close(18)
	  signal(1)=0.
	  ssignl(1)=0.
      nfit=1
	  symax=0.
	  ssymax=0.
	  do k=1,9
        do i=1,NHy
		  if(ytemp(i,k).gt.20.d0) then
		    if(nfit.lt.npt) nfit=nfit+1
		    signal(nfit)=ytemp(i,k)
		    ssignl(nfit)=sqrt(ytemp2(i,k))
		    if(signal(nfit).gt.symax) symax=signal(nfit)
		    if(ssignl(nfit).gt.ssymax) ssymax=ssignl(nfit)
		  endif
		enddo
	  enddo
	  xmax=(int(symax/100.)+1)*100.
	  ymax=(int(ssymax/10.)+1)*10.
	  bguess=ssymax**2/symax
C
C  Open Topdraw file
C
        write(fname,900) ifile
900     format('sigmay_',I5.5,'.top')
        open(18,file=fname,status='unknown',err=11111)
        write(18,710) xmax, ymax
	   do j=1,nfit
		  write(18,720) signal(j), ssignl(j)
       enddo
	   write(18,730)
c
c     Clear old parameter settings
c
        arglis(1)=0.
        call mnexcm(fcn,'CLEAR',arglis,0,ierflg,dummy)
        call mnseti('Pixel uncertainty fit, y-entrance')
        call mnparm(1,'a',0.d0,1.0d0,0.0d0,0.0d0,ier(1))
        call mnparm(2,'b',bguess,1.0d0,0.0d0,0.0d0,ier(2))
		call mnparm(3,'c',0.d0,1.0d0,0.0d0,0.0d0,ier(3))
		call mnparm(4,'d',0.d0,1.0d0,0.0d0,0.0d0,ier(4))
		call mnparm(5,'e',0.d0,0.0d0,0.0d0,0.0d0,ier(5))
c
c   Do the fit
c
        arglis(1)=2.
        call mnexcm(fcn,'SET STR',arglis,1,ierflg,dummy)
        arglis(1)=0.
        call mnexcm(fcn,'SET NOGRADIENT',arglis,0,ierflg,dummy)
        arglis(1)=10000.
        call mnexcm(fcn,'MIGRAD',arglis,1,ierflg,dummy)
        call mnexcm(fcn,'HESSE',arglis,1,ierflg,dummy)
C
C  Get the fit status
C
        CALL MNSTAT(FMIN,FEDM,ERRDEF,NPARI,NPARX,istat)
        arglis(1)=0.
        call mnexcm(fcn,'RETURN',arglis,0,ierflg,dummy)
        do i=1,nprm
          xpar(i,3)=par(i)
        enddo
C
C  Plot the fit
C
      xstep=xmax/30.
      do i=1,31
         x=(i-1)*xstep
	    write(18,720) x, fit(xpar(1,3),x)
	  enddo
	  write(18,750)
	  signal(1)=0.
	  ssignl(1)=0.
      nfit=1
	  do k=1,9
        do i=NHy,Ny
		  if(ytemp(i,k).gt.20.d0) then
		    if(nfit.lt.npt) nfit=nfit+1
		    signal(nfit)=ytemp(i,k)
		    ssignl(nfit)=sqrt(ytemp2(i,k))
		  endif
		enddo
	  enddo
	   do j=1,nfit
		  write(18,720) signal(j), ssignl(j)
       enddo
	   write(18,730)
c
c     Clear old parameter settings
c
        arglis(1)=0.
        call mnexcm(fcn,'CLEAR',arglis,0,ierflg,dummy)
        call mnseti('Pixel uncertainty fit, y-exit')
        call mnparm(1,'a',0.d0,1.0d0,0.0d0,0.0d0,ier(1))
        call mnparm(2,'b',bguess,1.0d0,0.0d0,0.0d0,ier(2))
		call mnparm(3,'c',0.d0,1.0d0,0.0d0,0.0d0,ier(3))
		call mnparm(4,'d',0.d0,1.0d0,0.0d0,0.0d0,ier(4))
		call mnparm(5,'e',0.d0,0.0d0,0.0d0,0.0d0,ier(5))
c
c   Do the fit
c
        arglis(1)=2.
        call mnexcm(fcn,'SET STR',arglis,1,ierflg,dummy)
        arglis(1)=0.
        call mnexcm(fcn,'SET NOGRADIENT',arglis,0,ierflg,dummy)
        arglis(1)=10000.
        call mnexcm(fcn,'MIGRAD',arglis,1,ierflg,dummy)
        call mnexcm(fcn,'HESSE',arglis,1,ierflg,dummy)
C
C  Get the fit status
C
        CALL MNSTAT(FMIN,FEDM,ERRDEF,NPARI,NPARX,istat)
        arglis(1)=0.
        call mnexcm(fcn,'RETURN',arglis,0,ierflg,dummy)
        do i=1,nprm
          xpar(i,4)=par(i)
        enddo
C
C  Plot the fit
C
      do i=1,31
         x=(i-1)*xstep
	    write(18,720) x, fit(xpar(1,4),x)
	  enddo
	  write(18,850) -xmax/10.,ymax/2., header
	  close(18)
C
C  Open Summary files
C
      write(fname,1000) ifile
1000  format('ztemp_',I5.5,'.txt')
      open(16,file=fname,status='unknown',err=11111)
      write(16,1005) cosx, cosy, cosz
1005  format(3(1x,f9.6))
      if(clslnx.ge.(0.8*xsize).or.clslny.gt.(0.4*ysize)) then
	write(16,1010) qavg, sxmax, pixmax
1010    format(3(1x,f8.1))
        do k=1,2
          write(16,1020) (xpar(i,k), i=1,nprm)
1020      format(5(1x,D15.8))
          do i=1,nprm
	    spar(i,k) = xpar(i,k)
	  enddo
        enddo
        spxmax = sxmax
        storep=.true.
      else
C
C If the cluster length is smaller than a pixel, there aren't enough points for a reliable
C fit.  Use the last (1 pixel length cluster) params instead for the z-template
C
	if(storep) then
	  write(16,1010) qavg, spxmax, pixmax
          do k=1,2
            write(16,1020) (spar(i,k), i=1,nprm)
	  enddo
        else
	  write(16,1010) qavg, sxmax, pixmax
          do k=1,2
            write(16,1020) (xpar(i,k), i=1,nprm)
          enddo
	endif
      endif
      do k=1,9
	write(16,1030) k, (k*0.125-0.625)*xsize
1030    format(1x,'bin ',i2,', xcenter = ',f8.2,' um')
	write(16,1040) (xtemp(i,k),i=1,Nx)
1040    format(23(1x,f8.1))
      enddo
      write(fname,1100) ifile
1100  format('ptemp_',I5.5,'.txt')
      open(17,file=fname,status='unknown',err=11111)
      write(17,1005) cosx, cosy, cosz
      write(17,1010) qavg, symax, pixmax
      do k=3,4
        write(17,1020) (xpar(i,k), i=1,nprm)
      enddo
	  do l=1,9
	    write(17,1110) l, (l*0.125-0.625)*ysize
1110    format(1x,'bin ',i2,', ycenter = ',f8.2,' um')
	    write(17,1040) (ytemp(i,l),i=1,Ny)
	  enddo
	  close(16)
	  close(17)
	  enddo
C	  
C  Open output file for the Lorentz Widths
C
      write(fname,1200) nfile
1200  format('lorentz_widths_xy',I5.5,'.out')
      open(21,file=fname,status='unknown',err=11111)
      write(21,1210) 2.*lorxw1, 2.*dx1sig, nx1max, cotbmn,
     > 2.*loryw1, 2.*dy1sig, ny1max, cotamn
1210    format(2(1x,f9.2,1x,f9.2,1x,i5,1x,f7.4))
      close(21)
      write(fname,1300) nfile
1300  format('lorentz_widths_new',I5.5,'.out')
      open(22,file=fname,status='unknown',err=11111)
      write(22,1310) 2.*lorxw1,lorxw1,2.*loryw1,loryw1
1310    format(4(1x,f8.2))
      close(22)      
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
      SUBROUTINE FCN(NPAR,GPAR,F,XPAR,IFLAG)
C
C
C  ************************************************************
C  * This routine reads the experimental data from a file and *
C  * calculates a chisquared.                                 *
C  ************************************************************
C
C
      parameter (npt = 200, nprm=5, maxarg=5)
      implicit double precision (a-h,o-z)
	  common /fitdat/ signal(npt), ssignl(npt), nfit
      COMMON / UNIT / ISYSRD,ISYSWR,ISYSSA,ISUMMR
      common /fparams/ par(nprm),parerr(nprm), 
     >     emat(nprm,nprm),f0,ipar(nprm),npr,istat,chnam(nprm)
      character*10 chnam
C
C  Dimension the parameter arrays
C
      DIMENSION XPAR(*), GPAR(*)
      REAL*4 CHI2
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
      DO i=1,nfit
	    wgt = 1.d0
		if(signal(i).le.0.1d0) wgt = 10.d0
        CHISQD=CHISQD+wgt*(ssignl(i)-fit(xpar,signal(i)))**2 
      ENDDO
C
C  Load the answer into the parameter list
C
      F=CHISQD
      IF(IFLAG.NE.3) GO TO 9999
C
C  Termination Section...print and save any results
C
      F0=CHISQD
      DO I=1,NPRM
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
      CHI2=CHISQD
9999  RETURN
      END
C
C
      double precision FUNCTION FIT(xpar,signal)
C
C *****************************************************
C * Pixel response function                           *
C * Parameters: xpar(5) - fit parameters              *
C *              signal - input signal size           *
C *****************************************************
C
      implicit double precision (a-h,o-z)
      real*8 xpar(5), signal
C
C  The function
C
	  arg=xpar(1)+xpar(2)*signal+xpar(3)*signal**2
     >   +xpar(4)*signal**3+xpar(5)*signal**4
	  aarg=abs(arg)
	  if(aarg.gt.0.) then
	    fit=arg/aarg*sqrt(aarg)
	  else
	    fit=0.d0
	  endif
      return
      end
