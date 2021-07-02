//! \file template_code6d.cpp
//!
//! Template hit reconstruction algorithms, add FPix to templates, add double pixels to templates, 
//! change dp calling seq (Add PSI46v2 response function and allow for zero-suppressed ROC output)
//! Change to Root historgraams for long-term compatibility
//! Add angle vs resolution for templates and "Standard Algorithms"
//! Tune CMSSW simulation for template 1 reconstruction
//! Change standard algorithm to always use edge method for y-reconstruction
//! Add Estar template number 4
//! Do cosmics
//! Add clustering algorithm
//! Update to include new cluster container
//! Update to include Phase 1 FPix

#include "PixelGeneric2D.h"



int SiPixelTemplateReco::PixelGeneric2D(int ID, float cotalpha, float cotbeta, float locBz, float locBx,  ClusMatrix & cluster, SiPixelGenError& gtempl,
				float& yrec, float& sigmay, float& xrec, float& sigmax, int& nypix, int& nxpix, float& yfrac, float& xfrac, bool IrradiationBiasCorrection = true)			
{
    // Local variables 
	int i, j;
	int nclusx, nclusy;
	float ypix, ypitch, fy, ly, yc, fypitch, lypitch, ytotal, loryloc, lorybias;
	float syone, dyone, sytwo, dytwo, sxone, dxone, sxtwo, dxtwo, pixmx, delx, dely, sigx, sigy, sigma, delta;
	float xpix, xpitch, fx, lx, xc, fxpitch, lxpitch, xtotal, lorxloc, lorxbias, effwidth;
//	const float ysize={150.}, xsize={100.}, zsize[2]={285.,270.}, lorwidth[4] = {0., 18.3, 145.2, 38.3};
	float xsize, ysize, zsize;
//	const float lorwidth[10][4] = {0., 13.5, 82.3, 34.0, 0., 10.21, 113.73, 25.05, 0., 0., 0., 0.,0., 8.40, 79.5, 21.2, 0., 13.5, 122.0, 33.7, 0., 8.93, 114.9, 22.45, 0., 3.40, 113.26, 17.92,
//	                               0., 6.31, 105.81, 19.38,0., 9.96, 134.00, 29.10,0., 8.40, 46.0, 21.2}; //Default
//	const float lorwidth[9][4] = {0., 13.5, 82.3, 34.0, 0., 10.21, 113.73, 25.05, 0., 0., 0., 0.,0., 8.40, 79.5, 21.2, 0., 13.5, 122.0, 33.7, 0., 8.93, 114.9, 22.45, 0., 3.40, 113.26, 17.92,0., 6.31, 105.81, 19.38,0., 9.96, 134.00, 29.10}; //Default
//	const float lorwidth[7][4] = {0., 13.5, 82.3, 34.0, 0., 16.55, 114.1, 28.45, 0., 0., 0., 0.,0., 8.40, 75.3, 21.2, 0., 13.5, 121.0, 34.0, 0., 8.93, 114.9, 22.45, 0., 3.40, 113.26, 17.92}; //100V bias, T=283K
//	const float lorwidth[7][4] = {0., 13.5, 82.3, 34.0, 0., 10.20, 114.1, 25.05, 0., 0., 0., 0.,0., 8.40, 75.3, 21.2, 0., 13.5, 121.0, 34.0, 0., 8.93, 114.9, 22.45, 0., 3.40, 113.26, 17.92}; //150V bias, T=283K
//	const float lorwidth[7][4] = {0., 13.5, 82.3, 34.0, 0., 5.90, 114.1, 20.95, 0., 0., 0., 0.,0., 8.40, 75.3, 21.2, 0., 13.5, 121.0, 34.0, 0., 8.93, 114.9, 22.45, 0., 3.40, 113.26, 17.92}; //200V bias, T=283K
//	const float lorwidth[7][4] = {0., 13.5, 82.3, 34.0, 0., 0.50, 114.1, 14.80, 0., 0., 0., 0.,0., 8.40, 75.3, 21.2, 0., 13.5, 121.0, 34.0, 0., 8.93, 114.9, 22.45, 0., 3.40, 113.26, 17.92}; //300V bias, T=283K
	
//	assert((ID >= 0 && ID < 5) || (ID >= 10 && ID < 17));
//	if(ID < 5) {idloc = ID; 
//	} else {
//		if(ID == 12 || ID == 14) {idloc = 2;} else {if(ID == 13) {idloc = 7;} else {if(ID == 15) {idloc = 8;} else {idloc = ID - 5;}}}
//	}
//	if(ID == 16) {idloc = 9;}
	    
// Get charge scaling factor

	
// Check that the cluster container is (up to) a 7x21 matrix and matches the dimensions of the double pixel flags

// Check that the cluster container is (up to) a 7x21 matrix and matches the dimensions of the double pixel flags

	nclusx = cluster.mrow;
	nclusy = (int)cluster.mcol;
	auto const xdouble = cluster.xdouble;
	auto const ydouble = cluster.ydouble;
	
// First, compute total charge

    float qtotal = 0.f;
	for(i=0; i<nclusx*nclusy; ++i) {
       qtotal +=cluster.matrix[i];
	}

   
	// Get LA and error info from the template object
	
	gtempl.qbin(ID, cotalpha, cotbeta, locBz, locBx, qtotal, pixmx, sigy, dely, sigx, delx, syone, dyone, sytwo, dytwo, sxone, dxone, sxtwo, dxtwo);

    if(!IrradiationBiasCorrection){
        dely = delx = 0.;
    }
	
	xsize = gtempl.xsize();
	ysize = gtempl.ysize();
	zsize = gtempl.zsize();
	

   lorxloc = gtempl.lorxwidth();
   loryloc = gtempl.lorywidth();
   lorxbias = gtempl.lorxbias();
   lorybias = gtempl.lorybias();


	
   
//   printf("\n ************ \n");
//   printf("xsize, ysize, zsize = %f, %f, %f \n", xsize, ysize, zsize);
//   printf("qtotal, pixmax, sigy, dely, sigx, delx = %f, %f, %f, %f, %f, %f \n", qtotal, pixmx, sigy, dely, sigx, delx);
//   printf("loryloc, lorybias, lorxloc, lorybias = %f, %f, %f, %f \n", loryloc, lorybias, lorxloc, lorxbias);
   
	
// Next, apply truncation and make the y-projection of the cluster       
	
// Next, sum the total charge and "decapitate" big pixels         

	   
// Next, apply truncation and make the y-projection of the cluster       
	
	std::vector<float> ysum(TYSIZE, 0.);
    for(i=0; i<nclusy; ++i) {
		for(j=0; j<nclusx; ++j) {
			if(cluster(j,i) > pixmx) {cluster(j,i) = pixmx;}
			ysum[i] += cluster(j,i);
		}
	}
	
//	printf("sigx/delx/sxone/dxone/sxtwo/dxtwo = %f/%f/%f/%f/%f/%f\n", sigx,delx,sxone,dxone,sxtwo,dxtwo);
		
// Next, make x-projection of the cluster and copy the double pixel flags into an 11 element container         

	std::vector<float> xsum(TXSIZE, 0.);
    for(j=0; j<nclusx; ++j) {
	   for(i=0; i<nclusy; ++i) {
		  xsum[j] += cluster(j,i);
	   }
    }
    

	        
// next, identify the y-cluster ends, count total pixels, nypix, and logical pixels, logypx   

    int fypix=-1;
	nypix=0;
	int lypix=0;
	fy = 0.; ly = 0.; ypix = 0.;
	fypitch = 0.; lypitch = 0.;
	for(i=0; i<TYSIZE; ++i) {
	   if(ydouble[i]) {ypitch = 2.*ysize;} else {ypitch = ysize;}
	   if(i == 0) {ypix = ypitch/2.;} else {ypix+=ypitch;}
	   if(ysum[i] > 0.) {
	      if(fypix == -1) {fypix = i; fy=ypix; fypitch=ypitch;}
		  ++nypix;
		  lypix = i;
		  ly = ypix-ypitch;
		  lypitch = ypitch;
		}
	}
	
    //Turn off strict length checking (allow gaps)
    /*
	if((lypix-fypix+1) != nypix) { 
	   if (theVerboseLevel > 1) {
         std::cout <<
           "ysum[0-9] = " << ysum[0] << ", " << ysum[1] << ", " << ysum[2] << ", " << ysum[3] << ", " << ysum[4] << ", "
		                  << ysum[5] << ", " << ysum[6] << ", " << ysum[7] << ", " << ysum[8] << ", " << ysum[9] << std::endl;
          std::cout <<
           "ysum[10-19] = " << ysum[10] << ", " << ysum[11] << ", " << ysum[12] << ", " << ysum[13] << ", " << ysum[14] << ", "
		                  << ysum[15] << ", " << ysum[16] << ", " << ysum[17] << ", " << ysum[18] << ", " << ysum[19] << std::endl;
          std::cout <<
           "ysum[20] = " << ysum[20] << std::endl;
       }
	
	   return 1; 
	}
    */
        	        
// next, identify the x-cluster ends, count total pixels, nxpix, and logical pixels, logxpx   

    int fxpix=-1;
	nxpix=0;
	int lxpix=0;
	fx = 0.; lx = 0.; xpix = 0.;
	fxpitch = 0.; lxpitch = 0.;
	for(i=0; i<TXSIZE; ++i) {
	   if(xdouble[i]) {xpitch = 2.*xsize;} else {xpitch = xsize;}
	   if(i == 0) {xpix = xpitch/2.;} else {xpix+=xpitch;}
	   if(xsum[i] > 0.) {
	      if(fxpix == -1) {fxpix = i; fx=xpix; fxpitch=xpitch;}
		  ++nxpix;
		  lxpix = i;
		  lx = xpix-xpitch;
		  lxpitch=xpitch;
		}
	}
    //Turn off strict length checking (allow gaps)
    /*
	if((lxpix-fxpix+1) != nxpix) { 
	
	   if (theVerboseLevel > 1) {
          std::cout <<
           "xsum[0-10] = " << xsum[0] << ", " << xsum[1] << ", " << xsum[2] << ", " << xsum[3] << ", " << xsum[4] << ", "
		                  << xsum[5] << ", " << xsum[6] << std::endl;
       }

	   return 2; 
	}
    */
                
	if (theVerboseLevel > 9) {
       std::cout <<
         " cot(alpha) = " << cotalpha << " cot(beta) = " << cotbeta << 
         " nclusx = " << nclusx << " nclusy = " << nclusy << std::endl;
       std::cout <<
         " cot(alpha) = " << cotalpha << " cot(beta) = " << cotbeta << 
         " nclusx = " << nclusx << " nclusy = " << nclusy << std::endl;
    }
		
// Do the y-reconstruction first, calculate the geometrical center of the cluster
	
	yc = (fy+ly)/2.;
	ytotal = fypitch+lypitch;
	
	if(nypix == 1) {
	
//  Get error and offset for this cluster

	   if(ydouble[fypix]) {
	      sigma = sytwo;
		  delta = dytwo;
	   } else {
          sigma = syone;
		  delta = dyone;
	   }
	   
	   if(sigma <= 0.) {
	      sigmay = 43.3;
	   } else {
          sigmay = sigma;
	   }
		   
// If the cluster size is one, use the geometrical center		   
		   
	   yrec = yc - delta;
	   yfrac = -1.;
	   	   
	} else {
	   
// For 4 pix > cluster > 1 pix
      if(nypix < 3) {
	     
		 effwidth = fabs(zsize*cotbeta + loryloc) - (ly-fy);
		 if(effwidth < 0. || effwidth > ytotal) {effwidth = ytotal/2.;}
		 yrec = yc+(ysum[lypix]-ysum[fypix])/(ysum[lypix]+ysum[fypix])*effwidth/2. - dely - lorybias;
//		 yrec = yc+(ysum[lypix]-ysum[fypix])/(ysum[lypix]+ysum[fypix])*effwidth/2. - lorybias;
		 sigmay = sigy;
		 yfrac = ysum[fypix]/(ysum[lypix]+ysum[fypix]);
		 
	  } else {
	  
	     yrec = yc+(ysum[lypix]-ysum[fypix])/(ysum[lypix]+ysum[fypix])*ytotal/4. - dely - lorybias;
//	     yrec = yc+(ysum[lypix]-ysum[fypix])/(ysum[lypix]+ysum[fypix])*ytotal/4. - lorybias;
		 sigmay = sigy;
		 yfrac = ysum[fypix]/(ysum[lypix]+ysum[fypix]);
	  
	  }	 	 
	}
//	printf("cota/cotb = %f/%f, nypix = %d, dely/lorybias = %f/%f \n", cotalpha, cotbeta, nypix, dely, lorybias);
	
// Do the x-reconstruction next, calculate the geometrical center of the cluster
	
	xc = (fx+lx)/2.;
	xtotal = fxpitch+lxpitch;
	
	if(nxpix == 1) {
		   
//  Get error and offset for this cluster

	   if(xdouble[fxpix]) {
	      sigma = sxtwo;
		  delta = dxtwo;
	   } else {
          sigma = sxone;
		  delta = dxone;
	   }
	   
	   if(sigma <= 0.) {
	      sigmax = 28.9;
	   } else {
          sigmax = sigma;
	   }
	   
// If the cluster size is one, use the geometrical center		   
		   
	   xrec = xc - delta;
	   xfrac = -1.;
	   	   
	} else {
	   
	   effwidth = fabs(zsize*cotalpha + lorxloc) - (lx-fx);
	   if(effwidth < 0. || effwidth > xtotal) {effwidth = xtotal/2.;}
	   xrec = xc+(xsum[lxpix]-xsum[fxpix])/(xsum[lxpix]+xsum[fxpix])*effwidth/2. - delx - lorxbias;
//	   xrec = xc+(xsum[lxpix]-xsum[fxpix])/(xsum[lxpix]+xsum[fxpix])*effwidth/2. - lorxbias;
	   sigmax = sigx;
	   xfrac = xsum[fxpix]/(xsum[lxpix]+xsum[fxpix]);	 
	}
   
//   printf("xrec/sigmax = %f/%f, yrec/sigmay = %f/%f \n", xrec, sigmax, yrec, sigmay);

	
    return 0;
} // PixelGeneric2D 

int SiPixelTemplateReco::PixelGeneric2D(int ID, float cotalpha, float cotbeta, ClusMatrix & cluster, SiPixelGenError& templ,
										float& yrec, float& sigmay, float& xrec, float& sigmax, int& nypix, int& nxpix, float& yfrac, float& xfrac, bool IrradiationBiasCorrection=true)			
{
    // Local variables 
	float locBx, locBz;
	locBx = 1.;
    if(cotbeta < 0.) locBx = -1.;
    locBz = locBx;
    if(cotalpha < 0.) locBz = -locBx;
    
	return SiPixelTemplateReco::PixelGeneric2D(ID, cotalpha, cotbeta, locBz, locBx, cluster, templ, 
												yrec, sigmay, xrec, sigmax, nypix, nxpix, yfrac, xfrac, IrradiationBiasCorrection);
	
} // PixelGeneric2D


