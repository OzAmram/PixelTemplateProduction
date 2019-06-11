
/* *******************************************************************   
 * Original Author: Morris Swartz
 * Translator to C++: Cristina
 *
 *
 * This program reads and aanalyzes 21x13 pixel summary files.     *   
 * Create cluster length templates for x-/y-fitting procedure      *
 * Implement 2-pass appproach (3/6/06)                             *
 * Store templates at Nx x Ny                                       *
 * Assume that data are corrected for ROC response.                *
 * Use input thickness (zsize).                                    *
 * Add (0,0) point to each fit (8 Feb 2007).                       *
 * Fake z error parameters for short clusters if possible (8 Feb)  *
 * Reverse run processing direction to make faking work with cotb  *
 * reflection scheme (cotbeta > 0 only is stored temp) [1 Mar 07]  *
 * Add pixmax information to both output files.       [21 Oct 08]  *
 * Make adjustable array and template sizes.          [21 Oct 08]  *
 * Add digits to the error parameterization.          [22 Jul 09]  *
 * Modify small cluster error parameterization.       [22 Jul 09]  *
 * Add single pixel processing to estimate Lor widths [31 Oct 09]  *
 * Update for new response modeling input files       [28 Jan 10]  *
 *
 * Translated from fortran to C++                      May 2019
 *******************************************************************
 */
#include "template_utils.h"


Double_t fit_fn(Double_t *xs, Double_t *par){
    //sqrt of 4th order polynomial, avoid negative
    Float_t x = xs[0];
    Double_t arg = par[0] + par[1]*x + par[2]*x*x + par[3] * x*x*x + par[4]*x*x*x*x;
    Double_t abs_arg = TMath::Abs(arg);
    if(abs_arg > 0.) return (arg/abs_arg)* TMath::Sqrt(abs_arg);
    else return 0.;
}


void gen_xy_template2d(const int nevents = 30000, const int npt = 200, const int nprm = 5, 
		       const int maxarg = 5){


    //size of the x and y arrays
    const int Nx = 21;
    const int Ny = 13;


    double  pixelt[Nx][Ny],  pixel[Nx][Ny],
    xsum[Nx], ysum[Ny],
    xpar[nprm][4], spar[nprm][4],
    xytemp[Nx][Ny][7][7], xytmp2[Nx][Ny][7][7]; 

    int nxytry[7][7];

    //declare larger arrays dynamically to avoid stack overflow
    //should probably make some of these actually 2d or 3d arrays rather than
    //doing the index gymnastics currently in the code

    double *pixev = new double[Nx*Ny*nevents];
    double *qsum = new double[nevents];
    double *xh = new double[nevents];
    double *yh = new double[nevents];
    int *i2d = new int[nevents];
    int *j2d = new int[nevents];

    bool storep = false;


    //size of pixel (in pixelav coordinates)
    float xsize, ysize, zsize;


    //define middle of dimension ranges
    const int NHx = Nx/2;
    const int NHy = Ny/2;

    FILE *f_config = fopen("pix_2t.proc", "r");

    if(f_config ==0){
        printf("Can't find pix_2t.proc! Exiting \n");
        exit(1);
    }

    int file_start, num_files, non_linear;
    float q100, q101, q100_frac, noise, common_frac, gain_frac, readout_noise;
	fscanf(f_config,"%d %d %f %f %f %f %f %f %f %d", &file_start, &num_files, &noise, &q100, &q101, &q100_frac, &common_frac, 
            &gain_frac, &readout_noise, &non_linear);
	fclose(f_config);
	printf("processing %d files starting from %d, noise = %f, threshold0 = %f, threshold1 = %f," 
            "rms threshold frac = %f, common_frac = %f, gain fraction = %f, readout noise = %f, nonlinear_resp = %d \n", 
            num_files, file_start, noise, q100, q101, q100_frac, common_frac, gain_frac, readout_noise, non_linear);

    fclose(f_config);

    //define RMS of noise
    double rten = 10.;
    double thr = q100;
    double thr10 = 0.1*thr;

    //search for best measured one-pixel offsets
    float lorxw1 = 0.;
    float loryw1 = 0.;
    float dx1sig=1.10;
    float dy1sig=1.10;
    float cotamn=1.10;
    float cotbmn=1.10;



    char fname[100];
    char header[120];
    char buffer[120];

    //loop over files backwards
    for(int iFile = file_start + num_files - 1; iFile >= file_start; iFile--){
        sprintf(fname, "template_events_d%i.out", iFile);
        FILE *f_evts = fopen(fname, "r");

        if(f_evts ==0) {
            printf("Cant find file %s! Exiting \n", fname);
            exit(1);
        }

        fgets(header, 120, f_evts);

        fgets(buffer, 120, f_evts); 

        sscanf(buffer, "%f %f %f", &xsize, &ysize, &zsize);
        printf("xsize = %f,  ysize = %f, zsize = %f \n", xsize, ysize, zsize);


        float z_center = zsize/2.0;

        float cosx, cosy, cosz, x1,y1,z1;
        float clslnx, clslny, cota, cotb;

        double qavg=0.;


        int nelec;

        //loop over events 
        int n;
        for(n=0; n < nevents; n++){
            if(fgets(buffer, 120, f_evts) ==0)
                //end of file
                break;

            sscanf(buffer, "%f %f %f %f %f %f %i",
                    &x1, &y1, &z1, &cosx, &cosy, &cosz, &nelec);

            if(n == 0){
                clslnx = fabs(zsize * cosx/cosz);
                clslny = fabs(zsize * cosy/cosz);
                cota = cosy/cosz;
                cotb = cosx/cosz;
            }



            float bixin[Nx][Ny];
            memset(bixin, 0., sizeof(bixin));
            //note the confusing index switch. input arrays are 13x21, we swap
            //to 21x13
            for(int j=0; j<Ny; j++){
                for (int i=0; i<Nx; i++){
                    fscanf(f_evts, " %f ", &bixin[i][j]);
                }
            }


            //  Pixelav gives hit position at face of pixel, translate to
            //  3d center of the pixel
            //  propagate from edge to center of pixel
            float xhit = x1 + (z_center - z1) * cosx/cosz;
            float yhit = y1 + (z_center - z1) * cosy/cosz;

            xh[n] = xhit;
            yh[n] = yhit;

            //hit position bins (9 bins within a pixel)
            //last bin is the same as first, translated 1 pixel over
            //divide into 2d bins (4x4) for max charge
            int i2dn  = int(xhit/xsize*6 + 3.5);

            //restrict from 0 to 6;
            i2dn = std::min(6, std::max(i2dn, 0));
            i2d[n] = i2dn;

            int j2dn  = int(yhit/ysize*6 + 3.5);
            //restrict from 0 to 6;
            j2dn = std::min(6, std::max(j2dn, 0));
            j2d[n] = j2dn;

            qsum[n] = 0.;
            memset(ysum, 0., sizeof(ysum));
            memset(pixelt, 0., sizeof(pixelt));

            for(int i=0; i<Nx; i++){
                for(int j=0; j<Ny; j++){
                    double qin = bixin[i][j] * rten;
                    if(qin < 0.) qin = 0.;
                    pixev[i*(Ny *nevents) + j*nevents + n] = qin;
		    qsum[n] += qin;
                }
            }
            qavg += qsum[n];

        } // finish loop over events

        int nmc = n;

        fclose(f_evts);
        qavg /= double(nmc);
        printf("finish reading run %i read %i events, qavg = %.1f \n", iFile, nmc, qavg);

        //zero things for this run
        memset(nxytry, 0, sizeof(nxytry));
        memset(xytemp, 0., sizeof(xytemp));
        memset(xytmp2, 0., sizeof(xytmp2));

        int mx1=0;
        float sx1 = 0.;
        float sx12 = 0.;

        int my1=0;
        float sy1 = 0.;
        float sy12 = 0.;

        for(n=0; n<nmc; n++){

            //count number of pixels with less than average charge
            //We make templates with only events with less than average charge
            //to avoid long tails (charge is truncated in reco anyway)
            if(qsum[n] < qavg){

                //fill 2d 7x7 2d template for average charge determination
                int k = i2d[n];
                int l = j2d[n];
                nxytry[k][l] +=1;

		for(int i=0; i < Nx; i++){
		  for(int j=0; j < Ny; j++){
		    xytemp[i][j][k][l] += pixev[i*(Ny *nevents) + j*nevents + n]; // pixev(i,j,n)?                                                                                                     
		    xytmp2[i][j][k][l] += pixev[i*(Ny *nevents) + j*nevents + n]**2;

		    // 0th bin of new pixel same as last bin of previous pixel 
		    // im so confused here...
                    if(k==0 && i>0){
		      xytemp[i][j][6][l] += pixev[(i-1)*(Ny *nevents) + j*nevents + n];
		      xytmp2[i][j][6][l] += pixev[(i-1)*(Ny *nevents) + j*nevents + n]**2;
 		    }
                    else if(k==6 && i>0){
		      xytemp[i-1][j][0][l] += pixev[i*(Ny *nevents) + j*nevents + n];
                      xytmp2[i-1][j][0][l] += pixev[i*(Ny *nevents) + j*nevents + n]**2;
		    }

                    if(l==0 && j>0){
                      xytemp[i][j][k][6] += pixev[i*(Ny *nevents) + (j-1)*nevents + n];
                      xytmp2[i][j][k][6] += pixev[i*(Ny *nevents) + (j-1)*nevents + n]**2;
                    }
                    else if(l==6 && j>0){
                      xytemp[i][j-1][k][0] += pixev[i*(Ny *nevents) + j*nevents + n];
                      xytmp2[i][j-1][k][0] += pixev[i*(Ny *nevents) + j*nevents + n]**2;
                    }

                    if((k==0 && l==0) && (i>0 && j>0)){
		      xytemp[i][j][6][6] += pixev[(i-1)*(Ny *nevents) + (j-1)*nevents + n];
                      xytmp2[i][j][6][6] += pixev[(i-1)*(Ny *nevents) + (j-1)*nevents + n]**2;
		    }
		    else if((k==6 && l==6) && (i>0 && j>0)){
		      xytemp[i-1][j-1][0][0] += pixev[i*(Ny *nevents) + j*nevents + n];
                      xytmp2[i-1][j-1][0][0] += pixev[i*(Ny *nevents) + j*nevents + n]**2;
                    }
                    else if((k==0 && l==6) && (i>0 && j>0)){
                      xytemp[i][j-1][0][0] +=pixev[(i-1)*(Ny *nevents) + j*nevents + n];
                      xytmp2[i][j-1][0][0] +=pixev[(i-1)*(Ny *nevents) + j*nevents + n]**2;
		    }
                    else if((k==6 && l==0) && (i>0 && j>0)){
                      xytemp[i-1][j][0][0] +=pixev[i*(Ny *nevents) + (j-1)*nevents + n];
                      xytmp2[i-1][j][0][0] +=pixev[i*(Ny *nevents) + (j-1)*nevents + n]**2;
                    }

		  }
                }

                if(k < 0 || k > 6){
                    printf("Problem with event %i, index (k) is =%i \n", n, k);
                }
                else{
                    nxytry[k][l] +=1;
                }

                if(l < 0 || l > 8){
                    printf("Problem with event %i, index (l) is =%i \n", n, l);
                }
                else{
		  nxytry[k][l] +=1; //why is it adding again!!!
                }

            } 
        }//end loop over events

        /*
         * print number of template entries (debug)
	 */
        for(int i=0; i<9; i++){
            printf("number of template entries %i ", nxytry[i]);
        }
        printf("\n");

	//renorm the distributions, find the maximum average pixel signal in 7x7 grid
        float pixmax = 0.;
	int imin = Nx;
	int imax = 1;
	int jmin = Ny;
	int jmax = 1;
        for(int i =0; i<Nx; i++){
            for(int j =0; j<Ny; j++){
                for(int k =0; k<7; k++){
                    for(int l =0; l<7; l++){
                        float sigxy = xytemp[i][j][k][l] / float(nxytry[k][l]);
			xytemp[i][j][k][l] = sigxy
                        if(sigxy > pixmax) pixmax = sigxy;
                        xytmp2[i][j][k][l] = xytmp2[i][j][k][l] / float(nxytry[k][l]) - sigxy**2; // is -sigxy**2  correct here?
			if(sigxy > thr10) {
			  if(i < imin) imin=i;
                          if(i > imax) imax=i;
			  if(j < jmin) jmin=j;
                          if(j > jmax) jmax=j;
			}
                    }
                }
            }
        }

        //analyze variance of charge vs charge for entry and exit sides of cluster
        char plot_title[100];

        //x projection fits
        double sxmax = 0.;
        double ssxmax = 0.;
        std::vector<float> xsignal1, xssignal1, xsignal2, xssignal2;
        std::vector<float> xssignal1err, xssignal2err;

        //start with zeros to get better fit
        xsignal1.push_back(0.);
        xssignal1.push_back(0.);
        xsignal2.push_back(0.);
        xssignal2.push_back(0.);

        //add small error to intercept to anchor fit
        float zero_err = 1./pow(10, 0.5);
        xssignal1err.push_back(zero_err);
        xssignal2err.push_back(zero_err);

        for (int k=0; k<9; k++){
            for(int i=0; i<Nx; i++){
                //separate two sides of cluster
                if(xtemp[i][k] > 20.){
                    if(i<= NHx){
                        xsignal1.push_back(xtemp[i][k]);
                        xssignal1.push_back(sqrt(xtemp2[i][k]));
                        xssignal1err.push_back(1.);

                        sxmax = std::max(sxmax, xtemp[i][k]);
                        ssxmax = std::max(ssxmax, sqrt(xtemp2[i][k]));
                    }
                    if(i>=NHx){
                        xsignal2.push_back(xtemp[i][k]);
                        xssignal2.push_back(sqrt(xtemp2[i][k]));
                        xssignal2err.push_back(1.);
                    }

                }
            }
        }
        float guess = ssxmax*ssxmax/sxmax;



        TGraphErrors *g_xtemp1 = new TGraphErrors(xsignal1.size(), xsignal1.data(), xssignal1.data(), nullptr, xssignal1err.data());
        g_xtemp1->SetTitle("X projections: Charge Variance vs. Charge");
        TGraphErrors *g_xtemp2 = new TGraphErrors(xsignal2.size(), xsignal2.data(), xssignal2.data(), nullptr, xssignal2err.data());
        TCanvas *c_xtemp = new TCanvas("c_xtemp", "", 800, 800);


        TF1 *f_xtemp1 = new TF1("f_xtemp1", fit_fn,
                -0.01, sxmax*1.01, nprm);
        f_xtemp1->SetParameter(0, 0.);
        f_xtemp1->SetParameter(1, guess);
        f_xtemp1->SetParameter(2, 0.);
        f_xtemp1->SetParameter(3, 0.);
        f_xtemp1->FixParameter(4, 0.);

        f_xtemp1->SetParError(0, 1.);
        f_xtemp1->SetParError(1, 1.);
        f_xtemp1->SetParError(2, 1.);
        f_xtemp1->SetParError(3, 1.);
        f_xtemp1->SetParError(4, 0.);
        f_xtemp1->SetLineColor(kRed);

        TF1 *f_xtemp2 = new TF1("f_xtemp2", fit_fn,
                -0.01, sxmax*1.01, nprm);
        f_xtemp2->SetParameter(0, 0.);
        f_xtemp2->SetParameter(1, guess);
        f_xtemp2->SetParameter(2, 0.);
        f_xtemp2->SetParameter(3, 0.);
        f_xtemp2->FixParameter(4, 0.);

        f_xtemp2->SetParError(0, 1.);
        f_xtemp2->SetParError(1, 1.);
        f_xtemp2->SetParError(2, 1.);
        f_xtemp2->SetParError(3, 1.);
        f_xtemp2->SetParError(4, 0.);
        f_xtemp2->SetLineColor(kBlue);


        //set minuit options for all fits
        //
        ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
        ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.1);
        ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);

        g_xtemp1->Fit(f_xtemp1, "VR EX0");
        g_xtemp1->Draw("AP");
        g_xtemp1->SetMarkerSize(1.4);
        g_xtemp1->SetMarkerStyle(20);


        g_xtemp2->Fit(f_xtemp2, "VR EX0");
        g_xtemp2->Draw("P same");
        g_xtemp2->SetMarkerSize(1.4);
        g_xtemp2->SetMarkerStyle(21);

        sprintf(plot_title, "sigmax_%i.png", iFile);
        c_xtemp->Print(plot_title);


        //y projection fits
        double symax = 0.;
        double ssymax = 0.;
        std::vector<float> ysignal1, yssignal1, ysignal2, yssignal2;
        std::vector<float>  yssignal1err,  yssignal2err;

        //start with zeros to anchor fit
        ysignal1.push_back(0.);
        yssignal1.push_back(0.);
        ysignal2.push_back(0.);
        yssignal2.push_back(0.);

        //add small error to intercept for stability
        yssignal1err.push_back(zero_err);
        yssignal2err.push_back(zero_err);

        for (int k=0; k<9; k++){
            for(int i=0; i<Ny; i++){
                //separate two sides of sensor
                if(ytemp[i][k] > 20.){
                    if(i<= NHy){
                        ysignal1.push_back(ytemp[i][k]);
                        yssignal1.push_back(sqrt(ytemp2[i][k]));
                        yssignal1err.push_back(1.);

                        symax = std::max(symax, ytemp[i][k]);
                        ssymax = std::max(ssymax, sqrt(ytemp2[i][k]));
                    }
                    if(i>= NHy){
                        ysignal2.push_back(ytemp[i][k]);
                        yssignal2.push_back(sqrt(ytemp2[i][k]));
                        yssignal2err.push_back(1.);
                    }

                }
            }
        }
        guess = ssymax*ssymax/symax;


        TGraph *g_ytemp1 = new TGraphErrors(ysignal1.size(), ysignal1.data(), yssignal1.data(), nullptr, yssignal1err.data());
        g_ytemp1->SetTitle("Y projections: Charge Variance vs. Charge");
        TGraph *g_ytemp2 = new TGraphErrors(ysignal2.size(), ysignal2.data(), yssignal2.data(), nullptr, yssignal2err.data());
        TCanvas *c_ytemp = new TCanvas("c_ytemp", "", 800, 800);

        TF1 *f_ytemp1 = new TF1("f_ytemp1", fit_fn,
                -0.01, symax*1.01, nprm);
        f_ytemp1->SetParameter(0, 0.);
        f_ytemp1->SetParameter(1, guess);
        f_ytemp1->SetParameter(2, 0.);
        f_ytemp1->SetParameter(3, 0.);
        f_ytemp1->FixParameter(4, 0.);

        f_ytemp1->SetParError(0, 1.);
        f_ytemp1->SetParError(1, 1.);
        f_ytemp1->SetParError(2, 1.);
        f_ytemp1->SetParError(3, 1.);
        f_ytemp1->SetParError(4, 0.);
        f_ytemp1->SetLineColor(kRed);

        TF1 *f_ytemp2 = new TF1("f_ytemp2", fit_fn,
                -0.01, symax*1.01, nprm);
        f_ytemp2->SetParameter(0, 0.);
        f_ytemp2->SetParameter(1, guess);
        f_ytemp2->SetParameter(2, 0.);
        f_ytemp2->SetParameter(3, 0.);
        f_ytemp2->FixParameter(4, 0.);

        f_ytemp2->SetParError(0, 1.);
        f_ytemp2->SetParError(1, 1.);
        f_ytemp2->SetParError(2, 1.);
        f_ytemp2->SetParError(3, 1.);
        f_ytemp2->SetParError(4, 0.);
        f_ytemp2->SetLineColor(kBlue);

        g_ytemp1->Fit(f_ytemp1, "VR EX0");
        g_ytemp1->Draw("AP");
        g_ytemp1->SetMarkerSize(1.4);
        g_ytemp1->SetMarkerStyle(20);


        g_ytemp2->Fit(f_ytemp2, "VR EX0");
        g_ytemp2->Draw("P same");
        g_ytemp2->SetMarkerSize(1.4);
        g_ytemp2->SetMarkerStyle(21);

        sprintf(plot_title, "sigmay_%i.png", iFile);
        c_ytemp->Print(plot_title);

        //get fit params
        for(int p=0; p<nprm; p++){
            xpar[p][0] = f_xtemp1->GetParameter(p);
            xpar[p][1] = f_xtemp2->GetParameter(p);
            xpar[p][2] = f_ytemp1->GetParameter(p);
            xpar[p][3] = f_ytemp2->GetParameter(p);
        }

        /*
        //print out charges and variances that are being fit (for debugging)
        printf("y1s \n");
        for(int i=0; i<ysignal1.size(); i++){
            printf("%.1f %.1f %.1f \n", ysignal1[i], yssignal1[i], yssignal1err[i]);
        }
        */
            


        //open summary files
        char file_out[100];

        //x template file
        sprintf(file_out, "zptemp_%i.txt", iFile);
        FILE *ztemp_file = fopen(file_out, "w+");
        fprintf(ztemp_file, "%9.6f %9.6f %9.6f \n", cosx, cosy, cosz);
	qavg,pixmax,imin,imax,jmin,jmax
	  fprintf(ztemp_file, "%8.1f %8.1f %i %i %i %i\n", qavg,pixmax,imin,imax,jmin,jmax); // format i?
	  for(int k=2; k <4; k++){
	    for(int p=0; p<nprm; p++){
	      fprintf(ztemp_file, "%15.8E ", xpar[p][k]);
	      spar[p][k] = xpar[p][k];
	    }
	    fprintf(ztemp_file, "\n");
	  }
	spxmax = sxmax;
	storep = true;

        // If the cluster length is smaller than a pixel, there aren't enough points for a reliable
        // fit.  Use the one from the previous file
        else if(storep){
            fprintf(ztemp_file, "%8.1f %8.1f %8.1f \n", qavg, spxmax, pixmax);
            for(int k=0; k <2; k++){
                for(int p=0; p<nprm; p++){
                    fprintf(ztemp_file,"%15.8E ", spar[p][k]);
                }
                fprintf(ztemp_file, "\n");
            }

        }
        //we are the first file so there are no previous to reference
        else{
            fprintf(ztemp_file, "%8.1f %8.1f %8.1f \n", qavg, sxmax, pixmax);
            for(int k=0; k <2; k++){
                for(int p=0; p<nprm; p++){
                    fprintf(ztemp_file, "%15.8E ", xpar[p][k]);
                }
                fprintf(ztemp_file, "\n");
            }
        }



        for(int k = 0; k<9; k++){
            float xcenter = ((k+1)*0.125 - 0.625) *xsize;
            fprintf(ztemp_file, "bin %2i,  xcenter = %8.2f  um \n", k+1, xcenter);
            for(int i=0; i<Nx; i++){
                fprintf(ztemp_file, "%8.1f ", xtemp[i][k] );
            }
            fprintf(ztemp_file, "\n");
        }

        fclose(ztemp_file);

        //y template file
        sprintf(file_out, "ptemp_%i.txt", iFile);
        FILE *ptemp_file = fopen(file_out, "w+");
        fprintf(ptemp_file, "%9.6f %9.6f %9.6f \n", cosx, cosy, cosz);
        fprintf(ptemp_file, "%8.1f %8.1f %8.1f \n", qavg, symax, pixmax);

        for(int k=0; k <2; k++){
            for(int p=0; p<nprm; p++){
                fprintf(ptemp_file, "%15.8E ", xpar[p][k]);
            }
            fprintf(ptemp_file, "\n");
        }

        for(int l = 0; l<9; l++){
            float ycenter = ((l+1)*0.125 - 0.625) *ysize;
            fprintf(ptemp_file, "bin %2i  ycenter = %8.2f  um \n", l+1, ycenter);
            for(int j=0; j<Ny; j++){
                fprintf(ptemp_file, "%8.1f ", ytemp[j][l] );
            }
            fprintf(ptemp_file, "\n");
        }

        fclose(ptemp_file);


    } //end loop over files


    //output lorentz widths
    char  lorz_out[100];
    sprintf(lorz_out, "lorentz_widths_xy%i.out", file_start);
    FILE *f_lor_width = fopen(lorz_out, "w+");
    fprintf(f_lor_width,  "%9.2f %9.2f %5i %7.4f"
            "%9.2f %9.2f %5i %7.4f",
            2.*lorxw1, 2.*dx1sig, nx1max, cotbmn,
            2.*loryw1, 2.*dy1sig, ny1max, cotamn);
    fclose(f_lor_width);

    sprintf(lorz_out, "lorentz_widths_new%i.out", file_start);
    FILE *f_lor_new = fopen(lorz_out, "w+");
    fprintf(f_lor_new, "%8.2f %8.2f %8.2f %8.2f \n", 2.*lorxw1, lorxw1, 2.*loryw1, loryw1);
    fclose(f_lor_new);


    delete[] xproj;
    delete[] yproj; 
    delete[] xprojt; 
    delete[] yprojt; 
    delete[] pixev;
    delete[] qsum; 
    delete[] xh; 
    delete[] yh; 
    delete[] iproj; 
    delete[] jproj; 
    delete[] i2d; 
    delete[] j2d;
}


int main(){
    gen_xy_template2d();
    return 0;
}

