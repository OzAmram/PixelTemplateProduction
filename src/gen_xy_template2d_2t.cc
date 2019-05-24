
/* *******************************************************************   
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
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <TGraph.h>
#include <TCanvas.h>
#include <TF1.h>


Double_t fit_fn(Double_t *x, Double_t *par){
    //sqrt of 4th order polynomial, avoid negative
    Float_t xx = x[0];
    Double_t arg = par[0] + par[1]*xx + par[2]*xx*xx + par[3] * xx*xx*xx + par[4]*xx*xx*xx*xx;
    if(arg > 0.) return TMath::Sqrt(arg);
    else return 0.;
}


void gen_xy_template(const int nevents = 30000, const int npt = 200, const int nprm = 5, 
		     const int nqpt = 150000,
		     const int maxarg = 5){


    //size of the x and y arrays
    const int Nx = 21;
    const int Ny = 13;
    const int NOPT = NOPT;

    double  pixelt[Nx][Ny], bixin[Nx][Ny], pixel[Nx][Ny],
    xsum[Nx], ysum[Ny],
    ier[nprm], xpar[nprm][4], spar[nprm][4],
    xtemp[Nx][9], xtemp2[Nx][9],
    ytemp[Ny][9], ytemp2[Ny][9],
    xytemp[Nx][Ny][4][4];

    int nxntry[9], nyntry[9], nxytry[4][4], jy[Ny];

    //declare larger arrays dynamically to avoid stack overflow

    double *pixev = new double[Nx*Ny*nevents];
    double *qsum = new double[nevents];
    double *xh = new double[nevents];
    double *yh = new double[nevents];
    int *i2d = new int[nevents];
    int *j2d = new int[nevents];

    bool storep = false;

    /*
    */


    //size of pixel (in pixelav coordinates)
    float xsize, ysize, zsize;


    //define middle of dimension ranges
    const int NHx = Nx/2; // why no +1?
    const int NHy = Ny/2;

    FILE *f_config = fopen("pix_2t.proc", "r");

    if(f_config ==0){
        printf("Can't find pix_2t.proc! Exiting \n");
        exit(1);
    }

    int file_start, num_files, linear;
    float q100, q101, fq100, noise, fcal, fgain, rnoise;
    fscanf(f_config, "%i %i %f %f %f %f %f %f %f %i", 
            &file_start, &num_files, &noise, &q100, &q101, 
            &fq100, &fcal, &fgain, &rnoise, &linear);
    printf("processing %i files starting at file %i \n"
            "thr0 = %f, thr1 = %f rms thres = %f \n"
            "preamp noise = %f cluster noise frac = %f gain noise frac %f \n"
            "readout noise %f linear response %i \n ", 
            num_files, file_start, q100, q101, 
            fq100, noise, fcal, fgain, rnoise, linear);

    fclose(f_config);

    //define RMS of noise
    double rten = 10.;
    thr= q100;
    thr10=0.1*thr;


    bool tstsig = true;
    if (linear == 2) tstsig = false;


    char fname[100];
    char header[120];
    char buffer[120];

    //loop over files
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

        float half_xsize = xsize/2.0;
        float half_ysize = ysize/2.0;

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



            float bixin[Ny][Nx];
            //note the confusing index switch. input arrays are 13x21, we swap
            //to 21x13
            for(int j=0; j<Ny; j++){
                for (int i=0; i<Nx; i++){
                    fscanf(f_evts, " %f ", &bixin[i][j]);
                }
            }


            //  Pixelav gives hit position at face of pixel, translate to
            //  center of the pixel
            //  propagate from edge to center of pixel
            float xhit = x1 + (z_center - z1) * cosx/cosz;
            float yhit = y1 + (z_center - z1) * cosy/cosz;

            xh[n] = xhit;
            yh[n] = yhit;

            //hit position bins (9 bins within a pixel)
            //last bin is the same as first, translated 1 pixel over
            //divide into 2d bins for max charge
            int i2dn  = int(xhit/xsize*6 + 3.5);
            //restrict from 0 to 6;
            i2dn = std::min(6, std::max(i2dn, 0));
            i2d[n] = i2dn;

            int j2dn  = int(yhit/ysize*6 + 3.5);
            //restrict from 0 to 6;
            j2dn = std::min(6, std::max(j2dn, 0));
            j2d[n] = j2dn;

            qsum[n] = 0.;
            memset(xsum, 0., sizeof(xsum));
            memset(ysum, 0., sizeof(ysum));
            memset(pixel, 0., sizeof(pixel));
            memset(pixelt, 0., sizeof(pixelt));

            for(int i=0; i<Nx; i++){
                for(int j=0; j<Ny; j++){
                    double qin = bixin[i][j] * rten;
                    if(qin < 0.) qin = 0.;

                    pixev[i*(Ny *nevents) + j*nevents + n] = qin;

		    qsum[n] += qin;
                }
            }

            /*
               printf("print cluster \n");
               for(int i=0; i<Nx; i++){
               for(int j=0; j<Ny; j++){
               printf("%.1f ", pixelt[i][j]);
               }
               printf("\n");
               }
               */

            qavg += qsum[n];

        }

        int nmc = n;

        fclose(f_evts);
        qavg /= double(nmc);
        printf("finish reading run %i read %i events, qavg = %.1f \n", iFile, nmc, qavg);

        //zero things for this run
        memset(nxytry, 0, sizeof(nxytry));

        memset(xytemp, 0., sizeof(xytemp));
        memset(xytmp2, 0., sizeof(xytemp));

        int mx1=0;
        float sx1 = 0.;
        float sx12 = 0.;

        int my1=0;
        float sy1 = 0.;
        float sy12 = 0.;


        for(n=0; n<nmc; n++){

            //count number of pixels with less than average charge
	    // make 2-d template shapes
            if(qsum[n] < qavg){
                int k = i2d[n];
                int l = j2d[n];
                nxytry[k][l] +=1;

                for(int i=0; i < Nx; i++){
                    for(int j=0; j < Ny; j++){
		        xytemp[i][j][k][l] += pixev[i*(Ny *nevents) + j*nevents + n]; // is this pixev(i,j,n)?
			xytmp2[i][j][k][l] += pixev[i*(Ny *nevents) + j*nevents + n]**2;
                    }
                }

		if(k == 1){
		  nxytry[6][0] = nxytry[6][0]+1;
		  for(int i=1; i < Nx; i++){
                    for(int j=0; j < Ny; j++){
		      xytemp[i][j][6][0] = xytemp[i][j][6][0] + pixev[(i-1)*(Ny *nevents) + j*nevents + n];
                      xytmp2[i][j][6][0] = xytmp2[i][j][6][0] + pixev[(i-1)*(Ny *nevents) + j*nevents + n]**2;
		    }
		  }
		}

                if(k == 7){
		  nxytry[0][0] = nxytry[0][0]+1;
                  for(int i=1; i < Nx; i++){
                    for(int j=0; j < Ny; j++){
                      xytemp[i-1][j][0][0] = xytemp[i-1][j][0][0] + pixev[i*(Ny *nevents) + j*nevents + n];
                      xytmp2[i-1][j][0][0] = xytmp2[i-1][j][0][0] + pixev[i*(Ny *nevents) + j*nevents + n]**2;
                    }
                  }
		}

		if(l == 1){
                  nxytry[0][6] = nxytry[0][6]+1;
                  for(int i=1; i < Nx; i++){
                    for(int j=0; j < Ny; j++){
                      xytemp[i][j][6][0] = xytemp[i][j][6][0] + pixev[(i-1)*(Ny *nevents) + j*nevents + n];
                      xytmp2[i][j][6][0] = xytmp2[i][j][6][0] + pixev[(i-1)*(Ny *nevents) + j*nevents + n]**2;
                    }
                  }
                }

                //hit position bin
                k = iproj[n];

                //renorm distributions and print them
                if(k < 0 || k > 8){
                    printf("Problem with event %i, index (k) is =%i \n", n, k);
                }
                else if(k==0 || k==8){
                    nxntry[0] += 1;
                    nxntry[8] +=1;
                }
                else{
                    nxntry[k] +=1;
                }
                for(int i=0; i<Nx; i++){
                    xtemp[i][k] += xproj[i*nevents+n];
                    xtemp2[i][k] += xproj[i*nevents+n]*xproj[i*nevents+n];
                    //0th bin of new pixel same as last bin of previous pixel
                    if(k==0 && i>0){
                        xtemp[i][8] += xproj[(i-1)*nevents+n];
                        xtemp2[i][8] += xproj[(i-1)*nevents+n] * xproj[(i-1)*nevents+n];
                    }
                    else if(k==8 && i>0){
                        xtemp[i-1][0] += xproj[i*nevents+n];
                        xtemp2[i-1][0] += xproj[i*nevents+n] * xproj[i*nevents+n];
                    }
                }

                l = jproj[n];


                if(l < 0 || l > 8){
                    printf("Problem with event %i, index (l) is =%i \n", n, l);
                }
                else if(l==0 || l==8){
                    nyntry[0] += 1;
                    nyntry[8] +=1;
                }
                else{
                    nyntry[l] +=1;
                }
                for(int j=0; j<Ny; j++){
                    ytemp[j][l] += yproj[j*nevents+n];
                    ytemp2[j][l] += yproj[j*nevents+n]*yproj[j*nevents+n];
                    if(l==0 && j>0){

                        ytemp[j][8] += yproj[(j-1)*nevents+n];
                        ytemp2[j][8] += yproj[(j-1)*nevents+n] * yproj[(j-1)*nevents+n];
                    }
                    else if(l==8 && j>0){
                        ytemp[j-1][0] += yproj[j*nevents+n];
                        ytemp2[j-1][0] += yproj[j*nevents+n] * yproj[j*nevents+n];
                    }
                }
            } //end if over less than avg pixel charge
        }//loop over events
        for(int i=0; i<9; i++){
            printf("%i ", nxntry[i]);
        }
        printf("\n");


        //find maximum average pixel signal
        float pixmax = 0.;

        for(int i =0; i<Nx; i++){
            for(int j =0; j<Ny; j++){
                for(int k =0; k<4; k++){
                    for(int l =0; l<4; l++){
                        float sigxy = xytemp[i][j][k][l] / float(nxytry[k][l]);
                        if(sigxy > pixmax) pixmax = sigxy;
                    }
                }
            }
        }

        if(mx1 > 10){
            sx1 /= float(mx1);
            sx12 /= float(mx1);
            sx12 = sqrt((sx12 - sx1*sx1)/mx1);
        }
        printf("Number of 1 x-clusters %i, dx %7.2f +/- %7.2f \n", mx1, sx1, sx12);

        if(mx1 > 10 && fabs(cotb) < cotbmn){
            lorxw1 = sx1;
            dx1sig = sx12;
            nx1max = mx1;
            cotbmn = fabs(cotb);
        }

        if(my1 > 10){
            sy1 /= float(my1);
            sy12 /= float(my1);
            sy12 = sqrt((sy12 - sy1*sy1)/float(my1));
        }
        printf("Number of 1 y-clusters %i, dy %7.2f +/- %7.2f \n", my1, sy1, sy12);

        //why extra condition for y? (ymax) 
        if(my1 > 10 && fabs(cota) < cotamn && my1 > ny1max){
            loryw1 = sy1;
            dy1sig = sy12;
            ny1max = my1;
            cotamn = fabs(cota);
        }


        //calculate average and variance of charge
        for(int k=0; k < 9; k++){
            for(int i=0; i<Nx; i++){
                xtemp[i][k] /= float(nxntry[k]);
                xtemp2[i][k] = xtemp2[i][k]/float(nxntry[k]) - xtemp[i][k]*xtemp[i][k];
                xtemp[i][k] = std::max(0., xtemp[i][k]);
            }
        }
        for(int l=0; l < 9; l++){
            for(int j=0; j<Ny; j++){
                ytemp[j][l] /= float(nyntry[l]);
                ytemp2[j][l] = ytemp2[j][l]/float(nyntry[l]) - ytemp[j][l]*ytemp[j][l];
                ytemp[j][l] = std::max(0., ytemp[j][l]);
            }
        }


        //analyze error vs signals for entry and exit sides of cluster
        char plot_title[100];

        //x projection fits
        double sxmax = 0.;
        double ssxmax = 0.;
        std::vector<float> xsignal1, xssignal1, xsignal2, xssignal2;

        //start with zeros to get better fit
        xsignal1.push_back(0.);
        xssignal1.push_back(0.);
        xsignal2.push_back(0.);
        xssignal2.push_back(0.);

        for (int k=0; k<9; k++){
            for(int i=0; i<Nx; i++){
                //separate two sides of cluster
                if(xtemp[i][k] > 20.){
                    if(i<= NHx){
                        xsignal1.push_back(xtemp[i][k]);
                        xssignal1.push_back(sqrt(xtemp2[i][k]));

                        //why sxmax from first side only?
                        sxmax = std::max(sxmax, xtemp[i][k]);
                        ssxmax = std::max(ssxmax, sqrt(xtemp2[i][k]));
                    }
                    if(i>=NHx){
                        xsignal2.push_back(xtemp[i][k]);
                        xssignal2.push_back(sqrt(xtemp2[i][k]));
                    }

                }
            }
        }
        float guess = ssxmax*ssxmax/sxmax;



        TGraph *g_xtemp1 = new TGraph(xsignal1.size(), xsignal1.data(), xssignal1.data());
        g_xtemp1->SetTitle("X projections: Charge Variance vs. Charge");
        TGraph *g_xtemp2 = new TGraph(xsignal2.size(), xsignal2.data(), xssignal2.data());
        TCanvas *c_xtemp = new TCanvas("c_xtemp", "", 800, 800);

        TF1 *f_xtemp1 = new TF1("f_xtemp1", fit_fn,
                0.00, sxmax*1.01, nprm);
        f_xtemp1->FixParameter(0, 0.);
        f_xtemp1->SetParameter(1, guess);
        f_xtemp1->SetParameter(2, 0.);
        f_xtemp1->SetParameter(3, 0.);
        f_xtemp1->FixParameter(4, 0.);

        f_xtemp1->SetParError(0, 0.);
        f_xtemp1->SetParError(1, 1.);
        f_xtemp1->SetParError(2, 1.);
        f_xtemp1->SetParError(3, 1.);
        f_xtemp1->SetParError(4, 0.);
        f_xtemp1->SetLineColor(kRed);

        TF1 *f_xtemp2 = new TF1("f_xtemp2", fit_fn,
                0.00, sxmax*1.01, nprm);
        f_xtemp2->FixParameter(0, 0.);
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

        g_xtemp1->Fit(f_xtemp1, "VWR");
        g_xtemp1->Draw("AP");
        g_xtemp1->SetMarkerSize(1.4);
        g_xtemp1->SetMarkerStyle(20);


        g_xtemp2->Fit(f_xtemp2, "VWR");
        g_xtemp2->Draw("P same");
        g_xtemp2->SetMarkerSize(1.4);
        g_xtemp2->SetMarkerStyle(21);

        sprintf(plot_title, "sigmax_%i.png", iFile);
        c_xtemp->Print(plot_title);


        //y projection fits
        double symax = 0.;
        double ssymax = 0.;
        std::vector<float> ysignal1, yssignal1, ysignal2, yssignal2;

        //start with zeros to get better fit
        ysignal1.push_back(0.);
        yssignal1.push_back(0.);
        ysignal2.push_back(0.);
        yssignal2.push_back(0.);


        for (int k=0; k<9; k++){
            for(int i=0; i<Ny; i++){
                //separate two sides of cluster
                if(ytemp[i][k] > 20.){
                    if(i<= NHy){
                        ysignal1.push_back(ytemp[i][k]);
                        yssignal1.push_back(sqrt(ytemp2[i][k]));

                        symax = std::max(symax, ytemp[i][k]);
                        ssymax = std::max(ssymax, sqrt(ytemp2[i][k]));
                    }
                    if(i>= NHy){
                        ysignal2.push_back(ytemp[i][k]);
                        yssignal2.push_back(sqrt(ytemp2[i][k]));
                    }

                }
            }
        }
        guess = ssymax*ssymax/symax;


        TGraph *g_ytemp1 = new TGraph(ysignal1.size(), ysignal1.data(), yssignal1.data());
        g_ytemp1->SetTitle("Y projections: Charge Variance vs. Charge");
        TGraph *g_ytemp2 = new TGraph(ysignal2.size(), ysignal2.data(), yssignal2.data());
        TCanvas *c_ytemp = new TCanvas("c_ytemp", "", 800, 800);

        TF1 *f_ytemp1 = new TF1("f_ytemp1", fit_fn,
                0.0, symax*1.01, nprm);
        f_ytemp1->FixParameter(0, 0.);
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
                0.0, symax*1.01, nprm);
        f_ytemp2->FixParameter(0, 0.);
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

        g_ytemp1->Fit(f_ytemp1, "VWR");
        g_ytemp1->Draw("AP");
        g_ytemp1->SetMarkerSize(1.4);
        g_ytemp1->SetMarkerStyle(20);


        g_ytemp2->Fit(f_ytemp2, "VWR");
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
        for(int i=0; i<ysignal2.size(); i++){
            printf("%.1f %.1f \n", ysignal2[i], yssignal2[i]);
        }
        */
            


        //open summary files
        char file_out[100];

        //x template file
        sprintf(file_out, "ztemp_%i.txt", iFile);
        FILE *ztemp_file = fopen(file_out, "w+");
        fprintf(ztemp_file, "%9.6f %9.6f %9.6f \n", cosx, cosy, cosz);
        if(clslnx > (0.8*xsize) || clslny > (0.4*ysize)){
            fprintf(ztemp_file, "%8.1f %8.1f %8.1f \n", qavg, sxmax, pixmax);

            for(int k=0; k <2; k++){
                for(int p=0; p<nprm; p++){
                    fprintf(ztemp_file, "%15.8E ", xpar[p][k]);
                    spar[p][k] = xpar[p][k];
                }
                fprintf(ztemp_file, "\n");
            }
            spxmax = sxmax;
            storep = true;
        }
        // If the cluster length is smaller than a pixel, there aren't enough points for a reliable
        // fit.  Use the last (1 pixel length cluster) params instead for the z-template
        else if(storep){
            fprintf(ztemp_file, "%8.1f %8.1f %8.1f \n", qavg, spxmax, pixmax);
            for(int k=0; k <2; k++){
                for(int p=0; p<nprm; p++){
                    fprintf(ztemp_file,"%15.8E ", xpar[p][k]);
                    spar[p][k] = xpar[p][k];
                }
                fprintf(ztemp_file, "\n");
            }

        }
        else{
            fprintf(ztemp_file, "%8.1f %8.1f %8.1f \n", qavg, sxmax, pixmax);
            for(int k=0; k <2; k++){
                for(int p=0; p<nprm; p++){
                    fprintf(ztemp_file, "%15.8E ", xpar[p][k]);
                    spar[p][k] = xpar[p][k];
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
                spar[p][k] = xpar[p][k];
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


    delete[] xproj, yproj, qsum, pixev,
        xh, yh, xprojt, yproj;
    delete[] i2d, j2d;
}


int main(){
    gen_xy_template();
    return 0;
}

