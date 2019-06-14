/*
 * gen_zp_template2d.cc
 * Original by Morris Swartz
 * Updated by Oz Amram          May 2019
 *
 * Create template from previous generated files.
 * Make error objects by calling 2d reco with created templates.

 *   Change notes:
 *   Template hit reconstruction algorithms, add FPix to templates, add double pixels to templates, 
 *   change dp calling seq (Add PSI46v2 response function and allow for zero-suppressed ROC output)
 *   Change to Root historgraams for long-term compatibility
 *   Add angle vs resolution for templates and "Standard Algorithms"
 *   Tune CMSSW simulation for template 1 reconstruction
 *   Change standard algorithm to always use edge method for y-reconstruction
 *   Add Estar template number 4
 *   Do cosmics
 *   Add clustering  algorithm
 *   Change response function, parameters to that used for real analysis
 *   Same as test_code_v9 but with weighting for flat eta distributions in the resolution histograms
 *   Add 2-threshold simulation
 *   Work on templates with non-zero scale factors, ask for template IDs to use
 *   Low seed threshold verion for thin pixels
 *   Remove generic reco to avoid needing generror objects
 *   Use CMSSW data structures from VI
 *   Write integers for the templates
 */

#include "template_utils.h"



// Main program  

int main(int argc, char *argv[])
{

    // Local variables 
    std::vector<float> pvec(6), wgauss(TYSIZE), vgauss(TYSIZE), xgauss(TYSIZE), ygauss(TYSIZE), zgauss(TYSIZE);
    static bool fpix;
    float pixin[TXSIZE][TYSIZE];
    bool ydouble[TYSIZE], xdouble[TXSIZE];
    float ytemp[9][TYSIZE], xtemp[9][TXSIZE], xpar[2][5], ypar[2][5];
    float sxmax, symax, sxmaxx, symaxx, cosx, cosy, cosz;
    static float thick, xsize, ysize, noise, zcen, gain_frac, q100_frac, common_frac, readout_noise, qscale, qperbit;
    static float qavg,  clslnx, clslny, fbin[3] = {1.5f, 1.0f, 0.85f};
    static float xhit, yhit, xrec, yrec, sigmax, sigmay, probx, proby, probQ,  signal, cotalpha, cotbeta, qclust, locBz, locBx,  pixmax;
    static float pixmaxy, pixmaxx;
    static int startfile, neh, nevent, tempID, nbad, ngood, non_linear, icol, ndcol, numrun; 
    int  id,NTy, NTyx,NTxx,IDtype;

    int imsort[60];
    float qmsort[60];


    static float Bfield,Vbias,temp,fluenc;
    static vector<int> nbin(5,0);
    float deltay;
    int  ierr, qbin, qb, jmin, jmax, imin, imax, numadd, idcol;
    int mrow = TXSIZE, mcol = TYSIZE;
    const int TXSHIFT = (TXSIZE - T2XSIZE)/2;
    double dx, dy, tote, bade, adc;
    float scaley, scalex, scalx[4], scaly[4], delyavg, delysig, offsetx[4], offsety[4];
    static int numbits;
    static float q100, q101, q50, q51, q10, qmax; 




    const double gain = 3.19;
    const double ped = 16.46;
    const double p0 = 0.01218;
    const double p1 = 0.711;
    const double p2 = 203.;
    const double p3 = 148.;	
    static double vcal = 47.;	
    static double vcaloffst = 60.;
    const float fmax = 0.5f;
    int write_temp_header, use_l1_offset;

    const int nvers = 21;

    float qin;
    static char infile[120], label[160], header[120], outfile0[120], outfile1[120], outfile2[120];
    int triplg(std::vector<float>&);
    //	int random(void);

    float cluster[TXSIZE][TYSIZE], clust[TXSIZE][TYSIZE], rclust[TXSIZE][TYSIZE];
    bool bclust[TXSIZE][TYSIZE];
    std::pair<int, int> pixel, max;

    FILE *output_file;

    struct timeval now0, now1;
    struct timezone timz;
    long deltas, deltaus;
    double deltat;

    TF1 *cfunc0 = new TF1("chisquare0",chisquare0,0.,100.,3);
    cfunc0->SetParLimits(1,0.01,50);
    cfunc0->SetParNames("norm","mean","scale");

    TF1 *cfunc1 = new TF1("chisquare1",chisquare1,0.,50.,3);
    cfunc1->SetParNames("norm","mean","scale");

    TF1 *vfunc = new TF1("vavilov",vavilov,-10.,50.,4);
    vfunc->SetParNames("norm","mean","sigma","kappa");

    TCanvas* c1 = new TCanvas("c1", header, 800, 800);
    c1->SetFillStyle(4000);


    //  Read which data and inputs to use (use c file i/o which is much faster than c++ i/o) 

    FILE *config_file = fopen("pix_2t.proc", "r");
    if (config_file==NULL) {
        printf("no pixel initialization file found \n");
        return 0;
    }


    char extra[80];
    fscanf(config_file,"%d %d %f %f %f %f %f %f %f %d %f %s", &startfile, &numrun, &noise, &q100, 
            &q101, &q100_frac, &common_frac, &gain_frac, &readout_noise, &non_linear, &qscale, &extra[0]);
    printf("processing %d files starting from %d, noise = %f, threshold0 = %f, threshold1 = %f," 
            "rms threshold frac = %f, common_frac = %f, gain fraction = %f, readout noise = %f, nonlinear_resp = %d \n", 
            numrun, startfile, noise, q100, q101, q100_frac, common_frac, gain_frac, readout_noise, non_linear);


    fscanf(config_file, " %d %d ", &use_l1_offset, &write_temp_header);
    fscanf(config_file, " %d %d %d %d %d %f %f %f %f %f",  &id, &NTy, &NTyx,&NTxx, &IDtype, &Bfield, &Vbias, &temp, &fluenc, &qscale);
    printf("Using params: Use_l1_offset=%d, write_temp_header=%d, ID=%d NTy=%d NTyx=%d NTxx=%d Dtype=%d Bfield=%.2f Bias Voltage = %.1f temparature = %.0f fluence = %.2f q-scale = %.4f \n",
            use_l1_offset, write_temp_header, id, NTy, NTyx, NTxx, IDtype, Bfield, Vbias, temp, fluenc, qscale);

    fclose(config_file);

    //  Calculate 50% of threshold in q units and enc noise in adc units

    q50=0.5*q100;
    q51=0.5*q101;
    q10=0.2*q50;

    //  Open template output file

    sprintf(infile,"template_summary_zp%5.5d.out",startfile);
    /*
       output_file = fopen(infile, "w");
       if (output_file==NULL) {
       printf("couldn't open template output file/n");
       return 0;
       }
       */

    //  Open Lorentz summary file and read stored quantities

    sprintf(infile,"lorentz_widths_new%5.5d.out",startfile);
    FILE *lorw_file = fopen(infile, "r");
    if (lorw_file==NULL) {
        printf("couldn't find Lorentz summary file/n");
        return 0;
    }
    float lorwdy, lorbsy, lorwdx, lorbsx;
    fscanf(lorw_file,"%f %f %f %f", &lorwdy, &lorbsy, &lorwdx, &lorbsx);
    fclose(lorw_file);
    //  Flip the signs to convert from pixelav to cmssw conventions
    lorwdy = -lorwdy;
    lorbsy = -lorbsy;
    lorwdx = -lorwdx;
    lorbsx = -lorbsx;

    if(use_l1_offset) {
        printf("using L1 parameters \n");
        vcal = 50.;	
        vcaloffst = 670.;
    }

    // Define the histograms to be used at each angle pair

    double  halfys=200.;
    double  halfxs=100.;
    int nx=100;	
    gStyle->SetOptFit(101);
    gStyle->SetHistLineWidth(2);
    static vector<TH1F*> hp(26);
    hp[0] = new TH1F("h101","dy_temp (all sig); #Deltay (#mum)",nx,-halfys,halfys);
    hp[1] = new TH1F("h102","dy_temp (signal @> 1.5mn); #Deltay (#mum)",nx,-halfys,halfys);      
    hp[2] = new TH1F("h103","dy_temp (1.5mn @> signal @> 1.0mn); #Deltay (#mum)",nx,-halfys,halfys);      
    hp[3] = new TH1F("h104","dy_temp (1.0mn @> signal @> 0.85mn); #Deltay (#mum)",nx,-halfys,halfys);     
    hp[4] = new TH1F("h105","dy_temp (0.85mn @> signal); #Deltay (#mum)",nx,-halfys,halfys);      
    hp[5] = new TH1F("h201","pully_temp (all sig); #Deltay/rmsy",100,-10.,10.);
    hp[6] = new TH1F("h202","pully_temp (signal @> 1.5mn); #Deltay/rmsy",100,-10.,10.);
    hp[7] = new TH1F("h203","pully_temp (1.5mn @> signal @> 1.0mn); #Deltay/rmsy",100,-10.,10.);
    hp[8] = new TH1F("h204","pully_temp (1.0mn @> signal @> 0.85mn); #Deltay/rmsy",100,-10.,10.);
    hp[9] = new TH1F("h205","pully_temp (0.85mn @> signal); #Deltay/rmsy",100,-10.,10.);
    hp[10] = new TH1F("h106","dx_temp (all sig); #Deltax (#mum)",nx,-halfxs,halfxs);
    hp[11] = new TH1F("h107","dx_temp (signal @> 1.5mn); #Deltax (#mum)",nx,-halfxs,halfxs);      
    hp[12] = new TH1F("h108","dx_temp (1.5mn @> signal @> 1.0mn); #Deltax (#mum)",nx,-halfxs,halfxs);      
    hp[13] = new TH1F("h109","dx_temp (1.0mn @> signal @> 0.85mn); #Deltax (#mum)",nx,-halfxs,halfxs);      
    hp[14] = new TH1F("h110","dx_temp (0.85mn @> signal); #Deltax (#mum)",nx,-halfxs,halfxs);  
    hp[15] = new TH1F("h206","pullx_temp (all sig); #Deltax/rmsx",100,-10.,10.);  
    hp[16] = new TH1F("h207","pullx_temp (signal @> 1.5mn); #Deltax/rmsx",100,-10.,10.);  
    hp[17] = new TH1F("h208","pullx_temp (1.5mn @> signal @> 1.0mn); #Deltax/rmsx",100,-10.,10.);  
    hp[18] = new TH1F("h209","pullx_temp (1.0mn @> signal @> 0.85mn); #Deltax/rmsx",100,-10.,10.);  
    hp[19] = new TH1F("h210","pullx_temp (0.85mn @> signal); #Deltax/rmsx",100,-10.,10.);  
    hp[20] = new TH1F("h401","deltay [template - cluster]; #Deltay (pixels)",nx,-6.,6.);
    hp[21] = new TH1F("h405","deltax [template - cluster]; #Deltax (pixels)",nx,-6.,6.);
    hp[22] = new TH1F("h501","chisq [all]",400,0.,200.);
    hp[23] = new TH1F("h502","chisq/pixel",200,0.,50.);
    hp[24] = new TH1F("h503","npixels",50,-0.5,49.5);	
    hp[25] = new TH1F("h115","Cluster Charge; Q_clus (e)",250,0.,500000.);



    // Set style for the the histograms	

    for(int i=0; i<26; ++i) {
        hp[i]->SetLineColor(2);
        hp[i]->SetFillColor(38);
    }

    std::vector<std::pair<int, int> > pixlst;

    // Create template object

    std::vector< SiPixelTemplateStore > thePixelTemp_;
    SiPixelTemplate templ(thePixelTemp_);

    //  Set the ID to zero to flag the special reco mode

    tempID = 0;

    //  Determine current time

    gettimeofday(&now0, &timz);

    // Loop over angle pair runs

    int lfile = startfile+numrun;

    for(int ifile = startfile; ifile < lfile; ++ifile) {

        //  Read in 1D z template information first

        sprintf(infile,"./ztemp_%5.5d.txt",ifile);

        //  Open input file and read header info 

        FILE *ztemp_file = fopen(infile, "r");
        if (ztemp_file==NULL) {
            printf("no z-template file %s \n", infile);
            return 0;
        }

        fscanf(ztemp_file,"%f  %f  %f", &cosx, &cosy, &cosz);
        //	   printf("cosx/cosy/cosz = %f/%f/%f \n", cosx, cosy, cosz);

        fscanf(ztemp_file,"%f  %f  %f", &qavg, &symax, &pixmaxy);
        printf("qavg/symax/pixmaxy = %f/%f/%f \n", qavg, symax, pixmaxy);

        symaxx = fmax*symax;

        //flip to match cmssw coords
        for(int i = 1; i > -1; --i) {
            fscanf(ztemp_file,"%f %f %f %f %f", &ypar[i][0], &ypar[i][1], &ypar[i][2], &ypar[i][3], &ypar[i][4]);
        }

        for (int k=0; k < 9; ++k) {

            // Skip labels   
            get_label(ztemp_file, label, 160);
            printf("%d %s\n", k, label);
            //read in template
            for(int i=0; i<TYSIZE; i++){
                fscanf(ztemp_file, " %f ", &ytemp[k][i]);
            }
        }
       
        fclose(ztemp_file);

        // Calculate the mean cluster size in pixels

        float pzfrst=0.f, pzlast=0.f;
        for(int i=0; i<TYSIZE; ++i) {
            for (int k=6; k > -1; --k) {
                if(ytemp[k][i] > symaxx) {
                    float dzsig = ytemp[k][i] - ytemp[k+1][i];
                    float frac;
                    if(dzsig > 0.f) {
                        frac = (ytemp[k][i] - symaxx)/dzsig;
                        if(frac > 1.f) frac = 1.f;
                        if(frac < 0.f) frac = 0.f;
                    } else {
                        frac = 0.f;
                    }
                    pzfrst = i-(k+frac-4.f)/8.f;
                    goto firstz;
                }
            }
        }
firstz: ;
        for(int i=TYSIZE-1; i>-1; --i) {
            for (int k=1; k < 9; ++k) {
                if(ytemp[k][i] > symaxx) {
                    float dzsig = ytemp[k][i] - ytemp[k-1][i];
                    float frac;
                    if(dzsig > 0.f) {
                        frac = (ytemp[k][i] - symaxx)/dzsig;
                        if(frac > 1.f) frac = 1.f;
                        if(frac < 0.f) frac = 0.f;
                    } else {
                        frac = 0.f;
                    }
                    pzlast = i-(k+frac-4.f)/8.f;
                    goto secondz;
                }
            }
        }
secondz: clslny = pzlast-pzfrst;
         if(clslny < 0.f) clslny = 0.f;

         //  Read in 1D p template information

         sprintf(infile,"./ptemp_%5.5d.txt",ifile);


         FILE *ptemp_file = fopen(infile, "r");
         if (ptemp_file==NULL) {
             printf("no p-template file %s \n", infile);
             return 0;
         }

         fscanf(ptemp_file,"%f  %f  %f", &cosx, &cosy, &cosz);
         //	   printf("cosx/cosy/cosz = %f/%f/%f \n", cosx, cosy, cosz);

         fscanf(ptemp_file,"%f  %f  %f", &qavg, &sxmax, &pixmaxx);
         printf("qavg/sxmax/pixmaxx = %f/%f/%f \n", qavg, sxmax, pixmaxx);


         sxmaxx = fmax*sxmax;

         for(int i = 1; i > -1; --i) {
             fscanf(ptemp_file,"%f %f %f %f %f", &xpar[i][0], &xpar[i][1], &xpar[i][2], &xpar[i][3], &xpar[i][4]);
         }

         for (int k=0; k < 9; ++k) {
             // Skip labels   
             get_label(ptemp_file, label, 160);
             printf("%s\n", label);
             //read in template
             for(int i=0; i<TXSIZE; i++){
                 fscanf(ptemp_file, " %f ", &xtemp[k][i]);
             }
         }
         fclose(ptemp_file);

         // Calculate the mean cluster size in pixels

         // Calculates something like a full-width half max for the cluster
         // length
         float ppfrst=0.f, pplast=0.f;
         for(int i=0; i<TXSIZE; ++i) {
             for (int k=7; k > -1; --k) {
                 if(xtemp[k][i] > sxmaxx) {
                     float dpsig = xtemp[k][i] - xtemp[k+1][i];
                     float frac;
                     if(dpsig > 0.f) {
                         frac = (xtemp[k][i] - sxmaxx)/dpsig;
                         if(frac > 1.f) frac = 1.f;
                         if(frac < 0.f) frac = 0.f;
                     } else {
                         frac = 0.f;
                     }
                     ppfrst = i-(k+frac-4.f)/8.f;
                     goto firstp;
                 }
             }
         }
firstp: ;
        for(int i=TXSIZE-1; i>-1; --i) {
            for (int k=1; k < 9; ++k) {
                if(xtemp[k][i] > sxmaxx) {
                    float dpsig = xtemp[k][i] - xtemp[k-1][i];
                    float frac;
                    if(dpsig > 0.f) {
                        frac = (xtemp[k][i] - sxmaxx)/dpsig;
                        if(frac > 1.f) frac = 1.f;
                        if(frac < 0.f) frac = 0.f;
                    } else {
                        frac = 0.f;
                    }
                    pplast = i-(k+frac-4.f)/8.f;
                    goto secondp;
                }
            }
        }
secondp: clslnx = pplast-ppfrst;
         if(clslnx < 0.f) clslnx = 0.f;







         //  Open input file and read header info 



         SiPixelTemplateEntry * slice = new SiPixelTemplateEntry;

         // Copy info into the slice and reformat from pixelav coordinates to CMSSW local coordinates

         slice->runnum = ifile;

         for(int i = 0; i < 2; ++i) {
             for(int j=0; j<5; ++j) {
                 slice->xpar[i][j] = xpar[i][j];
                 slice->ypar[i][j] = ypar[i][j];
             }
         }

         slice->clslenx = clslnx;
         slice->clsleny = clslny;

         //for now hard code
         slice->sxone = 30.;
         slice->dxone = 0.;
         slice->sxtwo = 30.;
         slice->dxtwo = 0.;

         slice->syone = 30.;
         slice->dyone = 0.;
         slice->sytwo = 30.;
         slice->dytwo = 0.;


         slice->costrk[0] = -cosy;
         slice->costrk[1] = -cosx;
         slice->costrk[2] = -cosz;
         slice->cotalpha = cosy/cosz;
         slice->cotbeta = cosx/cosz;

         pixmax = std::max(pixmaxx, pixmaxy);
         printf("qavg/sxmax/pixmax = %f/%f/%f \n", qavg, sxmax, pixmax);

         slice->qavg = qavg;
         slice->pixmax = pixmax;
         slice->sxmax = sxmax;
         slice->symax = symax;

         //fill templates into slice
         for(int k = 0; k < 9; k++){
             //read in templates reversed because of coordinate difference
             //between cmssw and pixelav
             for(int i=0; i<TXSIZE; i++){
                 slice->xtemp[k][TXSIZE - 1 - i] = xtemp[k][i];
             }
             for(int j=0; j<TYSIZE; j++){
                 slice->ytemp[k][TYSIZE - 1 - j] = ytemp[k][j];
             }
         }


         //  Create an input filename for this run 


         sprintf(infile,"template_events_d%d.out",ifile);

         printf("opening file %s to get pixel events \n", infile);

         //  Open input file and read header info 

         FILE *events_file = fopen(infile, "r");
         if (events_file==NULL) {
             printf("no pixel data file found: %s \n", infile);
             return 0;
         }

         // Read-in a header string first and print it    

         get_label(events_file, header, 80);
         printf("Header: %s\n", header);


         fscanf(events_file,"%f  %f  %f", &ysize, &xsize, &thick);
         zcen = thick/2.;
         printf("xsize/ysize/thick = %f/%f/%f \n", xsize, ysize, thick);
         fpix = false;
         if(thick > 285.) {fpix = true;}

         nevent=0;
         nbad = 0;
         ngood = 0;
         tote = 0.;
         bade = 0.;

         if(write_temp_header && ifile==startfile) {


             /*
                fprintf(output_file,"%s", header);
                fprintf(output_file,"%d %d %4.2f %d %d %d %d %5.1f %5.1f %4.2f %5.3f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %f %4.2f %4.2f %4.2f \n",
                id,nvers,Bfield,NTy,NTyx,NTxx,IDtype,Vbias, temp,fluenc,qscale,q50,lorwdy,
                lorwdx,ysize,xsize,thick,q51,lorbsy,lorbsx,1.5f,1.0f,0.85f);
                */
         }

         for(int i=0; i<26; ++i) { hp[i]->Reset();}

         qperbit = pixmax/(pow(2.,(double)(numbits)-1.));

         float rnelec = 0.f;

         ndcol = TYSIZE/2 +1;
         std::vector<int> ndhit(ndcol, 0);

         // Loop until end of input file 

         while(fscanf(events_file,"%f %f %f %f %f %f %d", &pvec[0], &pvec[1], &pvec[2], &pvec[3], &pvec[4], &pvec[5], &neh) != EOF) {

             // read the input cluster 
             read_cluster(events_file, pixin);
             ++nevent;

             //	   if(nevent > 100) break;

             // Add noise and analog response to cluster, reformat for flipped barrel coordinate system 

             triplg(vgauss);
             pixlst.clear();
             for(int i=0; i<ndcol; ++i) {ndhit[i] = 0;}
             icol = 0;
             if(vgauss[1] < 0.) {icol = 1;}
             for(int j=0; j<TXSIZE; ++j) {
                 triplg(wgauss);
                 triplg(xgauss);
                 triplg(ygauss);
                 triplg(zgauss);
                 for(int i=0; i<TYSIZE; ++i) {
                     bclust[j][i] = false;
                     qin = (10.*pixin[j][i] + xgauss[i]*noise);
                     rclust[TXSIZE-1-j][TYSIZE-1-i] = qin;
                     if(qin < q100*(1.+wgauss[i]*q100_frac)) {
                         clust[TXSIZE-1-j][TYSIZE-1-i] = 0.;
                     } else {
                         idcol = (TYSIZE-1-i+icol)/2;
                         ++ndhit[idcol];
                         if(non_linear == 0) {
                             qin *= (1.+gain_frac*ygauss[i]);
                             signal = qperbit*((int)((qin + zgauss[i]*readout_noise)/qscale/qperbit) + 0.5);
                         } else {
                             adc = (double)((int)(p3+p2*tanh(p0*(qin + vcaloffst)/(7.0*vcal) - p1)));
                             // The simulated clusters assume that qscale = 1 [no charge scaling here]
                             //						 signal = ((float)((1.+gain_frac*ygauss[i])*(vcal*gain*(adc-ped))) - vcaloffst + zgauss[i]*readout_noise)/qscale;		
                             signal = ((float)((1.+gain_frac*ygauss[i])*(vcal*gain*(adc-ped))) - vcaloffst + zgauss[i]*readout_noise);							 
                         }	 
                         clust[TXSIZE-1-j][TYSIZE-1-i] = (1.+vgauss[0]*common_frac)*signal;
                     }
                 }
             }

             // Simulate the second, higher threshold in single dcol hits

             for(int j=0; j<TXSIZE; ++j) {
                 for(int i=0; i<TYSIZE; ++i) {
                     if(clust[j][i] > 0.) {
                         idcol = (i+icol)/2;
                         if(ndhit[idcol] == 1) {
                             // Apply higher threshold on single hits in double columns
                             if(rclust[j][i] < q101*(1.+wgauss[i]*q100_frac)) {
                                 clust[j][i] = 0.;
                             }
                         }
                     }
                 }
             }


             // Simulate the seed finding

             qmax = 0.;
             for(int j=0; j<TXSIZE; ++j) {
                 for(int i=0; i<TYSIZE; ++i) {
                     if(clust[j][i] > qmax) {
                         qmax = clust[j][i];
                         max.first = j; max.second = i;

                     }
                 }
             }

             if(qmax < 1500.) continue;


             // Simulate clustering around maximum signal (seed)

             pixlst.push_back(max);
             bclust[max.first][max.second] = true;

             std::vector<std::pair<int, int> > pixlst_copy;



             numadd = 1;
             //iterively find all non zero pixels near our seed
             while(numadd > 0){
                 //use copy of vector to avoid modifying vector as we loop through it
                 pixlst_copy = pixlst;
                 numadd = 0;
                 for ( auto pixIter = pixlst_copy.begin(); pixIter != pixlst_copy.end(); ++pixIter ) {
                     //max's are +2 because we are doing <max in the loop
                     jmin = pixIter->first-1; 
                     jmax = pixIter->first+2;
                     imin = pixIter->second-1;
                     imax = pixIter->second+2;
                     if(jmin < 0) {jmin = 0;}
                     if(jmax > TXSIZE) {jmax = TXSIZE;}
                     if(imin < 0) {imin = 0;}
                     if(imax > TYSIZE) {imax = TYSIZE;}
                     for(int j=jmin; j<jmax; ++j) {
                         for(int i=imin; i<imax; ++i) {
                             if(clust[j][i] > 0.) {
                                 if(!bclust[j][i]) {
                                     bclust[j][i] = true;
                                     pixel.first = j; pixel.second = i;
                                     pixlst.push_back(pixel);
                                     ++numadd;
                                 }
                             }
                         }
                     }
                 }
             }
             //only add pixels we 'found' when clustering to final cluster
             memset(cluster, 0., sizeof(cluster));

             qclust=0.;
             imin = TYSIZE;
             imax = 0;
             for (auto pixIter = pixlst.begin() ; pixIter != pixlst.end(); ++pixIter ) {
                 int j = pixIter->first; 
                 int i = pixIter->second;
                 cluster[j][i] = clust[j][i];
                 qclust += clust[j][i];
                 if(cluster[j][i] > 0.f) {
                     if(i < imin) imin = i;
                     if(i > imax) imax = i;
                 }
             }
             hp[25]->Fill((double)qclust, 1.);

             rnelec += qclust;

             // Calculate the hit coordinates in the flipped coordinate system 

             yhit = -(pvec[0] + (zcen-pvec[2])*pvec[3]/pvec[5]);
             xhit = -(pvec[1] + (zcen-pvec[2])*pvec[4]/pvec[5]);

             // Do the template analysis on the cluster 

             cotalpha = pvec[4]/pvec[5];
             cotbeta = pvec[3]/pvec[5];

             // No double pixels 

             for(int i=0; i<TYSIZE; ++i) {
                 ydouble[i] = false;
             }

             for(int j=0; j<TXSIZE; ++j) {
                 xdouble[j] = false;
             }




             //        if(fabs(cotbeta) < 2.1) continue;
             // Do the template analysis on the cluster 
             SiPixelTemplateReco::ClusMatrix clusterPayload{&cluster[0][0], xdouble, ydouble, mrow,mcol};
             locBx = 1.;
             if(cotbeta < 0.) locBx = -1.;
             locBz = locBx;
             if(cotalpha < 0.) locBz = -locBx;

             //  Sideload this template slice

             int speed = 0;
             templ.sideload(slice, IDtype, locBx, locBz, lorwdy, lorwdx, q50, fbin, xsize, ysize, thick);
             ierr = PixelTempReco1D(tempID, cotalpha, cotbeta, locBz, locBx,  clusterPayload, templ, yrec, sigmay, proby, xrec, sigmax, probx,  qbin, speed, probQ);
             if(ierr != 0) {
                 ++nbad; 
                 printf("reconstruction failed with error %d \n", ierr);
             } else {
                 int k= int(xhit/xsize * 8. + 4.5);
                 int l= int(yhit/ysize * 8. + 4.5);
                 ngood++;
                 qb = qbin;
                 ++nbin[qb];
                 dy = yrec - (TYSIZE/2)*ysize - yhit;
                 hp[0]->Fill(dy);
                 if(sigmay > 0.f) hp[5]->Fill(dy/sigmay);
                 hp[1+qbin]->Fill(dy);
                 if(sigmay > 0.f) hp[6+qbin]->Fill(dy/sigmay);
                 dx = xrec - (TXSIZE/2)*xsize - xhit;
                 //printf("Hit bins %i %i \n", k,l);
                 //printf(" %.3f %.3f \n", dx,dy);
                 hp[10]->Fill(dx);
                 if(sigmax > 0.f) hp[15]->Fill(dx/sigmax);
                 hp[11+qbin]->Fill(dx);
                 if(sigmax > 0.f) hp[16+qbin]->Fill(dx/sigmax);

             }

         }

         fclose(events_file);
         printf(" low q failures = %d, malformed clusters = %d, successful fits = %d \n", nbin[4], nbad, ngood);	   

         /*
          * Histograms plotting
          */


         for(int i=0; i<20; ++i) {hp[i]->Fit("gaus");}

         scaley = hp[5]->GetRMS();
         scalex = hp[15]->GetRMS();   
         delyavg = (float)hp[20]->GetMean();
         delysig = (float)hp[20]->GetRMS();    
         for(int i=0; i<4; ++i) {
             scaly[i] = hp[6+i]->GetRMS();
             scalx[i] = hp[16+i]->GetRMS();
             offsety[i] = hp[1+i]->GetMean();
             offsetx[i] = hp[11+i]->GetMean();
         }    

         //create a function with 2 parameters in the range [0.,100.]

         Double_t cp01 = hp[22]->GetMean()/2.;
         Double_t cp02 = 0.75;
         Double_t cp00 = hp[22]->GetEntries()*2./cp01; 
         cfunc0->SetParameters(cp00,cp01,cp02);
         hp[22]->Fit("chisquare0");
         Double_t cpar0[3];
         cfunc0->GetParameters(cpar0);
         printf("cpar0 = %lf/%lf/%lf \n", cpar0[0], cpar0[1], cpar0[2]);

         double pixels = hp[24]->GetMean();
         Double_t cp11 = cpar0[1]/pixels;
         Double_t cp12 = cpar0[2];
         Double_t cp10 = hp[23]->GetEntries()*2./cp11; 
         cfunc1->SetParameters(cp10,cp11,cp12);
         cfunc1->SetParLimits(1,cp11,cp11);
         cfunc1->SetParLimits(2,cp12,cp12);
         hp[23]->Fit("chisquare1");

         //create a function with 4 parameters in the range [-10,50]

         Double_t fp1 = hp[25]->GetMean();
         Double_t fp2 = 0.1*p1;
         Double_t fp3 = 0.02*p1/20000.;
         Double_t fp0 = hp[25]->GetEntries()*20000./fp1; 
         vfunc->SetParameters(fp0,fp1,fp2,fp3);
         hp[25]->Fit("vavilov");
         Double_t par[4];
         vfunc->GetParameters(par);

         //  Create an output filename for this run 

         sprintf(outfile0,"template_histos%5.5d_newResp.ps[",ifile);
         sprintf(outfile1,"template_histos%5.5d_newResp.ps",ifile);
         sprintf(outfile2,"template_histos%5.5d_newResp.ps]",ifile);
         c1->Clear();
         c1->Print(outfile0);
         for(int i=0; i<26; ++i) {
             hp[i]->Draw();
             c1->Print(outfile1);
             c1->Clear();
         }
         c1->Print(outfile2);
         c1->Clear();
         // Write this template entry to the output file

         /*
            fprintf(output_file,"%d %8.6f %8.6f %8.6f \n", ifile, -cosy, -cosx, -cosz);
            rnelec /=float(nevent);
            fprintf(output_file,"%7.1f %5.1f %5.1f %d %d %d %d \n", rnelec, pixmax, symax, slice->iymin, slice->iymax, slice->jxmin, slice->jxmax);
            for(i = 0; i < 2; ++i) {
            fprintf(output_file,"%e %e %e %e %e \n", slice->xypar[i][0], slice->xypar[i][1], slice->xypar[i][2], 
            slice->xypar[i][3], slice->xypar[i][4]);
            }
            for(i = 0; i < 2; ++i) {
            fprintf(output_file,"%e %e %e %e %e \n", slice->lanpar[i][0], slice->lanpar[i][1], slice->lanpar[i][2], 
            slice->lanpar[i][3], slice->lanpar[i][4]);
            }
            for(l=0; l<7; ++l) {
            for(k=0; k<7; ++k) {
            for(j=0; j<T2XSIZE; ++j) {
            fprintf(output_file," %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d \n", 
            slice->xytemp[k][l][0][j],slice->xytemp[k][l][1][j],slice->xytemp[k][l][2][j],
            slice->xytemp[k][l][3][j],slice->xytemp[k][l][4][j],slice->xytemp[k][l][5][j],
            slice->xytemp[k][l][6][j],slice->xytemp[k][l][7][j],slice->xytemp[k][l][8][j],
            slice->xytemp[k][l][9][j],slice->xytemp[k][l][10][j],slice->xytemp[k][l][11][j],
            slice->xytemp[k][l][12][j],slice->xytemp[k][l][13][j],slice->xytemp[k][l][14][j],
            slice->xytemp[k][l][15][j],slice->xytemp[k][l][16][j],slice->xytemp[k][l][17][j],
            slice->xytemp[k][l][18][j],slice->xytemp[k][l][19][j],slice->xytemp[k][l][20][j]);
            }
            }
            }

            fprintf(output_file,"%7.5lf %7.5lf %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f\n", cp11,cp12,offsetx[0],offsetx[1],
            offsetx[2],offsetx[3],offsety[0],offsety[1],offsety[2],offsety[3]);
            fprintf(output_file,"%3.1f %3.1f %3.1lf %3.1lf %6.3lf %6.3f %6.3f %4.2f %4.2f %3.1f\n", clslny,clslnx,par[1],par[2],par[3],scalex,scaley,delyavg,delysig,1.);
            fprintf(output_file,"%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %3.1f %3.1f\n", scalx[0],scalx[1],scalx[2],scalx[3],scaly[0],scaly[1],scaly[2],scaly[3],1.,qavg);
            */
    }
    // Close output file   

    //fclose(output_file);  

    /*  Determine current time */

    gettimeofday(&now1, &timz);
    deltas = now1.tv_sec - now0.tv_sec;
    deltaus = now1.tv_usec - now0.tv_usec;
    deltat = ((double)deltaus)/1000000.;
    deltat += (double)deltas;
    printf("ellapsed time = %f seconds \n", deltat);

    return 0;
} // MAIN__ 





