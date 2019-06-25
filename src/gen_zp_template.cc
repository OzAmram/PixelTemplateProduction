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
    float ytemp[9][TYSIZE], xtemp[9][TXSIZE], xpar[2][5], ypar[2][5];
    float sxmax, symax, sxmaxx, symaxx, cosx, cosy, cosz;
    static float thick, xsize, ysize, noise, zcen, gain_frac, q100_frac, common_frac, readout_noise, qscale, qperbit;
    static float qavg_raw,  clslnx, clslny; 
    static float xrec, yrec, sigmax, sigmay, probx, proby, probQ,  signal, locBz, locBx,  pixmax;
    static float pixmaxy, pixmaxx;
    static int startfile,  nbad, ngood, non_linear, numrun; 
    int  id,NTy, NTyx,NTxx,IDtype;


    std::multiset<float> qmsort;
    static float fbin[] = {1.5, 1.0, 0.85}; //bins of charge in Q_cluster / Q_avg

    const int nevents = 30000;

    float xhit[nevents], yhit[nevents], cotalpha, cotbeta;

    float  qsmear[nevents], qmerge[nevents], npix[nevents], xfrac[nevents][3], yfrac[nevents][3],
    qtotal[nevents];

    int nelec[nevents];

    float ***cluster = setup_3d_array(nevents, TXSIZE, TYSIZE); //Cluster as found by seeding + clustering algo
    float **xsum1 = setup_2d_array(nevents, TXSIZE);
    float **xsum2 = setup_2d_array(nevents, TXSIZE);

    float **ysum1 = setup_2d_array(nevents, TYSIZE);
    float **ysum2 = setup_2d_array(nevents, TYSIZE);
    static float Bfield,Vbias,temp,fluenc;
    float nqbin[5];
    const int TXSHIFT = (TXSIZE - T2XSIZE)/2;
    double dx, dy, tote, bade, adc;
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
    //	int random(void);

    float clust[TXSIZE][TYSIZE], rclust[TXSIZE][TYSIZE], sigraw[TXSIZE+2][TYSIZE+2];
    bool bclust[TXSIZE][TYSIZE];
    std::pair<int, int> pixel, max;

    FILE *temp_output_file, *generr_output_file;

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
    temp_output_file = fopen(infile, "w");
    if (temp_output_file==NULL) {
        printf("couldn't open template output file/n");
        return 0;
    }

    sprintf(infile,"generror_summary_zp%5.5d.out",startfile);
    generr_output_file = fopen(infile, "w");
    if (generr_output_file==NULL) {
        printf("couldn't open generr output file/n");
        return 0;
    }

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

    double halfxs=300.;
    int nx=120;
    double halfys=300.;
    int ny=120;
    double chiymx=48.;
    double chixmx=150.;


    const int y_temp_idx =0;
    const int x_temp_idx =10;
    const int y_generic_idx =20;
    const int x_generic_idx =30;
    const int charge_idx = 40;

    const int n_hists = 60;


    gStyle->SetOptFit(101);
    gStyle->SetHistLineWidth(2);
    static vector<TH1F*> hp(n_hists);
    hp[y_temp_idx + 0] = new TH1F("h101","dy_temp (all sig); #Deltay (#mum)",ny,-halfys,halfys);
    hp[y_temp_idx + 1] = new TH1F("h102","dy_temp (signal @> 1.5mn); #Deltay (#mum)",ny,-halfys,halfys);      
    hp[y_temp_idx + 2] = new TH1F("h103","dy_temp (1.5mn @> signal @> 1.0mn); #Deltay (#mum)",ny,-halfys,halfys);      
    hp[y_temp_idx + 3] = new TH1F("h104","dy_temp (1.0mn @> signal @> 0.85mn); #Deltay (#mum)",ny,-halfys,halfys);     
    hp[y_temp_idx + 4] = new TH1F("h105","dy_temp (0.85mn @> signal); #Deltay (#mum)",ny,-halfys,halfys);      
    hp[y_temp_idx + 5] = new TH1F("h106","chi2y_temp (single pix)",ny,0.,chiymx);
    hp[y_temp_idx + 6] = new TH1F("h107","chi2y_temp (signal @> 1.5mn)",ny,0.,chiymx);
    hp[y_temp_idx + 7] = new TH1F("h108","chi2y_temp (1.5mn @> signal @> 1.0mn)",ny, 0., chiymx);
    hp[y_temp_idx + 8] = new TH1F("h109","chi2y_temp (1.0mn @> signal @> 0.85mn)",ny,0., chiymx);
    hp[y_temp_idx + 9]= new TH1F("h110","chi2y_temp (0.85mn @> signal)",ny,0.,chiymx);

    hp[x_temp_idx + 0] = new TH1F("h201","dx_temp (all sig); #Deltax (#mum)",nx,-halfxs,halfxs);
    hp[x_temp_idx + 1] = new TH1F("h202","dx_temp (signal @> 1.5mn); #Deltax (#mum)",nx,-halfxs,halfxs);      
    hp[x_temp_idx + 2] = new TH1F("h203","dx_temp (1.5mn @> signal @> 1.0mn); #Deltax (#mum)",nx,-halfxs,halfxs);      
    hp[x_temp_idx + 3] = new TH1F("h204","dx_temp (1.0mn @> signal @> 0.85mn); #Deltax (#mum)",nx,-halfxs,halfxs);     
    hp[x_temp_idx + 4] = new TH1F("h205","dx_temp (0.85mn @> signal); #Deltax (#mum)",nx,-halfxs,halfxs);      
    hp[x_temp_idx + 5] = new TH1F("h206","chi2x_temp (single pix)",nx,0.,chixmx);
    hp[x_temp_idx + 6] = new TH1F("h207","chi2x_temp (signal @> 1.5mn)",nx,0.,chixmx);
    hp[x_temp_idx + 7] = new TH1F("h208","chi2x_temp (1.5mn @> signal @> 1.0mn)",nx, 0., chixmx);
    hp[x_temp_idx + 8] = new TH1F("h209","chi2x_temp (1.0mn @> signal @> 0.85mn)",nx,0., chixmx);
    hp[x_temp_idx + 9]= new TH1F("h210","chi2x_temp (0.85mn @> signal)",nx,0.,chixmx);

    hp[y_generic_idx + 0] = new TH1F("h301","dy_generic (all sig); #Deltay (#mum)",ny,-halfys,halfys);
    hp[y_generic_idx + 1] = new TH1F("h302","dy_generic (signal @> 1.5mn); #Deltay (#mum)",ny,-halfys,halfys);      
    hp[y_generic_idx + 2] = new TH1F("h303","dy_generic (1.5mn @> signal @> 1.0mn); #Deltay (#mum)",ny,-halfys,halfys);      
    hp[y_generic_idx + 3] = new TH1F("h304","dy_generic (1.0mn @> signal @> 0.85mn); #Deltay (#mum)",ny,-halfys,halfys);     
    hp[y_generic_idx + 4] = new TH1F("h305","dy_generic (0.85mn @> signal); #Deltay (#mum)",ny,-halfys,halfys);      

    hp[x_generic_idx + 0] = new TH1F("h401","dx_generic (all sig); #Deltax (#mum)",nx,-halfxs,halfxs);
    hp[x_generic_idx + 1] = new TH1F("h402","dx_generic (signal @> 1.5mn); #Deltax (#mum)",nx,-halfxs,halfxs);      
    hp[x_generic_idx + 2] = new TH1F("h403","dx_generic (1.5mn @> signal @> 1.0mn); #Deltax (#mum)",nx,-halfxs,halfxs);      
    hp[x_generic_idx + 3] = new TH1F("h404","dx_generic (1.0mn @> signal @> 0.85mn); #Deltax (#mum)",nx,-halfxs,halfxs);     
    hp[x_generic_idx + 4] = new TH1F("h405","dx_generic (0.85mn @> signal); #Deltax (#mum)",nx,-halfxs,halfxs);      


    hp[charge_idx + 0] = new TH1F("h100","Number generated e",150,0.,5000000.);	
    hp[charge_idx + 1] = new TH1F ("h500","Cluster Charge",250,0.,500000.);
    hp[charge_idx + 2] = new TH1F ("h501","npix(signal @> 1.5mn)",40,0.5,40.5);
    hp[charge_idx + 3] = new TH1F ("h502","npix(1.5mn @> signal @> 1.0mn)",40,0.5,40.5);
    hp[charge_idx + 4] = new TH1F ("h503","npix(1.0mn @> signal @> 0.85mn)",40,0.5,40.5);
    hp[charge_idx + 5] = new TH1F ("h504","npix(0.85mn @> signal)",40,0.5,40.5);
    hp[charge_idx + 6] = new TH1F ("h505","2 Cluster Merged Charge",500,0.,1000000.);
    hp[charge_idx + 7] = new TH1F ("h506","chi2y_2cls (signal @> 1.5mn)",nx,0.,chiymx);
    hp[charge_idx + 8] = new TH1F ("h507","chi2y_2cls (1.5mn @> signal @> 1.0mn)",nx,0.,  chiymx);
    hp[charge_idx + 9] = new TH1F ("h508","chi2y_2cls (1.0mn @> signal @> 0.85mn)",nx,0., chiymx);
    hp[charge_idx + 10] = new TH1F ("h509","chi2y_2cls (0.85mn @> signal)",nx,0.,chiymx);
    hp[charge_idx + 11] = new TH1F ("h510","chi2x_2cls (signal @> 1.5mn)",nx,0.,chiymx);
    hp[charge_idx + 12] = new TH1F ("h511","chi2x_2cls (1.5mn @> signal @> 1.0mn)",nx,0., chiymx);
    hp[charge_idx + 13] = new TH1F ("h512","chi2x_2cls (1.0mn @> signal @> 0.85mn)",nx,0.,   chiymx);
    hp[charge_idx + 14] = new TH1F ("h513","chi2x_2cls (0.85mn @> signal)",nx,0.,chiymx);
    hp[charge_idx + 15] = new TH1F ("h606","measured Q/generated Q",300,0.,1.5);



    // Set style for the the histograms	

    for(int i=0; i<hp.size(); ++i) {
        if(hp[i] == NULL) continue;
        hp[i]->SetLineColor(2);
        hp[i]->SetFillColor(38);
    }

    std::vector<std::pair<int, int> > pixlst;

    // Create template object

    std::vector< SiPixelTemplateStore > thePixelTemp_;
    SiPixelTemplate templ(thePixelTemp_);

    //  Set the ID to zero to flag the special reco mode

    int tempID = 0;

    //  Determine current time

    gettimeofday(&now0, &timz);

    // Loop over angle pair runs

    int lfile = startfile+numrun;

    for(int ifile = startfile; ifile < lfile; ++ifile) {

        memset(nqbin, 0., sizeof(nqbin));
        memset(qtotal, 0., sizeof(qtotal));
        memset(qsmear, 0., sizeof(qsmear));
        memset(qmerge, 0., sizeof(qmerge));
        memset(npix, 0., sizeof(npix));
        memset(xfrac, 0., sizeof(xfrac));
        memset(yfrac, 0., sizeof(yfrac));
        memset(nelec, 0., sizeof(nelec));
        zero_3d_array(cluster, nevents, TXSIZE, TYSIZE);
        zero_2d_array(xsum1, nevents, TXSIZE);
        zero_2d_array(xsum2, nevents, TXSIZE);
        zero_2d_array(ysum1, nevents, TYSIZE);
        zero_2d_array(ysum2, nevents, TYSIZE);

        //  Read in 1D z template information first

        sprintf(infile,"./ztemp_%5.5d.txt",ifile);

        //  Open input file and read header info 

        FILE *ztemp_file = fopen(infile, "r");
        if (ztemp_file==NULL) {
            printf("no z-template file %s \n", infile);
            return 0;
        }

        fscanf(ztemp_file,"%f  %f  %f", &cosy, &cosx, &cosz);
        //	   printf("cosx/cosy/cosz = %f/%f/%f \n", cosx, cosy, cosz);

        fscanf(ztemp_file,"%f  %f  %f", &qavg_raw, &symax, &pixmaxy);
        printf("qavg_raw/symax/pixmaxy = %f/%f/%f \n", qavg_raw, symax, pixmaxy);

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

         fscanf(ptemp_file,"%f  %f  %f", &cosy, &cosx, &cosz);
         //	   printf("cosx/cosy/cosz = %f/%f/%f \n", cosx, cosy, cosz);

         fscanf(ptemp_file,"%f  %f  %f", &qavg_raw, &sxmax, &pixmaxx);
         printf("qavg_raw/sxmax/pixmaxx = %f/%f/%f \n", qavg_raw, sxmax, pixmaxx);


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

         float qavg = 0.f; //average charge after threshholding effects

         // loop over all events once to get single pixel avgs.
         int read_events = 0;
         qmsort.clear();


         for(int n=0; n<nevents; n++){

             float x1,y1,z1;
             //x and y flipped order in input files
             if(fscanf(events_file,"%f %f %f %f %f %f %d", &y1, &x1, &z1, &cosy, &cosx, &cosz, &nelec[n]) == EOF){
                 printf("File %s ended early!! \n\n", infile);
                 break;
             }
             read_events++;

             // read the input cluster 
             read_cluster(events_file, pixin);


             cotalpha = cosx/cosz;
             cotbeta = cosy/cosz;

             //  Pixelav gives hit position at face of pixel, translate to
             //  3d center of the pixel
             //  propagate from edge to center of pixel
             //  negative sign to convert to cmssw coords
             //  these hit coords are centered at 0
             xhit[n] = -(x1 + (zcen - z1) * cotalpha);
             yhit[n] = -(y1 + (zcen - z1) * cotbeta);

             int ndcol = TYSIZE/2 +1;
             std::vector<int> ndhit(ndcol, 0);
             int idcol;

             //sigraw is zero padded to allow overflow in double col
             //projections
             memset(sigraw, 0., sizeof(sigraw));
             memset(clust, 0., sizeof(clust));

             // Add noise and analog response to cluster, reformat for flipped barrel coordinate system 


             triplg(vgauss);
             qsmear[n] = (1.+vgauss[0]*common_frac);
             pixlst.clear();
             for(int i=0; i<ndcol; ++i) {ndhit[i] = 0;}
             int icol = 0;
             if(vgauss[1] < 0.) {icol = 1;}
             for(int j=0; j<TXSIZE; ++j) {
                 triplg(wgauss);
                 triplg(xgauss);
                 triplg(ygauss);
                 triplg(zgauss);

                 //do main cluster
                 for(int i=0; i<TYSIZE; ++i) {
                     sigraw[TXSIZE-1-j][TYSIZE-1-i] = pixin[j][i];
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
                         clust[TXSIZE-1-j][TYSIZE-1-i] = qsmear[n]*signal;
                     }
                 }


             }

             // Simulate the second, higher threshold in single double col hits
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
             for(int i=0; i<TXSIZE; ++i) {
                 for(int j=0; j<TYSIZE; ++j) {
                     if(clust[i][j] > qmax) {
                         qmax = clust[i][j];
                         max.first = i; max.second = j;

                     }
                 }
             }

             if(qmax < 1500.) continue;


             // Simulate clustering around maximum signal (seed)

             pixlst.clear();
             pixlst.push_back(max);
             memset(bclust, false, sizeof(bclust));
             bclust[max.first][max.second] = true;

             std::vector<std::pair<int, int> > pixlst_copy;



             int numadd = 1;
             //iterively find all non zero pixels near our seed
             while(numadd > 0){
                 //use copy of vector to avoid modifying vector as we loop through it
                 pixlst_copy = pixlst;
                 numadd = 0;
                 for ( auto pixIter = pixlst_copy.begin(); pixIter != pixlst_copy.end(); ++pixIter ) {
                     //max's are +2 because we are doing <max in the loop
                     int imin = pixIter->first-1; 
                     int imax = pixIter->first+2;
                     int jmin = pixIter->second-1;
                     int jmax = pixIter->second+2;
                     if(imin < 0) {imin = 0;}
                     if(imax > TXSIZE) {imax = TXSIZE;}
                     if(jmin < 0) {jmin = 0;}
                     if(jmax > TYSIZE) {jmax = TYSIZE;}
                     for(int i=imin; i<imax; ++i) {
                         for(int j=jmin; j<jmax; ++j) {
                             if(clust[i][j] > 0.) {
                                 if(!bclust[i][j]) {
                                     bclust[i][j] = true;
                                     pixel.first = i; pixel.second = j;
                                     pixlst.push_back(pixel);
                                     ++numadd;
                                 }
                             }
                         }
                     }
                 }
             }

             float qmeas=0.;
             npix[n] = 0.;
             for (auto pixIter = pixlst.begin() ; pixIter != pixlst.end(); ++pixIter ) {
                 int i = pixIter->first; 
                 int j = pixIter->second;
                 cluster[n][i][j] = clust[i][j];
                 qmeas += clust[i][j];
                 npix[n] += 1.;
             }
             qtotal[n] = qmeas;


             //keep 60 smallest charges
             if(qmsort.size() < 60){
                 qmsort.insert(qmeas);
             }
             else if(qmeas < *(qmsort.end())){
                 qmsort.erase(qmsort.end());
                 qmsort.insert(qmeas);
             }


             qavg += qmeas;




             //do x double  row projections
             for(int j=0; j<TXSIZE/2; j++){
                 xsum1[n][j] = xsum2[n][j] = 0.;
                 for(int i=0; i<TYSIZE; i++){
                     int j1 = 2*j;
                     int j2 = 2*j+1;

                     //sigraw padded with extra 0's to allow overflow
                     qin = 10. * (sigraw[j1][i] + sigraw[j1+1][i]);
                     qin += xgauss[i]*noise;
                     if(qin > q101*(1.+wgauss[i]*q100_frac)) {
                         if(non_linear == 0) {
                             qin *= (1.+gain_frac*ygauss[i]);
                             signal = qperbit*((int)((qin + zgauss[i]*readout_noise)/qscale/qperbit) + 0.5);
                         } else {
                             adc = (double)((int)(p3+p2*tanh(p0*(qin + vcaloffst)/(7.0*vcal) - p1)));
                             // The simulated clusters assume that qscale = 1 [no charge scaling here]
                             //						 signal = ((float)((1.+gain_frac*ygauss[i])*(vcal*gain*(adc-ped))) - vcaloffst + zgauss[i]*readout_noise)/qscale;		
                             signal = ((float)((1.+gain_frac*ygauss[i])*(vcal*gain*(adc-ped))) - vcaloffst + zgauss[i]*readout_noise);							 
                         }	 
                         xsum1[n][j] += qsmear[n]*signal;
                     }

                     qin = 10. * (sigraw[j2][i] + sigraw[j2+1][i]);
                     qin += xgauss[i]*noise;
                     if(qin > q101*(1.+wgauss[i]*q100_frac)) {
                         if(non_linear == 0) {
                             qin *= (1.+gain_frac*ygauss[i]);
                             signal = qperbit*((int)((qin + zgauss[i]*readout_noise)/qscale/qperbit) + 0.5);
                         } else {
                             adc = (double)((int)(p3+p2*tanh(p0*(qin + vcaloffst)/(7.0*vcal) - p1)));
                             // The simulated clusters assume that qscale = 1 [no charge scaling here]
                             //						 signal = ((float)((1.+gain_frac*ygauss[i])*(vcal*gain*(adc-ped))) - vcaloffst + zgauss[i]*readout_noise)/qscale;		
                             signal = ((float)((1.+gain_frac*ygauss[i])*(vcal*gain*(adc-ped))) - vcaloffst + zgauss[i]*readout_noise);							 
                         }	 
                         xsum2[n][j] += qsmear[n]*signal;
                     }
                 }
             }


             //do y double col projections 
             for(int i=0; i<TYSIZE/2; ++i) {
                 ysum1[n][i] = ysum2[n][i] = 0.;
                 for(int j=0; j<TXSIZE; j++){
                     int i1 = 2*i;
                     int i2 = 2*i+1;

                     //sigraw padded with extra 0's to allow overflow
                     qin = 10. * (sigraw[j][i1] + sigraw[j][i1+1]);
                     qin += xgauss[i1]*noise;
                     if(qin > q101*(1.+wgauss[i1]*q100_frac)) {
                         if(non_linear == 0) {
                             qin *= (1.+gain_frac*ygauss[i1]);
                             signal = qperbit*((int)((qin + zgauss[i1]*readout_noise)/qscale/qperbit) + 0.5);
                         } else {
                             adc = (double)((int)(p3+p2*tanh(p0*(qin + vcaloffst)/(7.0*vcal) - p1)));
                             // The simulated clusters assume that qscale = 1 [no charge scaling here]
                             //						 signal = ((float)((1.+gain_frac*ygauss[i])*(vcal*gain*(adc-ped))) - vcaloffst + zgauss[i]*readout_noise)/qscale;		
                             signal = ((float)((1.+gain_frac*ygauss[i1])*(vcal*gain*(adc-ped))) - vcaloffst + zgauss[i1]*readout_noise);							 
                         }	 
                         ysum1[n][i] += qsmear[n]*signal;
                     }

                     //sigraw padded with extra 0's to allow overflow
                     qin = 10. * (sigraw[j][i2] + sigraw[j][i2+1]);
                     qin += xgauss[i2]*noise;
                     if(qin > q101*(1.+wgauss[i2]*q100_frac)) {
                         if(non_linear == 0) {
                             qin *= (1.+gain_frac*ygauss[i2]);
                             signal = qperbit*((int)((qin + zgauss[i2]*readout_noise)/qscale/qperbit) + 0.5);
                         } else {
                             adc = (double)((int)(p3+p2*tanh(p0*(qin + vcaloffst)/(7.0*vcal) - p1)));
                             // The simulated clusters assume that qscale = 1 [no charge scaling here]
                             //						 signal = ((float)((1.+gain_frac*ygauss[i])*(vcal*gain*(adc-ped))) - vcaloffst + zgauss[i]*readout_noise)/qscale;		
                             signal = ((float)((1.+gain_frac*ygauss[i2])*(vcal*gain*(adc-ped))) - vcaloffst + zgauss[i2]*readout_noise);							 
                         }	 
                         ysum2[n][i] += qsmear[n]*signal;
                     }
                 }
             }


         }

         qavg /= read_events;
         fclose(events_file);

         int nxone=0;
         int nyone=0;
         int nxtwo=0;
         int nytwo=0;

         float dxone=0;
         float dyone=0;
         float dxtwo=0;
         float dytwo=0;

         float sxone=0.;
         float syone=0.;
         float sxtwo=0.;
         float sytwo=0.;


         //compute averages over all clusters 
         for(int n=0; n<read_events; n++){

             float qmeas = qtotal[n];
             float rcorr = qmeas/float(nelec[n]); //ratio of measured charge to generated charge


             hp[charge_idx]->Fill(float(nelec[n]));
             hp[charge_idx + 1] ->Fill(qmeas);
             hp[charge_idx + 15]->Fill(rcorr);

             if(n>0){
                 qmerge[n] = qtotal[n] + qsmear[n]*qtotal[n-1]/qsmear[n-1];
                 hp[charge_idx+6] ->Fill(qmerge[n]);
             }


             float xsum[TXSIZE], ysum[TYSIZE];

             memset(xsum, 0., sizeof(xsum));
             memset(ysum, 0., sizeof(ysum));
             for(int i=0; i<TXSIZE; i++){
                 for(int j=0; j<TYSIZE; j++){
                     xsum[i] += cluster[n][i][j];
                     ysum[j] += cluster[n][i][j];
                 }
             }

             int xwidth=0;
             int ywidth=0;
             int xstart=0;
             int ystart=0;
             int xend=0;
             int yend=0;

             int xw1(0), xw2(0), yw1(0), yw2(0);
             int xc1(0), xc2(0), yc1(0), yc2(0);

             for(int j=0; j<TXSIZE; j++){
                 if(xsum[j] >0.){
                     if(xstart==0) xstart = j;
                     xwidth++;
                     xend = j;
                 }
             }
             for(int i=0; i<TYSIZE; i++){
                 if(ysum[i] >0.){
                     if(ystart==0) ystart = i;
                     ywidth++;
                     yend = i;
                 }
             }


             //do double col version
             for(int j=0; j<TXSIZE/2; j++){
                 if(xsum1[n][j] >0.){
                     if(xc1 == 0) xc1 = j;
                     xw1++;
                 }
                 if(xsum2[n][j] >0.){
                     if(xc2 == 0) xc2 = j;
                     xw2++;
                 }
             }
             for(int i=0; i<TYSIZE/2; i++){
                 if(ysum1[n][i] >0.){
                     if(yc1 == 0) yc1 = i;
                     yw1++;
                 }
                 if(ysum2[n][i] >0.){
                     if(yc2 == 0) yc2 = i;
                     yw2++;
                 }
             }

             //compute avg shift and variance of single pixel clusters
             if(xwidth ==1){
                 nxone++;
                 float x0 = (xstart - TXSIZE/2) *xsize;
                 float deltax = x0-xhit[n];
                 //printf("1pix x: %.1f %.1f \n", x0, xhit[n]);
                 dxone += deltax;
                 sxone += deltax*deltax;
             }
             if(ywidth ==1){
                 nyone++;
                 float y0 = (ystart - TYSIZE/2) *ysize;
                 float deltay = y0-yhit[n];
                 //printf("1pix y: %.1f %.1f \n", y0, yhit[n]);
                 dyone += deltay;
                 syone += deltay*deltay;
             }

             //do the same for double sized single pixel clusters
             if(xw1 ==1){
                 nxtwo++;
                 //want middle of our double sized
                 float x0 = (xc1*2 + 1 - float(TXSIZE)/2) *xsize;
                 float deltax = x0-xhit[n];
                 //printf("2pix x: %.1f %.1f \n", x0, yhit[n]);
                 dxtwo += deltax;
                 sxtwo += deltax*deltax;
             }
             if(xw2 ==1){
                 nxtwo++;
                 float x0 = (xc2*2 +2 - float(TXSIZE)/2) *xsize;
                 float deltax = x0-xhit[n];
                 //printf("2pix x: %.1f %.1f \n", x0, yhit[n]);
                 dxtwo += deltax;
                 sxtwo += deltax*deltax;
             }
             if(yw1 ==1){
                 nytwo++;
                 //want middle of our double sized
                 float y0 = (yc1*2 +1 - float(TYSIZE)/2) *ysize;
                 float deltay = y0-yhit[n];
                 //printf("2piy y: %.1f %.1f \n", y0, yhit[n]);
                 dytwo += deltay;
                 sytwo += deltay*deltay;
             }
             if(yw2 ==1){
                 nytwo++;
                 float y0 = (yc2*2 +2 - float(TYSIZE)/2) *ysize;
                 float deltay = y0-yhit[n];
                 //printf("2piy y: %.1f %.1f \n", y0, yhit[n]);
                 dytwo += deltay;
                 sytwo += deltay*deltay;
             }

         }


         //compute avgs and std devs of single pixel residuals if there are enough events
         if(nyone <= 10){
             dyone=lorbsy;
             syone = ysize/sqrt(12);
         }
         else{
             dyone /= float(nyone);
             syone = syone/float(nyone) - dyone*dyone;
             if(syone < 0.) syone = 0.;
             syone = sqrt(syone);
         }

         if(nytwo <= 10){
             dytwo=lorbsy;
             sytwo = 2.*ysize/sqrt(12);
         }
         else{
             dytwo /= float(nytwo);
             sytwo = sytwo/float(nytwo) - dytwo*dytwo;
             if(sytwo < 0.) sytwo = 0.;
             sytwo = sqrt(sytwo);
         }

         if(nxone <= 10){
             dxone=lorbsx;
             sxone = 2.*xsize/sqrt(12);
         }
         else{
             dxone /= float(nxone);
             sxone = sxone/float(nxone) - dxone*dxone;
             if(sxone < 0.) sxone = 0.;
             sxone = sqrt(sxone);
         }

         if(nxtwo <= 10){
             dxtwo=lorbsx;
             sxtwo = xsize/sqrt(12);
         }
         else{
             dxtwo /= float(nxtwo);
             sxtwo = sxtwo/float(nxtwo) - dxtwo*dxtwo;
             if(sxtwo < 0.) sxtwo = 0.;
             sxtwo = sqrt(sxtwo);
         }

         printf("%i Single pixel y-clusters avg offset: %.1f, std dev %.1f \n", nyone, dyone, syone);
         printf("%i Single big pixel y-clusters  avg offset: %.1f, std dev %.1f \n", nytwo, dytwo, sytwo);
         printf("%i Single pixel x-clusters avg offset: %.1f, std dev %.1f \n", nxone, dxone, sxone);
         printf("%i Single big pixel x-clusters  avg offset: %.1f, std dev %.1f \n", nxtwo, dxtwo, sxtwo);



         // Copy info into the slice and reformat from pixelav coordinates to CMSSW local coordinates
         SiPixelTemplateEntry * slice = new SiPixelTemplateEntry;


         slice->runnum = ifile;

         for(int i = 0; i < 2; ++i) {
             for(int j=0; j<5; ++j) {
                 slice->xpar[i][j] = xpar[i][j];
                 slice->ypar[i][j] = ypar[i][j];
             }
         }

         slice->clslenx = clslnx;
         slice->clsleny = clslny;




         slice->sxone = sxone;
         slice->dxone = dxone;
         slice->sxtwo = sxtwo;
         slice->dxtwo = dxtwo;

         slice->syone = syone;
         slice->dyone = dyone;
         slice->sytwo = sytwo;
         slice->dytwo = dytwo;



         slice->costrk[0] = -cosx;
         slice->costrk[1] = -cosy;
         slice->costrk[2] = -cosz;
         slice->cotalpha = cosx/cosz;
         slice->cotbeta = cosy/cosz;

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
                 slice->xtemp[9-k-1][TXSIZE - 1 - i] = xtemp[k][i];
             }
             for(int j=0; j<TYSIZE; j++){
                 slice->ytemp[9-k-1][TYSIZE - 1 - j] = ytemp[k][j];
             }
         }


         locBx = 1.;
         if(cotbeta < 0.) locBx = -1.;
         locBz = locBx;
         if(cotalpha < 0.) locBz = -locBx;

         bool flip_x = false;
         bool flip_y = false;

         if(locBx*locBz < 0.f) {
             flip_x = true;
         }
         if(locBx < 0.f) {
             flip_y = true;
         }

         templ.sideload(slice, IDtype, locBx, locBz, lorwdy, lorwdx, q50, fbin, xsize, ysize, thick);


         nbad = 0;
         ngood = 0;
         tote = 0.;
         bade = 0.;

         if(write_temp_header && ifile==startfile) {


             fprintf(temp_output_file,"%s", header);
             fprintf(temp_output_file,"%d %d %4.2f %d %d %d %d %5.4f %5.4f %4.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %f %4.4f %4.4f %4.4f \n",
                     id,nvers,Bfield,NTy,NTyx,NTxx,IDtype,Vbias, temp,fluenc,qscale,q50,lorwdy,
                     lorwdx,ysize,xsize,thick,q51,lorbsy,lorbsx,fbin[0], fbin[1], fbin[2]);

             fprintf(generr_output_file,"%s", header);
             fprintf(temp_output_file,"%d %d %4.2f %d %d %d %d %5.4f %5.4f %4.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %f %4.4f %4.4f %4.4f \n",
                     id,nvers,Bfield,NTy,NTyx,NTxx,IDtype,Vbias, temp,fluenc,qscale,q50,lorwdy,
                     lorwdx,ysize,xsize,thick,q51,lorbsy,lorbsx,fbin[0], fbin[1], fbin[2]);
         }

         for(int i=0; i<hp.size(); ++i) { 
             if(hp[i] == NULL) continue;
             hp[i]->Reset();
         }

         qperbit = pixmax/(pow(2.,(double)(numbits)-1.));



         // Loop over all clusters and apply generic and 1d reco
         for(int n=0; n<read_events; n++){



             // Do generic reco on the cluster

             float xsum[TXSIZE], ysum[TYSIZE];

             memset(xsum, 0., sizeof(xsum));
             memset(ysum, 0., sizeof(ysum));
             for(int j=0; j<TXSIZE; j++){
                 for(int i=0; i<TYSIZE; i++){
                     xsum[j] += cluster[n][j][i];
                     ysum[i] += cluster[n][j][i];
                 }
             }

             int xwidth=0;
             int ywidth=0;
             int xstart=0;
             int ystart=0;
             int xend=0;
             int yend=0;

             for(int j=0; j<TXSIZE; j++){
                 if(xsum[j] >0.){
                     if(xstart==0) xstart = j;
                     xwidth++;
                     xend = j;
                 }
             }
             for(int i=0; i<TYSIZE; i++){
                 if(ysum[i] >0.){
                     if(ystart==0) ystart = i;
                     ywidth++;
                     yend = i;
                 }
             }

             //charges of first and last
             float Q_f_x = xsum[xstart];
             float Q_l_x = xsum[xend];
             float Q_f_y = ysum[ystart];
             float Q_l_y = ysum[yend];

             //edges of cluster
             //
             //f is lower edge of first pixel
             //l is upper edge of last pixel
             float e_f_x  = xstart *xsize;
             float e_l_x  = (xend+1)*xsize; 
             float e_f_y  = ystart *ysize;
             float e_l_y  = (yend+1) *ysize;

             int size_x = xend-xstart + 1;
             int size_y = yend-ystart + 1;

             bool isBigPix = false;

             //taken from CPEGeneric config as of June 2019
             float eff_charge_cut_lowX = 0.;
             float eff_charge_cut_lowY = 0.;
             float eff_charge_cut_highX = 1.0;
             float eff_charge_cut_highY = 1.0;
             float size_cutX = 3.0;
             float size_cutY = 3.0;


             float xrec_gen = SiPixelUtils::generic_position_formula(size_x, Q_f_x, Q_l_x, e_f_x, e_l_x,
                     lorbsx, thick, cotalpha,
                     xsize, isBigPix, isBigPix,
                     eff_charge_cut_lowX, eff_charge_cut_highX,
                     size_cutX);
             float yrec_gen = SiPixelUtils::generic_position_formula(size_y, Q_f_y, Q_l_y, e_f_y, e_l_y,
                     lorbsy, thick, cotbeta,
                     ysize, isBigPix, isBigPix,
                     eff_charge_cut_lowY, eff_charge_cut_highY,
                     size_cutY);

             float dx_gen = xrec - (TXSIZE/2)*xsize - xhit[n];
             float dy_gen = yrec - (TYSIZE/2)*ysize - yhit[n];




             // Do the template reco on the cluster 


             // No double pixels 
            bool ydouble[TYSIZE], xdouble[TXSIZE];
            memset(ydouble, false, sizeof(ydouble));
            memset(xdouble, false, sizeof(xdouble));



             float cluster_local[TXSIZE][TYSIZE];
             memset(cluster_local, 0., sizeof(cluster_local));
             for(int i=0; i<TXSIZE; i++){
                 for(int j=0; j<TYSIZE; j++){
                     cluster_local[i][j] = cluster[n][i][j];
                 }
             }


            int mrow = TXSIZE, mcol = TYSIZE;


             //        if(fabs(cotbeta) < 2.1) continue;
             // Do the template analysis on the cluster 
             SiPixelTemplateReco::ClusMatrix clusterPayload{&cluster_local[0][0], xdouble, ydouble, mrow,mcol};


             int speed = 0;
             int qbin;


             int ierr = PixelTempReco1D(tempID, cotalpha, cotbeta, locBz, locBx,  clusterPayload, templ, yrec, sigmay, proby, xrec, sigmax, probx,  qbin, speed, probQ);
             if(ierr != 0) {
                 ++nbad; 
                 printf("reconstruction failed with error %d \n", ierr);
             } else {
                 ngood++;
                 nqbin[qbin] += 1.;
                 dy = yrec - (TYSIZE/2)*ysize - yhit[n];
                 dx = xrec - (TXSIZE/2)*xsize - xhit[n];


                 //Fill template reco and generic reco residuals
                 hp[y_temp_idx]->Fill(dy);
                 hp[y_temp_idx+1+qbin]->Fill(dy);

                 hp[x_temp_idx]->Fill(dx);
                 hp[x_temp_idx+1+qbin]->Fill(dx);

                 hp[y_generic_idx]->Fill(dy_gen);
                 hp[y_generic_idx+1 +qbin]->Fill(dy_gen);

                 hp[x_generic_idx]->Fill(dx_gen);
                 hp[x_generic_idx+1 +qbin]->Fill(dx_gen);


                 //fill npix qbins here
                 hp[charge_idx+2+qbin]->Fill(npix[n]);

                 // fill Chisq from template fit here


                 /*
                 int k= int(xhit[n]/xsize * 8. + 4.5);
                 int l= int(yhit[n]/ysize * 8. + 4.5);
                 printf("Hit bins %i %i \n", k,l);
                 printf("template dx dy %.3f %.3f \n", dx,dy);
                 printf("generic dx dy %.3f %.3f \n", dx_gen,dy_gen);
                 */

             }

         }

         printf(" low q failures = %.0f, failed fits = %d, successful fits = %d, total read events %d \n", nqbin[4], nbad, ngood, read_events);	   

         /*
          * Histograms plotting
          */


         for(int i=0; i<5; ++i) {
             hp[y_temp_idx + i]->Fit("gaus");
             hp[x_temp_idx + i]->Fit("gaus");
             hp[y_generic_idx + i]->Fit("gaus");
             hp[x_generic_idx + i]->Fit("gaus");

         }


         /*
            float scaley, scalex, scalx[4], scaly[4], delyavg, delysig, offsetx[4], offsety[4];
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
         */

         //  Create an output filename for this run 

         sprintf(outfile0,"template_histos%5.5d_newResp.ps[",ifile);
         sprintf(outfile1,"template_histos%5.5d_newResp.ps",ifile);
         sprintf(outfile2,"template_histos%5.5d_newResp.ps]",ifile);
         c1->Clear();
         c1->Print(outfile0);
         for(int i=0; i<hp.size(); ++i) {
             if(hp[i] == NULL) continue;
             hp[i]->Draw();
             c1->Print(outfile1);
             c1->Clear();
         }
         c1->Print(outfile2);
         c1->Clear();



         // Write this template entry to the output file
         //
         fprintf(temp_output_file, "%i %8.6f %8.6f %8.6f \n", ifile, slice->costrk[0], slice->costrk[1], slice->costrk[2]);

         fprintf(temp_output_file, "%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f \n",
                 qavg, pixmax, symax, dyone, syone, sxmax, dyone, syone);
         //grab 30th smallest charge
         auto it = std::next(qmsort.begin(), 29);
         float qmin30 = *it;
         //grab 60th smallest charge
         it = std::next(qmsort.begin(), 59);
         float qmin60 = *it;

         fprintf(temp_output_file, "%.2f %.2f %.2f %.2f %.2f %.2f %.2f \n",
                 dytwo, sytwo, dxtwo, sytwo, qmin30, clslny, clslnx );

         //y charge variance fit
         for(int i=0; i<2; i++){
             for(int j=0; j<5; j++){
                 fprintf(temp_output_file, "%15.8E ", slice->ypar[i][j]);
             }
             fprintf(temp_output_file, "\n");
         }
         //y template
         for(int k = 0; k < 9; k++){
             for(int j=0; j<TYSIZE; j++){
                 fprintf(temp_output_file, "%8.1f ", slice->ytemp[k][j]);
             }
             fprintf(temp_output_file, "\n");
         }

         //x charge variance fit
         for(int i=0; i<2; i++){
             for(int j=0; j<5; j++){
                 fprintf(temp_output_file, "%15.8E ", slice->ypar[i][j]);
             }
             fprintf(temp_output_file, "\n");
         }
         //x template
         for(int k = 0; k < 9; k++){
             for(int j=0; j<TYSIZE; j++){
                 fprintf(temp_output_file, "%8.1f ", slice->ytemp[k][j]);
             }
             fprintf(temp_output_file, "\n");
         }


         //fitted stuff
         // MISSING




         float tqbin = nqbin[0] + nqbin[1] + nqbin[2] + nqbin[3]; //exclude low-q failures
         for(int i=0; i<4; i++){
             nqbin[i]/tqbin;
         }

         float fyone = nyone/nevents;
         float fxone = nyone/nevents;
         float fytwo = nyone/2./nevents;
         float fxtwo = nyone/2./nevents;

         //output to gen_errors

         fprintf(generr_output_file, "%i %8.6f %8.6f %8.6f \n", ifile, slice->costrk[0], slice->costrk[1], slice->costrk[2]);

         fprintf(temp_output_file, "%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f \n",
                 qavg, pixmax, symax, dyone, syone, sxmax, dyone, syone);

         fprintf(temp_output_file, "%.2f %.2f %.2f %.2f %.2f %.2f \n",
                 dytwo, sytwo, dxtwo, sytwo, qmin30, qmin60);


         //more fitted stuff


    }
    // Close output files

    fclose(temp_output_file);  
    fclose(generr_output_file);  

    /*  Determine current time */

    gettimeofday(&now1, &timz);
    deltas = now1.tv_sec - now0.tv_sec;
    deltaus = now1.tv_usec - now0.tv_usec;
    deltat = ((double)deltaus)/1000000.;
    deltat += (double)deltas;
    printf("ellapsed time = %f seconds \n", deltat);

    delete_3d_array(cluster, nevents, TXSIZE, TYSIZE);
    delete_2d_array(xsum1, nevents, TXSIZE);
    delete_2d_array(xsum2, nevents, TXSIZE);
    delete_2d_array(ysum1, nevents, TYSIZE);
    delete_2d_array(ysum2, nevents, TYSIZE);

    return 0;
} // MAIN__ 





