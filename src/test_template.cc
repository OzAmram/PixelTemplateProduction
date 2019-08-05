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
//! Add clustering  algorithm
//! Change response function, parameters to that used for real analysis
//! Same as test_code_v9 but with weighting for flat eta distributions in the resolution histograms
//! Add 2-threshold simulation
//! Work on templates with non-zero scale factors, ask for template IDs to use
//! Version to calibrate the LAs for CPEGeneric
//! Add charge centroid information 
//! Add charge collection efficiency information


#define TEMPL_DEBUG
#include "template_utils.h"
#include "PixelGeneric2D.cc"

// Global definitions 

// Main program  

int main(int argc, char *argv[])
{

    // Local variables 
    std::vector<float> pvec(6), wgauss(TYSIZE), vgauss(TYSIZE), xgauss(TYSIZE), ygauss(TYSIZE), zgauss(TYSIZE);
    float pixin[TXSIZE][TYSIZE];
    static bool fpix, noqbin4;
    bool ydouble[TYSIZE], xdouble[TXSIZE];
    int mrow = TXSIZE, mcol = TYSIZE;
    static float thick, xsize, ysize, noise, zcen, xcmssw, ycmssw, sxcmssw, sycmssw, gain_frac, q100_frac, common_frac, readout_noise, qscale, thrseed;
    static float xhit, yhit, xrec, yrec, sigmax, sigmay, probx, proby, signal, cotalpha, cotbeta, qclust, locBz, locBx, probxy, probQ, xclust, yclust;  
    static int ndata, nfile, neh, nevent, ID, nbad, non_linear, icol, ndcol, layer, nclust; 
    static vector<int> nbin(5,0);
    float zvtx, zdet;
    int izdet;
    static vector<double> fluence(4,0.);
    int i, j, k, ierr, qbin, qb, nypix, nxpix, etabin, jmin, jmax, imin, imax, numadd, idcol;
    double dx, dy, eta, dxc, dyc, log10probxy, log10probQ, tote, bade, weight, alpha, qnorm, dxclust, dyclust;
    double qpixels, qfrac;
    static int iyd, ixd, speed;	
    static float q100, q101, q50, q10, qmax; 

    float sigtmp, qin, yfrac, xfrac;
    static char infile[80], header[80], c, outfile0[80], outfile1[80], outfile2[80];
    int triplg(std::vector<float>&);
    //	int random(void);
    float cluster[TXSIZE][TYSIZE], clust[TXSIZE][TYSIZE], rclust[TXSIZE][TYSIZE];
    bool bclust[TXSIZE][TYSIZE];

    std::pair<int, int> pixel, max;

    FILE *ifp;

    struct timeval now0, now1;
    struct timezone timz;
    long deltas, deltaus;
    double deltat;


    //  Read which data and inputs to use (use c file i/o which is much faster than c++ i/o) 

    ifp = fopen("q_dist_2t.txt", "r");
    if (ifp==NULL) {
        printf("no pixel initialization file found \n");
        return 0;
    }


    fscanf(ifp,"%d %d %f %f %f %f %f %f %f %d", &ndata, &nfile, &noise, &q100, &q101, &q100_frac, &common_frac, &gain_frac, &readout_noise, &non_linear);
    fclose(ifp);
    printf("data file %d, mc file %d noise = %f, threshold0 = %f, threshold1 = %f, rms threshold frac = %f, common_frac = %f, gain fraction = %f, readout noise = %f, nonlinear_resp = %d \n", ndata, nfile, noise, q100, q101, q100_frac, common_frac, gain_frac, readout_noise, non_linear);

    FrontEndModel frontEnd;
    frontEnd.fe_type       = non_linear;
    frontEnd.gain_frac     = gain_frac;
    frontEnd.readout_noise = readout_noise;

    //  Calculate 50% of threshold in q units and enc noise in adc units

    q50=0.5*q100;
    q10=0.2*q50;

    //  Create an input filename for this run 


    sprintf(infile,"template_events_d%i.out",nfile);

    //  Open input file and read header info 

    ifp = fopen(infile, "r");
    if (ifp==NULL) {
        printf("Can't find pixel data file %s  \n", infile);
        return 0;
    }

    // Read-in a header string first and print it    

    for (i=0; (c=getc(ifp)) != '\n'; ++i) {
        if(i < 79) {header[i] = c;}
    }
    if(i > 78) {i=78;}
    header[i+1] ='\0';
    printf("%s\n", header);

    double  halfxs=300.;
    int nx=120;	
    gStyle->SetOptFit(101);
    gStyle->SetHistLineWidth(2);
    static vector<TH1F*> hp(41);
    hp[0] = new TH1F("h101","dy_temp (all sig); #Deltay (#mum)",nx,-halfxs,halfxs);
    hp[1] = new TH1F("h102","dy_temp (signal @> 1.5mn); #Deltay (#mum)",nx,-halfxs,halfxs);      
    hp[2] = new TH1F("h103","dy_temp (1.5mn @> signal @> 1.0mn); #Deltay (#mum)",nx,-halfxs,halfxs);      
    hp[3] = new TH1F("h104","dy_temp (1.0mn @> signal @> 0.85mn); #Deltay (#mum)",nx,-halfxs,halfxs);     
    hp[4] = new TH1F("h105","dy_temp (0.85mn @> signal); #Deltay (#mum)",nx,-halfxs,halfxs);      
    hp[5] = new TH1F("h201","log10(probxy) (all sig)",nx,-12.,0.);
    hp[6] = new TH1F("h202","log10(probxy) (signal @> 1.5mn)",nx,-12.,0.);      
    hp[7] = new TH1F("h203","log10(probxy) (1.5mn @> signal @> 1.0mn)",nx,-12.,0.);      
    hp[8] = new TH1F("h204","log10(probxy) (1.0mn @> signal @> 0.85mn)",nx,-12.,0.);     
    hp[9] = new TH1F("h205","log10(probxy) (0.85mn @> signal)",nx,-12.,0.);      
    hp[10] = new TH1F("h106","dx_temp (all sig); #Deltax (#mum)",nx,-halfxs,halfxs);
    hp[11] = new TH1F("h107","dx_temp (signal @> 1.5mn); #Deltax (#mum)",nx,-halfxs,halfxs);      
    hp[12] = new TH1F("h108","dx_temp (1.5mn @> signal @> 1.0mn); #Deltax (#mum)",nx,-halfxs,halfxs);      
    hp[13] = new TH1F("h109","dx_temp (1.0mn @> signal @> 0.85mn); #Deltax (#mum)",nx,-halfxs,halfxs);      
    hp[14] = new TH1F("h110","dx_temp (0.85mn @> signal); #Deltax (#mum)",nx,-halfxs,halfxs);    
    hp[15] = new TH1F("h206","log10(probQ) (all sig)",nx,-12.,0.);
    hp[16] = new TH1F("h207","log10(probQ) (signal @> 1.5mn)",nx,-12.,0.);      
    hp[17] = new TH1F("h208","log10(probQ) (1.5mn @> signal @> 1.0mn)",nx,-12.,0.);      
    hp[18] = new TH1F("h209","log10(probQ) (1.0mn @> signal @> 0.85mn)",nx,-12.,0.);     
    hp[19] = new TH1F("h210","log10(probQ) (0.85mn @> signal)",nx,-12.,0.);      
    hp[20] = new TH1F("h300","cotbeta (probxy<10-2)",nx,-10.,10.);      
    hp[21] = new TH1F("h301","cotalpha (probxy<10-2)",nx,-0.25,0.25);   
    hp[22] = new TH1F("h401","dy_cmssw (all sig); #Deltay (#mum)",nx,-halfxs,halfxs);
    hp[23] = new TH1F("h406","dx_cmssw (all sig); #Deltax (#mum)",nx,-halfxs,halfxs);
    hp[24] = new TH1F("h111","dy_temp (one pixel clust); #Deltay (#mum)",nx,-halfxs,halfxs);
    hp[25] = new TH1F("h112","dx_temp (one pixel clust); #Deltax (#mum)",nx,-halfxs,halfxs);
    hp[26] = new TH1F("h501","xy Probability; probxy",101,0,1.01);
    hp[27] = new TH1F("h502","Q Probability; probQ",101,0,1.01);
    hp[28] = new TH1F("h302","cotbeta (probQ<10-2)",nx,-10.,10.);      
    hp[29] = new TH1F("h303","cotalpha (probQ<10-2)",nx,-0.25,0.25);   
    hp[30] = new TH1F("h113","dy_temp (>one pixel clust); #Deltay (#mum)",nx,-halfxs,halfxs);
    hp[31] = new TH1F("h114","dx_temp (>one pixel clust); #Deltax (#mum)",nx,-halfxs,halfxs);
    hp[32] = new TH1F("h115","Cluster Charge; Q_clus (e)",nx,0.,120000.);
    hp[33] = new TH1F("h116","Normalized Cluster Charge; Q_clus (e)",nx,0.,60000.);
    hp[34] = new TH1F("h117","dy_temp (cot#beta > 0); #Deltay (#mum)",nx,-halfxs,halfxs);
    hp[35] = new TH1F("h118","dy_temp (cot#beta < 0); #Deltay (#mum)",nx,-halfxs,halfxs);
    hp[36] = new TH1F("h119","dy_cmssw (>one pixel clust); #Deltay (#mum)",nx,-halfxs,halfxs);
    hp[37] = new TH1F("h120","dx_cmssw (>one pixel clust); #Deltax (#mum)",nx,-halfxs,halfxs);
    hp[38] = new TH1F("h411","dy_baryc (all sig); #Deltay (#mum)",nx,-halfxs,halfxs);
    hp[39] = new TH1F("h412","dx_baryc (all sig); #Deltax (#mum)",nx,-halfxs,halfxs);
    hp[40] = new TH1F("h121","CCE; CCE",110,0.,1.10);


    // Set style for the the histograms	

    for(i=0; i<41; ++i) {
        hp[i]->SetLineColor(2);
        hp[i]->SetFillColor(38);
    }

    // Make some profile Histograms

    static vector<TProfile*> pp(53);
    pp[0] = new TProfile("pqhytmp","1.5>Q/Q_avg>1, y_temp",12,0,3.00,"s"); 
    pp[1] = new TProfile("pqhycms","1.5>Q/Q_avg>1, y_cmssw",12,0,3.00,"s"); 
    pp[2] = new TProfile("pqlytmp","1>Q/Q_avg, y_temp",12,0,3.00,"s"); 
    pp[3] = new TProfile("pqlycms","1>Q/Q_avg, y_cmssw",12,0,3.00,"s"); 
    pp[4] = new TProfile("pqhxtmp","1.5>Q/Q_avg>1, x_temp",12,0,3.00,"s"); 
    pp[5] = new TProfile("pqhxcms","1.5>Q/Q_avg>1, x_cmssw",12,0,3.00,"s"); 
    pp[6] = new TProfile("pqlxtmp","1>Q/Q_avg, x_temp",12,0,3.00,"s"); 
    pp[7] = new TProfile("pqlxcms","1>Q/Q_avg, x_cmssw",12,0,3.00,"s"); 
    pp[8] = new TProfile("pycfrac2","dyc vs yfrac, 2 pix clust; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
    pp[9] = new TProfile("pycfrac3","dyc vs yfrac, 3 pix clust; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
    pp[10] = new TProfile("pycfrac4","dyc vs yfrac, 4 pix clust; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
    pp[11] = new TProfile("pxcfrac2","dxc vs xfrac, 2 pix clust; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
    pp[12] = new TProfile("pxcfrac3","dxc vs xfrac, 3 pix clust; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
    pp[13] = new TProfile("pxcfrac4","dxc vs xfrac, 4 pix clust; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
    pp[14] = new TProfile("pyfrac2","dy vs yfrac, 2 pix clust; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
    pp[15] = new TProfile("pyfrac3","dy vs yfrac, 3 pix clust; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
    pp[16] = new TProfile("pyfrac4","dy vs yfrac, 4 pix clust; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
    pp[17] = new TProfile("pxfrac2","dx vs xfrac, 2 pix clust; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
    pp[18] = new TProfile("pxfrac3","dx vs xfrac, 3 pix clust; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
    pp[19] = new TProfile("pxfrac4","dx vs xfrac, 4 pix clust; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
    pp[20] = new TProfile("pycfrac2q1","dyc vs yfrac, 2 pix clust, 1.5>Q/Q_avg>1; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
    pp[21] = new TProfile("pycfrac3q1","dyc vs yfrac, 3 pix clust, 1.5>Q/Q_avg>1; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
    pp[22] = new TProfile("pycfrac4q1","dyc vs yfrac, 4 pix clust, 1.5>Q/Q_avg>1; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
    pp[23] = new TProfile("pxcfrac2q1","dxc vs xfrac, 2 pix clust, 1.5>Q/Q_avg>1; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
    pp[24] = new TProfile("pxcfrac3q1","dxc vs xfrac, 3 pix clust, 1.5>Q/Q_avg>1; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
    pp[25] = new TProfile("pxcfrac4q1","dxc vs xfrac, 4 pix clust, 1.5>Q/Q_avg>1; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
    pp[26] = new TProfile("pyfrac2q1","dy vs yfrac, 2 pix clust, 1.5>Q/Q_avg>1; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
    pp[27] = new TProfile("pyfrac3q1","dy vs yfrac, 3 pix clust, 1.5>Q/Q_avg>1; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
    pp[28] = new TProfile("pyfrac4q1","dy vs yfrac, 4 pix clust, 1.5>Q/Q_avg>1; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
    pp[29] = new TProfile("pxfrac2q1","dx vs xfrac, 2 pix clust, 1.5>Q/Q_avg>1; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
    pp[30] = new TProfile("pxfrac3q1","dx vs xfrac, 3 pix clust, 1.5>Q/Q_avg>1; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
    pp[31] = new TProfile("pxfrac4q1","dx vs xfrac, 4 pix clust, 1.5>Q/Q_avg>1; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
    pp[32] = new TProfile("pycfrac2q2","dyc vs yfrac, 2 pix clust, 1>Q/Q_avg; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
    pp[33] = new TProfile("pycfrac3q2","dyc vs yfrac, 3 pix clust, 1>Q/Q_avg; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
    pp[34] = new TProfile("pycfrac4q2","dyc vs yfrac, 4 pix clust, 1>Q/Q_avg; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
    pp[35] = new TProfile("pxcfrac2q2","dxc vs xfrac, 2 pix clust, 1>Q/Q_avg; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
    pp[36] = new TProfile("pxcfrac3q2","dxc vs xfrac, 3 pix clust, 1>Q/Q_avg; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
    pp[37] = new TProfile("pxcfrac4q2","dxc vs xfrac, 4 pix clust, 1>Q/Q_avg; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
    pp[38] = new TProfile("pyfrac2q2","dy vs yfrac, 2 pix clust, 1>Q/Q_avg; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
    pp[39] = new TProfile("pyfrac3q2","dy vs yfrac, 3 pix clust, 1>Q/Q_avg; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
    pp[40] = new TProfile("pyfrac4q2","dy vs yfrac, 4 pix clust, 1>Q/Q_avg; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
    pp[41] = new TProfile("pxfrac2q2","dx vs xfrac, 2 pix clust, 1>Q/Q_avg; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
    pp[42] = new TProfile("pxfrac3q2","dx vs xfrac, 3 pix clust, 1>Q/Q_avg; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
    pp[43] = new TProfile("pxfrac4q2","dx vs xfrac, 4 pix clust, 1>Q/Q_avg; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
    pp[44] = new TProfile("pqhxtmpvsa","1.5>Q/Q_avg>1, x_temp",25,-1.125,1.125,"s"); 
    pp[45] = new TProfile("pqhxcmsvsa","1.5>Q/Q_avg>1, x_cmssw",25,-1.125,1.125,"s"); 
    pp[46] = new TProfile("pqlxtmpvsa","1>Q/Q_avg, x_temp",25,-1.125,1.125,"s"); 
    pp[47] = new TProfile("pqlxcmsvsa","1>Q/Q_avg, x_cmssw",25,-1.125,1.125,"s"); 
    pp[48] = new TProfile("datnxcotb","nx(pix); cot(#beta)",50,-4.5,4.5," "); 
    pp[49] = new TProfile("ybarycvscotb","dy-centroid; cot(#beta)",18,-4.5,4.5," "); 
    pp[50] = new TProfile("ybarycvscota","dy-centroid; cot(#alpha)",10,-0.5,0.5," "); 
    pp[51] = new TProfile("xbarycvscotb","dx-centroid; cot(#beta)",18,-4.5,4.5," "); 
    pp[52] = new TProfile("xbarycvscota","dx-centroid; cot(#alpha)",10,-0.3,0.3," "); 



    // Set style for the the profile histograms	

    for(i=0; i<6; ++i) {
        pp[i+8]->SetLineColor(2);
        pp[i+8]->SetStats(kFALSE);
        pp[i+20]->SetLineColor(2);
        pp[i+20]->SetStats(kFALSE);
        pp[i+32]->SetLineColor(2);
        pp[i+32]->SetStats(kFALSE);
    }
    for(i=6; i<12; ++i) {
        pp[i+8]->SetLineColor(4);
        pp[i+8]->SetStats(kFALSE);
        pp[i+20]->SetLineColor(4);
        pp[i+20]->SetStats(kFALSE);
        pp[i+32]->SetLineColor(4);
        pp[i+32]->SetStats(kFALSE);
    }
    pp[48]->SetLineColor(2);
    pp[48]->SetLineStyle(2);


    fscanf(ifp,"%f  %f  %f", &ysize, &xsize, &thick);
    zcen = thick/2.;
    printf("xsize/ysize/thick = %f/%f/%f \n", xsize, ysize, thick);
    fpix = false;
    if(thick > 286.) {fpix = true;}

    nevent=0;
    nclust=0;
    nbad = 0;
    tote = 0.;
    bade = 0.;

    static vector<float> sx(5,0.), sx2(5,0.), scx(5,0.), sy(5,0.), sy2(5,0.), scy(5,0.); 

    static vector<float> sxp(5,0.), sxp2(5,0.), scxp(5,0.), syp(5,0.), syp2(5,0.), scyp(5,0.); 

    static vector<float> sxc(5,0.), scxc(5,0.), sxc2(5,0.), syc(5,0.), scyc(5,0.), syc2(5,0.), nt(12,0.), nb(12,0.); 

    std::vector<std::pair<int, int> > pixlst;


    printf("Enter > 0 for no qbin = 4 clusters in histos \n");
    scanf("%d", &i);
    noqbin4 = false;
    if(i > 0) noqbin4 = true;

    // Decide if this file corresponds to a CMSSW run 

    printf("Enter ID of Template and layer number [=1 changes vcal] \n");
    scanf("%d %d", &ID, &layer);

    if(layer == 1) {
      printf("using L1 parameters \n");
      frontEnd.vcal = 50.;
      frontEnd.vcaloffst = 670.;
    } 


    if(fpix) {printf("ID = %d, fpix \n", ID);} else {printf("ID = %d, barrel, vcal = %lf, offset = %lf \n", ID, frontEnd.vcal, frontEnd.vcaloffst);}


    // Initialize template store 

    std::vector< SiPixelTemplateStore > thePixelTemp_;
    SiPixelTemplate templ(thePixelTemp_);

    // Initialize template store, Pixelav 100V/300V simulation, +20C as thePixelTemp[6] 

    templ.pushfile(ID,thePixelTemp_);
    templ.interpolate(ID, 0.f, 0.f, -1.f);
    qscale = templ.qscale();	

    // Initialize GenError store 

    std::vector< SiPixelGenErrorStore > thePixelGenErr_;
    SiPixelGenError gtempl(thePixelGenErr_);

    gtempl.pushfile(ID,thePixelGenErr_);	


    // Ask for speed info

    printf("enter seeding threshold (4000) \n");
    scanf("%f", &thrseed);

    printf("seeding threshold = %f\n", thrseed);

    speed = -2;

    printf("speed = %d\n", speed);

    // Ask for external input into which pixels to join 

    //	printf("enter the first y- and x-columns to double (-1 -1 = none)\n");
    //	scanf("%d %d", &iyd, &ixd);

    iyd = -1; ixd = -1;
    printf("y/x double = %d/%d\n", iyd,ixd);

    int totale=0;  int goodx=0; int goody=0;

    ndcol = TYSIZE/2 + 1;
    std::vector<int> ndhit(ndcol, 0);

    /*  Determine current time */

    gettimeofday(&now0, &timz);

    // Loop until end of input file 

    while(fscanf(ifp,"%f %f %f %f %f %f %d", &pvec[0], &pvec[1], &pvec[2], &pvec[3], &pvec[4], &pvec[5], &neh) != EOF) {

        // read the input cluster 

        for (k=0; k < TXSIZE; ++k) {
            fscanf(ifp,
                    "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", 
                    &pixin[k][0],&pixin[k][1],&pixin[k][2],&pixin[k][3],&pixin[k][4],&pixin[k][5],&pixin[k][6],&pixin[k][7],&pixin[k][8],&pixin[k][9],
                    &pixin[k][10],&pixin[k][11],&pixin[k][12],&pixin[k][13],&pixin[k][14],&pixin[k][15],&pixin[k][16],&pixin[k][17],&pixin[k][18], 
                    &pixin[k][19],&pixin[k][20]);
        }
        ++nevent;

        // Add noise and analog response to cluster, reformat for flipped barrel coordinate system 

        triplg(vgauss);
        pixlst.resize(0);
        for(i=0; i<ndcol; ++i) {ndhit[i] = 0;}
        icol = 0;
        if(vgauss[1] < 0.) {icol = 1;}
        qpixels = 0.;
        for(j=0; j<TXSIZE; ++j) {
            triplg(wgauss);
            triplg(xgauss);
            triplg(ygauss);
            triplg(zgauss);
            for(i=0; i<TYSIZE; ++i) {
                bclust[j][i] = false;
                qin = (10.*pixin[j][i] + xgauss[i]*noise);
                rclust[TXSIZE-1-j][TYSIZE-1-i] = qin;
                if(qin < q100*(1.+wgauss[i]*q100_frac)) {
                    clust[TXSIZE-1-j][TYSIZE-1-i] = 0.;
                } else {
                    qpixels += qin;
                    idcol = (TYSIZE-1-i+icol)/2;
                    ++ndhit[idcol];
		    signal = frontEnd.apply_model( qin, ygauss[i], zgauss[i] );
                    clust[TXSIZE-1-j][TYSIZE-1-i] = (1.+vgauss[0]*common_frac)*signal;
                }
            }
        }



        // Simulate the second, higher threshold in single dcol hits

        for(j=0; j<TXSIZE; ++j) {
            for(i=0; i<TYSIZE; ++i) {
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

        memset(cluster, 0., sizeof(cluster));

        qclust=0.;
        xclust=0.;
        yclust=0.;
        for ( auto pixIter = pixlst.begin(); pixIter != pixlst.end(); ++pixIter ) {
            j = pixIter->first; 
            i = pixIter->second;
            cluster[j][i] = clust[j][i];
            qclust += clust[j][i];
            xclust += ((float)j)*xsize*clust[j][i];
            yclust += ((float)i)*ysize*clust[j][i];
        }
        xclust /= qclust;
        yclust /= qclust;
        hp[32]->Fill((double)qclust);

        // Calculate the hit coordinates in the flipped coordinate system 

        yhit = -(pvec[0] + (zcen-pvec[2])*pvec[3]/pvec[5]);
        xhit = -(pvec[1] + (zcen-pvec[2])*pvec[4]/pvec[5]);

        // Do the template analysis on the cluster 

        cotalpha = pvec[4]/pvec[5];
        alpha = atan((double)cotalpha);
        cotbeta = pvec[3]/pvec[5];
        qnorm = (double)(qclust/sqrt(1.+cotalpha*cotalpha+cotbeta*cotbeta));
        hp[33]->Fill(qnorm);

        eta = fabs(-log((double)(-cotbeta+sqrt((double)(1.+cotbeta*cotbeta)))));
        weight = 1./cosh((double)eta);
        tote += weight;
        etabin = (int)(eta/0.25);
        if(etabin > 11) {etabin = 11;}
        ++nt[etabin];
        zvtx = vgauss[2]*6.;
        zdet = abs(zvtx + cotbeta*4.);
        izdet = (int)(zdet/6.24);
        if(izdet < 4) {
            fluence[izdet] += weight*sqrt((double)(1.+cotbeta*cotbeta));
        }

        // Combine two single pixels into a double pixel 

        for(i=0; i<TYSIZE; ++i) {
            ydouble[i] = false;
        }
        if(iyd >= 0 && iyd < TYSIZE) {
            ydouble[iyd] = true;
            for(j=0; j<TXSIZE; ++j) {
                sigtmp = cluster[j][iyd]+cluster[j][iyd+1];
                cluster[j][iyd] = sigtmp;
            }
            for(i=iyd+1; i<TYSIZE-1; ++i) {
                for(j=0; j<TXSIZE; ++j) {
                    cluster[j][i] = cluster[j][i+1];
                }
            }    
            for(j=0; j<TXSIZE; ++j) {
                cluster[j][TYSIZE-1] = 0.;
            }
        }

        for(j=0; j<TXSIZE; ++j) {
            xdouble[j] = false;
        }
        if(ixd >= 0 && ixd < TXSIZE-1) {
            xdouble[ixd] = true;
            for(i=0; i<TYSIZE; ++i) {
                sigtmp = cluster[ixd][i]+cluster[ixd+1][i];
                cluster[ixd][i] = sigtmp;
            }
            for(j=ixd+1; j<TXSIZE-1; ++j) {
                for(i=0; i<TYSIZE; ++i) {
                    cluster[j][i] = cluster[j+1][i];
                }
            }    
            for(i=0; i<TYSIZE; ++i) {
                cluster[TXSIZE-1][i] = 0.;
            }
        }
        qfrac = (double)qclust/(double)neh;
        hp[40]->Fill(qfrac, weight);

        // Do the template analysis on the cluster 
        // Do the template analysis on the cluster 
        SiPixelTemplateReco::ClusMatrix clusterPayload{&cluster[0][0], xdouble, ydouble, mrow,mcol};
        SiPixelTemplateReco::ClusMatrix clusterPayloadC = clusterPayload;
        locBx = 1.;
        if(cotbeta < 0.) locBx = -1.;
        locBz = locBx;
        if(cotalpha < 0.) locBz = -locBx;
	std::cout << "cotalpha " << cotalpha << " cotbeta " << cotbeta << " Bz " << locBz << " Bx "<< locBx << std::endl;
	ierr = PixelTempReco1D(ID, cotalpha, cotbeta, locBz, locBx,  clusterPayload, templ, yrec, sigmay, proby, xrec, sigmax, probx,  qbin, speed, probQ);
	std::cout << " ierr " << ierr << std::endl;
        if(ierr != 0) {
            ++nbad; ++nb[etabin]; bade +=weight;
            printf("reconstruction failed with error %d \n", ierr);
        } else {
            qb = qbin;
            // Check resolution and weights
            if(qbin > 3) {
                qbin = 3; qb = 4;
                printf(" qbin = %d \n", qb);
            }
            ++nbin[qb];
            if(qb > 3 && noqbin4) continue;
            if(iyd != 0) {
                dy = yrec - (TYSIZE/2)*ysize - yhit;
                dyclust = yclust - (TYSIZE/2)*ysize - yhit;
            } else {
                dy = yrec - ((TYSIZE/2)-0.5)*ysize - yhit;
                dyclust = yclust - ((TYSIZE/2)-0.5)*ysize - yhit;
            }
            ++totale;
            probxy = probx*proby*(1.-log((double)(probx*proby)));
            if(probxy > 1.e-2) {++goodx;}else{hp[20]->Fill((double)cotbeta);hp[21]->Fill(alpha);}
            if(probQ < 1.e-2) {hp[28]->Fill((double)cotbeta);hp[29]->Fill(alpha);}
            if(probxy > 1.e-2) {++goody;}

            log10probxy = log10((double)probxy); log10probQ = log10((double)probQ);
            sy[qb] += dy; sy2[qb] += dy*dy; scy[qb] += dy*dy/(sigmay*sigmay);
            if(probxy > 1.e-2) {syp[qb] += dy; syp2[qb] += dy*dy; scyp[qb] += dy*dy/(sigmay*sigmay);}
            hp[0]->Fill(dy, weight);
            hp[1+qbin]->Fill(dy, weight);
            hp[38]->Fill(dyclust, weight);
            pp[49]->Fill((double)cotbeta, dyclust);
            pp[50]->Fill((double)cotalpha, dyclust);
            if(cotbeta > 0.) {hp[34]->Fill(dy, weight);} else {hp[35]->Fill(dy, weight);}
            hp[5]->Fill(log10probxy);
            hp[6+qbin]->Fill(log10probxy);
            if(qbin == 1) {pp[0]->Fill(eta, dy);}
            if(qbin > 1) {pp[2]->Fill(eta, dy);}
            if(ixd != 0) {
                dx = xrec - (TXSIZE/2)*xsize - xhit;
                dxclust = xclust - (TXSIZE/2)*xsize - xhit;
            } else {
                dx = xrec - ((TXSIZE/2)-0.5)*xsize - xhit;
                dxclust = xclust - ((TXSIZE/2)-0.5)*xsize - xhit;

            }
            sx[qb] += dx; sx2[qb] += dx*dx; scx[qb] += dx*dx/(sigmax*sigmax);
            if(probx > 1.e-2) {sxp[qb] += dx; sxp2[qb] += dx*dx; scxp[qb] += dx*dx/(sigmax*sigmax);}
            hp[10]->Fill(dx, weight);
            hp[11+qbin]->Fill(dx, weight);
            hp[39]->Fill(dxclust, weight);
            pp[51]->Fill((double)cotbeta, dxclust);
            pp[52]->Fill((double)cotalpha, dxclust);
            hp[15]->Fill(log10probQ);
            hp[16+qbin]->Fill(log10probQ);
            hp[26]->Fill((double)probxy);
            hp[27]->Fill((double)probQ);
            if(qbin == 1) {pp[4]->Fill(eta, dx);}
            if(qbin < 4 && qbin > 1) {pp[6]->Fill(eta, dx);}
            if(qbin == 1) {pp[44]->Fill(alpha, dx);}
            if(qbin < 4 && qbin > 1) {pp[46]->Fill(alpha, dx);}

            locBx = 1.;
            if(cotbeta < 0.) locBx = -1.;
            locBz = locBx;
            if(cotalpha < 0.) locBz = -locBx;            
            ierr = PixelGeneric2D(ID, cotalpha, cotbeta, locBz, locBx, clusterPayloadC, gtempl, ycmssw, sycmssw, xcmssw, sxcmssw, nypix, nxpix, yfrac, xfrac);
            if(iyd != 0) {
                dyc = ycmssw - (TYSIZE/2)*ysize - yhit;
            } else {
                dyc = ycmssw - ((TYSIZE/2)-0.5)*ysize - yhit;
            }
            syc[qb] += dyc; syc2[qb] += dyc*dyc; scyc[qb] += dyc*dyc/(sycmssw*sycmssw);
            if(ixd != 0) {
                dxc = xcmssw - (TXSIZE/2)*xsize - xhit;
            } else {
                dxc = xcmssw - ((TXSIZE/2)-0.5)*xsize - xhit;
            }
            sxc[qb] += dxc; sxc2[qb] += dxc*dxc; scxc[qb] += dxc*dxc/(sxcmssw*sxcmssw);
            hp[22]->Fill(dyc, weight);
            hp[23]->Fill(dxc, weight);
            if(nypix > 1) {hp[36]->Fill(dyc);}
            if(nxpix > 1) {hp[37]->Fill(dxc);}
            pp[48]->Fill((double)cotbeta,(double)nxpix);
            if(nypix == 1) {hp[24]->Fill(dy);} else {hp[30]->Fill(dy);}
            if(nxpix == 1) {hp[25]->Fill(dx);} else {hp[31]->Fill(dx);}
            if(qbin == 1) {pp[1]->Fill(eta, dyc);}
            if(qb < 4 && qb > 1) {pp[3]->Fill(eta, dyc);}
            if(qbin == 1) {pp[5]->Fill(eta, dxc);}
            if(qb < 4 && qb > 1) {pp[7]->Fill(eta, dxc);}
            if(qbin == 1) {pp[45]->Fill(alpha, dxc);}
            if(qb < 4 && qb > 1) {pp[47]->Fill(alpha, dxc);}
            if(nypix == 2 && yfrac >= 0.) {pp[8]->Fill(yfrac, dyc);}
            if(nypix == 3 && yfrac >= 0.) {pp[9]->Fill(yfrac, dyc);}
            if(nypix == 4 && yfrac >= 0.) {pp[10]->Fill(yfrac, dyc);}
            if(nxpix == 2 && xfrac >= 0.) {pp[11]->Fill(xfrac, dxc);}
            if(nxpix == 3 && xfrac >= 0.) {pp[12]->Fill(xfrac, dxc);}
            if(nxpix == 4 && xfrac >= 0.) {pp[13]->Fill(xfrac, dxc);}
            if(nypix == 2 && yfrac >= 0.) {pp[14]->Fill(yfrac, dy);}
            if(nypix == 3 && yfrac >= 0.) {pp[15]->Fill(yfrac, dy);}
            if(nypix == 4 && yfrac >= 0.) {pp[16]->Fill(yfrac, dy);}
            if(nxpix == 2 && xfrac >= 0.) {pp[17]->Fill(xfrac, dx);}
            if(nxpix == 3 && xfrac >= 0.) {pp[18]->Fill(xfrac, dx);}
            if(nxpix == 4 && xfrac >= 0.) {pp[19]->Fill(xfrac, dx);}
            if(qbin==1) {
                if(nypix == 2 && yfrac >= 0.) {pp[20]->Fill(yfrac, dyc);}
                if(nypix == 3 && yfrac >= 0.) {pp[21]->Fill(yfrac, dyc);}
                if(nypix == 4 && yfrac >= 0.) {pp[22]->Fill(yfrac, dyc);}
                if(nxpix == 2 && xfrac >= 0.) {pp[23]->Fill(xfrac, dxc);}
                if(nxpix == 3 && xfrac >= 0.) {pp[24]->Fill(xfrac, dxc);}
                if(nxpix == 4 && xfrac >= 0.) {pp[25]->Fill(xfrac, dxc);}
                if(nypix == 2 && yfrac >= 0.) {pp[26]->Fill(yfrac, dy);}
                if(nypix == 3 && yfrac >= 0.) {pp[27]->Fill(yfrac, dy);}
                if(nypix == 4 && yfrac >= 0.) {pp[28]->Fill(yfrac, dy);}
                if(nxpix == 2 && xfrac >= 0.) {pp[29]->Fill(xfrac, dx);}
                if(nxpix == 3 && xfrac >= 0.) {pp[30]->Fill(xfrac, dx);}
                if(nxpix == 4 && xfrac >= 0.) {pp[31]->Fill(xfrac, dx);}
            }
            if(qbin > 1) {
                if(nypix == 2 && yfrac >= 0.) {pp[32]->Fill(yfrac, dyc);}
                if(nypix == 3 && yfrac >= 0.) {pp[33]->Fill(yfrac, dyc);}
                if(nypix == 4 && yfrac >= 0.) {pp[34]->Fill(yfrac, dyc);}
                if(nxpix == 2 && xfrac >= 0.) {pp[35]->Fill(xfrac, dxc);}
                if(nxpix == 3 && xfrac >= 0.) {pp[36]->Fill(xfrac, dxc);}
                if(nxpix == 4 && xfrac >= 0.) {pp[37]->Fill(xfrac, dxc);}
                if(nypix == 2 && yfrac >= 0.) {pp[38]->Fill(yfrac, dy);}
                if(nypix == 3 && yfrac >= 0.) {pp[39]->Fill(yfrac, dy);}
                if(nypix == 4 && yfrac >= 0.) {pp[40]->Fill(yfrac, dy);}
                if(nxpix == 2 && xfrac >= 0.) {pp[41]->Fill(xfrac, dx);}
                if(nxpix == 3 && xfrac >= 0.) {pp[42]->Fill(xfrac, dx);}
                if(nxpix == 4 && xfrac >= 0.) {pp[43]->Fill(xfrac, dx);}
            }
        }

    }


    /*  Determine current time */

    gettimeofday(&now1, &timz);
    deltas = now1.tv_sec - now0.tv_sec;
    deltaus = now1.tv_usec - now0.tv_usec;
    deltat = ((double)deltaus)/1000000.;
    deltat += (double)deltas;
    printf("ellapsed time = %f seconds \n", deltat);

    printf(" total events = %d, probx > 10^{-3} = %d, proby > 10^{-3} = %d \n", totale, goodx, goody);
    printf(" low q failures = %d, malformed clusters = %d \n", nbin[4], nbad);
    printf(" detector efficiency [>= 1 pixel above seeding threshold = %f \n", float(nclust)/float(nevent));

    double flutot = 0.;
    for(j=0; j<4; ++j) { flutot += fluence[j];}
    printf("\n fluence for mods 0-3 = %lf, %lf, %lf, %lf \n", fluence[0]/flutot, fluence[1]/flutot,fluence[2]/flutot,fluence[3]/flutot);
    float nbtot = (float)(nbin[0]+nbin[1]+nbin[2]+nbin[3]+nbin[4]);
    printf("\n bin fractions 0-4 = %f, %f, %f, %f, %f \n", (float)nbin[0]/nbtot,nbin[1]/nbtot,nbin[2]/nbtot,(float)nbin[3]/nbtot,nbin[4]/nbtot);



    for(j=0; j<5; ++j) {
        sy[j] /= (float)nbin[j]; sy2[j] /= (float)nbin[j]; scy[j] /= (float)nbin[j];
        sy2[j] = sqrt((double)(sy2[j] - sy[j]*sy[j]));
        printf(" avg y residual[%1d] = %f +- %f, avg y chi^2 = %f \n", j, sy[j], sy2[j], scy[j]);       
        sx[j] /= (float)nbin[j]; sx2[j] /= (float)nbin[j]; scx[j] /= (float)nbin[j];
        sx2[j] = sqrt((double)(sx2[j] - sx[j]*sx[j]));
        printf(" avg x residual[%1d] = %f +- %f, avg x chi^2 = %f \n", j, sx[j], sx2[j], scx[j]); 
    }
    printf(" After 10^{-2}(x) 10^{-3}(y) probability cuts: \n");

    for(j=0; j<5; ++j) {
        syp[j] /= (float)nbin[j]; syp2[j] /= (float)nbin[j]; scyp[j] /= (float)nbin[j];
        syp2[j] = sqrt((double)(syp2[j] - syp[j]*syp[j]));
        printf(" avg y residual[%1d] = %f +- %f, avg y chi^2 = %f \n", j, syp[j], syp2[j], scyp[j]);       
        sxp[j] /= (float)nbin[j]; sxp2[j] /= (float)nbin[j]; scxp[j] /= (float)nbin[j];
        sxp2[j] = sqrt((double)(sxp2[j] - sxp[j]*sxp[j]));
        printf(" avg x residual[%1d] = %f +- %f, avg x chi^2 = %f \n", j, sxp[j], sxp2[j], scxp[j]); 
    }

    printf(" CMSSW algorithm \n");

    for(j=0; j<5; ++j) {
        syc[j] /= (float)nbin[j]; syc2[j] /= (float)nbin[j]; scyc[j] /= (float)nbin[j];
        syc2[j] = sqrt((double)(syc2[j] - syc[j]*syc[j]));
        printf(" avg y residual[%1d] = %f +- %f, avg y chi^2 = %f \n", j, syc[j], syc2[j], scyc[j]);       
        sxc[j] /= (float)nbin[j]; sxc2[j] /= (float)nbin[j]; scxc[j] /= (float)nbin[j];
        sxc2[j] = sqrt((double)(sxc2[j] - sxc[j]*sxc[j]));
        printf(" avg x residual[%1d] = %f +- %f, avg x chi^2 = %f \n", j, sxc[j], sxc2[j], scxc[j]); 
    }

    // Make resolution graphs from profile errors

    float binx, biny, bino;
    printf("\n Fraction malformed clusters = %lf\n", bade/tote);
    for(j=0; j<11; ++j) {
        binx=0.125+j*0.25; biny=nb[j]/nt[j];
        printf("%f  %f \n",binx, biny);
    }
    for(i=0; i<8; ++i) {
        printf("\n Profile errors %d, %s \n",i, pp[i]->GetTitle());
        for(j=0; j<12; ++j) {
            binx=pp[i]->GetBinCenter(j+1); biny=pp[i]->GetBinError(j+1);
            printf("%f  %f \n",binx, biny);
        }
    }

    for(i=44; i<48; ++i) {
        printf("\n Profile errors %d, %s \n",i, pp[i]->GetTitle());
        for(j=0; j<25; ++j) {
            binx=pp[i]->GetBinCenter(j+1); biny=pp[i]->GetBinError(j+1), bino=pp[i]->GetBinContent(j+1);
            printf("%f  %f %f \n",binx, biny, bino);
        }
    }

    /*
     * Histograms plotting
     */
    for(i=0; i<5; ++i) {hp[i]->Fit("gaus"); hp[i+10]->Fit("gaus");}
    for(i=22; i<26; ++i) {hp[i]->Fit("gaus");}
    for(i=30; i<32; ++i) {hp[i]->Fit("gaus");}
    for(i=34; i<41; ++i) {hp[i]->Fit("gaus");}

    //  Create an output filename for this run 

    sprintf(outfile0,"pixel_histos%5.5d.pdf[",nfile);
    sprintf(outfile1,"pixel_histos%5.5d.pdf",nfile);
    sprintf(outfile2,"pixel_histos%5.5d.pdf]",nfile);
    TCanvas c1("c1", header);
    c1.SetFillStyle(4000);
    c1.Print(outfile0);
    for(i=0; i<41; ++i) {
        hp[i]->Draw();
        c1.Print(outfile1);
    }
    for(i=0; i<8; ++i) {
        pp[i]->Draw();
        c1.Print(outfile1);
    }
    for(i=8; i<14; ++i) {
        pp[i]->Draw();
        pp[i+6]->Draw("same");
        c1.Print(outfile1);
    }
    for(i=20; i<26; ++i) {
        pp[i]->Draw();
        pp[i+6]->Draw("same");
        c1.Print(outfile1);
    }
    for(i=32; i<38; ++i) {
        pp[i]->Draw();
        pp[i+6]->Draw("same");
        c1.Print(outfile1);
    }
    pp[48]->Draw();
    c1.Print(outfile1);

    for(i=49; i<53; ++i) {
        pp[i]->Draw();
        c1.Print(outfile1);
    }
    c1.Print(outfile2);

    return 0;
} // MAIN__ 



