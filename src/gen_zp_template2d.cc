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

#define TEMPL_DEBUG
#include "template_utils.h"



// Main program  

int main(int argc, char *argv[])
{

    // Local variables 
    std::vector<float> pvec(6), wgauss(TYSIZE), vgauss(TYSIZE), xgauss(TYSIZE), ygauss(TYSIZE), zgauss(TYSIZE);
    float pixin[TXSIZE][TYSIZE];
    bool ydouble[TYSIZE], xdouble[TXSIZE];
    float ztemp[9][TYSIZE], ptemp[9][TXSIZE], xpar[4][5];
    float dummy[T2YSIZE];
    float sxmax, symax, sxmaxx, symaxx, cosx, cosy, cosz;
    static float thick, xsize, ysize, noise, zcen, gain_frac, q100_frac, common_frac, readout_noise, qscale;
    static float qavg, qxyavg, clslnx, clslny, fbin[3] = {1.5f, 1.0f, 0.85f};
    static float xhit, yhit, xrec, yrec, sigmax, sigmay, signal, cotalpha, cotbeta, qclust, locBz, locBx, probxy, probQ, pixmax;
    static float pixmaxy, pixmaxx;
    static int startfile, neh, nevent, tempID, nbad, fe_model_type, icol, ndcol, numrun; 
    int  id,NTy, NTyx,NTxx,IDtype;
    static float Bfield,Vbias,temp,fluenc;
    static vector<int> nbin(5,0);
    float deltay;
    int ierr, qbin, qb, jmin, jmax, imin, imax, idcol, edgeflagx, edgeflagy, npixels;
    int mrow = TXSIZE, mcol = TYSIZE;
    const int TXSHIFT = (TXSIZE - T2XSIZE)/2;
    double dx, dy;  
    float scaley, scalex, scalx[4], scaly[4], delyavg, delysig, offsetx[4], offsety[4];
    static float q100, q101, q50, q51,  qmax; 
    const float fmax = 0.5f;
    int write_temp_header, use_l1_offset, do_cluster_healing;

    const int nvers = 21;

    float qin;
    static char infile[120], label[160], header[120], dtitle[80], outfile0[120], outfile1[120], outfile2[120], histStore_outfile[120];
    unsigned int detType(1);
    static std::vector<double> cotalphaEdges;
    static std::vector<double> cotbetaEdges;
    int triplg(std::vector<float>&);
    //	int random(void);

    float cluster[TXSIZE][TYSIZE], clust[TXSIZE][TYSIZE], rclust[TXSIZE][TYSIZE];
    std::pair<int, int> pixel, max;

    FILE *output_file;

    struct timeval now0, now1;
    struct timezone timz;
    long deltas, deltaus;
    double deltat;
    float xtalk_frac, xtalk_noise;




    TCanvas* c1 = new TCanvas("c1", header, 800, 800);
    c1->SetFillStyle(4000);


    //  Read which data and inputs to use (use c file i/o which is much faster than c++ i/o) 

    FILE *config_file = fopen("pix_2t.proc", "r");
    if (config_file==NULL) {
        printf("no pixel initialization file found \n");
        return 0;
    }


    char extra[80];
    char line[160];
    fgets(line, 160, config_file);

    int num_read = sscanf(line,"%d %d %f %f %f %f %f %f %f %d %s", &startfile, &numrun, &noise, &q100, 
            &q101, &q100_frac, &common_frac, &gain_frac, &readout_noise, &fe_model_type, &extra[0]);
    printf("processing %d files starting from %d, noise = %f, threshold0 = %f, threshold1 = %f," 
            "rms threshold frac = %f, common_frac = %f, gain fraction = %f, readout noise = %f, front end model type = %d, extra = %s \n", 
            numrun, startfile, noise, q100, q101, q100_frac, common_frac, gain_frac, readout_noise, fe_model_type, extra);
    if(num_read < 10){
        printf("Error reading config file !. Only read %i params \n", num_read);
        return 0;
    }


    fgets(line, 160, config_file);
    num_read = sscanf(line, " %d %d %f %f %d", &use_l1_offset, &write_temp_header, &xtalk_frac, &xtalk_noise, &do_cluster_healing);
    if(num_read != 5){
        printf("Error reading config file !\n");
        printf("Line was %s \n", line);
        return 0;
    }
    fgets(line, 160, config_file);
    num_read = sscanf(line, " %d %d %d %d %d %f %f %f %f %f",  &id, &NTy, &NTyx,&NTxx, &IDtype, &Bfield, &Vbias, &temp, &fluenc, &qscale);
    NTy =0 ;
    printf("Using params: Use_l1_offset=%d, write_temp_header=%d, ID=%d NTy=%d NTyx=%d NTxx=%d Dtype=%d Bfield=%.2f "
            "Bias Voltage = %.1f temparature = %.0f fluence = %.2f q-scale = %.4f xtalk_frac = %.2f xtalk_noise %.2f \n",
            use_l1_offset, write_temp_header, id, NTy, NTyx, NTxx, IDtype, Bfield, Vbias, temp, fluenc, qscale, xtalk_frac, xtalk_noise );
    
    bool do_reso_hists = false;
    auto is_space = [](unsigned char const c) { return std::isspace(c); }; 
    if (fgets(line, 160, config_file) != NULL){
        std::string s_line = string(line);
        if(!std::all_of(s_line.begin(), s_line.end(), is_space)){ //read binning for reso histograms if line not blank
            do_reso_hists = true;


            num_read = sscanf(line, "%s", dtitle);
            if(num_read != 1){
                printf("Error reading config file !\n");
                printf("Line was %s \n", line);
                return 0;
            }
            
            fgets(line, 160, config_file);
            std::string s_binedge;
            istringstream linestream (line);
            while (getline(linestream, s_binedge, ' ')) {
                if (s_binedge == " ") {continue;}
                else {
                    cotbetaEdges.push_back(std::stod(s_binedge,0));
                }
            }

            fgets(line, 160, config_file);
            linestream.clear(); linestream.str(line);
            while (getline(linestream, s_binedge, ' ')) {
                if (s_binedge == " ") {continue;}
                else {
                    cotalphaEdges.push_back(std::stod(s_binedge,0));
                }
            }
            printf("Using cotbeta binning:\n\t");
            for (float e: cotbetaEdges) {
                printf("%f, ",e);
            }
            printf("\nUsing cotalpha binning:\n\t");
            for (float e: cotalphaEdges) {
                printf("%f, ",e);
            }
        }
    }


    fclose(config_file);
    //  Calculate 50% of threshold in q units and enc noise in adc units

    q50=0.5*q100;
    q51=0.5*q101;

    //  Open template output file

    sprintf(infile,"template_summary2D_zp%4.4d.out",id);
    output_file = fopen(infile, "w");
    if (output_file==NULL) {
        printf("couldn't open template output file/n");
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

    FrontEndModel frontEnd;
    frontEnd.fe_type       = fe_model_type;
    frontEnd.gain_frac     = gain_frac;
    frontEnd.readout_noise = readout_noise;
    frontEnd.threshold = q100;
    if(use_l1_offset) {
        printf("using L1 parameters \n");
        frontEnd.vcal = 50.;	
        frontEnd.vcaloffst = 670.;
    }

    sprintf(histStore_outfile,"pixel_histos%5.5d.root",id);
    PixelResolutionHistograms * fastSimResHistoStore;
    if(do_reso_hists){
        fastSimResHistoStore =  new PixelResolutionHistograms( histStore_outfile,                                // File name for histograms
			"",                                     // No subdirectory
			dtitle,                                 // Descriptive title	     
			detType, // unsigned int detType,             // Do we need this?
            cotbetaEdges,
            cotalphaEdges
            );
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


    // Create template object

    std::vector< SiPixelTemplateStore2D > thePixelTemp_;
    SiPixelTemplate2D templ(thePixelTemp_);

    //  Set the ID to -1 to flag the special reco mode

    tempID = -1;

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

        //why flip?
        for(int i = 1; i > -1; --i) {
            fscanf(ztemp_file,"%f %f %f %f %f", &xpar[i][0], &xpar[i][1], &xpar[i][2], &xpar[i][3], &xpar[i][4]);
        }

        for (int k=0; k < 9; ++k) {

            // Skip labels   
            get_label(ztemp_file, label, 160);
            printf("%d %s\n", k, label);
            for(int i=0; i<TYSIZE; i++){
                fscanf(ztemp_file, " %f ", &ztemp[k][i]);
            }
        }

        fclose(ztemp_file);


        // Calculate the mean cluster size in pixels
        clslny = get_clust_len(ztemp, TYSIZE, symaxx);


        //  Read in 1D p template information first

        sprintf(infile,"./ptemp_%5.5d.txt",ifile);

        //  Open input file and read header info 

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

        for(int i = 3; i > 1; --i) {
            fscanf(ptemp_file,"%f %f %f %f %f", &xpar[i][0], &xpar[i][1], &xpar[i][2], &xpar[i][3], &xpar[i][4]);
        }

        for (int k=0; k < 9; ++k) {
            // Skip labels   
            get_label(ptemp_file, label, 160);
            printf("%s\n", label);
            for(int i=0; i<TXSIZE; i++){
                fscanf(ptemp_file, " %f ", &ptemp[k][i]);
            }
        }
        fclose(ptemp_file);

        // Calculate the mean cluster size in pixels
        clslnx = get_clust_len(ptemp, TXSIZE, sxmaxx);


        //  Read in 2D template information

        sprintf(infile,"./zptemp_%5.5d.txt",ifile);

        //  Open input file and read header info 

        FILE *ztemp2d_file = fopen(infile, "r");
        if (ztemp2d_file==NULL) {
            printf("no zp-template file %s \n", infile);
            return 0;
        }


        SiPixelTemplateEntry2D * slice = new SiPixelTemplateEntry2D;

        // Copy info into the slice and reformat from pixelav coordinates to CMSSW local coordinates

        slice->runnum = ifile;

        for(int i = 0; i < 2; ++i) {
            for(int j=0; j<5; ++j) {slice->xypar[i][j] = xpar[i][j];}
        }

        slice->clslenx = clslnx;
        slice->clsleny = clslny;

        fscanf(ztemp2d_file,"%f  %f  %f", &cosx, &cosy, &cosz);
        //	   printf("cosx/cosy/cosz = %f/%f/%f \n", cosx, cosy, cosz);

        slice->costrk[0] = -cosy;
        slice->costrk[1] = -cosx;
        slice->costrk[2] = -cosz;
        slice->cotalpha = cosy/cosz;
        slice->cotbeta = cosx/cosz;

        fscanf(ztemp2d_file,"%f %f %d %d %d %d", &qxyavg, &pixmax, &imin, &imax, &jmin, &jmax);
        printf("qxyavg/pixmax/imin/imax/jmin/jmax = %f/%f/%d/%d/%d/%d \n", qxyavg, pixmax, imin, imax, jmin, jmax);

        slice->qavg = qxyavg;
        slice->pixmax = pixmax;
        slice->sxymax = symax;

        slice->iymin = TYSIZE - imax;
        slice->iymax = TYSIZE - imin;
        slice->jxmin = TXSIZE - jmax - TXSHIFT;
        slice->jxmax = TXSIZE - jmin - TXSHIFT;

        for(int i = 1; i > -1; --i) {
            fscanf(ztemp2d_file,"%f %f %f %f %f", &slice->lanpar[i][0], &slice->lanpar[i][1], &slice->lanpar[i][2], 
                    &slice->lanpar[i][3], &slice->lanpar[i][4]);
        }

        for(int l=6; l>-1; --l) {
            for(int k=6; k>-1; --k) {
                // Skip labels   
                get_label(ztemp2d_file, label, 160);
                printf("%d %s\n", k, label);
                for(int j=0; j<TXSHIFT; ++j) {
                    fscanf(ztemp2d_file,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", 
                            &dummy[20],&dummy[19],&dummy[18],&dummy[17],&dummy[16],&dummy[15],
                            &dummy[14],&dummy[13],&dummy[12],&dummy[11],&dummy[10],&dummy[9],
                            &dummy[8],&dummy[7],&dummy[6],&dummy[5],&dummy[4],&dummy[3],&dummy[2],&dummy[1],&dummy[0]);
                }
                for(int j=T2XSIZE-1; j>-1; --j) {
                    fscanf(ztemp2d_file,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", 
                            &dummy[20],&dummy[19],&dummy[18],&dummy[17],&dummy[16],&dummy[15],
                            &dummy[14],&dummy[13],&dummy[12],&dummy[11],&dummy[10],&dummy[9],
                            &dummy[8],&dummy[7],&dummy[6],&dummy[5],&dummy[4],&dummy[3],&dummy[2],&dummy[1],&dummy[0]);
                    for(int i=0; i<21; ++i) {slice->xytemp[k][l][i][j] = (short int)dummy[i];}
                }
                for(int j=0; j<TXSHIFT; ++j) {
                    fscanf(ztemp2d_file,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", 
                            &dummy[20],&dummy[19],&dummy[18],&dummy[17],&dummy[16],&dummy[15],
                            &dummy[14],&dummy[13],&dummy[12],&dummy[11],&dummy[10],&dummy[9],
                            &dummy[8],&dummy[7],&dummy[6],&dummy[5],&dummy[4],&dummy[3],&dummy[2],&dummy[1],&dummy[0]);
                }
            }
        }

        slice->scalexavg = 1.f;
        slice->scaleyavg = 1.f;

        fclose(ztemp2d_file);

        //  Create an input filename for this run 


        sprintf(infile,"template_events_d%05i.out",ifile);




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

        nevent=0;
        nbad = 0;

        if(write_temp_header && ifile==startfile) {
            fprintf(output_file,"%s", header);
            printf("Using params: ID=%d NTyx=%d NTxx=%d Dtype=%d Bfield=%.2f Bias Voltage = %.1f temparature = %.0f fluence = %.2f q-scale = %.4f \n",
                    id, NTyx, NTxx, IDtype, Bfield, Vbias, temp, fluenc, qscale);


            fprintf(output_file,"%d %d %4.2f %d %d %d %d %5.1f %5.1f %4.2f %5.3f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %f %4.2f %4.2f %4.2f \n",
                    id,nvers,Bfield,NTy,NTyx,NTxx,IDtype,Vbias, temp,fluenc,qscale,q50,lorwdy,
                    lorwdx,ysize,xsize,thick,q51,lorbsy,lorbsx,1.5f,1.0f,0.85f);
        }
        else {
            for(int i=0; i<26; ++i) { hp[i]->Reset();}
        }  


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
            for(int i=0; i<ndcol; ++i) {ndhit[i] = 0;}
            icol = 0;
            if(vgauss[1] < 0.) {icol = 1;}

            int xtalk_row_start = 0;
            int xtalk_unfold_row = 1;
            if(vgauss[2] < 0.) {
                xtalk_row_start = 1;
                xtalk_unfold_row = 0;
            }

            float xtalk_apply = xtalk_frac + xtalk_noise * vgauss[3];

            if (xtalk_frac > 0.) apply_xtalk(pixin, xtalk_row_start, xtalk_apply);



            for(int j=0; j<TXSIZE; ++j) {
                triplg(wgauss);
                triplg(xgauss);
                triplg(ygauss);
                triplg(zgauss);
                for(int i=0; i<TYSIZE; ++i) {
                    qin = (10.*pixin[j][i] + xgauss[i]*noise);
                    rclust[TXSIZE-1-j][TYSIZE-1-i] = qin;
                    if(qin < q100*(1.+wgauss[i]*q100_frac)) {
                        clust[TXSIZE-1-j][TYSIZE-1-i] = 0.;
                    } else {
                        idcol = (TYSIZE-1-i+icol)/2;
                        ++ndhit[idcol];
                        signal = frontEnd.apply_model( qin, ygauss[i], zgauss[i] );
                        clust[TXSIZE-1-j][TYSIZE-1-i] = (1.+vgauss[0]*common_frac)*signal;
                    }
                }
            }
            if(xtalk_frac > 0.) unfold_xtalk(clust, xtalk_unfold_row, xtalk_frac);

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

            if(qmax < clustering_thresh) continue;


            // Simulate clustering around maximum signal (seed)

            auto pixlst = clusterizer(clust, q100, max, do_cluster_healing);

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

            int xwidth(0), ywidth(0);
            for(int i=0; i<TXSIZE; i++){
                for(int j=0; j<TYSIZE; j++){
                    if(cluster[i][j] >0.){
                        xwidth++;
                        break;
                    }
                }
            }

            for(int j=0; j<TYSIZE; j++){
                for(int i=0; i<TXSIZE; i++){
                    if(cluster[i][j] >0.){
                        ywidth++;
                        break;
                    }
                }
            }

            //print_cluster(cluster);
            //printf("xwidth %i ywidth %i \n", xwidth, ywidth);




            //bool is_split = check_is_split(cluster);


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

            // No dead columns or module edges

            edgeflagy = 0; //code for no gap OR split cluster
            edgeflagx = 0;



            //        if(fabs(cotbeta) < 2.1) continue;
            // Do the template analysis on the cluster 
            SiPixelTemplateReco2D::ClusMatrix clusterPayload{&cluster[0][0], xdouble, ydouble, mrow,mcol};
            locBx = 1.;
            if(cotbeta < 0.) locBx = -1.;
            locBz = locBx;
            if(cotalpha < 0.) locBz = -locBx;

            //  Sideload this template slice

            templ.sideload(slice, IDtype, locBx, locBz, lorwdy, lorwdx, q50, fbin, xsize, ysize, thick);
            ierr = PixelTempReco2D(tempID, cotalpha, cotbeta, locBz, locBx, edgeflagy, edgeflagx, clusterPayload, templ, yrec, sigmay, xrec, sigmax, probxy, probQ, qbin, deltay, npixels);
            if(ierr != 0) {
                ++nbad; 
                //	      printf("reconstruction failed with error %d \n", ierr);
            } else {
                qb = qbin;
                ++nbin[qb];
                dy = yrec - (TYSIZE/2)*ysize - yhit;
                hp[0]->Fill(dy);
                if(sigmay > 0.f) hp[5]->Fill(dy/sigmay);
                hp[20]->Fill((double)deltay);
                hp[1+qbin]->Fill(dy, 1.);
                if(sigmay > 0.f) hp[6+qbin]->Fill(dy/sigmay);
                dx = xrec - (TXSIZE/2)*xsize - xhit;
                hp[10]->Fill(dx);
                if(sigmax > 0.f) hp[15]->Fill(dx/sigmax);
                hp[11+qbin]->Fill(dx);
                if(sigmax > 0.f) hp[16+qbin]->Fill(dx/sigmax);
                hp[22]->Fill((double)probxy);
                hp[24]->Fill((double)npixels);
                hp[23]->Fill((double)(probxy/npixels));

                // Fill the FastSim histograms
                if(do_reso_hists) fastSimResHistoStore->Fill( dx, dy, (double)cotalpha, (double)cotbeta, qbin, xwidth, ywidth );

            }

        }

        fclose(events_file);
        printf(" low q failures = %d, malformed clusters = %d \n", nbin[4], nbad);	   

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

        TF1 *cfunc0 = new TF1("chisquare0",chisquare0,0.,100.,3);
        cfunc0->SetParLimits(1,0.01,50);
        cfunc0->SetParNames("norm","mean","scale");
        Double_t cp01 = hp[22]->GetMean()/2.;
        Double_t cp02 = 0.75;
        Double_t cp00 = hp[22]->GetEntries()*2./cp01; 
        cfunc0->SetParameters(cp00,cp01,cp02);
        hp[22]->Fit("chisquare0");
        Double_t cpar0[3];
        cfunc0->GetParameters(cpar0);
        printf("cpar0 = %lf/%lf/%lf \n", cpar0[0], cpar0[1], cpar0[2]);


        TF1 *cfunc1 = new TF1("chisquare1",chisquare1,0.,50.,3);
        cfunc1->SetParNames("norm","mean","scale");
        double pixels = hp[24]->GetMean();
        Double_t cp11 = cpar0[1]/pixels;
        Double_t cp12 = cpar0[2];
        Double_t cp10 = hp[23]->GetEntries()*2./cp11; 
        cfunc1->SetParameters(cp10,cp11,cp12);
        cfunc1->SetParLimits(1,cp11,cp11);
        cfunc1->SetParLimits(2,cp12,cp12);
        hp[23]->Fit("chisquare1");

        //create a function with 4 parameters in the range [-10,50]

        TF1 *vfunc = new TF1("vavilov",vavilov,-10.,50.,4);
        vfunc->SetParNames("norm","mean","sigma","kappa");
        Double_t p1 = hp[25]->GetMean();
        Double_t p2 = 0.1*p1;
        Double_t p3 = 0.02*p1/20000.;
        Double_t p0 = hp[25]->GetEntries()*20000./p1; 
        vfunc->SetParameters(p0,p1,p2,p3);
        hp[25]->Fit("vavilov");
        Double_t par[4];
        vfunc->GetParameters(par);

        //  Create an output filename for this run 

        sprintf(outfile0,"template2d_histos%5.5d.pdf[",ifile);
        sprintf(outfile1,"template2d_histos%5.5d.pdf",ifile);
        sprintf(outfile2,"template2d_histos%5.5d.pdf]",ifile);
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

        fprintf(output_file,"%d %8.6f %8.6f %8.6f \n", ifile, -cosy, -cosx, -cosz);
        rnelec /=float(nevent);
        fprintf(output_file,"%7.1f %5.1f %5.1f %d %d %d %d \n", rnelec, pixmax, symax, slice->iymin, slice->iymax, slice->jxmin, slice->jxmax);
        for(int i = 0; i < 2; ++i) {
            fprintf(output_file,"%e %e %e %e %e \n", slice->xypar[i][0], slice->xypar[i][1], slice->xypar[i][2], 
                    slice->xypar[i][3], slice->xypar[i][4]);
        }
        for(int i = 0; i < 2; ++i) {
            fprintf(output_file,"%e %e %e %e %e \n", slice->lanpar[i][0], slice->lanpar[i][1], slice->lanpar[i][2], 
                    slice->lanpar[i][3], slice->lanpar[i][4]);
        }
        for(int l=0; l<7; ++l) {
            for(int k=0; k<7; ++k) {
                for(int j=0; j<T2XSIZE; ++j) {
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
    }
    // Close output file   

    fclose(output_file);  
    if(do_reso_hists) delete fastSimResHistoStore;

    /*  Determine current time */

    gettimeofday(&now1, &timz);
    deltas = now1.tv_sec - now0.tv_sec;
    deltaus = now1.tv_usec - now0.tv_usec;
    deltat = ((double)deltaus)/1000000.;
    deltat += (double)deltas;
    printf("ellapsed time = %f seconds \n", deltat);

    return 0;
} // MAIN__ 





