//template_utils.h
//common functions for template making
//by Oz Amram              May 2019

#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <boost/multi_array.hpp>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <set>
#include <sys/time.h>

#define SI_PIXEL_TEMPLATE_STANDALONE
#include "../cmssw_code/SiPixelTemplate.cc"
#include "../cmssw_code/SiPixelTemplate2D.cc"
#include "../cmssw_code/VVIObjF.cc"
#include "../cmssw_code/SiPixelTemplateReco.cc"
#include "../cmssw_code/SiPixelTemplateReco2D.cc"
#include "../cmssw_code/SiPixelUtils.cc"
#include "../cmssw_code/SiPixelGenError.cc"


using namespace std;

#include <TF1.h>
#include "Math/MinimizerOptions.h"
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TObject.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "Math/DistFunc.h"

#ifdef TEMPL_DEBUG
#include  "ranlux.c"
#endif

float clustering_thresh = 1000.;
float minErrX = 1.25;
float minErrY = 3.;



//-----------------------------------------------------------------------------
//  A simple struct to collect the information about the electronics model.
//  The kind of front-end response that will be modeled is selected by fe_type.
//  
//  Note: this is a struct and all variables are public.
//-----------------------------------------------------------------------------
struct FrontEndModel
{
    //--- Variables, read from the config file.
    int   fe_type       = 0;
    float gain_frac     = 0.f;
    float readout_noise = 0.f;

    //--- Variables we can change, but we start with good default values
    double vcal = 47.0;	
    double vcaloffst = 60.0;

    //--- PhaseII - initial guess
    double threshold = 1000; // threshold in e-
    double qperTOT = 1500; // e- per TOT
    int nbitsTOT = 4; // fixed and carved in stone?
    int tot_max = pow(2, nbitsTOT);
    int tot = 0;

    //--- Constants (could be made variables later)
    const double gain  = 3.19;
    const double ped   = 16.46;
    const double p0    = 0.01218;
    const double p1    = 0.711;
    const double p2    = 203.;
    const double p3    = 148.;	

    //-----------------------------------------------------------------------------
    //  A function to simulate the electronics response: different gains, TOT, etc.
    //-----------------------------------------------------------------------------
    float apply_model( float qin, float ygauss, float zgauss )
    {
        float signal = 0.f;

        switch ( fe_type ) {
            case 0:
                signal = qin * (1. + gain_frac*ygauss) + zgauss*readout_noise;
                break;

            case 1:
                {
                    double adc = (double)((int)(p3+p2*tanh(p0*(qin + vcaloffst)/(7.0*vcal) - p1)));
                    signal = ((float)((1.+gain_frac*ygauss)*(vcal*gain*(adc-ped))) - vcaloffst + zgauss*readout_noise);
                }
                break;

            case 2:
                if(qin < threshold) signal = 0;
                else{
                    tot = int((qin - threshold)/qperTOT);
                    tot = std::min(tot, tot_max);
                    double step = qperTOT/2;
                    signal = (qperTOT*tot + threshold)+step;
                }
                break;

            default:
                std::cout << "ElectronicModel::apply_model: illegal fe_type = " << fe_type << std::endl;
                assert(0);
        }

        return signal;
    }

};



// Calculates something like a full-width half max for the cluster
// length
template <typename TwoD> //template to allow arrays of different sizes
float get_clust_len(TwoD &temp, int size, float max){
    float pfrst=0.f, plast=0.f, clsln=0.f;
    for(int i=0; i<size; ++i) {
        for (int k=7; k > -1; --k) {
            if(temp[k][i] > max) {
                float dzsig = temp[k][i] - temp[k+1][i];
                float frac;
                if(dzsig > 0.f) {
                    frac = (temp[k][i] - max)/dzsig;
                    if(frac > 1.f) frac = 1.f;
                    if(frac < 0.f) frac = 0.f;
                } else {
                    frac = 0.f;
                }
                pfrst = i-(k+frac-4.f)/8.f;
                goto first;
            }
        }
    }
first: ;
       for(int i=size; i>-1; --i) {
           for (int k=1; k < 9; ++k) {
               if(temp[k][i] > max) {
                   float dzsig = temp[k][i] - temp[k-1][i];
                   float frac;
                   if(dzsig > 0.f) {
                       frac = (temp[k][i] - max)/dzsig;
                       if(frac > 1.f) frac = 1.f;
                       if(frac < 0.f) frac = 0.f;
                   } else {
                       frac = 0.f;
                   }
                   plast = i-(k-frac-4.f)/8.f;
                   goto second;
               }
           }
       }
second: clsln = plast-pfrst;
        if(clsln < 0.f) clsln = 0.f;
        return clsln;
}




//vavilov distribution (to fit to)
Double_t vavilov(Double_t *v, Double_t *par)
{
    Float_t arg = 0.;
    if (par[2] != 0) arg = (float)(v[0] - par[1])/fabs(par[2]);

    float beta2 = 1.f;   
    float kappa = (float)par[3];
    VVIObjF vvidist(kappa, beta2, 0);
    float xl, xu;
    vvidist.limits(xl,xu);
    if(arg < xl) arg = xl;
    if(arg > xu) arg = xu;
    Double_t fitval = par[0]*(double)vvidist.fcn(arg);
    return fitval;
}


Double_t chisquare0(Double_t *v, Double_t *par)
{
    Double_t fitval = par[0]*ROOT::Math::chisquared_pdf(par[2]*v[0], par[1]);
    return fitval;
}

Double_t chisquare1(Double_t *v, Double_t *par)
{
    Double_t fitval = par[0]*ROOT::Math::chisquared_pdf(par[2]*v[0], par[1]);
    return fitval;
}

int is_empty_str(char *s) {
    while (*s != '\0') {
        if (!isspace((unsigned char)*s))
            return 0;
        s++;
    }
    return 1;
}

void get_label(FILE *ifp, char *label, unsigned int size){
    memset(label, ' ', sizeof(char) * size);
    bool got_label = false;
    int i =0;
    while(!got_label){
        fgets(label, size, ifp);
        if(!is_empty_str(label)) got_label = true;
        i++;
    }
    return;
}

void read_cluster(FILE *ifp, float pixin[TXSIZE][TYSIZE]){
    for (int i=0; i < TXSIZE; i++) {
        for(int  j=0; j < TYSIZE; j++){

            fscanf(ifp, " %f ", &pixin[i][j]);
        }
    }
    return;
}

void apply_xtalk(float pixin[TXSIZE][TYSIZE], int irow_start, float xtfrac = 0.1){
    //matrix looks like | 1-x  x | 
    //                  | x  1-x |
   if(xtfrac == 0.) return;
   float m11 = 1. - xtfrac;
   float m22 = m11;
   float m12 = xtfrac;
   float m21 = m12;


   for(int i=irow_start; i<TXSIZE-1; i += 2) {
       for(int j=0; j<TYSIZE; j++) {
           float old_ij = pixin[i][j];
           float old_i1j = pixin[i+1][j];
           pixin[i][j] = old_ij*m11 +  old_i1j*m12;
           pixin[i+1][j] = old_ij*m21 +  old_i1j*m22;
       }
   }


}


void unfold_xtalk(float clust[TXSIZE][TYSIZE], int irow_start, float xtfrac=0.1){
    //unfold cross talk. Remember cluster is flipped compared to before, 
    // inverse matrix is (1/(2-x)) * |1 - x  -x |
    //                               | -x  1 -x |
   if(xtfrac == 0.) return;
   float m11 = 1.-xtfrac;
   float minv11 = m11/(1. - 2.*xtfrac);
   float minv22 = minv11;
   float minv12 = -xtfrac/(1. - 2.*xtfrac);
   float minv21 = minv12;       

   for(int i=irow_start; i<TXSIZE-1; i += 2) {
       for(int j=0; j<TYSIZE; j++) {
           float old_ij = clust[i][j];
           float old_i1j = clust[i+1][j];
           clust[i][j] = old_ij*minv11 + old_i1j*minv12;
           clust[i+1][j] = old_ij*minv21 + old_i1j*minv22;
       }
   }      		 
}



template <typename TwoD> //template to allow arrays of different sizes
void print_cluster(TwoD &clust){
    for (int i=0; i < TXSIZE; i++) {
        for(int  j=0; j < TYSIZE; j++){

            printf("%.1f ", clust[i][j]);
        }
        printf("\n");
    }
}

//see https://stackoverflow.com/questions/936687/how-do-i-declare-a-2d-array-in-c-using-new
float** setup_2d_array(int size1, int size2){
    float **a = new float*[size1];
    for(int i=0; i< size1; i++){
        a[i] = new float[size2];

        for(int j=0; j<size2; j++){
            a[i][j] = 0.;
        }

    }
    return a;
}

//return 4 parameters for template output
//first two are always mean and std dev
//second two are fitted mean and sigma of gaussian if fit went well
std::vector<float> get_gaussian_pars(TH1F *h, double min_std = 3.){
    std::vector<float> pars;
    double h_mean = h->GetMean();
    double h_std = std::max(h->GetStdDev(), min_std);
    pars.push_back(h_mean);
    pars.push_back(h_std);
    if(h->Integral() > 50.){
        h->Fit("gaus");
        TF1 *fit = h->GetFunction("gaus");
        double mean = fit->GetParameter(1);
        double sigma = std::max(fit->GetParameter(2), min_std);
        if(fabs(mean) < 300. && abs(sigma) < 175. ){
            pars.push_back(mean);
            pars.push_back(sigma);
        }
        else{
            pars.push_back(h_mean);
            pars.push_back(h_std);
        }
    }
    else{
        pars.push_back(h_mean);
        pars.push_back(h_std);
    }
    return pars;
}
std::vector<float> fit_pol5(TProfile *h){
    std::vector<float> pars;
    Double_t stats[6];
    h->GetStats(stats);
    Double_t nentries = stats[0];
    //printf("integral %.2f \n", nentries);
    if(nentries > 10.){
        h->Fit("pol5");
        h->SetStats(false);
        TF1 *fit = h->GetFunction("pol5");
        Double_t chi2 = fit->GetChisquare();
        for(int i=0; i<6; i++){
            if(chi2 < 1e6){
                pars.push_back(fit->GetParameter(i));
            }
            else{
                pars.push_back(0.);
            }

        }
    }
    else{
        for(int i=0; i<6; i++){
            pars.push_back(0.);
        }
    }


    return pars;
}

std::vector<float> get_chi2_pars(TH1F *h, double h_min){
    std::vector<float> pars;
    if(h->Integral() < 20.){
        pars.push_back(0.10);
        pars.push_back(0.16);
        return pars;
    }
    double h_mean = h->GetMean();
    if(h_min < 0.) h_min = 0.;
    if(h_min > h_mean) h_min = 0.;
    double dmean = h_mean - h_min;
    if(dmean < 0.1) dmean = h_mean;
    pars.push_back(dmean);
    pars.push_back(h_min);
    return pars;
}

std::vector<float> get_vavilov_pars(TH1F *h){

    TF1 *vfunc = new TF1("vavilov",vavilov,-10.,50.,4);
    vfunc->SetParNames("norm","mean","sigma","kappa");

    float par_guess [4];
    float mean = h->GetMean();
    par_guess[0] = h->GetEntries()*20000./mean; 
    par_guess[1] = mean;
    par_guess[2] = 0.1*mean;
    par_guess[3] = 0.02*mean/20000.;
    vfunc->SetParameters(par_guess[0], par_guess[1], par_guess[2], par_guess[3]);
    vfunc->SetParLimits(1, 0., 1e9);
    vfunc->SetParLimits(2, 0., 1e9);
    vfunc->SetParLimits(3, 0.01, 10.);
    h->Fit("vavilov");
    std::vector<float> pars;
    //we don't care about norm
    for(int i=1; i<4; i++){
        float par = vfunc->GetParameter(i);
        if((i==1 || i == 2) && par < 1e-4) par = par_guess[i];
        pars.push_back(par);
    }
    return pars;

}





void delete_2d_array(float **a, int size1, int size2 ){
    for(int i=0; i< size1; i++){
        delete[] a[i];
    }
}
float*** setup_3d_array(int size1, int size2, int size3){
    float ***a = new float**[size1];
    for(int i=0; i< size1; i++){
        a[i] = new float*[size2];
        for(int j=0; j<size2; j++){
            a[i][j] = new float[size3];

        }
    }
    for(int i=0; i< size1; i++){
        for(int j=0; j<size2; j++){
            for(int k=0; k<size3; k++){
                a[i][j][k] = 0.;
            }
        }
    }

    return a;
}

void zero_3d_array(float ***a, int size1, int size2, int size3){
    for(int i=0; i< size1; i++){
        for(int j=0; j<size2; j++){
            for(int k=0; k<size3; k++){
                a[i][j][k] = 0.;
            }
        }
    }
}

void zero_2d_array(float **a, int size1, int size2){
    for(int i=0; i< size1; i++){
        for(int j=0; j<size2; j++){
            a[i][j] = 0.;
        }
    }
}


void delete_3d_array(float ***a, int size1, int size2, int size3){
    for(int i=0; i< size1; i++){
        for(int j=0; j<size2; j++){
            delete[] a[i][j];
        }
        delete[] a[i];
    }
    delete[] a;
}



// ***************************************************************** 
//! Calculate 21 gaussianly-distributed random numbers.
//! \param x - a vector holding 21 random numbers generated with rms = 1.            
// ***************************************************************** 
int triplg(std::vector<float>& x)
{
    // Initialized data 

    static int fcall = -1;

#ifdef TEMPL_DEBUG
    static int lux = 3;
    static int seed = 1234567;
    static int k1=0;
    static int k2=0;
    static float rlux_out[TYTEN];
    static int size = TYTEN;
#endif




    // Local variables 
    static float r1, r2;
    static int i__;
    static float r__;
    static int ibase;
    static std::vector<float> rbuff(TYTEN);
    static float twopi;
    double arg, phi;



    // Function Body 

    //  Initalize the parameters 

    if (fcall) {
        twopi = 2.*acos((double)-1.);
        ibase = TYTEN;
        fcall = 0;
#ifdef TEMPL_DEBUG
        rluxgo_(&lux, &seed, &k1, &k2);
#endif
    }

    //  If all random numbers used up, generate 210 more 

    if (ibase == TYTEN) {

#ifdef TEMPL_DEBUG
        ranlux_(rlux_out, &size);
#endif
        for (i__ = 0; i__ < TYTEN-1; i__ += 2) {
#ifdef TEMPL_DEBUG
            r1 = rlux_out[i__];
            r2 = rlux_out[i__+1];
#else
            r1 = ((float)random())/((float)RAND_MAX);
            r2 = ((float)random())/((float)RAND_MAX);
#endif
            arg = (double)(1. - r1);
            if (arg < 1.e-30) {arg = 1.e-30;}
            r__ = sqrt(log(arg) * (-2.));
            phi = twopi * r2;
            rbuff[i__] = r__ * cos(phi);
            rbuff[i__+1] = r__ * sin(phi);
        }
        ibase = 0;
    }
    for (i__ = 0; i__ < TYSIZE; ++i__) {
        x[i__] = rbuff[ibase + i__];
    }
    ibase += TYSIZE;
    return 0;
} // triplg 
