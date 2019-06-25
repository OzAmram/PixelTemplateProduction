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
#include "TProfile.h"
#include "TStyle.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "Math/DistFunc.h"



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

void get_label(FILE *ifp, char *label, unsigned int size){
    memset(label, ' ', sizeof(char) * size);
    bool got_label = false;
    int i =0;
    while(!got_label){
        fgets(label, size, ifp);
        if(label[0] != '\n') got_label = true;
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
    }

//  If all random numbers used up, generate 210 more 

    if (ibase == TYTEN) {
	   for (i__ = 0; i__ < TYTEN-1; i__ += 2) {
	      r1 = ((float)random())/((float)RAND_MAX);
	      r2 = ((float)random())/((float)RAND_MAX);
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
