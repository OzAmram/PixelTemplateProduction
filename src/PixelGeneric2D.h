//
//   PixelGenericReco2D.h (v2.00)
//
//
//    Created by Morris Swartz on 04/03/07.
//  Copyright 2007 __TheJohnsHopkinsUniversity__. All rights reserved.
//
//
 
#ifndef SI_PIXEL_TEMPLATE_STANDALONE
#include "CondFormats/SiPixelObjects/interface/SiPixelGenError.h"
#else
#include "../cmssw_code/SiPixelGenError.h"
#endif

#include <vector>

#ifndef SiPixelTemplateClusMatrix
#define SiPixelTemplateClusMatrix 1

namespace SiPixelTemplateReco {

     struct ClusMatrix {
      float & operator()(int x, int y) { return matrix[mcol*x+y];} 
      float operator()(int x, int y) const { return matrix[mcol*x+y];} 
      float * matrix;
      bool const * xdouble;
      bool const * ydouble;
      int mrow, mcol;
     };
#endif


namespace SiPixelTemplateReco
{
   
   typedef boost::multi_array<float, 2> array_2d;
   
   int PixelGeneric2D(int ID, float cotalpha, float cotbeta, float locBz, float locBx, ClusMatrix & cluster,
                      SiPixelGenError& templ,
                      float& yrec, float& sigmay, float& xrec, float& sigmax, int& nypix, int& nxpix, float& yfrac, float& xfrac);
   int PixelGeneric2D(int ID, float cotalpha, float cotbeta, ClusMatrix & cluster,
                      SiPixelGenError& templ,
                      float& yrec, float& sigmay, float& xrec, float& sigmax, int& nypix, int& nxpix, float& yfrac, float& xfrac);
}


