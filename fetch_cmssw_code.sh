#!/bin/sh

mkdir cmssw_code
cd cmssw_code/
wget https://raw.githubusercontent.com/cms-sw/cmssw/master/RecoLocalTracker/SiPixelRecHits/src/SiPixelTemplateReco.cc 
wget https://raw.githubusercontent.com/cms-sw/cmssw/master/RecoLocalTracker/SiPixelRecHits/src/SiPixelTemplateReco2D.cc 
wget https://raw.githubusercontent.com/cms-sw/cmssw/master/RecoLocalTracker/SiPixelRecHits/src/VVIObj.cc                        
wget https://raw.githubusercontent.com/cms-sw/cmssw/master/RecoLocalTracker/SiPixelRecHits/src/VVIObjF.cc                        
wget https://raw.githubusercontent.com/cms-sw/cmssw/master/RecoLocalTracker/SiPixelRecHits/src/sicif.h                       

wget https://raw.githubusercontent.com/cms-sw/cmssw/master/RecoLocalTracker/SiPixelRecHits/interface/SiPixelTemplateReco.h 
wget https://raw.githubusercontent.com/cms-sw/cmssw/master/RecoLocalTracker/SiPixelRecHits/interface/SiPixelTemplateReco2D.h 
wget https://raw.githubusercontent.com/cms-sw/cmssw/master/RecoLocalTracker/SiPixelRecHits/interface/VVIObj.h                          
wget https://raw.githubusercontent.com/cms-sw/cmssw/master/RecoLocalTracker/SiPixelRecHits/interface/VVIObjF.h                       

wget https://raw.githubusercontent.com/cms-sw/cmssw/master/CondFormats/SiPixelTransient/src/SiPixelGenError.cc
wget https://raw.githubusercontent.com/cms-sw/cmssw/master/CondFormats/SiPixelTransient/src/SiPixelTemplate.cc
wget https://raw.githubusercontent.com/cms-sw/cmssw/master/CondFormats/SiPixelTransient/src/SiPixelTemplate2D.cc

wget https://raw.githubusercontent.com/cms-sw/cmssw/master/CondFormats/SiPixelTransient/interface/SiPixelGenError.h
wget https://raw.githubusercontent.com/cms-sw/cmssw/master/CondFormats/SiPixelTransient/interface/SiPixelTemplate.h
wget https://raw.githubusercontent.com/cms-sw/cmssw/master/CondFormats/SiPixelTransient/interface/SiPixelTemplate2D.h
wget https://raw.githubusercontent.com/cms-sw/cmssw/master/CondFormats/SiPixelTransient/interface/SiPixelTemplateDefs.h
