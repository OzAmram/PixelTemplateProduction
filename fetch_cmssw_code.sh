#!/bin/sh
#Get the need code from cmssw that is used for template generation.
#You can set the branch variable if you want to use a different branch than cmssw master

base=https://raw.githubusercontent.com/
#branch=cms-sw/cmssw/master/
branch=pmaksim1/cmssw/sideload-for-SiPixelTemplate/
mkdir cmssw_code
cd cmssw_code/
rm *
wget ${base}${branch}RecoLocalTracker/SiPixelRecHits/src/SiPixelTemplateReco.cc 
wget ${base}${branch}RecoLocalTracker/SiPixelRecHits/src/SiPixelTemplateReco2D.cc 
wget ${base}${branch}RecoLocalTracker/SiPixelRecHits/src/VVIObj.cc                        
wget ${base}${branch}RecoLocalTracker/SiPixelRecHits/src/VVIObjF.cc                        
wget ${base}${branch}RecoLocalTracker/SiPixelRecHits/src/sicif.h                       

wget ${base}${branch}RecoLocalTracker/SiPixelRecHits/interface/SiPixelTemplateReco.h 
wget ${base}${branch}RecoLocalTracker/SiPixelRecHits/interface/SiPixelTemplateReco2D.h 
wget ${base}${branch}RecoLocalTracker/SiPixelRecHits/interface/VVIObj.h                          
wget ${base}${branch}RecoLocalTracker/SiPixelRecHits/interface/VVIObjF.h                       

wget ${base}${branch}CondFormats/SiPixelTransient/src/SiPixelGenError.cc
wget ${base}${branch}CondFormats/SiPixelTransient/src/SiPixelTemplate.cc
wget ${base}${branch}CondFormats/SiPixelTransient/src/SiPixelTemplate2D.cc

wget ${base}${branch}CondFormats/SiPixelTransient/interface/SiPixelGenError.h
wget ${base}${branch}CondFormats/SiPixelTransient/interface/SiPixelTemplate.h
wget ${base}${branch}CondFormats/SiPixelTransient/interface/SiPixelTemplate2D.h
wget ${base}${branch}CondFormats/SiPixelTransient/interface/SiPixelTemplateDefs.h
