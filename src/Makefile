root_flags= `root-config --cflags --ldflags --libs --glibs`
compiler_flags=-std=c++11 -Wall -g
vdt_dir=/cvmfs/cms.cern.ch/slc7_amd64_gcc820/cms/vdt/0.4.0-nmpfii/include
includes= -I. -I../cmssw_code/ -I$(vdt_dir) -I/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/boost/1.67.0-pafccj/include
DEPS= template_utils.h ranlux.c ../cmssw_code/*
EXE= compare_templates gen_zp_template gen_zp_template2d gen_xy_template gen_xy_template2d test_template
ODIR=../bin
BOOST_FLAG=-DSI_PIXEL_TEMPLATE_USE_BOOST

all: $(EXE)


gen_zp_template: gen_zp_template.cc $(DEPS)
	g++ gen_zp_template.cc -o $(ODIR)/gen_zp_template $(compiler_flags)  $(root_flags) $(includes)

gen_zp_template2d: gen_zp_template2d.cc $(DEPS)
	g++ gen_zp_template2d.cc -o $(ODIR)/gen_zp_template2d $(compiler_flags)  $(root_flags) $(includes)

gen_xy_template: gen_xy_template.cc $(DEPS)
	g++ gen_xy_template.cc   -o $(ODIR)/gen_xy_template $(compiler_flags) $(root_flags) $(includes)

gen_xy_template2d: gen_xy_template2d.cc $(DEPS)
	g++ gen_xy_template2d.cc  -o $(ODIR)/gen_xy_template2d $(compiler_flags) $(root_flags) $(includes)

compare_templates: compare_templates.cc $(DEPS)
	g++ compare_templates.cc -o $(ODIR)/compare_templates $(compiler_flags) $(root_flags) $(includes)

test_template: test_template.cc $(DEPS)
	g++ test_template.cc -o $(ODIR)/test_template $(compiler_flags) $(root_flags) $(includes)

load_template_test: load_template_test.cc $(DEPS)
	g++ load_template_test.cc -o $(ODIR)/load_template_test $(compiler_flags) $(root_flags) $(includes)

clean:
	rm -f $(ODIR)/*

.PHONY: all clean

