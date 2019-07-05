root_flags= `root-config --cflags --ldflags --libs --glibs`
compiler_flags=-std=c++11 -Wall -O0 -g
vdt_dir=/usr/local/include
includes= -I. -I../cmssw_code/ -I$(vdt_dir)
DEPS= template_utils.h ranlux.c ../cmssw_code/*
EXE= compare_templates gen_zp_template gen_zp_template2d gen_zp_template2d_old gen_xy_template gen_xy_template2d
ODIR=../bin

all: $(EXE)


gen_zp_template: gen_zp_template.cc $(DEPS)
	g++ gen_zp_template.cc -o $(ODIR)/gen_zp_template $(compiler_flags)  $(root_flags) $(includes)

gen_zp_template2d: gen_zp_template2d.cc $(DEPS)
	g++ gen_zp_template2d.cc -o $(ODIR)/gen_zp_template2d $(compiler_flags)  $(root_flags) $(includes)

gen_zp_template2d_old: gen_zp_template2d_2018_2t.cxx
	g++ gen_zp_template2d_2018_2t.cxx -o $(ODIR)/gen_zp_template2d_old $(compiler_flags)  $(root_flags) $(includes)

gen_xy_template: gen_xy_template.cc $(DEPS)
	g++ gen_xy_template.cc   -o $(ODIR)/gen_xy_template $(compiler_flags) $(root_flags) $(includes)

gen_xy_template2d: gen_xy_template2d_2t.cc $(DEPS)
	g++ gen_xy_template2d_2t.cc  -o $(ODIR)/gen_xy_template2d $(compiler_flags) $(root_flags) $(includes)

compare_templates: compare_templates.cc $(DEPS)
	g++ compare_templates.cc -o $(ODIR)/compare_templates $(compiler_flags) $(root_flags) $(includes)

load_template_test: load_template_test.cc $(DEPS)
	g++ load_template_test.cc -o $(ODIR)/load_template_test $(compiler_flags) $(root_flags) $(includes)

clean:
	rm -f $(ODIR)/*

.PHONY: all clean

