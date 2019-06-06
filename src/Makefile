root_flags= `root-config --cflags --ldflags --libs --glibs`
compiler_flags=-std=c++11 -Wall -O0 -g
vdt_dir=/usr/local/include
includes= -I. -I../cmssw_code/ -I$(vdt_dir)
execs= compare_templates gen_zp_template gen_zp_template2d gen_xy_template #gen_xy_template2d
bin_dir=../bin

makeAll: $(execs)


gen_zp_template: gen_zp_template.cc template_utils.h
	g++ gen_zp_template.cc -o $(bin_dir)/gen_zp_template $(compiler_flags)  $(root_flags) $(includes)

gen_zp_template2d: gen_zp_template2d.cc template_utils.h
	g++ gen_zp_template2d.cc -o $(bin_dir)/gen_zp_template2d $(compiler_flags)  $(root_flags) $(includes)

gen_zp_template2d_old: gen_zp_template2d_2018_2t.cxx
	g++ gen_zp_template2d_2018_2t.cxx -o $(bin_dir)/gen_zp_template2d_old $(compiler_flags)  $(root_flags) $(includes)

gen_xy_template: gen_xy_templateqd4_2t.cc template_utils.h
	g++ gen_xy_templateqd4_2t.cc   -o $(bin_dir)/gen_xy_template $(compiler_flags) $(root_flags) $(includes)

gen_xy_template2d: gen_xy_template2d_2t.cc template_utils.h
	g++ gen_xy_template2d_2t.cc  -o $(bin_dir)/gen_xy_template2d $(compiler_flags) $(root_flags) $(includes)

compare_templates: compare_templates.cc template_utils.h
	g++ compare_templates.cc -o $(bin_dir)/compare_templates $(compiler_flags) $(root_flags) $(includes)

.PHONY: clean

clean:
	rm -f $(bin_dir)/*
