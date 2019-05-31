root_flags= `root-config --cflags --ldflags --libs --glibs`
compiler_flags=-std=c++11 -Wall -O0 -g
includes= -I. -I../cmssw_code/ -I/usr/local/include
execs= gen_zp_template gen_zp_template2d gen_xy_template #gen_xy_template2d

makeAll: $(execs)


gen_zp_template: gen_zp_template.cc
	g++ gen_zp_template.cc -o gen_zp_template $(compiler_flags)  $(root_flags) $(includes)

gen_zp_template2d: gen_zp_template2d.cc
	g++ gen_zp_template2d.cc -o gen_zp_template2d $(compiler_flags)  $(root_flags) $(includes)

gen_xy_template: gen_xy_templateqd4_2t.cc
	g++ gen_xy_templateqd4_2t.cc   -o gen_xy_template $(compiler_flags) $(root_flags) $(includes)

gen_xy_template2d: gen_xy_template2d_2t.cc
	g++ gen_xy_template2d_2t.cc  -o gen_xy_template2d $(compiler_flags) $(root_flags) $(includes)

.PHONY: clean

clean:
	rm -f $(execs)
