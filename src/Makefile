root_flags= `root-config --cflags --ldflags --libs --glibs`
includes= -I. -I../cmssw_code/ -I/usr/local/include
execs= gen_zp_template2d gen_xy_template gen_xy_template2d

makeAll: $(execs)

gen_zp_template2d: gen_zp_template2d.cc
	g++ gen_zp_template2d.cc  -o gen_zp_template2d   $(root_flags) $(includes)

gen_xy_template: gen_xy_templateqd4_2t.cc
	g++ gen_xy_templateqd4_2t.cc  -o gen_xy_template $(root_flags) $(includes)

gen_xy_template2d: gen_xy_template2d_2t.cc
	g++ gen_xy_template2d_2t.cc  -o gen_xy_template2d $(root_flags) $(includes)

.PHONY: clean

clean:
	rm -rf *.so *.o gen_zp_template2d gen_xy_template
