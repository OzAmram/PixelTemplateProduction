
#!/bin/sh
# A simple script to do both steps of making 1D templates
# Should be run in a directory with pixelav events (named template_events_dXXXXX.out) 
# Requires a config file named pix_2t.proc in the same directory

# Note that first two steps (gen_xy_template and gen_xy_template2d) can actually be run in parallel 
# So if you are pressed for time you can do that

exec_dir=../bin/

./${exec_dir}gen_xy_template
./${exec_dir}gen_xy_template2d
./${exec_dir}gen_zp_template2d
