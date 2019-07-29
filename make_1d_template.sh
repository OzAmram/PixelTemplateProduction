#!/bin/sh
# A simple script to do both steps of making 1D templates
# Should be run in a directory with pixelav events (named template_events_dXXXXX.out) 
# Requires a config file named pix_2t.proc in the same directory

exec_dir=../bin/

./${exec_dir}gen_xy_template
./${exec_dir}gen_zp_template
