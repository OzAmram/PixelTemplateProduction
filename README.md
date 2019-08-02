# Pixel Template Production
Code for production of pixel templates for CMS

## Compiling

Everything is run standalone from the cmssw environment but makes use of some of the code in
it. 

Run the script **fetch\_cmssw\_code.sh** grabs the latest version of all the needed pixel code from the [cmssw github](https://github.com/cms-sw/cmssw).
If you want to grab from a branch other than the cmssw master (eg to test some
changes), you can change the `branch` variable in the script to point to a different branch.

The code shared with CMSSW uses the [vdt](https://github.com/dpiparo/vdt) library so you must have it installed.
For instructions to install vdt see [their github](https://github.com/dpiparo/vdt).
The Makefile assumes it installed to /usr/local/, if it is installed to some
other location, change the `vdt_dir` variable in the Makefile to point to the correct location. 



All of the source code is in the src/ directory which contains a Makefile. So you should be able to compile by simply changing to the src/ directory and running `make`. All the compiled executables are put in the bin/ directory. 




## Making Templates

There are two simple bash scripts that run the necessary executables to make templates: **make\_1d\_templates.sh** and **make\_2d\_templates.sh**. 

They should be run inside a directory containing pixelav events and a config file named `pix_2t.proc`. They take 1 argument, which is the location of the bin/ directory. 

Many files are produced when making templates but for the 1d production the real output is one template file named `template_summary_zpXXXX.out`and one gen errors file named `generror_summary_zpXXXX.out`. For 2d production it is one 2d template file named `template_summary2D_zpXXXX.out`. (The XXXX will be the starting file number in your config). 

The format `pix_2t.proc` is as follows:

> start\_file nfiles noise thresh1 thresh2 thresh1\_noise_frac common\_noise\_frac gain\_noise\_frac readout\_noise frontend\_type

> use\_l1\_offset write\_header


> id NTy NTyx NTxx DType Bfield VBias temp fluenc qscale 

Note that NTy is not used by the 2D templates so its value doesn't matter, but to keep the format consistent something must be there. 
Extra parameters on any of the lines will be ignored. 

An example config for 1D barrel templates is: 

> 58401 205 250. 1600. 1600. 0.073 0.080 0.080 350. 1

> 0 1

> 900 60 5 29 0 3.8 125. 253. 0. 1.




## Description of Executables
**gen\_xy\_template** : Makes 1D projections of average charge distributions from pixelav events. Also does charge variance fits. Makes files like `ptemp_XXXX.txt` and `ztemp_XXXX.txt`.

**gen\_xy\_template2d** : Makes 2d projections of average charge distributions from pixelav events.  Makes files like `zptemp_XXXX.txt`.

**gen\_zp\_template**: Uses the 1D projections and pixelav events. Runs the generic and template reco (using CMSSW code) on pixelav events to get resolution different algorithms and compute corrections to be saved in templates. Outputs one template file named `template_summary_zpXXXX.out`and one gen errors file named `generror_summary_zpXXXX.out`.

**gen\_zp\_template2d**: Uses the 1D and 2D projections and pixelav events. Runs the 2D template reco (using CMSSW code) on pixelav events to get resolution. Outputs one 2d template file named `template_summary2D_zpXXXX.out`.

**compare\_templates**: Takes in the file names of two templates and checks that all numerical values are the same within some threshold (default is 10^-5). It lists any discrepancies with the line number for investigation. Useful for testing changes. 

**test_template**: Uses pre-made 1D templates to run local version of CMSSW 1D template reco and makes various plots. Useful for test a new set of templates. 
