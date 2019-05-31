### Pixel Template Production
Code for production of pixel templates for CMS


In order to compile properly you must have [vdt](https://github.com/dpiparo/vdt) installed.
To instructions to install vdt see [their github]((https://github.com/dpiparo/vdt).
The Makefile assumes it installed to /usr/local/, if it is installed to some
other location, change the `vdt_dir` variable to point to the correct location. 


Everything is run standalone from the cmssw environment but makes use of some of the code in
it. 
The script **fetch\_cmssw\_code.sh** grabs the latest version of all the needed pixel code from the [cmssw github](https://github.com/cms-sw/cmssw).
If you want to grab from a branch other than the cmssw master (eg to test some
changes), you can change the `branch` variable in the script to point to a different branch
