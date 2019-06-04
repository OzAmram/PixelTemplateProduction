/* compare_templates
 * Author Oz Amram June 2019
 *
 * Compares all the floats in two files to see if they disagree beyond some
 * tolerance.
 * Ignores lines with anything that cannot be interpretted as a float
 * Very useful for comparing two sets of templates, to see if something has
 * changed when changing the code
 *
 * Requires 2 file names to compare.
 * Takes optional 3rd argument, tolerance, which sets the threshhold to be
 * considerered 'disagreeing'. Default is 10^-5.
 */

#include "template_utils.h"
#include <set>
#include <exception>

int main(int argc, char *argv[]){
    if (argc < 3){
        printf("You must give two files to compare. Optional 3rd argument is tolerance \n");
        return 1;
    }

    FILE *f1 = fopen(argv[1], "r");
    FILE *f2 = fopen(argv[2], "r");

    float tolerance = 1e-5;

    if(argc > 3){
        sscanf(argv[3], " %f ", &tolerance);
    }
    float val1, val2;
    int nargs1, nargs2;
    char b1[180], b2[180];
    int checked_vals = 0;
    int line_num=0;
    int fails =0;
    int success = 0;
    std::string s1, s2;
    std::set<int> bad_lines;
    
    //pointers to substrings
    std::string::size_type ps1, ps2, new1, new2;
    char *l1;
    char *l2;


    //get line by line
    while((l1 = fgets(b1,180, f1)) != NULL && (l2 = fgets(b2, 180, f2)) != NULL){
        //printf("b1 %s, b2 %s \n", b1, b2);
        line_num++;
        s1 = std::string(b1);
        s2 = std::string(b2);
        ps1 = 0;
        ps2 = 0;
        //loop over all substrings in the line
        while(s1.substr(ps1).length() > 0 && s2.substr(ps2).length() > 0){

            try{
                //try to read a float from the string
                val1 = std::stof(s1.substr(ps1), &new1);
                val2 = std::stof(s2.substr(ps2), &new2);
                ps1+=new1;
                ps2+=new2;
                

                float diff = std::fabs((val1 - val2)/val1);
                if(diff > tolerance){
                    printf("Difference between files found on line %i \n"
                           "%s has %f. " 
                           "%s has %f.\n \n", 
                           
                           line_num, argv[1], val1, argv[2], val2);
                    
                    bad_lines.insert(line_num);
                    fails++;
                }
                else{
                    success++;
                
                }
            }
            catch(const std::invalid_argument& ia){
                // Couldn't convert to a float, likely a string
                //printf("Not floats. Strings are %s %s", s1.substr(ps1).c_str(), s2.substr(ps2).c_str());
                break;
            }


        checked_vals++;
        }
    }
    //check if f2 is finished as well (can lazy evaluate the conditional) 
    l2 = fgets(b2, 180, f2);
            
   
    if(l1 == NULL && l2 == NULL)
        printf("End of both files reached. \n"); 
    else{
        printf("only one of the files is done \n"
               "File 1 is done (1=yes, 0=no): %i"
               " File 2 is done (1=yes, 0=no): %i \n",
               l1 == NULL, l2 == NULL);
    }


    printf("Read %i lines. Fails: %i. Success %i \n", line_num, fails, success);
    if(bad_lines.size() > 0){
        printf("Lines with discrepancies are: ");
        for(auto iter = bad_lines.begin(); iter != bad_lines.end(); iter++){
            printf("%i ", *iter);
        }
        printf("\n");
    }


}
