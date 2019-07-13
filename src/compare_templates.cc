/* compare_templates
 * Author Oz Amram June 2019
 *
 * Compares all the doubles in two files to see if they disagree beyond some
 * tolerance.
 * Ignores lines with anything that cannot be interpretted as a double
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

    double tolerance = 1e-5;

    if(argc > 3){
        tolerance = atof(argv[3]);
    }
    double val1, val2;
    const int max_line = 300;
    char b1[max_line], b2[max_line];
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
    while((l1 = fgets(b1,max_line, f1)) != NULL && (l2 = fgets(b2, max_line, f2)) != NULL){
        //printf("b1 %s, b2 %s \n", b1, b2);
        line_num++;
        s1 = std::string(b1);
        s2 = std::string(b2);
        ps1 = 0;
        ps2 = 0;
        //loop over all substrings in the line
        while(s1.substr(ps1).length() > 0 && s2.substr(ps2).length() > 0){

            try{
                //try to read a double from the string
                val1 = std::stod(s1.substr(ps1), &new1);
                val2 = std::stod(s2.substr(ps2), &new2);
                ps1+=new1;
                ps2+=new2;
                

                double diff = std::fabs((val1 - val2)/val1);
                //ignore numbers less than certain amount because too noisy
                double min_val = 1e-8;
                bool too_small = std::fabs(val1) < min_val && std::fabs(val2) < min_val;
                if(diff > tolerance && !too_small){
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
                // Couldn't convert to a double, likely a string
                //printf("Not doubles. Strings are %s %s", s1.substr(ps1).c_str(), s2.substr(ps2).c_str());
                ps1 += s1.substr(ps1).length();
                ps2 += s2.substr(ps2).length();
                continue;
            }
            catch(const std::out_of_range& ia){
                printf("Out of range?. Line is %i. %s %s \n", line_num, s1.c_str(), s2.c_str());
                printf("ps1, new1, ps2, new2: %d %d %d %d \n", (int) ps1,  (int)new1,  (int)new2,  (int)ps2);
                break;
            }


        checked_vals++;
        }
    }
    //check if f2 is finished as well (can lazy evaluate the conditional) 
    l2 = fgets(b2, 180, f2);
            
    fclose(f1);
    fclose(f2);
   
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
