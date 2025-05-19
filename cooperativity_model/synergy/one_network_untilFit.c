#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>


#define N 1 //number of promoter sequences

// #define NT 19 //number of TFs
// #define LC 8 //length of consensus sequences
// #define LX 256 //length of regulatory sequences
// #define NE 2 //number of environments; one TF with on and off mode = 2 environments


#include "../init.h"
#include  "../regnet.h"
#include "../evolve_untilFit.h"


int main(int argc, char *argv[]){

     //read evolutionary parameters
    double T  = atof(argv[1]);
    double alpha = atof(argv[2]);
    double Npop = atof(argv[3]);
    double d0 = atof(argv[4]);
    double epsilon = atof(argv[5]);
    double offset = atof(argv[6]);
    double eta = atof(argv[7]);
    double cmax = atof(argv[8]);
    double R = atof(argv[9]);
    double fitLevel = atof(argv[10]); //only run sims up to first hit of high fitness
    int ss = atoi(argv[11]);
    int seed = atoi(argv[12]);

    //initial things

    srand(seed);
    srand48(seed);
    
    // create input/output tables instead of reading in
    // ONLY WORKS FOR 2 ENV

    double cmin = cmax*pow(10.0, -R);
    
    //cognate TF
    input_table[0][0] = cmin;
    input_table[0][1] = cmax;

    //noncognate TF
    for(int i=1; i<NT; i++){
        input_table[i][0] = cmax;
        input_table[i][1] = cmin;

    }
    
    //works only for one gene!!
    desired_output[0][0] = 0.0;
    desired_output[0][1] =1.0;

    // read_input_table("../inputs/input_table_2env_xtalkingTFs18_coop.dat");

    // read_desired_output("../inputs/output_table_2env_xtalkingTFs18.dat");

    //argv[12] file name for saving output

    int size = strlen(argv[13])+100; 
    char temp[size];
    strcpy(temp, argv[13]);

    if (ss == 0){
        char f0[size];
        strcpy(f0, temp);
        strcat(f0, "/params.txt");
        FILE *par = fopen(f0, "w");
        fprintf(par, "T\t%lf\n alpha\t%lf\n Npop\t%lf\n d0\t%lf\n eps\t%lf\n offset\t%lf\n eta\t%lf\n R\t%lf\n motif length \t%d\n num of TFs \t%d\n num of env \t%d\n promoter length \t%d\n cmin cmax \n %.15lf %.15lf", T, alpha, Npop, d0, epsilon, offset,eta, R, LC, NT, NE, LX, cmin, cmax);
        // fprintf(par, "%s \t %d",temp, size);
        fclose(par);
    }


    char number[20];
    sprintf(number,"%d.dat",ss);
    
    //initialize_sequences();
    //read in initial sequences so that we have the same starting set

    char fInitSeq[200];
    sprintf(fInitSeq, "../inputs/initialSequences_l%d_L%d_xtalkingTFs%d/sequencesStart_", LC, LX, NT-1);
    strcat(fInitSeq, number);

    read_initial_sequences(fInitSeq);
    
    char f1[size];
    strcpy(f1,temp);

    // char f2[size];
    // strcpy(f2,temp);

    char f3[size];
    strcpy(f3,temp);

    char f4[size];
    strcpy(f4,temp);

    char f5[size];
    strcpy(f5,temp);

    strcat(f1,"/fitness/out_fitness"); strcat(f1,number);
    // strcat(f2,"/rho/out_rho"); strcat(f2,number);
    strcat(f3,"/gene_exp/out_output_table"); strcat(f3,number);
    strcat(f4,"/selection/out_selection"); strcat(f4,number);
    strcat(f5,"/sequences/out_sequence"); strcat(f5,number);
    
    //do not save rho data, upon saving put back here f2 as an argument!
    evolve(T,alpha, Npop, d0, epsilon, offset, eta, fitLevel, f1,f3,f4,f5);

    return 0;


}
