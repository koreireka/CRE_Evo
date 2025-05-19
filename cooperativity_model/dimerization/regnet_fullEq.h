#ifndef REGNET
#define REGNET

#include <stdbool.h>

//////////////////////////////////////////////////////////////////////////////////
///     Functions for calculating the gene expression of the promoter          ///
//////////////////////////////////////////////////////////////////////////////////



double calc_rho(int indexTF, double epsilon, double ec, double c){  //X[N][LX], C[NT][LC] global variables
    // Calculating rho
    // Optionally index for gene as well but for now we have only one promoter

    //Scan over promoter here and calculate contributions by index
    
    double sum = 0.0;
    int LCper2 = LC/2;

    
    for(int i=0; i<=LX-LC; i++){

        //first monomer contribution
        int k1 = 0;

        for(int j=0; j<LCper2; j++){

            if(X[i+j] != C[indexTF][j]) k1++;
        }

        //second monomer contribution
        int k2 = 0;

        for(int j=LCper2; j<LC; j++){

            if(X[i+j] != C[indexTF][j]) k2++;
        }

    //dimer cooperativity only at eta=1 l=12 with full partition function
    sum +=  ( (c*c) /( (c*c) + exp(epsilon*(k1+k2)+ec) + c*exp(epsilon*k1+ec) + c*exp(epsilon*k2+ec) ) );
        
    }
    

    return sum;

}


double sigmoid(double sum, double offset, double d0 ){
    // Calculating the sigmoid function
    return 1.0/(1.0+exp(-d0*sum+d0*offset));

}


void run_network(double rho[NE][NT], double g[NE], double d0, double epsilon, double offset, double ec){ //X[LX], C[NT][LC] global variables

    int j, k;

    for(j=0; j<NE; j++) g[j] =0.0;

    // Even though there is a for loop here over genes, rho structure only good for 1 gene now! +1 dimension needed for more
    // input_table[k][j] --> concentration of TF k in environment j (either cmin or cmax)
    
    for(j=0; j<NE; j++){

        double sum = 0.0;
        
        for(k=0; k<NT; k++){ 
            rho[j][k] = calc_rho(k, epsilon, ec, input_table[k][j]);
            sum += rho[j][k];
        }

        g[j] = sigmoid(sum, offset, d0);
    }


    return;
}

#endif
