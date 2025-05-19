#ifndef REGNET
#define REGNET

#include <stdbool.h>

//////////////////////////////////////////////////////////////////////////////////
///     Functions for calculating the gene expression of the promoter          ///
//////////////////////////////////////////////////////////////////////////////////


int *mismatches(char *x, char *c, int *d){
    // Calculating mismatches
    int j;
    for(int i=0; i<=LX-LC; i++){

        int n = 0;

        for(j=0; j<LC; j++){

            if(x[i+j] != c[j]) n++;

        }

        d[n]++;
    }


    return d;
}


double calc_rho(int indexTF, double epsilon, double ec, double c){  //X[N][LX], C[NT][LC] global variables
    // Calculating rho
    // Optionally index for gene as well but for now we have only one promoter
    
    int *d = (int*)calloc(LC+1, sizeof(int));

    for(int ii=0; ii<LC+1; ii++) d[ii] = 0;

    d = mismatches(X,C[indexTF],d);

    double sum = 0.0;
   
    //dimer cooperativity only at eta=1 l=12
    for(int n=0; n<LC+1; n++) sum +=  d[n] * ( (c*c) /( (c*c)+exp(epsilon*n+ec) ) ); // if + sign --> ec<0
            

    free(d);

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
