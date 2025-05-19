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

    //printf("%d ", n);
    }


    return d;
}


double calc_rho(int indexTF, double epsilon, double eta, double c){  //X[N][LX], C[NT][LC] global variables
    // Calculating rho
    // Optionally index for gene as well but for now we have only one promoter
    
    int *d = (int*)calloc(LC+1, sizeof(int));

    for(int ii=0; ii<LC+1; ii++) d[ii] = 0;

    d = mismatches(X[0],C[indexTF],d);

    double sum = 0.0;
    for(int n=0; n<LC+1; n++) sum +=  d[n] * ( c/(c+exp(epsilon*n)) );
            
    double rho = pow(sum, eta); // (-): repression (+): activation

    free(d);

    return rho;

}


double sigmoid(double sum, double offset, double d0 ){
    // Calculating the sigmoid function
    return 1.0/(1.0+exp(-d0*sum+d0*offset));

}


void run_network(double rho[NE][NT], double g[N][NE], double d0, double epsilon, double offset, double eta){ //X[N][LX], C[NT][LC] global variables

    int i,j, k;

    for(i=0; i<N; i++) {
    for(j=0; j<NE; j++) g[i][j] =0.0;
    }

    // Even though there is a for loop here over genes, rho structure only good for 1 gene now! +1 dimension needed for more
    for(i=0; i<N; i++){
        for(j=0; j<NE; j++){

            double sum = 0.0;
            
            for(k=0; k<NT; k++){ 
                rho[j][k] = calc_rho(k, epsilon, eta, input_table[k][j]);
                sum += rho[j][k];
            }

            g[i][j] = sigmoid(sum, offset, d0);
        }
    }

    return;
}

#endif
