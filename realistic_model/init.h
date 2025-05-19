#ifndef INIT
#define INIT

//Contains all stuff for initialization of the code
//Contains also all static variables

double input_table[NT][NE]; //one TF with on and off mode = 2 environments, NE=number of environments
double desired_output[N][NE]; //desired gene expression in each environment
char C[NT][LC]; //consensus sequence for TF
char X[N][LX]; //N promoter sequence



//////////////////////////////////////////////////////////////////////////////////
///     Function for one_network.c, these get called there as initialization    ///
//////////////////////////////////////////////////////////////////////////////////

int numberOfMatches(char *c1, char *c2){

    int res = 0;

    for(int i=0; i<LC; i++) 
        if(c1[i] == c2[i]) res++;

    return res;
}

void initialize_sequences(){

    int j;

    for (int i=0; i<N; i++){
        for(j=0; j<LX; j++) X[i][j] = rand()%4;
    }

    // generalized approach for arbitrary number of TFs, they should differ at least in 3 position (naive approach)
    for (int i=0; i<NT; i++)
    {

        int success=0;
        while(!success)
        {
            char test[LC];

            for(j=0; j<LC; j++) test[j]= rand()%4;

            int valid = 1;
            for(int ii = 0; ii<i; ii++)
            {
                if(numberOfMatches(C[ii], test) >= 4 ) {valid = 0; break;}
            }
            if(valid)
            {
                success=1;
                for(j=0; j<LC; j++) C[i][j] = test[j];
            }
        }
    }

    return;

}


void read_initial_sequences(const char *fname){

    FILE *f = fopen(fname, "r");

    for(int i=0; i<NT; i++){
        for(int j=0; j<LC; j++) fscanf(f, "%1d ", &C[i][j]);
    }
    
    
    for (int i=0; i<N; i++)
        for(int j=0; j<LX; j++) fscanf(f, "%1d ", &X[i][j]);

    fclose(f);

}


void read_input_table(const char *fname){

    FILE *f = fopen(fname, "r");

    for(int i=0; i<NT; i++){

        for(int j=0; j<NE; j++) fscanf(f,"%lf ",&input_table[i][j]);
        fscanf(f,"\n");
    }

    fclose(f);
}


void read_desired_output(const char *fname){

    FILE *f = fopen(fname, "r");

     for(int i=0; i<N; i++){

        for(int j=0; j<NE; j++) fscanf(f,"%lf ",&desired_output[i][j]);
        fscanf(f,"\n");
    }

    fclose(f);
}


//////////////////////////////////////////////////////////////////////////////////
///     Print out functions during the evolution of the promoters              ///
//////////////////////////////////////////////////////////////////////////////////

void print_rho(FILE *out, double t, double rho[NE][NT]){

    fprintf(out,"%.10lf ",t);

    for(int i=0; i<NE; i++)
        for(int j=0; j<NT; j++)
            fprintf(out,"%.10lf ",rho[i][j]);

    fprintf(out,"\n");
}

void print_gene_expression(FILE *out, double t, double g[N][NE]){

    fprintf(out,"%.10lf ",t);
    for(int i=0; i<N; i++){
        for(int j=0; j<NE; j++) fprintf(out,"%.10lf ",g[i][j]);
    }
    fprintf(out,"\n");
}

void print_fitness(FILE *out, double t, double f){

    fprintf(out,"%.10lf %.15lf",t, f);
    fprintf(out,"\n");
}

void print_selection(FILE *out, double t, long double s){

    fprintf(out,"%lf %.15Lf",t, s);
    fprintf(out,"\n");
}


void print_sequences(FILE *out){

    // TF motifs
    for(int i=0; i<NT; i++){
        for(int j=0; j<LC; j++) fprintf(out, "%d",C[i][j]);
    fprintf(out, " ");
    }

     //Promoters
   for(int i=0; i<N; i++){
        for(int j=0; j<LX; j++) fprintf(out, "%d", X[i][j]);
        fprintf(out, " ");
    }


    fprintf(out,"\n");


}


#endif
