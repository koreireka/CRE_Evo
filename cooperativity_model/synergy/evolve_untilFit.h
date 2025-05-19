#ifndef EVOLVE
#define EVOLVE

//////////////////////////////////////////////////////////////////////////////////
///     Functions for the evolution of the promoter                            ///
//////////////////////////////////////////////////////////////////////////////////

double fitness(double alpha, double g[N][NE]){

    double f = 0.0;

    for(int i=0; i<N; i++){
        for(int j=0; j<NE; j++) f += pow(desired_output[i][j]-g[i][j],2);
    }

    f *= alpha;
    f /= (double)(NE);

    f = exp(-f);

    return f;
}


int *propose_mutation(int *m){

    int gene = rand()%N;
    int loc = rand()%LX;
    int new_letter = rand()%4;

    while (new_letter == X[gene][loc]) new_letter = rand()%4;

    m[0] = gene;            // which promoter
    m[1] = loc;             // where
    m[2] = X[gene][loc];    //previous letter

    X[gene][loc] = (char)new_letter;

    return m;

}


void reject_mutation(int *m){

    X[m[0]][m[1]] = m[2];

    return;
}


double rand_exp(double lambda){
    double u = drand48(); //change here !!!!
    return -log(u)/lambda;
}


void evolve(double T, double alpha, double Npop, double d0, double epsilon, double offset, double eta, double fitLevel,
            const char *fitness_out_filename,
            // const char *rho_out_filename,  do not save rho data -- takes up too much space, could calculate later if needed form sequences.
            const char *gene_exp_out_filename,
            const char *selection_out_filename,
            const char *sequences_out_filename){

    //X[N][LX], C[NT][LC] global variables

    FILE *f_fitness_out = fopen(fitness_out_filename,"w");
    // FILE *f_rho_out = fopen(rho_out_filename,"w");
    FILE *f_g_out = fopen(gene_exp_out_filename,"w");
    FILE *f_selection_out = fopen(selection_out_filename, "w");
    FILE *f_sequences_out = fopen(sequences_out_filename, "w");


    int mutation[3]; // which promoter, where on promoter, old number
    long double s, f0, f1, pfix;
    double rho[NE][NT];     // current interaction strength depends on the concentration now!!!
    double g[N][NE]; //current gene expression

    double tau;
    double lambda = Npop*LX; //1.0; This is now rescaled such that time units = 1/u 

    // initial steps
    for(int i=0; i<N; i++){
        for(int j=0; j<NE; j++) g[i][j] = 0.0;
    }

    tau = 0.0;

    //run network for all possible inputs
    run_network(rho, g, d0, epsilon, offset, eta);

    //calculate fitness
    f0 = fitness(alpha, g);


    // print_rho(f_rho_out, tau, rho);
    print_fitness(f_fitness_out, tau, f0);
    print_sequences(f_sequences_out);
    print_gene_expression(f_g_out, tau, g);

     //run evolution until T and first fitness hit
    while((tau < T) && (f0 < fitLevel)){

        tau += rand_exp(lambda);

        //propose mutation
        propose_mutation(mutation);

        //calculate mutated fitness
        run_network(rho, g, d0, epsilon, offset, eta);

        f1 = fitness(alpha, g);

        //calculate fixation probability

         if (f0 == 0.0) f0 = 0.000000000000001;
         s = (f1-f0)/f0;  //division by f0=0 cause problems!


//        //with very very strong selection
//        if( s > 0.0) pfix = 1.0;
//        else pfix = 0.0;

        if(s == 0) {
                pfix = 1.0/Npop; //corrected this too !
        }

        else {
                pfix = (1.0-exp(-2.0*s))/(1.0-exp(-2.0*Npop*s)); //corrected version!
        }

        //throw random number
        double u = drand48(); //(double)rand()/(double)RAND_MAX;

        if(u < pfix) {

            f0 = f1; //old fitness became the new one

            //only print things out if they change
            print_gene_expression(f_g_out, tau, g);
            // print_rho(f_rho_out,tau, rho);
            print_fitness(f_fitness_out, tau, f0);
            print_sequences(f_sequences_out);
            print_selection(f_selection_out, tau, Npop*s);


        }

        else {

            reject_mutation(mutation);
            // print out full selection trajectories
            // print_selection(f_selection_out, tau, Npop*s);


        }
    }


    fclose(f_fitness_out);
    // fclose(f_rho_out);
    fclose(f_g_out);
    fclose(f_selection_out);
    fclose(f_sequences_out);

    return;
}

#endif
