#include "AbcSmc.h"
#include "Cholera_Sim.h"
#include <unistd.h>

vector<double> simulator(vector<double> args, const unsigned long int rng_seed, const unsigned long int serial, const ABC::MPI_par* mp) {
    const double b       = 0.02*1/12;//birth rate
    const double beta    = atof(args[0]);//10//transmission rate
    const double C       = atof(args[1]);//0.0002;//prob of symptomatic infection
    
    const double epsilon = 30.0/583;//waning rate from symptomatic infection
    const double gamma   = 30.0/14;//recovery rate from symptomatic infection
    const double mu      = 0.02*1/12;//death rate
    const double rho     = 30.0/102.2;//recovery from asymptomatic infection
    const double perturbation = 1;//makes denominator of rain function nonzero
    
    //const double N       = 790590.0/4;
    const double S       = 33509.20478;
    const double I       = 0.4518533614;
    const double Y       = 164119.6172;
    const double R       = 18.22614054;
    
    Cholera_Sim sim(b, beta, C, epsilon, gamma, mu, rho, perturbation);
    sim.initialize(S, I, Y, R);
    
    const double meanRain = sim.readRain(); //read in rain from txt file and get mean of all rainfall
    //const double weeklyMeanRain = meanRain/4.33; //uncomment for weekly
    
    const int max_time   = 168;
    //variables to initialize second simulation
    double S_2;
    double I_2;
    double Y_2;
    double R_2;
    
    //burn-in simulation
    for (unsigned int i = 0; i < max_time; ++i) {
        if(i == max_time - 1){
            S_2 = sim.getCompartment()[0];
            I_2 = sim.getCompartment()[1];
            Y_2 = sim.getCompartment()[2];
            R_2 = sim.getCompartment()[3];
        }
        sim.step_simulation(1);
        //const double rain = meanRain/(sim.getRain(i)+perturbation);
        //cerr << rain << endl;
    }
    sim.reset_time();
    sim.initialize(S_2,I_2,Y_2,R_2);
    
    for(unsigned int i = 0; i < max_time; ++i){
        sim.printX();
        sim.step_simulation(1);
    }
}

void usage(){
    cout<<"\n\t Usage: Input model parameters <beta> <c>";
    exit(-1);
}

int main(int argc, char* argv[]) {
    
    if (argc!=3) usage();
    
    vector<double> initialValues(argc-1);
    for(int i = 1; i<argc; i++){
        initialValues[i-1] = atof(argv[i]);
    }
    
    bool process_db = false;
    bool simulate_db = false;
    int buffer_size = -1;
    
    for (int i=2; i < argc;  i++ ) {
        if ( strcmp(argv[i], "--process") == 0  ) {
            process_db = true;
        } else if ( strcmp(argv[i], "--simulate") == 0  ) {
            simulate_db = true;
            buffer_size = buffer_size == -1 ? 1 : buffer_size;
        } else if ( strcmp(argv[i], "-n" ) == 0 ) {
            buffer_size = atoi(argv[++i]);
        } else {
            usage();
            exit(101);
        }
    }
    AbcSmc* abc = new AbcSmc();
    abc->parse_config(string(argv[1]));
    if (process_db) {
        gsl_rng_set(RNG, time(NULL) * getpid()); // seed the rng using sys time and the process id
        abc->process_database(RNG);
    }
    
    if (simulate_db) {
        abc->set_simulator(simulator);
        abc->simulate_next_particles(buffer_size);
    }
    
    return 0;
}

/*int main() {
    const double b       = 0.02*1/12;//birth rate
    const double beta    = 10.0;//transmission rate
    const double C       = 0.0002;//prob of symptomatic infection
    const double epsilon = 30.0/583;//waning rate from symptomatic infection
    const double gamma   = 30.0/14;//recovery rate from symptomatic infection
    const double mu      = 0.02*1/12;//death rate
    const double rho     = 30.0/102.2;//recovery from asymptomatic infection
    const double perturbation = 1;//makes denominator of rain function nonzero

    //const double N       = 790590.0/4;
    const double S       = 33509.20478;
    const double I       = 0.4518533614;
    const double Y       = 164119.6172;
    const double R       = 18.22614054;

    Cholera_Sim sim(b, beta, C, epsilon, gamma, mu, rho, perturbation);
    sim.initialize(S, I, Y, R);
    
    const double meanRain = sim.readRain(); //read in rain from txt file and get mean of all rainfall

    const int max_time   = 168;
    //variables to initialize second simulation
    double S_2;
    double I_2;
    double Y_2;
    double R_2;

    //burn-in simulation
    for (int i = 0; i < max_time; ++i) {
        if(i == max_time - 1){
            S_2 = sim.getCompartment()[0];
            I_2 = sim.getCompartment()[1];
            Y_2 = sim.getCompartment()[2];
            R_2 = sim.getCompartment()[3];
        }
        sim.step_simulation(1);
        //const double rain = meanRain/(sim.getRain(i)+perturbation);
        //cerr << rain << endl;
    }
    sim.reset_time();
    sim.initialize(S_2,I_2,Y_2,R_2);
   
    for(int i = 0; i < max_time; ++i){
        sim.printX();
        sim.step_simulation(1);
    }


    return 0;
}*/
