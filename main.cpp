#include "AbcSmc.h"
#include "Cholera_Sim.h"
#include <fstream>
#include <vector>

const gsl_rng* RNG = gsl_rng_alloc(gsl_rng_taus2);

vector<double> simulator(vector<double> args, const unsigned long int rng_seed, const unsigned long int serial, const ABC::MPI_par* mp) {
    const double b       = 0.02*1/12;//birth rate
    const double beta    = (double)(args[0]);//10//transmission rate
    const double C       = (double)(args[1]);//0.0002;//prob of symptomatic infection
    
    const double epsilon = 30.0/583;//waning rate from symptomatic infection
    const double gamma   = 30.0/14;//recovery rate from symptomatic infection
    const double mu      = 0.02*1/12;//death rate
    const double rho     = 30.0/102.2;//recovery from asymptomatic infection
    
    //const double N       = 790590.0/4;
    const double S       = 33509.20478;
    const double I       = 0.4518533614;
    const double Y       = 164119.6172;
    const double R       = 18.22614054;
    
    Cholera_Sim sim(b, beta, C, epsilon, gamma, mu, rho);
    sim.initialize(S, I, Y, R);
    
    const double meanRain = sim.readRain(); //read in rain from txt file and get mean of all rainfall
    //const double weeklyMeanRain = meanRain/4.33; //uncomment for weekly
    vector <double> monthlyCases;
    ifstream cases;
    cases.open("Vellore_monthly_cases.txt");
    double num;
    
    while(cases){
        cases >> num;
        monthlyCases.push_back(num);
    }
 
    const int max_time   = 168;
    vector<double> metrics(168);
    
    //burn-in simulation -- get rid --> make rainfall vector twice as long
    for (unsigned int i = 0; i < 2*max_time; ++i) {
        if(i >= max_time){
            //sim.printX();
            metrics[i-168] = gamma*sim.getCompartment()[1];
        }
        sim.step_simulation(1);
        //const double rain = meanRain/(sim.getRain(i)+perturbation);
        //cerr << sim.getRain(i) << endl;
    }
    return metrics;
}

void usage(){
    cout<<"\n\t Usage: Input model parameters <beta> <c>";
    exit(-1);
}

int main(int argc, char* argv[]) {
    
    if (not(argc == 3 or argc == 5 or argc == 6)) usage();
    
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

    //const double N       = 790590.0/4;
    const double S       = 33509.20478;
    const double I       = 0.4518533614;
    const double Y       = 164119.6172;
    const double R       = 18.22614054;

    Cholera_Sim sim(b, beta, C, epsilon, gamma, mu, rho);
    sim.initialize(S, I, Y, R);
    
    const double meanRain = sim.readRain(); //read in rain from txt file and get mean of all rainfall

    vector <double> monthlyCases;
    ifstream cases;
    cases.open("Vellore_monthly_cases.txt");
    double num;
        while(cases){
            cases >> num;
            monthlyCases.push_back(num);
        }
    cout<<"length of monthly cases "<<monthlyCases.size()<<"\n";
    for(int i = 0; i < monthlyCases.size();i++){
        cout<<"monthly cases["<<i<<"] "<<monthlyCases[i]<<"\n";
    }
    
    const int max_time   = 168;

    for (unsigned int i = 0; i < 2*max_time; ++i) {
        if(i >= max_time){
            
            //cerr<<sim.getCompartment()[1]<<endl;
        }
        sim.step_simulation(1);
        //const double rain = meanRain/(sim.getRain(i)+perturbation);
        //cerr << sim.getRain(i) << endl;
    }


    return 0;
}*/
