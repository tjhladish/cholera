#include "Cholera_Sim.h"

int main() { 
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
}
