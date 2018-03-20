#include <fstream>
#include <vector>
#include <cstring>

#ifdef USE_LL
#include "Cholera_Sim_LL.h"
#endif

#ifdef USE_KING
#include "Cholera_Sim.h"
#endif

using namespace std;

int main() {
    const double b       = 0.02*1/12;//birth rate
    const double mu      = 0.02*1/12;//death rate
    const double N       = 790590.0/4; //total population
    
    const int max_time   = 168;
    
    double recoveryRate;
    
#ifdef USE_LL

        const double theta      = 0.0014805;
        const double gamma_A    = 30.0/1;//asymptomatic recovery rate
        const double gamma_S    = 30.0/14; //symptomatic recovery rate
        const double rho_A      = 30.0/102.2; //asympomatic waning rate
        const double rho_S      = 30.0/583; //symptomatic waning rate
    
        recoveryRate = gamma_S;
    
        //set population distribution
        const double S         = 64912.89976;
        const double I         = 6.782346297;
        const double RA        = 689.2424262;
        const double RS        = 273.5754724;
        
        Cholera_Sim_LL sim(b,theta,gamma_A,gamma_S,mu,rho_A,rho_S);
        sim.initialize(S,I,RA,RS);

#endif
#ifdef USE_KING
        const double beta    = 608.67;//transmission rate
        const double C       = 0.000047823;//prob of symptomatic infection
        const double epsilon = 30.0/583;//waning rate from symptomatic infection
        const double gamma   = 30.0/14;//recovery rate from symptomatic infection
        const double rho     = 30.0/102.2;//recovery from asymptomatic infection
    
        recoveryRate = gamma;
        
        //set population distribution
        const double S       = 33509.20478;
        const double I       = 0.4518533614;
        const double Y       = 164119.6172;
        const double R       = 18.22614054;
        
        Cholera_Sim sim(b, beta, C, epsilon, gamma, mu, rho);
        sim.initialize(S, I, Y, R);
#endif
    
    const double meanRain = sim.readRain(); //read in rain from txt file and get mean of all rainfall
    
    vector <double> monthlyCases;
    ifstream cases;
    cases.open("Vellore_monthly_cases.txt");
    double num;
    while(cases){
        cases >> num;
        monthlyCases.push_back(num);
    }
    
    for (unsigned int i = 0; i < 2*max_time; ++i) {
        if(i >= max_time){
            
            cout<<recoveryRate*sim.getCompartment()[1]<<endl;
        }
        sim.step_simulation(1);
    }


    return 0;
}
