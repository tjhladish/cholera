#ifndef CHOLERA_SIR_H
#define CHOLERA_SIR_H

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <math.h>
#include "DiffEq_Sim.h"
#include <array>
#include <fstream>

using namespace std;

//const vector<double> RAIN = {(267.029+45.824+30.495), 26.68579, 169.00367, 275.74283, 418.06171, 91.20199, 140.94853, 140.94853, 456.00995, 39.24, 98.1, 215.82, 307.38, 18.086, 280.333, 253.204, 352.677, 39.537, 237.222, 263.58, 777.561, 68.464, 136.928, 179.718, 470.69, 11.488, 183.808, 333.152, 620.352, 164.065, 94.985, 198.605, 405.845, 30.348, 106.218, 227.61, 394.524, 10.091, 171.547, 272.457, 544.914, 40.552, 172.346, 263.588, 537.314, 16.198, 121.485, 251.069, 421.148, 110.964, 206.076, 467.634, 515.19, 13.714, 143.997, 198.853, 322.279};

//const vector<double> RAIN ={14.83620868,7.927745098,3.963872549,94.67992717,41.33752801,34.08930392	92.30160364	55.26770868	125.7113866	188.0008123	149.3813683	81.99553501

vector<double> RAIN;


//const vector<double> MEAN_RAIN = {8.407,14.214,28.9,46.386,76.279,43.7,60.8,94.15,102.729,207.964,196.021,74.892};//assuming rainfall starts in January

class Cholera_Sim : public DiffEq_Sim {

    private:
        const double b;
        const double beta;
        const double C;
        const double epsilon;
        const double gamma;
        const double mu;
        const double rho;
        const double perturbation;
        double meanRain;
        const double avgWeeksinMonth = 4.33;

    public:
        Cholera_Sim() : b(0.0), beta(0.0), C(0.0), epsilon(0.0), gamma(0.0), mu(0.0), rho(0.0), perturbation(0.0) { nbins=4; }
        Cholera_Sim(double _b, double _beta, double _C, double _epsilon, double _gamma, double _mu, double _rho, double _perturbation):
                    b(_b), beta(_beta), C(_C), epsilon(_epsilon), gamma(_gamma), mu(_mu), rho(_rho), perturbation(_perturbation) { nbins=4; }
        ~Cholera_Sim() {};

        void initialize(double S, double I, double Y, double R) {
            x = new double[nbins];
            x[0] = S; 
            x[1] = I; 
            x[2] = Y;
            x[3] = R;
        }
    
    vector<double> getCompartment(){
        vector<double> comp;
        for(unsigned int i = 0; i < nbins - 1; i++){
            comp.push_back(x[i]);
        }
        return comp;
    }
    
    double readRain(){
        double sumRain = 0.0;
        if(ifstream in {"Vellore_monthly_rainfall.txt"}){
            double rain;
            while(in >> rain){
                RAIN.push_back(rain);
                sumRain+=rain;
            }
        }
        return meanRain = sumRain/RAIN.size();
    }
    
    double getRain(double t){
        int rain_idx = (int)t;
        return RAIN[rain_idx];
    }
    
    double getWeeklyRain(double t){
        int rain_idx = t/avgWeeksinMonth;
        return RAIN[rain_idx]/avgWeeksinMonth;
    }
    

        void derivative(double const x[], double dxdt[]) {
            const double S = x[0];
            const double I = x[1];
            const double Y = x[2];
            const double R = x[3];
            const double N = S+I+Y+R;
            const double _t = get_time();
            const double rain = getRain(_t);
            //const double weekly_rain = getWeeklyRain(_t); //uncomment for weekly
            const double waterTerm = meanRain/(rain+perturbation);
            //const double waterTerm = (meanRain/avgWeeksinMonth)/(weekly_rain+perturbation); //uncomment for weekly

            dxdt[0] = (b*N) - waterTerm*(beta*S*(I+Y)/N) - (mu*S) + (rho*Y) + (epsilon*R);
            dxdt[1] = (waterTerm*C*beta*S*(I+Y)/N) - (gamma*I) - (mu*I);
            dxdt[2] = (waterTerm*(1-C)*beta*S*(I+Y)/N) - (rho*Y) - (mu*Y);
            dxdt[3] = (gamma*I) - (epsilon*R) - (mu*R);
        }

};

#endif
