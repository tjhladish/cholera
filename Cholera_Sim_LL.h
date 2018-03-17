#ifndef CHOLERA_SIR_H
#define CHOLERA_SIR_H

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <math.h>
#include "DiffEq_Sim.h"
#include <vector>
#include <fstream>

using namespace std;

//const vector<double> RAIN = {(267.029+45.824+30.495), 26.68579, 169.00367, 275.74283, 418.06171, 91.20199, 140.94853, 140.94853, 456.00995, 39.24, 98.1, 215.82, 307.38, 18.086, 280.333, 253.204, 352.677, 39.537, 237.222, 263.58, 777.561, 68.464, 136.928, 179.718, 470.69, 11.488, 183.808, 333.152, 620.352, 164.065, 94.985, 198.605, 405.845, 30.348, 106.218, 227.61, 394.524, 10.091, 171.547, 272.457, 544.914, 40.552, 172.346, 263.588, 537.314, 16.198, 121.485, 251.069, 421.148, 110.964, 206.076, 467.634, 515.19, 13.714, 143.997, 198.853, 322.279};

//const vector<double> RAIN ={14.83620868,7.927745098,3.963872549,94.67992717,41.33752801,34.08930392	92.30160364	55.26770868	125.7113866	188.0008123	149.3813683	81.99553501

vector<double> RAIN;


//const vector<double> MEAN_RAIN = {8.407,14.214,28.9,46.386,76.279,43.7,60.8,94.15,102.729,207.964,196.021,74.892};//assuming rainfall starts in January

class Cholera_Sim_LL : public DiffEq_Sim {

    private:
        const double b;
        const double theta;
        const double gamma_A;
        const double gamma_S;
        const double mu;
        const double rho_A;
        const double rho_S;
        double meanRain;
        const double avgWeeksinMonth = 4.33;

    public:
        Cholera_Sim_LL() : b(0.0), theta(0.0), gamma_A(0.0), gamma_S(0.0), mu(0.0), rho_A(0.0), rho_S(0.0) { nbins=4; }
        Cholera_Sim_LL(double _b, double _theta, double _gamma_A, double _gamma_S, double _mu, double _rho_A, double _rho_S):
                    b(_b), theta(_theta), gamma_A(_gamma_A), gamma_S(_gamma_S), mu(_mu), rho_A(_rho_A), rho_S(_rho_S) { nbins=4; }
        ~Cholera_Sim_LL() {};

        void initialize(double S, double I, double RA, double RS) {
            x = new double[nbins];
            x[0] = S; 
            x[1] = I; 
            x[2] = RA;
            x[3] = RS;
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
        RAIN.insert(RAIN.end(), RAIN.begin(), RAIN.end());
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
            const double RA = x[2];
            const double RS = x[3];
            const double N = S+I+RA+RS;
            const double _t = get_time();
            const double rain = getRain(_t);
            //const double weekly_rain = getWeeklyRain(_t); //uncomment for weekly
            const double waterTerm = meanRain/(rain+meanRain);
            const double choleraDeath = 0;
            //const double waterTerm = (meanRain/avgWeeksinMonth)/(weekly_rain+perturbation); //uncomment for weekly

            dxdt[0] = (b*N) - (waterTerm*theta*S) - (mu*S) + (rho_A*RA) + (rho_S*RS);
            dxdt[1] = (waterTerm*theta*S) - (gamma_A*I) - (gamma_S*I) - (choleraDeath*I) - (mu*I);
            dxdt[2] = (gamma_A*I)- (rho_A*RA) - (mu*RA);
            dxdt[3] = (gamma_S*I) - (rho_S*RS) - (mu*RS);
        }

};

#endif
