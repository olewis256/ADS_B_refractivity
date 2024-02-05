#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstdio>

#include "tracer.h"
#include "integrate.h"
#include "inputs.h"
#include "adjoint.h"

#include "spline.h"

double n_lev = 30;
double pi = acos(-1);

std::vector<double> n_h(n_lev, 0);
std::vector<double> n(n_lev, 0);
std::vector<double> n_init(n_lev, 0);
std::vector<double> n_target(n_lev, 0);


int main()
{
    Atmosphere atmosphere_input("Refractivity_data/22Sep_12z_Watnall_profile.txt");
    atmosphere_input.process();

    const std::vector<double>& n_test = atmosphere_input.N();
    const std::vector<double>& h_test = atmosphere_input.H();

    ADSB adsb_input("ADS_B_data/sep_NE_paperII_input_5.000_central10.txt", 1000);
    adsb_input.process();

    const std::vector<double> sin_obsAoA_test = adsb_input.sin_obsAoA();
    const std::vector<double> obsAoA_test = adsb_input.obsAoA();
    const std::vector<double> H_test = adsb_input.rH();
    const std::vector<double> d_test = adsb_input.rD();
   
    Tracer rayTracer;

    tk::spline interp_n(h_test,n_test,tk::spline::cspline,true);

    double h0 = 0.577;

    for(int i(0); i < n_lev; i++)
    {

        n_h[i] = h0 + exp(2.7*i/n_lev) - 1.0;

        n[i] = log(1.000 + interp_n(n_h[0])*exp(-(n_h[i] - h0)/8.0)/1e6);

        n_init[i] = log(1.000 + interp_n(n_h[0])*exp(-(n_h[i] - h0)/8.0)/1e6);
        std::cout << interp_n(n_h[i]) << std::endl;
        n_target[i] = log(1.000 + interp_n(n_h[i])/1e6);

    }

    double dr = 0.1;

    Adjoint populate(H_test, sin_obsAoA_test, d_test, n, n_h);

    std::vector<double> result = populate.retrieve_synthetic(h0,10, 5e-11, dr, n_target);

    for(int i(0); i < n_lev; i++){
        std::cout << (exp(result[i]) - 1.0)*1e6 << ' ' << (exp(n_target[i]) - 1.0)*1e6 << std::endl;
    }

    std::ostringstream file3;
    file3 << "Profiles/PAPERII_retrieve_NE_RK3_NEW.txt";
    std::ofstream myfile3(file3.str());

    for(int i(0); i < n_lev; i++)
    {
    myfile3 << (exp(result[i]) - 1)*1e6 << ' ' << n_h[i] << ' '  << (exp(n_target[i]) - 1)*1e6 << ' ' << (exp(n_init[i]) - 1)*1e6 << std::endl;// ' ' << N_dry[i] << ' ' << Temp[i] << ' ' << Rh[i] << endl;
    }

    myfile3.close();

    return 0;
}

