#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstdio>

#include "tracer.h"
#include "inputs.h"
#include "adjoint.h"
#include "constants.h"

#include "../externals/spline.h"

double n_lev = 30;

std::vector<double> n_h(n_lev, 0);
std::vector<double> logn(n_lev, 0);
std::vector<double> logndry(n_lev, 0);
std::vector<double> logn_init(n_lev, 0);
std::vector<double> logn_target(n_lev, 0);


int main()
{

    // Read in refractivity data

    Atmosphere atmosphere_input("../refractivity/22Sep_12z_Watnall_profile.txt");
    atmosphere_input.process();

    const std::vector<double>& N_profile = atmosphere_input.N();
    const std::vector<double>& NDRY_profile = atmosphere_input.NDRY();
    const std::vector<double>& h_profile = atmosphere_input.H();

    // Read in ADS-B data

    ADSB adsb_input("../ADS_B_data/sep_NE_paperII_input_5.000_central10.txt", 5000);
    adsb_input.process();

    const std::vector<double> u_adsb = adsb_input.sin_obsAoA();
    const std::vector<double> o_adsb = adsb_input.obsAoA();
    const std::vector<double> h_adsb = adsb_input.rH();
    const std::vector<double> d_adsb = adsb_input.rD();
   
    Tracer rayTracer;

    tk::spline interp_n(h_profile,N_profile,tk::spline::cspline,true);
    tk::spline interp_ndry(h_profile,NDRY_profile,tk::spline::cspline,true);

    for(int i(0); i < n_lev; i++)
    {

        n_h[i] = OBSERVER_H + exp(2.7*i/n_lev) - 1.0;

        logn[i] = log(1.000 + interp_n(n_h[0])*exp(-(n_h[i] - OBSERVER_H)/8)/1e6);
        logndry[i] = log(1.000 + interp_ndry(n_h[i])/1e6);
        logn_init[i] = log(1.000 + interp_n(n_h[0])*exp(-(n_h[i] - OBSERVER_H)/8)/1e6);

        logn_target[i] = log(1.000 + interp_n(n_h[i])/1e6);

    }

    Adjoint Profile(h_adsb, u_adsb, d_adsb, logn, logndry, n_h);

    // Tracing and learning parameters
    double dr = 0.1;
    double learn_rate =50e-9;
    int iterations =20;
    double noise_std = 0.05;

    std::vector<double> retrieval = Profile.retrieve_synthetic(OBSERVER_H, iterations, learn_rate, dr, logn_target, noise_std);

    std::ostringstream file;
    file << "../Retrievals/PAPERII_retrieve_NE_RK3_NEW_" << noise_std << ".txt";
    std::ofstream rfile(file.str());

    for(int i(0); i < n_lev; i++)
    {
    rfile << (exp(retrieval[i]) - 1.0)*1e6 << ' ' << n_h[i] << ' '  << (exp(logn_target[i]) - 1.0)*1e6 << ' ' << (exp(logn_init[i]) - 1.0)*1e6 << std::endl;
    }

    rfile.close();

    return 0;
}