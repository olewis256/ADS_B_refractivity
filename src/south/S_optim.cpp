#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstdio>
#include <string.h> 

#include "../tracer.h"
#include "../inputs.h"
#include "../adjoint.h"
#include "../constants.h"

#include "../../externals/spline.h"

const double n_lev = 30;
std::string input_atmo;

double learn_rate;
int iterations;
int t;
const double dr = 0.1;

std::vector<double> n_h(n_lev, 0);
std::vector<double> logn(n_lev, 0);
std::vector<double> logndry(n_lev, 0);
std::vector<double> logn_init(n_lev, 0);
std::vector<double> logn_target(n_lev, 0);

std::vector<double> u_adsb;
std::vector<double> h_adsb;
std::vector<double> d_adsb;
std::vector<double> r_adsb;
std::vector<double> azim_adsb;
std::vector<double> arc_adsb;
std::vector<double> lat_adsb;

int main()
{

    std::cout << "time start t0 (t0 + 15mins slot): ";
    std::cin >> t;
    std::cout << "learning rate (normally Xe-9): ";
    std::cin >> learn_rate;
    std::cout << "iterations: ";
    std::cin >> iterations;

    // Read in refractivity data

    Atmosphere atmosphere_input("../../refractivity/Oct_25_N_16.00_(km)_updated_orog.txt");
    atmosphere_input.process();

    const std::vector<double>& N_profile = atmosphere_input.N();
    const std::vector<double>& NDRY_profile = atmosphere_input.NDRY();
    const std::vector<double>& h_profile = atmosphere_input.H();

    // Read in ADS-B data

    std::string adsbfile = "../../ADS_B_data/oct_S_paperII_input_central5_t" + std::to_string(t) + "_" + std::to_string(t+900) + ".txt";

    ADSB adsb_input(adsbfile, 5000);
    adsb_input.process();

    std::vector<double> u_adsb_I = adsb_input.sin_obsAoA();
    const std::vector<double> o_adsb_I = adsb_input.obsAoA();
    const std::vector<double> r_adsb_I = adsb_input.repAoA();
    const std::vector<double> azim_adsb_I = adsb_input.azim();
    const std::vector<double> h_adsb_I = adsb_input.rH();
    const std::vector<double> d_adsb_I = adsb_input.rD();
    const std::vector<double> arc_adsb_I = adsb_input.arc();
    const std::vector<double> lat_adsb_I = adsb_input.lat();

    Tracer rayTracer;

    tk::spline interp_n(h_profile,N_profile,tk::spline::cspline,true);
    tk::spline interp_ndry(h_profile,NDRY_profile,tk::spline::cspline,true);

    for(int i(0); i < n_lev; i++)
    {

        n_h[i] = OBSERVER_H + exp(2.7*i/n_lev) - 1.0;

        logn[i] = log(1.000 + interp_n(n_h[i])/1e6);//log(1.000 + interp_n(n_h[0])*exp(-(n_h[i] - OBSERVER_H)/9)/1e6);
        logndry[i] = 0.0;//log(1.000 + interp_ndry(n_h[i])/1e6);
        logn_init[i] = log(1.000 + interp_n(n_h[i])/1e6);//log(1.000 + interp_n(n_h[0])*exp(-(n_h[i] - OBSERVER_H)/9)/1e6);

        logn_target[i] = log(1.000 + interp_n(n_h[i])/1e6);

    }

    for(int j(0); j < (int) o_adsb_I.size(); j++)
    {
        if(d_adsb_I[j] < 600)
        {
            u_adsb.push_back(u_adsb_I[j]);
            d_adsb.push_back(EARTH_R * arc_adsb_I[j]);
            h_adsb.push_back(h_adsb_I[j]);
            azim_adsb.push_back(azim_adsb_I[j]);
            r_adsb.push_back(r_adsb_I[j]);
            
        }
        
    
    }

    std::cout << "Size of updated array: " << u_adsb.size() << std::endl;

    Adjoint Profile(h_adsb, u_adsb, d_adsb, logn, logndry, n_h, 0.0);

    std::vector<std::vector<double>> flightpath = Profile.retrieve_paths(OBSERVER_H, iterations, learn_rate, dr);

    std::ostringstream file;
    file << "../../flightpaths/PAPERII_retrieve_paths_oct_S_t" << t << "_" << t+900 << ".txt";
    std::ofstream rfile(file.str());

    for(int i(0); i < (int) u_adsb.size(); i++)
    {
    rfile << asin(u_adsb[i])*180.0/PI << ' ' << h_adsb[i] << ' '  << d_adsb[i] << ' ' << flightpath[0][i] << ' ' << flightpath[1][i] << ' ' << r_adsb[i] << ' ' << azim_adsb[i] <<  std::endl;
    }

    rfile.close();

    std::ostringstream file2;
    file2 << "../../retrievals/PAPERII_retrieve_oct_S_TRUE_t" << t << "_" << t+900 << ".txt";
    std::ofstream rfile2(file2.str());

    for(int i(0); i < n_lev; i++)
    {
    rfile2 << (exp(flightpath[2][i]) - 1.0)*1e6 << ' ' << n_h[i] << ' '  << (exp(logn_target[i]) - 1.0)*1e6 << ' ' << (exp(logn_init[i]) - 1.0)*1e6 << ' ' << (exp(logndry[i]) - 1.0)*1e6 << std::endl;
    }

    rfile2.close();
    
    return 0;
}