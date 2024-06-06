
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
double h;

std::vector<double> n_h(n_lev, 0);
std::vector<double> logn(n_lev, 0);

int main()
{

    // Read in refractivity data

    Atmosphere atmosphere_input("../refractivity/Sep21_N_16.00_(km)_updated_orog.txt");
    atmosphere_input.process();

    const std::vector<double>& N_profile = atmosphere_input.N();
    const std::vector<double>& NDRY_profile = atmosphere_input.NDRY();
    const std::vector<double>& h_profile = atmosphere_input.H();

    // Read in ADS-B data

    ADSB adsb_input("../ADS_B_data/sep_NE_paperII_input_all_within02.txt", 263513);
    adsb_input.process();

    const std::vector<double> u_adsb = adsb_input.sin_obsAoA();
    const std::vector<double> o_adsb = adsb_input.obsAoA();
    const std::vector<double> r_adsb = adsb_input.repAoA();
    const std::vector<double> azim_adsb = adsb_input.azim();
    const std::vector<double> h_adsb = adsb_input.rH();
    const std::vector<double> d_adsb = adsb_input.rD();
    const std::vector<double> t_adsb = adsb_input.rTime();
   
    Tracer rayTracer;

    tk::spline interp_n(h_profile,N_profile,tk::spline::cspline);

    for(int i(0); i < n_lev; i++)
    {


        n_h[i] = exp(2.7*i/n_lev) - 1.0;

        if(n_h[i] < h_profile[0])
        {

            logn[i] = log(1.000 + interp_n(h_profile[0])*exp((h_profile[0]-n_h[i])/7)/1e6);
        
        }

        else
        {
            logn[i] = log(1.000 + interp_n(n_h[i])/1e6);
        }

        //logn[i] = log(1.000 + 300*exp(-n_h[i]/8)/1e6);

        std::cout << n_h[i] << ' ' << (exp(logn[i]) - 1)*1e6 << std::endl;

    }

    // Tracing 
    double dr = 0.1;

    std::ostringstream file;
    file << "../flightpaths/PAPERII_synobs_all.txt";
    std::ofstream rfile(file.str());

    for(int i(0); i < (int) u_adsb.size(); i++)

    {
    std::cout << 100.0*i/(u_adsb.size()) << std::endl;

    h = rayTracer.trace(OBSERVER_H, u_adsb[i], d_adsb[i], dr, logn, n_h)[0];
    
    rfile << asin(u_adsb[i])*180.0/PI << ' ' << h_adsb[i] << ' ' << h << ' ' << d_adsb[i] << ' ' << t_adsb[i] << ' ' << r_adsb[i] << ' ' << azim_adsb[i] <<  std::endl;
    }

    rfile.close();



    return 0;
}

