#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstdio>
#include <string.h> 

#include "../tracer.h"
#include "../adjoint.h"
#include "../inputs.h"
#include "../constants.h"

#include "../../externals/spline.h"

const double n_lev = 30;
std::string input_atmo;

double learn_rate;
int iterations;
int t;
int num_batch;
const double dr = 0.1;

double max_d, min_d;
double max_h, min_h;

double n_err = 0.0;
double h_err = 0.0;
int h_err_input = 0;
double zero_point = 0.0;

double N_f, h_f;

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

std::vector<double> N;
std::vector<double> h;   

std::vector<double> forward;
double ue, he, se;

int main(int argc, char* argv[])
{

        std::cout << "------------------------------------------------------------" << std::endl;
        std::cout << "Retrieves one profile for selected time sample." << std::endl;
        std::cout << "------------------------------------------------------------ \n \n" << std::endl;
        std::cout << "time start t0 (t0 + 15mins slot): ";
        std::cin >> t;

        // Read in refractivity data

        std::ifstream input_f("ropp_profiles/ropp_profile.txt");
    
        while (input_f >> N_f >> h_f)
        {
                
                N.push_back(N_f);
                h.push_back(h_f);
        }     

        // Read in ADS-B data

        std::string adsbfile = "../../ADS_B_data/sep_NE_paperII_input_central10_t" + std::to_string(t) + "_" + std::to_string(t+900) + ".txt";

        ADSB adsb_input(adsbfile);
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

        tk::spline interp_n(h,N,tk::spline::cspline,true);

        for(int i(0); i < n_lev; i++)
        {

            n_h[i] = OBSERVER_H + exp(2.7*i/n_lev) - 1.0;

            logn[i] = log(1.000 + interp_n(n_h[i])/1e6);

        }

    
        for(int j(0); j < 2000; j++)
        {   
            u_adsb.push_back(u_adsb_I[j]);
            d_adsb.push_back(EARTH_R * arc_adsb_I[j]);
            h_adsb.push_back(h_adsb_I[j]);
            azim_adsb.push_back(azim_adsb_I[j]);
            r_adsb.push_back(r_adsb_I[j]);     
        }    

        std::ostringstream file;
        file << "ROPP_adsb.txt";
        std::ofstream myfile(file.str());
        
        for(int j(0); j < 2000; j++)
        {
    
            forward = rayTracer.trace(OBSERVER_H, u_adsb[j], d_adsb[j], dr, logn, n_h);
            ue = forward[1];
            he = forward[0];
            se = forward[2];

            myfile << ue << ' ' << he  << ' ' << se << ' ' << u_adsb[j] << std::endl;
     

        }

        myfile.close();
    

    return 0;
}