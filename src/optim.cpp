#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstdio>
#include <string.h> 

#include "tracer.h"
#include "inputs.h"
#include "adjoint.h"
#include "constants.h"

#include "../externals/spline.h"

const double n_lev = 30;
std::string input_atmo;

double learn_rate;
int iterations;
int t;
int num_batch;
const double dr = 0.1;

double max_d, min_d;
double max_h, min_h;

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

int main(int argc, char* argv[])
{

    if(!strcmp(argv[1], "real"))
    {
        std::cout << "------------------------------------------------------------" << std::endl;
        std::cout << "Retrieves one profile for selected time sample." << std::endl;
        std::cout << "------------------------------------------------------------ \n \n" << std::endl;
        std::cout << "time start t0 (t0 + 15mins slot): ";
        std::cin >> t;
        std::cout << "learning rate (normally Xe-9): ";
        std::cin >> learn_rate;
        std::cout << "iterations: ";
        std::cin >> iterations;

        // Read in refractivity data

        Atmosphere atmosphere_input("refractivity/Sep21_N_14.00_(km)_updated_orog.txt");
        atmosphere_input.process();

        const std::vector<double>& N_profile = atmosphere_input.N();
        const std::vector<double>& NDRY_profile = atmosphere_input.NDRY();
        const std::vector<double>& h_profile = atmosphere_input.H();

        // Read in ADS-B data

        std::string adsbfile = "ADS_B_data/sep_NE_paperII_input_central10_t" + std::to_string(t) + "_" + std::to_string(t+900) + ".txt";

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

        tk::spline interp_n(h_profile,N_profile,tk::spline::cspline,true);
        tk::spline interp_ndry(h_profile,NDRY_profile,tk::spline::cspline,true);

        for(int i(0); i < n_lev; i++)
        {

            n_h[i] = OBSERVER_H + exp(2.7*i/n_lev) - 1.0;

            logn[i] = log(1.000 + interp_n(n_h[0])*exp(-(n_h[i] - OBSERVER_H)/9)/1e6);
            logndry[i] = log(1.000 + interp_ndry(n_h[i])/1e6);
            logn_init[i] = log(1.000 + interp_n(n_h[0])*exp(-(n_h[i] - OBSERVER_H)/9)/1e6);

            logn_target[i] = log(1.000 + interp_n(n_h[i])/1e6);

        }

        if((argc >= 3) && (!strcmp(argv[2], "block")))
        {
                
            std::cout << "Max height of blocked region: ";
            std::cin >> max_h;
            std::cout << "Min height of blocked region: ";
            std::cin >> min_h;
            std::cout << "Max distance of blocked region: ";
            std::cin >> max_d;
            std::cout << "Min distance of blocked region: ";
            std::cin >> min_d;

            for(int j(0); j < (int) o_adsb_I.size(); j++)
            {   
                if( ((h_adsb_I[j] < min_h) || (h_adsb_I[j] > max_h)) && ((d_adsb_I[j] < min_d) || (d_adsb_I[j] > max_d)) )
                {
                    u_adsb.push_back(u_adsb_I[j]);
                    d_adsb.push_back(EARTH_R * arc_adsb_I[j]);
                    h_adsb.push_back(h_adsb_I[j]);
                    azim_adsb.push_back(azim_adsb_I[j]);
                    r_adsb.push_back(r_adsb_I[j]);    
                }
            }
        }
    
        else
        {
            for(int j(0); j < (int) o_adsb_I.size(); j++)
            {   
                u_adsb.push_back(sin(asin(u_adsb_I[j]) -0.03*PI/180.0));
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
        file << "flightpaths/PAPERII_retrieve_paths_t" << t << "_" << t+900 << ".txt";
        std::ofstream rfile(file.str());

        for(int i(0); i < (int) u_adsb.size(); i++)
        {
        rfile << asin(u_adsb[i])*180.0/PI << ' ' << h_adsb[i] << ' '  << d_adsb[i] << ' ' << flightpath[0][i] << ' ' << flightpath[1][i] << ' ' << r_adsb[i] << ' ' << azim_adsb[i] <<  std::endl;
        }

        rfile.close();

        std::ostringstream file2;
        file2 << "retrievals/PAPERII_retrieve_NE_TRUE_t" << t << "_" << t+900 << ".txt";
        std::ofstream rfile2(file2.str());

        for(int i(0); i < n_lev; i++)
        {
        rfile2 << (exp(flightpath[2][i]) - 1.0)*1e6 << ' ' << n_h[i] << ' '  << (exp(logn_target[i]) - 1.0)*1e6 << ' ' << (exp(logn_init[i]) - 1.0)*1e6 << ' ' << (exp(logndry[i]) - 1.0)*1e6 << std::endl;
        }

        rfile2.close();

    }

    else if(!strcmp(argv[1], "real_batch"))
    {
        std::cout << "------------------------------------------------------------" << std::endl;
        std::cout << "Retrieves profiles for all four time slots and all batches." << std::endl;
        std::cout << "------------------------------------------------------------ \n \n" << std::endl;
        std::cout << "Number of batches: ";
        std::cin >> num_batch;
        std::cout << "learning rate (normally Xe-9): ";
        std::cin >> learn_rate;
        std::cout << "iterations: ";
        std::cin >> iterations;

        // Read in refractivity data

        Atmosphere atmosphere_input("refractivity/Sep21_N_14.00_(km)_updated_orog.txt");
        atmosphere_input.process();

        const std::vector<double>& N_profile = atmosphere_input.N();
        const std::vector<double>& NDRY_profile = atmosphere_input.NDRY();
        const std::vector<double>& h_profile = atmosphere_input.H();

        Tracer rayTracer;

        tk::spline interp_n(h_profile,N_profile,tk::spline::cspline,true);
        tk::spline interp_ndry(h_profile,NDRY_profile,tk::spline::cspline,true);

        for(int i(0); i < n_lev; i++)
        {

            n_h[i] = OBSERVER_H + exp(2.7*i/n_lev) - 1.0;

            logn[i] = log(1.000 + interp_n(n_h[0])*exp(-(n_h[i] - OBSERVER_H)/9)/1e6);
            logndry[i] = 0.0;//log(1.000 + interp_ndry(n_h[i])/1e6);
            logn_init[i] = log(1.000 + interp_n(n_h[0])*exp(-(n_h[i] - OBSERVER_H)/9)/1e6);

            logn_target[i] = log(1.000 + interp_n(n_h[i])/1e6);

        }

        // Read in ADS-B data

        for(int j(0); j < 4; j++)
        {
                
            t = (int) j * 900;

            for(int i(0); i < num_batch; i++)
            {

                std::string adsbfile = "ADS_B_data/sep_NE_paperII_input_central5_t" + std::to_string(t) + "_" + std::to_string(t+900) + "_sample_" + std::to_string(i) + ".txt";

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
                file << "flightpaths/PAPERII_retrieve_paths_low_t" << t << "_" << t+900 << "_" << i << ".txt";
                std::ofstream rfile(file.str());

                for(int i(0); i < (int) u_adsb.size(); i++)
                {
                rfile << asin(u_adsb[i])*180.0/PI << ' ' << h_adsb[i] << ' '  << d_adsb[i] << ' ' << flightpath[0][i] << ' ' << flightpath[1][i] << ' ' << r_adsb[i] << ' ' << azim_adsb[i] <<  std::endl;
                }

                rfile.close();

                std::ostringstream file2;
                file2 << "retrievals/PAPERII_retrieve_NE_TRUE_low_t" << t << "_" << t+900 << "_" << i << ".txt";
                std::ofstream rfile2(file2.str());

                for(int i(0); i < n_lev; i++)
                {
                rfile2 << (exp(flightpath[2][i]) - 1.0)*1e6 << ' ' << n_h[i] << ' '  << (exp(logn_target[i]) - 1.0)*1e6 << ' ' << (exp(logn_init[i]) - 1.0)*1e6 << ' ' << (exp(logndry[i]) - 1.0)*1e6 << std::endl;
                }

                rfile2.close();

                u_adsb.clear();
                d_adsb.clear();
                h_adsb.clear();
                azim_adsb.clear();
                r_adsb.clear();
            }
        }
    }

    if(!strcmp(argv[1], "real_paths"))
    {

        Tracer tracer;

        std::cout << "------------------------------------------------------------" << std::endl;
        std::cout << "Retrieves one paths for selected time sample." << std::endl;
        std::cout << "------------------------------------------------------------ \n \n" << std::endl;
        std::cout << "time start t0 (t0 + 15mins slot): ";
        std::cin >> t;

        // Read in refractivity data

        Atmosphere atmosphere_input("refractivity/Sep21_N_14.00_(km)_updated_orog.txt");
        atmosphere_input.process();

        const std::vector<double>& N_profile = atmosphere_input.N();
        const std::vector<double>& NDRY_profile = atmosphere_input.NDRY();
        const std::vector<double>& h_profile = atmosphere_input.H();

        // Read in ADS-B data

        std::string adsbfile = "ADS_B_data/sep_NE_paperII_input_central10_t" + std::to_string(t) + "_" + std::to_string(t+900) + ".txt";

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

        tk::spline interp_n(h_profile,N_profile,tk::spline::cspline,true);
        tk::spline interp_ndry(h_profile,NDRY_profile,tk::spline::cspline,true);

        for(int i(0); i < n_lev; i++)
        {

            n_h[i] = OBSERVER_H + exp(2.7*i/n_lev) - 1.0;

            logn[i] = log(1.000 + interp_n(n_h[0])*exp(-(n_h[i] - OBSERVER_H)/9)/1e6);
            logndry[i] = log(1.000 + interp_ndry(n_h[i])/1e6);
            logn_init[i] = log(1.000 + interp_n(n_h[0])*exp(-(n_h[i] - OBSERVER_H)/9)/1e6);

            logn_target[i] = log(1.000 + interp_n(n_h[i])/1e6);

        }

        if((argc >= 3) && (!strcmp(argv[2], "block")))
        {
                
            std::cout << "Max height of blocked region: ";
            std::cin >> max_h;
            std::cout << "Min height of blocked region: ";
            std::cin >> min_h;
            std::cout << "Max distance of blocked region: ";
            std::cin >> max_d;
            std::cout << "Min distance of blocked region: ";
            std::cin >> min_d;

            for(int j(0); j < (int) o_adsb_I.size(); j++)
            {   
                if( ((h_adsb_I[j] < min_h) || (h_adsb_I[j] > max_h)) && ((d_adsb_I[j] < min_d) || (d_adsb_I[j] > max_d)) )
                {
                    u_adsb.push_back(u_adsb_I[j]);
                    d_adsb.push_back(EARTH_R * arc_adsb_I[j]);
                    h_adsb.push_back(h_adsb_I[j]);
                    azim_adsb.push_back(azim_adsb_I[j]);
                    r_adsb.push_back(r_adsb_I[j]);    
                }
            }
        }
    
        else
        {
            for(int j(0); j < (int) o_adsb_I.size(); j++)
            {   
                u_adsb.push_back(u_adsb_I[j]);
                d_adsb.push_back(EARTH_R * arc_adsb_I[j]);
                h_adsb.push_back(h_adsb_I[j]);
                azim_adsb.push_back(azim_adsb_I[j]);
                r_adsb.push_back(r_adsb_I[j]);     
            }    
        
        }

        std::cout << "Size of updated array: " << u_adsb.size() << std::endl;

        std::vector<std::vector<double>> path;

        std::ostringstream file;
        file << "flightpaths/PAPERII_ray_paths" << t << "_" << t+900 << ".txt";
        std::ofstream rfile(file.str());

        for(int j(0); j < (int) u_adsb.size(); j++)
        {
            path = tracer.trace_paths(OBSERVER_H, u_adsb[j], d_adsb[j], dr, logn_target, n_h);
            for(int i(0); i < (int) path.size(); i++)
            {
                rfile << j << ' ' << path[i][0] << ' ' << path[i][1] << '\n';
                std::cout << j << ' ' << i << std::endl;
            }
        }
        rfile.close();

    }

    else
    {
    
        std::cout << "learning rate: (normally X*1e-9) ";
        std::cin >> learn_rate;
        std::cout << "iterations: ";
        std::cin >> iterations;

        // Read in refractivity data
        if(!strcmp(argv[1], "1"))
        {
            input_atmo = "refractivity/22Sep_12z_Watnall_profile.txt";
        }
        if(!strcmp(argv[1], "2"))
        {
            input_atmo = "refractivity/18Jul_12z_2022_Watnall_profile_RH.txt";
        }
        if(!strcmp(argv[1], "3"))
        {
            input_atmo = "refractivity/15Dec_12z_2022_Watnall_profile_RH.txt";
        }
        
        Atmosphere atmosphere_input(input_atmo);
        atmosphere_input.process();

        const std::vector<double>& N_profile = atmosphere_input.N();
        const std::vector<double>& NDRY_profile = atmosphere_input.NDRY();
        const std::vector<double>& h_profile = atmosphere_input.H();

        // Read in ADS-B data

        ADSB adsb_input("ADS_B_data/sep_NE_paperII_input_5.000_central10.txt", 5000);
        adsb_input.process();

        const std::vector<double> u_adsb = adsb_input.sin_obsAoA();
        const std::vector<double> o_adsb = adsb_input.obsAoA();
        const std::vector<double> h_adsb = adsb_input.rH();
        const std::vector<double> d_adsb = adsb_input.rD();
    
        Tracer rayTracer;

        tk::spline interp_n(h_profile,N_profile,tk::spline::cspline,true);
        tk::spline interp_ndry(h_profile,NDRY_profile,tk::spline::cspline,true);

        double noise[3] = {0.00, 0.01, 0.05};

        for(int i(0); i < 3; i++)
        {

            for(int i(0); i < n_lev; i++)
            {

                n_h[i] = OBSERVER_H + exp(2.7*i/n_lev) - 1.0;

                logn[i] = log(1.000 + interp_n(n_h[0])*exp(-(n_h[i] - OBSERVER_H)/8)/1e6);
                logndry[i] = log(1.000 + interp_ndry(n_h[i])/1e6);
                logn_init[i] = log(1.000 + interp_n(n_h[0])*exp(-(n_h[i] - OBSERVER_H)/8)/1e6);

                logn_target[i] = log(1.000 + interp_n(n_h[i])/1e6);

            }

            Adjoint Profile(h_adsb, u_adsb, d_adsb, logn, logndry, n_h, std::stoi(argv[1]));

            // Tracing and learning parameters
            double noise_std = noise[i];

            std::vector<double> retrieval = Profile.retrieve_synthetic(OBSERVER_H, iterations, learn_rate, dr, logn_target, noise_std);

            std::ostringstream file;
            file << "retrievals/PAPERII_retrieve_NE_RK3_" << argv[1] << "_" << noise_std << ".txt";
            std::ofstream rfile(file.str());

            for(int i(0); i < n_lev; i++)
            {
            rfile << (exp(retrieval[i]) - 1.0)*1e6 << ' ' << n_h[i] << ' '  << (exp(logn_target[i]) - 1.0)*1e6 << ' ' << (exp(logn_init[i]) - 1.0)*1e6 << std::endl;
            }

            rfile.close();
        }
    }

    return 0;
}