#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstdio>

#include "tracer.h"
#include "inputs.h"
#include "adjoint.h"
#include "constants.h"

#include "spline.h"

int n_lev = 120;
double h0 = 0.577;
double dr = 0.1;
double ue, he, se;

std::vector<double> n_h(n_lev, 0);
std::vector<double> n(n_lev, 0);
std::vector<double> n_init(n_lev, 0);
std::vector<double> n_target(n_lev, 0);

std::vector<double> forward;
std::vector<double> backward;

Tracer rayTracer;

int main()
{

    Atmosphere atmosphere_input("Refractivity_data/22Sep_12z_Watnall_profile.txt");
    atmosphere_input.process();

    const std::vector<double>& n_test = atmosphere_input.N();
    const std::vector<double>& h_test = atmosphere_input.H();

    tk::spline interp_n(h_test,n_test,tk::spline::cspline,true);

    for(int i(0); i < n_lev; i++)
    {

        n_h[i] = h0 + exp(2.7*i/n_lev) - 1.0;

        n[i] = log(1.000 + interp_n(n_h[0])*exp(-(n_h[i] - h0)/8.0)/1e6);

        n_init[i] = log(1.000 + interp_n(n_h[0])*exp(-(n_h[i] - h0)/8.0)/1e6);
        
        n_target[i] = log(1.000 + interp_n(n_h[i])/1e6);

    }

    forward = rayTracer.trace(h0, sin(1.0*PI/180.0), 300, dr, n_target, n_h);
    ue = -forward[1];
    he = forward[0];
    se = forward[2];

    //backward = rayTracer.trace(he, ue, 299, dr, n_target, n_h, false);

    std::ostringstream file2;
    file2 << "PAPERII_paths1.txt";
    std::ofstream myfile2(file2.str());
  
    std::vector<std::vector<double>> paths = rayTracer.trace_paths(h0, sin(1.0*PI/180.0), 300, dr, n_target, n_h);

    std::vector<std::vector<double>> pathsb = rayTracer.trace_paths(he, ue, se, dr, n_target, n_h, false);
    std::cout << paths.size() << ' ' << pathsb.size() <<std::endl;
    for(int i(0); i < pathsb.size(); i++)
    {
    myfile2 << paths[i][0] << ' ' << paths[i][1]  << ' ' << paths[i][2] << ' ' << pathsb[i][0] << ' ' << pathsb[i][1]  << ' ' << pathsb[i][2] << std::endl;// ' ' << N_dry[i] << ' ' << Temp[i] << ' ' << Rh[i] << endl;
    }

    myfile2.close();


    return 0;
}