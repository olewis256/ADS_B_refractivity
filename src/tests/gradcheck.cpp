#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstdio>

#include "../tracer.h"
#include "../inputs.h"
#include "../adjoint.h"
#include "../constants.h"

int n_lev = 20;
double dr = 0.1;

double u0, h0, s0;
double ue, he, se;

int index_n;
double delta_n;
double adj_n;
double lam;

double hpos, hneg;
double error;

std::vector<double> n_h(n_lev, 0);
std::vector<double> n(n_lev, 0);
std::vector<double> n1(n_lev, 0);
std::vector<double> n0(n_lev, 0);
std::vector<double> n_init(n_lev, 0);
std::vector<double> n_target(n_lev, 0);

std::vector<double> forward;
double forward_target;
std::vector<double> backward;

Tracer rayTracer;

int main()
{      

    std::ostringstream file;
    file << "../../Gradient/PAPERII_grad_" << dr << "km.txt";
    std::ofstream rfile(file.str());

    for (int k(0); k < n_lev; k++)
    {
        
        index_n = k;
        adj_n = 0.0;

        for(int i(0); i < n_lev; i++)
        {

            n_h[i] = exp(2.7*i/n_lev) - 1.0;

            n[i] = log(1.000 + 320*exp(-n_h[i]/8.0)/1e6);
            n1[i] = log(1.000 + 320*exp(-n_h[i]/8.0)/1e6);
            n0[i] = log(1.000 + 320*exp(-n_h[i]/8.0)/1e6);

            n_init[i] = log(1.000 + 320*exp(-n_h[i]/8.0)/1e6);
            
            n_target[i] = log(1.000 + 315*exp(-n_h[i]/8)/1e6);

        }

        delta_n = 1e-6;

        n1[index_n] = log(1.000 + ((exp(n1[index_n]) - 1.0)*1e6 + delta_n)/1e6 );
        n0[index_n] = log(1.000 + ((exp(n0[index_n]) - 1.0)*1e6 - delta_n)/1e6 );

        h0 = 0.01;
        u0 = sin(0.2*PI/180.0);
        s0 = 350;

        forward_target = rayTracer.trace(h0, u0, s0, dr, n_target, n_h)[0];

        forward = rayTracer.trace(h0, u0, s0, dr, n, n_h);

        ue = -forward[1];
        he = forward[0];
        se = forward[2];

        lam = 2*(he - forward_target);

        rayTracer.backprop(he, ue, se, dr, n, n[0], n_h, lam, 0.0, 0.0, &index_n, &adj_n, &h0);

        hpos = rayTracer.trace(h0, u0, s0, dr, n1, n_h)[0] - forward_target;
        hneg = rayTracer.trace(h0, u0, s0, dr, n0, n_h)[0] - forward_target;

        error = (hpos*hpos - hneg*hneg) / (2*delta_n*1e-6);

        std::cout << adj_n << ' ' << error << ' ' << 100*(adj_n - error)/error << std::endl;
        
        rfile << adj_n << ' ' << error << ' ' << n_h[k] << ' ' << k << std::endl;
        
    }
    rfile.close();
    return 0;
}