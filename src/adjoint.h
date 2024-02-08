#ifndef ADJOINT_H
#define ADJOINT_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <thread>
#include <chrono>
#include <random>

#include "tracer.h"

std::default_random_engine gen(time(0));

void progressBackProp(int progress, int total, int width = 50) {

    float percentage = (float) progress / total;
    int bar = (int) (percentage * width);

    std::cout << "[";

    for (int i = 0; i < bar; ++i) {
        std::cout << "=";
    }

    for (int i = bar; i < width; ++i) {
        std::cout << " ";
    }

    std::cout << "] " << std::setw(3) << (int) (percentage * 100) << "%\r backpropagating... ";
    
    std::cout.flush();
}


class Adjoint
{
    private:

        std::vector<double> h, u, d;
        std::vector<double> n_optim, n_optim_h, ndry;

        double s, t;
        double u_i;
        int steps;
        int size;

        double loss;
        double n_rms;
        double h_target;

        std::vector<double> pos;

        std::vector<double> init_pos;
        std::vector<double> target_pos;

        std::vector<std::vector<double>> forward;

        double h_end, u_end, d_end;

        double lam, mu;

        std::vector<double>::iterator i_low, i_up;

        double dw0, dw1;
        double dn0, dn1;

        double perturb;
        double AoA;

    public:

        Adjoint(const std::vector<double>& h0, const std::vector<double>& u0, const std::vector<double>& dmax, std::vector<double>& n, std::vector<double>& ndry, std::vector<double>& n_h)

        : h(h0), u(u0), d(dmax), size(h.size()), n_optim(n), n_optim_h(n_h), ndry(ndry)
        {}

        std::vector<double> retrieve(double obs_height, int n_iter, double lrate, double dr);
        std::vector<double> retrieve_synthetic(double obs_height, int n_iter, double lrate, double dr, std::vector<double>& n_target, double noise);


};

std::vector<double> Adjoint::retrieve(double obs_height, int n_iter, double lrate, double dr)
{   

    Tracer tracer;

    for (int i(0); i < n_iter; i++)
    {   
     
        loss = 0.0;

        for (int j(0); j < size; j++)
        {

            init_pos = tracer.trace(obs_height, u[j], d[j], dr, n_optim, n_optim_h);

            loss += pow((init_pos[0] - h[j]), 2);

            lam = 2*(init_pos[0] - h[j]);
            mu = 0.0;

            h_end = init_pos[0];
            u_end = -init_pos[1];
            d_end = init_pos[2];

        }

        //n_optim = trace.backprop(h_end, u_end, d_end, dr, n_optim, n_optim, n_optim_h, lam, mu, i, lrate);        

        std::cout << "\n " << loss << std::endl;

    }

    return n_optim;


};

std::vector<double> Adjoint::retrieve_synthetic(double obs_height, int n_iter, double lrate, double dr, std::vector<double>& n_target, double noise)
{   

    Tracer tracer;

    std::cout << "Noise standard deviation: " << noise << "\n";

    std::cout << "Obs Height" << ' ' << obs_height << std::endl;

    target_pos.clear();

    std::normal_distribution<double> distribution (0.0,noise);

    for (int k(0); k < size; k++)
    {
        target_pos.push_back(tracer.trace(obs_height, u[k], d[k], dr, n_target, n_optim_h)[0]);
    }

    if (noise > 0.0)
    {   

        std::cout << "Noise added" << std::endl;

        std::ostringstream filenoise;
        filenoise << "../Noise/PAPERII_noise_NE_RK3_NEW_" << noise << ".txt";
        std::ofstream rfilenoise(filenoise.str());

        for(int k(0); k < size; k++)
        {
            perturb = distribution(gen);

            AoA = asin(u[k]);
            AoA += perturb*PI/180.0;
            u[k] = sin(AoA);

            rfilenoise << perturb << std::endl; 

        }

        rfilenoise.close();

    }
    else if (noise < 0.0)
    {
        std::cerr << "Noise should be above or equal to 0.0 \n";
    }

    std::vector<double> m(n_optim.size(), 0.0);
    std::vector<double> v(n_optim.size(), 0.0);

    std::ostringstream filerms;
    filerms << "../RMS/PAPERII_RMS_NE_RK3_NEW_" << noise << ".txt";
    std::ofstream rfilerms(filerms.str());

    for (int i(0); i < n_iter; i++)
    {   

        //lrate = lrate + (0.5e-9 - lrate)*1.0*i/n_iter;

        loss = 0.0;
        n_rms = 0.0;
        for (int jj(0); jj < size; jj++)
        {
            loss += pow((tracer.trace(obs_height, u[jj], d[jj], dr, n_optim, n_optim_h)[0] - target_pos[jj]), 2);
        }

        for (int jjj(0); jjj < n_optim.size(); jjj++)
        {
            n_rms += pow((exp(n_optim[jjj]) - exp(n_target[jjj]))*1e6, 2);
        }

        for (int j(0); j < size; j++)
        {

            progressBackProp(j, size);

            init_pos = tracer.trace(obs_height, u[j], d[j], dr, n_optim, n_optim_h);

            lam = 2*(init_pos[0] - target_pos[j]);
            mu = 0.0;

            h_end = init_pos[0];
            u_end = -init_pos[1];
            d_end = init_pos[2];

            tracer.backprop(h_end, u_end, d_end, dr, n_optim, ndry, n_target[0], n_optim_h, lam, mu, lrate, i, m, v);
        

        }

        

        std::cout << "\n " << "iteration: " << (i+1) << ' ' << "loss: " << loss << " RMS N: " << sqrt(n_rms / n_optim.size()) << "\n" << std::endl;;            
        rfilerms << i << ' ' << loss << ' ' << sqrt(n_rms / n_optim.size()) << std::endl;
    }

    rfilerms.close();

    return n_optim;

};



#endif
