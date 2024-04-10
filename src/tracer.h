#ifndef TRACER_H
#define TRACER_H

#include <vector>

class Tracer
{
    private:

        double h, u;
        double s, t;
        int steps;
        int size;
        double r;

        double lam, mu;

        double dw0_1, dw1_1, dw0_2, dw1_2, dw0_3, dw1_3;

        double dn0, dn1;

        std::vector<double>::iterator i_lev1, i_lev2, i_lev3;

        std::vector<double> final_pos;
        std::vector<std::vector<double> > paths;

        double h_obs;

        //-----------------------
        // Runge Kutta gradients
        //-----------------------

        double n1, n2, n3;
        double k1h, k2h, k3h;
        double k1u, k2u, k3u;
        double k1lam, k2lam, k3lam;
        double k1mu, k2mu, k3mu;

        //---------------------------------
        // Adam optimisier hyperparameters
        //---------------------------------

        const double beta1 = 0.9;
        const double beta2 = 0.999;
        const double epsilon = 1e-8;

        double m_est = 0.0;
        double v_est = 0.0;

        double gamma = 0.2;

       
    public:

        //-------------------
        // Defining functions
        //--------------------

        std::vector<double> trace(const double h0, const double u0, const double d, const double dr_i, std::vector<double>& n, std::vector<double>& n_h, bool forward = true);

        std::vector<std::vector<double> > trace_paths(const double h0, const double u0, const double d, const double dr_i, std::vector<double>& n, std::vector<double>& n_h, bool forward = true);

        std::vector<double> backprop(const double h0, const double u0, const double d, const double dr_i, std::vector<double>& n, std::vector<double>& ndry, const double n_surface_true,
                      std::vector<double>& n_h, const double lam0, const double mu0, const double lrate, int iter, std::vector<double>& m, std::vector<double>& v,
                      int* index_n = nullptr, double* dn_adj = nullptr, double* obs_height = nullptr);

};

#endif