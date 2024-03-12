#ifndef ADJOINT_H
#define ADJOINT_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <thread>
#include <chrono>
#include <random>
#include <algorithm>

#include "constants.h"

class Adjoint
{

    //-----------------------------------------------------------------------------------------------------------------------
    //
    // Adjoint class.  Contains function algorithms for retrieving refractivity profiles using raw ADS-B observations.  
    // During each iteration in range n_iter, all observations are simulated using the observed AoA and known position of the
    // aircraft.  The loss function is evaluated using tracer.trace and backpropagated through the refractivity
    // field using tracer.backprop for each individual ray. The gradients of the loss function with respect to
    // each refractivity value in n is accumulated and then n updated using gradient descent.  The retrieval is
    // terminated after n_iter iterations.
    //
    // Each function has slight variations:
    //  +  .retrieve                   : retrieves refractivity profiles using raw ADS-B observations
    //  +  .retrieve_synthetic         : retrieves refractivity profiles using synthetic observations (input n field)
    //  +  .retrieve_paths             : initial and retrieved flightpaths and retrieved n field using raw ADS-B obs.
    //
    // Common inputs:
    //  +  obs_height                  : height of observer (initial height of ray)  
    //  +  n_iter                      : number of iterations to run the retrieval algorithm
    //  +  lrate                       : learning rate - controls speed of convergence (careful!)
    //
    //------------------------------------------------------------------------------------------------------------------------

    private:

        int argi;

        std::vector<double> h, u, d;
        std::vector<double> n_optim, n_optim_h, ndry;

        double s, t;
        double u_i;
        int steps;
        int size;

        double loss, loss_prev;
        double n_loss, n_sum, n_grad_sum;
        double n_rms;
        double h_target;

        std::vector<double> pos;
        std::vector<double> init_pos;
        std::vector<double> target_pos;

        std::vector<std::vector<double> > forward;

        double h_end, u_end, d_end;

        double lam, mu;

        std::vector<double>::iterator i_low, i_up;

        double dw0, dw1;
        double dn0, dn1;

        double perturb;
        double AoA;

        //---------------------------------
        // Adam optimisier hyperparameters
        //---------------------------------

        const double beta1 = 0.9;
        const double beta2 = 0.999;
        const double epsilon = 1e-8;

        double m_est = 0.0;
        double v_est = 0.0;


    public:

        //-------------------
        // Defining functions
        //--------------------

        Adjoint(const std::vector<double>& h0, const std::vector<double>& u0, const std::vector<double>& dmax,
                std::vector<double>& n, std::vector<double>& ndry, std::vector<double>& n_h, int argi);

        std::vector<double> retrieve(double obs_height, int n_iter, double lrate, double dr);

        std::vector<double> retrieve_synthetic(double obs_height, int n_iter, double lrate, double dr, std::vector<double>& n_target, double noise);

        std::vector<std::vector<double> > retrieve_paths(double obs_height, int n_iter, double lrate, double dr);

        void progressBackProp(int prog, int tot);


};



#endif
