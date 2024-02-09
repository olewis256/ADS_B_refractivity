#ifndef TRACER_H
#define TRACER_H

#include <vector>
#include <algorithm>
#include <cmath>

#include "constants.h"

class Tracer
{
    private:

        double h, u;
        double s, t;
        double u_i;
        int steps;
        int size;
        double r;

        double lam, mu;

        double dw0_1, dw1_1, dw0_2, dw1_2, dw0_3, dw1_3;

        double dn0, dn1;

        //-----------------------
        // Runge Kutta gradients
        //-----------------------

        double n1, n2, n3;
        double k1h, k2h, k3h;
        double k1u, k2u, k3u;
        double k1lam, k2lam, k3lam;
        double k1mu, k2mu, k3mu;

        std::vector<double>::iterator i_lev1, i_lev2, i_lev3;

        const double gamma = 0.8;

        std::vector<double> final_pos;
        std::vector<std::vector<double>> paths;

        double h_obs;

        //---------------------------------
        // Adam optimisier hyperparameters
        //---------------------------------

        const double beta1 = 0.9;
        const double beta2 = 0.999;
        const double epsilon = 1e-8;
        double m_est0 = 0.0;
        double m_est1 = 0.0;
        double v_est0 = 0.0;
        double v_est1 = 0.0;
       
    public:

        std::vector<double> trace(const double h0, const double u0, const double d, const double dr_i, std::vector<double>& n, std::vector<double>& n_h, bool forward);

        std::vector<std::vector<double>> trace_paths(double h0, double u0, double d, double dr_i, std::vector<double>& n, std::vector<double>& n_h, bool forward);

        void backprop(const double h0, const double u0, const double d, const double dr_i, std::vector<double>& n, std::vector<double>& ndry, const double n_surface_true, std::vector<double>& n_h,
                      const double lam0, const double mu0, const double lrate, int iter, std::vector<double>& m, std::vector<double>& v, int* index_n = nullptr, double* dn_adj = nullptr, double* obs_height = nullptr);


};

std::vector<double> Tracer::trace(const double h0, const double u0, const double dmax, const double dr,
                                                 std::vector<double>& n, std::vector<double>& n_h, bool forward = true)
{
    steps = (int) 600/dr;


    h = h0;
    u = u0;

    if (forward)
    {
        s = 0.0;
        t = 1.0;
    }
    else
    {
        s = dmax;
        t = -1.0;
    }

    for(int i(0); i < steps; i++)
    { 
        r = h + EARTH_R;

        i_lev1 = std::upper_bound(n_h.begin(), n_h.end(), h);

        dw0_1 = - 1  / ( n_h[i_lev1 - n_h.begin()] - n_h[i_lev1 - 1 - n_h.begin()] );
        dw1_1 = -1 * dw0_1;
            
        n1 = dw0_1*n[i_lev1 - 1 - n_h.begin()] + dw1_1*n[i_lev1 - n_h.begin()]; 


        k1h = dr*u;
        k1u = dr*( (1 - u*u) * (n1 + (1 / r)) );

        i_lev2 = std::upper_bound(n_h.begin(), n_h.end(), (h + k1h/2));

        dw0_2 = - 1  / ( n_h[i_lev2 - n_h.begin()] - n_h[i_lev2 - 1 - n_h.begin()] );
        dw1_2 = -1 * dw0_2;

        n2 = dw0_2*n[i_lev2 - 1 - n_h.begin()] + dw1_2*n[i_lev2 - n_h.begin()]; 

        
        k2h = dr*(u + k1u/2);
        k2u = dr*( (1 - (u + k1u/2)*(u + k1u/2))  * ( n2 + 1 / (r + k1h/2)));

        i_lev3 = std::upper_bound(n_h.begin(), n_h.end(), (h - k1h + 2*k2h));

        dw0_3 = - 1  / ( n_h[i_lev3 - n_h.begin()] - n_h[i_lev3 - 1 - n_h.begin()] );
        dw1_3 = -1 * dw0_3;

        n3 = dw0_3*n[i_lev3 - 1 - n_h.begin()] + dw1_3*n[i_lev3 - n_h.begin()];


        k3h = dr*(u - k1u + 2*k2u);
        k3u = dr*( (1 - (u - k1u + 2*k2u)*(u - k1u + 2*k2u)) * ( n3 + 1 / (r - k1h + 2*k2h)));

        u_i = u;
        
        h = h + (k1h + 4*k2h + k3h)/6;
        u = u + (k1u + 4*k2u + k3u)/6;

        s = s + t*EARTH_R*asin(cos(asin(u_i))*dr / (EARTH_R + h));



        if(s >= dmax || s <= 0)
        {
            final_pos = {h, u, s};

            break;
        }
    }

    return final_pos;

}

void Tracer::backprop(const double h0, const double u0, const double dmax, const double dr, std::vector<double>& n, std::vector<double>& ndry, const double n_init, std::vector<double>& n_h,
                      const double lam0, const double mu0, const double lrate, int iter, std::vector<double>& m, std::vector<double>& v, int* index_n, double* dn_adj, double* obs_height)

// -----------------------------------------------------------------------------------------------------
//
// Adjoint-based gradient descent.  Analagous to backpropagation in machine learning, the gradient of 
// the loss function is propagated through the refractivity field.  Initally the loss function is 
// differentiatedwith respect to the state variables h and u, before the full d L / d n is evaluated
// through the backwards ray tracing.  The gradients d L / d n1 and d L / d n0 are evaluated at each step,
// incrementing the refractivity values immediately above and below the ray position.  A third order Runge-
// Kutta scheme is used to integrate the rays.
//
// mu = d L / d u (= 0)
//
// lam = d L / d h
//
// then d L / d n_x (eval at r_i) = mu * (1 - u*u) * dw_x/dh where w is the interpolation weight:
//
// n(h, N) = w0*N0 + w1*N1
//
// Inputs:
//  +  h0                 : initial height of ray (aircraft altitude)         
//  +  u0                 : initial direction of ray (evaluated by -dir end point)
//  +  dmax               : initial distance of ray (aircraft distance across surface)
//  +  dr                 : ray step size
//  +  (vector) n         : refractive index field to be optimised
//  +  n_surface_true     : fixed field at surface
//  +  (vector) n_h       : height of refractive index levels
//  +  lam0               : initial adjoint associated with h
//  +  mu0                : initial adjoint associated with u
//  +  lrate              : learning rate - controls who far along gradient direction to go (careful!)
//
// Some variables:
//  + dw0, etc            : the gradient of the weight used in the linear interpolation scheme. E.g. 
//                          n_eval = N0*(y1 - y)/(y1 - y0) + N1*(y - y0)/(y1 - y0), there
//                          dndh_eval = N0 * -1/(y1 - y0) + N1 * 1/(y1 - y0)
//  + k1,2,etc            : used for Runge-Kutta integration, y = y + (k1 + 4*k2 + k3)/6
//
//------------------------------------------------------------------------------------------------------

{
    steps = (int) 600/dr;

    h = h0;
    u = u0;
    s = dmax;

    mu = mu0;
    lam = lam0;

    if(obs_height != nullptr){
        h_obs = *obs_height;
    }
    else{
        h_obs = OBSERVER_H;
    }


    for(int i(0); i < steps; i++)
    { 
        r = h + EARTH_R;

        // k1 
        {
            i_lev1 = std::upper_bound(n_h.begin(), n_h.end(), h);

            dw0_1 = - 1  / ( n_h[i_lev1 - n_h.begin()] - n_h[i_lev1 - 1 - n_h.begin()] );
            dw1_1 = -1 * dw0_1;
                
            n1 = dw0_1*n[i_lev1 - 1 - n_h.begin()] + dw1_1*n[i_lev1 - n_h.begin()]; 
        }

        k1h = dr*u;
        k1u = dr*( (1 - u*u) * (n1 + (1 / r)) );

        // k2
        {
            i_lev2 = std::upper_bound(n_h.begin(), n_h.end(), (h + k1h/2));

            dw0_2 = - 1  / ( n_h[i_lev2 - n_h.begin()] - n_h[i_lev2 - 1 - n_h.begin()] );
            dw1_2 = -1 * dw0_2;

            n2 = dw0_2*n[i_lev2 - 1 - n_h.begin()] + dw1_2*n[i_lev2 - n_h.begin()]; 
        }

        k2h = dr*(u + k1u/2);
        k2u = dr*( (1 - (u + k1u/2)*(u + k1u/2) )  * ( n2 +  1 / (r + k1h/2)  ) );

        // k3
        {
            i_lev3 = std::upper_bound(n_h.begin(), n_h.end(), (h - k1h + 2*k2h));

            dw0_3 = - 1  / ( n_h[i_lev3 - n_h.begin()] - n_h[i_lev3 - 1 - n_h.begin()] );
            dw1_3 = -1 * dw0_3;

            n3 = dw0_3*n[i_lev3 - 1 - n_h.begin()] + dw1_3*n[i_lev3 - n_h.begin()];
        }

        k3h = dr*(u - k1u + 2*k2u);
        k3u = dr*( (1 - (u - k1u + 2*k2u)*(u - k1u + 2*k2u)) * ( n3 + 1 / (r - k1h + 2*k2h)));

        // k lam and mu


        {
            k1lam = dr*( mu * (1 - u*u) * ( n1*n1 + 1 / (r*r) ) );
            k1mu = dr*( -lam + 2 * mu * u * (n1 + (1 / r)));

            k2lam = dr*( (mu + k1mu/2) * (1 - (u + k1u/2)*(u + k1u/2) ) * ( n2*n2 + 1 / ( (r + k1h/2)*(r + k1h/2) ) ) );
            k2mu = dr*( -(lam + k1lam/2) + 2 * (mu + k1mu/2) * (u + k1u/2) * ( n2 + 1 / (r + k1h/2) ) );

            k3lam = dr*( (mu - k1mu + 2*k2mu) * (1 - (u - k1u + 2*k2u)*(u - k1u + 2*k2u) )  * ( n3*n3 + 1 / ( (r - k1h + 2*k2h)*(r - k1h + 2*k2h) ) ) );
            k3mu = dr*( -(lam - k1lam + 2*k2lam) + 2 * (mu - k1mu + 2*k2mu) * (u - k1u + 2*k2u) * (n3 + 1 / (r - k1h + 2*k2h) ) );
        }

        // Evaluate gradients for upper and lower refractivity values
        //-----------------------------------------------------------

        dn0 = mu*(u*u - 1)*dw0_1;
        dn1 = mu*(u*u - 1)*dw1_1;

        if (index_n != nullptr && dn_adj != nullptr)
        {

            if(n[i_lev1 - 1 - n_h.begin()] == n[*index_n])
            {
                *dn_adj += dn0*dr / exp(n[i_lev1 - 1 - n_h.begin()]);
            }

            if(n[i_lev1 - n_h.begin()] == n[*index_n])
            {
                *dn_adj += dn1*dr / exp(n[i_lev1 - n_h.begin()]);
            }

        }

        // Update respective refractivity values using gradient descent
        //-------------------------------------------------------------


        u_i = u;

        h = h + (k1h + 4*k2h + k3h)/6;
        u = u + (k1u + 4*k2u + k3u)/6;

        s = s - EARTH_R*asin(cos(asin(u_i))*dr / (EARTH_R + h)); // reverse ray tracing, so negative sign used

        mu = mu + (k1mu + 4*k2mu + k3mu)/6;
        lam = lam + (k1lam + 4*k2lam + k3lam)/6;

        if( s <= 1e-6 || h <= h_obs)
        {
            break;
        }

        if (index_n == nullptr && dn_adj == nullptr)
        {      

            //------------------
            //  Adam optimiser
            //------------------


            // m[i_lev1 - 1 - n_h.begin()] = beta1*m[i_lev1 - 1 - n_h.begin()] + (1.0-beta1)*dn0;
            // m[i_lev1 - n_h.begin()] = beta1*m[i_lev1 - n_h.begin()] + (1.0-beta1)*dn1;

            // v[i_lev1 - 1 - n_h.begin()] = beta2*v[i_lev1 - 1 - n_h.begin()] + (1.0-beta2)*(dn0*dn0);
            // v[i_lev1 - n_h.begin()] = beta2*v[i_lev1 - n_h.begin()] + (1.0-beta2)*(dn1*dn1);

            // m_est0 = m[i_lev1 - 1 - n_h.begin()] / (1.0 - pow(beta1, (iter+1)));
            // m_est1 = m[i_lev1 - n_h.begin()] / (1.0 - pow(beta1, (iter+1)));

            // v_est0 = v[i_lev1 - 1 - n_h.begin()] / (1.0 - pow(beta2, (iter+1)));
            // v_est1 = v[i_lev1 - n_h.begin()] / (1.0 - pow(beta2, (iter+1)));

            // n[i_lev1 - 1 - n_h.begin()] = n[i_lev1 - 1 - n_h.begin()] - lrate*m_est0 / sqrt(v_est0 + epsilon);
            // n[i_lev1 - n_h.begin()] = n[i_lev1 - n_h.begin()] - lrate*m_est1 / sqrt(v_est1 + epsilon);

            n[i_lev1 - 1 - n_h.begin()] = n[i_lev1 - 1 - n_h.begin()]*(1.0 - gamma) + gamma*std::max(n[i_lev1 - 1 - n_h.begin()] - lrate*dn0, 0.0);
            n[i_lev1 - n_h.begin()] = n[i_lev1 - n_h.begin()]*(1.0 - gamma) + gamma*std::max(n[i_lev1 - n_h.begin()] - lrate*dn1, 0.0);

            if( h > 0)
            {
                if (n[i_lev1 - 1 - n_h.begin()] < ndry[i_lev1 - 1 - n_h.begin()]) { n[i_lev1 - 1 - n_h.begin()] = ndry[i_lev1 - 1 - n_h.begin()]; }

                if (n[i_lev1 - n_h.begin()] < ndry[i_lev1 - n_h.begin()]) { n[i_lev1 - n_h.begin()] = ndry[i_lev1 - n_h.begin()]; }

                //n[i_lev1 - n_h.begin()] = ndry[i_lev1 - n_h.begin()];
                //n[i_lev1 - 1 - n_h.begin()] = ndry[i_lev1 - 1 - n_h.begin()];
            }
            
            n[0] = n_init; // Clamp surface refractivity

        }


        

        
    }

};

std::vector<std::vector<double>> Tracer::trace_paths(double h0, double u0, double dmax, double dr,
                                                 std::vector<double>& n, std::vector<double>& n_h, bool forward = true)
{
    steps = (int) 600/dr;


    h = h0;
    u = u0;

    if (forward)
    {
        s = 0.0;
        t = 1.0;
    }
    else
    {
        s = dmax;
        t = -1.0;
    }

    std::vector<double> vals;
    paths.clear();

    for(int i(0); i < steps; i++)
    { 

        vals = {s, h, u};
        paths.push_back(vals);
    
        r = h + EARTH_R;

        i_lev1 = std::upper_bound(n_h.begin(), n_h.end(), h);

        dw0_1 = - 1  / ( n_h[i_lev1 - n_h.begin()] - n_h[i_lev1 - 1 - n_h.begin()] );
        dw1_1 = -1 * dw0_1;
            
        n1 = dw0_1*n[i_lev1 - 1 - n_h.begin()] + dw1_1*n[i_lev1 - n_h.begin()]; 


        k1h = dr*u;
        k1u = dr*( (1 - u*u) * (n1 + (1 / r)) );

        i_lev2 = std::upper_bound(n_h.begin(), n_h.end(), (h + k1h/2));

        dw0_2 = - 1  / ( n_h[i_lev2 - n_h.begin()] - n_h[i_lev2 - 1 - n_h.begin()] );
        dw1_2 = -1 * dw0_2;

        n2 = dw0_2*n[i_lev2 - 1 - n_h.begin()] + dw1_2*n[i_lev2 - n_h.begin()]; 

        
        k2h = dr*(u + k1u/2);
        k2u = dr*( (1 - (u + k1u/2)*(u + k1u/2))  * ( n2 + 1 / (r + k1h/2)));

        i_lev3 = std::upper_bound(n_h.begin(), n_h.end(), (h - k1h + 2*k2h));

        dw0_3 = - 1  / ( n_h[i_lev3 - n_h.begin()] - n_h[i_lev3 - 1 - n_h.begin()] );
        dw1_3 = -1 * dw0_3;

        n3 = dw0_3*n[i_lev3 - 1 - n_h.begin()] + dw1_3*n[i_lev3 - n_h.begin()];


        k3h = dr*(u - k1u + 2*k2u);
        k3u = dr*( (1 - (u - k1u + 2*k2u)*(u - k1u + 2*k2u)) * ( n3 + 1 / (r - k1h + 2*k2h)));

        h = h + (k1h + 4*k2h + k3h)/6;
        u = u + (k1u + 4*k2u + k3u)/6;

        s = s + t*EARTH_R*asin(cos(asin(u_i))*dr / (EARTH_R + h));

        

        if(s > dmax || s < 0)
        {
            break;
        }
    }

    return paths;

}


#endif
