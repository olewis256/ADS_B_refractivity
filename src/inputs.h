#ifndef INPUTS_H
#define INPUTS_H

#include <iostream>
#include <fstream>
#include <vector>

class Atmosphere 
{
    private:

        std::string file;
        std::vector<double> H_input;
        std::vector<double> N_input;
        std::vector<double> NDRY_input;
        std::vector<double> NWET_input;

    public:

        Atmosphere(const std::string& file);

        void process();

        const std::vector<double>& H() const;
        const std::vector<double>& N() const;
        const std::vector<double>& NDRY() const;
        const std::vector<double>& NWET() const;

};

class ADSB 
{
    private:

        int length;

        std::string file;
        std::vector<double> obsAOA_input;
        std::vector<double> sin_obsAOA_input;
        std::vector<double> rH_input;
        std::vector<double> rD_input;
        std::vector<double> repAOA_input;
        std::vector<double> azim_input;
        std::vector<double> time_input;
        std::vector<double> arc_input;

    public:

        ADSB(const std::string& file, int length = 0);

        void process();

        const std::vector<double>& obsAoA() const;
        const std::vector<double>& sin_obsAoA() const;
        const std::vector<double>& rH() const;
        const std::vector<double>& rD() const;
        const std::vector<double>& repAoA() const;
        const std::vector<double>& azim() const;
        const std::vector<double>& rTime() const;
        const std::vector<double>& arc() const;
};

#endif