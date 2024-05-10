#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "inputs.h"
#include "constants.h"

Atmosphere::Atmosphere(const std::string& file)
    : file(file) { 
}

void Atmosphere::process() 
{

    std::ifstream inputFile(file);

    if (!inputFile.is_open())
    {
        std::cerr << "Error: Unable to open refractivity profile file " << file << std::endl;
        return;
    }
    else
    {
        std::cout << "Reading in atmosphere file: " << file << std::endl;
    }

    double h, n, ndry, nwet;

    while (inputFile >> h >> n >> ndry >> nwet)
    {   
        H_input.push_back(h*1e-3);
        N_input.push_back(n);
        NDRY_input.push_back(ndry);
        NWET_input.push_back(nwet);

    }

    inputFile.close();
}

const std::vector<double>& Atmosphere::H() const {
    return H_input;
}
const std::vector<double>& Atmosphere::N() const {
    return N_input;
}
const std::vector<double>& Atmosphere::NDRY() const {
    return NDRY_input;
}
const std::vector<double>& Atmosphere::NWET() const {
    return NWET_input;
}

ADSB::ADSB(const std::string& file, int length)
    : file(file), length(length) {
}

void ADSB::process() 
{

    std::ifstream inputFile(file);

    if (!inputFile.is_open())
    {
        std::cerr << "Error: Unable to open ADS-B data file " << file << std::endl;
        return;
    }
    else
    {
        std::cout << "Reading in: " << file << std::endl;
    }

    double obsaoa, h, d, repaoa, azi, t, arc, lat, lon;
    std::string icao;

    while (inputFile >> obsaoa >> h >> d >> repaoa >> azi >> t >> arc >> lat >> lon >> icao && (length == 0 || obsAOA_input.size() < length))
    {
        obsAOA_input.push_back(obsaoa);
        sin_obsAOA_input.push_back(sin(obsaoa*PI/180.0));
        rH_input.push_back(h);
        rD_input.push_back(d);
        repAOA_input.push_back(repaoa);
        azim_input.push_back(azi);
        time_input.push_back(t);
        arc_input.push_back(arc);
        lat_input.push_back(lat);

    }

    std::cout << "Number of observations: " << obsAOA_input.size() << std::endl;


    // while (inputFile >> obsaoa >> repaoa >> h >> d && (sin_obsAOA_input.size() < length))
    // {
    //     //obsAOA_input.push_back(obsaoa);
    //     sin_obsAOA_input.push_back(sin(obsaoa*pi/180.0));
    //     rH_input.push_back(h);
    //     rD_input.push_back(d);
    //     repAOA_input.push_back(repaoa);
    //     //azim_input.push_back(azi);
    //     //time_input.push_back(t);

    // }

    inputFile.close();
} 

const std::vector<double>& ADSB::obsAoA() const {
    return obsAOA_input;
}
const std::vector<double>& ADSB::sin_obsAoA() const {
    return sin_obsAOA_input;
}
const std::vector<double>& ADSB::rH() const {
    return rH_input;
}
const std::vector<double>& ADSB::rD() const {
    return rD_input;
}
const std::vector<double>& ADSB::repAoA() const {
    return repAOA_input;
}
const std::vector<double>& ADSB::azim() const {
    return azim_input;
}
const std::vector<double>& ADSB::rTime() const {
    return time_input;
}
const std::vector<double>& ADSB::arc() const {
    return arc_input;
}
const std::vector<double>& ADSB::lat() const {
    return lat_input;
}