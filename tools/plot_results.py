import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import smplotlib
import sys
import pandas as pd
import argparse
import os
import statistics as stat

r0 = 6391.579965656317

rep = lambda h, d : np.arctan( ((r0+h)*np.cos(d/r0) - (r0+0.560)) / ((r0+h)*np.sin(d/r0))) * 180/np.pi

parser = argparse.ArgumentParser(
                    prog='Refraction plotting',
                    description='Plotting tools ADS-B interferometer.\nIncludes plotting for real and synthetic data, as well as loss and apparent flightpath statistics.',
                    epilog='')

parser.add_argument('--index_profile', type=int,
                    help='Choose number for synthetic retrieved profile\n1: Sep \'23, 2: July \'22, 3: Dec \'22')

parser.add_argument('--difference', action='store_true',
                    help='Plot percentage difference between retrieved and synthetic profile')

parser.add_argument('--profile', action='store_true',
                    help='Plot retrieved profile')

parser.add_argument('--loss', action='store_true',
                    help='Plot loss and RMS statistics')

parser.add_argument('--true', action='store_true',
                    help='Plot true profile using raw ADS-B observations')

parser.add_argument('--paths', action='store_true',
                    help='Plot true and retrieved flightpaths')

parser.add_argument('--pdf', action='store_true',
                    help='Plot how errors depend on a variety of parameters')

parser.add_argument('--south', action='store_true',
                    help='Plot south side data')

parser.add_argument('--fill', action='store_true'
                    )

parser.add_argument('--heights', action='store_true',
                    help='Save the plot (jpeg, dpi=400)')

parser.add_argument('--error_func', action='store_true',
                    help='Save the plot (jpeg, dpi=400)')

parser.add_argument('--both', action='store_true',
                    help='Plot both profile and profile differences together')

parser.add_argument('--ray_path', action='store_true',
                    help='Plot both profile and profile differences together')

parser.add_argument('--sensitivity', action='store_true',
                    help='Plot both profile and profile differences together')

parser.add_argument('--n_err', type=int,
                    help='Magnitude of refractivity uncertainty in sensitivity')

parser.add_argument('--h_err', type=int,
                    help='Magnitude of altitude uncertainty in sensitivity')

parser.add_argument('--may', action='store_true',
                    help='Plot May data')     

parser.add_argument('--model', action='store_true',
                    help='Modelled refraction')      

parser.add_argument('--igarss', action='store_true',
                    help='Modelled refraction')      

parser.add_argument('--igarss_paths', action='store_true',
                    help='Modelled refraction')        


args = parser.parse_args()

data = [0]*3
datar = [0]*4
noises = [0, 0.01, 0.05]

if (args.index_profile in [1,2,3]):

    for i in range(3):
        
        data[i] = pd.read_csv("../retrievals/PAPERII_retrieve_NE_RK3_{}_{}.txt".format(args.index_profile, noises[i]), sep=' ', skiprows=1, header=None)
        data[i].columns = ["retrieve", "h", "target", "init"]

if (args.loss):
    for i in range(3):
        
        data[i] = pd.read_csv("../RMS/PAPERII_RMS_NE_RK3_{}.txt".format(noises[i]), sep=' ', header=None)
        data[i].columns = ["iter", "loss", "Nrms"]


if (args.true and args.profile):

    if(args.south):
        datar[0] = pd.read_csv("../retrievals/PAPERII_retrieve_oct_S_TRUE_t0_900.txt", sep=' ', header=None)
        datar[1] = pd.read_csv("../retrievals/PAPERII_retrieve_oct_S_TRUE_t900_1800.txt", sep=' ', header=None)
        datar[2] = pd.read_csv("../retrievals/PAPERII_retrieve_oct_S_TRUE_t1800_2700.txt", sep=' ', header=None)
        datar[3] = pd.read_csv("../retrievals/PAPERII_retrieve_oct_S_TRUE_t2700_3600.txt", sep=' ', header=None)

    else:
        datar[0] = pd.read_csv("../retrievals/PAPERII_retrieve_NE_TRUE_t0_900.txt", sep=' ', header=None)
        datar[1] = pd.read_csv("../retrievals/PAPERII_retrieve_NE_TRUE_t900_1800.txt", sep=' ', header=None)
        datar[2] = pd.read_csv("../retrievals/PAPERII_retrieve_NE_TRUE_t1800_2700.txt", sep=' ', header=None)
        datar[3] = pd.read_csv("../retrievals/PAPERII_retrieve_NE_TRUE_t2700_3600.txt", sep=' ', header=None)
        
    
    datar[0].columns = ["retrieve", "h", "target", "init", "ndry"]
    datar[1].columns = ["retrieve", "h", "target", "init", "ndry"]
    datar[2].columns = ["retrieve", "h", "target", "init", "ndry"]
    datar[3].columns = ["retrieve", "h", "target", "init", "ndry"]

    fig, ax1 = plt.subplots(figsize=(8, 8))

    ax1.plot(datar[0]['target'], datar[0]['h'], color='black', label='UKV')
    # ax1.plot(datar[0]['ndry'], datar[0]['h'], color='black', linestyle=':', label='UKV Ndry')
    ax1.plot(datar[0]['init'], datar[0]['h'], linewidth=0.5, linestyle='--', color='black', label='initial')
    ax1.plot(datar[0]['retrieve'], datar[0]['h'], color='green', label='0-15 mins')
    ax1.plot(datar[1]['retrieve'], datar[1]['h'], color='red', label='15-30 mins')
    ax1.plot(datar[2]['retrieve'], datar[2]['h'], color='blue', label='30-45 mins')
    ax1.plot(datar[3]['retrieve'], datar[3]['h'], color='deepskyblue', label='45-60 mins')

    plt.legend()

    ax1.set_ylabel("Altitude (km)")
    ax1.set_xlabel("Refractivity (ppm)")

    if(args.south):
        path = "../plots/PAPERII_retrieve_profiles_oct_S_true.jpeg"
    else:
        path = "../plots/PAPERII_retrieve_profiles_true.jpeg"

    if(os.path.exists(path)):
    
        ftrue = input("Plot already saved as {}, do you want to overwrite? [y/n]".format(path))

        if(ftrue == 'y'):
        
            os.remove(path)

            plt.savefig(path, dpi=400)

    else:

        print("Saving image file to: {}".format(path))

        plt.savefig(path, dpi=400)
    
    plt.show()

if (args.true and args.difference):


    if(args.south):
        datar[0] = pd.read_csv("../retrievals/PAPERII_retrieve_oct_S_TRUE_t0_900.txt", sep=' ', header=None)
        datar[1] = pd.read_csv("../retrievals/PAPERII_retrieve_oct_S_TRUE_t900_1800.txt", sep=' ', header=None)
        datar[2] = pd.read_csv("../retrievals/PAPERII_retrieve_oct_S_TRUE_t1800_2700.txt", sep=' ', header=None)
        datar[3] = pd.read_csv("../retrievals/PAPERII_retrieve_oct_S_TRUE_t2700_3600.txt", sep=' ', header=None)

    else:
        datar[0] = pd.read_csv("../retrievals/PAPERII_retrieve_NE_TRUE_t0_900.txt", sep=' ', header=None)
        datar[1] = pd.read_csv("../retrievals/PAPERII_retrieve_NE_TRUE_t900_1800.txt", sep=' ', header=None)
        datar[2] = pd.read_csv("../retrievals/PAPERII_retrieve_NE_TRUE_t1800_2700.txt", sep=' ', header=None)
        datar[3] = pd.read_csv("../retrievals/PAPERII_retrieve_NE_TRUE_t2700_3600.txt", sep=' ', header=None)
    
    datar[0].columns = ["retrieve", "h", "target", "init", "ndry"]
    datar[1].columns = ["retrieve", "h", "target", "init", "ndry"]
    datar[2].columns = ["retrieve", "h", "target", "init", "ndry"]
    datar[3].columns = ["retrieve", "h", "target", "init", "ndry"]

    fig, ax1 = plt.subplots(figsize=(6, 8))


    ax1.plot(100*(datar[0]['init'] - datar[0]['target'])/datar[0]['target'], datar[0]['h'], linewidth=0.5, linestyle='--', color='black', label='initial')
    ax1.plot(100*(datar[0]['retrieve'] - datar[0]['target'])/datar[0]['target'], datar[0]['h'], color='green', label='0-15 mins')
    ax1.plot(100*(datar[1]['retrieve'] - datar[0]['target'])/datar[0]['target'], datar[0]['h'], color='red', label='15-30 mins')
    ax1.plot(100*(datar[2]['retrieve'] - datar[0]['target'])/datar[0]['target'], datar[0]['h'], color='blue', label='30-45 mins')
    ax1.plot(100*(datar[2]['retrieve'] - datar[0]['target'])/datar[0]['target'], datar[0]['h'], color='deepskyblue', label='45-60 mins')

    plt.axvline(0.0, linestyle=':', color='black', label='UKV')

    ax1.set_xlim(-25, 25)
    ax1.set_ylabel("Altitude (km)")
    ax1.set_xlabel(r"$\frac{N_{retrieved} - N_{UKV}}{N_{UKV}}$"+" (%)")

    plt.legend()

    fig.tight_layout()

    path = "../plots/PAPERII_retrieve_difference_real.jpeg"

    if(os.path.exists(path)):
    
        ftrue = input("Plot already saved as {}, do you want to overwrite? [y/n]".format(path))

        if(ftrue == 'y'):
        
            os.remove(path)

            plt.savefig(path, dpi=400)

    else:

        print("Saving image file to: {}".format(path))

        plt.savefig(path, dpi=400)

    plt.show()



if (args.paths) and not (args.may):

    if(args.south):
        datar[0] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t0_900.txt", sep=' ', header=None)
        datar[1] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t900_1800.txt", sep=' ', header=None)
        datar[2] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t1800_2700.txt", sep=' ', header=None)
        datar[3] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t2700_3600.txt", sep=' ', header=None)
    else:
        datar[0] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t0_900.txt", sep=' ', header=None)
        datar[1] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t900_1800.txt", sep=' ', header=None)
        datar[2] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t1800_2700.txt", sep=' ', header=None)
        datar[3] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t2700_3600.txt", sep=' ', header=None)

    datar[0].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    datar[1].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    datar[2].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    datar[3].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]

    obsAoA1, obsAoA2, obsAoA3, obsAoA4 = datar[0]['obsAoA'], datar[1]['obsAoA'], datar[2]['obsAoA'], datar[3]['obsAoA']
    repAoA1, repAoA2, repAoA3, repAoA4 = datar[0]['repAoA'], datar[1]['repAoA'], datar[2]['repAoA'], datar[3]['repAoA']
    h1, h2, h3, h4 = datar[0]['h'], datar[1]['h'], datar[2]['h'], datar[3]['h']
    d1, d2, d3, d4 = datar[0]['d'], datar[1]['d'], datar[2]['d'], datar[3]['d']
    hinit1, hinit2, hinit3, hinit4  = datar[0]['hinit'], datar[1]['hinit'], datar[2]['hinit'], datar[3]['hinit']
    hret1, hret2, hret3, hret4 = datar[0]['hret'], datar[1]['hret'], datar[2]['hret'], datar[3]['hret']
    azim1, azim2, azim3, azim4 = datar[0]['azim'], datar[1]['azim'], datar[2]['azim'], datar[3]['azim']

    repAoA_sph1, repAoA_sph2, repAoA_sph3, repAoA_sph4 = rep(h1, d1), rep(h2, d2), rep(h3, d3), rep(h4, d4)
    repAoA_init1, repAoA_init2, repAoA_init3, repAoA_init4 = rep(hinit1, d1), rep(hinit2, d2), rep(hinit3, d3), rep(hinit4, d4)
    repAoA_ret1, repAoA_ret2, repAoA_ret3, repAoA_ret4 = rep(hret1, d1), rep(hret2, d2), rep(hret3, d3), rep(hret4, d4)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2,nrows=2,sharey=True,sharex=True,figsize=(10,10))

    ax1.scatter(-azim1, repAoA1, color='black', label='obs', s=1.5)
    ax1.scatter(-azim1, repAoA_init1, linewidth=0.5, color='blue', label='initial',s=1.5)
    ax1.scatter(-azim1, repAoA_ret1, color='green', label='retrieved',s=1.5)
    # ax1.scatter(-azim1, repAoA_sph1, color='skyblue', label='retrieved',s=1.5)
    ax1.set_ylabel("Reported AoA (deg.)")
    ax1.set_xlabel("(a)")

    print("Maximum difference between ellipsoid and spheroid repAoA for plot 1:", max(repAoA1 - repAoA_sph1))

    ax1.set_title('$0 \leq t < 15 $ (mins.)')

    ax2.scatter(-azim2, repAoA2, color='black',s=1.5)
    ax2.scatter(-azim2, repAoA_init2, linewidth=0.5, color='blue',s=1.5)
    ax2.scatter(-azim2, repAoA_ret2, color='green',s=1.5)
    # ax2.scatter(-azim2, repAoA_sph2, color='skyblue', label='retrieved',s=1.5)

    ax2.set_xlabel("(b)")

    ax2.sharey(ax1)

    ax2.set_title('$15 \leq t < 30$ (mins.)')

    ax3.scatter(-azim3, repAoA3, color='black', s=1.5)
    ax3.scatter(-azim3, repAoA_init3, linewidth=0.5, color='blue',s=1.5)
    ax3.scatter(-azim3, repAoA_ret3, color='green',s=1.5)
    # ax3.scatter(-azim3, repAoA_sph3, color='skyblue', label='retrieved',s=1.5)

    ax3.set_xlabel("Horizontal angle (deg.) \n (c)")
    ax4.set_xlabel("Horizontal angle (deg.) \n (d)")
    ax3.set_ylabel("Reported AoA (deg.)")

    ax1.sharex(ax3)

    ax3.set_title('$30 \leq t < 45$ (mins.)')

    ax4.scatter(-azim4, repAoA4, color='black', s=1.5)
    ax4.scatter(-azim4, repAoA_init4, linewidth=0.5, color='blue',s=1.5)
    ax4.scatter(-azim4, repAoA_ret4, color='green',s=1.5)
    # ax4.scatter(-azim4, repAoA_sph4, color='skyblue', label='retrieved',s=1.5)

    ax4.set_title('$45 \leq t < 60$ (mins.)')

    fig.legend(bbox_to_anchor=(0.3, 0.68), markerscale=2)
    fig.tight_layout()

    if(args.south):
        path = "../plots/PAPERII_retrieve_oct_S_flightpaths.jpeg"
    else:
        path = "../plots/PAPERII_retrieve_flightpaths.jpeg"

    if(os.path.exists(path)):
    
        ftrue = input("Plot already saved as {}, do you want to overwrite? [y/n]".format(path))

        if(ftrue == 'y'):
        
            os.remove(path)

            fig.savefig(path, dpi=400)

    else:

        print("Saving image file to: {}".format(path))

        fig.savefig(path, dpi=400)
    
    plt.show()

if (args.true and args.loss):

    if(args.south):
        datar[0] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t0_900.txt", sep=' ', header=None)
        datar[1] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t900_1800.txt", sep=' ', header=None)
        datar[2] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t1800_2700.txt", sep=' ', header=None)
        datar[3] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t2700_3600.txt", sep=' ', header=None)
    else:
        datar[0] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t0_900.txt", sep=' ', header=None)
        datar[1] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t900_1800.txt", sep=' ', header=None)
        datar[2] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t1800_2700.txt", sep=' ', header=None)
        datar[3] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t2700_3600.txt", sep=' ', header=None)

    datar[0].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    datar[1].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    datar[2].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    datar[3].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]

    obsAoA1, obsAoA2, obsAoA3, obsAoA4 = datar[0]['obsAoA'], datar[1]['obsAoA'], datar[2]['obsAoA'], datar[3]['obsAoA']
    repAoA1, repAoA2, repAoA3, repAoA4 = datar[0]['repAoA'], datar[1]['repAoA'], datar[2]['repAoA'], datar[3]['repAoA']
    h1, h2, h3, h4 = datar[0]['h'], datar[1]['h'], datar[2]['h'], datar[3]['h']
    d1, d2, d3, d4 = datar[0]['d'], datar[1]['d'], datar[2]['d'], datar[3]['d']
    hinit1, hinit2, hinit3, hinit4  = datar[0]['hinit'], datar[1]['hinit'], datar[2]['hinit'], datar[3]['hinit']
    hret1, hret2, hret3, hret4 = datar[0]['hret'], datar[1]['hret'], datar[2]['hret'], datar[3]['hret']
    azim1, azim2, azim3, azim4 = datar[0]['azim'], datar[1]['azim'], datar[2]['azim'], datar[3]['azim']

    repAoA1_sph, repAoA2_sph, repAoA3_sph, repAoA4_sph = rep(h1, d1), rep(h2, d2), rep(h3, d3), rep(h4, d4)

    error1, error1init = rep(hret1, d1) - repAoA1, rep(hinit1, d1) - repAoA1
    error2, error2init = rep(hret2, d2) - repAoA2, rep(hinit2, d2) - repAoA2
    error3, error3init = rep(hret3, d3) - repAoA3, rep(hinit3, d3) - repAoA3
    error4, error4init = rep(hret4, d4) - repAoA4, rep(hinit4, d4) - repAoA4

    n = 100

    vals = np.linspace(0, 450, n)
    valsp = np.linspace(0, 450, n)
    mean_err1, mean_err2, mean_err3, mean_err4 = np.zeros(n),np.zeros(n),np.zeros(n),np.zeros(n)
    mean_err_init1, mean_err_init2, mean_err_init3, mean_err_init4 = np.zeros(n),np.zeros(n),np.zeros(n),np.zeros(n)

    for i in range(len(vals)-1):
        
        mask1 = (d1 >= vals[i]) & (d1 < vals[i+1])
        mask2 = (d2 >= vals[i]) & (d2 < vals[i+1])
        mask3 = (d3 >= vals[i]) & (d3 < vals[i+1])
        mask4 = (d4 >= vals[i]) & (d4 < vals[i+1])



        if (~np.isnan(np.mean(error1[mask1])) and ~np.isnan(np.mean(error1init[mask1]))):
            mean_err1[i] = np.mean(error1[mask1])
            mean_err_init1[i] = np.mean(error1init[mask1])
        if (~np.isnan(np.mean(error2[mask2])) and ~np.isnan(np.mean(error2init[mask2]))):
            mean_err2[i] = np.mean(error2[mask2])
            mean_err_init2[i] = np.mean(error2init[mask2])
        if (~np.isnan(np.mean(error3[mask3])) and ~np.isnan(np.mean(error3init[mask3]))):
            mean_err3[i] = np.mean(error3[mask3])
            mean_err_init3[i] = np.mean(error3init[mask3])
        if (~np.isnan(np.mean(error4[mask4])) and ~np.isnan(np.mean(error4init[mask4]))):
            mean_err4[i] = np.mean(error4[mask4])
            mean_err_init4[i] = np.mean(error4init[mask4])
           

    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2,nrows=2,sharey=True,sharex=True,figsize=(10,10))

    
    ax1.plot(valsp, mean_err1, color='green', label='retrieved')
    ax1.plot(valsp, mean_err_init1, color='black', label='initial')
    ax1.axhline(0.0, color='black')
   
    ax1.set_ylabel("Reported AoA error (deg.)")#\n(Aircraft altitude ($10^2$ km))")
    ax1.set_xlabel("(a)")

    ax1.set_title('$0 < t \leq 15 $ (mins.)')

    ax2.plot(valsp, mean_err2, color='green')
    ax2.plot(valsp, mean_err_init2, color='black')
    ax2.axhline(0.0, color='black')
    ax2.set_xlabel("(b)")

    ax2.set_title('$15 < t \leq 30$ (mins.)')

    ax3.plot(valsp, mean_err3, color='green')
    ax3.plot(valsp, mean_err_init3, color='black')
    ax3.axhline(0.0, color='black')

    ax3.set_xlabel("Distance (km) \n (c)")
    ax4.set_xlabel("Distance (km) \n (d)")
    ax3.set_ylabel("Reported AoA error (deg.)")

    ax1.sharex(ax3)

    ax3.set_title('$30 < t \leq 45$ (mins.)')

    ax4.plot(valsp, mean_err4, color='green')
    ax4.plot(valsp, mean_err_init4, color='black')
    ax4.axhline(0.0, color='black')
    ax4.set_title('$45 < t\leq 60$ (mins.)')

    fig.legend(bbox_to_anchor=(0.49, 0.9), markerscale=2)
    fig.tight_layout()

    if(args.south):
        path = "../plots/PAPERII_retrieve_flightpaths_oct_S_error.jpeg"
    else:
        path = "../plots/PAPERII_retrieve_flightpaths_error.jpeg"

    if(os.path.exists(path)):
    
        ftrue = input("Plot already saved as {}, do you want to overwrite? [y/n]".format(path))

        if(ftrue == 'y'):
        
            os.remove(path)

            fig.savefig(path, dpi=400)

    else:

        print("Saving image file to: {}".format(path))

        fig.savefig(path, dpi=400)

    plt.show()

if (args.true and args.pdf):

    if(args.south):
        datar[0] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t0_900.txt", sep=' ', header=None)
        datar[1] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t900_1800.txt", sep=' ', header=None)
        datar[2] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t1800_2700.txt", sep=' ', header=None)
        datar[3] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t2700_3600.txt", sep=' ', header=None)
    else:
        datar[0] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t0_900.txt", sep=' ', header=None)
        datar[1] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t900_1800.txt", sep=' ', header=None)
        datar[2] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t1800_2700.txt", sep=' ', header=None)
        datar[3] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t2700_3600.txt", sep=' ', header=None)

    datar[0].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    datar[1].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    datar[2].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    datar[3].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]

    obsAoA1, obsAoA2, obsAoA3, obsAoA4 = datar[0]['obsAoA'], datar[1]['obsAoA'], datar[2]['obsAoA'], datar[3]['obsAoA']
    repAoA1, repAoA2, repAoA3, repAoA4 = datar[0]['repAoA'], datar[1]['repAoA'], datar[2]['repAoA'], datar[3]['repAoA']
    h1, h2, h3, h4 = datar[0]['h'], datar[1]['h'], datar[2]['h'], datar[3]['h']
    d1, d2, d3, d4 = datar[0]['d'], datar[1]['d'], datar[2]['d'], datar[3]['d']
    hinit1, hinit2, hinit3, hinit4  = datar[0]['hinit'], datar[1]['hinit'], datar[2]['hinit'], datar[3]['hinit']
    hret1, hret2, hret3, hret4 = datar[0]['hret'], datar[1]['hret'], datar[2]['hret'], datar[3]['hret']
    azim1, azim2, azim3, azim4 = datar[0]['azim'], datar[1]['azim'], datar[2]['azim'], datar[3]['azim']

    repAoA1_sph, repAoA2_sph, repAoA3_sph, repAoA4_sph = rep(h1, d1), rep(h2, d2), rep(h3, d3), rep(h4, d4)

    error1, error1init = rep(hret1, d1) - repAoA1, rep(hinit1, d1) - repAoA1
    error2, error2init = rep(hret2, d2) - repAoA2, rep(hinit2, d2) - repAoA2
    error3, error3init = rep(hret3, d3) - repAoA3, rep(hinit3, d3) - repAoA3
    error4, error4init = rep(hret4, d4) - repAoA4, rep(hinit4, d4) - repAoA4

    error = np.concatenate([error1, error2, error3, error4])
    error_init = np.concatenate([error1init, error2init, error3init, error4init])
    print(stat.stdev(error_init), stat.stdev(error))
    print(stat.mean(error_init), stat.mean(error))
    

    
    
    fig, ax1 = plt.subplots(figsize=(8, 6))


    ax1.hist(error, bins=100, 
                   histtype=u'step', density=True, color='green', label='retrieved') 
    ax1.hist(error_init, bins=100, 
                   histtype=u'step', density=True, color='blue', label='initial') 

    ax1.set_xlim(-max(error_init), max(error_init))
    handle1 = Line2D([], [], c='r')
    handle2 = Line2D([], [], c='b')
    handles, labels = ax1.get_legend_handles_labels()
    new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]

    plt.legend(handles=new_handles, labels=labels)

    plt.axvline(0.0, linestyle='--', color='black')

    ax1.set_xlabel("Reported AoA error (deg.)")
    ax1.set_ylabel("Number of broadcasts")
    
    if(args.south):
        path = "../plots/PAPERII_retrieve_oct_S_STDEV.jpeg"
    else:
        path = "../plots/PAPERII_retrieve_STDEV.jpeg"

    if(os.path.exists(path)):
    
        ftrue = input("Plot already saved as {}, do you want to overwrite? [y/n]".format(path))

        if(ftrue == 'y'):
        
            os.remove(path)

            plt.savefig(path, dpi=400)

    else:

        print("Saving image file to: {}".format(path))

        plt.savefig(path, dpi=400)



    
    plt.show()

if (args.fill):

    upper, lower = np.zeros((4, 30)), np.zeros((4, 30))

    data = [[0]*3]*4

    for j in range(4):
        for i in range(3):

            data[j][i] = pd.read_csv("../retrievals/PAPERII_retrieve_NE_TRUE_t{}_{}_{}.txt".format(j*900, (j+1)*900, i), sep=' ', header=None)
            
            data[j][i].columns = ["retrieve", "h", "target", "init", "ndry"]

        for k in range(len(data[j][i])):

            upper[j][k] = max([data[j][0]["retrieve"][k], data[j][1]["retrieve"][k], data[j][2]["retrieve"][k]])
            lower[j][k] = min([data[j][0]["retrieve"][k], data[j][1]["retrieve"][k], data[j][2]["retrieve"][k]])

            print([data[j][0]["retrieve"][k], data[j][1]["retrieve"][k], data[j][2]["retrieve"][k]])

    fig, ax1 = plt.subplots(figsize=(8, 6))
    print(upper, lower)
    ax1.plot(data[0][0]['target'], data[0][0]['h'], color='black', label='UKV')
    ax1.plot(data[0][0]['init'], data[0][0]['h'], linewidth=0.5, linestyle='--', color='black', label='initial')
    ax1.fill_betweenx(data[0][0]['h'], lower[0], upper[0])
    ax1.fill_betweenx(data[1][0]['h'], lower[1], upper[1])
    ax1.fill_betweenx(data[2][0]['h'], lower[2], upper[2])
    ax1.fill_betweenx(data[3][0]['h'], lower[3], upper[3])

    plt.legend()

    ax1.set_ylabel("Altitude (km)")
    ax1.set_xlabel("Refractivity (ppm)")

    if(args.south):
        path = "../plots/PAPERII_retrieve_profiles_oct_S_true_fill.jpeg"
    else:
        path = "../plots/PAPERII_retrieve_profiles_true_fill.jpeg"

    if(os.path.exists(path)):
    
        ftrue = input("Plot already saved as {}, do you want to overwrite? [y/n]".format(path))

        if(ftrue == 'y'):
        
            os.remove(path)

            plt.savefig(path, dpi=400)

    else:

        print("Saving image file to: {}".format(path))

        plt.savefig(path, dpi=400)



    
    plt.show()



    # vals = np.linspace(0, 12, n)
    # valsp = np.linspace(0, 12, n)
    # mean_err1, mean_err2, mean_err3, mean_err4 = np.zeros(n),np.zeros(n),np.zeros(n),np.zeros(n)
    # mean_err_init1, mean_err_init2, mean_err_init3, mean_err_init4 = np.zeros(n),np.zeros(n),np.zeros(n),np.zeros(n)

    # for i in range(len(vals)-1):
        
    #     mask1 = (h1 >= vals[i]) & (h1 < vals[i+1])
    #     mask2 = (h2 >= vals[i]) & (h2 < vals[i+1])
    #     mask3 = (h3 >= vals[i]) & (h3 < vals[i+1])
    #     mask4 = (h4 >= vals[i]) & (h4 < vals[i+1])



    #     if (~np.isnan(np.mean(error1[mask1])) and ~np.isnan(np.mean(error1init[mask1]))):
    #         mean_err1[i] = np.mean(error1[mask1])
    #         mean_err_init1[i] = np.mean(error1init[mask1])
    #     if (~np.isnan(np.mean(error2[mask2])) and ~np.isnan(np.mean(error2init[mask2]))):
    #         mean_err2[i] = np.mean(error2[mask2])
    #         mean_err_init2[i] = np.mean(error2init[mask2])
    #     if (~np.isnan(np.mean(error3[mask3])) and ~np.isnan(np.mean(error3init[mask3]))):
    #         mean_err3[i] = np.mean(error3[mask3])
    #         mean_err_init3[i] = np.mean(error3init[mask3])
    #     if (~np.isnan(np.mean(error4[mask4])) and ~np.isnan(np.mean(error4init[mask4]))):
    #         mean_err4[i] = np.mean(error4[mask4])
    #         mean_err_init4[i] = np.mean(error4init[mask4])
           

    
    # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2,nrows=2,sharey=True,sharex=True,figsize=(10,10))

    
    # ax1.plot(mean_err1, valsp, color='green', label='retrieved')
    # ax1.plot(mean_err_init1, valsp, color='black', label='initial')
    # ax1.axvline(0.0, color='black')
   
    # ax1.set_ylabel("Reported AoA error (deg.)")#\n(Aircraft altitude ($10^2$ km))")
    # ax1.set_xlabel("(a)")

    # ax1.set_title('$0 < t \leq 15 $ (mins.)')

    # ax2.plot(mean_err2, valsp, color='green')
    # ax2.plot(mean_err_init2, valsp, color='black')
    # ax2.axvline(0.0, color='black')
    # ax2.set_xlabel("(b)")

    # ax2.set_title('$15 < t \leq 30$ (mins.)')

    # ax3.plot(mean_err3, valsp, color='green')
    # ax3.plot(mean_err_init3, valsp, color='black')
    # ax3.axvline(0.0, color='black')

    # ax3.set_xlabel("Distance (km) \n (c)")
    # ax4.set_xlabel("Distance (km) \n (d)")
    # ax3.set_ylabel("Reported AoA error (deg.)")

    # ax1.sharex(ax3)

    # ax3.set_title('$30 < t \leq 45$ (mins.)')

    # ax4.plot(mean_err4, valsp, color='green')
    # ax4.plot(mean_err_init4, valsp, color='black')
    # ax4.axvline(0.0, color='black')
    # ax4.set_title('$45 < t\leq 60$ (mins.)')

    # fig.legend(bbox_to_anchor=(0.49, 0.9), markerscale=2)
    # fig.tight_layout()

    # path = "../plots/PAPERII_retrieve_flightpaths_error_AOA.jpeg"

    # if(os.path.exists(path)):
    
    #     ftrue = input("Plot already saved as {}, do you want to overwrite? [y/n]".format(path))

    #     if(ftrue == 'y'):
        
    #         os.remove(path)

    #         fig.savefig(path, dpi=400)

    # else:

    #     print("Saving image file to: {}".format(path))

    #     fig.savefig(path, dpi=400)

    # plt.show()

if (args.index_profile in [1,2,3] and args.profile):

    fig, ax1 = plt.subplots(figsize=(8, 8))


    ax1.plot(data[0]['target'], data[0]['h'], color='black', label='radiosonde')
    ax1.plot(data[0]['init'], data[0]['h'], linewidth=0.5, linestyle='--', color='black', label='initial')
    ax1.plot(data[0]['retrieve'], data[0]['h'], color='green', label='0.00 RMS')
    ax1.plot(data[1]['retrieve'], data[0]['h'], color='blue', label='{} RMS'.format(noises[1]))
    ax1.plot(data[2]['retrieve'], data[0]['h'], color='red', label='{} RMS'.format(noises[2]))

    ax1.set_ylabel("Altitude (km)")
    ax1.set_xlabel("Refractivity (ppm)")

    fig.legend()

    fig.tight_layout()
    num = len(data[0])
    RMS1, RMS2, RMS3 = 0,0,0
    for i in range(num):
        RMS1 += (data[0]['retrieve'][i] - data[0]['target'][i])**2
        RMS2 += (data[1]['retrieve'][i] - data[0]['target'][i])**2
        RMS3 += (data[2]['retrieve'][i] - data[0]['target'][i])**2

    print(np.sqrt(RMS1 / num), np.sqrt(RMS2 / num), np.sqrt(RMS3 / num))

    path = "../plots/PAPERII_retrieve_profiles_synthetic_{}.jpeg".format(args.index_profile)

    if(os.path.exists(path)):
    
        ftrue = input("Plot already saved as {}, do you want to overwrite? [y/n]".format(path))

        if(ftrue == 'y'):
        
            os.remove(path)

            plt.savefig(path, dpi=400)

    else:

        print("Saving image file to: {}".format(path))

        plt.savefig(path, dpi=400)



    
    plt.show()

if (args.index_profile in [1,2,3] and args.difference):

    fig, ax1 = plt.subplots(figsize=(6, 8))


    ax1.plot(100*(data[0]['init'] - data[0]['target'])/data[0]['target'], data[0]['h'], color='black', linestyle='--', label='initial')
    ax1.plot(100*(data[0]['retrieve'] - data[0]['target'])/data[0]['target'], data[0]['h'], color='green', label='0.00 RMS')
    ax1.plot(100*(data[1]['retrieve'] - data[0]['target'])/data[0]['target'], data[0]['h'], color='blue', label='{} RMS'.format(noises[1]))
    ax1.plot(100*(data[2]['retrieve'] - data[0]['target'])/data[0]['target'], data[0]['h'], color='red', label='{} RMS'.format(noises[2]))

    plt.axvline(0.0, linestyle=':', color='black')

    ax1.set_xlim(-15, 15)
    ax1.set_ylabel("Altitude (km)")
    ax1.set_xlabel(r'$\frac{N_{retrieved} - N_{radiosonde}}{N_{radiosonde}}$'+' (%)')

    plt.legend()

    fig.tight_layout()

    path = "../plots/PAPERII_retrieve_difference_synthetic_{}.jpeg".format(args.index_profile)

    if(os.path.exists(path)):
    
        ftrue = input("Plot already saved as {}, do you want to overwrite? [y/n]".format(path))

        if(ftrue == 'y'):
        
            os.remove(path)

            plt.savefig(path, dpi=400)

    else:

        print("Saving image file to: {}".format(path))

        plt.savefig(path, dpi=400)



    
    plt.show()

if (args.index_profile in [1,2,3] and args.loss):

    fig, (ax1, ax2) = plt.subplots(ncols=2,figsize=(10, 6))

    ax1.plot(data[0]['iter'], data[0]['loss'], color='green', label='0.00 RMS')
    ax1.plot(data[0]['iter'], data[1]['loss'], color='blue', label='{} RMS')
    ax1.plot(data[0]['iter'], data[2]['loss'], color='red', label='{} RMS')

    ax1.set_ylabel("Loss")
    ax1.set_xlabel("Iteration number")
    ax1.set_yscale('log')

    

    ax2.plot(data[0]['iter'], data[0]['Nrms'], color='green', label='0.00 RMS')
    ax2.plot(data[0]['iter'], data[1]['Nrms'], color='blue', label='{} RMS'.format(noises[1]))
    ax2.plot(data[0]['iter'], data[2]['Nrms'], color='red', label='{} RMS'.format(noises[2]))

    ax2.set_ylabel("Refractivity RMS (ppm)")
    ax2.set_xlabel("Iteration number")
    
    plt.legend()
    fig.tight_layout()

    fig.savefig("../Plots/PAPERII_retrieve_loss.jpeg", dpi=700)

    plt.show()

if (args.heights):


    if(args.south):
        datar[0] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t0_900.txt", sep=' ', header=None)
        datar[1] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t900_1800.txt", sep=' ', header=None)
        datar[2] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t1800_2700.txt", sep=' ', header=None)
        datar[3] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t2700_3600.txt", sep=' ', header=None)
    else:
        datar[0] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t0_900.txt", sep=' ', header=None)
        datar[1] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t900_1800.txt", sep=' ', header=None)
        datar[2] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t1800_2700.txt", sep=' ', header=None)
        datar[3] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t2700_3600.txt", sep=' ', header=None)

    datar[0].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    datar[1].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    datar[2].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    datar[3].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]

    fig, ax1 = plt.subplots(figsize=(8, 6))

    ax1.hist(datar[0]['h'], bins=20, histtype=u'step', 
             color='green', label='retrieved',orientation='horizontal') 
    ax1.hist(datar[1]['h'], bins=20, histtype=u'step', 
             color='red', label='initial',orientation='horizontal') 
    ax1.hist(datar[2]['h'], bins=20, histtype=u'step', 
             color='blue', label='initial',orientation='horizontal') 
    ax1.hist(datar[3]['h'], bins=20, histtype=u'step',
              color='deepskyblue', label='initial',orientation='horizontal') 


    plt.legend()

    ax1.set_ylabel("Altitude (km)")
    ax1.set_xlabel("Refractivity (ppm)")

    if(args.south):
        path = "../plots/PAPERII_heights_oct_S.jpeg"
    else:
        path = "../plots/PAPERII_heights.jpeg"

    if(os.path.exists(path)):
    
        ftrue = input("Plot already saved as {}, do you want to overwrite? [y/n]".format(path))

        if(ftrue == 'y'):
        
            os.remove(path)

            plt.savefig(path, dpi=400)

    else:

        print("Saving image file to: {}".format(path))

        plt.savefig(path, dpi=400)



    
    plt.show()

if (args.error_func):

    if(args.south):
        datar[0] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t0_900.txt", sep=' ', header=None)
        datar[1] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t900_1800.txt", sep=' ', header=None)
        datar[2] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t1800_2700.txt", sep=' ', header=None)
        datar[3] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t2700_3600.txt", sep=' ', header=None)
    else:
        datar[0] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t0_900.txt", sep=' ', header=None)
        datar[1] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t900_1800.txt", sep=' ', header=None)
        datar[2] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t1800_2700.txt", sep=' ', header=None)
        datar[3] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t2700_3600.txt", sep=' ', header=None)

    datar[0].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    datar[1].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    datar[2].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    datar[3].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]

    obsAoA1, obsAoA2, obsAoA3, obsAoA4 = datar[0]['obsAoA'], datar[1]['obsAoA'], datar[2]['obsAoA'], datar[3]['obsAoA']
    repAoA1, repAoA2, repAoA3, repAoA4 = datar[0]['repAoA'], datar[1]['repAoA'], datar[2]['repAoA'], datar[3]['repAoA']
    h1, h2, h3, h4 = datar[0]['h'], datar[1]['h'], datar[2]['h'], datar[3]['h']
    d1, d2, d3, d4 = datar[0]['d'], datar[1]['d'], datar[2]['d'], datar[3]['d']
    hinit1, hinit2, hinit3, hinit4  = datar[0]['hinit'], datar[1]['hinit'], datar[2]['hinit'], datar[3]['hinit']
    hret1, hret2, hret3, hret4 = datar[0]['hret'], datar[1]['hret'], datar[2]['hret'], datar[3]['hret']
    azim1, azim2, azim3, azim4 = datar[0]['azim'], datar[1]['azim'], datar[2]['azim'], datar[3]['azim']

    repAoA1_sph, repAoA2_sph, repAoA3_sph, repAoA4_sph = rep(h1, d1), rep(h2, d2), rep(h3, d3), rep(h4, d4)

    error1, error1init = rep(hret1, d1) - repAoA1, rep(hinit1, d1) - repAoA1
    error2, error2init = rep(hret2, d2) - repAoA2, rep(hinit2, d2) - repAoA2
    error3, error3init = rep(hret3, d3) - repAoA3, rep(hinit3, d3) - repAoA3
    error4, error4init = rep(hret4, d4) - repAoA4, rep(hinit4, d4) - repAoA4

    error = np.concatenate([error1, error2, error3, error4])
    error_init = np.concatenate([error1init, error2init, error3init, error4init])
    print(stat.stdev(error_init), stat.stdev(error))
    print(stat.mean(error_init), stat.mean(error))
    

    
    
    fig, ax1 = plt.subplots(figsize=(8, 6))


    plot = ax1.scatter(d1, h1, c=error1, s=50, edgecolors="none")

   
    
    ax1.set_xlabel("Distance (km)")
    ax1.set_ylabel("Altitude (km)")

    cbar = fig.colorbar(plot, ax=ax1)
    cbar.set_label('Retrieved - true rep. AoA (deg.)',labelpad=10)
    
    if(args.south):
        path = "../plots/PAPERII_retrieve_oct_S_STDEV_dist.jpeg"
    else:
        path = "../plots/PAPERII_retrieve_STDEV_dist1.jpeg"

    if(os.path.exists(path)):
    
        ftrue = input("Plot already saved as {}, do you want to overwrite? [y/n]".format(path))

        if(ftrue == 'y'):
        
            os.remove(path)

            plt.savefig(path, dpi=400)

    else:

        print("Saving image file to: {}".format(path))

        plt.savefig(path, dpi=400)



    
    plt.show()


# if(typeplot==6):

#     data[0] = pd.read_csv("../flightpaths/PAPERII_synobs_all.txt", sep=' ', header=None)
#     data[0].columns = ["obsAoA", "h", "hm", "d", "t", "repAoA", "azim"]

#     obsAoA = data[0]['obsAoA']
#     repAoA = data[0]['repAoA']
#     h = data[0]['h']
#     hm = data[0]['hm']
#     t = data[0]['t']

#     print(len(h))

#     d = data[0]['d']
#     azim = data[0]['azim']

#     repAoA_sph = np.arctan( ((r0+h)*np.cos(d/r0) - (r0+0.515)) / ((r0+h)*np.sin(d/r0))) * 180/np.pi
#     repAoA_m = np.arctan( ((r0+hm)*np.cos(d/r0) - (r0+0.515)) / ((r0+hm)*np.sin(d/r0))) * 180/np.pi

#     fig, ax1 = plt.subplots(figsize=(8, 6))

#     o_m = repAoA - repAoA_m

#     ax1.scatter(repAoA-repAoA_m, obsAoA, s=1.5,c=t, edgecolors="none")
#     ax1.set_ylabel("Observed AoA (deg.)")
#     ax1.set_xlabel("(O - M) refracted angle (deg.)")
#     ax1.set_xlim(-0.5, 0.5)
#     plt.axvline(0, linestyle='--')
#     plt.savefig("../plots/test.jpeg")

if (args.index_profile in [1,2,3] and args.both):

    fig, (ax1, ax2) = plt.subplots(ncols=2,sharey=True,figsize=(10,10))

    ax1.plot(data[0]['target'], data[0]['h'], color='black', label='radiosonde')
    ax1.plot(data[0]['init'], data[0]['h'], linewidth=0.5, linestyle='--', color='black', label='initial')
    ax1.plot(data[0]['retrieve'], data[0]['h'], color='green', label='0.00 RMS')
    ax1.plot(data[1]['retrieve'], data[0]['h'], color='blue', label='{} RMS'.format(noises[1]))
    ax1.plot(data[2]['retrieve'], data[0]['h'], color='red', label='{} RMS'.format(noises[2]))
    ax1.text(100, 2, s="(a)")
    ax1.set_ylabel("Altitude (km)")
    ax1.set_xlabel("Refractivity (ppm)")
    ax1.legend()

    ax2.plot(100*(data[0]['init'] - data[0]['target'])/data[0]['target'], data[0]['h'], color='black', linestyle='--', label='initial')
    ax2.plot(100*(data[0]['retrieve'] - data[0]['target'])/data[0]['target'], data[0]['h'], color='green', label='0.00 RMS')
    ax2.plot(100*(data[1]['retrieve'] - data[0]['target'])/data[0]['target'], data[0]['h'], color='blue', label='{} RMS'.format(noises[1]))
    ax2.plot(100*(data[2]['retrieve'] - data[0]['target'])/data[0]['target'], data[0]['h'], color='red', label='{} RMS'.format(noises[2]))
    
    ax2.axvline(0.0, linestyle=':', color='black')
    ax2.text(-10, 2, s="(b)")
    ax2.set_xlim(-15, 15)
    ax2.set_xlabel(r'$\frac{N_{retrieved} - N_{radiosonde}}{N_{radiosonde}}$'+' (%)')

    

    fig.tight_layout()

    path = "../plots/PAPERII_retrieve_difference_synthetic_{}_both.jpeg".format(args.index_profile)

    if(os.path.exists(path)):
    
        ftrue = input("Plot already saved as {}, do you want to overwrite? [y/n]".format(path))

        if(ftrue == 'y'):
        
            os.remove(path)

            plt.savefig(path, dpi=400)

    else:

        print("Saving image file to: {}".format(path))

        plt.savefig(path, dpi=400)

    plt.show()

if(args.ray_path):
    data = pd.read_csv("../flightpaths/PAPERII_ray_paths0_900.txt", sep=' ', header=None)

    data.columns = ['j', 's', 'h']

    fig, ax = plt.subplots()

    ax.scatter(data['s'], data['h'], s=1.5)
    ax.set_ylabel("Altitude (km)")
    ax.set_xlabel("Distance (km)")
    ax.grid(which='both')

    alt_u = np.array([1.5, 2.0, 2.5, 3.0, 4.0, 4.5, 5.0, 6.0, 7.0, 7.5])
    alt = np.array([0.5, 1.0, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 3.5, 4.0])
    dist = np.array([20, 40, 60, 80, 100, 120, 140, 160, 180, 200])

    plt.scatter(dist, alt)
    plt.scatter(dist, alt_u)

    plt.show()

if (args.sensitivity and args.profile and isinstance(args.n_err, int) and isinstance(args.h_err, int)):

    

    fig, ax1 = plt.subplots(figsize=(8, 8))

    print(args.n_err, args.h_err)

    data0 = pd.read_csv("../retrievals/PAPERII_retrieve_sens_NE_RK3_0_{}m_{}ppm_0deg.txt".format(args.h_err, args.n_err), sep=' ', header=None)
    if (args.n_err == 0):
        data1 = pd.read_csv("../retrievals/PAPERII_retrieve_sens_NE_RK3_0_-{}m_{}ppm_0deg.txt".format(args.h_err, args.n_err), sep=' ', header=None)
    elif (args.h_err == 0):
        data1 = pd.read_csv("../retrievals/PAPERII_retrieve_sens_NE_RK3_0_{}m_-{}ppm_0deg.txt".format(args.h_err, args.n_err), sep=' ', header=None)
    else:
        data1 = pd.read_csv("../retrievals/PAPERII_retrieve_sens_NE_RK3_0_-{}m_-{}ppm_0deg.txt".format(args.h_err, args.n_err), sep=' ', header=None)
    data = pd.read_csv("../retrievals/PAPERII_retrieve_sens_NE_RK3_0_0m_0ppm_0deg.txt", sep=' ', header=None)

    data5 = pd.read_csv("../retrievals/PAPERII_retrieve_sens_NE_RK3_0_5m_0ppm_0deg.txt", sep=' ', header=None)
    data5m = pd.read_csv("../retrievals/PAPERII_retrieve_sens_NE_RK3_0_-5m_0ppm_0deg.txt", sep=' ', header=None)

    data0.columns = ["retrieve", "h", "target", "init"]
    data1.columns = ["retrieve", "h", "target", "init"]
    data.columns = ["retrieve", "h", "target", "init"]

    data5.columns = ["retrieve", "h", "target", "init"]
    data5m.columns = ["retrieve", "h", "target", "init"]


    ax1.plot(data['target'], data['h'], label='radiosonde', color='green')
    ax1.plot(data['init'], data['h'], linewidth=0.5, linestyle='--', color='black', label='initial')
    ax1.plot(data['retrieve'], data['h'], color='black', label='retrieved')
    # ax1.plot(data1['retrieve'], data['h'], color='black',linewidth=0.2)
    # ax1.plot(data5m['retrieve'], data['h'], color='black',linewidth=0.2)
    ax1.fill_betweenx(data['h'], data1['retrieve'], data0['retrieve'], alpha=0.3, color='blue', label='10 m')
    ax1.fill_betweenx(data['h'], data5m['retrieve'], data5['retrieve'], alpha=0.3, color='red', label='5 m')
    

    ax1.set_ylabel("Altitude (km)")
    ax1.set_xlabel("Refractivity (ppm)")

    fig.legend()

    fig.tight_layout()

    path = "../plots/PAPERII_retrieve_sense_synthetic_0.0{}m_{}ppm_0deg.jpeg".format(args.h_err, args.n_err)

    if(os.path.exists(path)):
    
        ftrue = input("Plot already saved as {}, do you want to overwrite? [y/n]".format(path))

        if(ftrue == 'y'):
        
            os.remove(path)

            plt.savefig(path, dpi=400)

    else:

        print("Saving image file to: {}".format(path))

        plt.savefig(path, dpi=400)



    
    plt.show()    

if (args.may) and (args.paths): 

    datar[0] = pd.read_csv("../flightpaths/May8_flightpath_-10_to_10_t0_1800.txt", sep=' ', header=None)
    datar[1] = pd.read_csv("../flightpaths/May8_flightpath_-10_to_10_t1800_3600.txt", sep=' ', header=None)
    datar[2] = pd.read_csv("../flightpaths/May8_flightpath_-10_to_10_t3600_5400.txt", sep=' ', header=None)
    datar[3] = pd.read_csv("../flightpaths/May8_flightpath_-10_to_10_t5400_7200.txt", sep=' ', header=None)

    datar[0].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    datar[1].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    datar[2].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    datar[3].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    
    obsAoA1, obsAoA2, obsAoA3, obsAoA4 = datar[0]['obsAoA'], datar[1]['obsAoA'], datar[2]['obsAoA'], datar[3]['obsAoA']
    repAoA1, repAoA2, repAoA3, repAoA4 = datar[0]['repAoA'], datar[1]['repAoA'], datar[2]['repAoA'], datar[3]['repAoA']
    h1, h2, h3, h4 = datar[0]['h'], datar[1]['h'], datar[2]['h'], datar[3]['h']
    d1, d2, d3, d4 = datar[0]['d'], datar[1]['d'], datar[2]['d'], datar[3]['d']
    hinit1, hinit2, hinit3, hinit4  = datar[0]['hinit'], datar[1]['hinit'], datar[2]['hinit'], datar[3]['hinit']
    hret1, hret2, hret3, hret4 = datar[0]['hret'], datar[1]['hret'], datar[2]['hret'], datar[3]['hret']
    azim1, azim2, azim3, azim4 = datar[0]['azim'], datar[1]['azim'], datar[2]['azim'], datar[3]['azim']

    repAoA_sph1, repAoA_sph2, repAoA_sph3, repAoA_sph4 = rep(h1, d1), rep(h2, d2), rep(h3, d3), rep(h4, d4)
    repAoA_init1, repAoA_init2, repAoA_init3, repAoA_init4 = rep(hinit1, d1), rep(hinit2, d2), rep(hinit3, d3), rep(hinit4, d4)
    repAoA_ret1, repAoA_ret2, repAoA_ret3, repAoA_ret4 = rep(hret1, d1), rep(hret2, d2), rep(hret3, d3), rep(hret4, d4)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2,nrows=2,sharey=True,sharex=True,figsize=(10,10))

    ax1.scatter(-azim1, repAoA1, color='black', label='obs', s=1.5)
    ax1.scatter(-azim1, repAoA_init1, linewidth=0.5, color='blue', label='initial',s=1.5)
    ax1.scatter(-azim1, repAoA_ret1, color='green', label='retrieved',s=1.5)
    # plot = ax1.scatter(-azim1, repAoA1, c=d1,s=5, edgecolors='none')
    # ax1.scatter(d1, h1, s=3, color='black')

    ax1.set_ylabel("Reported AoA (deg.)")
    # ax1.set_ylabel("Altitude (km)")
    ax1.set_xlabel("(a)")

    print("Maximum difference between ellipsoid and spheroid repAoA for plot 1:", max(repAoA1 - repAoA_sph1))
    print("Maximum difference between ellipsoid and spheroid repAoA for plot 2:", max(repAoA2 - repAoA_sph2))
    print("Maximum difference between ellipsoid and spheroid repAoA for plot 3:", max(repAoA3 - repAoA_sph3))
    print("Maximum difference between ellipsoid and spheroid repAoA for plot 4:", max(repAoA4 - repAoA_sph4))

    ax1.set_title('$0 \leq t < 30 $ (mins.)')

    ax2.scatter(-azim2, repAoA2, color='black',s=1.5)
    ax2.scatter(-azim2, repAoA_init2, linewidth=0.5, color='blue',s=1.5)
    ax2.scatter(-azim2, repAoA_ret2, color='green',s=1.5)
    # plot = ax2.scatter(-azim2, repAoA2, c=d2, s=5, edgecolors='none')
    # ax2.scatter(d2, h2, s=3, color='black')
    ax2.set_xlabel("(b)")

    ax2.sharey(ax1)

    ax2.set_title('$30 \leq t < 60$ (mins.)')

    ax3.scatter(-azim3, repAoA3, color='black', s=1.5)
    ax3.scatter(-azim3, repAoA_init3, linewidth=0.5, color='blue',s=1.5)
    ax3.scatter(-azim3, repAoA_ret3, color='green',s=1.5)
    # plot = ax3.scatter(-azim3, repAoA3, c=d3, s=5, edgecolors='none')
    # ax3.scatter(d3, h3, s=3, color='black')

    ax3.set_xlabel("Horizontal angle (deg.) \n (c)")
    ax4.set_xlabel("Horizontal angle (deg.) \n (d)")
    # ax3.set_xlabel("Distance (km) \n (c)")
    # ax4.set_xlabel("Distance (km) \n (d)")
    # ax3.set_ylabel("Altitude (km)")
    ax3.set_ylabel("Reported AoA (deg.)")

    ax1.sharex(ax3)

    ax3.set_title('$60 \leq t < 90$ (mins.)')

    ax4.scatter(-azim4, repAoA4, color='black', s=1.5)
    ax4.scatter(-azim4, repAoA_init4, linewidth=0.5, color='blue',s=1.5)
    ax4.scatter(-azim4, repAoA_ret4, color='green',s=1.5)
    # plot = ax4.scatter(-azim4, repAoA4, c=d4, s=5, edgecolors='none')
    # ax4.scatter(d4, h4, s=3, color='black')

    ax4.set_title('$90 \leq t < 120$ (mins.)')

    # cbar_ax = fig.add_axes([0.95, 0.15, 0.05, 0.7])
    # fig.colorbar(plot, cax=cbar_ax)

    fig.legend(bbox_to_anchor=(0.3, 0.68), markerscale=2)
    fig.tight_layout()

    
    path = "../plots/May8_retrieve_dist.jpeg"

    if(os.path.exists(path)):
    
        ftrue = input("Plot already saved as {}, do you want to overwrite? [y/n]".format(path))

        if(ftrue == 'y'):
        
            os.remove(path)

            fig.savefig(path, dpi=400)

    else:

        print("Saving image file to: {}".format(path))

        fig.savefig(path, dpi=400)
    
    plt.show()

if (args.may) and (args.profile): 

    data = pd.read_csv("../refractivity/9May_00z_2024_Watnall_profile_RH.txt", sep=' ', header=None)

    data.columns = ['h', 'N', 'Nd', 'Nw']

    data2 = pd.read_csv("../refractivity/8May_00z_2024_Watnall_profile_RH.txt", sep=' ', header=None)

    data2.columns = ['h', 'N', 'Nd', 'Nw']


    datar[0] = pd.read_csv("../retrievals/May8_retrieve_-10_to_10_t0_1800.txt", sep=' ', header=None)
    datar[1] = pd.read_csv("../retrievals/May8_retrieve_-10_to_10_t1800_3600.txt", sep=' ', header=None)
    datar[2] = pd.read_csv("../retrievals/May8_retrieve_-10_to_10_t3600_5400.txt", sep=' ', header=None)
    datar[3] = pd.read_csv("../retrievals/May8_retrieve_-10_to_10_t5400_7200.txt", sep=' ', header=None)
 
    datar[0].columns = ["retrieve", "h", "target", "init", "ndry"]
    datar[1].columns = ["retrieve", "h", "target", "init", "ndry"]
    datar[2].columns = ["retrieve", "h", "target", "init", "ndry"]
    datar[3].columns = ["retrieve", "h", "target", "init", "ndry"]


    fig, ax1 = plt.subplots(figsize=(8, 8))

    # ax1.plot(datar[0]['target'], datar[0]['h'], color='black', label='UKV')
    # ax1.plot(datar[0]['ndry'], datar[0]['h'], color='black', linestyle=':', label='UKV Ndry')
    ax1.plot(datar[0]['init'], datar[0]['h'], linewidth=0.5, linestyle='--', color='black', label='initial')
    ax1.plot(datar[0]['retrieve'], datar[0]['h'], color='green', label='0-30 mins')
    ax1.plot(datar[1]['retrieve'], datar[1]['h'], color='red', label='30-60 mins')
    ax1.plot(datar[2]['retrieve'], datar[2]['h'], color='blue', label='60-90 mins')
    ax1.plot(datar[3]['retrieve'], datar[3]['h'], color='yellow', label='90-120 mins')
    ax1.plot(data['N'], data['h']/1e3, label='radiosonde (00z 9 May)', color='black')
    ax1.plot(data2['N'], data2['h']/1e3, label='radiosonde (00z 8 May', color='orange')

    plt.legend()
    plt.ylim(0, 12)
    ax1.set_ylabel("Altitude (km)")
    ax1.set_xlabel("Refractivity (ppm)")

    
    path = "../plots/PAPERII_retrieve_profiles_May8_2.jpeg"


    if(os.path.exists(path)):
    
        ftrue = input("Plot already saved as {}, do you want to overwrite? [y/n]".format(path))

        if(ftrue == 'y'):
        
            os.remove(path)

            plt.savefig(path, dpi=400)

    else:

        print("Saving image file to: {}".format(path))

        plt.savefig(path, dpi=400)
    
    plt.show() 

if (args.may) and (args.model): 

    datar[0] = pd.read_csv("../flightpaths/May_NE_ray_paths.txt", sep=' ', header=None)
 
    datar[0].columns = ["j", "hm", "um", "sm", "h", "azim", "t"]

    repm = rep(datar[0]['hm'], datar[0]['sm'])
    rep = rep(datar[0]['h'], datar[0]['sm'])

    fig, ax1 = plt.subplots(figsize=(8, 8))

    # ax1.scatter(np.arcsin(datar[0]['um'])*180/np.pi-rep, np.arcsin(datar[0]['um'])*180/np.pi, color='black', s=1.5, label='reported')
    # ax1.scatter(np.arcsin(datar[0]['um'])*180/np.pi-repm, np.arcsin(datar[0]['um'])*180/np.pi, color='blue', s=1.5, label='model')
    plot = ax1.scatter(rep-repm, np.arcsin(datar[0]['um'])*180/np.pi, s=1.5, c=datar[0]['sm'], edgecolors='none')

    plt.legend()
    ax1.set_xlim(-0.5, 0.5)
    plt.axvline(0)
    ax1.set_ylabel("Altitude (km)")
    ax1.set_xlabel("Refractivity (ppm)")

    cbar = plt.colorbar(plot)
    cbar.set_label("Range (km)")
    
    path = "../plots/May_refrac_model_vs_rep_diff.jpeg"


    if(os.path.exists(path)):
    
        ftrue = input("Plot already saved as {}, do you want to overwrite? [y/n]".format(path))

        if(ftrue == 'y'):
        
            os.remove(path)

            plt.savefig(path, dpi=400)

    else:

        print("Saving image file to: {}".format(path))

        plt.savefig(path, dpi=400)
    
    plt.show() 

if (args.igarss):

    
    datar[0] = pd.read_csv("../retrievals/IGARSS_retrieve_0.0_0.25hr.txt", sep=' ', header=None)
    datar[1] = pd.read_csv("../retrievals/IGARSS_retrieve_0.5_0.75hr.txt", sep=' ', header=None)
    datar[2] = pd.read_csv("../retrievals/IGARSS_retrieve_1.0_1.25hr.txt", sep=' ', header=None)

    data2 = pd.read_csv("../refractivity/OCT25_N_15.00_(km)_updated_orog.txt", sep=' ', header=None)
    data3 = pd.read_csv("../refractivity/OCT25_N_17.00_(km)_updated_orog.txt", sep=' ', header=None)
    data4 = pd.read_csv("../refractivity/OCT25_N_16.00_(km)_updated_orog.txt", sep=' ', header=None)

    datar[0].columns = ["retrieve", "h", "target", "init", "ndry"]
    datar[1].columns = ["retrieve", "h", "target", "init", "ndry"]
    datar[2].columns = ["retrieve", "h", "target", "init", "ndry"]

    data2.columns = ['h', 'N', 'Nd', 'Nw']
    data3.columns = ['h', 'N', 'Nd', 'Nw']
    data4.columns = ['h', 'N', 'Nd', 'Nw']

    fig, ax1 = plt.subplots(figsize=(8, 8))

    # ax1.plot(datar[0]['target'], datar[0]['h'], color='black', label='UKV')
    # ax1.plot(datar[0]['init'], datar[0]['h'], linewidth=0.5, linestyle='--', color='black', label='initial')
    ax1.plot(datar[0]['retrieve'], datar[0]['h'], color='green', label='15:30-15:45')
    ax1.plot(datar[1]['retrieve'], datar[1]['h'], color='red', label='16:00-16:15')
    ax1.plot(datar[2]['retrieve'], datar[2]['h'], color='blue', label='16:30-16:45')
    ax1.plot(data2['N'], data2['h']/1e3, linewidth=0.8, linestyle='--', color='green', label='15:00 UKV')
    ax1.plot(data4['N'], data4['h']/1e3, linewidth=0.8, linestyle='--', color='red', label='16:00 UKV')
    ax1.plot(data3['N'], data3['h']/1e3, linewidth=0.8, linestyle='--', color='blue', label='17:00 UKV')
    plt.ylim(0, 12)
    plt.xlim(50, 350)
    plt.legend()

    ax1.set_ylabel("Altitude (km)")
    ax1.set_xlabel("Refractivity (ppm)")

    
    path = "../plots/IGARSS_profile.jpeg"
   

    if(os.path.exists(path)):
    
        ftrue = input("Plot already saved as {}, do you want to overwrite? [y/n]".format(path))

        if(ftrue == 'y'):
        
            os.remove(path)

            plt.savefig(path, dpi=400)

    else:

        print("Saving image file to: {}".format(path))

        plt.savefig(path, dpi=400)
    
    plt.show()

if (args.igarss_paths):


    datar[0] = pd.read_csv("../flightpaths/IGARSS_retrieve_paths_0.0_0.25hr.txt", sep=' ', header=None)
    datar[1] = pd.read_csv("../flightpaths/IGARSS_retrieve_paths_0.5_0.75hr.txt", sep=' ', header=None)
    datar[2] = pd.read_csv("../flightpaths/IGARSS_retrieve_paths_1.0_1.25hr.txt", sep=' ', header=None)
    datar[3] = pd.read_csv("../flightpaths/IGARSS_retrieve_paths_1.0_1.25hr.txt", sep=' ', header=None)

    datar[0].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    datar[1].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    datar[2].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]
    datar[3].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]


    obsAoA1, obsAoA2, obsAoA3, obsAoA4 = datar[0]['obsAoA'], datar[1]['obsAoA'], datar[2]['obsAoA'], datar[3]['obsAoA']
    repAoA1, repAoA2, repAoA3, repAoA4 = datar[0]['repAoA'], datar[1]['repAoA'], datar[2]['repAoA'], datar[3]['repAoA']
    h1, h2, h3, h4 = datar[0]['h'], datar[1]['h'], datar[2]['h'], datar[3]['h']
    d1, d2, d3, d4 = datar[0]['d'], datar[1]['d'], datar[2]['d'], datar[3]['d']
    hinit1, hinit2, hinit3, hinit4  = datar[0]['hinit'], datar[1]['hinit'], datar[2]['hinit'], datar[3]['hinit']
    hret1, hret2, hret3, hret4 = datar[0]['hret'], datar[1]['hret'], datar[2]['hret'], datar[3]['hret']
    azim1, azim2, azim3, azim4 = datar[0]['azim'], datar[1]['azim'], datar[2]['azim'], datar[3]['azim']

    repAoA_sph1, repAoA_sph2, repAoA_sph3, repAoA_sph4 = rep(h1, d1), rep(h2, d2), rep(h3, d3), rep(h4, d4)
    repAoA_init1, repAoA_init2, repAoA_init3, repAoA_init4 = rep(hinit1, d1), rep(hinit2, d2), rep(hinit3, d3), rep(hinit4, d4)
    repAoA_ret1, repAoA_ret2, repAoA_ret3, repAoA_ret4 = rep(hret1, d1), rep(hret2, d2), rep(hret3, d3), rep(hret4, d4)

    fig, ax1 = plt.subplots(figsize=(10,10))#ncols=2,nrows=2,sharey=True,sharex=True,

    ax1.scatter(-azim2, repAoA2, color='black', label='obs', s=1.5)
    ax1.scatter(-azim2, repAoA_init2, linewidth=0.5, color='blue', label='initial',s=1.5)
    ax1.scatter(-azim2, repAoA_ret2, color='green', label='retrieved',s=1.5)
    # ax1.scatter(-azim1, repAoA_sph1, color='skyblue', label='retrieved',s=1.5)
    ax1.set_ylabel("Reported AoA (deg.)")
    #ax1.set_xlabel("(a)")
    ax1.set_xlabel("Horizontal angle (deg.)")
    print("Maximum difference between ellipsoid and spheroid repAoA for plot 1:", max(repAoA1 - repAoA_sph1))

    ax1.set_title('$30 \leq t < 45 $ (mins.)')

    # ax2.scatter(-azim2, repAoA2, color='black',s=1.5)
    # ax2.scatter(-azim2, repAoA_init2, linewidth=0.5, color='blue',s=1.5)
    # ax2.scatter(-azim2, repAoA_ret2, color='green',s=1.5)
    # # ax2.scatter(-azim2, repAoA_sph2, color='skyblue', label='retrieved',s=1.5)

    # ax2.set_xlabel("(b)")

    # ax2.sharey(ax1)

    # ax2.set_title('$15 \leq t < 30$ (mins.)')

    # ax3.scatter(-azim3, repAoA3, color='black', s=1.5)
    # ax3.scatter(-azim3, repAoA_init3, linewidth=0.5, color='blue',s=1.5)
    # ax3.scatter(-azim3, repAoA_ret3, color='green',s=1.5)
    # # ax3.scatter(-azim3, repAoA_sph3, color='skyblue', label='retrieved',s=1.5)

    # ax3.set_xlabel("Horizontal angle (deg.) \n (c)")
    # ax4.set_xlabel("Horizontal angle (deg.) \n (d)")
    # ax3.set_ylabel("Reported AoA (deg.)")

    # ax1.sharex(ax3)

    # ax3.set_title('$30 \leq t < 45$ (mins.)')

    # ax4.scatter(-azim4, repAoA4, color='black', s=1.5)
    # ax4.scatter(-azim4, repAoA_init4, linewidth=0.5, color='blue',s=1.5)
    # ax4.scatter(-azim4, repAoA_ret4, color='green',s=1.5)
    # # ax4.scatter(-azim4, repAoA_sph4, color='skyblue', label='retrieved',s=1.5)

    # ax4.set_title('$45 \leq t < 60$ (mins.)')

    fig.legend(bbox_to_anchor=(0.3, 0.68), markerscale=2)
    fig.tight_layout()

    # if(args.south):
    #     path = "../plots/PAPERII_retrieve_oct_S_flightpaths.jpeg"
    # else:
    #     path = "../plots/PAPERII_retrieve_flightpaths.jpeg"

    # if(os.path.exists(path)):
    
    #     ftrue = input("Plot already saved as {}, do you want to overwrite? [y/n]".format(path))

    #     if(ftrue == 'y'):
        
    #         os.remove(path)

    #         fig.savefig(path, dpi=400)

    # else:

    #     print("Saving image file to: {}".format(path))

    #     fig.savefig(path, dpi=400)
    
    plt.show()