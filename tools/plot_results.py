import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import smplotlib
import sys
import pandas as pd
import argparse
import os
import statistics as stat

r0 = 6383.5713292490545
h0 = 0.575

rep = lambda h, d : np.arctan( ((r0+h)*np.cos(d/r0) - (r0+h0)) / ((r0+h)*np.sin(d/r0))) * 180/np.pi

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

parser.add_argument('--model', action='store_true',
                    help='Modelled refraction')      
   
parser.add_argument('--south', action='store_true',
                    help='South side data')  

args = parser.parse_args()

data = [0]*3
dataerr1 = [0]*20
dataerr2 = [0]*20
dataerr3 = [0]*20
datar = [0]*4
datar4 = [0]*2
noises = [0, 0.01, 0.05]

if (args.index_profile in [1,2,3]):

    for i in range(3):
        
        data[i] = pd.read_csv("../retrievals/PAPERII_retrieve_NE_RK3_{}_{}_575_test.txt".format(args.index_profile, noises[i]), sep=' ', skiprows=1, header=None)
        data[i].columns = ["retrieve", "h", "target", "init"]

if (args.loss):
    for i in range(3):
        
        data[i] = pd.read_csv("../RMS/PAPERII_RMS_NE_RK3_{}.txt".format(noises[i]), sep=' ', header=None)
        data[i].columns = ["iter", "loss", "Nrms"]


if (args.true and args.loss):

    if(args.south):
        datar[0] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t0_900.txt", sep=' ', header=None)
        datar[1] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t900_1800.txt", sep=' ', header=None)
        datar[2] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t1800_2700.txt", sep=' ', header=None)
        datar[3] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_oct_S_t2700_3600.txt", sep=' ', header=None)
    else:
        datar[0] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t0_900_updated_check.txt", sep=' ', header=None)
        datar[1] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t900_1800_updated_check.txt", sep=' ', header=None)
        datar[2] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t1800_2700_updated_check.txt", sep=' ', header=None)
        datar[3] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t2700_3600_updated_check.txt", sep=' ', header=None)

    datar[0].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim", "time"]
    datar[1].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim", "time"]
    datar[2].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim", "time"]
    datar[3].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim", "time"]

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

    n = 50

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
    ax1.plot(valsp, mean_err_init1, color='blue', label='initial')
    ax1.axhline(0.0, color='black')
   
    ax1.set_ylabel("LoS AoA difference (deg.)")#\n(Aircraft altitude ($10^2$ km))")
    ax1.set_xlabel("(a)")

    ax1.set_title('$0 \leq t < 15 $ (mins.)')

    ax2.plot(valsp, mean_err2, color='green')
    ax2.plot(valsp, mean_err_init2, color='blue')
    ax2.axhline(0.0, color='black')
    ax2.set_xlabel("(b)")

    ax2.set_title('$15 \leq t < 30$ (mins.)')

    ax3.plot(valsp, mean_err3, color='green')
    ax3.plot(valsp, mean_err_init3, color='blue')
    ax3.axhline(0.0, color='black')

    ax3.set_xlabel("Distance (km) \n (c)")
    ax4.set_xlabel("Distance (km) \n (d)")
    ax3.set_ylabel("LoS AoA difference (deg.)")

    ax1.sharex(ax3)

    ax3.set_title('$30 \leq t < 45$ (mins.)')

    ax4.plot(valsp, mean_err4, color='green')
    ax4.plot(valsp, mean_err_init4, color='blue')
    ax4.axhline(0.0, color='black')
    ax4.set_title('$45 \leq t < 60$ (mins.)')

    fig.legend(bbox_to_anchor=(0.49, 0.9), markerscale=2)
    fig.tight_layout()

    if(args.south):
        path = "../plots/PAPERII_retrieve_flightpaths_oct_S_error.jpeg"
    else:
        path = "../plots/PAPERII_retrieve_flightpaths_error_dist.jpeg"

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


    datar[0] = pd.read_csv("../flightpath/PAPERII_retrieve_paths_t0_900_updated_check.txt", sep=' ', header=None)
    datar[1] = pd.read_csv("../flightpath/PAPERII_retrieve_paths_t900_1800_updated_check.txt", sep=' ', header=None)
    datar[2] = pd.read_csv("../flightpath/PAPERII_retrieve_paths_t1800_2700_updated_check.txt", sep=' ', header=None)
    datar[3] = pd.read_csv("../flightpath/PAPERII_retrieve_paths_t2700_3600_updated_check.txt", sep=' ', header=None)

    datar[0].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim", "time"]
    datar[1].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim", "time"]
    datar[2].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim", "time"]
    datar[3].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim", "time"]

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
    print("Retrieved mean: ", stat.mean(error), "Retrieved std: ", stat.stdev(error))
    print("Initial mean: ", stat.mean(error_init), "Initial std: ", stat.stdev(error_init))
    

    
    
    fig, ax1 = plt.subplots(figsize=(8, 6))


    ax1.hist(error, bins=500, 
                   histtype=u'step', density=True, color='green', label='retrieved') 
    ax1.hist(error_init, bins=500, 
                   histtype=u'step', density=True, color='blue', label='initial') 

    ax1.set_xlim(-max(error_init), max(error_init))
    handle1 = Line2D([], [], c='r')
    handle2 = Line2D([], [], c='b')
    handles, labels = ax1.get_legend_handles_labels()
    new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]
    ax1.set_xlim(-0.15, 0.15)
    plt.legend(handles=new_handles, labels=labels)

    plt.axvline(0.0, linestyle='--', color='black')

    ax1.set_xlabel("LoS AoA difference (deg.)")
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


if (args.paths):

    datar[0] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t0_900_updated_check.txt", sep=' ', header=None)
    datar[1] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t900_1800_updated_check.txt", sep=' ', header=None)
    datar[2] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t1800_2700_updated_check.txt", sep=' ', header=None)
    datar[3] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths_t2700_3600_updated_check.txt", sep=' ', header=None)

    datar[0].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim", "time"]
    datar[1].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim", "time"]
    datar[2].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim", "time"]
    datar[3].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim", "time"]

    print("Number of obs in 1st sample: ", len(datar[0]))
    print("Number of obs in 2nd sample: ", len(datar[1]))
    print("Number of obs in 3rd sample: ", len(datar[2]))
    print("Number of obs in 4th sample: ", len(datar[3]))

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

    ax1.scatter(-azim1, repAoA1, color='black', label='AoA$_{LoS}$', s=1.5)
    ax1.scatter(-azim1, repAoA_init1, linewidth=0.5, color='blue', label='AoA$_{init\_model}$',s=1.5)
    ax1.scatter(-azim1, repAoA_ret1, color='green', label='AoA$_{ret\_model}$',s=1.5)
    # ax1.scatter(-azim1, repAoA_sph1, color='skyblue', label='retrieved',s=1.5)
    ax1.set_ylabel("LoS AoA (deg.)")
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

    ax3.set_xlabel("Azimuthal angle (deg.) \n (c)")
    ax4.set_xlabel("Azimuthal angle (deg.) \n (d)")
    ax3.set_ylabel("LoS AoA (deg.)")

    ax1.sharex(ax3)

    ax3.set_title('$30 \leq t < 45$ (mins.)')

    ax4.scatter(-azim4, repAoA4, color='black', s=1.5)
    ax4.scatter(-azim4, repAoA_init4, linewidth=0.5, color='blue',s=1.5)
    ax4.scatter(-azim4, repAoA_ret4, color='green',s=1.5)
    # ax4.scatter(-azim4, repAoA_sph4, color='skyblue', label='retrieved',s=1.5)

    ax4.set_title('$45 \leq t < 60$ (mins.)')

    fig.legend(bbox_to_anchor=(0.3, 0.68), markerscale=2)
    fig.tight_layout()

    path = "../plots/PAPERII_retrieve_flightpaths_updated.jpeg"

    if(os.path.exists(path)):
    
        ftrue = input("Plot already saved as {}, do you want to overwrite? [y/n]".format(path))

        if(ftrue == 'y'):
        
            os.remove(path)

            fig.savefig(path, dpi=400)

    else:
         
       print("Saving image file to: {}".format(path))

       fig.savefig(path, dpi=400)
    
    plt.show()

if (args.index_profile in [1,2,3] and args.both):

    data5 = pd.read_csv("../retrievals/PAPERII_retrieve_sens_NE_RK3_0_5m_0ppm_0deg_575.txt", sep=' ', header=None)
    data5m = pd.read_csv("../retrievals/PAPERII_retrieve_sens_NE_RK3_0_-5m_0ppm_0deg_575.txt", sep=' ', header=None)

    data10 = pd.read_csv("../retrievals/PAPERII_retrieve_sens_NE_RK3_0_10m_0ppm_0deg_575.txt", sep=' ', header=None)
    data10m = pd.read_csv("../retrievals/PAPERII_retrieve_sens_NE_RK3_0_-10m_0ppm_0deg_575.txt", sep=' ', header=None)

    for i in range(20):
        
        dataerr2[i] = pd.read_csv("../retrievals/PAPERII_retrieve_NE_RK3_{}_0.01_575_sample{}.txt".format(args.index_profile, i), sep=' ', skiprows=1, header=None)
        dataerr2[i].columns = ["retrieve", "h", "target", "init"]

        dataerr3[i] = pd.read_csv("../retrievals/PAPERII_retrieve_NE_RK3_{}_0.05_575_sample{}.txt".format(args.index_profile, i), sep=' ', skiprows=1, header=None)
        dataerr3[i].columns = ["retrieve", "h", "target", "init"]
    
    err_N2, err_N3 = [], []
    abs_N2, abs_N3 = [], []
    max_prof3, min_prof3 = [], []

    N2_mean = []
    N3_mean = []

    std_N2 = []
    std_N3 = []

    for j in range(len(dataerr2[0]["retrieve"])):

        data2_val = [dataerr2[k]["retrieve"][j] for k in range(20)]
        data3_val = [dataerr3[k]["retrieve"][j] for k in range(20)]

        N2_mean.append(np.mean(data2_val))
        N3_mean.append(np.mean(data3_val))
        
        err_N2.append([np.abs(-np.min(data2_val)+N2_mean[j]),
                       np.abs(np.max(data2_val)-N2_mean[j])])
        err_N3.append([np.abs(-np.min(data3_val)+N3_mean[j]),
                       np.abs(np.max(data3_val)-N3_mean[j])])
        
        std_N2.append(np.std(data2_val))
        std_N3.append(np.std(data3_val))
        
        abs_N2.append([100*np.abs(-np.min(data2_val) + N2_mean[j])/N2_mean[j],
                       100*np.abs((np.max(data2_val) - N2_mean[j])/N2_mean[j])])
        abs_N3.append([100*np.abs(-np.min(data3_val) + N3_mean[j])/N3_mean[j],
                       100*np.abs((np.max(data3_val) - N3_mean[j])/N3_mean[j])])             
        
        # N_mean.append(np.mean())
   
    data10.columns = ["retrieve", "h", "target", "init"]
    data10m.columns = ["retrieve", "h", "target", "init"]

    data5.columns = ["retrieve", "h", "target", "init"]
    data5m.columns = ["retrieve", "h", "target", "init"]

    fig, (ax1, ax2) = plt.subplots(ncols=2,sharey=True,figsize=(12,10))

    ax1.plot(data[0]['target'], data[0]['h'], color='black', label='radiosonde')
    ax1.plot(data[0]['init'], data[0]['h'], linewidth=0.5, linestyle='--', color='black', label='initial')
    ax1.errorbar(data[0]['retrieve'], data[0]['h'], color='green', label='0.00$\degree$ RMSE')
    ax1.errorbar(N2_mean, data[0]['h'], xerr=np.transpose(err_N2), color='blue', label='0.01$\degree$ RMSE')
    ax1.errorbar(N3_mean, data[0]['h'], xerr=np.transpose(err_N3), color='red', label='0.05$\degree$ RMSE')
    # ax1.plot(max_prof3, data[0]['h'], color='black')
    # ax1.plot(min_prof3, data[0]['h'], color='black')
    
    ax1.fill_betweenx(data[0]['h'], data10m['retrieve'], data10['retrieve'], alpha=0.2, color='blue', label='10 m')
    ax1.fill_betweenx(data[0]['h'], data5m['retrieve'], data5['retrieve'], alpha=0.2, color='red', label='5 m')
    
    
    ax1.legend()
    ax1.text(100, 2, s="(a)")

    
    ax1.set_ylabel("Altitude (km)")
    ax1.set_xlabel("Refractivity (ppm)")

    ax2.plot(100*(data[0]['init'] - data[0]['target'])/data[0]['target'], data[0]['h'], color='black', linestyle='--', linewidth=1)
    ax2.plot(100*(data[0]['retrieve'] - data[0]['target'])/data[0]['target'], data[0]['h'], color='green')
    ax2.plot(100*(N2_mean - data[0]['target'])/data[0]['target'], data[0]['h'], color='blue', label='0.01$\degree$ RMSE')#, xerr=np.transpose(abs_N2))
    ax2.plot(100*(N3_mean - data[0]['target'])/data[0]['target'], data[0]['h'], color='red', label='0.05$\degree$ RMSE')#, xerr=np.transpose(abs_N3))
    ax2.fill_betweenx(data[0]['h'], 100*(data10m['retrieve']-data[0]['target'])/data[0]['target'], 100*(data10['retrieve']-data[0]['target'])/data[0]['target'], alpha=0.2, color='blue')
    ax2.fill_betweenx(data[0]['h'], 100*(data5m['retrieve']-data[0]['target'])/data[0]['target'], 100*(data5['retrieve']-data[0]['target'])/data[0]['target'], alpha=0.2, color='red')
    
    ax3 = ax2.twiny()
    ax3.plot(100*(std_N2/data[0]['target']), data[0]['h'], linestyle=':', linewidth=1, color='blue', label='0.01$\degree$ RMSE SD')
    ax3.plot(100*(std_N3/data[0]['target']), data[0]['h'], linestyle=':', linewidth=1, color='red', label='0.05$\degree$ RMSE SD')
    
    ax3.set_xlim(-1.0, 1.0)
    ax3.set_xticks([0,0.5,1,1.0])
    ax3.legend(loc='upper left')
    ax3.set_xlabel('Refractivity error SD (%)')
    ax2.axvline(0.0, color='black')
    ax2.text(-10, 2, s="(b)")
    ax2.set_xlim(-15, 15)
    ax2.set_xlabel('Refractivity error (%)')

    print(np.std(abs_N2))
    print(np.std(abs_N3))

    num = len(data[0]['h'])
    RMS1, RMS2, RMS3, RMSinit, RMSE10, RMSEm10, RMSE5, RMSEm5 = 0,0,0,0,0,0,0,0
    for i in range(num):
        RMS1 += (data[0]['retrieve'][i] - data[0]['target'][i])**2
        RMS2 += (N2_mean[i] - data[0]['target'][i])**2
        RMS3 += (N3_mean[i] - data[0]['target'][i])**2
        RMSinit += (data[0]['init'][i] - data[0]['target'][i])**2
        RMSE10 += (data10['retrieve'][i]-data[0]['target'][i])**2
        RMSEm10 += (data10m['retrieve'][i]-data[0]['target'][i])**2
        RMSE5 += (data5['retrieve'][i]-data[0]['target'][i])**2
        RMSEm5 += (data5m['retrieve'][i]-data[0]['target'][i])**2

    print(np.sqrt(RMS1/num))
    print(np.sqrt(RMS2/num))
    print(np.sqrt(RMS3/num))
    print(np.sqrt(RMSinit/num))
    print("+10m:", np.sqrt(RMSE10/num))
    print("-10m:",np.sqrt(RMSEm10/num))
    print("+5m:",np.sqrt(RMSE5/num))
    print("-5m:",np.sqrt(RMSEm5/num))

    print("Start height: ", data[0]['h'][0])
    

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