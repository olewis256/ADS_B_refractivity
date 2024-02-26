import numpy as np
import matplotlib.pyplot as plt
import smplotlib
import sys
import pandas as pd

r0 = 6364.765

rep = lambda h, d : np.arctan( ((r0+h)*np.cos(d/r0) - (r0+0.527)) / ((r0+h)*np.sin(d/r0))) * 180/np.pi

typeplot = int(sys.argv[1])

if len(sys.argv) == 3:  
    typeplot2 = int(sys.argv[2])


data = [0]*3
datar = [0]*4
noises = [0, 0.01, 0.05]

if (typeplot == 1) or (typeplot == 2):

    for i in range(3):
        
        data[i] = pd.read_csv("../retrievals/PAPERII_retrieve_NE_RK3_{}_{}.txt".format(typeplot2, noises[i]), sep=' ', header=None)
        data[i].columns = ["retrieve", "h", "target", "init"]

if (typeplot == 3):
    for i in range(3):
        
        data[i] = pd.read_csv("../RMS/PAPERII_RMS_NE_RK3_{}.txt".format(noises[i]), sep=' ', header=None)
        data[i].columns = ["iter", "loss", "Nrms"]


if (typeplot == 4):

    datar[0] = pd.read_csv("../retrievals/PAPERII_retrieve_NE_TRUE_t0_900.txt", sep=' ', header=None)
    datar[1] = pd.read_csv("../retrievals/PAPERII_retrieve_NE_TRUE_t900_1800.txt", sep=' ', header=None)
    datar[2] = pd.read_csv("../retrievals/PAPERII_retrieve_NE_TRUE_t1800_2700.txt", sep=' ', header=None)
    datar[3] = pd.read_csv("../retrievals/PAPERII_retrieve_NE_TRUE_t2700_3600.txt", sep=' ', header=None)
    
    datar[0].columns = ["retrieve", "h", "target", "init", "ndry"]
    datar[1].columns = ["retrieve", "h", "target", "init", "ndry"]
    datar[2].columns = ["retrieve", "h", "target", "init", "ndry"]
    datar[3].columns = ["retrieve", "h", "target", "init", "ndry"]

    fig, ax1 = plt.subplots(figsize=(8, 6))


    ax1.plot(datar[0]['target'], datar[0]['h'], color='black', label='radiosonde')
    ax1.plot(datar[0]['init'], datar[0]['h'], linewidth=0.5, linestyle='--', color='black', label='initial')
    ax1.plot(datar[0]['retrieve'], datar[0]['h'], color='green', label='0-15 mins')
    ax1.plot(datar[1]['retrieve'], datar[1]['h'], color='red', label='15-30 mins')
    ax1.plot(datar[2]['retrieve'], datar[2]['h'], color='blue', label='30-45 mins')
    ax1.plot(datar[3]['retrieve'], datar[3]['h'], color='deepskyblue', label='45-60 mins')
    # ax1.plot(data[0]['ndry'], data[0]['h'], color='blue', label='dry refrac.')

    plt.legend()

    ax1.set_ylabel("Altitude (km)")
    ax1.set_xlabel("Refractivity (ppm)")

    plt.savefig("../plots/PAPERII_retrieve_profiles_above150.jpeg", dpi=700)
    
    plt.show()

if (typeplot == 45):

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
    # ax1.scatter(-azim1, repAoA_init1, linewidth=0.5, color='blue', label='initial',s=1.5)
    # ax1.scatter(-azim1, repAoA_ret1, color='green', label='retrieved',s=1.5)
    ax1.scatter(-azim1, repAoA_sph1, color='skyblue', label='retrieved',s=1.5)
    ax1.set_ylabel("Reported AoA (deg.)")
    ax1.set_xlabel("(a)")

    print(repAoA1 - repAoA_sph1)

    ax1.set_title('$0 < t \leq 15 $')

    print(azim2)
    ax2.scatter(-azim2, repAoA2, color='black',s=1.5)
    # ax2.scatter(-azim2, repAoA_init2, linewidth=0.5, color='blue',s=1.5)
    # ax2.scatter(-azim2, repAoA_ret2, color='green',s=1.5)
    ax2.scatter(-azim2, repAoA_sph2, color='skyblue', label='retrieved',s=1.5)

    ax2.set_xlabel("(b)")

    ax2.sharey(ax1)

    ax2.set_title('$15 < t \leq 30$')

    ax3.scatter(-azim3, repAoA3, color='black', s=1.5)
    # ax3.scatter(-azim3, repAoA_init3, linewidth=0.5, color='blue',s=1.5)
    # ax3.scatter(-azim3, repAoA_ret3, color='green',s=1.5)
    ax3.scatter(-azim3, repAoA_sph3, color='skyblue', label='retrieved',s=1.5)

    ax3.set_xlabel("Horizontal angle (deg.) \n (c)")
    ax4.set_xlabel("Horizontal angle (deg.) \n (d)")
    ax3.set_ylabel("Reported AoA (deg.)")

    ax1.sharex(ax3)

    ax3.set_title('$30 < t \leq 45$')

    ax4.scatter(-azim4, repAoA4, color='black', s=1.5)
    # ax4.scatter(-azim4, repAoA_init4, linewidth=0.5, color='blue',s=1.5)
    # ax4.scatter(-azim4, repAoA_ret4, color='green',s=1.5)
    ax4.scatter(-azim4, repAoA_sph4, color='skyblue', label='retrieved',s=1.5)

    ax4.set_title('$45 < t\leq 60$')

    fig.legend(bbox_to_anchor=(0.3, 0.68), markerscale=2)
    fig.tight_layout()
    fig.savefig("../Plots/PAPERII_retrieve_flightpaths.jpeg", dpi=700)

    # plt.scatter(-azim1, repAoA1, c=d1, edgecolors='none')
    # plt.ylabel("Reported AoA (deg.)")
    # plt.xlabel("Horizontal angle (deg.)")
    # cbar = plt.colorbar()
    # cbar.set_label("Distance (km)")
    
    plt.show()

if (typeplot == 46):

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
    error2, error2init = rep(hret2, d2) - repAoA2_sph, rep(hinit2, d2) - repAoA2_sph
    error3, error3init = rep(hret3, d3) - repAoA3_sph, rep(hinit3, d3) - repAoA3_sph
    error4, error4init = rep(hret4, d4) - repAoA4_sph, rep(hinit4, d4) - repAoA4_sph

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
    ax1.axhline(0.0, linestyle='--', color='black')
   
    ax1.set_ylabel("Reported AoA error (deg.)\n(Aircraft altitude ($10^2$ km))")
    ax1.set_xlabel("(a)")

    ax1.set_title('$0 < t \leq 15 $')

    ax2.plot(valsp, mean_err2, color='green')
    ax2.plot(valsp, mean_err_init2, color='black')
    ax2.axhline(0.0, linestyle='--', color='black')
    ax2.set_xlabel("(b)")

    ax2.set_title('$15 < t \leq 30$')

    ax3.plot(valsp, mean_err3, color='green')
    ax3.plot(valsp, mean_err_init3, color='black')
    ax3.axhline(0.0, linestyle='--', color='black')

    ax3.set_xlabel("Distance (km) \n (c)")
    ax4.set_xlabel("Distance (km) \n (d)")
    ax3.set_ylabel("Reported AoA error (deg.)")

    ax1.sharex(ax3)

    ax3.set_title('$30 < t \leq 45$')

    ax4.plot(valsp, mean_err4, color='green')
    ax4.plot(valsp, mean_err_init4, color='black')
    ax4.axhline(0.0, linestyle='--', color='black')
    ax4.set_title('$45 < t\leq 60$')

    fig.legend(bbox_to_anchor=(0.49, 0.9), markerscale=2)
    fig.tight_layout()

    fig.savefig("../Plots/PAPERII_retrieve_flightpaths_error.jpeg", dpi=700)
    
    plt.show()

if (typeplot == 1):

    fig, ax1 = plt.subplots(figsize=(8, 6))


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

    # fig.savefig("../Plots/PAPERII_retrieve_NE_RK3_Sep9.jpeg", dpi=700)

    plt.show()

if (typeplot == 2):

    fig, ax1 = plt.subplots(figsize=(6, 8))


    ax1.plot(100*(data[0]['init'] - data[0]['target'])/data[0]['target'], data[0]['h'], color='black', linestyle='--', label='initial')
    ax1.plot(100*(data[0]['retrieve'] - data[0]['target'])/data[0]['target'], data[0]['h'], color='green', label='0.00 RMS')
    ax1.plot(100*(data[1]['retrieve'] - data[0]['target'])/data[0]['target'], data[0]['h'], color='blue', label='{} RMS'.format(noises[1]))
    ax1.plot(100*(data[2]['retrieve'] - data[0]['target'])/data[0]['target'], data[0]['h'], color='red', label='{} RMS'.format(noises[2]))

    plt.axvline(0.0, linestyle=':', color='black')

    ax1.set_xlim(-20, 20)
    ax1.set_ylabel("Altitude (km)")
    ax1.set_xlabel(r'$\frac{N_{retrieved} - N_{radiosonde}}{N_{radiosonde}}$'+' (%)')

    plt.legend()

    fig.tight_layout()

    # fig.savefig("../Plots/PAPERII_retrieveError_NE_RK3_Sep9.jpeg", dpi=700)

    plt.show()

if (typeplot == 3):

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

if(typeplot==6):

    data[0] = pd.read_csv("../flightpaths/PAPERII_synobs_all.txt", sep=' ', header=None)
    data[0].columns = ["obsAoA", "h", "hm", "d", "t", "repAoA", "azim"]

    obsAoA = data[0]['obsAoA']
    repAoA = data[0]['repAoA']
    h = data[0]['h']
    hm = data[0]['hm']
    t = data[0]['t']

    print(len(h))

    d = data[0]['d']
    azim = data[0]['azim']

    repAoA_sph = np.arctan( ((r0+h)*np.cos(d/r0) - (r0+0.515)) / ((r0+h)*np.sin(d/r0))) * 180/np.pi
    repAoA_m = np.arctan( ((r0+hm)*np.cos(d/r0) - (r0+0.515)) / ((r0+hm)*np.sin(d/r0))) * 180/np.pi

    fig, ax1 = plt.subplots(figsize=(8, 6))

    o_m = repAoA - repAoA_m

    ax1.scatter(repAoA-repAoA_m, obsAoA, s=1.5,c=t, edgecolors="none")
    ax1.set_ylabel("Observed AoA (deg.)")
    ax1.set_xlabel("(O - M) refracted angle (deg.)")
    ax1.set_xlim(-0.5, 0.5)
    plt.axvline(0, linestyle='--')
    plt.savefig("../plots/test.jpeg")

