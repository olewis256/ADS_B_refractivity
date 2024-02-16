import numpy as np
import matplotlib.pyplot as plt
import smplotlib
import sys
import pandas as pd

r0 = 6364.539

typeplot = int(sys.argv[1])


data = [0]*3
noises = [0, 0.01, 0.03]

if (typeplot == 1) or (typeplot == 2):
    for i in range(3):
        
        data[i] = pd.read_csv("../retrievals/PAPERII_retrieve_NE_RK3_NEW_{}.txt".format(noises[i]), sep=' ', header=None)
        data[i].columns = ["retrieve", "h", "target", "init"]

if (typeplot == 3):
    for i in range(3):
        
        data[i] = pd.read_csv("../RMS/PAPERII_RMS_NE_RK3_{}.txt".format(noises[i]), sep=' ', header=None)
        data[i].columns = ["iter", "loss", "Nrms"]


if (typeplot == 4):

    data[0] = pd.read_csv("../retrievals/PAPERII_retrieve_NE_TRUE0.txt", sep=' ', header=None)
    data[0].columns = ["retrieve", "h", "target", "init"]

    fig, ax1 = plt.subplots(figsize=(8, 6))


    ax1.plot(data[0]['target'], data[0]['h'], color='black', label='radiosonde')
    ax1.plot(data[0]['init'], data[0]['h'], linewidth=0.5, linestyle='--', color='black', label='initial')
    ax1.plot(data[0]['retrieve'], data[0]['h'], color='green', label='retrieved (real)')

    plt.legend()

    ax1.set_ylabel("Altitude (km)")
    ax1.set_xlabel("Refractivity (ppm)")
    
    plt.show()

if (typeplot == 45):

    data[0] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths0.txt", sep=' ', header=None)
    data[0].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]

    obsAoA = data[0]['obsAoA']
    repAoA = data[0]['repAoA']
    h = data[0]['h']
    d = data[0]['d']
    hinit = data[0]['hinit']
    hret = data[0]['hret']
    azim = data[0]['azim']

    fig, ax1 = plt.subplots(figsize=(8, 6))

    repAoA_sph = np.arctan( ((r0+h)*np.cos(d/r0) - (r0+0.515)) / ((r0+h)*np.sin(d/r0))) * 180/np.pi
    repAoA_init = np.arctan( ((r0+hinit)*np.cos(d/r0) - (r0+0.515)) / ((r0+hinit)*np.sin(d/r0))) * 180/np.pi
    repAoA_ret = np.arctan( ((r0+hret)*np.cos(d/r0) - (r0+0.515)) / ((r0+hret)*np.sin(d/r0))) * 180/np.pi


    ax1.scatter(-azim, repAoA, color='black', label='obs', s=1.5)
    ax1.scatter(-azim, repAoA_init, linewidth=0.5, color='blue', label='initial',s=1.5)
    ax1.scatter(-azim, repAoA_ret, color='green', label='retrieved',s=1.5)

    ax1.set_ylabel("Reported AoA (deg.)")
    ax1.set_xlabel("Horizontal angle (deg.)")

    plt.legend()

    fig.savefig("../Plots/PAPERII_retrieve_flightpath_t900_1800.jpeg", dpi=700)
    
    plt.show()

if (typeplot == 46):

    data[0] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths0.txt", sep=' ', header=None)
    data[0].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]

    obsAoA = data[0]['obsAoA']
    repAoA = data[0]['repAoA']
    h = data[0]['h']
    d = data[0]['d']
    hinit = data[0]['hinit']
    hret = data[0]['hret']
    azim = data[0]['azim']

    fig, ax1 = plt.subplots(figsize=(8, 6))

    # repAoA = np.arctan( ((r0+h)*np.cos(d/r0) - (r0+0.577)) / ((r0+h)*np.sin(d/r0))) * 180/np.pi
    repAoA_init = np.arctan( ((r0+hinit)*np.cos(d/r0) - (r0+0.577)) / ((r0+hinit)*np.sin(d/r0))) * 180/np.pi
    repAoA_ret = np.arctan( ((r0+hret)*np.cos(d/r0) - (r0+0.577)) / ((r0+hret)*np.sin(d/r0))) * 180/np.pi

    error_init = repAoA_init - repAoA
    error_retrieve = repAoA_ret - repAoA

    ax1.scatter(d, abs(error_init), linewidth=0.5, color='blue', label='init',s=1.5)
    ax1.scatter(d, abs(error_retrieve), color='green', label='retrieved', s=1.5)

    ax1.set_ylabel("Absolute reported AoA error (deg.)")
    ax1.set_xlabel("Distance (km)")

    plt.legend()

    fig.savefig("../Plots/PAPERII_retrieve_error_distance_t900_1800.jpeg", dpi=700)
    
    plt.show()

    # data[0] = pd.read_csv("../flightpaths/PAPERII_retrieve_paths0.txt", sep=' ', header=None)
    # data[0].columns = ["obsAoA", "h", "d", "hinit", "hret", "repAoA", "azim"]

    # obsAoA = data[0]['obsAoA']
    # repAoA = data[0]['repAoA']
    # h = data[0]['h']
    # d = data[0]['d']
    # hinit = data[0]['hinit']
    # hret = data[0]['hret']
    # azim = data[0]['azim']

    # fig, ax1 = plt.subplots(figsize=(8, 6))

    # # repAoA = np.arctan( ((r0+h)*np.cos(d/r0) - (r0+0.577)) / ((r0+h)*np.sin(d/r0))) * 180/np.pi
    # repAoA_init = np.arctan( ((r0+hinit)*np.cos(d/r0) - (r0+0.515)) / ((r0+hinit)*np.sin(d/r0))) * 180/np.pi
    # repAoA_ret = np.arctan( ((r0+hret)*np.cos(d/r0) - (r0+0.515)) / ((r0+hret)*np.sin(d/r0))) * 180/np.pi

    # error_init = abs(repAoA_init - repAoA)
    # error_retrieve = abs(repAoA_ret - repAoA)

    # vals = np.linspace(0, 450, 20)
    # valsp = np.linspace(0, 450, 19)
    # mean_err = []
    # mean_err_init = []

    # for i in range(len(vals)-1):
    #     if np.isnan(np.mean(error_retrieve[(d >= vals[i]) & (d < vals[i+1])])):
    #         mean_err.append(0)
    #         mean_err_init.append(0)
    #     else: 
    #         mean_err.append(np.mean(error_retrieve[(d >= vals[i]) & (d < vals[i+1])]))
    #         mean_err_init.append(np.mean(error_init[(d >= vals[i]) & (d < vals[i+1])]))

    # print(mean_err)
    # ax1.plot(valsp, mean_err_init, linewidth=0.5, color='blue', label='init')
    # ax1.plot(valsp, mean_err, color='green', label='retrieved')

    # ax1.set_ylabel("Mean absolute repAoA error (deg.)")
    # ax1.set_xlabel("Distance (km)")

    # plt.legend()

    # fig.savefig("../Plots/PAPERII_retrieve_error_distance.jpeg", dpi=700)
    
    # plt.show()

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

    fig.savefig("../Plots/PAPERII_retrieveError_NE_RK3_Sep9.jpeg", dpi=700)

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

