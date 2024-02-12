import numpy as np
import matplotlib.pyplot as plt
import smplotlib
import sys
import pandas as pd

typeplot = int(sys.argv[1])


data = [0]*3
noises = [0, 0.01, 0.03]

if (typeplot == 1) or (typeplot == 2):
    for i in range(3):
        
        data[i] = pd.read_csv("../Retrievals/PAPERII_retrieve_NE_RK3_NEW_{}.txt".format(noises[i]), sep=' ', header=None)
        data[i].columns = ["retrieve", "h", "target", "init"]

if (typeplot == 3):
    for i in range(3):
        
        data[i] = pd.read_csv("../RMS/PAPERII_RMS_NE_RK3_{}.txt".format(noises[i]), sep=' ', header=None)
        data[i].columns = ["iter", "loss", "Nrms"]

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

    fig.savefig("../Plots/PAPERII_retrieve_NE_RK3.jpeg", dpi=700)

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

    fig.savefig("../Plots/PAPERII_retrieveError_NE_RK3.jpeg", dpi=700)

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
