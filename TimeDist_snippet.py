import uproot
import numpy as np
import random
# Open the ROOT file
with uproot.open("electron_2MeV.root") as fMyFile:
    # Access the TTree
    myTree = fMyFile["Hits"]
    # Define variables and constants
    posZ = 0.0  # translation to add “by hand”
    fibre_dt = 2.8  # decay time WLS fibre
    myEff = 100 * (0.108 * 0.4 * 0.9)  # in %
    TTS = 0.4  # fibre transit time spread / metre
    timeFibre = 6.26  # average time, in ns, for photon to travel 1m of fibre
    Att_leng = 5.0  # fibre attenuation length metres
    # Set the seed for the random number generator
    mySeed = 2
    random.seed(mySeed)
    # Loop over the events
    for data in myTree.iterate(["Hit_Z", "Time_ns"], cut="Event_Number==0", library="np", step_size="1GB"):
        Hit_Z = data["Hit_Z"]
        Time_ns = data["Time_ns"]
        nEvt = len(Hit_Z)
        for i in range(nEvt):
            detected = random.random()  # to choose if this hit makes a signal
            side = random.random()  # to choose if this hit goes to detector’s back or front side
            if side <= 0.5:  # front electronics
                distTravel = 2.0 - posZ - Hit_Z[i] / 1e3  # distance travelled from hit position to SiPM
                totalEff = myEff * np.exp(-distTravel / Att_leng)
                if detected < totalEff / 100.0:  # we have a hit!
                    WLStime = random.expovariate(1 / fibre_dt)  # decay time fibre
                    spreadT = (random.random() - 0.5) * distTravel * TTS  # time spread
                    hitTime = Time_ns[i] + WLStime + spreadT + distTravel * timeFibre
                    print(hitTime)
