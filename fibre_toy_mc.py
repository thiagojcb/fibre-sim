#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#Created on Tue Feb  7 15:03:58 2023
#
#@author: Thiago Bezerra
#
#Simple toy MC to get the number of light propagated in a WLS fibre

import random
import math
import numpy as np
import matplotlib.pyplot as plt

def propagate_rdm_ligth(n_exterior, nevt = 1000, save_output = False):

    #critical angle for second clad to exterior interface
    critic_ang3 = math.asin( n_exterior / n_clad_out )

    core_evt      = 0
    core_clad_evt = 0
    clad_clad_evt = 0
    loss_evt      = 0

    if save_output:
        fileName = "data/Fibre_toy_MC_"+str(n_exterior)+".txt"
        f = open(fileName,"w")


    for i in range(0,nevts):

        is_core       = 0
        is_clad_in    = 0
        is_clad_out   = 0
        is_loss       = 0
        new_angle1    = 0
        new_angle2    = 0        

        #light emission angle, wrt surface (0 ~ 90 degree)
        angle = math.pi * random.random() / 2

        if angle < critic_ang1:
            #refraction core-to-clad angle
            new_angle1 = math.asin( n_core * math.sin(angle) / n_clad_in) 

            if new_angle1 < critic_ang2:
                #refraction clad-to-clad angle
                new_angle2 = math.asin(
                    n_clad_in * math.sin(new_angle1)
                    / n_clad_out )

                if new_angle2 < critic_ang3:

                    loss_evt += 1
                    is_loss = 1

                else:
                    clad_clad_evt += 1
                    is_clad_out = 1

            else:
                core_clad_evt += 1
                is_clad_in = 1

        else:
            core_evt += 1
            is_core = 1

        if save_output:
            f.write(str(angle*180/PI)+" "
                    +str(new_angle1*180/PI)+" "
                    +str(is_core)+" "
                    +str(new_angle2*180/PI)+" "
                    +str(is_clad_in)+" "
                    +str(is_clad_out)+" "
                    +str(is_loss)+"\n")
            
    return core_evt, core_clad_evt, clad_clad_evt, loss_evt            
            
def scan_ref_index(events):
    loss_i = list()
    clad_i = list()
    n_i    = list()

    #Scan over a range of external index values
    #From air up to outer cladding
    min_n   = 1.00
    max_n   = n_clad_out
    n_width = 0.01

    for ni in np.arange(min_n,max_n,n_width):

        core_evt, core_clad_evt, clad_clad_evt, loss_evt \
            = propagate_rdm_ligth(ni,events)

        #print("core light intensity",100*core_evt/nevts,"%")
        #print("clad-in light intensity",100*core_clad_evt/nevts,"%")
        #print("clad-out light intensity",100*clad_clad_evt/nevts,"%")
        #print("light loss",100*loss_evt/nevts,"%")

        loss_i.append(loss_evt/nevts)
        clad_i.append(clad_clad_evt/nevts)
        n_i.append(ni)
        print(ni,loss_evt/nevts)

    plt.xlabel('external material refractive index')
    plt.plot(n_i,loss_i,'r-',label='loss')
    plt.plot(n_i,clad_i,'b-',label='outter clad reflection')
    plt.ylabel('fraction of light')
    plt.grid(linestyle='--', linewidth=0.8, alpha=0.5)
    plt.legend(loc='best')
    plt.show()    

##### MAIN

PI = math.pi

#SG fibre values
n_core     = 1.60
n_clad_in  = 1.49
n_clad_out = 1.42
n_exterior = 1.00

#critical angle for core-first clad interface
critic_ang1 = math.asin( n_clad_in  / n_core     )
#critical angle for first to sencond clad interface
critic_ang2 = math.asin( n_clad_out / n_clad_in  )

nevts = 1000000

scan_ref_index(nevts)





