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

PI = math.pi

#SG fibre values
n_core     = 1.60
n_clad_in  = 1.49
n_clad_out = 1.42
n_exterior = 1.00

critic_ang1 = math.asin( n_clad_in  / n_core     ) #critical angle for core-first clad interface
critic_ang2 = math.asin( n_clad_out / n_clad_in  ) #critical angle for first to sencond clad interface

nevts = 1000000

loss_i = list()
clad_i = list()
n_i    = list()

for ni in np.arange(1.00,1.42,0.01):
    n_exterior = ni
    critic_ang3 = math.asin( n_exterior / n_clad_out ) #critical angle for second clad to exterior interface

    core_evt      = 0
    core_clad_evt = 0
    clad_clad_evt = 0
    loss_evt      = 0

#    if ni == 1.00 or (ni > 1.33 and ni<1.34):
#        fileName = "Fibre_toy_MC_"+str(ni)+".txt"
#        f = open(fileName,"w")


    for i in range(0,nevts):

        is_core       = 0
        is_clad_in    = 0
        is_clad_out   = 0
        is_loss       = 0
        new_angle1    = 0
        new_angle2    = 0        

        angle = math.pi * random.random() / 2 #light emission angle, wrt surface (0 ~ 90 degree)

        if angle < critic_ang1:

            new_angle1 = math.asin( n_core * math.sin(angle) / n_clad_in) #refraction core-to-clad angle

            if new_angle1 < critic_ang2:

                new_angle2 = math.asin( n_clad_in * math.sin(new_angle1) / n_clad_out) #refraction clad-to-clad angle

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
            
        #if ni == 1.00 or (ni > 1.33 and ni<1.34):
        #    f.write(str(angle*180/PI)+" "+str(new_angle1*180/PI)+" "+str(new_angle2*180/PI)+" "+str(is_core)+" "+str(is_clad_in)+" "+str(is_clad_out)+" "+str(is_loss)+"\n")
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






