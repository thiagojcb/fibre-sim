#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#Created on Fri Feb  10 17:30:58 2023
#
#@author: Thiago Bezerra
#
#Analytical calculation of a WLS fibre trapping efficiency

import math
import numpy as np
import matplotlib.pyplot as plt

def get_trap_eff(n_exterior):

    core_eff      = 0
    core_clad_eff = 0
    clad_clad_eff = 0    

    #critical angle for second clad to exterior interface
    critic_ang3 = math.asin( n_exterior / n_clad_out )
    core_eff = 0.5 * (-math.cos(PI/2 - critic_ang1) + math.cos(0))

    #refraction core-to-clad angle of clad-to-clad crictical angle
    new_angle1 = math.asin( n_clad_in * math.sin(critic_ang2) / n_core)
    core_clad_eff = 0.5 * (-math.cos(PI/2 - new_angle1) + math.cos(0))

    #refraction clad-to-clad angle of clad-to-exterior crictical angle
    new_angle2 = math.asin( n_clad_out * math.sin(critic_ang3) / n_core )    
    clad_clad_eff = 0.5 * (-math.cos(PI/2 - new_angle2) + math.cos(0))

    #print(n_exterior,(critic_ang1)*180/PI,(new_angle1)*180/PI,(new_angle2)*180/PI)
    
    return core_eff, core_clad_eff, clad_clad_eff

def get_trap_eff_skew(n_exterior):

    core_eff      = 0
    core_clad_eff = 0
    clad_clad_eff = 0    

    #critical angle for second clad to exterior interface
    critic_ang3 = math.asin( n_exterior / n_clad_out )
    core_eff = (0.5 * (-math.cos(PI/2 - critic_ang1) + math.cos(0))
                    * math.cos(PI/2 - critic_ang1))

    #refraction core-to-clad angle of clad-to-clad crictical angle
    new_angle1 = math.asin( n_clad_in * math.sin(critic_ang2) / n_core)
    core_clad_eff = (0.5 * (-math.cos(PI/2 - new_angle1) + math.cos(0))
                         * math.cos(PI/2 - new_angle1))

    #refraction clad-to-clad angle of clad-to-exterior crictical angle
    new_angle2 = math.asin( n_clad_out * math.sin(critic_ang3) / n_core )    
    clad_clad_eff = (0.5 * (-math.cos(PI/2 - new_angle2) + math.cos(0))
                         * math.cos(PI/2 - new_angle2))

    #print(n_exterior,(critic_ang1)*180/PI,(new_angle1)*180/PI,(new_angle2)*180/PI)
    
    return core_eff, core_clad_eff, clad_clad_eff
            
def scan_ref_index(min_n, max_n, width_n):
    loss_i  = list()
    core_i  = list()
    clad1_i = list()
    clad2_i = list()
    n_i     = list()

    #Scan over a range of external index values

    for ni in np.arange(min_n,max_n,width_n):

        core_eff, core_clad_eff, clad_clad_eff \
            = get_trap_eff(ni)
        loss = 1 - clad_clad_eff
        loss_i.append(loss*100)
        clad1_i.append(core_clad_eff*100)
        clad2_i.append(clad_clad_eff*100)
        core_i.append(core_eff*100)
        n_i.append(ni)
        #print(ni,core_eff,core_clad_eff,clad_clad_eff,loss)
        print(ni,core_clad_eff/core_eff,clad_clad_eff/core_eff,
              clad_clad_eff/core_clad_eff)

    plt.xlabel('External environment refractive index')
#    plt.plot(n_i,loss_i,'r-',label='loss')
    plt.plot(n_i,core_i,'b-',label='core')
    plt.plot(n_i,clad1_i,'g--',label='1st clad')
    plt.plot(n_i,clad2_i,'r-.',label='2nd clad')
    plt.ylabel('Trapping efficiency (%)')
    plt.grid(linestyle='--', linewidth=0.8, alpha=0.5)
    plt.ylim([0, 20])
    plt.legend(loc='best')    
    plt.show()    

##### MAIN

PI = math.pi

#SG fibre values
n_core     = 1.59
n_clad_in  = 1.49
n_clad_out = 1.42
n_exterior = 1.00

#critical angle for core-first clad interface
critic_ang1 = math.asin( n_clad_in  / n_core     )
#critical angle for first to sencond clad interface
critic_ang2 = math.asin( n_clad_out / n_clad_in  )

#print(get_trap_eff(1.00))

scan_ref_index(1.00,1.42,0.01)





