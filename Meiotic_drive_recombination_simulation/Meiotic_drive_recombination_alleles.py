# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 15:45:38 2015

@author: Lucas hemmer

This code was used a final project for a Biological Simulations Course
at the University of Kansas in 2015 as co-taught by John K. Kelly and Michael Tourtellot.

It is a recreation of the results and conclusions of an academic paper:
    "Scrambling eggs: Meiotic drive and the evolution of female recombination rates"
    originally published in 2012 in Genetics by Yaniv Brandvain and Graham Coop. 
    The link for the original paper: doi:10.1534/genetics.111.136721

The results of the paper suggest alleles involved in meiotic drive, or non-Mendelian segregation
of chromosomes, can have an effect on recombination rates in females buy separating modifiers of
meiotic drive away from driver alleles. 

The original code was in Mathematica but I recreated the same theoritial results in Python.
The script was originally in Python 2.7 but I shifted it over to Python 3 since handing in the assignment.

The code simulates populations and the frequency of driving alleles in populations over time

While I created the script, please cite the original paper for the results and conclusions.
The recreation of the figures from the paper are noted in the script.

"""


#IMPORT Libraries -------------------------------------------------------

import math
import matplotlib.pyplot as plt

#INITIALIZE PARAMETERS --------------------------------------------------

Gens = 10000    # generations
r = 0.25       #frequency of recombination, between 0-0.5
D = 0.1            #freq of driver allele, between 0-1
M = 0.01     #freq of recombination modifier, additive, between 0-1
deltaR = 0.05   #increase in recombination for each Modifier allele 

alphaM1 = 1.0  #strength of distortion of first meiotic division

wDD = 0.0 # fitness of homozygous driver 
wDd = 1.0 # fitness of heterozygote
wdd = 1.0 # fitness in absence of driver

## ------------------------------------
def driveStrength(r, alphaM1):
    Xm1 = alphaM1*(1-r) + (r/2)
    return Xm1

Xm1 = driveStrength(r, alphaM1) #strength of drive in ith meiotic division, 1st division
#Xm2 = alphaM2*r + ((1-r)/2)  #second division

## ------------------------------------   
def Graph1locus(Title, Driver, driveLabel):

    plt.title(Title)
    plt.xlabel('Generations')
    plt.ylabel('Frequency')
    plt.axis([0, (len(Driver)+1)/2, 0.0, 1.0])
    plt.plot(Driver, label = driveLabel)
    plt.legend()

    plt.show() 

## ------------------------------------   
def GraphwithMod(Title, Driver, driveLabel, Modifier, modLabel):

    plt.title(Title)
    plt.xlabel('Generations')
    plt.ylabel('Frequency')
    plt.axis([0, (len(Driver)+1)/2, 0.0, 1.0])
    plt.plot(Driver, label = driveLabel)
    plt.plot(Modifier, label = modLabel)
    plt.legend()

    plt.show() 
 
    
#SINGLE LOCUS DRIVE -----------------------------------------------------
#equivalent drive fitness in males and females 

## ------------------------------------
## drive is equal between sexes

onelocus_equaldrive = []
onelocus_equaldrive.append(D)

for i in range(Gens):
    
    if i == 0:
        Dtemp = D
        #print i, Dtemp
    
    #genotype = [DD, Dd, dd]    
    g = [Dtemp*Dtemp, 2*Dtemp*(1-Dtemp), (1-Dtemp)*(1-Dtemp)]  #genotype freq
    meanW = g[2]*wdd + g[1]*wDd + g[0]*wDD   #mean population fitness
    Dtemp = (g[0]*wDD + Xm1*g[1]*wDd) / meanW  #new freq of D after selection
    onelocus_equaldrive.append(Dtemp) #keeping track
    #print i+1, Dtemp
    
    if onelocus_equaldrive[-1] - onelocus_equaldrive[-2] < 10**-10:
        break

#Not illustrated, check with predicted equilibrium equations
Graph1locus("1 Locus Meiotic Drive", onelocus_equaldrive, "Drive freq")
    
drive_eq = (-1 + r*wDd + 2*alphaM1*wDd - 2*r*alphaM1*wDd) / (-1 + 2*wDd - wDD) 
print("prediction =", drive_eq)

## ------------------------------------
## drive is female limited

onelocus_femaledrive = []
onelocus_femaledriver_male = []
onelocus_femaledriver_female = [] 
onelocus_femaledrive.append(D)
onelocus_femaledriver_male.append(D)
onelocus_femaledriver_female.append(D)

for i in range(Gens):
    
    if i == 0:
        Dtemp = male_D = female_D = D        
        #print i, Dtemp, male_D, female_D
    
    #genotype = [DD, Dd, dd]
    g = [male_D*female_D, male_D*(1-female_D) + female_D*(1-male_D), (1-male_D)*(1-female_D)]
    
    #mean population fitness 
    meanW = g[2]*wdd + g[1]*wDd + g[0]*wDD
    
    #number of new gametes from males and females and total after drive and selection
    male_D = (wDD*g[0] + ((wDd*g[1]) / 2)) / meanW
    female_D = (wDD*g[0] + wDd*g[1]*Xm1) / meanW  
    Dtemp = (male_D/2) + (female_D/2)
    
    #keeping track
    onelocus_femaledrive.append(Dtemp)
    onelocus_femaledriver_male.append(male_D)   
    onelocus_femaledriver_female.append(female_D)   
    #print i+1, Dtemp, male_D, female_D
    
    if onelocus_femaledrive[-1] - onelocus_femaledrive[-2] < 10**-10:
        break

#Not illustrated, check with predicted equilibrium equations 
Graph1locus("Female-Limited Meiotic Drive", onelocus_femaledrive, "Drive freq")

#predicted equilibrium values
maleDrive_eq = (1 + 2*female_D - math.sqrt((-1-2*female_D)**2 - 8*female_D**2)) / (4*female_D)      #D is actually male D freq, but equal to total D freq intitially
print("male prediction =", maleDrive_eq)

femaleDrive_eq = 0.5 * (2 - r - math.sqrt(2*r - r**2)) 
print("female prediction =", femaleDrive_eq)

print("driver freq prediction =", maleDrive_eq*0.5 + femaleDrive_eq*0.5)

#SINGLE LOCUS DRIVE - UNLINKED MODIFIER ---------------------------------

## ------------------------------------
## drive is equal between sexes

onelocus_modDrive = []
onelocus_Modifier = []
onelocus_modDrive.append(D)
onelocus_Modifier.append(M)
Xm1_mm = driveStrength(r, alphaM1)
Xm1_Mm = driveStrength(r + deltaR, alphaM1)
Xm1_MM = driveStrength(r + 2*deltaR, alphaM1)
    
for i in range(Gens):
    
    #intialize first generation    
    if i == 0:
        Dtemp = D
        Mtemp = M
        DM = Dtemp * Mtemp
        Dm = Dtemp * (1-Mtemp)
        dM = (1-Dtemp) * (Mtemp) 
    
    #genotype = [DDMM, DDMm, DDmm, DdMM, DdMm, Ddmm, ddMM, ddMm, ddmm]    
    g = [DM*DM, 2*DM*Dm, Dm*Dm, 2*DM*dM, 2*DM*(1-DM-Dm-dM)+2*Dm*dM, 2*Dm*(1-DM-Dm-dM), 
         dM*dM, 2*dM*(1-DM-Dm-dM), (1-DM-Dm-dM)*(1-DM-Dm-dM)]  #genotype freq
    
    #mean population fitness     
    meanW = (g[0]+g[1]+g[2])*wDD + (g[3]+g[4]+g[5])*wDd + (g[6]+g[7]+g[8])*wdd
    
    #gamete haplotype frequencies, fitness/selection and drive taken into account
    DM = ((g[0] + g[1]*0.5)*wDD + (g[3]*Xm1_MM + g[4]*Xm1_Mm*0.5)*wDd) / meanW
    Dm = ((g[1]*0.5 + g[2])*wDD + (g[4]*Xm1_Mm*0.5 + g[5]*Xm1_mm)*wDd) / meanW
    dM = ((g[3]*(1-Xm1_MM) + g[4]*(1-Xm1_Mm)*0.5)*wDd + (g[6] + g[7]*0.5)*wdd) / meanW
    dm = ((g[4]*(1-Xm1_Mm)*0.5 + g[5]*(1-Xm1_mm))*wDd + (g[7]*0.5 + g[8])*wdd) / meanW
    
    #LD produced by driver and modifier 
    LD = DM*(1 - DM - Dm - dM) - Dm*dM
    #print i, "D freq", DM + Dm, "M freq", DM+dM, LD 
    
    #keeping track    
    onelocus_modDrive.append(DM+Dm)
    onelocus_Modifier.append(DM+dM)
    
    if onelocus_modDrive[-1] - onelocus_modDrive[-2] < 10**-10 and onelocus_Modifier[-1] - onelocus_Modifier[-2] < 10**-10:
        break

#Figure 3A
GraphwithMod("Meiotic Drive with Unlinked Modifier", onelocus_modDrive, "Driver freq", onelocus_Modifier, "Modifier freq")

## ------------------------------------
## drive is female limited

onelocus_femMod = []
#onelocus_femMod_male = []
#onelocus_femMod_female = [] 
onelocus_femMod.append(D)
#onelocus_femMod_male.append(D)
#onelocus_femMod_female.append(D)
onelocus_femModifier = []
onelocus_femModifier.append(M)

for i in range(Gens):
    if i == 0:
        Dtemp = male_D = female_D = D  
        Mtemp = male_M = female_M = M
        DM = Dtemp * Mtemp
        Dm = Dtemp * (1-Mtemp)
        dM = (1-Dtemp) * (Mtemp)

    #genotype = [DDMM, DDMm, DDmm, DdMM, DdMm, Ddmm, ddMM, ddMm, ddmm] 
    g = [DM*DM, 2*DM*Dm, Dm*Dm, 2*DM*dM, 2*DM*(1-DM-Dm-dM)+2*Dm*dM, 2*Dm*(1-DM-Dm-dM), 
         dM*dM, 2*dM*(1-DM-Dm-dM), (1-DM-Dm-dM)*(1-DM-Dm-dM)]

    #mean population fitness
    meanW = (g[0]+g[1]+g[2])*wDD + (g[3]+g[4]+g[5])*wDd + (g[6]+g[7]+g[8])*wdd

    #male gamete freq after fitness/selection and drive
    DM_male = ((g[0] + g[1]*0.5)*wDD + (g[3]*0.5 + g[4]*0.25)*wDd) / meanW
    Dm_male = ((g[1]*0.5 + g[2])*wDD + (g[4]*0.25 + g[5]*0.5)*wDd) / meanW
    dM_male = ((g[3]*0.5 + g[4]*0.25)*wDd + (g[6] + g[7]*0.5)*wdd) / meanW
    dm_male = ((g[4]*0.25 + g[5]*0.5)*wDd + (g[7]*0.5 + g[8])*wdd) / meanW
    #print DM_male + Dm_male + dM_male + dm_male

    #female gamete freq after fitness/selection and drive
    DM_female = ((g[0] + g[1]*0.5)*wDD + (g[3]*Xm1_MM + g[4]*Xm1_Mm*0.5)*wDd) / meanW
    Dm_female = ((g[1]*0.5 + g[2])*wDD + (g[4]*Xm1_Mm*0.5 + g[5]*Xm1_mm)*wDd) / meanW
    dM_female = ((g[3]*(1-Xm1_MM) + g[4]*(1-Xm1_Mm)*0.5)*wDd + (g[6] + g[7]*0.5)*wdd) / meanW
    dm_female = ((g[4]*(1-Xm1_Mm)*0.5 + g[5]*(1-Xm1_mm))*wDd + (g[7]*0.5 + g[8])*wdd) / meanW
    #print DM_female + Dm_female + dM_female + dm_female

    #gamete freq
    DM = DM_male*0.5 + DM_female*0.5
    Dm = Dm_male*0.5 + Dm_female*0.5
    dM = dM_male*0.5 + dM_female*0.5
    dm = dm_male*0.5 + dm_female*0.5
    #print DM + Dm + dM + dm
    
    onelocus_femMod.append(DM+Dm)
    #onelocus_femMod_male.append(DM_male+Dm_male)
    #onelocus_femMod_female.append(DM_female+Dm_female) 
    onelocus_femModifier.append(DM+dM)
    
    #LD_female = -0.5*(-female_D + male_D*(-1+2*female_D))*(-1+(DM+dM))*(DM+dM)*wDd*(-1+2*alphaM1)*deltaR
    #print "LD pred", LD_female

    #print i, "D", (DM+Dm), "M", (DM+dM), "LD", LD_female

    if onelocus_femMod[-1] - onelocus_femMod[-2] < 10**-10 and onelocus_femModifier[-1] - onelocus_femModifier[-2] < 10**-10:
        break

#Figure 3B
GraphwithMod("Female-Limited Drive with Unlinked Modifier", onelocus_femMod, "Driver freq", onelocus_femModifier, "Modifier freq")

#SINGLE LOCUS DRIVE - LINKED MODIFIER -----------------------------------

## ------------------------------------
## drive is in both sexes
## recombination modifier linked to non-driver haplotype


onelocus_linkedD1 = []
onelocus_linkedMod1 = []
onelocus_linkedD1.append(D)
onelocus_linkedMod1.append(M)

deltaR = 0.05
Xm1_mm = driveStrength(r, alphaM1)
Xm1_Mm = driveStrength(r + deltaR, alphaM1)
Xm1_MM = driveStrength(r + 2*deltaR, alphaM1)
    
for i in range(Gens):
    
    #intialize first generation    
    if i == 0:
        Dtemp = D
        Mtemp = M
        DM = Dtemp * Mtemp
        Dm = Dtemp * (1-Mtemp)
        dM = (1-Dtemp) * (Mtemp) 
    
    #genotype = [DDMM, DDMm, DDmm, DdMM, DdMm, Ddmm, ddMM, ddMm, ddmm]    
    g = [DM*DM, 2*DM*Dm, Dm*Dm, 2*DM*dM, 2*DM*(1-DM-Dm-dM)+2*Dm*dM, 2*Dm*(1-DM-Dm-dM), 
         dM*dM, 2*dM*(1-DM-Dm-dM), (1-DM-Dm-dM)*(1-DM-Dm-dM)]  #genotype freq
    
    #mean population fitness     
    meanW = (g[0]+g[1]+g[2])*wDD + (g[3]+g[4]+g[5])*wDd + (g[6]+g[7]+g[8])*wdd
    
    #gamete haplotype frequencies, fitness/selection and drive taken into account
    DM = ((g[0] + g[1]*0.5)*wDD + (g[3]*Xm1_MM)*wDd) / meanW
    Dm = ((g[1]*0.5 + g[2])*wDD + (g[4]*Xm1_Mm + g[5]*Xm1_mm)*wDd) / meanW
    dM = ((g[3]*(1-Xm1_MM) + g[4]*(1-Xm1_Mm))*wDd + (g[6] + g[7]*0.5)*wdd) / meanW
    dm = ((g[5]*(1-Xm1_mm))*wDd + (g[7]*0.5 + g[8])*wdd) / meanW
    
    #LD produced by driver and modifier 
    LD = DM*(1 - DM - Dm - dM) - Dm*dM
    #print i, "D freq", DM + Dm, "M freq", DM+dM, LD 
    #print DM+Dm+dM+dm
    
    #keeping track    
    onelocus_linkedD1.append(DM+Dm)
    onelocus_linkedMod1.append(DM+dM)
    
    if onelocus_linkedD1[-1] - onelocus_linkedD1[-2] < 10**-10 and onelocus_linkedMod1[-1] - onelocus_linkedMod1[-2] < 10**-10:
        break

#Figure 4A
GraphwithMod("Meiotic Drive, Modifier Linked to Non-Driver", onelocus_linkedD1, "Driver freq", onelocus_linkedMod1, "Modifier freq")
    
## ------------------------------------
## drive is female limited
## recombination modifier linked to non-driver haplotype

onelocus_femlinkD1 = []
onelocus_femlinkMod1 = []
onelocus_femlinkD1.append(D)
onelocus_femlinkMod1.append(M)

deltaR = 0.05
Xm1_mm = driveStrength(r, alphaM1)
Xm1_Mm = driveStrength(r + deltaR, alphaM1)
Xm1_MM = driveStrength(r + 2*deltaR, alphaM1)
    
for i in range(Gens):
    
    #intialize first generation    
    if i == 0:
        Dtemp = D
        Mtemp = M
        DM = Dtemp * Mtemp
        Dm = Dtemp * (1-Mtemp)
        dM = (1-Dtemp) * (Mtemp) 
    
    #genotype = [DDMM, DDMm, DDmm, DdMM, DdMm, Ddmm, ddMM, ddMm, ddmm]    
    g = [DM*DM, 2*DM*Dm, Dm*Dm, 2*DM*dM, 2*DM*(1-DM-Dm-dM)+2*Dm*dM, 2*Dm*(1-DM-Dm-dM), 
         dM*dM, 2*dM*(1-DM-Dm-dM), (1-DM-Dm-dM)*(1-DM-Dm-dM)]  #genotype freq
    
    #mean population fitness     
    meanW = (g[0]+g[1]+g[2])*wDD + (g[3]+g[4]+g[5])*wDd + (g[6]+g[7]+g[8])*wdd
    
    DM_male = ((g[0] + g[1]*0.5)*wDD + (g[3]*0.5)*wDd) / meanW
    Dm_male = ((g[1]*0.5 + g[2])*wDD + (g[4]*0.5 + g[5]*0.5)*wDd) / meanW
    dM_male = ((g[3]*0.5 + g[4]*0.5)*wDd + (g[6] + g[7]*0.5)*wdd) / meanW
    dm_male = ((g[5]*0.5)*wDd + (g[7]*0.5 + g[8])*wdd) / meanW
    #print DM_male + Dm_male + dM_male + dm_male

    #female gamete freq after fitness/selection and drive
    DM_female = ((g[0] + g[1]*0.5)*wDD + (g[3]*Xm1_MM)*wDd) / meanW
    Dm_female = ((g[1]*0.5 + g[2])*wDD + (g[4]*Xm1_Mm + g[5]*Xm1_mm)*wDd) / meanW
    dM_female = ((g[3]*(1-Xm1_MM) + g[4]*(1-Xm1_Mm))*wDd + (g[6] + g[7]*0.5)*wdd) / meanW
    dm_female = ((g[5]*(1-Xm1_mm))*wDd + (g[7]*0.5 + g[8])*wdd) / meanW
    #print DM_female + Dm_female + dM_female + dm_female

    #gamete freq
    DM = DM_male*0.5 + DM_female*0.5
    Dm = Dm_male*0.5 + Dm_female*0.5
    dM = dM_male*0.5 + dM_female*0.5
    dm = dm_male*0.5 + dm_female*0.5  
    
    #LD produced by driver and modifier 
    LD = DM*(1 - DM - Dm - dM) - Dm*dM
    #print i, "D freq", DM + Dm, "M freq", DM+dM, LD 
    #print DM+Dm+dM+dm
    
    #keeping track    
    onelocus_femlinkD1.append(DM+Dm)
    onelocus_femlinkMod1.append(DM+dM)
    
    if onelocus_femlinkD1[-1] - onelocus_femlinkD1[-2] < 10**-10 and onelocus_femlinkMod1[-1] - onelocus_femlinkMod1[-2] < 10**-10:
        break

#Figure 4B
GraphwithMod("Female-Limited Drive, Modifier Linked to Non-Driver", onelocus_femlinkD1, "Driver freq", onelocus_femlinkMod1, "Modifier freq")
        
## ------------------------------------
## drive is in both sexes
## recombination modifier linked to driver haplotype

onelocus_linkedD2 = []
onelocus_linkedMod2 = []
onelocus_linkedD2.append(D)
onelocus_linkedMod2.append(M)

#recombination suppressors are favored in this model
deltaR = -0.05
Xm1_mm = driveStrength(r, alphaM1)
Xm1_Mm = driveStrength(r + deltaR, alphaM1)
Xm1_MM = driveStrength(r + 2*deltaR, alphaM1)
    
for i in range(Gens):
    
    #intialize first generation    
    if i == 0:
        Dtemp = D
        Mtemp = M
        DM = Dtemp * Mtemp
        #Dm = 0.0
        #dM = 0.0 
        Dm = Dtemp * (1-Mtemp)
        dM = (1-Dtemp) * (Mtemp)
    #genotype = [DDMM, DDMm, DDmm, DdMM, DdMm, Ddmm, ddMM, ddMm, ddmm]    
    g = [DM*DM, 2*DM*Dm, Dm*Dm, 2*DM*dM, 2*DM*(1-DM-Dm-dM)+2*Dm*dM, 2*Dm*(1-DM-Dm-dM), 
         dM*dM, 2*dM*(1-DM-Dm-dM), (1-DM-Dm-dM)*(1-DM-Dm-dM)]  #genotype freq
    
    #mean population fitness     
    meanW = (g[0]+g[1]+g[2])*wDD + (g[3]+g[4]+g[5])*wDd + (g[6]+g[7]+g[8])*wdd
    
    #gamete haplotype frequencies, fitness/selection and drive taken into account
    DM = ((g[0] + g[1]*0.5)*wDD + (g[3]*Xm1_MM + g[4]*Xm1_Mm)*wDd) / meanW
    Dm = ((g[1]*0.5 + g[2])*wDD + (g[5]*Xm1_mm)*wDd) / meanW
    dM = ((g[3]*(1-Xm1_MM)*wDd) + (g[6] + g[7]*0.5)*wdd) / meanW
    dm = ((g[4]*(1-Xm1_Mm) + g[5]*(1-Xm1_mm))*wDd + (g[7]*0.5 + g[8])*wdd) / meanW
    
    #LD produced by driver and modifier 
    LD = DM*(1 - DM - Dm - dM) - Dm*dM
    #print i, "D freq", DM + Dm, "M freq", DM+dM, LD 
    #print DM+Dm+dM+dm
    
    #keeping track    
    onelocus_linkedD2.append(DM+Dm)
    onelocus_linkedMod2.append(DM+dM)
    
    if onelocus_linkedD2[-1] - onelocus_linkedD2[-2] < 10**-10 and onelocus_linkedMod2[-1] - onelocus_linkedMod2[-2] < 10**-10:
        break

#Figure 4C
GraphwithMod("Meoitic Drive, Modifier Linked to Driver", onelocus_linkedD2, "Driver freq", onelocus_linkedMod2, "Modifier freq")
 
## ------------------------------------
## drive is female limited
## recombination modifier linked to driver haplotype

onelocus_femlinkD2 = []
onelocus_femlinkMod2 = []
onelocus_femlinkD2.append(D)
onelocus_femlinkMod2.append(M)

#recombination suppressors are favored in this model
deltaR = -0.05
Xm1_mm = driveStrength(r, alphaM1)
Xm1_Mm = driveStrength(r + deltaR, alphaM1)
Xm1_MM = driveStrength(r + 2*deltaR, alphaM1)
    
for i in range(Gens):
    
    #intialize first generation    
    if i == 0:
        Dtemp = D
        Mtemp = M
        DM = Dtemp * Mtemp
        Dm = Dtemp * (1-Mtemp)
        dM = (1-Dtemp) * (Mtemp) 
    
    #genotype = [DDMM, DDMm, DDmm, DdMM, DdMm, Ddmm, ddMM, ddMm, ddmm]    
    g = [DM*DM, 2*DM*Dm, Dm*Dm, 2*DM*dM, 2*DM*(1-DM-Dm-dM)+2*Dm*dM, 2*Dm*(1-DM-Dm-dM), 
         dM*dM, 2*dM*(1-DM-Dm-dM), (1-DM-Dm-dM)*(1-DM-Dm-dM)]  #genotype freq
    
    #mean population fitness     
    meanW = (g[0]+g[1]+g[2])*wDD + (g[3]+g[4]+g[5])*wDd + (g[6]+g[7]+g[8])*wdd
    
    DM_male = ((g[0] + g[1]*0.5)*wDD + (g[3]*0.5 + g[4]*0.5)*wDd) / meanW
    Dm_male = ((g[1]*0.5 + g[2])*wDD + (g[5]*0.5)*wDd) / meanW
    dM_male = ((g[3]*0.5)*wDd + (g[6] + g[7]*0.5)*wdd) / meanW
    dm_male = ((g[4]*0.5 + g[5]*0.5)*wDd + (g[7]*0.5 + g[8])*wdd) / meanW
    #print DM_male + Dm_male + dM_male + dm_male

    #female gamete freq after fitness/selection and drive
    DM_female = ((g[0] + g[1]*0.5)*wDD + (g[3]*Xm1_MM + g[4]*Xm1_Mm)*wDd) / meanW
    Dm_female = ((g[1]*0.5 + g[2])*wDD + (g[5]*Xm1_mm)*wDd) / meanW
    dM_female = ((g[3]*(1-Xm1_MM)*wDd + (g[6] + g[7]*0.5)*wdd)) / meanW
    dm_female = ((g[4]*(1-Xm1_Mm) + g[5]*(1-Xm1_mm))*wDd + (g[7]*0.5 + g[8])*wdd) / meanW
    #print DM_female + Dm_female + dM_female + dm_female

    #gamete freq
    DM = DM_male*0.5 + DM_female*0.5
    Dm = Dm_male*0.5 + Dm_female*0.5
    dM = dM_male*0.5 + dM_female*0.5
    dm = dm_male*0.5 + dm_female*0.5  
    
    #LD produced by driver and modifier 
    LD = DM*(1 - DM - Dm - dM) - Dm*dM
    #print i, "D freq", DM + Dm, "M freq", DM+dM, LD 
    #print DM+Dm+dM+dm
    
    #keeping track    
    onelocus_femlinkD2.append(DM+Dm)
    onelocus_femlinkMod2.append(DM+dM)
    
    if onelocus_femlinkD2[-1] - onelocus_femlinkD2[-2] < 10**-10 and onelocus_femlinkMod2[-1] - onelocus_femlinkMod2[-2] < 10**-10:
        break

#Figure 4D
GraphwithMod("Female-Limited Drive, Modifier Linked to Driver", onelocus_femlinkD2, "Driver freq", onelocus_femlinkMod2, "Modifier freq")
 
