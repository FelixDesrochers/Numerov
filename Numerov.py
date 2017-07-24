'''This script solves the 1 dimensional time-independant Schrodinger equation for any given potential. Its takes a potentiel as an entry and outputs the wanted energy level and a figure with
the potentiel and the wave fonctions corresponding to the energy level that have been specicied.

The one dimensionl time independent Schrodinger is as follow:

    \frac{-\hbar^2}{2*m}*\frac{d}{dx} \psi(x) +  V(x) \psi(x) = E * \psi(x)

To solve this diffenrential equation a Numerov method wil be used. This method aims to solve any equation of the type : \frac{d}{dx}y(x) + s(x)*x =0.

The script takes as inputs the wanted energy levels and the potnetial. Other parameters can be adjust in the script to modify the precision such as the tolerance for the wave fonction at the
ending point and the number of division'''

###################
# Importing modules
###################

#Append path to acces all modules in anaconda
import sys
sys.path.append('/Users/felix/anaconda/lib/python3.6/site-packages')

import math
import numpy as np
#import matplotlib.pyplot as plt

#Imports the Fct_Numerov module which defines many functions that will be used in this script
import Fct_Numerov

############################
# 1) Initializing parameters
############################

#Indication :Theses parameters determine the precision of the calculations and can be adjust as wanted

#Setting the range from wich we will evaluate the potential and teh number of division we will include
x_V_min = -50.0
x_V_max = 50.0
nbr_division_V = 500

#Setting the number of division from the initial point in the classical forbidden zone x_0 t the ending point x_max
nbr_division = 200

#Setting the initial augmentation after the point where the wave function will be set to zero
Initial_augmentation = 0.01

#Setting the tolerance for the wave fonction at the ending point (x_max) to accept the energy level as the wnated energy level
Tolerance = 0.01


###########################################################################
# 2) Entering the parameters concerning the energy levels and the potential
###########################################################################

# i) Energy levels
E_level = input('Which first energy levels do you want (enter an integer) : ')
E_lvl = E_level.split(',')

# ii) Potential
potential=input('Potential (as a fonction of x): ')

#Verify if the potential expression is correct (Syntax, bounadries value and "global concavity")
i=1
while i ==1:
    #Changes the expression to be sure it matches python mathematical syntax
    potential = Fct_Numerov.ModifyPotential(potential)

    #Verify if the potential expression has any syntax error
    potential = Fct_Numerov.VerifySyntaxPotential(potential)

    #Verify if the potential seems to respect the boundaries conditions
    potential = Fct_Numerov.VerifyLimitsPotential(potential)

    #Convert the potential into a numpy array (see the settings for this potential array in the "Initializing parameters section")
    EvaluatePotential = np.vectorize(Fct_Numerov.EvaluateOnePotential)
    DivisionPotential = (x_V_max - x_V_min) / nbr_division_V
    PositionPotential = np.arange(x_V_min,x_V_max,DivisionPotential)

    PotentialArray = EvaluatePotential(PositionPotential,potential)

    #Verify the concavity of the potential
    concavity = Fct_Numerov.VerifyConcavity(PotentialArray)

    #If it is correct exit the loop
    if concavity == 'positive':
        i = 0
    #Else ask for a new one or take this one anyway
    elif concavity == 'negative':
        potential2 = input('The concavityof the potential isn\'t correct enter a new one (or "O" to overule): ')

        if potential2 == 'O':
            i = 0
        else :
            potential = potential2


###################################
# 3) Numerov algorithm
###################################

#Initializing constants
MassElectron = 9.10938215 * (10**(-31))
HBar = 1.054572 * (10**(-34))

#Initializing paramaters for the while loop
EnergyLevelFound = {} #Defines energy levels that avec been found. Has the structure {0:E0, 1:E1, 2:E2, ...}
E_guess_try = {} #Defines the lowest and higest energy levels that have been used so far for each number of nodes. Has the structure {NbrNodes1:[Energy guessed min, Energy guessed max], ...}
E_First_Guess = (HBar * ((MassElectron)**(-1/2)))/2 #Takes as intial guess the energy level of the quantum harmonic oscillator with k=1
iteration = 1 #Defines the number of iterations

while not len(EnergyLevelFound) == E_level:

    #########################################################
    # i) Initial Energy guess

    E_guess = E_First_Guess if iteration == 1 else Fct_Numerov.E_Guess(EnergyLevelFound,E_guess_try)

    ##########################################################
    # ii) Setting the initial and final points (where \psi =0)

    #Gets the meeting points with the energy and the potential
    MeetingPoints = Fct_Numerov.MeetingPointsNumerov(E_guess, PotentialArray, PositionPotential)

    #Sets the minimum and maximum value for the position where the wave function equals zero
    Position_min,Position_max = Fct_Numerov.DetermineMinAndMax(MeetingPoints, x_V_min, x_V_max)

    ###############################################################
    # iii) Calculate the wave fonction for the guessed energy value

    WaveFunction = Fct_Numerov.WaveFunctionNumerov(potential, E_guess, nbr_division, Initial_augmentation, Position_min, Position_max)

    ##########################################################
    # iv) Determine the number of nodes in the wave fonction

    NumberOfNodes = Fct_Numerov.NumberNodes(WaveFunction)


    ####################################################################################
    # v) See if the wave fonction for this energy respects the restriction (if yes save)

    VerificationTolerance = Fct_Numerov.VerifyTolerance(WaveFunction,Tolerance)

    if VerificationTolerance == 'yes':
        EnergyLevelFound.update({NumberOfNodes:WaveFunction})

    ######################################################################################
    # vi) Saves Energy guess and the corresponding number of nodes (no matter if it fails)

    E_guess_try = Fct_Numerov.SaveEnergy(NumberOfNodes, E_guess, E_guess_try)

    iteration += 1

##############################################
#Temporary try to see what the script displays
for i,Energy in EnergyLevelFound.items():
    print(i, Energy)


######################################
# 4) Output (energy levels and figure)
######################################

# i) Figure

# ii) Energy levels
