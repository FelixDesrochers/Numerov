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

#Imports the Fct_Numerov module which defines many functions that will be used in this script
import Fct_Numerov

############################
# 1) Initializing parameters
############################

#Indication :Theses parameters determine the precision of the calculations and can be adjust as wanted

#Setting the range from wich we will evaluate the potential and teh number of division we will include
x_V_min = -5
x_V_max = 5
nbr_division_V = 1000000

#Setting the number of division from the initial point in the classical forbidden zone x_0 t the ending point x_max
nbr_division = 10000

#Setting the initial augmentation after the point where the wave function will be set to zero
Initial_augmentation = 0.00001

#Setting the tolerance for the wave fonction at the ending point (x_max) to accept the energy level as the wnated energy level
Tolerance = 0.00000001


###########################################################################
# 2) Entering the parameters concerning the energy levels and the potential
###########################################################################

# i) Energy levels
E_level = input('Which first energy levels do you want (enter an integer) : ')
E_level = int(E_level)
E_lvl = list(range(0,E_level))

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

    #Translate the potential
    trans_x,trans_y = Fct_Numerov.GetTranslation(potential)
    potential = Fct_Numerov.TranslatePotential(potential,trans_x, trans_y)

    #Convert the potential into a numpy array (see the settings for this potential array in the "Initializing parameters section")
    EvaluatePotential = np.vectorize(Fct_Numerov.EvaluateOnePotential)
    DivisionPotential = (x_V_max - x_V_min) / nbr_division_V
    PositionPotential = np.arange(x_V_min,x_V_max,DivisionPotential)

    PotentialArray = EvaluatePotential(PositionPotential,potential)

    #Recenters this new potential array for more accuracy
    #PotentialArray,PositionPotential = Fct_Numerov.RecenterPotential(PotentialArray,PositionPotential)

    #Defines the initial Energy guess that will be used to verify the concavity
    First_E_guess = Fct_Numerov.GetFirstEnergyGuess(PotentialArray)

    #Verify the concavity of the potential
    concavity,First_E_guess = Fct_Numerov.VerifyConcavity(PotentialArray, First_E_guess)

    #If it is correct exit the loop
    if concavity == 'positive':
        i = 0
    #Else ask for a new one or take this one anyway
    elif concavity == 'negative':
        potential2 = input('The concavity of the potential isn\'t correct enter a new one (or "O" to overule): ')

        if potential2 == 'O':
            i = 0
        else :
            potential = potential2

###################################
# 3) Numerov algorithm
###################################

#Initializing paramaters for the while loop
EnergyLevelFound = {} #Defines energy levels that avec been found. Has the structure {0:E0, 1:E1, 2:E2, ...}
WaveFunctionFound = {} #defines the wave functions that have been found. Has the structure {0:WaveFunction0, 1:WaveFunction1, ...}
E_guess_try = {} #Defines the lowest and higest energy levels that have been used so far for each number of nodes. Has the structure {NbrNodes1:[Energy guessed min, Energy guessed max], ...}
iteration = 1 #Defines the number of iterations

#Continue while we don't have the n first energy level
E_found = list()

while not E_found == list(range(E_level)):
#while not len(E_found) == E_level:

    #########################################################
    # i) Initial Energy guess

    E_guess = Fct_Numerov.E_Guess(EnergyLevelFound,E_guess_try,iteration, First_E_guess)
    print('E_guess: ', float(E_guess))

    ##########################################################
    # ii) Setting the initial and final points (where \psi =0)

    #Gets the meeting points with the energy and the potential
    MeetingPoints = Fct_Numerov.MeetingPointsPotential(E_guess, PotentialArray, PositionPotential)

    #Sets the minimum and maximum value for the position where the wave function equals zero
    Position_min,Position_max = Fct_Numerov.DetermineMinAndMax(MeetingPoints)

    ###############################################################
    # iii) Calculate the wave fonction for the guessed energy value

    WaveFunction = Fct_Numerov.WaveFunctionNumerov(potential, E_guess, nbr_division, Initial_augmentation, Position_min, Position_max)

    ###############################################################################
    # iv) Determine the number of nodes in the wave fonction and set the tolerance

    NumberOfNodes,PositionNodes,x_max = Fct_Numerov.NumberNodes(WaveFunction)
    print('NumberOfNodes:', NumberOfNodes)

    ####################################################################################
    # v) See if the wave fonction for this energy respects the restriction (if yes save)

    VerificationTolerance = Fct_Numerov.VerifyTolerance(WaveFunction,Tolerance,E_guess,E_guess_try,NumberOfNodes)

    if VerificationTolerance == 'yes':
        print('Niveau d\'energie trouve!!\n\n')
        NumberOfNodesCorrected = Fct_Numerov.CorrectNodeNumber(NumberOfNodes,PositionNodes,x_max,E_guess,E_guess_try)
        EnergyLevelFound.update({NumberOfNodesCorrected:E_guess})
        WaveFunctionFound.update({NumberOfNodesCorrected:WaveFunction})

    ######################################################################################
    # vi) Saves Energy guess and the corresponding number of nodes (no matter if it fails)

    E_guess_try = Fct_Numerov.SaveEnergy(NumberOfNodes, E_guess, E_guess_try)
    print('E_guess_try: ',E_guess_try)

    #Increments the iteration
    print('iterations:',iteration,'\n')
    iteration += 1

    ############################################################
    # vii) Verify if the condition is respected
    E_found = list()
    for i in EnergyLevelFound.keys():
        E_found.append(i)
    E_found.sort()
    print('Energy level found',EnergyLevelFound)
    print('E_found: ',E_found)



######################################
# 4) Output (energy levels and figure)
######################################

# i) Figure

#Draw all the wave functions
Fct_Numerov.DrawWaveFunction(WaveFunctionFound, EnergyLevelFound, PositionPotential, PotentialArray)

# ii) Energy levels
Fct_Numerov.OuputEnergy(EnergyLevelFound)

#Displays the figure
Fct_Numerov.DisplayFigure()
