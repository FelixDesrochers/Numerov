"""
Numerov Schrodinger equation solver

Description: This script solves the 1 dimensional time-independant Schrodinger equation for any given potential. Its takes a potentiel as an entry and outputs the wanted energy level and
             a figure with the potentiel and the wave fonctions corresponding to the energy level that have been specicied.

Indications: The script ask for the number of energy levels that are desired and the potential. The potential must be centered at x=0 (the programm will itself translate the potential in y so
             that all the values are positive). Also if the potential behaves weirdly or the desired number of energy level is big, you may encounter some problem with the MeetingPoints function.
             Usually changing the x_V_min and x_V_max fixes the problem.

author: Félix Desrochers
email: felix.desrochers@polymtl.ca

MIT License

Copyright (c) 2017 Félix Desrochers

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

###################
# Importing modules
###################

import numpy as np

#Imports the Fct_Numerov module which defines many functions that will be used in this script
import Fct_Numerov


############################
# 1) Initializing parameters
############################

#Indication :Theses parameters determine the precision of the calculations and can be adjust as wanted

#Setting the range from wich we will evaluate the potential and the number of division we will use to discretise the potential
x_V_min = -13
x_V_max = 13
nbr_division_V = 800000

#Setting the number of division from the initial point in the classical forbidden zone x_0 to the ending point x_max
nbr_division = 5000

#Setting the initial augmentation after the point where the wave function will be set to zero
Initial_augmentation = 0.00001

#Setting the tolerance for the wave fonction at the ending point (x_max) to accept the energy level as the wnated energy level
Tolerance = 0.00000001


###########################################################################
# 2) Entering the parameters concerning the energy levels and the potential
###########################################################################

# i) Energy levels
E_level = int(input('Which first energy levels do you want (enter an integer) : '))
E_lvl = list(range(0,E_level))

# ii) Potential
potential=input('Potential (as a fonction of x): ')

#Verify if the potential expression is correct (Syntax, bounadries value and "global concavity")
i=1
while i == 1:
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

    #Translate the potential
    PositionPotential,PotentialArray = Fct_Numerov.TranslationPotential(PositionPotential, PotentialArray)

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
WaveFunctionFound = {} #Defines the wave functions that have been found. Has the structure {0:WaveFunction0, 1:WaveFunction1, ...}
E_guess_try = {} #Defines the lowest and higest energy levels that have been used so far for each number of nodes. Has the structure {NbrNodes1:[Energy guessed min, Energy guessed max], ...}
iteration = 1 #Defines the number of iterations
E_found = list() #A list of the energy levels that have been found (ex: [0,2,3] if the energy levels 0,2 and 3 have been found)
# Note: All the wave function are a list of tuple and have the following structure: [(x0, \psi(x0), (x1, \psi(x1)), ...]

#Continue while we don't have the n first energy level
while not E_found == list(range(E_level)):

    #########################################################
    # i) Initial Energy guess

    E_guess = Fct_Numerov.E_Guess(EnergyLevelFound,E_guess_try,iteration, First_E_guess)
    print('E_guess: ', E_guess)

    ##########################################################
    # ii) Setting the initial and final points (where \psi=0)

    #Gets the meeting points with the energy and the potential
    MeetingPoints,end_program,E_guess = Fct_Numerov.MeetingPointsPotential(E_guess, PotentialArray, PositionPotential, E_guess_try)

    if end_program:
        break

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

    ###################################################################################
    # vii) Updates the Energy levels found list to verify if the condition is respected
    E_found = list()
    for i in EnergyLevelFound.keys():
        E_found.append(i)
    E_found.sort()
    print('Energy level found',EnergyLevelFound)


######################################
# 4) Output (energy levels and figure)
######################################

# i) Energy levels
Fct_Numerov.OuputEnergy(EnergyLevelFound)

# ii) Figure
#Get all the information about what to draw
y_max,min_x,max_x,WavPlot,WavLines,EnergyLines = Fct_Numerov.DefineWhatToPlot(WaveFunctionFound,EnergyLevelFound)

#Draw all the wave functions
Fct_Numerov.DrawWaveFunction(y_max, min_x, max_x, WavPlot, WavLines, EnergyLines, PositionPotential, PotentialArray)



