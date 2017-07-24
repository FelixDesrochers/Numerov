
#Importing the necessary modules
import sys
sys.path.append('/Users/anaconda/lib/python3.6/site-packages')

import math
import numpy as np

######################################################################################
# 1) Potential functions
# These functions are used to modify and be sure that the entered potential is correct
######################################################################################


def ModifyPotential(potential):
    '''This fonction replaces any mathematical expression that is usually used but that is incorect in python.
    For instance:
        x^2 -> x**2
        |x| -> math.fabs(x)'''

    #Relacing exponential
    potential = potential.replace('^','**')

    #Replacing absolute value
    pot_list = potential.rsplit('|')
    for i in [ i for i in range(1,(len(pot_list)-1)*2) if i%2==1 ]:
        insertion = 'math.fabs(' if i%4 ==1 else ')'
        pot_list.insert(i,insertion)

    potential=''.join(pot_list)

    return potential


def VerifySyntaxPotential(potential):
    ''' Verify if the potential entered has an invalid syntax and demands another potential untils there is no more syntax error '''
    i=0
    while i == 0:
        #Tries to evaluate the potential at x=0 and asks for a new one until there is a syntax error
        try:
            x=0
            eval(potential)
        except SyntaxError:
            potential = input('The potential is incorrect enter a new one: ')
            potential = ModifyPotential(potential)
        else:
            i=1

    return potential


def VerifyLimitsPotential(potential):
    '''Verify if the potential seems to verify the borders conditions (to allow bound states) if it doesn't it ask to the user if he is sure that the potential respects these conditions'''
    #Verify if the potential is bigger than V(x=0) for x=100 and x=-100
    i=1
    while i == 1:
        eval_pot = list()
        x=-100
        eval_pot.append(eval(potential))
        x=100
        eval_pot.append(eval(potential))
        eval_pot = np.array(eval_pot)
        x = 0
        #if it doesn't ask for a new potential
        if eval_pot[eval_pot < eval(potential)]:
            QuestionPotential = input('The potential doesn\'t seem to be correct. Are you it corresponds to a bound state (y/n)? ')

            if QuestionPotential == 'n':
                potential = input('Enter a new potential: ')
                #Check the syntax for the new potential
                potential = ModifyPotential(potential)
                potential = VerifySyntaxPotential(potential)
            elif QuestionPotential == 'y':
                i = 0

        #If it respects the condition, exit the while loop
        else :
                i = 0

    return potential


def VerifyConcavity(PotentialArray):
    '''Evaluates the concavity of the potential and returns its value positive if the concavity is correct or negative if it is incorrect'''
    i=1
    E_verify=2
    index_min=[]
    index_max=[]

    for i in range(0,len(PotentialArray)-2):
        #Gets all the points where the potential meets the E_verify value
        if PotentialArray[i] > E_verify and PotentialArray[i+1] < E_verify:
            index_min.append(i)

        elif PotentialArray[i] < E_verify and PotentialArray[i+1] > E_verify:
            index_max.append(i)

        elif PotentialArray[i] == E_verify:
            if PotentialArray[i-1]>E_verify and PotentialArray[i+1]<E_verify:
                index_min.append(i)

            elif PotentialArray[i-1]<E_verify and PotentialArray[i+1]>E_verify:
                index_max.append(i)

        #Gets the concavity value
        index_max = np.array(index_max)
        index_min = np.array(index_min)

        if (index_max.max() > index_min.max()) and (index_min.min() < index_min.min()):
            concavity = 'positive'
        else:
            concavity = 'negative'
    return concavity


def EvaluateOnePotential(position,potential):
    '''Defines a function that evaluate the potential at a certain point x. This function will be vectorized with np.vectorize to evaluate the potential on a list of position [x1,x2,...]'''
    x = position
    EvalPotential = eval(potential)

    return EvalPotential


##################################################
# 2) Numerov algorithm functions
# Defines the functions used in the Numerov method
##################################################

#########################
# i) Initial Energy guess

def E_Guess(EnergyLevelFound, E_guess_try):
    '''Defines the energy guess depending on the energy levels that have been found and on the energy that have already been guessed. '''

    # I) Define the energy level that we want to find E_level_guess (the lowest energy level that hasn't been found yet)
    #List for the energy that have been found
    Lvl_found = list(EnergyLevelFound.keys())
    Lvl_found.sort()
    #Gets the energy level that we want to find fo this guess
    E_level_missing = [Energy for index,Energy in enumerate(Lvl_found) if not Energy == index]
    if E_level_missing == None:
        E_level_guess = max(Lvl_found) +1
    else:
        E_level_missing.sort()
        E_level_guess = E_level_missing[0]

    # II) Defining the energy guess depending on the guess that have already been done (E_guess_try)
    #Finds the closest energy energy level (number of nodes) that has been guessed and that corresponds to a smaller or an equal number of nodes than E_level_guess
    try:
        E_level_smaller = max([ E for E in E_guess_try.keys() if E <= E_level_guess ])
    except ValueError:
        E_level_smaller = None
    #Finds the closest energy energy level (number of nodes) that has been guessed and that corresponds to a bigger number of nodes than E_level_guess
    try:
        E_level_bigger = min([ E for E in E_guess_try.key() if E > E_level_guess ])
    except ValueError:
        E_level_bigger = None

    #Define the energy guess
    #If the smaller and higher exist take the average
    if E_level_smaller and E_level_bigger:
        E_guess = ( E_guess_try[E_level_smaller][1] + E_guess_try[E_level_bigger][0] ) / 2

    #If only the higher exists take the half
    elif E_level_bigger:
        E_guess = E_guess_try[E_level_bigger][0]/2

    #If only the smaller exists take the double
    elif E_level_smaller:
        E_guess = E_guess_try[E_level_smaller][1] * 2

    return E_guess

##################################################################################
# ii) Setting the minimal and maximal points (where the wave function equals zero)

def MeetingPointsPotential(E_guess, PotentialArray, PositionPotential):
    '''Finds the minimal and maximal points where the energy that has been guessed is equal to the potential.
    Parameters:
        E_guess: the guessed energy
        PotentialArray: a Numpy array that contains the potential for certain points
        PositionPotential: a Numpy array that contains the positions that correspond to the potential array
    Returns:
        MeetingPoints: a tuple of the smallest and biggest meeting point'''

    MeetingPoints = [None,None] # a list containing the smallest and highest meeting point

    for i in range(0,len(PotentialArray)-2):
        #Gets all the meeting points
        if (PotentialArray[i] < E_guess and PotentialArray[i+1] > E_guess) or (PotentialArray[i] > E_guess and PotentialArray[i+1] < E_guess) or PotentialArray[i] == E_guess:
            #And filter them
            if (MeetingPoints[0] == None) or (PositionPotential[i] < MeetingPoints[0]):
                MeetingPoints[0] = PositionPotential[i]
            elif (MeetingPoints[1] == None) or (PositionPotential[i] > MeetingPoints[1]):
                MeetingPoints[1] = PositionPotential[i]

    MeetingPoints = tuple(MeetingPoints)
    return MeetingPoints

def DetermineMinAndMax(MeetingPoints,x_V_min,x_V_max):
    '''This function determines the minimal and maximal position where the wave function will be set to 0 depending on the points where the potential meets the guess energy and on
    the minimum and maximum that are initially set for the potential

    Parameter:
        MeetingPoints: the minimum and maximum point where the potentila meets the guessed energy
        x_V_min: The minimum value of the position for the potential
        x_V_max: The maximum value of the poisition for the potential

    Returns:
        Position_min: the minimum value where psi=0
        Position_max: the maximum value where psi=0'''

    #Sets the min and max as the medium value between the meeting with the potential and the ending
    Position_min = (x_V_min + MeetingPoints[0])/2
    Position_max = (x_V_max + MeetingPoints[1])/2

    return Position_min,Position_max

#######################################
# iii) Calculate the wave function

def WaveFunctionNumerov(potential, E_guess, nbr_division, Initial_augmentation, Position_min, Position_max):
    '''This function calculates the wave function values depending on the x coordinate by using the Numerov method. The function returns a list that contains tuple with the x coordinate and
    the wave function value. It has the general form: [(x0, psi(x0)), (x1, psi(x1)), ...]'''

    #Initializing the wave function
    WaveFunction = []

    #Setting the divisions
    Division = (Position_max - Position_min) / nbr_division

    #Setting constant
    HBar = 1.054571726 * (10**(-34))
    MassElectron = 9.10938215 * (10**(-31))

    #Setting the first values of the wave function
    WaveFunction.append((Position_min,0))
    WaveFunction.append((Position_min+Division, Initial_augmentation))

    #Defing an array and an index to use in the for loop
    index = 0
    PositionArray = np.arange(Position_min, Position_max, Division)

    #Calculating the wave function for other values
    for i in np.arange(Position_min + 2 * Division, Position_max + Division, Division):
        #Evaluating the potential
        #For V_i+1
        x = i
        V_plus1 = eval(potential)

        #For V_i
        x = PositionArray[index]
        V = eval(potential)

        #For V_i-1
        x = PositionArray[index + 1]
        V_minus1 = eval(potential)

        #Setting the k**2 values ( where k**2 = (2m/HBar)*(E-V(x)) )
        k_2_plus1 = ((2*MassElectron)/(HBar**2)) * (E_guess - V_plus1)
        k_2 = ((2*MassElectron)/(HBar**2)) * (E_guess - V)
        k_2_minus1 = ((2*MassElectron)/(HBar**2)) * (E_guess - V_minus1)

        #Calculating the wave function
        psi = ((2 * (1 - (5/12) * (Division**2) * (k_2)) * WaveFunction[-1][1]) - (1 + (1/12) * (Division**2) * k_2_minus1 ) * WaveFunction[-2][1] ) / (1 + (1/12) * (Division**2) * k_2_plus1)

        #Saving the wave function and the x coordinate
        WaveFunction.append((i,psi))

        #Incrementing the index
        index += 1

    return WaveFunction




########################################################
# iv) Determine the number of nodes in the wave function

def NumberNodes(WaveFunction):
    '''This function evaluates the number of nodes in the wavefunction. The number of nodes will allow us the determine the energy level to which a certain wave function corresponds'''

    #Initialize the number of nodes
    NumberOfNodes = 0

    #Calculate the number of nodes
    for i in range(0,len(WaveFunction)-2):

        if (WaveFunction[i][1] > 0 and WaveFunction[i+1][1] < 0) or (WaveFunction[i][1] > 0 and WaveFunction[i+1][1] < 0 ) or (WaveFunction[i][1] == 0):
            NumberOfNodes += 1

    return NumberOfNodes



#####################################################
# v) Verify if wave function respects the restriction

def VerifyTolerance(WaveFunction, Tolerance):
    '''See if the wave function for the given energy level respects the tolerance. Returns yes if it respects the tolerance and no if not.'''

    VerificationTolerance = 'yes' if WaveFunction[-1][1] < Tolerance else 'no'

    return VerificationTolerance

#######################################################
# vi) Saves energy and the correponding number of nodes

def SaveEnergy(NumberOfNodes, E_guess, E_guess_try):
    '''This function saves the guessed energy and the number of nodes corresponding to it. It return a dictionnary that has the general form: {NumberofNodes:[E_min, E_max], NumberOfnodes:[E_min, E_max]}'''

    #Checks if the key Number of Nodes exists. If it doesn't, define the two values in the list corresponding to the key NumberOfNodes as E_guess.
    try:
        E_guess_try[NumberOfNodes]

    except KeyError:
        E_guess_try[NumberOfNodes] = [E_guess, E_guess]
        return E_guess_try

    #Checks if the energy guess is smaller than the smallest value in the list
    if E_guess < E_guess_try[NumberOfNodes][0]:
        E_guess_try[NumberOfNodes][0] = E_guess

    #Checks if the energy guess is greater than the biggest value in the list
    elif E_guess > E_guess_try[NumberOfNodes][1]:
        E_guess_try[NumberOfNodes][1] = E_guess

    return E_guess_try
