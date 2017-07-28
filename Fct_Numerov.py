
#Importing the necessary modules
import sys
sys.path.append('/Users/anaconda/lib/python3.6/site-packages')

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

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
            potential = input('The potential\'s syntax is incorrect enter a new one: ')
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


def GetFirstEnergyGuess(PotentialArray):
    '''Defines the first energy level as a value between the the average potential and the minimum value. More explicitly: (1/50000)*((V_average + V_min)/2)'''

    First_E_guess = PotentialArray.min() +  (1/500000) * (PotentialArray.mean() + PotentialArray.min())

    return First_E_guess


def VerifyConcavity(PotentialArray, First_E_guess):
    '''Evaluates the concavity of the potential and returns its value positive if the concavity is correct or negative if it is incorrect'''

    i = 1

    while i == 1:
        print('First Energy guess:', First_E_guess)
        index_min=list()
        index_max=list()
        try:
            for i in range(0,len(PotentialArray)-2):
                #Gets all the points where the potential meets the E_verify value
                if PotentialArray[i] > First_E_guess and PotentialArray[i+1] < First_E_guess:
                    index_min.append(i)

                elif PotentialArray[i] < First_E_guess and PotentialArray[i+1] > First_E_guess:
                    index_max.append(i)

                elif PotentialArray[i] == First_E_guess:
                    if PotentialArray[i-1] > First_E_guess and PotentialArray[i+1] < First_E_guess:
                        index_min.append(i)

                    elif PotentialArray[i-1] < First_E_guess and PotentialArray[i+1] > First_E_guess:
                        index_max.append(i)

            #Gets the concavity value
            index_max = np.array(index_max)
            index_min = np.array(index_min)

            print('index max: ',index_max)
            print('index_min: ',index_min)

            if (index_max.max() > index_min.max()) and (index_max.min() > index_min.min()):
                concavity = 'positive'
            else:
                concavity = 'negative'

        except ValueError:
            First_E_guess = First_E_guess/2

        else:
            i = 0

    return concavity,First_E_guess


def EvaluateOnePotential(position,potential):
    '''Defines a function that evaluate the potential at a certain point x. This function will be vectorized with np.vectorize to evaluate the potential on a list of position [x1,x2,...]'''
    x = position
    EvalPotential = eval(potential)

    return EvalPotential

def GetTranslation(potential):
    '''Checks approximately where the minimum of the potential is and outputs the necessary translation in x and y to recenter the minimum at x=0 and y=0'''

    #i) Create a potential Array
    Nbr_Div = 50000
    x_min = -2
    x_max = 2
    Div = (x_max - x_min)/Nbr_Div
    y = list()

    for i in range(Nbr_Div):
        position = x_min + (Div * i)
        y.append(EvaluateOnePotential(position,potential))

    #ii) Gets the minimum value for the potential
    trans_y = min(y)
    index = y.index(trans_y)
    trans_x = x_min + (Div * index)

    print('trans_x; ',trans_x)
    print('trans_y; ',trans_y)

    return trans_x,trans_y

def TranslatePotential(potential,trans_x,trans_y):
    '''Modify the potential expression to center its minimum at x=0 and y=0'''
    #x translation
    #potential = potential.replace('x','(x+' + str(trans_x) + ')')

    #y translation
    potential = potential + '-' +  str(trans_y)

    print(potential)

    return potential

def RecenterPotential(PositionPotential,PotentialArray):
    '''Translate the potential array to rencenter it at (0,0) '''

    #Translate the potential array
    trans_y_2 = -min(PotentialArray)
    index = list(PotentialArray).index(-trans_y_2)
    #PotentialArray = PotentialArray + trans_y_2

    #Translate the position potential
    trans_x_2 = -PositionPotential[index]
    PositionPotential = PositionPotential + trans_x_2

    return PositionPotential

##################################################
# 2) Numerov algorithm functions
# Defines the functions used in the Numerov method
##################################################

#########################
# i) Initial Energy guess

def E_Guess(EnergyLevelFound, E_guess_try, iteration, First_E_guess):
    '''Defines the energy guess depending on the energy levels that have been found and on the energy that have already been guessed. '''

    print('Iteration: ',iteration)
    #If it is the first time, return the first energy level of the quantum harmonic oscillator
    if iteration == 1:
        E_guess = First_E_guess  #Takes as intial guess the First_E_guess that has previously been defined
        return E_guess

    # I) Define the energy level that we want to find E_level_guess (the lowest energy level that hasn't been found yet)
    #List for the energy that have been found
    Lvl_found = list(EnergyLevelFound.keys())
    Lvl_found.sort()
    #Gets the energy level that we want to find
    E_level_missing = [index for index,Energy in enumerate(Lvl_found) if not Energy <= index]
    if not E_level_missing:
        if not Lvl_found:
            E_level_guess = 0
        else:
            E_level_guess = max(Lvl_found) +1
    else:
        E_level_guess = min(E_level_missing)

    # II) Defining the energy guess depending on the guess that have already been done (E_guess_try)
    #Finds the closest energy energy level (number of nodes) that has been guessed and that corresponds to a smaller or an equal number of nodes than E_level_guess
    try:
        E_level_smaller = max([ E for E in E_guess_try.keys() if E <= E_level_guess ])
    except ValueError:
        E_level_smaller = None
    #Finds the closest energy energy level (number of nodes) that has been guessed and that corresponds to a bigger number of nodes than E_level_guess
    try:
        E_level_bigger = min([ E for E in E_guess_try.keys() if E > E_level_guess ])
    except ValueError:
        E_level_bigger = None

    #Define the energy guess
    #If the smaller and higher exist take the average
    if (not E_level_smaller == None) and (not E_level_bigger ==None):
        E_guess = ( E_guess_try[E_level_smaller][1] + E_guess_try[E_level_bigger][0] ) / 2

    #If only the higher exists take the half
    elif not E_level_bigger == None:
        E_guess = E_guess_try[E_level_bigger][0]/2

    #If only the smaller exists take the double
    elif not E_level_smaller == None:
        E_guess = E_guess_try[E_level_smaller][1] * 2

    print('E_level_guess:', E_level_guess )
    print('E_level_bigger: ', E_level_bigger)
    print('E_level_smaller: ', E_level_smaller)

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
                print('index rencontre min: ',i)
                MeetingPoints[0] = PositionPotential[i]
            elif (MeetingPoints[1] == None) or (PositionPotential[i] > MeetingPoints[1]):
                MeetingPoints[1] = PositionPotential[i]
                print('index renccontre max: ', i)

    MeetingPoints = tuple(MeetingPoints)
    return MeetingPoints

def DetermineMinAndMax(MeetingPoints):
    '''This function determines the minimal and maximal position where the wave function will be set to 0 depending on the points where the potential meets the guess energy and on
    the minimum and maximum that are initially set for the potential

    Parameter:
        MeetingPoints: the minimum and maximum point where the potentila meets the guessed energy
        x_V_min: The minimum value of the position for the potential
        x_V_max: The maximum value of the poisition for the potential

    Returns:
        Position_min: the minimum value where psi=0
        Position_max: the maximum value where psi=0'''

    #Sets the min and max as the half of the distance between the min and the max plus the min or the max


    Position_min = MeetingPoints[0] - (MeetingPoints[1] - MeetingPoints[0])/1
    Position_max =  MeetingPoints[1] + (MeetingPoints[1] - MeetingPoints[0])/1

    print('MeetingPoint: ', MeetingPoints)
    print('min:',Position_min)
    print('max:',Position_max)
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

    #Setting the first values of the wave function
    WaveFunction.append((float(Position_min),0))
    WaveFunction.append((float(Position_min+Division), Initial_augmentation))

    #Defing an array and an index to use in the for loop
    index = 0
    PositionArray = np.arange(Position_min, Position_max, Division)

    #Calculating the wave function for other values
    for i in np.arange(Position_min + (2 * Division), Position_max, Division):
        #Evaluating the potential
        #For V_i+1
        x = i
        V_plus1 = eval(potential)

        #For V_i
        x = PositionArray[index+1]
        V = eval(potential)

        #For V_i-1
        x = PositionArray[index]
        V_minus1 = eval(potential)

        #Setting the k**2 values ( where k**2 = (2m/HBar)*(E-V(x)) )
        k_2_plus1 = 2 * (E_guess - V_plus1)
        k_2 = 2 * (E_guess - V)
        k_2_minus1 = 2 * (E_guess - V_minus1)

        #Calculating the wave function
        psi = ((2 * (1 - (5/12) * (Division**2) * (k_2)) * (WaveFunction[-1][1])) - (1 + (1/12) * (Division**2) * k_2_minus1 ) * (WaveFunction[-2][1])) / (1 + (1/12) * (Division**2) * k_2_plus1)

        #Saving the wave function and the x coordinate
        WaveFunction.append((i,psi))

        #Incrementing the index
        index += 1

    return WaveFunction




########################################################
# iv) Determine the number of nodes in the wave function

def NumberNodes(WaveFunction):
    '''This function evaluates the number of nodes in the wavefunction. The number of nodes will allow us the determine the energy level to which a certain wave function corresponds'''

    #Initialize the number of nodes and their position
    NumberOfNodes = 0
    PositionNodes = list()

    #Calculate the number of nodes
    for i in range(1,len(WaveFunction)-1):
        if (WaveFunction[i][1] > 0 and WaveFunction[i+1][1] < 0) or (WaveFunction[i][1] < 0 and WaveFunction[i+1][1] > 0) or (WaveFunction[i][1] == 0):
            NumberOfNodes += 1
            PositionNodes.append(WaveFunction[i][0])


    #Gets the biggest position
    x = list()
    for position,wave in WaveFunction:
        x.append(position)
    x_max = max(x)

    print('PositionNodes:', PositionNodes)
    #print('Position max x :', x_max)

    return NumberOfNodes,PositionNodes,x_max



#####################################################
# v) Verify if wave function respects the restriction

def VerifyTolerance(WaveFunction, Tolerance, E_guess, E_guess_try, NumberOfNodes):
    '''See if the wave function for the given energy level respects the tolerance. Returns yes if it respects the tolerance and no if not.'''

    # i) Checks if the last value of the wave function respects the tolerance
    VerificationTolerance = 'yes' if math.fabs(WaveFunction[-1][1]) < Tolerance else 'no'
    print('Last value Wave Function: ', WaveFunction[-1][1])

    # ii) Checks if the energy guess doesn't change a lot
    try:
        E_minus = E_guess_try[NumberOfNodes][1]
        E_plus = E_guess_try[NumberOfNodes + 1][0]
    except KeyError:
        pass
    else:
        if (E_guess < E_plus and E_guess > E_minus) and ((E_minus/E_plus) > 0.9999999999) :
            VerificationTolerance = 'yes'

    return VerificationTolerance

def CorrectNodeNumber(NumberOfNodes,PositionNodes,x_max,E_guess,E_guess_try):
    ''' '''
    NumberOfNodesCorrected = NumberOfNodes
    #Correct the number of nodes if E_guess is between the lowest energy for this number of nodes and the maximum for the number of nodes - 1
    try:
        if (E_guess_try[NumberOfNodes][1] > E_guess) and (E_guess_try[NumberOfNodes - 1][1] < E_guess):
            NumberOfNodesCorrected -= 1
    #If the dictionnary E_guess_try doesn't contain these keys check if the Last number of nodes is close to the maximum value in x x_max
    except KeyError:
        if (PositionNodes/x_max) > 94:
            NumberOfNodesCorrected -= 1

    return NumberOfNodesCorrected

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



#####################################
# 3) Ouput (Energy levels and figure)
####################################

############################
# i) Draw the figure

#Draw the wave Functions, the energy levels and sets the axis limits
def DrawWaveFunction(WaveFunctionFound, EnergyLevelFound, PositionPotential, PotentialArray):

    #Determine the maximum energy to set the maximum value for the y axis
    E_max = EnergyLevelFound[max(EnergyLevelFound)]
    y_max = 1.15* E_max
    Y_by_E_level = (y_max/(max(EnergyLevelFound)+2))

    #Define a new figure with two subplot: the energy levels and the corresponding wave function
    f,(En,Wav) = plt.subplots(1,2,sharey=True)

    # i) Draw the wave functions
    for i in WaveFunctionFound.keys():
        x =[]
        y= []
        for j in range(2000,len(WaveFunctionFound[i])-2000):
            x.append(WaveFunctionFound[i][j][0])
            y.append(WaveFunctionFound[i][j][1])

        x = np.array(x)
        y = np.array(y)

        mult = (0.9 * Y_by_E_level)/(2 * y.max())
        y = mult * y + (Y_by_E_level * (i+1))
        Wav.plot(x,y,'b',label=r"$\psi(x)$",zorder=3)

    #Determines the min and max in x
    min_x = x.min()
    max_x = x.max()

    #Draw line to specify wher the wave function is centered
    for i in WaveFunctionFound.keys():
        for j in range(len(x)):
            y[j] = (Y_by_E_level * (i+1))
        Wav.plot(x,y,'k--',zorder=1)

    #Sets the axis limits
    Wav.axis([min_x, max_x, 0, y_max])

    #Draw the potential
    Wav.plot(PositionPotential, PotentialArray, 'r',label='Potential',zorder=2)

    # ii) Draw the Energy levels
    for i in WaveFunctionFound.keys():
        for j in range(len(x)):
            y[j] = EnergyLevelFound[i]
        PlotColor = cm.hot(i/len(WaveFunctionFound))
        En.plot(x,y,'--',color=PlotColor,label='E'+str(i))

    #Set the axis limit
    En.axis([min_x, max_x, 0, y_max])

    #Draw the potential
    En.plot(PositionPotential, PotentialArray, 'r',label='Potential')

    # iii) Sets differents esthetic components like the legend

    #For the wave function
    Wav.set_xlabel(r'x ($a_0$)')
    Wav.set_title('Wave Function')

    #Verify if the labels reappear multiple times
    handles, labels = plt.gca().get_legend_handles_labels()
    newLabels, newHandles = [], []
    for handle, label in zip(handles, labels):
        if label not in newLabels:
            newLabels.append(label)
            newHandles.append(handle)
    Wav.legend(newHandles, newLabels, loc='upper right')

    #For the energy levels
    En.set_xlabel(r'x ($a_0$)')
    En.set_ylabel('Energy (Hartree)')
    En.set_title('Energy levels')
    En.legend(loc='upper right')


#Draw the potential
def DrawPotential(PositionPotential, PotentialArray):
    plt.plot(PositionPotential, PotentialArray, 'r')

#Displays the figure
def DisplayFigure():
    plt.show()

#############################
# ii) ouput the energy levels
def OuputEnergy(EnergyLevelFound):
    for i,Energy in EnergyLevelFound.items():
        print('Energy level', i, ':', Energy)

