"""
Module for the Numerov Scrodinger equation solver

Description: This module defines all the necessary functions that are used in the main script Numerov.py

author: Félix Desrochers
email: felix.desrochers@polymtl.ca
license: copyleft
Feel free to modify and improve this code, but preserve any derivative of it open source and keep the information above. Thanks!
"""

####################################
#Importing the necessary modules
####################################
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import animation


######################################################################################
# 1) Potential functions
# These functions are used to modify and be sure that the entered potential is correct
######################################################################################


def ModifyPotential(potential):
    """This fonction replaces any mathematical expression that is usually used but that is incorect in python.

    Parameters:
    ----------
        potential (str) : a string that indicates the mathematical form of the potential

    Returns:
    --------
        potential (str) : a new potential that now has changed any mathematical expression that is usually used but that is incorrect in python
        For instance:
            x^2 -> x**2
            |x| -> math.fabs(x)

    """

    #Relacing exponential
    potential = potential.replace('^','**')

    #Replacing absolute value
    pot_list = potential.rsplit('|')
    for i in [ i for i in range(1,(len(pot_list)-1)*2) if i%2==1 ]:
        insertion = 'np.absolute(' if i%4 ==1 else ')'
        pot_list.insert(i,insertion)

    potential=''.join(pot_list)

    #Replacing trigonometric functions
    potential = potential.replace('cos','np.cos')
    potential = potential.replace('sin','np.sin')
    potential = potential.replace('tan','np.tan')

    return potential


def VerifySyntaxPotential(potential):
    """ Verify if the potential entered has an invalid syntax and demands another potential untils there is no more syntax error

    Parameters:
    -----------
        potential (str) : a string that indicates the mathematical form of the potential

    Returns:
    --------
        potential (str) : a new string with a valid python mathematical syntax

    """

    i=0
    while i == 0:
        #Tries to evaluate the potential at x=0 and asks for a new one until there is no more syntax error
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
    """Verify if the potential seems to verify the borders conditions (to allow bound states). If it doesn't it ask to the user if he is sure that the potential respects these conditions

    Parameters:
    -----------
        potential (str) : a string that indicates the mathematical form of the potential

    Returns:
    --------
        potential (str) : a new string with a valid python mathematical syntax and with value bigger than V(x=0) for x=-100 and x=100

    """

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

        #if it doesn't respect the condition ask for a new potential
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
    """Defines the first energy level as a value between the the average potential and the minimum value. More explicitly: (1/50000)*((V_average + V_min)/2)


    Parameters:
    -----------
        PotentialArray (numpy.ndarray) : a numpy array that contains the potential value between 'x_V_min' and 'x_V_max' at every interval of length 'Division'

    Returns:
    --------
        First_E_guess (float) : the first energy guess that will be used in the Numerov algorithm. It correponds to the average of the minimum value of the potential and the average of the
                                potential times 1/50000

    """

    First_E_guess = PotentialArray.min() +  (1/500000) * (PotentialArray.mean() + PotentialArray.min())

    return First_E_guess


def VerifyConcavity(PotentialArray, First_E_guess):
    """Evaluates the concavity of the potential and returns its value: positive if the concavity is correct or negative if it is incorrect. To be positive, the smallest meeting
    point of an energy guess with the potential needs to have a negative derivative and the greatest meeting point needs to have a positive one. If the function finds no meeting point
    then tries a smaller energy guess and restarts the process


    Parameters:
    -----------
        PotentialArray (numpy.ndarray) : a numpy array that contains the potential value between 'x_V_min' and 'x_V_max' at every interval of length 'Division'

        First_E_guess (float) : the first energy guess that will be used in the Numerov algorithm

    Returns:
    --------
        concavity (str) : a string that indicates the global concavity of the potential. It can either be positive if it respects the condition or negative if it doesn't

    """

    i = 1
    #Continue while it doesn't find meeting points
    while i == 1:
        print('First Energy guess:', First_E_guess)
        index_min=list()
        index_max=list()

        #Tries to find meeting points and to compare them
        try:
            for i in range(0,len(PotentialArray)-2):

                #Gets all the points where the potential meets the E_verify value and filters them depending on their derivatives
                if PotentialArray[i] > First_E_guess and PotentialArray[i+1] < First_E_guess:
                    index_min.append(i)

                elif PotentialArray[i] < First_E_guess and PotentialArray[i+1] > First_E_guess:
                    index_max.append(i)

                elif PotentialArray[i] == First_E_guess:
                    if PotentialArray[i-1] > First_E_guess and PotentialArray[i+1] < First_E_guess:
                        index_min.append(i)

                    elif PotentialArray[i-1] < First_E_guess and PotentialArray[i+1] > First_E_guess:
                        index_max.append(i)

            #Defines the concavity value depending on
            print('index max: ',index_max)
            print('index_min: ',index_min)

            if (max(index_max) > max(index_min)) and (min(index_max) > min(index_min)):
                concavity = 'positive'
            else:
                concavity = 'negative'

        #If we are not able to compare the potential, we define a new energy guess
        except ValueError:
            First_E_guess = First_E_guess/2

        #If it is able to compare them, exit the loop
        else:
            i = 0

    return concavity,First_E_guess


def EvaluateOnePotential(position,potential):
    """Defines a function that evaluate the potential at a certain point x. This function will be vectorized with np.vectorize to evaluate the potential on a list of position [x1,x2,...]

    Parameters:
    -----------
        position (float) : a float that defines the x position where we want to evaluate the potential

        potential (str) : a string that defines the mathematical expression of the potential

    Returns:
    --------
        EvalPotential (float) : the potential value at the x position

    """

    x = position
    EvalPotential = eval(potential)

    return EvalPotential

def TranslationPotential(PositionPotential, PotentialArray):
    """Checks approximately where the minimum of the potential is and outputs the necessary translation in x and y to recenter the minimum at x=0 and y=0

    Parameters:
    -----------
        potential (str) : a string that defines the mathematical expression of the potential

    Returns:
    --------
        trans_x (float) : the necessary x translation to replace the minimum of the potential at x=0

        trans_y (float) : the necessary y translation to be sure that all the potential values are positive

    """

    # i) Gets the minimum value for the potential and the translation in y
    trans_y = PotentialArray.min()
    #index = float(np.where(PotentialArray==trans_y)[0])

    # ii) Defines the necessary translation in x
    #trans_x = x_min + (Div * index)
    #trans_x = PositionPotential[index]

    # iii) Translates the potential
    PotentialArray = PotentialArray - trans_y
    #PositionPotential = PositionPotential - trans_x

    #print('trans_x; ',trans_x)
    print('trans_y; ',trans_y)

    return PositionPotential, PotentialArray

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
    """Finds the minimal and maximal points where the energy that has been guessed is equal to the potential.

    Parameters:
    -----------
        E_guess: the guessed energy
        PotentialArray: a Numpy array that contains the potential for certain points
        PositionPotential: a Numpy array that contains the positions that correspond to the potential array

    Returns:
    --------
        MeetingPoints: a tuple of the smallest and biggest meeting point"""

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

def DetermineMinAndMax(MeetingPoints, E_guess, E_guess_try, PotentialArray, PositionPotential):
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
    i = 1
    iteration = 1
    end_program = False

    while i == 1:
        try:
            Position_min = MeetingPoints[0] - (MeetingPoints[1] - MeetingPoints[0])/1
            Position_max =  MeetingPoints[1] + (MeetingPoints[1] - MeetingPoints[0])/1
        except:
            E_guess = (E_guess + max([k for j,k in E_guess_try.values() if k < E_guess]))/2
            iteration += 1
            MeetingPoints = MeetingPointsPotential(E_guess, PotentialArray, PositionPotential)
            if iteration > 10:
                end_program = True
        else:
            i = 0
            print('MeetingPoint: ', MeetingPoints)
            print('min:',Position_min)
            print('max:',Position_max)

    return Position_min,Position_max, end_program

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
    VerificationTolerance = 'yes' if np.absolute(WaveFunction[-1][1]) < Tolerance else 'no'
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

#Define the wave funcions to plot, the lines corresponding to these wave function and the energy lines
def DefineWhatToPlot(WaveFunctionFound, EnergyLevelFound):

    # i) Determine the maximum energy to set the maximum value for the y axis
    y_max = 1.1*EnergyLevelFound[max(EnergyLevelFound)]
    Y_by_E_level = (y_max/(max(EnergyLevelFound)+2))

    # ii) For the wave function
    WavPlot = []
    for i in WaveFunctionFound.keys():
        x=[]
        y=[]
        for j in range(800,len(WaveFunctionFound[i])-800):
            if not (j > 7500 and np.absolute(WaveFunctionFound[i][j][1]) > (max(y)*0.07)):
                x.append(WaveFunctionFound[i][j][0])
                y.append(WaveFunctionFound[i][j][1])
        x = np.array(x)
        y = np.array(y)

        mult = (0.9 * Y_by_E_level)/(2 * y.max())
        y = (mult * y) + (Y_by_E_level * (i+1))
        WavPlot.append((x,y))

    # iii) Determines the min and max in x
    min_x = x.min()
    max_x = x.max()

    # iv) Get lines to where the wave function is centered
    WavLines = []
    for i in WaveFunctionFound.keys():
        Wav_line_y=[]
        for j in range(len(x)):
            Wav_line_y.append(Y_by_E_level * (i+1))
        WavLines.append((x,Wav_line_y))


    # v) get lines for all the Energy levels
    EnergyLines = []
    for i in WaveFunctionFound.keys():
        En_y = []
        for j in range(len(x)):
            En_y.append(EnergyLevelFound[i])
        EnergyLines.append((x,En_y))

    return y_max, min_x, max_x, WavPlot, WavLines, EnergyLines


#Draw the wave Functions, the energy levels and sets the axis limits
def DrawWaveFunction(y_max, min_x, max_x, WavPlot, WavLines, EnergyLines, PositionPotential, PotentialArray):

    ###################################################################################################
    # i) Define a new figure with two subplot: the energy levels and the corresponding wave function
    f,(En,Wav) = plt.subplots(1,2,sharey=True)

    #Set figure title
    f.suptitle("Schrödinger equation solutions",fontsize=20,fontweight='bold')

    ################################
    # ii) Draw the wave functions
    lines = [Wav.plot(x,y,'b',label=r"$Re(\psi(x))$",zorder=3)[0] for x,y in WavPlot]
    lines2 =  [Wav.plot(x,y,'m',label=r"$Im(\psi(x))$",zorder=3)[0] for x,y in WavPlot]

    for x,y in WavLines:
        Wav.plot(x,y,'k--',zorder=1)

    #Sets the axis limits
    Wav.axis([min_x, max_x, 0, y_max])

    #Draw the potential
    Wav.plot(PositionPotential, PotentialArray, 'r',label='Potential',zorder=2)

    ################################
    # iii) Draw the Energy levels
    i = 1
    for x,y in EnergyLines:
        PlotColor = cm.viridis(i/len(EnergyLines))
        En.plot(x,y,'--',color=PlotColor,label='E'+str(i),zorder=2)
        i+=1

    #Set the axis limit
    En.axis([min_x, max_x, 0, y_max])

    #Draw the potential
    En.plot(PositionPotential, PotentialArray, 'r',label='Potential',zorder=1)

    ####################################################
    # iv) Sets differents esthetic components

    #For the wave function set the title and the axis title
    Wav.set_xlabel(r'x ($a_0$)')
    Wav.set_title('Wave Function',fontsize=14)

    #Verify if the labels reappear multiple times and set legend for the wave function
    handles, labels = plt.gca().get_legend_handles_labels()
    newLabels, newHandles = [], []
    for handle, label in zip(handles, labels):
        if label not in newLabels:
            newLabels.append(label)
            newHandles.append(handle)
    leg1 = Wav.legend(newHandles, newLabels, fancybox=True, loc='upper right')
    leg1.get_frame().set_alpha(1)

    #Identify each wave function
    for i in range(len(EnergyLines)):
        Wav.text(((max_x - min_x) * 0.04) + min_x, WavLines[i][1][0] - (0.2 * (y_max/(len(EnergyLines)+2))), r'$\Psi_{%s}(x)$'%(i))

    #For the energy levels set the title, the axis title and the legend
    En.set_xlabel(r'x ($a_0$)')
    En.set_ylabel('Energy (Hartree)')
    En.set_title('Energy levels',fontsize=14)
    leg2 = En.legend(fancybox=True, loc='upper right')
    leg2.get_frame().set_alpha(1)

    #################################
    # v) Animate the wave function

    def init():
        for line in lines:
            line.set_data([], [])
        return lines

    def UpdateData(t):
        for j,line in enumerate(lines):
            x = WavPlot[j][0]
            y = ((WavPlot[j][1] - (WavLines[j][1][0]))  * np.cos(EnergyLines[j][1][0]*t/20)) + (WavLines[j][1][0])
            line.set_data(x,y)
        for j,line in enumerate(lines2):
            x = WavPlot[j][0]
            y = ((WavPlot[j][1] - (WavLines[j][1][0]))  * np.sin(EnergyLines[j][1][0]*t/20)) + (WavLines[j][1][0])
            line.set_data(x,y)

        return lines,lines2

    anim = animation.FuncAnimation(f, UpdateData, init_func=init, interval=20, blit=False)
    plt.show()

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

