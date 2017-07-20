'''This script solves the 1 dimensional time-independant Schrodinger equation for any given potential. Its takes a potentiel as an entry and outputs the wanted energy level and a figure with the potentiel and the wave fonctions corresponding to the energy level that have been specicied.

The one dimensionl time independent Schrodinger is as follow:

    \frac{-\hbar^2}{2*m}*\frac{d}{dx} \psi(x) +  V(x) \psi(x) = E * \psi(x)

To solve this diffenrential equation a Numerov method wil be used. This method aims to solve any equation of the type : \frac{d}{dx}y(x) + s(x)*x =0
