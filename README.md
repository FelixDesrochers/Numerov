# Numerov

A python script that solves the one dimensional time-independent Schrodinger equation for bound states.  The script uses a Numerov method to solve the differential equation and displays the wanted energy levels and a figure with an approximate wave fonction for each of these energy levels.

<p align="center"><img alt="$$&#10;\frac{-\hbar}{2m}\frac{d^2}{dx^2}\psi(x) = E\psi(x) &#10;$$" src="https://rawgit.com/leegao/readme2tex (fetch/master/svgs/b39a5217b0c2dc2208c93980fc82ab24.svg?invert_in_darkmode" align=middle width="155.85124499999998pt" height="35.777445pt"/></p>


## Running

To run this code simply clone this repository and run the Numerov.py script with python (the numpy and matplotlib modules are required):
 
```
<img alt="$ git clone https://github.com/FelixDesrochers/Numerov/&#10;$" src="https://rawgit.com/leegao/readme2tex (fetch/master/svgs/03c2fb4726389b92a0c0c21969414b1f.svg?invert_in_darkmode" align=middle width="425.88430500000004pt" height="24.65759999999998pt"/> cd Numerov
<img alt="$ python Numerov.py&#10;```&#10;&#10;Then the program will ask you to enter the number of energy levels you want to display and the desired potential (make sure that the potential is centered at x=0):&#10;&#10;```&#10;$" src="https://rawgit.com/leegao/readme2tex (fetch/master/svgs/f1e30fdb082040289689d49153e14141.svg?invert_in_darkmode" align=middle width="700.2748499999999pt" height="118.35582pt"/> >> Which first energy levels do you want (enter an integer) : 4
<img alt="$ &gt;&gt; Potential (as a fonction of x): 3*(x^4)-2*(x^3)-6*(x^2)+x+5&#10;```&#10;&#10;Note: The programm may sometimes display less energy levels than what has been asked. To solve this problem modify the values of x_V_min and x_V_max in the parameters section of the Numerov.py script.&#10;&#10;&#10;## Examples&#10;&#10;### Harmonic Oscillator&#10;&#10;For instance, if we want the energy levels for the quantum harmonic oscillator we would run the following commands:&#10;&#10;```sh&#10;$" src="https://rawgit.com/leegao/readme2tex (fetch/master/svgs/7ea7f9b0c9b5b3b9b2515b8af8e3b08c.svg?invert_in_darkmode" align=middle width="756.9853499999999pt" height="276.16413pt"/> python Numerov.py
<img alt="$ &gt;&gt; Which first energy levels do you want (enter an integer) : 8&#10;$" src="https://rawgit.com/leegao/readme2tex (fetch/master/svgs/682b803e13a59383813b0c7a6b896df9.svg?invert_in_darkmode" align=middle width="433.391805pt" height="24.65759999999998pt"/> >> Potential (as a fonction of x): x**2
```

The program then displays the following figure:

<img src="/Examples/Harm_pot.gif?raw=true" width="1200" height="600" />

And the following energies (in hartree):

```
Energy level 0 : 0.707658207399
Energy level 1 : 2.12132034435
Energy level 2 : 3.5355339059
Energy level 3 : 4.94974746826
Energy level 4 : 6.36396103051
Energy level 5 : 7.77817459266
Energy level 6 : 9.19238815459
Energy level 7 : 10.6066017163
```

### Other examples

This program can also solve the Schrödinger equation for all sorts of unorthodox such as the double-well potential or the absolute value.

#### Double-Well Potential

For a "double-well potential" with the following form,

<p align="center"><img alt="$$&#10;V(x) = (x^4)+6\cdot(x^2)-9&#10;$$" src="https://rawgit.com/leegao/readme2tex (fetch/master/svgs/a2b6b57062ea23112daf3d3f3b53b4e0.svg?invert_in_darkmode" align=middle width="184.94189999999998pt" height="18.312359999999998pt"/></p>

we get the following results:

```sh
<img alt="$ python Numerov.py&#10;$" src="https://rawgit.com/leegao/readme2tex (fetch/master/svgs/19686c0a55a955512c1d0929696475da.svg?invert_in_darkmode" align=middle width="142.544325pt" height="22.831379999999992pt"/> >> Which first energy levels do you want (enter an integer) : 7
$ >> Potential (as a fonction of x): (x^4)+6*(x^2)-9
```

<img src="/Examples/Double_pot.gif?raw=true" width="1200" height="600" />

```
Energy level 0 : 2.35727297545
Energy level 1 : 2.35937176485
#Energy level 2 : 6.54881394689
Energy level 3 : 6.68442950475
Energy level 4 : 9.41561275062
Energy level 5 : 10.6784965326
Energy level 6 : 12.8773418447
```

#### Absolute Value

And for the absolute value (V(x)=|x|), we get:

<img src="/Examples/abs_value.gif?raw=true" width="1200" height="600" />

```
Energy level 0 : 0.81639999999
Energy level 1 : 1.85575743448
Energy level 2 : 2.57809582976
Energy level 3 : 3.24460762395
Energy level 4 : 3.82571482969
Energy level 5 : 4.38167123906
Energy level 6 : 4.89181971232
```

# Contributing

I am oppen to any improvement suggestion or contribution. If you wish to contribute to this repository just follow these simple steps:

1. Fork it (<https://github.com/yourname/yourproject/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request


# License
MIT - [http://alco.mit-license.org](http://alco.mit-license.org)

(See the LICENSE.md for more informations)
