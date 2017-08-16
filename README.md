# Numerov

A python script that solves the one dimensional time-independent Schrodinger equation for bound states.  The script uses a Numerov method to solve the differential equation and displays the wanted energy levels and a figure with an approximate wave fonction for each of these energy levels.



## Running

To run this code simply clone this repository and run the Numerov.py script with python (the numpy and matplotlib modules are required):
 
```
$ git clone https://github.com/FelixDesrochers/Numerov/
$ cd Numerov
$ python Numerov.py
```

Then the program will ask you to enter the number of energy levels you want to display and the desired potential (make sure that the potential is centered at x=0):

```
$ >> Which first energy levels do you want (enter an integer) : 4
$ >> Potential (as a fonction of x): 3*(x^4)-2*(x^3)-6*(x^2)+x+5
```

Note: The programm may sometimes display less energy levels than what has been asked. To solve this problem modify the values of x_V_min and x_V_max in the parameters section of the Numerov.py script.


## Examples

### Harmonic Oscillator

For instance, if we want the energy levels for the quantum harmonic oscillator we would run the following commands:

```sh
$ python Numerov.py
$ >> Which first energy levels do you want (enter an integer) : 8
$ >> Potential (as a fonction of x): x**2
```

The program then displays the following figure:

<img src="/Examples/Harm_pot.gif?raw=true" width="1200" height="600" />

And the following energies:

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

### Other Examples

This program can also solve the Schrödinger equation for all sorts of unorthodox such as the double-well potential or the absolute value.

#### Double-Well Potential

For a "double-well potential" with the following input,

```sh
$ python Numerov.py
$ >> Which first energy levels do you want (enter an integer) : 7
$ >> Potential (as a fonction of x): (x^4)-6*(x^2)+9
```

we get the following results:

<img src="/Examples/Double_pot.gif?raw=true" width="1200" height="600" />

```
Energy level 0 : 2.35727297545
Energy level 1 : 2.35937176485
Energy level 2 : 6.54881394689
Energy level 3 : 6.68442950475
Energy level 4 : 9.41561275062
Energy level 5 : 10.6784965326
Energy level 6 : 12.8773418447
```

### Absolute value potential 

And for the absolute value potential (V(x)=|x|), we get:

<img src="/Examples/abs_value.gif?raw=true" width="1200" height="600" />

```
Energy level 0 : 0.81639999999
Energy level 1 : 1.85575743448
Energy level 2 : 2.57809582976
Energy level 3 : 3.24460762395
Energy level 4 : 3.82571482969
Energy level 5 : 4.38167123906
Energy level 6 : 4.89181971232
Energy level 7 : 5.38661378006
Energy level 8 : 5.8513002713
Energy level 9 : 6.30526300457
```

# Algorithm

All informations about the used algorithm are described in the explain_algorithm.pdf file.


# Contributing

I am open to any improvement suggestion or contribution. If you wish to contribute to this repository just follow these simple steps:

1. Fork it (<https://github.com/yourname/yourproject/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request


# License
MIT - [http://alco.mit-license.org](http://alco.mit-license.org)

(See the LICENSE.md for more informations)
