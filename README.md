# Numerov

A python script that solves the one dimensional time-independent Schrodinger equation for bound states.  The script uses a Numerov method to solve the differential equation and displays the wanted energy levels and a figure with an approximate wave fonction for each of these energy levels.

## Description


## Examples

### Harmonic oscillator

```sh
$ pyhton Numerov.py
$ Which first energy levels do you want (enter an integer) : 8
$ Potential (as a fonction of x): x**2
```

<img src="/Examples/Harm_pot.gif?raw=true" width="1100" height="700" />

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

```sh
$ pyhton Numerov.py
$ Which first energy levels do you want (enter an integer) : 7
$ Potential (as a fonction of x): (x^4)+6*(x^2)-9
```

<img src="/Examples/Double_pot.gif?raw=true" width="1100" height="700" />

```
Energy level 0 : 2.35727297545
Energy level 1 : 2.35937176485
Energy level 2 : 6.54881394689
Energy level 3 : 6.68442950475
Energy level 4 : 9.41561275062
Energy level 5 : 10.6784965326
Energy level 6 : 12.8773418447
```
