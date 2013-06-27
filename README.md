ktransit
========
**A simple exoplanet transit modeling tool**

This package contains routines creating and optionally fitting a transiting planet model.
The model is based on the work of `Mandel & Agol (2002) <http://iopscience.iop.org/1538-4357/580/2/L171/fulltext/16756.text.html>`_. 
`Goodman & Weare (2010) <http://cims.nyu.edu/~weare/papers/d13.pdf>`_.


The basic module, **ktransit**

```
import ktransit
import matplotlib.pyplot as plt

M = ktransit.LCModel()
M.add_star()
M.add_planet()
M.add_data()

tmod = M.transitmodel
plt.plot(M.time,tmod)
```

##OR
```
time,earthlc = ktransit.give_me_earth()
plt.plot(time,earthlc)
```
