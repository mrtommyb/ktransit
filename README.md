ktransit
========
This package contains routines creating and optionally fitting a transiting planet model.
The model is based on the work of Mandel and Agol (2002) and is designed to be very simple
to use.

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
