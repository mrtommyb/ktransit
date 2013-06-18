ktransit
========

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
