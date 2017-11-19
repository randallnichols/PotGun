'''  file = potgun.py
Motion of projectile driven by hot gas injected into a smooth barrel.  
Output is a text table of velocity and position in the barrel
over specified time intervals terminating at the muzzle.

Calculations use metric units - meters, kilograms, seconds, and degrees K
Input parameters are converted to the above units

Calculation strategy:
Known at time t the number of moles n, volume V and Temperature T:
   Solve ideal gas equation: PV = nRT for pressure P (in Pascals)
Calculate force on projectile F = pressure * bore area (in Newtons)
Calculate acceleration A = force / mass of projectile (in meters/s^2)
Integrate(dA / dt) to get velocity in meters/sec (yint[t,0])
Integrate(dV / dt) to get position in barrel in meters (yint[t,1])

15Jul17 - initial load to git copied from 38T13.py
19Nov17 - complete code revamp and fix calculation error
'''

import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import lagrange

# input parameters
mw = 85.0   # Munition weight in grams
bdi = 1.5   # bore diameter in inches
bli = 40.0  # barrel length in inches
xi = 4.0    # munition starting point in inches (void space behind)
bfi = 0.5   # barrel friction in psig equivalent
gvl = 3.36  # total injected gas volume in liters at STP
gms = 30.0  # gas injection duration in milliseconds
gtf = 400.0 # gas injection temperature in deg F

# constants and units conversion factors
ps = 100000.0  # ambient pressure in Pascals = 14.5 psia = 0.9869 atm
ts = 273.15    # standard temperature in deg K = 32 deg F
ta = 298.15    # ambient temperature in deg K = 20 deg C = 68 deg F
R = 8.3144621  # ideal gas constant in m3*Pa/K*mol
cf1 = 12.0 / 0.3048    # meters to inches conversion factor = 39.37007874
cf2 = 1.0 / 0.0254**3  # in3 to m3 conversion = 61023.74409
cf3 = 6894.7573  # Pascal to psi conversion factor

# calculate parameters and convert to metric
mkg = mw/1000.0  # Munition weight in kilograms
bd = bdi / cf1   # bore diameter in meters
bs = (np.pi/4.0) * bd**2  # bore cross section area in m2
bl = bli/cf1     # barrel length in meters
gvm = gvl/22.4   # total moles of gas to be injected
gs1 = gms/1000.0 # gas injection time in seconds
gvr = gvm/gs1    # gas injection rate in moles per second
gtk = (gtf - 32.0)/1.8 + ts  # gas injection temp in deg K
x0 = xi/cf1      # empty barrel length behind projectile in meters 
v0 = bs*x0       # initial empty volume behind projectile in m3
n0 = ps*v0/(R*gtk) # initial moles air in empty space
bf = bfi*cf3     # barrel friction in Pascals equivalent

# setup ordinary differential equations (ode's)
def derv(yint, t):
    if t < gs1:              # if gas is still being injected
        nm = gvr*t + n0      # number of moles injected at time t
    else:                    # if gas is no longer being injected
        nm = gvm             # maximum moles already injected
    bpa = nm*R*gtk/(bs*yint[1]) # barrel pressure in Pa abs
    bpg = bpa - ps - bf      # effective gage barrel pressure in Pa
    ode0 = bpg*bs / mkg      # acceleration in m/s2 (a = force/mass)
    ode1 = yint[0]           # velocity in m/s - integral of acceleration
    return [ode0, ode1]      # return two ordinary differentials (ode's)

y0 = np.array([0.0, x0,])       # initial values of velocity and position
t = np.arange(0.0, 0.04, 0.001) # sequence of time points for solutions

yint = odeint(derv, y0, t)      # find solutions to ode's at time points t

# print output header
print('  time    gas in     Press     Accel     Vel      Pos')
print('   ms      mol/s      psig      m/s2     m/s       m')

n = 0              # initialize output table line counter
mpng = 0.0         # initialize barrel pressure
for x in yint:     # for every solution in yint (every time point in t)
    if x[1] > bl:  # if current position is beyond barrel length
        break      # stop iterating and go to find muzzle conditions
    if t[n] < gs1: # if gas injection is ON
        gr = gvr   # gas injection rate
        nm = gvr*t[n] + n0 # caclulate total moles in barrel at time t
    else:          # if gas injection is OFF
        gr = 0.0   # gas injection is zero
        nm = gvm   # total moles is the maximum that was injected

    bpna = nm*R*gtk/(bs*x[1])      # barrel pressure in Pa abs at time t
    bpng = (bpna - ps)             # barrel pressure in Pa gauge at time t
    if bpng > mpng:                # if barrel pressure is rising:
        mpng = bpng                # save pressure as maximum
    acc = bpng*bs/mkg              # acceleration in m/s2
    # format line of output table
    sg3 = "{:5.1f} {:9.3f} {:10.1f}".format(t[n]*1000.0, gr, bpng/cf3)
    sg4 = "{:10.1f} {:8.3f} {:8.3f}".format(acc, x[0], x[1])
    print(sg3 + sg4)               # print a line of output table
    n += 1                         # increment line counter

# after end of barrel has be detected, interpolate muzzle parameters
lgx = [yint[n-2,1], yint[n-1,1], yint[n,1]] # last 3 x positions
lgv = [yint[n-2,0], yint[n-1,0], yint[n,0]] # last 3 velocities
lgt = [t[n-2], t[n-1], t[n]] # last 3 times
pfv = lagrange(lgx, lgv)     # create vel vs x interpolation object
exv = pfv(bl)                # interpolate velocity at muzzle in m/s
pft = lagrange(lgx, lgt)     # create time vs x interpolation object
ext = pft(bl)*1000           # interpolate muzzle exit time in milliseconds

# format muzzle parameter output
sg5 = "{:6.2f} {:19.1f}".format(ext, mpng/cf3)
sg6 = "{:8.3f} {:8.3f}".format(exv, bl)
# print muzzle parameter header
print('  EXIT TIME          MAX P           VELOCITY     POS')
# print muzzle parameters
print(sg5 + '           ' + sg6)

'''  Example Output:

  time    gas in     Press     Accel     Vel      Pos
   ms      mol/s      psig      m/s2     m/s       m
  0.0     4.464        0.0       0.0    0.000    0.102
  1.0     4.464       22.1    2042.1    0.978    0.102
  2.0     4.464       42.9    3968.3    3.955    0.104
  3.0     4.464       60.1    5559.0    8.710    0.110
  4.0     4.464       71.5    6607.8   14.797    0.122
  5.0     4.464       76.4    7064.5   21.634    0.140
  6.0     4.464       76.2    7047.0   28.675    0.165
  7.0     4.464       72.9    6738.3   35.539    0.198
  8.0     4.464       68.1    6293.9   42.015    0.236
  9.0     4.464       62.9    5812.5   48.022    0.281
 10.0     4.464       57.8    5345.9   53.553    0.332
 11.0     4.464       53.2    4916.4   58.634    0.388
 12.0     4.464       49.0    4531.1   63.308    0.449
 13.0     4.464       45.3    4189.2   67.618    0.515
 14.0     4.464       42.0    3887.1   71.607    0.585
 15.0     4.464       39.1    3620.2   75.312    0.658
 16.0     4.464       36.6    3383.8   78.765    0.735
 17.0     4.464       34.3    3173.7   81.996    0.815
 18.0     4.464       32.3    2986.2   85.028    0.899
 19.0     4.464       30.5    2818.1   87.882    0.985
 20.0     0.000       27.5    2540.0   90.577    1.075
 21.0     0.000       24.2    2234.6   92.914    1.166
 22.0     0.000       21.3    1968.0   94.967    1.260
 23.0     0.000       18.8    1734.1   96.769    1.356
 24.0     0.000       16.5    1527.7   98.351    1.454
  EXIT TIME          MAX P           VELOCITY     POS
 24.71                76.4             99.357    1.524

https://www.sensorsone.com/pressure-and-area-to-force-calculator/
'''
