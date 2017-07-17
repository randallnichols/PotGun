'''  file = potgun.py
15Jul17 - initial load to git copied from 38T13.py

Motion of projectile driven by hot injected gas in a smooth barrel.  
Output is a text table of velocity and position in the barrel
over specified time intervals terminating at the muzzle. 
'''

import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import lagrange

mw = 85.0   # Munition weight in grams
bdi = 1.5   # bore diameter in inches
bli = 39.37*2.0 # barrel length in inches 39.37
xi = 1.75   # empty barrel length in inches
bf = 3400.0 # barrel friction - 0.5 psig backpressure equivalent
gvl = 0.448*10.0  # injected gas volume in liters at STP
gms = 40.0  # gas injection duration in milliseconds
gtf = 400.0 # gas injection temperature in deg F

gc = 9.80665   # acceleration of gravity in m/s2 9.80665
ps = 100000.0  # standard pressure in Pascals = 14.5 psia = 0.9869 atm
ts = 273.15    # standard temperature in deg K = 32 deg F
ta = 298.15    # ambient temperature in deg K = 20 deg C = 68 deg F
R = 8.3144621  # ideal gas constant in m3*Pa/K*mol
cf1 = 12.0 / 0.3048  # meters to inches conversion factor 39.37007874015748
cf2 = 0.0254**3  # in3 to m3 conversion 1.6387064e-05 1/x = 61023.74409473229
cf3 = 6894.7573  # Pascal to psi ratio

mkg = mw/1000.0  # Munition weight in kilograms 0.085
bd = bdi / cf1   # bore diameter in meters 0.0381
bs = (np.pi/4.0) * bd**2  # bore cross section area in m2 0.00114009182796937
bl = bli/cf1     # barrel length in meters 0.508
bv = bs*bl       # barrel volume in m3 0.00057916664860844 = 0.5791 liters
gvm = gvl/22.4   # moles of gas injected 0.0625
gs1 = gms/1000.0
gvr = gvm/gs1    # gas injection rate in moles per second 1.25
gtk = (gtf - 32.0)/1.8 + ts  # gas injection temp in deg K 477.5944
x0 = xi/cf1      # empty barrel length under projectile in meters 
v0 = bs*x0       # initial empty volume behind projectile 2.17187493228e-05
n0 = ps*v0/(R*gtk) # initial moles air in empty space 0.0005469421436
psig = []
x0

def derv(yint, t, xx):
    if t < gs1:
        nm = gvr*t + n0 # number of moles injected at time t
    else:
        nm = gvm
    bpa = nm*R*gtk/(bs*(yint[1]+x0)) # barrel pressure in Pa abs
    bpg = bpa - ps - bf # effective gage barrel pressure in Pa
    ode0 = bpg*bs / (mkg*gc) # acceleration in m/s2
    ode1 = yint[0]
    return [ode0, ode1]

y0 = np.array([0.0, 0.0,])
t = np.linspace(0, 0.07, 71)

yint = odeint(derv, y0, t, args=(22,))

print('  time    gas in     Press     Accel     Vel      Pos')
print('   ms      mol/s      psig      m/s2     m/s       m')

n = 0
mpng = 0.0
for x in yint:
    if x[1] > bl: break # if barrel length is exceeded
    if t[n] < gs1: # if gas injection is ON
        gr = gvr # gas injection rate
        nm = gvr*t[n] + n0 # total moles injected at time t
    else: # if gas injection is OFF
        gr = 0.0 # gas injection is zero
        nm = gvm # total moles is the maximum

    bpna = nm*R*gtk/(bs*(x[1]+x0)) # barrel pressure in Pa abs
    bpng = (bpna - ps) # barrel pressure in Psig
    if bpng > mpng: mpng = bpng
    acc = bpng*bs / (mkg*gc) # acceleration in m/s2
    sg3 = "{:5.1f} {:9.5f} {:10.3f}".format(t[n]*1000.0, gr, bpng/cf3)
    sg4 = "{:10.2f} {:8.4f} {:8.4f}".format(acc, x[0], x[1])
    print(sg3 + sg4)
    n += 1
    
lgx = [yint[n-2,1], yint[n-1,1], yint[n,1]] # last 3 x positions
lgv = [yint[n-2,0], yint[n-1,0], yint[n,0]] # last 3 velocities
lgt = [t[n-2], t[n-1], t[n]] # last 3 times 
pfv = lagrange(lgx, lgv) # velocity polynomial fit
exv = pfv(bl)  # projectile muzzle velocity in m/s
pft = lagrange(lgx, lgt) # muzzle exit time polynomial fit
ext = pft(bl)*1000  # projectile muzzle exit time in milliseconds
sg5 = "{:7.3f} {:7.5f} {:10.3f}".format(ext, gvl, mpng/cf3)
sg6 = "{:8.4f} {:8.4f}".format(exv, bl)
print(sg5 + '           ' + sg6)

'''  Example Output
  time    gas in     Press     Accel     Vel      Pos
   ms      mol/s      psig      m/s2     m/s       m
  0.0   3.12500      0.000      0.00   0.0000   0.0000
  1.0   3.12500     35.452    334.32   0.1673   0.0001
  2.0   3.12500     70.182    661.82   0.6665   0.0004
  3.0   3.12500    102.606    967.60   1.4838   0.0015
  4.0   3.12500    130.606   1231.63   2.5877   0.0035
  5.0   3.12500    152.255   1435.79   3.9270   0.0067
  6.0   3.12500    166.571   1570.80   5.4361   0.0114
  7.0   3.12500    173.810   1639.06   7.0462   0.0177
  8.0   3.12500    175.191   1652.08   8.6958   0.0255
  9.0   3.12500    172.335   1625.15  10.3371   0.0350
 10.0   3.12500    166.785   1572.81  11.9377   0.0462
 11.0   3.12500    159.759   1506.55  13.4781   0.0589
 12.0   3.12500    152.110   1434.42  14.9488   0.0731
 13.0   3.12500    144.382   1361.55  16.3467   0.0888
 14.0   3.12500    136.895   1290.94  17.6727   0.1058
 15.0   3.12500    129.816   1224.19  18.9299   0.1241
 16.0   3.12500    123.223   1162.01  20.1226   0.1436
 17.0   3.12500    117.132   1104.58  21.2555   0.1643
 18.0   3.12500    111.534   1051.78  22.3333   0.1861
 19.0   3.12500    106.398   1003.35  23.3605   0.2090
 20.0   0.00000     99.365    937.03  24.3414   0.2328
 21.0   0.00000     90.018    848.89  25.2333   0.2576
 22.0   0.00000     81.839    771.76  26.0428   0.2833
 23.0   0.00000     74.652    703.98  26.7799   0.3097
 24.0   0.00000     68.310    644.17  27.4534   0.3368
 25.0   0.00000     62.688    591.16  28.0706   0.3646
 26.0   0.00000     57.683    543.96  28.6377   0.3929
 27.0   0.00000     53.209    501.77  29.1601   0.4218
 28.0   0.00000     49.192    463.89  29.6426   0.4512
 29.0   0.00000     45.572    429.75  30.0892   0.4811
 29.888                               30.4587   0.5080
 
'''
