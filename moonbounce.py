
# -*- coding: utf-8 -*-
# Copyright (C) Michael Coughlin and Christopher Stubbs(2015)
#
# skybrightness is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# skybrightness is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with skybrightness.  If not, see <http://www.gnu.org/licenses/>.

"""This module calculates ROLO from Kieffer and Stone [2005]. 
"""

import datetime, time
import numpy as np
import matplotlib
#matplotlib.rc('text', usetex=True)
#matplotlib.use('Agg')
#matplotlib.rcParams.update({'font.size': 16})
import matplotlib.pyplot as plt

def smooth(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

file = 'data/kieffer_stone_table4.txt'
data_out = np.loadtxt(file)

params = ['wavelength','a0','a1','a2','a3','b1','b2','b3','d1','d2','d3'] 

data = {}
ii = 0
for param in params:
    data[params[ii]] = data_out[:,ii]
    ii = ii + 1

file = 'data/kieffer_stone_eq8.txt'
data_out = np.loadtxt(file)
data_stone = {}
bands = ["u","g","r","i","z","y"]
for ii,band in zip(range(len(bands)),bands):
    data_stone[band] = {}
    data_stone[band]["phase"] = [2.0,5.0,10.0,45.0,90.0]
    data_stone[band]["omegaM"] = 6.4177*1e-5 # steradians
    data_stone[band]["wavelength"] = data_out[ii,0] # nm
    data_stone[band]["lunar"] = data_out[ii,1:6] # microW / m^2
    data_stone[band]["solar"] = data_out[ii,6] * 1e6 # microW/m^2 nm

    data_stone[band]["A"] = data_stone[band]["lunar"] * np.pi/(data_stone[band]["omegaM"]*data_stone[band]["solar"])
    data_stone[band]["magA"] = -2.5 * np.log10(data_stone[band]["A"])

    #print data_stone[band]["magA"]    
    #print data_stone[band]["Ak"]

c1 = 0.00034115 
c2 = -0.0013425
c3 = 0.00095906
c4 = 0.00066229

p1=4.06054
p2=12.8802
p3=-30.5858
p4=16.7498

## spline lunar coefficients onto 1 nm spacing, with smoothing
nm = np.arange(350,1060)
rnm = data['wavelength']

smoothparam = 1.0
a0=np.interp(nm,rnm,smooth(data['a0'],smoothparam))
a1=np.interp(nm,rnm,smooth(data['a1'],smoothparam))
a2=np.interp(nm,rnm,smooth(data['a2'],smoothparam))
a3=np.interp(nm,rnm,smooth(data['a3'],smoothparam))
b1=np.interp(nm,rnm,smooth(data['b1'],smoothparam))
b2=np.interp(nm,rnm,smooth(data['b2'],smoothparam))
b3=np.interp(nm,rnm,smooth(data['b3'],smoothparam))
d1=np.interp(nm,rnm,smooth(data['d1'],smoothparam))
d2=np.interp(nm,rnm,smooth(data['d2'],smoothparam))
d3=np.interp(nm,rnm,smooth(data['d3'],smoothparam))

#a0[0:10] = data['a0'][0:10]

plotName = "plots/moonbounce.png"
plt.figure()
plt.plot(nm,a0,'bo-',label='a0')
plt.plot(nm,a1,'gx-',label='a1')
plt.plot(nm,a2,'r*-',label='a2')
plt.xlabel("Wavelength [microns]")
plt.ylabel("Coefficients")
plt.legend()
plt.savefig(plotName)
plotName = "plots/moonbounce.eps"
plt.savefig(plotName)
plt.close()

## make boolean masks for passbands
umask = np.zeros((len(nm),))
gmask = np.zeros((len(nm),))
rmask = np.zeros((len(nm),))
imask = np.zeros((len(nm),))
zmask = np.zeros((len(nm),))
ymask = np.zeros((len(nm),))

umask[np.intersect1d(np.where(nm>=320),np.where(nm<385))] = 1.0
gmask[np.intersect1d(np.where(nm>=401),np.where(nm<550))] = 1.0
rmask[np.intersect1d(np.where(nm>=562),np.where(nm<695))] = 1.0
imask[np.intersect1d(np.where(nm>=695),np.where(nm<844))] = 1.0
zmask[np.intersect1d(np.where(nm>=826),np.where(nm<920))] = 1.0
ymask[np.intersect1d(np.where(nm>=950),np.where(nm<1058))] = 1.0

#umask=(nm>=320)&(nm<385)
#gmask=(nm>=401)&(nm<550)
#rmask=(nm>=562)&(nm<695)
#imask=(nm>=695)&(nm<844)
#zmask=(nm>=826)&(nm<920)
#ymask=(nm>=950)&(nm<1058)

## compute reflectance spectra vs. lunar phase angle 
## note they have a really nasty combination of degrees and radians!
smallest_angle = 2.0
smallest_rad=smallest_angle*np.pi/180
phase=np.arange(smallest_rad,0.95*np.pi,0.01)
phasedeg=180*phase/np.pi

magmatrix = np.zeros((len(phase),len(nm)))
Rmatrix = np.zeros((len(phase),len(nm)))

plotName = "plots/Amag.png"
plt.figure()

for loopcount in xrange(len(phase)): ## evaluate this one phase angle at a time
    g=phase[loopcount]
    gdeg=phasedeg[loopcount] # convert to degrees for their exponents!
    ## next expression comes from lunar irradiance paper, uses their fit
    ## values, smoothed and splined
    lnA = a0+a1*g+a2*g*g+a3*g*g*g
    lnA = lnA + b1*(-g)+b2*(-g)*(-g)*(-g)+b3*(-g)*(-g)*(-g)*(-g)*(-g)
    lnA = lnA + d1*np.exp(-gdeg/p1)+d2*np.exp(-gdeg/p2)+d3*np.cos((np.pi/180)*(gdeg*180/np.pi-p3)/p4)
    ## now convert to magnitudes
    Rmag=-2.5*np.log10(np.exp(lnA))
    ## this is an implicit loop over wavelengths
    magmatrix[loopcount,:]=Rmag #% this is reflection in magnitudes
    Rmatrix[loopcount,:]=np.exp(lnA) #% this is reflection matrix in linear units. First index is angle, second is lambda.

    if loopcount in [0, 50, 100]:
        plt.plot(nm,-2.5*np.log10(np.exp(lnA)),label='%.0f'%gdeg)
plt.xlabel("Wavelength [microns]")
plt.ylabel("-2.5 log10 (A)")
plt.legend()
plt.savefig(plotName)
plotName = "plots/Amag.eps"
plt.savefig(plotName)
plt.close()

## compute differences from full moon by adding up A(lambda) across passbands.

## first, figure out color of moon compared to sun, at full moon. This is
## the difference in total reflection, compared to no wavelength dep.

Fullmoon=Rmatrix[0,:] # reflectance(lambda) at smallest angle we consider
ufull=np.sum(Fullmoon*umask)/np.sum(umask)
gfull=np.sum(Fullmoon*gmask)/np.sum(gmask)
rfull=np.sum(Fullmoon*rmask)/np.sum(rmask)
ifull=np.sum(Fullmoon*imask)/np.sum(imask)
zfull=np.sum(Fullmoon*zmask)/np.sum(zmask)
yfull=np.sum(Fullmoon*ymask)/np.sum(ymask)

du = np.zeros((len(phase),))
dg = np.zeros((len(phase),))
dr = np.zeros((len(phase),))
di = np.zeros((len(phase),))
dz = np.zeros((len(phase),))
dy = np.zeros((len(phase),))

## compute change in magnitude as a function of phase, from full moon value
for loopcount in xrange(len(phase)):
   utemp=sum(Rmatrix[loopcount,:]*umask)/np.sum(umask)
   du[loopcount]=-2.5*np.log10(utemp)
   
   gtemp=sum(Rmatrix[loopcount,:]*gmask/np.sum(gmask))
   dg[loopcount]=-2.5*np.log10(gtemp)
   
   rtemp=sum(Rmatrix[loopcount,:]*rmask/np.sum(rmask))
   dr[loopcount]=-2.5*np.log10(rtemp)
   
   itemp=sum(Rmatrix[loopcount,:]*imask/np.sum(imask))
   di[loopcount]=-2.5*np.log10(itemp)
   
   ztemp=sum(Rmatrix[loopcount,:]*zmask/np.sum(zmask))
   dz[loopcount]=-2.5*np.log10(ztemp)
   
   ytemp=sum(Rmatrix[loopcount,:]*ymask/np.sum(ymask))
   dy[loopcount]=-2.5*np.log10(ytemp)

plotName = "plots/filtermag.png"
plt.figure()
plt.plot(phasedeg,du,'b-',label='u')
plt.plot(phasedeg,dg,'g-',label='g')
plt.plot(phasedeg,dr,'r-',label='r')
plt.plot(phasedeg,di,'c-',label='i')
plt.plot(phasedeg,dz,'m-',label='z')
plt.plot(phasedeg,dy,'k-',label='y')
plt.xlabel("Phase [degrees]")
plt.ylabel("Magnitudes")
plt.legend(loc=4)
plt.savefig(plotName)
plotName = "plots/filtermag.eps"
plt.savefig(plotName)
plt.close()

print "ROLO Model"
phases = [2, 5, 10, 45, 90]
for ii in xrange(len(phases)):
   index = np.argmin(np.abs(phasedeg - phases[ii]))
   print '%5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f'%(phases[ii],du[index],dg[index],dr[index],di[index],dz[index],dy[index])

print "Lunar irradiances"
phases = [2, 5, 10, 45, 90]
for ii in xrange(len(phases)):
   print '%5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f'%(phases[ii],data_stone["u"]["magA"][ii],data_stone["g"]["magA"][ii],data_stone["r"]["magA"][ii],data_stone["i"]["magA"][ii],data_stone["z"]["magA"][ii],data_stone["y"]["magA"][ii])
