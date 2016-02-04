# -*- coding: utf-8 -*-
# Copyright (C) Michael Coughlin and Christopher Stubbs (2015)
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

"""This module provides example methods to calculate the lunar contribution
to sky brightness using calculations from Coughlin, Stubbs, and Claver [2015].
"""

import datetime, time
import numpy as np
import matplotlib
matplotlib.rc('text', usetex=True)
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 16})
import matplotlib.pyplot as plt
import ephem

# Inputs
latitude = '-30:15:06.37'
longitude = '-70:44:17.50'
elevation = 2552.0
passband1 = 'u'
passband2 = 'g'
passband3 = 'r'
ra_obs = 12.0 # Hours
dec_obs = -60.0 # Degrees
time_of_observation = "2015/04/30 18:00:00"

deltats = np.arange(0,7*86400,60)

tts = []
total_mags_1 = []
total_mags_2 = []
total_mags_3 = []


for deltat in deltats:

    utc_date = datetime.datetime.strptime(time_of_observation, "%Y/%m/%d %H:%M:%S")
    utc_date = utc_date + datetime.timedelta(seconds=deltat)

    print utc_date

    # Where is the moon?
    moon = ephem.Moon()
    moon.compute(utc_date)
    ra_moon = (24/(2*np.pi))*float(repr(moon.ra))
    dec_moon = (180/np.pi)*float(repr(moon.dec))
    
    # Where is the sun?
    sun = ephem.Sun()
    sun.compute(utc_date)
    #ra_sun = (24/(2*np.pi))*float(repr(sun.ra))
    #dec_sun = (180/np.pi)*float(repr(sun.dec))

    # Coverting both target and moon ra and dec to radians
    ra1 = float(repr(moon.ra))
    ra2 = ra_obs * ((2*np.pi)/24)
    d1 = float(repr(moon.dec))
    d2 = dec_obs * ((2*np.pi)/360)
    
    # Calculate angle between target and moon
    cosA = np.sin(d1)*np.sin(d2) + np.cos(d1)*np.cos(d2)*np.cos(ra1-ra2)
    angle = np.arccos(cosA)*(360/(2*np.pi))
    print "Angle between moon and target: %.5f"%(angle)
    
    # Establish the location of the telescope
    telescope = ephem.Observer()
    telescope.lat = latitude
    telescope.long = longitude
    telescope.elevation = elevation
    telescope.date = utc_date
    
    # Determine altitude and azimuth of the target
    star = ephem.FixedBody()
    star._ra  = ra2
    star._dec = d2
    star.compute(telescope)
    alt_target = float(repr(star.alt)) * (360/(2*np.pi))
    az_target = float(repr(star.az)) * (360/(2*np.pi))
    print "Altitude / Azimuth of target: %.5f / %.5f"%(alt_target,az_target)
    
    # Determine altitude and azimuth of the moon
    star._ra  = ra1
    star._dec = d1
    star.compute(telescope)
    alt_moon = float(repr(star.alt)) * (360/(2*np.pi))
    az_moon = float(repr(star.az)) * (360/(2*np.pi))
    print "Altitude / Azimuth of moon: %.5f / %.5f"%(alt_moon,az_moon)
    
    # Determine altitude and azimuth of the sun
    star._ra  = sun.ra
    star._dec = sun.dec
    star.compute(telescope)
    alt_sun = float(repr(star.alt)) * (360/(2*np.pi))
    az_sun = float(repr(star.az)) * (360/(2*np.pi))
    print "Altitude / Azimuth of sun: %.5f / %.5f"%(alt_sun,az_sun)

    # Moon phase data (from Coughlin, Stubbs, and Claver Table 2) 
    moon_phases = [2,10,45,90]
    moon_data = {'u':[2.60,3.05,4.06,5.52],
                 'g':[2.36,2.78,3.77,5.19],
                 'r':[2.10,2.50,3.45,4.84],
                 'i':[1.92,2.31,3.23,4.59],
                 'z':[2.17,2.57,3.53,4.92],
                 'y':[1.72,2.12,2.91,4.31]}
    
    # Determine moon data for this phase
    moon_data_passband1 = moon_data[passband1]
    delta_mag_1 = np.interp(moon.moon_phase,moon_phases,moon_data_passband1)
    moon_data_passband2 = moon_data[passband2]
    delta_mag_2 = np.interp(moon.moon_phase,moon_phases,moon_data_passband2)
    moon_data_passband3 = moon_data[passband3]
    delta_mag_3 = np.interp(moon.moon_phase,moon_phases,moon_data_passband3)    

    # Fits to solar sky brightness (from Coughlin, Stubbs, and Claver Table 4) 
    sun_data = {'u':[88.5,-0.5,-0.5,0.4],
                'g':[386.5,-2.2,-2.4,0.8],
                'r':[189.0,-1.4,-1.1,0.8],
                'i':[164.8,-1.5,-0.7,0.6],
                'z':[231.2,-2.8,-0.7,1.4],
                'zs':[131.1,-1.4,-0.5,0.2],
                'y':[92.0,-1.3,-0.2,0.9]}
    
    # Determine sun data for this phase
    sun_data_passband = sun_data[passband1]
    flux = sun_data_passband[0] + sun_data_passband[1]*angle +\
           sun_data_passband[2]*alt_target + sun_data_passband[3]*alt_moon
    flux = flux* (10**11)
    #flux_mag = -2.5 * np.log10(flux)
    flux_mag_1 = flux

    sun_data_passband = sun_data[passband2]
    flux = sun_data_passband[0] + sun_data_passband[1]*angle +\
           sun_data_passband[2]*alt_target + sun_data_passband[3]*alt_moon
    flux = flux* (10**11)
    #flux_mag = -2.5 * np.log10(flux)
    flux_mag_2 = flux

    sun_data_passband = sun_data[passband3]
    flux = sun_data_passband[0] + sun_data_passband[1]*angle +\
           sun_data_passband[2]*alt_target + sun_data_passband[3]*alt_moon
    flux = flux* (10**11)
    #flux_mag = -2.5 * np.log10(flux)
    flux_mag_3 = flux

    # Determine total magnitude contribution
    total_mag_1 = delta_mag_1 + flux_mag_1
   
    print "Sun-> Moon conversion: %.5f"%delta_mag_1
    print "Sky brightness contribution: %.5f"%flux_mag_1
    print "Total magnitude reduction: %.5f"%total_mag_1
    
    print "" 

    if alt_moon < 0:
        total_mag_1 = 0
        total_mag_2 = 0
        total_mag_3 = 0
    if alt_target < 0:
        continue
    if alt_sun > 0:
        continue

    tts.append(deltat)
    total_mags_1.append(flux_mag_1)
    total_mags_2.append(flux_mag_2)
    total_mags_3.append(flux_mag_3)

tts = np.array(tts)

plotName = "plots/delta_mag.png"
plt.figure()
plt.plot(tts/86400.0,total_mags_1,'bo',label='u')
plt.plot(tts/86400.0,total_mags_2,'gx',label='g')
plt.plot(tts/86400.0,total_mags_3,'r*',label='r')
plt.xlabel(r"%s + time [days]"%time_of_observation)
#plt.ylabel(r"\Delta M")
plt.ylabel(r"Flux [photons/s]")
plt.ylim([0,5e13])
plt.legend()
plt.savefig(plotName)
plotName = "plots/delta_mag.eps"
plt.savefig(plotName)
plt.close()
