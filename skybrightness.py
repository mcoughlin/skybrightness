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

"""This module provides example methods to calculate the lunar contribution
to sky brightness using calculations from Coughlin, Stubbs, and Claver [2015].
"""

import datetime, time
import numpy as np
import matplotlib
#matplotlib.rc('text', usetex=True)
#matplotlib.use('Agg')
#matplotlib.rcParams.update({'font.size': 16})
import matplotlib.pyplot as plt

import ephem

# Inputs
latitude = '-30:15:06.37'
longitude = '-70:44:17.50'
elevation = 2552.0
passband = 'r'
ra_obs = 12.0 # Hours
dec_obs = -60.0 # Degrees
time_of_observation = "2015/05/05 6:00:00"

# Where is the moon?
moon = ephem.Moon()
utc_date = datetime.datetime.strptime(time_of_observation, "%Y/%m/%d %H:%M:%S")
moon.compute(utc_date)
ra_moon = (24/(2*np.pi))*float(repr(moon.ra))
dec_moon = (180/np.pi)*float(repr(moon.dec))

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

# Moon phase data (from Coughlin, Stubbs, and Claver Table 2) 
moon_phases = [2,10,45,90]
moon_data = {'u':[2.7,3.1,4.2,5.7],
             'g':[2.4,2.8,3.8,5.2],
             'r':[2.1,2.5,3.4,4.9],
             'i':[1.9,2.3,3.3,4.7],
             'z':[1.9,2.2,3.2,4.6],
             'y':[1.8,2.2,3.1,4.5]}

# Determine moon data for this phase
moon_data_passband = moon_data[passband]
delta_mag = np.interp(moon.moon_phase,moon_phases,moon_data_passband)
delta_mag_error = 0.1*delta_mag

# Fits to solar sky brightness (from Coughlin, Stubbs, and Claver Table 4) 
sun_data = {'u':[88.5,-0.5,-0.5,0.4],
            'g':[386.5,-2.2,-2.4,0.8],
            'r':[189.0,-1.4,-1.1,0.8],
            'i':[164.8,-1.5,-0.7,0.6],
            'z':[231.2,-2.8,-0.7,1.4],
            'zs':[131.1,-1.4,-0.5,0.2],
            'y':[92.0,-1.3,-0.2,0.9]}

sun_data_error = {'u':[6.2,0.1,0.1,0.1],
            'g':[34.0,0.2,0.2,0.5],
            'r':[32.7,0.2,0.2,0.5],
            'i':[33.1,0.2,0.2,0.5],
            'z':[62.3,0.3,0.4,0.9],
            'zs':[45.6,0.2,0.3,0.6],
            'y':[32.7,0.2,0.2,0.5]}

# Determine sun data for this phase
sun_data_passband = sun_data[passband]
flux = sun_data_passband[0] + sun_data_passband[1]*angle +\
       sun_data_passband[2]*alt_target + sun_data_passband[3]*alt_moon
flux = flux* (10**11)
flux_mag = -2.5 * np.log10(flux)

sun_data_passband_error = sun_data_error[passband]
flux_error = np.sqrt(sun_data_passband_error[0]**2 + sun_data_passband_error[1]**2 * angle**2 +\
       sun_data_passband_error[2]**2 * alt_target**2 + sun_data_passband_error[3]**2 * alt_moon**2)
flux_error = flux_error * (10**11)
flux_mag_error = 1.08574 * flux_error / flux

# Determine total magnitude contribution
total_mag = delta_mag + flux_mag
total_mag_error = np.sqrt(delta_mag_error**2 + flux_mag_error**2)

print "Sun-> Moon conversion: %.5f +- %.5f"%(delta_mag,delta_mag_error)
print "Sky brightness contribution: %.5f +- %.5f"%(flux_mag,flux_mag_error)
print "Total magnitude reduction: %.5f +- %.5f"%(total_mag,total_mag_error)


