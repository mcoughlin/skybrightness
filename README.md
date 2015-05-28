# skybrightness
Calculate lunar sky brightness based on Coughlin, Stubbs, and Claver [2015]

HOW-TO:

Requirements:
datetime, numpy, pyephem

Inputs:
skybrightness.py takes 7 main arguments:
Latitude, Longitude, and Elevation of site
Passband of your instrument
Right ascension (in hours) and Declination (in degrees) of target object
Time of observation (%Y/%m/%d %H:%M:%S string) 

Outputs:
Angles between the moon and the target
Altitude and azimuth of the target
Altitude and azimuth of the moon
Magnitude contribution from sun->moon conversion
Magnitude contribution from moon

What it does:
It computes the angles described above and then uses the data given in Tables 2 and 4 to contribute the contributions.

