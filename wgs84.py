# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 11:41:14 2013

@author: Rafal Samborski (rafal@samborski.co)

Based on java class written by Christopher Dorner (http://www.avionixsl.com/)
Formulars adapted from Book: Server-side Gps and Assisted-Gps in Java (by Neil Harper)
Can be found on Google Books, Chapter 2

"""

import numpy as np
    
A_SEMIMAJOR_AXIS = 6378137.0;

F_FLATTENING = 1.0 / 298.257223563;

B_SEMIMINOR_AXIS = 6356752.3142;

# second eccentricity squared
ESQ_DASH = (A_SEMIMAJOR_AXIS**2 - B_SEMIMINOR_AXIS**2) / B_SEMIMINOR_AXIS**2;

# eccentricity squared
ESQ = 2 * F_FLATTENING - F_FLATTENING**2;    
 
 
def calculateNlat( lat ):
    return A_SEMIMAJOR_AXIS / np.sqrt( 1 - ( ESQ * np.sin( np.radians( lat ) )**2 ) );


def lla2ecef( lat, lon, alt ):
    
    #converts ellipsoidal coordinates into cartesian ecef
    #lat in decimal degrees
    #lon in decimal degrees
    #alt in meters
    
    # check for bad values
    if lat < -90 or lat > 90 or lon < -180 or lon > 180:
        return None

    N = calculateNlat( lat );
    x = (N + alt) * np.cos( np.radians( lat ) ) * np.cos( np.radians( lon ) );
    y = (N + alt) * np.cos( np.radians( lat ) ) * np.sin( np.radians( lon ) );
    z = ((N * (1 - ESQ)) + alt) * np.sin( np.radians( (lat) ) );
    return x, y, z
    
def ecef2lla( x, y, z ):
    #converts ecef to lla
    #x in meters
    #y in meters
    #z in meters

    lat = 0.0;
    lon = 0.0;
    alt = 0.0;

    # if it is one of the poles, no need to compute
    if x == 0 and y == 0:
        lon = 0.0;
        if z > 0:
            lat = 90;
            alt = z - B_SEMIMINOR_AXIS;
        else:
            lat = -90;
            alt = z + B_SEMIMINOR_AXIS;
        return lat, lon, alt

    p = np.sqrt( x**2 + y**2 );
    
    # parametric latitude of point p
    theta = np.arctan2( (z * A_SEMIMAJOR_AXIS), (p * B_SEMIMINOR_AXIS) );
    lat = np.rad2deg( np.arctan2( (z + (ESQ_DASH * B_SEMIMINOR_AXIS * \
        np.power( np.sin( theta ), 3 ))), (p - (ESQ * A_SEMIMAJOR_AXIS * \
        np.power( np.cos( theta ), 3 ))) ) );
    lon = np.rad2deg( np.arctan2( y, x ) );
    alt = (p / np.cos( np.radians( lat ) )) - calculateNlat( lat );
    
    return lat, lon, alt;