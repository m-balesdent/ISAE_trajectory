# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 11:24:54 2023

@author: mbalesde
"""

import numpy as np



def compute_atmos(altitude):
    """
    Function calculating the atmosphere parameters as functions of the altitude (US 76 model)
        
    :altitude: altitude in meters
    
    """
	
	# Definition of constant values
	#Altitude
    altitudeTable = np.zeros(8);
    altitudeTable[0] = 0.0;
    altitudeTable[1] = 11000.0; 
    altitudeTable[2] = 20000.0;
    altitudeTable[3] = 32000.0;
    altitudeTable[4] = 47000.0;
    altitudeTable[5] = 51000.0; 
    altitudeTable[6] = 71000.0;
    altitudeTable[7] = 84852.0;
    
    #Gradient of the temperature on the atmosphere section
    gradientTable = np.zeros(8);
    gradientTable[0] = -6.5e-3;
    gradientTable[1] =  0.0; 
    gradientTable[2] =  1.0e-3;
    gradientTable[3] =  2.8e-3;
    gradientTable[4] =  0.0;
    gradientTable[5] = -2.8e-3; 
    gradientTable[6] = -2.0e-3;
    gradientTable[7] =  0.0;
    temperatureSeaLevel= 288.15;
    pressureSeaLevel=101325.0;
    gravity0 = 9.80665;
    molecularWeightAir = 28.9644e-3;
    gazConstant = 8.31432;
    earthRadius = 6356766.0;
    ratioCpCv=1.40;    
    temperatureTable = np.zeros(8);
    pressureTable= np.zeros(8);
    
    temperatureTable[0] = temperatureSeaLevel;
    for i in range(1,8):
    	temperatureTable[i] = temperatureTable[i-1] + gradientTable[i-1]* (altitudeTable[i] - altitudeTable[i-1]);
    
    #Pressure at the beginning of the section
    pressureTable[0] = pressureSeaLevel;
    gmr = gravity0 * molecularWeightAir / gazConstant;
    for i in range(1,8):
    	pressureTable[i] = computePressure(temperatureTable[i-1],
    									temperatureTable[i],
    									gradientTable[i-1],
    									pressureTable[i-1], gmr,
    									altitudeTable[i-1],
    									altitudeTable[i]);
    
    

    
    imin=0;
    imax=7;
    
    gmr = gravity0 * molecularWeightAir / gazConstant;
    
    #Calculation of altitude
    h = altitude * earthRadius / (altitude + earthRadius);
    
    # 
    if (h < altitudeTable[0]):
    	temp = temperatureTable[0] + gradientTable[0] * (h-altitudeTable[0]);
    	press = computePressure(temperatureTable[0], temp, gradientTable[0],
    						pressureTable[0], gmr, altitudeTable[0], h);
    else:
    	if (h > altitudeTable[imax]):
    		temp = temperatureTable[imax];
    		press = computePressure(temperatureTable[imax], temp, 
    							gradientTable[imax], pressureTable[imax], 
    							gmr, altitudeTable[imax], h);
    	else:
    
    		
    		while((imax - imin) > 1):
    	
    			k = np.int((imin + imax) / 2);
    			if (h >= altitudeTable[k]):
    				imin = k;
    			else:
    				imax = k;
    
    		temp = temperatureTable[imin] + gradientTable[imin] * (h-altitudeTable[imin]);
    
    		press = computePressure(temperatureTable[imin], temp, 
    							gradientTable[imin], pressureTable[imin], 
    							gmr, altitudeTable[imin], h);
    

    dens = press * molecularWeightAir / (gazConstant * temp);
    sound = np.sqrt ((temp * ratioCpCv * gazConstant) / molecularWeightAir);
    
    return(temp,press,dens,sound)



    
def computePressure(temp0,temp,grad0,press0,gmr,h0,h):
	if (grad0 == 0.0):
		return press0 * np.exp(-gmr * (h-h0) / temp0)
	else:
		return press0 * pow(temp0/temp, gmr/grad0)
    