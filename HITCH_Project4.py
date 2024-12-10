#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import scipy
import scipy.stats as stats
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits


# In[2]:


def pause():
    
    # pauses until a keyboard entry (e.g. carrage return)
    
    print('\r')     
    dummy = input('Pause')     
    print('\r')


# # A.i #

# In[3]:


hubbleData = np.genfromtxt('hubbleoriginal.csv',dtype=float,delimiter=',',usecols=(1,2),skip_header=1,unpack='True') # extract the Hubble data from the csv file


# In[4]:


rhub, vhub = hubbleData                                        # extract the distance and velocity columns from the Hubble data

plt.scatter(rhub, vhub)                                        # scatter plot velocities vs. distances

plt.title('Galaxy Recessional Velocity vs. Distance to Earth') # title plot and label axes
plt.xlabel('Distance (Mega Parsecs)')
plt.ylabel('Recessional Velocity (km/s)')

pause()


# # A.ii #

# In[35]:


values, variance = np.polyfit(rhub, vhub, 1, cov = True) # calculate slope/intercept/variance of the distances and velocities 
x = np.linspace(np.min(rhub), np.max(rhub), len(rhub))   # define a "time" array for easy plotting

m = values[0]                                            # extract the best slope value
b = values[1]                                            # extract the best intercept value

mVar = np.sqrt(variance[0, 0])                           # extract the variances in m and b
bVar = np.sqrt(variance[1, 1])


# In[44]:


plt.scatter(rhub, vhub)                                        # scatter plot velocities vs. distances
plt.plot(x, m*x+b, c = 'r')                                    # plot the line of best fit

plt.title('Galaxy Recessional Velocity vs. Distance to Earth') # title plot and label axes
plt.xlabel('Distance (Mega Parsecs)')
plt.ylabel('Recessional Velocity (km/s)')

pause()


# # A.iii #

# In[45]:


plt.scatter(rhub, vhub)                                                            # scatter plot velocities vs. distances
plt.plot(x, m*x+b, label = 'Best Fit', c = 'r')                                    # plot the line of best fit

plt.plot(x, (m + mVar)*x + b + bVar, ls = '--', c = 'g', label = '68% Confidence') # plot the uncertainty lines
plt.plot(x, (m - mVar)*x + b - bVar, ls = '--', c = 'g')

plt.title('Galaxy Recessional Velocity vs. Distance to Earth')                     # title plot and label axes
plt.xlabel('Distance (Mega Parsecs)')
plt.ylabel('Recessional Velocity (km/s)')

plt.legend()

pause()


# # A.iv #

# In[8]:


scipy.stats.pearsonr(rhub, vhub) # calculate the Pearson correlation coefficient and the p-value of the original Hubble data


# # B.i #

# In[9]:


slowGalaxyData = np.genfromtxt('ned1dlevel5.csv',dtype=float,delimiter=',',usecols=(3, 11),skip_header=2,unpack='True') # extract galaxy data for v < c/8 from csv file
rhubSlow, vhubSlow = slowGalaxyData                                                                                                  # extract the distance and velocity columns from the slow galaxy data                                 


# In[10]:


fastGalaxyData = np.genfromtxt('ned4dlevel5.csv',dtype=float,delimiter=',',usecols=(3, 10),skip_header=2,unpack='True') # extract galaxy data for v > c/8 from csv file
rhubFast, vhubFast = fastGalaxyData                                                                                                  # extract the distance and velocity columns from the fast galaxy data                    


# In[48]:


plt.scatter(rhubFast, vhubFast)                       # scatter plot the fast galaxy data
plt.scatter(rhubSlow, vhubSlow, c = 'g')                       # scatter plot the slow galaxy data
plt.scatter(rhub, vhub, c = 'r')                               # scatter plot the Hubble data



plt.title('Galaxy Recessional Velocity vs. Distance to Earth') # title plot and label axes
plt.xlabel('Distance (Mega Parsecs)')
plt.ylabel('Recessional Velocity (km/s)')

pause()


# # B.ii #

# In[12]:


c = scipy.constants.c                                                 # define the speed of light in m/s

velocities = np.append(np.append(vhub, vhubSlow), vhubFast) * 1000    # append all the velocity data to one array and multiply by 1000 to convert to m/s
distances = np.append(np.append(rhub, rhubSlow), rhubFast)            # append all the distance data to one array for later

z = np.sqrt(1 + (velocities / c)) / np.sqrt(1 - (velocities / c)) - 1 # calculate z using the formula from the Project 4 file

vEffective = c * z                                                    # mulitply z by the speed of light for a normalized velocity array

vhubEff = vEffective[0:len(vhub)]                                     # splice the effective orignal Hubble velocites

vSlowEff = vEffective[len(vhub):len(vhub) + len(vhubSlow)]            # splice the effective slow galaxy velocites

vFastEff = vEffective[len(vhub) + len(vhubSlow):]                     # splice the effective slow galaxy velocites


# In[38]:


plt.scatter(distances, vEffective)                                                 # scatter plot all the effective velocites vs. distances

plt.title('Effective Galaxy Recessional Velocity vs. Distance to Earth', y = 1.05) # title plot and label axes
plt.xlabel('Distance (Mega Parsecs)')
plt.ylabel('Effective Recessional Velocity (m/s)')

pause()


# # B.iii #

# In[14]:


hubEffValues, hubEffVariance = np.polyfit(rhub, vhubEff, 1, cov = True) # calculate slope/intercept/variance of the distances and effective velocities 

mEff = hubEffValues[0]                                                  # extract the effective slope
bEff = hubEffValues[1]                                                  # extract the effective intercept


# In[15]:


slowEffValues, slowEffVariance = np.polyfit(rhubSlow, vSlowEff, 1, cov = True) # calculate slope/intercept/variance of the distances and effective slow velocities 

mSlowEff = slowEffValues[0]                                                    # extract the effective slope
bSlowEff = slowEffValues[1]                                                    # extract the effective intercept

xSlow = np.linspace(np.min(rhubSlow), np.max(rhubSlow), len(rhubSlow))         # define a "time" array for easy plotting


# In[16]:


fastEffValues, fastEffVariance = np.polyfit(rhubFast, vFastEff, 1, cov = True) # calculate slope/intercept/variance of the distances and effective fast velocities

mFastEff = fastEffValues[0]                                                    # extract the effective slope
bFastEff = fastEffValues[1]                                                    # extract the effective intercept

xFast = np.linspace(np.min(rhubFast), np.max(rhubFast), len(rhubFast))         # define a "time" array for easy plotting


# In[37]:


plt.scatter(distances, vEffective)                                                 # scatter plot effective velocities vs. distance

plt.plot(x, mEff*x + bEff, c = 'r', lw = 3, label = 'Hubble Best Fit')             # plot the line of best fit for Hubble galaxy data

plt.plot(xSlow, mSlowEff*xSlow + bSlowEff, c = 'g', label = 'Slow Best Fit')       # plot the line of best fit for slow galaxy data

plt.plot(xFast, mFastEff*xFast + bFastEff, c = 'orange', label = 'Fast Best Fit')  # plot the line of best fit for fast galaxy data

plt.title('Effective Galaxy Recessional Velocity vs. Distance to Earth', y = 1.05) # title plot and label axes
plt.xlabel('Distance (Mega Parsecs)')
plt.ylabel('Recessional Velocity (m/s)')

plt.legend()

pause()


# # B.iv #

# In[18]:


scipy.stats.pearsonr(rhubSlow, vhubSlow) # calculate the Pearson correlation coefficient and the p-value of the slow galaxy data


# In[19]:


scipy.stats.pearsonr(rhubFast, vhubFast) # calculate the Pearson correlation coefficient and the p-value of the fast galaxy data


# # B.vi #

# In[20]:


slowMethods = np.genfromtxt('ned1dlevel5.csv',dtype=str,delimiter=',',usecols=(4),skip_header=2,unpack='True') # extract method data from the slow galaxy csv file

slowSupernovas = slowMethods == 'SNIa'                                                                                      # find elements that are Supernova Type Ia

vSlowSN = slowSupernovas * vSlowEff                                                                                         # use only SNIa effective velocities

rSlowSN = slowSupernovas * rhubSlow                                                                                         # use only SNIa distances


# In[21]:


fastMethods = np.genfromtxt('ned4dlevel5.csv',dtype=str,delimiter=',',usecols=(4),skip_header=2,unpack='True') # extract method data from the fast galaxy csv file

fastSupernovas = fastMethods == 'SNIa'                                                                                      # find elements that are Supernova Type Ia

vFastSN = fastSupernovas * vFastEff                                                                                         # use only SNIa effective velocities

rFastSN = fastSupernovas * rhubFast                                                                                         # use only SNIa distances


# In[22]:


vSupernovas = np.append(vSlowSN, vFastSN) # append the Supernova velocites to 1 array

rSupernovas = np.append(rSlowSN, rFastSN) # append the Supernova distances to 1 array


# In[42]:


mSN, bSN = np.polyfit(rSupernovas, vSupernovas/c, 1)                          # calculate slope/intercept of the Supernova distances and velocities 

xSN = np.linspace(np.min(rSupernovas), np.max(rSupernovas), len(rSupernovas)) # define a "time" array for easy plotting


# In[47]:


plt.scatter(rSupernovas, vSupernovas/c) # scatter plot Redshift vs. Distance for SNIa's

plt.plot(xSN, mSN*xSN + bSN, c = 'r')   # plot the best-fit line

plt.title('Redshift vs. Distance')      # title plot and label axes
plt.xlabel('Distance (Mega Parsecs)')
plt.ylabel('Redshift')

pause()


# # C.i #

# In[24]:


raHours = np.genfromtxt('hubbleoriginal.csv',dtype=float,delimiter=',',usecols=(3),skip_header=1,unpack='True')    # extract Right Ascension data
decDegrees = np.genfromtxt('hubbleoriginal.csv',dtype=float,delimiter=',',usecols=(4),skip_header=1,unpack='True') # extract the Declination data


# # C.ii #

# In[25]:


ra = raHours * np.pi/12      # convert RA from hours to radians
dec = decDegrees * np.pi/180 # convert DEC from degrees to radians


# # C.iii #

# In[26]:


def fit_correction(x,y,r,d):
    """
    This function will fit Hubble's correction to his original dataset.
    This correction takes the form:
    v = Ho*r + Xcos(a)cos(d) + Ysin(a)cos(d) + Zsin(d)
    where v is the recessional velocity of a galaxy, r is the distance
    to the galaxy, a is the right ascension and d is the declination.
    Input:
    x = The distances to a set of galaxies
    y = The corresponding recessional velocities
    r = The corresponding right ascensions
    d = The corresponding declinations
    Output:
    fit = The linear fit to the corrected velocities and distances
    correction = The correction to the original velocities provided
    ---------------------------------------------------------------
    Example:
    V, D, ra, dec = np.genfromtxt(datafile).T
    fit, dV = fit_correction(D,V,ra,dec)
    """
    # Generating the X matrix for the correction model above
    X = np.vstack(( x,
                    np.cos(r)*np.cos(d),
                    np.sin(r)*np.cos(d),
                    np.sin(d) ))
    
    # Fitting and reporting the coefficients of the fit
    A = np.dot( np.linalg.inv(np.dot(X,X.T)), np.dot(y,X.T) )
    print("Ho: {}, X: {}, Y: {}, Z: {}".format(*tuple(A)))
    
    # Finding the fit and correction
    fit = x * A[0]
    correction = np.dot(A[1:],X[1:])
    
    # Finding and reporting the chi square and p values
    x2 = ( (y-fit)**2 / np.var(y-fit) ).sum()
    p = 1 - stats.chi2.cdf(x2, x.size - 2)
    print("X2: {}, p-val: {}".format(x2,p))
    return fit, correction


# In[27]:


vhubEffKMS = vhubEff / 1000                                       # convert the effective velocites back to km/s

hubFit, hubCorrection = fit_correction(rhub, vhubEffKMS, ra, dec) # use the fit_correction function on the effective data


# # C.iv #

# In[49]:


correctedVelocities = 465.914209474116*rhub - 68.77896102613333*np.cos(ra)*np.cos(dec) + 236.47980365532118*np.sin(ra)*np.cos(dec) - 200.27641364430508*np.sin(dec) # manually define the corrected velocities

plt.scatter(rhub, correctedVelocities)                                                                                                                              # scatter plot the corrected velocities vs. distance

plt.plot(rhub, hubFit, c = 'r')                                                                                                                                     # plot the line of best fit from fit_correction

plt.title('Corrected Recessional Velocities vs. Distance')                                                                                                          # title plot and label axes
plt.xlabel('Distance (Mega Parsecs)')
plt.ylabel('Corrected Recessional Velocities (km/s)')

pause()


# # C.v #

# In[29]:


scipy.stats.pearsonr(rhub, correctedVelocities) # calculate the Pearson correlation coefficient and the p-value of the corrected velocity data


# In[30]:


np.polyfit(rhub, correctedVelocities, 1, cov = True) # calculate slope/intercept/variance of the distances and corrected velocities 

