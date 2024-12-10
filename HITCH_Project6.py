#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import scipy
import matplotlib


# In[2]:


def pause():
    
    # pauses until a keyboard entry (e.g. carrage return)
    
    print('\r')     
    dummy = input('Pause')     
    print('\r')


# # A.i #

# In[3]:


data = np.genfromtxt('wasp_4b.tsv', skip_header = 36, usecols = (0, 1, 2), delimiter = ';') # import the data from tsv file


# In[4]:


time = data[:, 0] # extract the Julian time information from the data

flux = data[:, 1] # extract the relative flux information from the data


# # A.ii #

# In[5]:


plt.plot(time, flux) # plot flux vs. time 

plt.title('Relative Flux vs. Time') # title plot and label axes
plt.xlabel('Barycentric Julian Date')
plt.ylabel('Relative Flux')

pause()


# # A.iii #

# In[6]:


plt.plot(np.diff(time)) # plot the differences between neighboring values of the time array to find the largest jumps in time
plt.ylim(0, 10)

pause()


# In[41]:


np.where(np.diff(time)>1)[0] # find the element locations of the different starting times for the 6 different transits


# In[34]:


#for i range(len(transits) - 1):
    #transits[i] = np.array([time[np.where(np.diff(time)>1)[0][i]:np.where(np.diff(time)>1)[0][i + 1]], flux[np.where(np.diff(time)>1)[0][i]:np.where(np.diff(time)>1)[0][i + 1]]])
    
# I tried to make for loop but I ended up having to hard code this part:

transit1 = np.array([time[0:299], flux[0:299]])

transit2 = np.array([time[300:625], flux[300:625]])

transit3 = np.array([time[626:994], flux[626:994]])

transit4 = np.array([time[995:1400], flux[995:1400]])

transit5 = np.array([time[1401:1765], flux[1401:1765]])

transit6 = np.array([time[1766:len(time)], flux[1766:len(time)]])


# # A.iv #

# In[37]:


#for i in range(len(transits)):
    #transits[i] = np.array((transits[i][0, :] - np.min(transits[i][0, :])), transits[i][1, :])
    
# I tried to make for loop again but also ended up having to hard code this part:

transit1 = np.array([transit1[0, :] - np.min(transit1[0, :]), transit1[1, :]]) # subtract the smallest time value from each time array, making them all start at t = 0

transit2 = np.array([transit2[0, :] - np.min(transit2[0, :]), transit2[1, :]])

transit3 = np.array([transit3[0, :] - np.min(transit3[0, :]), transit3[1, :]])

transit4 = np.array([transit4[0, :] - np.min(transit4[0, :]), transit4[1, :]])

transit5 = np.array([transit5[0, :] - np.min(transit5[0, :]), transit5[1, :]])

transit6 = np.array([transit6[0, :] - np.min(transit6[0, :]), transit6[1, :]])


# In[38]:


transits = [transit1, transit2, transit3, transit4, transit5, transit6] # define a list of the six transits

for i in transits:
    plt.plot(i[0, :], i[1, :])
    
plt.title('Relative Flux vs. Time') # title plot and label axes
plt.xlabel('Relative Time')
plt.ylabel('Relative Flux')

pause()


# In[11]:


colors = ['b', 'g', 'r', 'orange', 'purple', 'pink'] # define a list of colors

for i in range(len(transits)): # for loop to scatter plot all the transits starting at t = 0
    plt.scatter(transits[i][0, :] - np.min(transits[i][0, :]), transits[i][1, :], c = colors[i])
    
plt.title('Relative Flux vs. Time') # title plot and label axes
plt.xlabel('Relative Time')
plt.ylabel('Relative Flux')

pause()


# # B.i.a #

# In[12]:


newTime = np.arange(0, np.max(transit4[0, :]), np.max(np.diff(transit4[0, :]))) # uniform time array using the longest transit (transit 4)

dt = np.max(np.diff(transit4[0, :])) # define a uniform change in time


# In[28]:


newFluxes = [newFlux1, newFlux2, newFlux3, newFlux4, newFlux5, newFlux6] # define a list of interpolated fluxes

for i in range(len(transits)): # for loop that interpolates each flux data with the new time array
    newFluxes[i] = np.interp(newTime, transits[i][0, :], transits[i][1, :])


# In[26]:


for i in newFluxes:
    plt.plot(newTime, i)

plt.title('Relative Flux vs. Time') # title plot and label axes
plt.xlabel('Relative Time')
plt.ylabel('Relative Flux')

pause()


# In[15]:


physShifts = []

for i in newFluxes: # for loop that interpolates each flux time-series to the reference, shifting each with respect to the reference one time-step at a time, multiplying, and summing
    
    crossCor = np.zeros(len(newFlux2)) # empty array of length 35

    for shift in range(len(newFlux2)): # for loop to calculate a singular maximum cross correlation, from the Coding Tips
        
        crossCor[shift] = np.sum(newFlux2 * np.roll(i, shift))
        
        relShift = 35 - np.where(crossCor == np.max(crossCor))[0][0]
        
    physShifts = np.append(physShifts, relShift)
        
physShifts


# # B.i.b #

# In[16]:


def findShift(fluxes):
    '''
    Takes: 
        Array or list of different flux arrays
        
    Parameters:
        fluxes : array or list full of multiple transit fluxes
        
    Returns:
        list of time-shifts to maximize the cross correlation between the reference flux data and the other transit fluxes
    '''
    
    shifts = []
    
    for i in fluxes:
        
        crossCor = np.fft.ifft(np.fft.fft(i)*np.conj(np.fft.fft(newFlux2)))
        
        shift = np.where(crossCor == np.max(crossCor))
        
        shifts.append(shift)
        
    return shifts


# In[17]:


specShifts = findShift(newFluxes)
specShifts


# # B.ii #

# In[18]:


for i in range(len(transits)): # for loop that shifts each transit by my calculated shifts in part B.i
    plt.scatter(transits[i][0, :] - dt*(specShifts[i][0]+1), transits[i][1, :], c = colors[i])
    
plt.title('Relative Flux vs. Time') # title plot and label axes
plt.xlabel('Relative Time')
plt.ylabel('Relative Flux')

pause()


# # C.i #

# In[19]:


correctedTimes = [] # empty flux and time lists for appending
correctedFluxes = []

for i in range(len(transits)): # for loop that appends all the corrected time and flux data to their own lists
    correctedTimes = np.append(correctedTimes, transits[i][0, :] - dt*(specShifts[i][0]+1))
    correctedFluxes = np.append(correctedFluxes, transits[i][1, :])


# In[20]:


correctedTransits = np.array([correctedTimes, correctedFluxes])


# In[21]:


plt.scatter(correctedTransits[0, :], correctedTransits[1, :]) # scatter plot the corrected data
plt.axvline(0.055) # estimate the x-value of the lowest point in the transits
plt.axhline(0.973) # estimate the y-value of the lowest point in the transits

plt.title('Relative Flux vs. Time') # title plot and label axes
plt.xlabel('Relative Time')
plt.ylabel('Relative Flux')

pause()


# # C.ii #

# In[22]:


plt.scatter(correctedTransits[0, :], correctedTransits[1, :]) # scatter plot the corrected data
plt.axhline(0.975) # estimate upper error bound for transit
plt.axhline(0.971) # estimate lower error bound for transit

plt.axhline(1.002) # estimate upper error bound for non-transit
plt.axhline(0.998) # estimate lower error bound for non-transit

transUncertainty = 0.975 - 0.971 # spread of values
nonTransUncertainty = 1.002 - 0.998 # spread of values
print(transUncertainty, nonTransUncertainty)

plt.title('Relative Flux vs. Time') # title plot and label axes
plt.xlabel('Relative Time')
plt.ylabel('Relative Flux')

pause()


# # C.iii #

# In[23]:


# the uncertainty in the transit depth is +/- 0.002 and the uncertainty in the radius of WASP-4 is +/- 0.04, so:
rPlanet = np.sqrt((1-0.973)*0.87**2) # radius of the planet in units of Solar Radii

rPlanetError = rPlanet * np.sqrt((0.002 / 1)**2 + (0.04 / 0.87)**2) # using error propagation methods to find the uncertainty in the radius of the planet

print(f'So the radius of the planet WASP-4b is {rPlanet:.3f} +/- {rPlanetError:.3f} Rsun')


# # C.iv #

# In[24]:


# using rPlanet = 0.143 +/- 0.007 Rsun, and mPlanet = 1.21 +/- 0.12 mJupiter
rSun = 7e8 # radius of the Sun in meters
mJupiter = 1.898e27 # mass of Jupiter in kilograms

rPlanet_meters = rPlanet * rSun # radius of WASP-4b in meters
rPlanetError_meters = rPlanetError * rSun # uncertainty in radius of WASP-4b in meters

mPlanet = 1.21 * mJupiter # mass of WASP-4b in kilograms
mPlanetError = 0.12 * mJupiter # uncertainty in mass of WASP-4b in kilograms

vPlanet = 4/3 * np.pi * rPlanet_meters**3 # volume of WASP-4b in cubic meters
vPlanetError = vPlanet * np.sqrt(3*(rPlanetError_meters/rPlanet_meters)**2)

density_Planet = mPlanet / vPlanet # density of WASP-4b in kilograms per cubic meter
density_PlanetError = density_Planet * np.sqrt((mPlanetError/mPlanet)**2 + (vPlanetError/vPlanet)**2) # uncertainty in density of WASP-4b in kilograms per cubic meter

print(f'So the density of the planet WASP-4b is {density_Planet:.1f} +/- {density_PlanetError:.1f} kg/m^3')

