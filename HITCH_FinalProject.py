#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from matplotlib.pyplot import figure


# In[2]:


def pause():
    
    # pauses until a keyboard entry (e.g. carrage return)
    
    print('\r')     
    dummy = input('Pause')     
    print('\r')


# In[3]:


data = pyfits.getdata('spec-0429-51820-0056.fits.gz') # read in the spectrum data from the project files
flux = np.array([i[0] for i in data]) # extract the flux data from the spectrum data file
wavelength = np.array([10**i[1] for i in data]) # extract the wavelength data from the spectrum data file


# In[4]:


data


# In[5]:


plt.plot(wavelength, flux) # plot flux vs. wavelength

plt.title('Flux vs. Wavelength')  # title plot and label axes
plt.ylabel('Flux')
plt.xlabel('Wavelength')

pause()


# In[6]:


fit = np.polyfit(wavelength, flux, 4) # define a polyfit to correct the spectrum data


# In[7]:


correctedFlux = (flux - (fit[0]*wavelength**4 + fit[1]*wavelength**3 + fit[2]*wavelength**2 + fit[3]*wavelength + fit[4])) # corrected flux data array


# In[8]:


plt.plot(wavelength, correctedFlux) # plot corrected flux vs. wavelength

plt.title('Corrected Flux vs. Wavelength')  # title plot and label axes
plt.ylabel('Corrected Flux')
plt.xlabel('Wavelength')

pause()


# In[9]:


def findShift(fluxes, reference):
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
        
        crossCor = np.fft.ifft(np.fft.fft(i)*np.conj(np.fft.fft(reference)))
        
        shift = np.where(crossCor == np.max(crossCor))
        
        shifts.append(shift)
        
    return shifts


# In[10]:


lineData = np.genfromtxt('linelist-0429-51820-0056.csv',dtype=float,delimiter=',',usecols=(1),skip_header=1,unpack='True') # read in the rest wavelength data from the project files


# In[11]:


lineData # show the data


# In[12]:


restarray = np.zeros(len(wavelength)) - 100 # define a rest wavelength array
points = [] # empty list for appending
for i in range(len(lineData)): # for loop to determine where the wavelength data is greater than or equal to rest wavelengths
    wh = np.where(wavelength >= lineData[i])
    points.append(np.min(wh)) # append the minimum wavelength element to the points array
for i in range(len(points)): # for loop to create vertical lines at each rest wavelength point
    restarray[points[i]] = 150
restarray[0] = -100


# In[13]:


figure(figsize = (20, 5)) # define figure
plt.plot(wavelength, correctedFlux, label='Flux') # plot corrected flux vs. wavelength
plt.plot(wavelength, restarray, label='Rest Wavelengths') # plot rest array vs. wavelength

plt.title('Corrected Flux vs. Wavelength')  # title plot and label axes
plt.ylabel('Corrected Flux')
plt.xlabel('Wavelength')

plt.legend()

pause()


# In[14]:


figure(figsize = (20, 5))
plt.plot(wavelength[2400:2600], correctedFlux[2400:2600], label='Flux') # plot corrected flux vs. wavelength on specific interval
plt.plot(wavelength[2400:2600], restarray[2400:2600], label='Rest Wavelengths') # plot rest array vs. wavelength on specific interval

plt.axvline(wavelength[2426], c='g', label='Intervals') # section off each local maxima
plt.axvline(wavelength[2436], c='g')
plt.axvline(wavelength[2448], c='g')
plt.axvline(wavelength[2460], c='g')

plt.axvline(wavelength[2537], c='g')
plt.axvline(wavelength[2547], c='g')
plt.axvline(wavelength[2555], c='g')

plt.title('Corrected Flux vs. Wavelength')  # title plot and label axes
plt.ylabel('Corrected Flux')
plt.xlabel('Wavelength')

plt.legend()

pause()


# In[15]:


redShifts = [] # define empty list for appending

redShifts.append(wavelength[np.where(correctedFlux == np.max(correctedFlux[2426:2436]))[0][0]]/lineData[-6]) # find the maximum value's location in each interval, and divide by the rest wavelength to calculate the redshift

redShifts.append(wavelength[np.where(correctedFlux == np.max(correctedFlux[2436:2448]))[0][0]]/lineData[-5]) # repeat for 4 other lines

redShifts.append(wavelength[np.where(correctedFlux == np.max(correctedFlux[2448:2460]))[0][0]]/lineData[-4])

redShifts.append(wavelength[np.where(correctedFlux == np.max(correctedFlux[2537:2547]))[0][0]]/lineData[-3])

redShifts.append(wavelength[np.where(correctedFlux == np.max(correctedFlux[2547:2555]))[0][0]]/lineData[-2])


# In[16]:


np.array(redShifts) - 1


# In[17]:


avgRedshift = np.average(np.array(redShifts) - 1) # average redshift 
avgRedshift


# In[18]:


error = [] # define empty list for appending
for i in range(5): # for loop to run through all my redShift data
    error = np.append(error, np.std(redShifts[:i+1])) # calculate the standard deviation for different number of spectrum lines


# In[19]:


error


# In[20]:


print(f'So the redshift is {avgRedshift*100:.2f}% +/- {error[4]*100:.2f}%')

