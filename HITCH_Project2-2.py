#!/usr/bin/env python
# coding: utf-8

# # A.i #

# In[1]:


from astropy.io import fits as pyfits

import matplotlib.pyplot as plt

import numpy as np

import pandas as pd

import scipy.stats


# In[2]:


def pause():
    
    # pauses until a keyboard entry (e.g. carrage return)
    
    print('\r')     
    dummy = input('Pause')     
    print('\r')


# In[3]:


hdulist = pyfits.open('frame-g-000094-1-0131.fits') 

im = hdulist[0].data

hdulist.close()


# In[4]:


plt.imshow(im, vmin=np.percentile(im, 1), vmax=np.percentile(im, 99.9), cmap='gray', extent = [-51.030819, 759.977, -296.13902388, 294.2791194])

plt.title('Sloan Digital Sky Survey (SDSS) Image')

plt.xlabel('Right Ascension (arcseconds)')

plt.ylabel('Declination (arcseconds)')

pause()


# This grey-scale image looks very similar to the one in the Project 2 file, however it appears to have a little less detail, and is also upside down (due to how I plotted it with imshow)

# # A.ii #

# In[5]:


csvData = pd.read_csv('frame-g-000094-1-0131.csv')


# In[6]:


x = np.array(csvData[' Column'])

y = np.array(csvData['Row'])

gBand = np.array(csvData[' g-band'])


# In[7]:


empty1 = []                                          # define an empty list for appending

for i in range(len(x)):                              # loop through all the stars
    
    starPositions = np.array([x[i], y[i], gBand[i]]) # define array for star data of single star
    
    empty1.append(starPositions)                     # append data to empty list
    
twoDarray = np.array(empty1)                         # convert to a numpy array

unique2D = np.unique(twoDarray, axis=0)              # remove duplicate entries


# In[8]:


uniqueX = unique2D[:, 0]

uniqueY = unique2D[:, 1]

uniqueGband = unique2D[:, 2] 


# In[9]:


plt.scatter(uniqueX, uniqueY, marker='*')

plt.title('SDSS Scatter Plot')

plt.xlabel('Pixels')

plt.ylabel('Pixels')

pause()


# The largest difference that's apparent between my grey-scale image and this scatter plot is the number of visible stars. Many stars in the grey-scale image hardly show up and are so faint that you can't even see them. This scatter plot helps see all the stars in the image, but makes it difficult to tell you're even looking at the same data.

# # A.iii #

# In[10]:


gBandWeighted = np.array(2**(uniqueGband)/1000000)


# In[11]:


plt.scatter(uniqueX, uniqueY, marker='*', s = gBandWeighted)

plt.title('SDSS Weighted Scatter Plot')

plt.xlabel('Pixels')

plt.ylabel('Pixels')

pause()


# The differences between my scatter plot and my grey-scale image are mostly accounted for, but it's still kinda hard to tell they're the same image. 

# # B.i #

# In[12]:


np.random.seed(3)

randomX = np.random.uniform(0, 2048, len(uniqueX))


# In[13]:


np.random.seed(5)

randomY = np.random.uniform(0, 1498, len(uniqueY))


# In[14]:


plt.scatter(randomX, randomY, marker='*')

plt.title('Randomly Star Scatter Plot')

plt.xlabel('Pixels')

plt.ylabel('Pixels')

pause()


# Other than the physical locations of the stars, this scatter plot is qualitatively AND quantitatively indistinguishable from my scatter plot in part A.ii

# # B.ii #

# In[15]:


def neighbor(data):
    '''
    Takes a 2-D numpy array of x and y star positions, and returns the distance of the closest star, for each star
    
    Parameters:
        data : 2-D numpy array full of star positions
        
    Returns:
        Numpy array with length (# of stars) full of closest star distances for each star
    '''
    
    x, y = data.shape                                          # extract the shape of input array
    
    emptyList = []                                             # define an empty list for appending

    for i in range(x):                                         # loop through all the stars
    
        emptyList2 = []                                        # define another empty list for appending
    
        for j in range(x-1):                                   # loop through all the stars (excluding reference star)
        
            if j+1 != i:
        
                xDistance = (data[i, 0] - data[j+1, 0])**2     # find the change in x and square it
        
                yDistance = (data[i, 1] - data[j+1, 1])**2     # find the change in y and square it
        
                totalDistance = np.sqrt(xDistance + yDistance) # add the distances and take the square root
        
            emptyList2.append(totalDistance)                   # append distances to empty list
    
        nearest = np.min(np.array(emptyList2))                 # extract only the smallest distance
    
        emptyList.append(nearest)                              # append shortest neighboring star distances to empty list
    
    nearestNeighbors = np.array(emptyList)                     # convert to numpy array
    
    return nearestNeighbors


# In[16]:


empty2 = []                                              # use same technique as in part A.ii to create a 2D array of unique positions

for i in range(len(uniqueX)):
    
    uniquePositions = np.array([uniqueX[i], uniqueY[i]])
    
    empty2.append(uniquePositions)
    
uniquePositions2D = np.array(empty2)


# In[17]:


nearestStars = neighbor(uniquePositions2D)


# In[18]:


empty3 = []                                              # use same technique as in part A.ii to create a 2D array of random positions
for i in range(len(randomX)):
    
    randomPositions = np.array([randomX[i], randomY[i]])
    
    empty3.append(randomPositions)
    
randomPositions2D = np.array(empty3)


# In[19]:


nearestRandom = neighbor(randomPositions2D)


# # B.ii.a. #

# In[20]:


plt.hist(nearestRandom)

plt.title('Histogram A (random)')

plt.xlabel('Nearest Neighbor Distance (Pixels)')

plt.ylabel('# of Stars')

pause()


# In[25]:


plt.hist(nearestStars);

plt.title('Histogram B (SDSS Data)')

plt.xlabel('Nearest Neighbor Distance (Pixels)')

plt.ylabel('# of Stars')

pause()


# There is a large visual difference between these two histograms; The real star data from the image has many more stars that are closer together than my random star locations, which is evident by the huge count in bin 1 of the SDSS data histogram.

# # B.ii.b. #

# In[22]:


scipy.stats.skew(nearestRandom)


# In[23]:


scipy.stats.skew(nearestStars)


# In light of the visual differences between the two histograms, these values of skew make sense because the histogram for the star data has many more values on the left side, which results in a higher value of skewness. It's higher compared to my histogram for the random star locations, which is still skewed, but much less than the real star data. 

# # B.ii.c #

# In[24]:


scipy.stats.chisquare(nearestStars, nearestRandom)


# With a chi-squared value of nearly 20,000 and a p-value of 0, this implies that the nearest star distance for the real data is completely uncorrelated to the random distribution I created. 
