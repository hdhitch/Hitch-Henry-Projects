#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def pause():
    
    # pauses until a keyboard entry (e.g. carrage return)
    
    print('\r')     
    dummy = input('Pause')     
    print('\r')


# In[ ]:


def nuemann(x, y):
    '''
    Takes:
        A 6-digit seed number [x], and returns a list length [y + 1] of random numbers
    
    Parameters:
        x : 6-digit seed
        y : desired length of list - 1
        
    Returns:
        List length [y + 1] full of psuedo-random numbers using Jon Von Nuemanns original technique
    '''
    
    first = str(x**2)             # First seed gets squared and converted into a string for splicing
    
    list1 = [x]                   # Define a list that starts with our input seed
    
    for i in range(y):
        newSeed = int(first[3:9]) # First string is spliced to take the middle 6 digits and turn them back into an integer
        
        first = str(newSeed**2)   # Redefine the string as the new seed squared
        
        list1.append(newSeed)     # Append the new seed to list1
            
    return list1                  # Return list 1


# In[3]:


def nuemannUpdated(x, y):
    '''
    Takes: 
        A 6-digit seed number [x], and returns a list length [y + 1] of pseudo-random numbers
    
    Parameters:
        x : 6-digit seed
        y : desired length of list
        
    Returns:
        List length [y + 1] full of psuedo-random numbers using an updated Nuemann technique which doesnt generate infinite loops or 0s as often
    '''
    
    
    first = str(x**2)                           # First seed gets squared and converted into a string for splicing
    
    list1 = [x]                                 # Define a list that starts with our input seed
    
    for i in range(y):
        
        newSeed = int(first[3:9])
        
        if newSeed%100 == 0 or newSeed < 100000: # If the new seed is evenly divisble by 100 or less than 6 digits, does something different with the number
            
            first = str((newSeed+123456)**2)     # Takes the current seed and adds 123456 to it
            
            list1.append(newSeed)                # Appends back to the list
            
        else:                                    # Otherwise, applies Jon Von Nuemann's original techinque 
            first = str(newSeed**2)
            list1.append(newSeed)
            
    return list1


# In[4]:


def chiHist(observed, bins):
    '''
    Takes:
        An observed data set, an expected data set, and a number of bins
        
    Parameters:
        observed : 1-D numpy.array of data
        expected : 1-D numpy.array of expected data, len(observed)
        bins : histogram-like bins that splits up the data, ideally a minimum of 5 numbers per bin
        
    Returns:
        The chi squared value between the obvserved and expected data sets
        The correlated p-value
    '''
    plt.hist(observed, bins)                                     # Create a histogram with the observed data
    
    axis = plt.gca()                                             # Assign an axis so we can extract the bin heights of the histogram

    patches = axis.patches

    heights = np.array([patch.get_height() for patch in patches])
        
    expected = len(observed)/bins                                 # Expected = Probability * Number of events
    
    chi = np.sum(((heights-expected)**2)/expected)                # Define chi squared according to the equation, and sum up over all the data
    
    dof = bins - 1                                                # Only have 1 constraint, so the # of degrees of freedom is equal to the # of bins - 1
    
    p = 1 - stats.chi2.cdf(chi, dof)                              # 1 - CDF = P-value
    
    return chi, p                                                 # Return chi squared and p-value


# In[5]:


def chiSquared(observed, expected):
    chi = np.sum(((observed-expected)**2)/expected) 
    
    return chi


# In[6]:


def serialTest(observed, bins):
    same = 0
    neighboring = 0
    
    for i in range(len(observed)-1):                                                            # For loop that adds 1 if two values in the observed data are in the same bin
        for j in range(bins):
            if j/bins <= observed[i] < (j+1)/bins and j/bins <= observed[i+1] < (j+1)/bins:
                same += 1
                
    for i in range(len(observed)-1):                                                            # For loop that adds 1 if two values in the observed data are in neighboring bins
        for j in range(bins):
            if j/bins <= observed[i] < (j+1)/bins and (j-1)/bins <= observed[i+1] < (j+2)/bins: 
                neighboring += 1
                
    neighboring -= same                                                                         # Account for the neighboring loop counting values in the same bin too
    
    return same, neighboring

