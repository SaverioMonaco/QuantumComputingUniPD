#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Usual imports
import matplotlib.pyplot as plt               # for plotting
import numpy as np                            # matrix handling

from scipy.optimize import curve_fit          # fitting a curve

import warnings
warnings.filterwarnings('ignore')             # Many warnings appear when trying to plt.savefig a figure that
                                              # has transparency


# In[2]:


s_1 = np.genfromtxt('s_1.csv')   # contains the normalized spacings from random hermitian matrices
s_2 = np.genfromtxt('s_2.csv')   # contains the normalized spacings from random real diagonal matrices


# In[3]:


plt.subplots(nrows=1, ncols=2, figsize=(15, 3))

plt.subplot(1, 2, 1)
# Plot the histogram of normalized spacing for random hermitian matrices
count_1, bin_1, _ = plt.hist(s_1[:,1:-1].flatten(),bins=100,color='forestgreen', density=True)
plt.title('Histogram of normalized spacing for random hermitian matrices')

plt.subplot(1, 2, 2)
# Plot the histogram of normalized spacing for random diagonal matrices
count_2, bin_2, _ = plt.hist(s_2[:,1:-1].flatten(),bins=100, density=True)
plt.title('Histogram of normalized spacing for random diagonal matrices')

plt.savefig('imgs/dists.svg', format='svg')

plt.show()


# In[4]:


# Plot the histogram in one figure
plt.hist(s_2[:,1:-1].flatten(), bins=100, alpha=0.5, label='Random Real-diagonal Matrices', density=True)
plt.hist(s_1[:,1:-1].flatten(), bins=100, alpha=0.5, label='Random Hermitian Matrices', density=True,color='forestgreen')
plt.legend(loc='upper right')
plt.xlim(-0.2,7)
plt.title("Comparison spacings generated")

plt.savefig('imgs/comparisondist.svg', format='svg')

plt.show()


# The function to fit to the histogram is:
# $$P(s)=a\,s^b\,\exp(-c\,s^d)$$
# $a,b,c,d$ to be determined

# In[5]:


# Function P(x) for fitting
def P(x, a, b, c, d):
    return a * np.power(x,b) * np.exp(-c * np.power(x, d) )

# Function log(P(x)) to avoid overflow in fitting (not used)
def Plog(x, a, b, c, d):
    return np.log(a) + b*np.log(x) - c * np.power(x,d)


# In[6]:


# Find the 4 parameters for each histogram
param_1 = curve_fit(P, bin_1[1:], count_1)
param_2 = curve_fit(P, bin_2[1:], count_2)


# In[7]:


#################################################################
# Plot histogram and fitted curve for random hermitian matrices #
#################################################################
plt.hist(s_1[:,1:-1].flatten(),bins=100,color='forestgreen',label='data', alpha=0.8,density=True)
plt.title('Normalized spacings of random hermitian matrices', fontsize=16)
plt.plot(bin_1[1:],P(bin_1[1:], param_1[0][0], param_1[0][1], param_1[0][2], param_1[0][3]),label='fit',color='red')

plt.text(4,.2, r'$P(s)=a\,s^b\,\exp(-c\,s^d)$' +
             '\n\na ='+str(round(param_1[0][0],4))+'\nb ='+str(round(param_1[0][1],4))+
               '\nc ='+str(round(param_1[0][2],4))+'\nd ='+str(round(param_1[0][3],4)),
         fontdict = {'fontsize' : 15})

plt.xlim(-.1,8)
plt.legend(prop={'size': 12})
plt.ylabel('Counts')

plt.savefig('imgs/fit_randomhermitian.svg', format='svg')

plt.show()


# In[8]:


################################################################
# Plot histogram and fitted curve for random diagonal matrices #
################################################################
plt.hist(s_2[:,1:-1].flatten(),bins=100,label='data', alpha=0.8,density=True,color='#1f77b4')
plt.title('Normalized spacings of random diagonal matrices', fontsize=16)
plt.plot(bin_2[1:],P(bin_2[1:], param_2[0][0], param_2[0][1], param_2[0][2], param_2[0][3]),label='fit',color='red')

plt.text(4,.2, r'$P(s)=a\,s^b\,\exp(-c\,s^d)$' +
             '\n\na ='+str(round(param_2[0][0],4))+'\nb ='+str(round(param_2[0][1],4))+
               '\nc ='+str(round(param_2[0][2],4))+'\nd ='+str(round(param_2[0][3],4)),
         fontdict = {'fontsize' : 15})

plt.xlim(-.1,8)
plt.legend(prop={'size': 12})
plt.ylabel('Counts')

plt.savefig('imgs/fit_randomdiagonal.svg', format='svg')

plt.show()


# In[ ]:




