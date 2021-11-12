#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np
import numpy.polynomial.polynomial as poly


# In[2]:


# Load the data points from files
loop     = np.genfromtxt('loop.csv', delimiter=',')
loop2     = np.genfromtxt('loop2.csv', delimiter=',')
internal = np.genfromtxt('internal.csv', delimiter=',')


# In[3]:


# Generic dense x vector for plotting lines
x = np.arange(loop[0,1],loop[-1,1],1)


# In[4]:


###########################################
#        LOOP MATRIX MULTIPLICATION       #
###########################################
plt.figure()
# Plot the scatter of timings
plt.plot(loop[:,1],loop[:,0], 'o', ms=3, label='Data points')

# Fit a polynomial
coefs = poly.polyfit(loop[:,1],loop[:,0], 3)
fit_loop = poly.polyval(x, coefs)

# This is for printing the equation
eqstr = "f(x) = "
for i in range(len(coefs)):
    if i == (len(coefs) - 1):
        eqstr += str("{:.2e}".format(coefs[i]))
    else:
        eqstr += "+"+str("{:.2e}".format(coefs[i]))+"*x^"+str(len(coefs)-1-i)
    if i == 1:
        eqstr +="+\n          "
plt.text(0, 4, eqstr, fontdict = {'fontsize' : 12})

plt.plot(x,fit_loop,label='Fit')
plt.title("Time Scaling (Loop)",fontsize=16)
plt.legend(prop={'size': 12})
plt.ylabel('t[s]     ').set_rotation(0)
plt.xlabel('N')
plt.savefig('imgs/loop.svg', format='svg')


# In[5]:


###########################################
#        LOOP2 MATRIX MULTIPLICATION      #
###########################################
plt.figure()
# Plot the scatter of timings
plt.plot(loop2[:,1],loop2[:,0], 'o', ms=3, label='Data points')

# Fit a polynomial
coefs2 = poly.polyfit(loop2[:,1],loop2[:,0], 3)
fit_loop2 = poly.polyval(x, coefs2)

# This is for printing the equation
eqstr = "f(x) = "
for i in range(len(coefs)):
    if i == (len(coefs) - 1):
        eqstr += str("{:.2e}".format(coefs2[i]))
    else:
        eqstr += "+"+str("{:.2e}".format(coefs2[i]))+"*x^"+str(len(coefs2)-1-i)
    if i == 1:
        eqstr +="+\n          "
plt.text(0, 4, eqstr, fontdict = {'fontsize' : 12})

plt.plot(x,fit_loop2,label='Fit')
plt.title("Time Scaling (Loop2)",fontsize=16)
plt.legend(prop={'size': 12})
plt.ylabel('t[s]     ').set_rotation(0)
plt.xlabel('N')
plt.savefig('imgs/loop2.svg', format='svg')


# In[6]:


###########################################
#     INTERNAL MATRIX MULTIPLICATION      #
###########################################
plt.figure()
# Plot the scatter of timings
plt.plot(internal[:,1],internal[:,0], 'o', ms=3, label='Data points')

# Fit a polynomial
coefs = poly.polyfit(internal[:,1],internal[:,0], 3)
fit_int = poly.polyval(x, coefs)

# This is for printing the equation
eqstr = "f(x) = "
for i in range(len(coefs)):
    if i == (len(coefs) - 1):
        eqstr += str("{:.2e}".format(coefs[i]))
    else:
        eqstr += "+"+str("{:.2e}".format(coefs[i]))+"*x^"+str(len(coefs)-1-i)
    if i == 1:
        eqstr +="+\n          "
plt.text(0, 0.06, eqstr, fontdict = {'fontsize' : 12})

plt.plot(x,fit_int, label='Fit',color='red')
plt.title("Time Scaling (Internal)",fontsize=16)
plt.legend(prop={'size': 12})
plt.ylabel('t[s]     ').set_rotation(0)
plt.xlabel('N')
plt.savefig('imgs/internal.svg', format='svg')


# In[7]:


###################################################
#               LOOPS SIDE BY SIDE                #
###################################################

plt.subplots(nrows=1, ncols=2, figsize=(15, 3))

plt.subplot(1, 2, 1)
# Plot the histogram of normalized spacing for random hermitian matrices
plt.plot(loop[:,1],loop[:,0], 'o', ms=3, label='Data points')
plt.plot(x,fit_loop,label='Fit')
plt.title("Time Scaling (Loop)",fontsize=16)
plt.legend(prop={'size': 12})
plt.ylabel('t[s]     ').set_rotation(0)
plt.xlabel('N')
plt.ylim(0,8)

plt.subplot(1, 2, 2)
# Plot the histogram of normalized spacing for random diagonal matrices
plt.plot(loop2[:,1],loop2[:,0], 'o', ms=3, label='Data points')
plt.plot(x,fit_loop2,label='Fit',color='magenta')
plt.title("Time Scaling (Loop2)",fontsize=16)
plt.legend(prop={'size': 12})
plt.ylabel('t[s]     ').set_rotation(0)
plt.xlabel('N')
plt.ylim(0,8)

plt.savefig('imgs/loopssidebyside.svg', format='svg')

plt.show()


# In[8]:


###########################################
#                 COMPARE                 #
###########################################
#            Internal vs Loop             #
###########################################
plt.figure()
plt.plot(x,fit_int,label='Internal multiplication',color='red')
plt.plot(x,fit_loop,label='Loop multiplication',color='orange')
plt.legend(prop={'size': 12})
plt.title("Time Scaling",fontsize=16)
plt.ylabel('t[s]     ').set_rotation(0)
plt.xlabel('N')
plt.savefig('imgs/compare.svg', format='svg')


# In[9]:


###########################################
#                 COMPARE                 #
###########################################
#            Loop 1 vs Loop 2             #
###########################################
plt.figure()
plt.plot(x,fit_loop,label='Loop1',color='orange')
plt.plot(x,fit_loop2,label='Loop2',color='magenta')
plt.legend(prop={'size': 12})
plt.title("Time Scaling",fontsize=16)
plt.ylabel('t[s]     ').set_rotation(0)
plt.xlabel('N')
plt.savefig('imgs/compareloops.svg', format='svg')


# In[ ]:




