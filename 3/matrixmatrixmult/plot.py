import matplotlib.pyplot as plt
import numpy as np
import numpy.polynomial.polynomial as poly

# Load the data points from files
loop     = np.genfromtxt('loop.csv', delimiter=',')
internal = np.genfromtxt('internal.csv', delimiter=',')

# Generic dense x vector for plotting lines
x = np.arange(10,1000,1)

###########################################
#        LOOP MATRIX MULTIPLICATION       #
###########################################
plt.figure()
# Plot the scatter of timings
plt.plot(loop[:,1],loop[:,0], 'o', ms=3)

# Fit a polynomial
coefs = poly.polyfit(loop[:,1],loop[:,0], 3)
fit_loop = poly.polyval(x, coefs)

# This is for printing the equation
eqstr = ""
for i in range(len(coefs)):
    if i == (len(coefs) - 1):
        eqstr += str("{:.2e}".format(coefs[i]))
    else:
        eqstr += "+"+str("{:.2e}".format(coefs[i]))+"*x^"+str(len(coefs)-1-i)
    if i == 2:
        eqstr +="\n"
plt.text(0, 6, eqstr)

plt.plot(x,fit_loop)
plt.title("Time Scaling (Loop)")
plt.ylabel('t[s]').set_rotation(0)
plt.xlabel('N')
plt.savefig('loop.png')

###########################################
#     INTERNAL MATRIX MULTIPLICATION      #
###########################################
plt.figure()
# Plot the scatter of timings
plt.plot(internal[:,1],internal[:,0], 'o', ms=3)

# Fit a polynomial
coefs = poly.polyfit(internal[:,1],internal[:,0], 3)
fit_int = poly.polyval(x, coefs)

# This is for printing the equation
eqstr = ""
for i in range(len(coefs)):
    if i == (len(coefs) - 1):
        eqstr += str("{:.2e}".format(coefs[i]))
    else:
        eqstr += "+"+str("{:.2e}".format(coefs[i]))+"*x^"+str(len(coefs)-1-i)
    if i == 2:
        eqstr +="\n"
plt.text(0, .08, eqstr)

plt.plot(x,fit_int)
plt.title("Time Scaling (Internal)")
plt.ylabel('t[s]').set_rotation(0)
plt.xlabel('N')
plt.savefig('internal.png')


###########################################
#                 COMPARE                 #
###########################################
plt.figure()
plt.plot(x,fit_int,label='Internal multiplication')
plt.plot(x,fit_loop,label='Loop multiplication')
plt.legend()
plt.title("Time Scaling")
plt.ylabel('t[s]').set_rotation(0)
plt.xlabel('N')
plt.savefig('compare.png')
