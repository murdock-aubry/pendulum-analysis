import matplotlib.pyplot as plt
import os
from matplotlib import rc
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2

save=True # if True then we save images as files
font = {'size'   : 4}
rc('font', **font)
mydpi = 160
myalpha = 1

font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 13}
rc('font', **font)

fontsize = 12

############################################ CONSTANTS AND UNCERTAINTIES ############################################


reading_unc = 0.003 #reding uncertainty for x, y values. Propagated under assumption that reading is accurate up to diameter of tracking tape

D = 0.004 #distance from rope attachment to centre of mass, value is constant over all sized of mass

L_vals = np.array([0.075 + D, 0.10 + D, 0.125 + D, 0.150 + D, 0.20 + D, 0.175 + D, 0.19 + D]) #distance from pivot point to COM
L_vals_unc = [0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002]

masses = np.array([0.05, 0.1, 0.2, 0.4, 0.5, 0.15, 0.25, 0.3, 0.14])


############################################ IMPORTING DATA ############################################

interval_min = 10
interval_max = 100

P2_times = (np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta1.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta1.txt", skiprows= 2, usecols= 0, unpack= True)[0])[interval_min: interval_max]
P2_theta = (np.arctan(np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta1.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta1.txt", skiprows= 2, usecols= 2, unpack= True) ))[interval_min: interval_max]
P2_x = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta1.txt", skiprows= 2, usecols= 1, unpack= True)[interval_min: interval_max]
P2_y = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta1.txt", skiprows= 2, usecols= 2, unpack= True)[interval_min: interval_max]
P2_unc = np.sqrt(abs((reading_unc * (P2_y / (P2_x ** 2 + P2_y**2)))**2  -  (reading_unc * (P2_x / (P2_x**2  + P2_y**2)))**2))/2




def theta(t, A, T, tau, phi):
    return A * np.exp(-t/tau) * np.cos((2 * np.pi * t / T) + phi)




popt_P2, pcov_P2 = curve_fit(theta, P2_times, P2_theta, sigma = P2_unc, p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
chisquared_1 = np.sum(((P2_theta - theta(P2_times, *popt_P2))/P2_unc )**2) / (len(P2_theta) - len(popt_P2))
chisquared_2 = np.sum(((P2_theta + theta(P2_times - popt_P2[1]/2, *popt_P2))/P2_unc )**2) / (len(P2_theta) - len(popt_P2))
print(chisquared_1, chisquared_2)




period = popt_P2[1]


times_trans = P2_times + period/2

times = np.linspace(min(P2_times), max(P2_times[interval_min: interval_max]), 2000)


#Original Fit


plt.plot(times, theta(times, *popt_P2), c = 'blue', label = 'Functional Form')
plt.scatter(P2_times, P2_theta, c = 'red', s = 5, label = 'Collected Data')
plt.text(2.5, 0.11, r'$\chi^2/DOF = $ %3.2f'%(chisquared_1), fontsize=fontsize, c = 'red', alpha = myalpha)
plt.xlabel(r'Time [s]')
plt.ylabel(r'Angular Amplitude [rads]')
plt.legend(loc = 'upper right')
plt.tight_layout()
if (save): plt.savefig('symmetry_plot.png',dpi=mydpi)
plt.show()


#Shifted fit


plt.plot(times, -theta(times + popt_P2[1]/2, popt_P2[0], popt_P2[1], popt_P2[2], popt_P2[3]), c = 'orange', label = 'Shifted Functional Form')
plt.scatter(P2_times, P2_theta, c = 'red', s = 5, label = 'Collected Data')
plt.text(2.5, 0.11, r'$\chi^2/DOF = $ %3.2f'%(chisquared_2), fontsize=fontsize, c = 'red', alpha = myalpha)
plt.xlabel(r'Time [s]')
plt.ylabel(r'Angular Amplitude [rads]')
plt.legend(loc = 'upper right')
plt.tight_layout()
if (save): plt.savefig('symmetry_shifted.png',dpi=mydpi)
plt.show()


average1 = np.average(abs(P2_theta - theta(P2_times, *popt_P2)))

average2 = np.average(abs(P2_theta + theta(P2_times + popt_P2[1]/2, popt_P2[0], popt_P2[1], popt_P2[2], popt_P2[3])))




plt.plot(times, theta(times, 0, 1, 1, 1), c = 'blue', label = 'Functional Form')
plt.errorbar(P2_times, P2_theta - theta(P2_times, *popt_P2), yerr = P2_unc, capsize = 2, c = 'green', fmt='o', label = 'Collected Data Residual')
plt.xlabel(r'Time [s]')
plt.ylabel(r'Angular Amplitude [rads]')
plt.legend(loc = 'lower right')
plt.tight_layout()
if (save): plt.savefig('symmetry_plot_residual.png',dpi=mydpi)
plt.show()


plt.plot(times, theta(times, 0, 1, 1, 1), c = 'orange', label = 'Shifted Functional Form')
plt.errorbar(P2_times, P2_theta + theta(P2_times - popt_P2[1]/2, popt_P2[0], popt_P2[1], popt_P2[2], popt_P2[3]), yerr = P2_unc, capsize = 2, c = 'green', fmt='o', label = 'Collected Data Residual')
plt.xlabel(r'Time [s]')
plt.ylabel(r'Angular Amplitude [rads]')
plt.legend(loc = 'lower right')
plt.tight_layout()
if (save): plt.savefig('symmetry_shifted_residual.png',dpi=mydpi)
plt.show()


