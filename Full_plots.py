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


reading_unc = 0.009 #reding uncertainty for x, y values. Propagated under assumption that reading is accurate up to diameter of tracking tape

D = 0.004 #distance from rope attachment to centre of mass, value is constant over all sized of mass

L_vals = np.array([0.075 + D, 0.10 + D, 0.125 + D, 0.150 + D, 0.20 + D]) #distance from pivot point to COM
L_vals_unc = [0.001, 0.001, 0.001, 0.001, 0.001]

masses = np.array([0.05, 0.1, 0.2, 0.4, 0.5])


############################################ IMPORTING DATA ############################################



############################################ VARYING L

L1_times = np.loadtxt("./Data 2/Varying L (m = 100g)/75mm.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying L (m = 100g)/75mm.txt", skiprows= 2, usecols= 0, unpack= True)[0]
L1_theta = np.arctan(np.loadtxt("./Data 2/Varying L (m = 100g)/75mm.txt",skiprows= 2, usecols= 1, unpack= True)/np.loadtxt("./Data 2/Varying L (m = 100g)/75mm.txt", skiprows= 2, usecols= 2, unpack= True))
L1_x = np.loadtxt("./Data 2/Varying L (m = 100g)/75mm.txt",skiprows= 2, usecols= 1, unpack= True)
L1_y = np.loadtxt("./Data 2/Varying L (m = 100g)/75mm.txt", skiprows= 2, usecols= 2, unpack= True)
L1_unc = np.sqrt(abs((reading_unc * (L1_y / (L1_x ** 2 + L1_y**2)))**2  -  (reading_unc * (L1_x / (L1_x**2  + L1_y**2)))**2))



L2_times = np.loadtxt("./Data 2/Varying L (m = 100g)/100mm.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying L (m = 100g)/100mm.txt", skiprows= 2, usecols= 0, unpack= True)[0]
L2_theta = np.arctan(np.loadtxt("./Data 2/Varying L (m = 100g)/100mm.txt",skiprows= 2, usecols= 1, unpack= True)/np.loadtxt("./Data 2/Varying L (m = 100g)/100mm.txt", skiprows= 2, usecols= 2, unpack= True))
L2_x = np.loadtxt("./Data 2/Varying L (m = 100g)/100mm.txt",skiprows= 2, usecols= 1, unpack= True)
L2_y = np.loadtxt("./Data 2/Varying L (m = 100g)/100mm.txt", skiprows= 2, usecols= 2, unpack= True)
L2_unc = np.sqrt(abs((reading_unc * (abs(L2_y) / (L2_x ** 2 + L2_y**2)))**2  -  (reading_unc * (abs(L2_x) / (L2_x**2  + L2_y**2)))**2))



L3_times = np.loadtxt("./Data 2/Varying L (m = 100g)/125mm.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying L (m = 100g)/125mm.txt", skiprows= 2, usecols= 0, unpack= True)[0]
L3_theta = np.arctan(np.loadtxt("./Data 2/Varying L (m = 100g)/125mm.txt",skiprows= 2, usecols= 1, unpack= True)/np.loadtxt("./Data 2/Varying L (m = 100g)/125mm.txt", skiprows= 2, usecols= 2, unpack= True))
L3_x = np.loadtxt("./Data 2/Varying L (m = 100g)/125mm.txt",skiprows= 2, usecols= 1, unpack= True)
L3_y = np.loadtxt("./Data 2/Varying L (m = 100g)/125mm.txt", skiprows= 2, usecols= 2, unpack= True)
L3_unc = np.sqrt(abs((reading_unc * (abs(L3_y) / (L3_x ** 2 + L3_y**2)))**2  -  (reading_unc * (abs(L3_x) / (L3_x**2  + L3_y**2)))**2))



L4_times = np.loadtxt("./Data 2/Varying L (m = 100g)/150mm.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying L (m = 100g)/150mm.txt", skiprows= 2, usecols= 0, unpack= True)[0]
L4_theta = np.arctan(np.loadtxt("./Data 2/Varying L (m = 100g)/150mm.txt",skiprows= 2, usecols= 1, unpack= True)/np.loadtxt("./Data 2/Varying L (m = 100g)/150mm.txt", skiprows= 2, usecols= 2, unpack= True))
L4_x = np.loadtxt("./Data 2/Varying L (m = 100g)/150mm.txt",skiprows= 2, usecols= 1, unpack= True)
L4_y = np.loadtxt("./Data 2/Varying L (m = 100g)/150mm.txt", skiprows= 2, usecols= 2, unpack= True)
L4_unc = np.sqrt(abs((reading_unc * (abs(L4_y) / (L4_x ** 2 + L4_y**2)))**2  -  (reading_unc * (abs(L4_x) / (L4_x**2  + L4_y**2)))**2))


L5_times = np.loadtxt("./Data 2/Varying L (m = 100g)/200mm.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying L (m = 100g)/200mm.txt", skiprows= 2, usecols= 0, unpack= True)[0]
L5_theta = np.arctan(np.loadtxt("./Data 2/Varying L (m = 100g)/200mm.txt",skiprows= 2, usecols= 1, unpack= True)/np.loadtxt("./Data 2/Varying L (m = 100g)/200mm.txt", skiprows= 2, usecols= 2, unpack= True))
L5_x = np.loadtxt("./Data 2/Varying L (m = 100g)/200mm.txt",skiprows= 2, usecols= 1, unpack= True)
L5_y = np.loadtxt("./Data 2/Varying L (m = 100g)/200mm.txt", skiprows= 2, usecols= 2, unpack= True)
L5_unc = np.sqrt(abs((reading_unc * (abs(L5_y) / (L5_x ** 2 + L5_y**2)))**2  -  (reading_unc * (abs(L5_x) / (L5_x**2  + L5_y**2)))**2))



############################################ VARYING M

M1_times = np.loadtxt("./Data 2/Varying m (L = 100mm)/50g.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying m (L = 100mm)/50g.txt", skiprows= 2, usecols= 0, unpack= True)[0]
M1_theta = np.arctan(np.loadtxt("./Data 2/Varying m (L = 100mm)/50g.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying m (L = 100mm)/50g.txt", skiprows= 2, usecols= 2, unpack= True) )
M1_x = np.loadtxt("./Data 2/Varying m (L = 100mm)/50g.txt",skiprows= 2, usecols= 1, unpack= True)
M1_y = np.loadtxt("./Data 2/Varying m (L = 100mm)/50g.txt", skiprows= 2, usecols= 2, unpack= True)
M1_unc = np.sqrt(abs((reading_unc * (M1_y / (M1_x ** 2 + M1_y**2)))**2  -  (reading_unc * (M1_x / (M1_x**2  + M1_y**2)))**2))



M2_times = np.loadtxt("./Data 2/Varying m (L = 100mm)/100g.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying m (L = 100mm)/100g.txt", skiprows= 2, usecols= 0, unpack= True)[0]
M2_theta = np.arctan(np.loadtxt("./Data 2/Varying m (L = 100mm)/100g.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying m (L = 100mm)/100g.txt", skiprows= 2, usecols= 2, unpack= True) )
M2_x = np.loadtxt("./Data 2/Varying m (L = 100mm)/100g.txt",skiprows= 2, usecols= 1, unpack= True)
M2_y = np.loadtxt("./Data 2/Varying m (L = 100mm)/100g.txt", skiprows= 2, usecols= 2, unpack= True)
M2_unc = np.sqrt(abs((reading_unc * (M2_y / (M2_x ** 2 + M2_y**2)))**2  -  (reading_unc * (M2_x / (M2_x**2  + M2_y**2)))**2))


M3_times = np.loadtxt("./Data 2/Varying m (L = 100mm)/200g.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying m (L = 100mm)/200g.txt", skiprows= 2, usecols= 0, unpack= True)[0]
M3_theta = np.arctan(np.loadtxt("./Data 2/Varying m (L = 100mm)/200g.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying m (L = 100mm)/200g.txt", skiprows= 2, usecols= 2, unpack= True) )
M3_x = np.loadtxt("./Data 2/Varying m (L = 100mm)/200g.txt",skiprows= 2, usecols= 1, unpack= True)
M3_y = np.loadtxt("./Data 2/Varying m (L = 100mm)/200g.txt", skiprows= 2, usecols= 2, unpack= True)
M3_unc = np.sqrt(abs((reading_unc * (M3_y / (M3_x ** 2 + M3_y**2)))**2  -  (reading_unc * (M3_x / (M3_x**2  + M3_y**2)))**2))


M4_times = np.loadtxt("./Data 2/Varying m (L = 100mm)/400g.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying m (L = 100mm)/400g.txt", skiprows= 2, usecols= 0, unpack= True)[0]
M4_theta = np.arctan(np.loadtxt("./Data 2/Varying m (L = 100mm)/400g.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying m (L = 100mm)/400g.txt", skiprows= 2, usecols= 2, unpack= True) )
M4_x = np.loadtxt("./Data 2/Varying m (L = 100mm)/400g.txt",skiprows= 2, usecols= 1, unpack= True)
M4_y = np.loadtxt("./Data 2/Varying m (L = 100mm)/400g.txt", skiprows= 2, usecols= 2, unpack= True)
M4_unc = np.sqrt(abs((reading_unc * (M4_y / (M4_x ** 2 + M4_y**2)))**2  -  (reading_unc * (M4_x / (M4_x**2  + M4_y**2)))**2))


M5_times = np.loadtxt("./Data 2/Varying m (L = 100mm)/500g.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying m (L = 100mm)/500g.txt", skiprows= 2, usecols= 0, unpack= True)[0]
M5_theta = np.arctan(np.loadtxt("./Data 2/Varying m (L = 100mm)/500g.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying m (L = 100mm)/500g.txt", skiprows= 2, usecols= 2, unpack= True) )
M5_x = np.loadtxt("./Data 2/Varying m (L = 100mm)/500g.txt",skiprows= 2, usecols= 1, unpack= True)
M5_y = np.loadtxt("./Data 2/Varying m (L = 100mm)/500g.txt", skiprows= 2, usecols= 2, unpack= True)
M5_unc = np.sqrt(abs((reading_unc * (M5_y / (M5_x ** 2 + M5_y**2)))**2  -  (reading_unc * (M5_x / (M5_x**2  + M5_y**2)))**2))



M6_times = np.loadtxt("./Data 2/Varying m (L = 100mm)/150g.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying m (L = 100mm)/150g.txt", skiprows= 2, usecols= 0, unpack= True)[0]
M6_theta = np.arctan(np.loadtxt("./Data 2/Varying m (L = 100mm)/150g.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying m (L = 100mm)/150g.txt", skiprows= 2, usecols= 2, unpack= True) )
M6_x = np.loadtxt("./Data 2/Varying m (L = 100mm)/150g.txt",skiprows= 2, usecols= 1, unpack= True)
M6_y = np.loadtxt("./Data 2/Varying m (L = 100mm)/150g.txt", skiprows= 2, usecols= 2, unpack= True)
M6_unc = np.sqrt(abs((reading_unc * (M6_y / (M6_x ** 2 + M6_y**2)))**2  -  (reading_unc * (M6_x / (M6_x**2  + M6_y**2)))**2))


M7_times = np.loadtxt("./Data 2/Varying m (L = 100mm)/250g.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying m (L = 100mm)/250g.txt", skiprows= 2, usecols= 0, unpack= True)[0]
M7_theta = np.arctan(np.loadtxt("./Data 2/Varying m (L = 100mm)/250g.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying m (L = 100mm)/250g.txt", skiprows= 2, usecols= 2, unpack= True) )
M7_x = np.loadtxt("./Data 2/Varying m (L = 100mm)/250g.txt",skiprows= 2, usecols= 1, unpack= True)
M7_y = np.loadtxt("./Data 2/Varying m (L = 100mm)/250g.txt", skiprows= 2, usecols= 2, unpack= True)
M7_unc = np.sqrt(abs((reading_unc * (M7_y / (M7_x ** 2 + M7_y**2)))**2  -  (reading_unc * (M7_x / (M7_x**2  + M7_y**2)))**2))


M8_times = np.loadtxt("./Data 2/Varying m (L = 100mm)/300g.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying m (L = 100mm)/300g.txt", skiprows= 2, usecols= 0, unpack= True)[0]
M8_theta = np.arctan(np.loadtxt("./Data 2/Varying m (L = 100mm)/300g.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying m (L = 100mm)/300g.txt", skiprows= 2, usecols= 2, unpack= True) )
M8_x = np.loadtxt("./Data 2/Varying m (L = 100mm)/300g.txt",skiprows= 2, usecols= 1, unpack= True)
M8_y = np.loadtxt("./Data 2/Varying m (L = 100mm)/300g.txt", skiprows= 2, usecols= 2, unpack= True)
M8_unc = np.sqrt(abs((reading_unc * (M8_y / (M8_x ** 2 + M8_y**2)))**2  -  (reading_unc * (M8_x / (M8_x**2  + M8_y**2)))**2))

M9_times = np.loadtxt("./Data 2/Varying m (L = 100mm)/140g.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying m (L = 100mm)/140g.txt", skiprows= 2, usecols= 0, unpack= True)[0]
M9_theta = np.arctan(np.loadtxt("./Data 2/Varying m (L = 100mm)/140g.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying m (L = 100mm)/140g.txt", skiprows= 2, usecols= 2, unpack= True) )
M9_x = np.loadtxt("./Data 2/Varying m (L = 100mm)/140g.txt",skiprows= 2, usecols= 1, unpack= True)
M9_y = np.loadtxt("./Data 2/Varying m (L = 100mm)/140g.txt", skiprows= 2, usecols= 2, unpack= True)
M9_unc = np.sqrt(abs((reading_unc * (M9_y / (M9_x ** 2 + M9_y**2)))**2  -  (reading_unc * (M9_x / (M9_x**2  + M9_y**2)))**2))


############################################ VARYING THETA


P1_times = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta0.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta0.txt", skiprows= 2, usecols= 0, unpack= True)[0]
P1_theta = np.arctan(np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta0.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta0.txt", skiprows= 2, usecols= 2, unpack= True) )
P1_x = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta0.txt", skiprows= 2, usecols= 1, unpack= True)
P1_y = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta0.txt", skiprows= 2, usecols= 2, unpack= True) 
P1_unc = np.sqrt(abs((reading_unc * (P1_y / (P1_x ** 2 + P1_y**2)))**2  -  (reading_unc * (P1_x / (P1_x**2  + P1_y**2)))**2))



P3_times = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta1.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta1.txt", skiprows= 2, usecols= 0, unpack= True)[0]
P3_theta = np.arctan(np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta1.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta1.txt", skiprows= 2, usecols= 2, unpack= True) )
P3_x = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta1.txt", skiprows= 2, usecols= 1, unpack= True)
P3_y = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta1.txt", skiprows= 2, usecols= 2, unpack= True) 
P3_unc = np.sqrt(abs((reading_unc * (P3_y / (P3_x ** 2 + P3_y**2)))**2  -  (reading_unc * (P3_x / (P3_x**2  + P3_y**2)))**2))



P2_times = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta2.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta2.txt", skiprows= 2, usecols= 0, unpack= True)[0]
P2_theta = np.arctan(np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta2.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta2.txt", skiprows= 2, usecols= 2, unpack= True) )
P2_x = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta2.txt", skiprows= 2, usecols= 1, unpack= True)
P2_y = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta2.txt", skiprows= 2, usecols= 2, unpack= True) 
P2_unc = np.sqrt(abs((reading_unc * (P2_y / (P2_x ** 2 + P2_y**2)))**2  -  (reading_unc * (P2_x / (P2_x**2  + P2_y**2)))**2))



P4_times = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta3.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta3.txt", skiprows= 2, usecols= 0, unpack= True)[0]
P4_theta = np.arctan(np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta3.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta3.txt", skiprows= 2, usecols= 2, unpack= True) )
P4_x = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta3.txt", skiprows= 2, usecols= 1, unpack= True)
P4_y = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta3.txt", skiprows= 2, usecols= 2, unpack= True) 
P4_unc = np.sqrt((reading_unc / P4_x) ** 2 + (reading_unc / P4_y) ** 2)
P4_unc = np.sqrt(abs((reading_unc * (P4_y / (P4_x ** 2 + P4_y**2)))**2  -  (reading_unc * (P4_x / (P4_x**2  + P4_y**2)))**2))


P5_times = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta4.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta4.txt", skiprows= 2, usecols= 0, unpack= True)[0]
P5_theta = np.arctan(np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta4.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta4.txt", skiprows= 2, usecols= 2, unpack= True) )
P5_x = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta4.txt", skiprows= 2, usecols= 1, unpack= True)
P5_y = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta4.txt", skiprows= 2, usecols= 2, unpack= True) 
P5_unc = np.sqrt(abs((reading_unc * (P5_y / (P5_x ** 2 + P5_y**2)))**2  -  (reading_unc * (P5_x / (P5_x**2  + P5_y**2)))**2))



P6_times = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta5.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta5.txt", skiprows= 2, usecols= 0, unpack= True)[0]
P6_theta = np.arctan(np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta5.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta5.txt", skiprows= 2, usecols= 2, unpack= True) )
P6_x = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta5.txt", skiprows= 2, usecols= 1, unpack= True)
P6_y = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta5.txt", skiprows= 2, usecols= 2, unpack= True) 
P6_unc = np.sqrt(abs((reading_unc * (P6_y / (P6_x ** 2 + P6_y**2)))**2  +  (reading_unc * (P6_x / (P6_x**2  + P6_y**2)))**2))



############################################ FUNCTIONAL FORMS ############################################

thetas = np.array([P1_theta[0], P2_theta[0], P3_theta[0], P4_theta[0], P5_theta[0], P6_theta[0]])

def theta(t, A, T, tau, phi):
    return A * np.exp(-t/tau) * np.cos((2 * np.pi * t / T) + phi)

def linear(x, A, B):
    return A*x + B

def quad(x, A, B, C):
    return A*x**2 + B*x + C

def sqrt(x, A, B):
    return A * np.sqrt(x) + B



############################################ CURVE FITTING DATA ############################################


popt_L1, pcov_L1 = curve_fit(theta, L1_times, L1_theta, sigma = L1_unc, p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_L2, pcov_L2 = curve_fit(theta, L2_times, L2_theta, sigma = L2_unc, p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_L3, pcov_L3 = curve_fit(theta, L3_times, L3_theta, sigma = L3_unc, p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_L4, pcov_L4 = curve_fit(theta, L4_times, L4_theta, sigma = L4_unc, p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_L5, pcov_L5 = curve_fit(theta, L5_times, L5_theta, sigma = L5_unc, p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)


popt_M1, pcov_M1 = curve_fit(theta, M1_times, M1_theta, sigma = M1_unc, p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_M2, pcov_M2 = curve_fit(theta, M2_times, M2_theta, sigma = M2_unc, p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_M3, pcov_M3 = curve_fit(theta, M3_times, M3_theta, sigma = M3_unc, p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_M4, pcov_M4 = curve_fit(theta, M4_times, M4_theta, sigma = M4_unc, p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_M5, pcov_M5 = curve_fit(theta, M5_times, M5_theta, sigma = M5_unc, p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)



popt_P1, pcov_P1 = curve_fit(theta, P1_times, P1_theta, sigma = P1_unc, p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_P2, pcov_P2 = curve_fit(theta, P2_times, P2_theta, sigma = P2_unc, p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_P3, pcov_P3 = curve_fit(theta, P3_times, P3_theta, sigma = P3_unc, p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_P4, pcov_P4 = curve_fit(theta, P4_times, P4_theta, sigma = P4_unc, p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_P5, pcov_P5 = curve_fit(theta, P5_times, P5_theta, sigma = P5_unc, p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_P6, pcov_P6 = curve_fit(theta, P6_times, P6_theta, sigma = P6_unc, p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)


print(popt_L1[2]/popt_L1[1], popt_L2[2]/popt_L2[1], popt_L3[2]/popt_L3[1], popt_L4[2]/popt_L4[1], popt_L5[2]/popt_L5[1])
print(np.sqrt(np.diag(pcov_P1)[2]), np.sqrt(np.diag(pcov_P2)[2]), np.sqrt(np.diag(pcov_P3)[2]), np.sqrt(np.diag(pcov_P4)[2]), np.sqrt(np.diag(pcov_P5)[2]), np.sqrt(np.diag(pcov_P6)[2]))

average = np.average([popt_L1[2]/popt_L1[1], popt_L2[2]/popt_L2[1], popt_L3[2]/popt_L3[1], popt_L4[2]/popt_L4[1], popt_L5[2]/popt_L5[1],\
                     popt_M1[2]/popt_M1[1], popt_M2[2]/popt_M2[1], popt_M3[2]/popt_M3[1], popt_M4[2]/popt_M4[1], popt_M5[2]/popt_M5[1],\
                        popt_P1[2]/popt_M1[1], popt_P2[2]/popt_P2[1], popt_P3[2]/popt_P3[1]])

print(average)

chisquared_L1 = np.sum(((L1_theta - theta(L1_times, *popt_L1))/L1_unc )**2) / (len(L1_theta) - len(popt_L1))
chisquared_L2 = np.sum(((L2_theta - theta(L2_times, *popt_L2))/L2_unc )**2) / (len(L2_theta) - len(popt_L2))
chisquared_L3 = np.sum(((L3_theta - theta(L3_times, *popt_L3))/L3_unc )**2) / (len(L3_theta) - len(popt_L3))
chisquared_L4 = np.sum(((L4_theta - theta(L4_times, *popt_L4))/L4_unc )**2) / (len(L4_theta) - len(popt_L4))
chisquared_L5 = np.sum(((L5_theta - theta(L5_times, *popt_L5))/L5_unc )**2) / (len(L5_theta) - len(popt_L5))

chisquared_M1 = np.sum(((M1_theta - theta(M1_times, *popt_M1))/M1_unc )**2) / (len(M1_theta) - len(popt_M1))
chisquared_M2 = np.sum(((M2_theta - theta(M2_times, *popt_M2))/M2_unc )**2) / (len(M2_theta) - len(popt_M2))
chisquared_M3 = np.sum(((M3_theta - theta(M3_times, *popt_M3))/M3_unc )**2) / (len(M3_theta) - len(popt_M3))
chisquared_M4 = np.sum(((M4_theta - theta(M4_times, *popt_M4))/M4_unc )**2) / (len(M4_theta) - len(popt_M4))
chisquared_M5 = np.sum(((M5_theta - theta(M5_times, *popt_M5))/M5_unc )**2) / (len(M5_theta) - len(popt_M5))

chisquared_P1 = np.sum(((P1_theta - theta(P1_times, *popt_P1))/P1_unc )**2) / (len(P1_theta) - len(popt_P1))
chisquared_P2 = np.sum(((P2_theta - theta(P2_times, *popt_P2))/P2_unc )**2) / (len(P2_theta) - len(popt_P2))
chisquared_P3 = np.sum(((P3_theta - theta(P3_times, *popt_P3))/P3_unc )**2) / (len(P3_theta) - len(popt_P3))
chisquared_P4 = np.sum(((P4_theta - theta(P4_times, *popt_P4))/P4_unc )**2) / (len(P4_theta) - len(popt_P4))
chisquared_P5 = np.sum(((P5_theta - theta(P5_times, *popt_P5))/P5_unc )**2) / (len(P5_theta) - len(popt_P5))
chisquared_P6 = np.sum(((P6_theta - theta(P6_times, *popt_P6))/P6_unc )**2) / (len(P6_theta) - len(popt_P6))


############################################ PLOTTING ############################################

zeros = np.zeros(5000)

############################################ VARYING LENGTH

'''

plt.plot(L1_times, theta(L1_times, *popt_L1), c = 'red', label = 'Functional Fit')
plt.scatter(L1_times, L1_theta, s = 2, c = 'green', label = 'Collected Data Residual, L1')
plt.text(50, 0.15, r'$\chi^2/DOF = $ %3.2f'%(chisquared_L1), fontsize=fontsize, c = 'red', alpha = myalpha)
plt.xlabel('Time (s)')
plt.ylabel('Angle (rads)')
plt.legend(loc = 'upper right')
plt.tight_layout()
if (save): plt.savefig('L1.png',dpi=mydpi)
plt.show()



plt.plot(L2_times, theta(L2_times, *popt_L2), c = 'red', label = 'Functional Fit')
plt.scatter(L2_times, L2_theta, s = 2, c = 'green', label = 'Collected Data, L2')
plt.text(55, 0.14, r'$\chi^2/DOF = $ %3.2f'%(chisquared_L2), fontsize=fontsize, c = 'red', alpha = myalpha)
plt.xlabel('Time (s)')
plt.ylabel('Angle (rads)')
plt.legend(loc = 'upper right')
plt.tight_layout()
if (save): plt.savefig('L2.png',dpi=mydpi)
plt.show()


plt.plot(L3_times, theta(L3_times, *popt_L3), c = 'red', label = 'Functional Fit')
plt.scatter(L3_times, L3_theta, s = 2, c = 'green', label = 'Collected Data, L3')
plt.text(50, -0.2, r'$\chi^2/DOF = $ %3.2f'%(chisquared_L3), fontsize=fontsize, c = 'red', alpha = myalpha)
plt.xlabel('Time (s)')
plt.ylabel('Angle (rads)')
plt.legend(loc = 'upper right')
plt.tight_layout()
if (save): plt.savefig('L3.png',dpi=mydpi)
plt.show()


plt.plot(L4_times, theta(L4_times, *popt_L4), c = 'red', label = 'Functional Fit')
plt.scatter(L4_times, L4_theta, s = 2, c = 'green', label = 'Collected Data, L4')
plt.text(45, -0.2, r'$\chi^2/DOF = $ %3.2f'%(chisquared_L4), fontsize=fontsize, c = 'red', alpha = myalpha)
plt.xlabel('Time (s)')
plt.ylabel('Angle (rads)')
plt.legend(loc = 'upper right')
plt.tight_layout()
if (save): plt.savefig('L4.png',dpi=mydpi)
plt.show()


plt.plot(L5_times, theta(L5_times, *popt_L5), c = 'red', label = 'Functional Fit')
plt.scatter(L5_times, L5_theta, s = 2, c = 'green', label = 'Collected Data, L5')
plt.text(55, -0.3, r'$\chi^2/DOF = $ %3.2f'%(chisquared_L5), fontsize=fontsize, c = 'red', alpha = myalpha)
plt.xlabel('Time (s)')
plt.ylabel('Angle (rads)')
plt.legend(loc = 'upper right')
plt.tight_layout()
if (save): plt.savefig('L5.png',dpi=mydpi)
plt.show()



'''


############################################ VARYING MASS



'''


plt.plot(M1_times, theta(M1_times, *popt_M1), c = 'red', label = 'Functional Fit')
plt.scatter(M1_times, M1_theta, s = 2, c = 'green', label = 'Collected Data, M1')
plt.text(45, 0.13, r'$\chi^2/DOF = $ %3.2f'%(chisquared_M1), fontsize=fontsize, c = 'red', alpha = myalpha)
plt.xlabel('Time (s)')
plt.ylabel('Angle (rads)')
plt.legend(loc = 'upper right')
plt.tight_layout()
if (save): plt.savefig('M1.png',dpi=mydpi)
plt.show()


plt.plot(M2_times, theta(M2_times, *popt_M2), c = 'red', label = 'Functional Fit')
plt.scatter(M2_times, M2_theta, s = 2, c = 'green', label = 'Collected Data, M2')
plt.text(50, 0.13, r'$\chi^2/DOF = $ %3.2f'%(chisquared_M2), fontsize=fontsize, c = 'red', alpha = myalpha)
plt.xlabel('Time (s)')
plt.ylabel('Angle (rads)')
plt.legend(loc = 'upper right')
plt.tight_layout()
if (save): plt.savefig('M2.png',dpi=mydpi)
plt.show()


plt.plot(M3_times, theta(M3_times, *popt_M3), c = 'red', label = 'Functional Fit')
plt.scatter(M3_times, M3_theta, s = 2, c = 'green', label = 'Collected Data, M3')
plt.text(45, -0.2, r'$\chi^2/DOF = $ %3.2f'%(chisquared_M3), fontsize=fontsize, c = 'red', alpha = myalpha)
plt.xlabel('Time (s)')
plt.ylabel('Angle (rads)')
plt.legend(loc = 'upper right')
plt.tight_layout()
if (save): plt.savefig('M3.png',dpi=mydpi)
plt.show()


plt.plot(M4_times, theta(M4_times, *popt_M4), c = 'red', label = 'Functional Fit')
plt.scatter(M4_times, M4_theta, s = 2, c = 'green', label = 'Collected Data, M4')
plt.text(50, -0.22, r'$\chi^2/DOF = $ %3.2f'%(chisquared_M4), fontsize=fontsize, c = 'red', alpha = myalpha)
plt.xlabel('Time (s)')
plt.ylabel('Angle (rads)')
plt.legend(loc = 'upper right')
plt.tight_layout()
if (save): plt.savefig('M4.png',dpi=mydpi)
plt.show()


plt.plot(M5_times, theta(M5_times, *popt_M5), c = 'red', label = 'Functional Fit')
plt.scatter(M5_times, M5_theta, s = 2, c = 'green', label = 'Collected Data, M5')
plt.text(50, -0.2, r'$\chi^2/DOF = $ %3.2f'%(chisquared_M5), fontsize=fontsize, c = 'red', alpha = myalpha)
plt.xlabel('Time (s)')
plt.ylabel('Angle (rads)')
plt.legend(loc = 'upper right')
plt.tight_layout()
if (save): plt.savefig('M5.png',dpi=mydpi)
plt.show()


'''


############################################ VARYING THETA

'''


plt.plot(P1_times, theta(P1_times, *popt_P1), c = 'red', label = 'Functional Fit')
plt.scatter(P1_times, P1_theta, s = 2, c = 'green', label = 'Collected Data, P1')
plt.text(50, 0.08, r'$\chi^2/DOF = $ %3.2f'%(chisquared_P1), fontsize=fontsize, c = 'red', alpha = myalpha)
plt.xlabel('Time (s)')
plt.ylabel('Angle (rads)')
plt.legend(loc = 'upper right')
plt.tight_layout()
if (save): plt.savefig('P1.png',dpi=mydpi)
plt.show()


plt.plot(P2_times, theta(P2_times, *popt_P2), c = 'red', label = 'Functional Fit')
plt.scatter(P2_times, P2_theta, s = 2, c = 'green', label = 'Collected Data, P2')
plt.text(50, 0.1, r'$\chi^2/DOF = $ %3.2f'%(chisquared_P2), fontsize=fontsize, c = 'red', alpha = myalpha)
plt.xlabel('Time (s)')
plt.ylabel('Angle (rads)')
plt.legend(loc = 'upper right')
plt.tight_layout()
if (save): plt.savefig('P2.png',dpi=mydpi)
plt.show()


plt.plot(P3_times, theta(P3_times, *popt_P3), c = 'red', label = 'Functional Fit')
plt.scatter(P3_times, P3_theta, s = 2, c = 'green', label = 'Collected Data, P3')
plt.text(45, 0.13, r'$\chi^2/DOF = $ %3.2f'%(chisquared_P3), fontsize=fontsize, c = 'red', alpha = myalpha)
plt.xlabel('Time (s)')
plt.ylabel('Angle (rads)')
plt.legend(loc = 'upper right')
plt.tight_layout()
if (save): plt.savefig('P3.png',dpi=mydpi)
plt.show()


plt.plot(P4_times, theta(P4_times, *popt_P4), c = 'red', label = 'Functional Fit')
plt.scatter(P4_times, P4_theta, s = 2, c = 'green', label = 'Collected Data, P4')
plt.text(50, 0.25, r'$\chi^2/DOF = $ %3.2f'%(chisquared_P4), fontsize=fontsize, c = 'red', alpha = myalpha)
plt.xlabel('Time (s)')
plt.ylabel('Angle (rads)')
plt.legend(loc = 'upper right')
plt.tight_layout()
if (save): plt.savefig('P4.png',dpi=mydpi)
plt.show()


plt.plot(P5_times, theta(P5_times, *popt_P5), c = 'red', label = 'Functional Fit, P5')
plt.scatter(P5_times, P5_theta, s = 2, c = 'green', label = 'Collected Data')
plt.text(50, 0.4, r'$\chi^2/DOF = $ %3.2f'%(chisquared_P5), fontsize=fontsize, c = 'red', alpha = myalpha)
plt.xlabel('Time (s)')
plt.ylabel('Angle (rads)')
plt.legend(loc = 'upper right')
plt.tight_layout()
if (save): plt.savefig('P5.png',dpi=mydpi)
plt.show()





plt.plot(P6_times, theta(P6_times, *popt_P6), c = 'red', label = 'Functional Fit, P6')
plt.scatter(P6_times, P6_theta, s = 2, c = 'green', label = 'Collected Data')
plt.text(45, 0.6, r'$\chi^2/DOF = $ %3.2f'%(chisquared_P6), fontsize=fontsize, c = 'red', alpha = myalpha)
plt.xlabel('Time (s)')
plt.ylabel('Angle (rads)')
plt.legend(loc = 'upper right')
plt.tight_layout()
if (save): plt.savefig('P6.png',dpi=mydpi)
plt.show()


'''