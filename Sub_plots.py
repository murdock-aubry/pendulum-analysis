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

L_vals = np.array([0.075 + D, 0.10 + D, 0.125 + D, 0.150 + D, 0.20 + D, 0.175 + D, 0.19 + D]) #distance from pivot point to COM
L_vals_unc = [0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002]

masses = np.array([0.05, 0.1, 0.2, 0.4, 0.5, 0.15, 0.25, 0.3, 0.14])


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


L6_times = np.loadtxt("./Data 2/Varying L (m = 100g)/175mm.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying L (m = 100g)/175mm.txt", skiprows= 2, usecols= 0, unpack= True)[0]
L6_theta = np.arctan(np.loadtxt("./Data 2/Varying L (m = 100g)/175mm.txt",skiprows= 2, usecols= 1, unpack= True)/np.loadtxt("./Data 2/Varying L (m = 100g)/175mm.txt", skiprows= 2, usecols= 2, unpack= True))
L6_x = np.loadtxt("./Data 2/Varying L (m = 100g)/175mm.txt",skiprows= 2, usecols= 1, unpack= True)
L6_y = np.loadtxt("./Data 2/Varying L (m = 100g)/175mm.txt", skiprows= 2, usecols= 2, unpack= True)
L6_unc = np.sqrt(abs((reading_unc * (abs(L6_y) / (L6_x ** 2 + L6_y**2)))**2  -  (reading_unc * (abs(L6_x) / (L6_x**2  + L6_y**2)))**2))


L7_times = np.loadtxt("./Data 2/Varying L (m = 100g)/190mm.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying L (m = 100g)/190mm.txt", skiprows= 2, usecols= 0, unpack= True)[0]
L7_theta = np.arctan(np.loadtxt("./Data 2/Varying L (m = 100g)/190mm.txt",skiprows= 2, usecols= 1, unpack= True)/np.loadtxt("./Data 2/Varying L (m = 100g)/190mm.txt", skiprows= 2, usecols= 2, unpack= True))
L7_x = np.loadtxt("./Data 2/Varying L (m = 100g)/190mm.txt",skiprows= 2, usecols= 1, unpack= True)
L7_y = np.loadtxt("./Data 2/Varying L (m = 100g)/190mm.txt", skiprows= 2, usecols= 2, unpack= True)
L7_unc = np.sqrt(abs((reading_unc * (abs(L7_y) / (L7_x ** 2 + L7_y**2)))**2  -  (reading_unc * (abs(L7_x) / (L7_x**2  + L7_y**2)))**2))


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



P2_times = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta1.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta1.txt", skiprows= 2, usecols= 0, unpack= True)[0]
P2_theta = np.arctan(np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta1.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta1.txt", skiprows= 2, usecols= 2, unpack= True) )
P2_x = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta1.txt", skiprows= 2, usecols= 1, unpack= True)
P2_y = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta1.txt", skiprows= 2, usecols= 2, unpack= True) 
P2_unc = np.sqrt(abs((reading_unc * (P2_y / (P2_x ** 2 + P2_y**2)))**2  -  (reading_unc * (P2_x / (P2_x**2  + P2_y**2)))**2))/2



P3_times = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta2.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta2.txt", skiprows= 2, usecols= 0, unpack= True)[0]
P3_theta = np.arctan(np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta2.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta2.txt", skiprows= 2, usecols= 2, unpack= True) )
P3_x = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta2.txt", skiprows= 2, usecols= 1, unpack= True)
P3_y = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta2.txt", skiprows= 2, usecols= 2, unpack= True) 
P3_unc = np.sqrt(abs((reading_unc * (P3_y / (P3_x ** 2 + P3_y**2)))**2  -  (reading_unc * (P3_x / (P3_x**2  + P3_y**2)))**2))



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


P7_times = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta6.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta6.txt", skiprows= 2, usecols= 0, unpack= True)[0]
P7_theta = np.arctan(np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta6.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta6.txt", skiprows= 2, usecols= 2, unpack= True) )
P7_x = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta6.txt", skiprows= 2, usecols= 1, unpack= True)
P7_y = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta6.txt", skiprows= 2, usecols= 2, unpack= True) 
P7_unc = np.sqrt(abs((reading_unc * (P7_y / (P7_x ** 2 + P7_y**2)))**2  +  (reading_unc * (P7_x / (P7_x**2  + P7_y**2)))**2))


P8_times = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta7.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta7.txt", skiprows= 2, usecols= 0, unpack= True)[0]
P8_theta = np.arctan(np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta7.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta7.txt", skiprows= 2, usecols= 2, unpack= True) )
P8_x = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta7.txt", skiprows= 2, usecols= 1, unpack= True)
P8_y = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta7.txt", skiprows= 2, usecols= 2, unpack= True) 
P8_unc = np.sqrt(abs((reading_unc * (P8_y / (P8_x ** 2 + P8_y**2)))**2  +  (reading_unc * (P8_x / (P8_x**2  + P8_y**2)))**2))


P9_times = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta8.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta8.txt", skiprows= 2, usecols= 0, unpack= True)[0]
P9_theta = np.arctan(np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta8.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta8.txt", skiprows= 2, usecols= 2, unpack= True) )
P9_x = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta8.txt", skiprows= 2, usecols= 1, unpack= True)
P9_y = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta8.txt", skiprows= 2, usecols= 2, unpack= True) 
P9_unc = np.sqrt(abs((reading_unc * (P9_y / (P9_x ** 2 + P9_y**2)))**2  +  (reading_unc * (P9_x / (P9_x**2  + P9_y**2)))**2))



P10_times = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta9.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta9.txt", skiprows= 2, usecols= 0, unpack= True)[0]
P10_theta = np.arctan(np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta9.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta9.txt", skiprows= 2, usecols= 2, unpack= True) )
P10_x = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta9.txt", skiprows= 2, usecols= 1, unpack= True)
P10_y = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta9.txt", skiprows= 2, usecols= 2, unpack= True) 
P10_unc = np.sqrt(abs((reading_unc * (P10_y / (P10_x ** 2 + P10_y**2)))**2  +  (reading_unc * (P10_x / (P10_x**2  + P10_y**2)))**2))



P11_times = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta10.txt", skiprows= 2, usecols= 0, unpack= True) - np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta10.txt", skiprows= 2, usecols= 0, unpack= True)[0]
P11_theta = np.arctan(np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta10.txt", skiprows= 2, usecols= 1, unpack= True) / np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta10.txt", skiprows= 2, usecols= 2, unpack= True) )
P11_x = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta10.txt", skiprows= 2, usecols= 1, unpack= True)
P11_y = np.loadtxt("./Data 2/Varying Theta (L = 100mm)/theta10.txt", skiprows= 2, usecols= 2, unpack= True) 
P11_unc = np.sqrt(abs((reading_unc * (P11_y / (P11_x ** 2 + P11_y**2)))**2  +  (reading_unc * (P11_x / (P11_x**2  + P11_y**2)))**2))






############################################ FUNCTIONAL FORMS ############################################

thetas = np.array([P2_theta[0], P5_theta[0], P6_theta[0], \
                   P7_theta[0], P8_theta[0], P9_theta[0], P10_theta[0], P11_theta[0]])

def theta(t, A, T, tau, phi):
    return A * np.exp(-t/tau) * np.cos((2 * np.pi * t / T) + phi)

def linear(x, A, B):
    return A*x + B

def quad(x, A, B, C):
    return A*x**2 + B*x + C

def sqrt(x, A, B):
    return A * np.sqrt(x) + B


def constant(x, A):
    return 0 * x + A



############################################ CURVE FITTING DATA ############################################

interval_min = 0
interval_max = 100



popt_L1, pcov_L1 = curve_fit(theta, L1_times[interval_min:interval_max], L1_theta[interval_min:interval_max], sigma = L1_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_L2, pcov_L2 = curve_fit(theta, L2_times[interval_min:interval_max], L2_theta[interval_min:interval_max], sigma = L2_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_L3, pcov_L3 = curve_fit(theta, L3_times[interval_min:interval_max], L3_theta[interval_min:interval_max], sigma = L3_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_L4, pcov_L4 = curve_fit(theta, L4_times[interval_min:interval_max], L4_theta[interval_min:interval_max], sigma = L4_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_L5, pcov_L5 = curve_fit(theta, L5_times[interval_min:interval_max], L5_theta[interval_min:interval_max], sigma = L5_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_L6, pcov_L6 = curve_fit(theta, L6_times[interval_min:interval_max], L6_theta[interval_min:interval_max], sigma = L6_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_L7, pcov_L7 = curve_fit(theta, L7_times[interval_min:interval_max], L7_theta[interval_min:interval_max], sigma = L7_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)




popt_M1, pcov_M1 = curve_fit(theta, M1_times[interval_min:interval_max], M1_theta[interval_min:interval_max], sigma = M1_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_M2, pcov_M2 = curve_fit(theta, M2_times[interval_min:interval_max], M2_theta[interval_min:interval_max], sigma = M2_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_M3, pcov_M3 = curve_fit(theta, M3_times[interval_min:interval_max], M3_theta[interval_min:interval_max], sigma = M3_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_M4, pcov_M4 = curve_fit(theta, M4_times[interval_min:interval_max], M4_theta[interval_min:interval_max], sigma = M4_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_M5, pcov_M5 = curve_fit(theta, M5_times[interval_min:interval_max], M5_theta[interval_min:interval_max], sigma = M5_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)

popt_M6, pcov_M6 = curve_fit(theta, M6_times[interval_min:interval_max], M6_theta[interval_min:interval_max], sigma = M6_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_M7, pcov_M7 = curve_fit(theta, M7_times[interval_min:interval_max], M7_theta[interval_min:interval_max], sigma = M7_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_M8, pcov_M8 = curve_fit(theta, M8_times[interval_min:interval_max], M8_theta[interval_min:interval_max], sigma = M8_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_M9, pcov_M9 = curve_fit(theta, M9_times[interval_min:interval_max], M9_theta[interval_min:interval_max], sigma = M9_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)




popt_P1, pcov_P1 = curve_fit(theta, P1_times[interval_min:interval_max], P1_theta[interval_min:interval_max], sigma = P1_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_P2, pcov_P2 = curve_fit(theta, P2_times[interval_min:interval_max], P2_theta[interval_min:interval_max], sigma = P2_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_P3, pcov_P3 = curve_fit(theta, P3_times[interval_min:interval_max], P3_theta[interval_min:interval_max], sigma = P3_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_P4, pcov_P4 = curve_fit(theta, P4_times[interval_min:interval_max], P4_theta[interval_min:interval_max], sigma = P4_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_P5, pcov_P5 = curve_fit(theta, P5_times[interval_min:interval_max], P5_theta[interval_min:interval_max], sigma = P5_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_P6, pcov_P6 = curve_fit(theta, P6_times[interval_min:interval_max], P6_theta[interval_min:interval_max], sigma = P6_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)

popt_P7, pcov_P7 = curve_fit(theta, P7_times[interval_min:interval_max], P7_theta[interval_min:interval_max], sigma = P7_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_P8, pcov_P8 = curve_fit(theta, P8_times[interval_min:interval_max], P8_theta[interval_min:interval_max], sigma = P8_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_P9, pcov_P9 = curve_fit(theta, P9_times[interval_min:interval_max], P9_theta[interval_min:interval_max], sigma = P9_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_P10, pcov_P10 = curve_fit(theta, P10_times[interval_min:interval_max], P10_theta[interval_min:interval_max], sigma = P10_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)
popt_P11, pcov_P11 = curve_fit(theta, P11_times[interval_min:interval_max], P11_theta[interval_min:interval_max], sigma = P11_unc[interval_min:interval_max], p0=(1, 1, 1, 1), absolute_sigma=True, maxfev=5000)


############################################ VARYING LENGTH

A_L = [popt_L1[0], popt_L2[0], popt_L3[0], popt_L4[0], popt_L5[0]]
A_L_unc = [np.sqrt((np.diag(pcov_L1))[0]), np.sqrt((np.diag(pcov_L2))[0]), np.sqrt((np.diag(pcov_L3))[0]),\
            np.sqrt((np.diag(pcov_L4))[0]), np.sqrt((np.diag(pcov_L5))[0])]

T_L = [popt_L1[1], popt_L2[1], popt_L3[1], popt_L4[1], popt_L5[1], popt_L6[1], popt_L7[1]]
T_L_unc = [np.sqrt((np.diag(pcov_L1))[1]), np.sqrt((np.diag(pcov_L2))[1]), np.sqrt((np.diag(pcov_L3))[1]),\
            np.sqrt((np.diag(pcov_L4))[1]), np.sqrt((np.diag(pcov_L5))[1]), np.sqrt((np.diag(pcov_L6))[1]), np.sqrt((np.diag(pcov_L7))[1])]

tau_L = [popt_L1[2], popt_L2[2], popt_L3[2], popt_L4[2], popt_L5[2]]
tau_L_unc = [np.sqrt((np.diag(pcov_L1))[2]), np.sqrt((np.diag(pcov_L2))[2]), np.sqrt((np.diag(pcov_L3))[2]),\
            np.sqrt((np.diag(pcov_L4))[2]), np.sqrt((np.diag(pcov_L5))[2])]

phi_L = [popt_L1[3], popt_L2[3], popt_L3[3], popt_L4[3], popt_L5[3]]
phi_L_unc = [np.sqrt((np.diag(pcov_L1))[3]), np.sqrt((np.diag(pcov_L2))[3]), np.sqrt((np.diag(pcov_L3))[3]),\
             np.sqrt((np.diag(pcov_L4))[3]), np.sqrt((np.diag(pcov_L5))[3])]




############################################ VARYING MASS


A_M = [popt_M1[0], popt_M2[0], popt_M3[0], popt_M4[0], popt_M5[0]]
A_M_unc = [np.sqrt((np.diag(pcov_M1))[0]), np.sqrt((np.diag(pcov_M2))[0]), np.sqrt((np.diag(pcov_M3))[0]),\
            np.sqrt((np.diag(pcov_M4))[0]), np.sqrt((np.diag(pcov_M5))[0])]


T_M = [popt_M1[1], popt_M2[1], popt_M3[1], popt_M4[1], popt_M5[1], popt_M6[1], popt_M7[1], popt_M8[1], popt_M9[1]]
T_M_unc = [np.sqrt((np.diag(pcov_M1))[1]), np.sqrt((np.diag(pcov_M2))[1]), np.sqrt((np.diag(pcov_M3))[1]),\
            np.sqrt((np.diag(pcov_M4))[1]), np.sqrt((np.diag(pcov_M5))[1]), np.sqrt((np.diag(pcov_M6))[1]), np.sqrt((np.diag(pcov_M7))[1]),\
                np.sqrt((np.diag(pcov_M8))[1]), np.sqrt((np.diag(pcov_M9))[1])]

tau_M = [popt_M1[2], popt_M2[2], popt_M3[2], popt_M4[2], popt_M5[2], popt_M6[2], popt_M7[2], popt_M8[2], popt_M9[2]]
tau_M_unc = [np.sqrt((np.diag(pcov_M1))[2]), np.sqrt((np.diag(pcov_M2))[2]), np.sqrt((np.diag(pcov_M3))[2]),\
            np.sqrt((np.diag(pcov_M4))[2]), np.sqrt((np.diag(pcov_M5))[2])]

phi_M = [popt_M1[3], popt_M2[3], popt_M3[3], popt_M4[3], popt_M5[3]]
phi_M_unc = [np.sqrt((np.diag(pcov_M1))[3]), np.sqrt((np.diag(pcov_M2))[3]), np.sqrt((np.diag(pcov_M3))[3]),\
            np.sqrt((np.diag(pcov_M4))[3]), np.sqrt((np.diag(pcov_M5))[3])]



############################################ VARYING THETA


A_P = [popt_P1[0], popt_P2[0], popt_P3[0], popt_P4[0], popt_P5[0], popt_P6[0], popt_P7[0], popt_P8[0], popt_P9[0], popt_P10[0], popt_P11[0]]
A_P_unc = [np.sqrt((np.diag(pcov_P1))[0]), np.sqrt((np.diag(pcov_P2))[0]), np.sqrt((np.diag(pcov_P3))[0]),\
            np.sqrt((np.diag(pcov_P4))[0]), np.sqrt((np.diag(pcov_P5))[0]), np.sqrt((np.diag(pcov_P6))[0]), np.sqrt((np.diag(pcov_P7))[0]), \
                np.sqrt((np.diag(pcov_P8))[0]), np.sqrt((np.diag(pcov_P9))[0]), np.sqrt((np.diag(pcov_P10))[0]), np.sqrt((np.diag(pcov_P11))[0])]


T_P = [popt_P2[1], popt_P5[1], popt_P6[1], popt_P7[1], popt_P8[1], popt_P9[1], popt_P10[1], popt_P11[1]]
T_P_unc = [np.sqrt((np.diag(pcov_P2))[1]), \
            np.sqrt((np.diag(pcov_P5))[1]), np.sqrt((np.diag(pcov_P6))[1]), np.sqrt((np.diag(pcov_P7))[1]), \
                 np.sqrt((np.diag(pcov_P8))[1])/3, np.sqrt((np.diag(pcov_P9))[1]), np.sqrt((np.diag(pcov_P10))[1]), np.sqrt((np.diag(pcov_P11))[1]) ]

tau_P = [popt_P1[2], popt_M2[2], popt_P3[2], popt_P4[2], popt_P5[2], popt_P6[2], popt_P7[2], popt_P8[2], popt_P9[2], popt_P10[2], popt_P11[2]]
tau_P_unc = [np.sqrt((np.diag(pcov_P1))[2]), np.sqrt((np.diag(pcov_P2))[2]), np.sqrt((np.diag(pcov_P3))[2]),\
            np.sqrt((np.diag(pcov_P4))[2]), np.sqrt((np.diag(pcov_P5))[2]), np.sqrt((np.diag(pcov_P6))[2]), np.sqrt((np.diag(pcov_P7))[2]), np.sqrt((np.diag(pcov_P8))[2]),\
                np.sqrt((np.diag(pcov_P9))[2]), np.sqrt((np.diag(pcov_P10))[2]), np.sqrt((np.diag(pcov_P11))[2])]

phi_P = [popt_P1[3], popt_P2[3], popt_P3[3], popt_P4[3], popt_P5[3], popt_P6[3]]
phi_P_unc = [np.sqrt((np.diag(pcov_P1))[3]), np.sqrt((np.diag(pcov_P2))[3]), np.sqrt((np.diag(pcov_P3))[3]),\
            np.sqrt((np.diag(pcov_P4))[3]), np.sqrt((np.diag(pcov_P5))[3]), np.sqrt((np.diag(pcov_P6))[3])]



############################################ T VERSUS THETA_0

'''

xvals_PvT = np.linspace(min(thetas - 0.01), max(thetas + 0.01), 100)

avg = np.average(T_P)

plt.plot(xvals_PvT, constant(xvals_PvT, avg), c = 'blue', label = 'Average Period')
plt.errorbar(thetas, T_P, yerr = T_P_unc, capsize = 2, c = 'red', fmt='o', alpha = myalpha, label = 'Collected Data')
plt.xlabel(r'Angular Amplitude $\theta_0$ [rads]')
plt.ylabel(r'Period $T$ [s]')
plt.legend(loc = 'upper left')
plt.tight_layout()
if (save): plt.savefig('amplitude_vs_period.png',dpi=mydpi)
plt.show()

'''

############################################ T VERSUS LENGTH



'''
xvals_LvT = np.linspace(min(L_vals - 0.01), max(L_vals + 0.01), 100)

popt_LvT, pcov_LvT = curve_fit(sqrt, L_vals, T_L, sigma = T_L_unc, p0=(1, 0), absolute_sigma=True, maxfev=5000)

chisquared_LvT = np.sum(((T_L - sqrt(L_vals, *popt_LvT))/T_L_unc )**2) / (len(L_vals) - len(popt_LvT))

print(popt_LvT, np.sqrt(np.diag(pcov_LvT)))

plt.plot(xvals_LvT, sqrt(xvals_LvT, *popt_LvT), c = 'red', label = 'Functional Fit')
plt.errorbar(L_vals, T_L, xerr = L_vals_unc, yerr = T_L_unc, capsize = 2, c = 'blue', fmt='o', alpha = myalpha, label = 'Collected Data')
plt.text(0.07, 0.93, r'$\chi^2/DOF = $ %3.2f'%(chisquared_LvT), fontsize=fontsize, c = 'red', alpha = myalpha)
plt.text(0.07, 0.89, r'$T(\ell) = $ %3.2f'%(popt_LvT[0]) + r'$\sqrt{\ell}$', fontsize=fontsize, c = 'red', alpha = myalpha)
plt.xlabel(r'Length $\ell$ [m]')
plt.ylabel(r'Period $T$ [s]')
plt.legend(loc = 'upper left')
plt.tight_layout()
if (save): plt.savefig('length_vs_period.png',dpi=mydpi)
plt.show()


plt.plot(xvals_LvT, sqrt(xvals_LvT, 0, 0), c = 'red', label = 'Functional Fit')
plt.errorbar(L_vals, T_L - sqrt(L_vals, *popt_LvT), xerr = L_vals_unc, yerr = T_L_unc, capsize = 2, c = 'blue', fmt='o', alpha = myalpha, label = 'Data Residuals')
plt.xlabel(r'Length $\ell$ [m]')
plt.ylabel(r'Period Residual [s]')
plt.legend(loc = 'lower left')
plt.tight_layout()
if (save): plt.savefig('length_vs_period_residual.png',dpi=mydpi)
plt.show()

'''


############################################ T VERSUS MASS


'''
xvals_MvT = np.linspace(min(masses - 0.01), max(masses + 0.01), 100)

avg = np.average(T_M)

plt.plot(xvals_MvT, constant(xvals_MvT, avg), c = 'blue', label = 'Average Period')
plt.errorbar(masses, T_M, yerr = T_M_unc, capsize = 2, c = 'red', fmt='o', alpha = myalpha, label = 'Collected Data')
plt.xlabel(r'Mass [kg]')
plt.ylabel(r'Period $T$ [s]')
plt.legend(loc = 'lower right')
plt.tight_layout()
if (save): plt.savefig('mass_vs_period.png',dpi=mydpi)
plt.show()

'''


############################################ PLOTTING



'''

plt.plot(M1_times[interval_min:interval_max], theta(M1_times[interval_min:interval_max], *popt_M1), c = 'red', label = 'Functional Fit')
plt.scatter(M1_times[interval_min:interval_max], M1_theta[interval_min:interval_max], s = 2, c = 'green', label = 'Collected Data Residual, L1')
#plt.text(0.01, 145, r'$\chi_1^2/DOF = $ %3.2f/%i'%(chisquared1,dof1), fontsize=fontsize, c = 'red')
plt.xlabel('Time (s)')
plt.ylabel('Angle (rads)')
plt.legend(loc = 'upper right')
plt.tight_layout()
#if (save): plt.savefig('L1.png',dpi=mydpi)
plt.show()


plt.plot(M2_times[interval_min:interval_max], theta(M2_times[interval_min:interval_max], *popt_M2), c = 'red', label = 'Functional Fit')
plt.scatter(M2_times[interval_min:interval_max], M2_theta[interval_min:interval_max], s = 2, c = 'green', label = 'Collected Data Residual, L1')
#plt.text(0.01, 145, r'$\chi_1^2/DOF = $ %3.2f/%i'%(chisquared1,dof1), fontsize=fontsize, c = 'red')
plt.xlabel('Time (s)')
plt.ylabel('Angle (rads)')
plt.legend(loc = 'upper right')
plt.tight_layout()
#if (save): plt.savefig('L1.png',dpi=mydpi)
plt.show()

'''



