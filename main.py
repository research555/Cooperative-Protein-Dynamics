import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from dps import draw_phase_space
from tgd import toggle_derivative


#############################
##### System Parameters #####
#############################

# production rates (basal/induced) for LacI/TetR mRNA
kappa_L_m0= 0.0082
kappa_T_m0= 0.0149
kappa_L_m= 1
kappa_T_m= 0.3865*3 #making a switch

# Hill function parameters for LacI/TetR regulation
theta_L= 600
theta_T= 500
eta_L= 4 
eta_T= eta_L

# dilution rates for LacI/TetR mRNA
gamma_L_m = 0.04
gamma_T_m = gamma_L_m

# production rates for LacI/TetR proteins from mRNA
kappa_L_p = 0.1
kappa_T_p = kappa_L_p

# dilution rates for LacI/TetR proteins
gamma_L_p = 0.002
gamma_T_p = gamma_L_p

# saving those in a list of length 14
p = np.zeros((14))
p[0]= kappa_L_m0
p[1]= kappa_T_m0
p[2]= kappa_L_m
p[3]= kappa_T_m
p[4]= theta_L
p[5]= theta_T
p[6]= eta_L
p[7]= eta_T
p[8]= gamma_L_m
p[9]= gamma_T_m
p[10]= kappa_L_p
p[11]= kappa_T_p
p[12]= gamma_L_p
p[13]= gamma_T_p

################################
##### Numerical simulation #####
################################

# Y[0] is Lacl mRNA concentration
# Y[1] is TetR mRNA concentration
# Y[2] is LacI protein concentration
# Y[3] is TetR protein concentration
Y1 = np.array([100, 100, 2000, 500])  # assume no mRNA or protein is present
Y0 = np.array([100, 100, 500, 2000])  # assume no mRNA or protein is present

Tend = 6000

# this is our ODE solver. Many options are available
# the one that interest us is "args", that allow to pass
# arguments to functions. We will use it to ensure
# toggle_derivative has our "p" list with all parameters.
#    NOTE : args expects a "tuple"

sol = solve_ivp(toggle_derivative, (0, Tend), Y0, args=(p,))
sol1 = solve_ivp(toggle_derivative, (0, Tend), Y1, args=(p,))



#depending on who "wins", we chose a different color for the plot
if sol.y[2][-1] > 1000 and sol.y[3][-1] < 1000: # LacI wins
    trajectory_color= 'g' # green
    
elif sol.y[2][-1] < 1000 and sol.y[3][-1] > 1000 : # TetR wins
    trajectory_color= 'c' # cyan
    
else:
    trajectory_color= 'k' # black

# depending on who "wins", we chose a different color for the plot
if sol1.y[2][-1] > 1000 and sol1.y[3][-1] < 1000:  # LacI wins
    trajectory_color1 = 'g'  # green
elif sol1.y[2][-1] < 1000 and sol1.y[3][-1] > 1000:  # TetR wins
    trajectory_color1 = 'c'  # cya
else:
    trajectory_color1 = 'k'  # black


####################
##### Plotting #####
####################

# plotting the protein concentrations as a function of time
plt.figure()
plt.plot(sol.t, sol.y[0], "g")
plt.plot(sol.t, sol.y[1], "c")
plt.plot(sol.t, sol.y[2], "g", linewidth=3)
plt.plot(sol.t, sol.y[3], "c", linewidth=3)
plt.title('LacI & TetR evolution over time')
plt.xlabel('Timesteps')
plt.ylabel('Concentration')
plt.axis([0, Tend, 0, 2000])
plt.legend(('LacI mRNA', 'TetR mRNA', 'LacI', 'TetR'))


######################################################
##### plotting the trajectory in the state space #####
######################################################

plt.figure()
plt.xlabel('[LacI]')
plt.ylabel('[TetR]')
plt.axis ([0, 2000, 0, 2000])
plt.title('Superimposed Plot of State Space and Protein Trajectory For LacI & TetR')
plt.plot(sol.y[2], sol.y[3], '.', linewidth=2, color=trajectory_color)
plt.plot(sol1.y[2], sol1.y[3], '.', linewidth=2, color=trajectory_color1)

plt.show()



######################################
###### Plotting the vector field #####
######################################

plt.figure()
draw_phase_space(p)
plt.show()


######################################
##### Plotting many trajectories #####
######################################

# we want to plot 121 trajectories coming from 11x11 starting points
LacI_vector= np.linspace(0, 2000, 11)
TetR_vector= np.linspace(0, 2000, 11)

# now for each 121 starting points, we want to draw a curve
for i in range(len(LacI_vector)):
    for j in range(len(TetR_vector)):
        # Now, assuming that mRNAs are at steady state, we can compute them
        mRNAt= (p[1] + p[3]/(1+(LacI_vector[i]/p[4])**p[6]))/p[9]
        mRNAl= (p[0] + p[2]/(1+(TetR_vector[j]/p[5])**p[7]))/p[8]
        
        # We then have the 4 values we need for each starting point
        Y0 = [mRNAl,mRNAt,LacI_vector[i],TetR_vector[j]]
        
        # we compute the curve
        sol = solve_ivp(toggle_derivative,[0,Tend],Y0,args=(p,))
        
        # we chose the color based on who "wins"
        if sol.y[2][-1] > 1000 and sol.y[3][-1] <1000 : # LacI wins
            trajectory_color= 'g' # green
    
        elif sol.y[2][-1] < 1000 and sol.y[3][-1] > 1000 : # TetR wins
            trajectory_color= 'c' # blue; cyan actually
        
        else:
            trajectory_color= 'k' # black
         
        # we plot the trajectory
        plt.plot(sol.y[2], sol.y[3],'-',linewidth=2, color=trajectory_color)

plt.show()