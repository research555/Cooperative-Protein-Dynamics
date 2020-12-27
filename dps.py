import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

# we import our own functions from our files
# place them all in the same folder for them to be found
from tgd import toggle_derivative
from intersections import intersection

def draw_phase_space(p):

    #######################################
    ##### draw LacI/TetR vector field #####
    #######################################
    
    # we want a 35x35 vector field ranging with values ranging from 0 to 2000, this can be changed
    LacI_vector= np.linspace(0, 2000, 35)
    TetR_vector= np.linspace(0, 2000, 35)
    
    # we will need an array of 35 values for each mRNA
    mRNAt = np.zeros((35))
    mRNAl = np.zeros((35))
    
    # a 35x35x2 matrix to store the lacI/tetR values
    # a 35x35x2 matrix to store the derivative values
    grid = np.zeros((35, 35, 2))
    deriv = np.zeros((35, 35, 2))
      
    # Now, assuming that mRNAs are at steady state
    # we compute the tetR mRNA values for each of the 35 lacI values
    # and the lacI mRNA values for each of the 35 tetR values
    for i in range(len(LacI_vector)):
        lacIprot = LacI_vector[i]
        mRNAt[i] = (p[1] + p[3]/(1+(lacIprot/p[4])**p[6]) + p[1])/p[9] #Artur's code
    for j in range(len(TetR_vector)):
        tetRprot = TetR_vector[j]
        mRNAl[j] = (p[0] + p[2]/(1+(TetR_vector[j]/p[5])**p[7]))/p[8]
    
    # we need 35x35 arrows : for each arrow we need the
    # current values (grid), and the derivative values (deriv)
    for i in range(len(LacI_vector)):
        for j in range(len(TetR_vector)):
            grid[i][j] = [LacI_vector[i],TetR_vector[j]]
            dydt = toggle_derivative(0, [mRNAl[j], mRNAt[i], LacI_vector[i], TetR_vector[j]], p)
            deriv[i][j] = [dydt[2], dydt[3]]

    plt.quiver(grid[:,:,0],grid[:,:,1],deriv[:,:,0],deriv[:,:,1])

    return ()
    
    ###############################################################################
    ##### draw LacI/TetR nullclines with steady state approximation for mRNAs #####
    ###############################################################################
    
    # we will take 300 points between 0 and 2000
    LacI_vector = np.linspace(0,2000,300)
    TetR_vector = np.linspace(0,2000,300)
    
    # we also need vectors to store nullclines:
    nLacI_vector = np.zeros(300)
    nTetR_vector = np.zeros(300)
    
    # for every LacI value we compute a TetR nullcline value
    # and we do the same for TetR
    for i in range(len(LacI_vector)):
        mRNAt = (p[1] + p[3]/(1+(LacI_vector[i]/p[4])**p[6]))/p[9]
        nTetR_vector[i]= p[11]*mRNAt/p[13]

    for j in range(len(TetR_vector)):
        mRNAl = (p[0] + p[2]/(1+(TetR_vector[j]/p[5])**p[7]))/p[8]                 
        nLacI_vector[j]= p[10]*mRNAl/p[12]

    
    # then we use our "intersection" function to compute nullcline intersections
    # thus giving us the steady states
    plt.plot(nLacI_vector, TetR_vector, 'g', linewidth=4, label='LacI nullcline')
    plt.plot(LacI_vector, nTetR_vector, 'c', linewidth=4, label='TetR nullcline')
    plt.legend()