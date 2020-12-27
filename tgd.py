def toggle_derivative(t,y,p):

    # Y[0] is Lacl mRNA concentration
    # Y[1] is TetR mRNA concentration
    # Y[2] is LacI protein concentration
    # Y[3] is TetR protein concentration

    ####################################################
    ##### Deterministic Model of the Toggle Switch #####
    ####################################################

    dydt= [0, 0, 0, 0]

    # LacI mRNA
    dydt[0] = -p[8]*y[0] + p[2]/(1 + (y[3] / p[5]) ** p[7]) + p[0]
    # TetR mRNA
    dydt[1] = - p[9] * y[1] + p[3] / (1 + (y[2] / p[4]) ** p[6]) + p[1]
    # LacI
    dydt[2] = p[10]*y[0] - p[12]*y[2] 
    # TetR
    dydt[3] = p[11]*y[1] - p[13]*y[3]
    
    return dydt
