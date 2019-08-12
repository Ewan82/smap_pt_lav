import numpy as np
#from SALib.sample import saltelli
#from SALib.analyze import sobol


def head(suction, g=9.80665, rho_water=1000.):
    """
    Converts suction (KPa) to head (m)
    :param suction: value of suction in KPa (float)
    :param g: gravity (float)
    :param rho_water: density of water (float)
    :return: head in m (float)
    """
    return (-1000.*suction)/(g*rho_water)


def cosby(f_clay, f_sand, f_silt, a=3.10, b=15.70, c=0.3, d=0.505, e=0.037, f=0.142, g=2.17, h=0.63, i=1.58, j=5.55, k=0.64,
          l=1.26):
    """
    Function calculating JULES soil hydrological parameters for Brooks & Corey model from soil texture fraction
    :param f_clay: fraction of clay (float)
    :param f_sand: fraction of sand (float)
    :param f_silt: fraction of silt (float)
    :param a: defined fn parameter (float)
    :param b: defined fn parameter (float)
    :param c: defined fn parameter (float)
    :param d: defined fn parameter (float)
    :param e: defined fn parameter (float)
    :param f: defined fn parameter (float)
    :param g: defined fn parameter (float)
    :param h: defined fn parameter (float)
    :param i: defined fn parameter (float)
    :param j: defined fn parameter (float)
    :param k: defined fn parameter (float)
    :param l: defined fn parameter (float)
    :return: b exponent, theta_sat, saturated soil moisture pressure (m),
             saturate soil hydraulic conductivity (kg/m2/s), theta_wilt, theta_crit
    """
    # Hydrological parameters
    psi_c = -33.
    psi_w = -1500.
    b = a + b*f_clay - c*f_sand
    v_sat = d - e*f_clay - f*f_sand

    psi_s = 0.01*np.exp(g - h*f_clay - i*f_sand)
    k_s = np.exp(-j - k*f_clay + l*f_sand)
    head_c = head(psi_c)
    head_w = head(psi_w)
    v_wilt = v_sat*(psi_s/head_w)**(1/b)
    v_crit = v_sat*(psi_s/head_c)**(1/b)
    # Heat capacity, combines heat capapcities of clay, sand and silt linearly
    c_s = (1 - v_sat)*(f_clay*2.373e6 + f_sand*2.133e6 + f_silt*2.133e6)
    # Thermal conductivity
    lam = (0.025**(v_sat))*(1.16**(f_clay*(1-v_sat)))*(1.57**(f_sand*(1-v_sat)))*(1.57**(f_silt*(1-v_sat)))
    return np.array([b, v_sat, psi_s, k_s, v_wilt, v_crit, c_s, lam])


def cosby_wrap(X):
    return cosby(0.2, 0.3, X[:,0], X[:,1], X[:,2], X[:,3], X[:,4], X[:,5], X[:,6], X[:,7], X[:,8], X[:,9], X[:,10],
                 X[:,11])


def sensitivity():
    param_means = [3.10, 15.70, 0.3, 0.505, 0.037, 0.142, 2.17, 0.63, 1.58, 5.55,
                   0.64, 1.26]
    problem = {'num_vars': 12,
               'names': ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l'],
               'bounds': [[i-0.2*i, i+0.2*i] for i in param_means]}
    param_values = saltelli.sample(problem, 1000000)

    Y = cosby_wrap(param_values)
    return Y