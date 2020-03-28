import streamlit as st
import numpy as np
import pandas as pd
import math

#Function to solve Kepler's Equation for the Eccentric Anomaly. 
#Parameters: M_e defined as mean anomaly, e defined as eccentricity
def kepler_ellip(M_e, e):
    
    #Iteration Tolerance
    tol = 1E-8

    #Initial Guess at Eccentric Anomaly
    if M_e < math.pi:
        E = M_e + (e / 2)
    if M_e > math.pi:
        E = M_e - (e / 2)
    
    #Initial Conditions
    f = E - e*math.sin(E) - M_e
    f_prime = 1 - e*math.cos(E)
    ratio = f / f_prime

    #Numerical interation for ratio compared to tolerance 
    while abs(ratio) > tol:
        f = E - e*math.sin(E) - M_e
        f_prime = 1 - e*math.cos(E)
        ratio = f / f_prime

        if abs(ratio) > tol:
            E = E - ratio
        if abs(ratio) < tol:
            break
    
    #Returning E as the Eccentric Anomaly
    return E

def kepler_hyper(M_h, e):
    #Iteration Tolerance
    tol = 1E-8

    #Initial Guess at Hyperbolic Eccentric Anomaly
    F_h = M_h
    
    #Initial Conditions
    f = e*math.sinh(F_h) - F_h - M_h
    f_prime = e* math.cosh(F_h) - 1
    ratio = f / f_prime

    #Numerical interation for ratio compared to tolerance 
    while abs(ratio) > tol:
        f = e*math.sinh(F_h) - F_h - M_h
        f_prime = e* math.cosh(F_h) - 1
        ratio = f / f_prime

        if abs(ratio) > tol:
            F_h = F_h - ratio
        if abs(ratio) < tol:
            break
    
    #Returning F_h as Hyperbolic Eccentric Anomaly
    return F_h

def orbit_elements(r, v):
    #Block of code for encapsulation
    #Mostly vector calculations 
    #################################################################
    #################################################################
    #Cross Product for two vectors a x b
    def cross(a, b):
        c = [a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]]

        return c

    #Dot Product for two vectors a . b 
    def dot(a, b):
        return sum(x*y for x,y in zip(a,b))

    #Magnitude of given vector a
    def mag(a):
        magnitude = math.sqrt(sum(v**2 for v in a))
        
        return magnitude

    #Multiplying a vector (b) by coefficient (a)
    def coeff(a, b):
        newVector = [0,0,0]
        newVector[0] = a * b[0]
        newVector[1] = a * b[1]
        newVector[2] = a * b[2]

        return newVector 

    #Subtracting two vectors (a - b)
    def subtract_vectors(a, b):
        newVector = [0,0,0]
        newVector[0] = a[0] - b[0]
        newVector[1] = a[1] - b[1]
        newVector[2] = a[2] - b[2]

        return newVector
    #################################################################
    #################################################################
    #End of encapsulation
    
    #Earth Gravitational Constant
    mu = 3.986E5

    #Magnitudes of state vectors
    r_mag = mag(r)
    v_mag = mag(v)

    #Radial Component of velocity
    v_r = (dot(r,v)) / r_mag
    
    #Specific Angular Momentum Vector & Magnitude
    h = cross(r, v)
    h_mag = mag(h)

    #Inclination (degrees)
    i = math.acos(h[2] / h_mag) 

    #Node Line and Node Magnitude
    k_hat = [0, 0, 1]
    node = cross(k_hat, h)
    node_mag = mag(node)

    #Right Ascension of the Ascending Node
    omega = math.acos(node[0] / node_mag) 
    if node[1] < 0:
        omega = 360 - omega

    #Eccentricity Vector and Magnitude
    #Equation:  e = (1/mu)(V x h - mu * (r / r_mag))
    x = cross(v,h)
    c = mu * (1/r_mag)
    newr = coeff(c, r)
    e = subtract_vectors(x, newr)
    inv = 1 / mu

    e = coeff(inv, e)
    e_mag = mag(e)
    
    #Argument of Perigee
    node_unit = coeff(1/node_mag, node)
    eccentricity_unit = coeff(1/e_mag, e)
    w = math.acos(dot(node_unit, eccentricity_unit)) 
    if e[2] < 0:
        w = 2*math.pi - w
    
    #True Anomaly
    eccentricity_unit = coeff(1/e_mag, e)
    position_unit = coeff(1/r_mag, r)
    theta = math.acos(dot(eccentricity_unit, position_unit)) 
    if v_r < 0:
        theta = 2*math.pi - theta
    #List elements: 
    #[Eccentricity, Specific Angular Momentum, Inclination, Right Ascension of the Ascending Node, Arument of Perigee, True Anomaly]
    return [e_mag, h_mag, i * (180 / math.pi), omega * (180 / math.pi), w * (180 / math.pi), theta * (180 / math.pi)]

if __name__ == "__main__":
    #Heading title of website
    st.title('Basic Orbital Dynamics Scripts')

    #List of possible scripts to choose from in primary drop down display
    scripts = ['Select Script', 'Orbital Elements From State Vectors', "Kepler's Equation (Elliptical)", "Kepler's Equation (Hyperbolic)"]
    page = st.selectbox('', options=scripts)

    #Director page for Kepler's Equation (Elliptical) Script
    if page == "Kepler's Equation (Elliptical)":
        #User sidebar numerical input
        mean_anomaly = st.sidebar.number_input('Mean Anomaly')
        eccentricity = st.sidebar.number_input('Eccentricity')
        
        #Function Call for Solving Equation: Returned as a float (hopefully)
        eccentric_anomaly = kepler_ellip(mean_anomaly, eccentricity)

        #Output of two given parameters (1/2) and desired output (3)
        output1 = st.write('Eccentricity:', eccentricity)
        output2 = st.write('Mean Anomaly:', mean_anomaly)
        output3 = st.write('Eccentric Anomaly:', eccentric_anomaly)

    if page == "Kepler's Equation (Hyperbolic)":
        #User sidebar numerical input
        hyper_mean_anomaly = st.sidebar.number_input('Hyperbolic Mean Anomaly')
        eccentricity = st.sidebar.number_input('Eccentricity')
        
        #Function Call for Solving Equation: Returned as a float (hopefully)
        hyper_eccentric_anomaly = kepler_hyper(hyper_mean_anomaly, eccentricity)

        #Output of two given parameters (1/2) and desired output (3)
        output1 = st.write('Eccentricity:', eccentricity)
        output2 = st.write('Hyperbolic Mean Anomaly:', hyper_mean_anomaly)
        output3 = st.write('Eccentric Anomaly:', hyper_eccentric_anomaly)
    
    if page == 'Orbital Elements From State Vectors':
        
        #Position Vector 
        velocity_title = st.sidebar.markdown('**Position Vector in Perifocal Frame**')
        r_x = st.sidebar.number_input('Position Vector X-Component (km)')
        r_y = st.sidebar.number_input('Position Vector Y-Component (km)')
        r_z = st.sidebar.number_input('Position Vector Z-Component (km)')
        
        #Compile components into list(vector)
        r_vector = [r_x, r_y, r_z]

        #Formatter
        empty = st.sidebar.markdown('')

        #Velocity Vector
        velocity_title = st.sidebar.markdown('**Velocity Vector in Perifocal Frame**')
        v_x = st.sidebar.number_input('Velocity Vector X-Component (km/s)')
        v_y = st.sidebar.number_input('Velocity Vector Y-Component (km/s)')
        v_z = st.sidebar.number_input('Velocity Vector Z-Component (km/s)')

        #Compile components into list(vector)
        v_vector = [v_x, v_y, v_z]

        #Function calls to grab orbital elements also to avoid dividing by zero
        if r_x != 0.0 and r_y != 0.0 and r_z != 0.0 and v_x != 0.0 and v_y != 0.0 and v_z != 0.0:
            orbital_parameters = orbit_elements(r_vector, v_vector)
            #Output of all orbital elements
            #[Eccentricity, Specific Angular Momentum, Inclination, Right Ascension of the Ascending Node, Arument of Perigee, True Anomaly]
            e = orbital_parameters[0]
            st.write('Eccentricity:', e)

            sam = orbital_parameters[1]
            st.write('Specific Angular Momentum:', sam)

            inclination = orbital_parameters[2]
            st.write('Inclination:', inclination)

            right_ascension_of_ascending_node = orbital_parameters[3]
            st.write('Right Ascension of the Ascending Node:', right_ascension_of_ascending_node)

            arg_perigee = orbital_parameters[4]
            st.write('Arument of Perigee:', arg_perigee)

            true_anomaly = orbital_parameters[5]
            st.write('True Anomaly:', true_anomaly)
            
        if r_x == 0.0 and r_y == 0.0 and r_z == 0.0 and v_x == 0.0 and v_y == 0.0 and v_z == 0.0:
            st.write('Please enter state vectors.')

        
       




        
