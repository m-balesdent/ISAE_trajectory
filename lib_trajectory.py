# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 08:59:26 2023

@author: mbalesde
"""

import numpy as np
import time
import copy
import scipy
import constants as Cst
from modele_atmos import compute_atmos
import matplotlib.pyplot as plt

def flight_dynamics(t,x,Parameters,Rocket_model):
    """
    Function calculating the right member of the equation of motions
        
    :t: time 
    
    :x: the state vector (radius, relative velocity, flight path angle ,longitude,mass)
    
    :Parameters: dictionary of parameters for the simulation (current number of stage, simulation mode)
    
    :Rocket_model: dictionary composed of rocket parameters
    
    """
    # State variables
    r = x[0] 
    V = x[1]
    gamma = x[2]
    longi = x[3]
    m = x[4]

    #Parameters for simulation
    Number_stage  = Parameters['Number_stage']
    Integration = Parameters['Mode_simu']


    #Launcher characteristics for drag calculation
    
    #Compute current atmosphere
    alt = r-Cst.RT
    (Temp,Pa,rho,c) = compute_atmos(alt)
    Mach = V/c
    
    #Compute the current gravity term
    g_current = Cst.mu/((alt+Cst.RT)**2)


    #Control command to get the command associated to the vehicle state 
    if Number_stage == 1:
        #Pitch angle
        theta= control_first_stage(t,x,Parameters,Rocket_model)
        #Angle of attack
        alpha = theta - gamma
        
        # Thrust
        thrust= Rocket_model['First_stage']['Mass_flow_rate']* Cst.g0*Rocket_model['First_stage']['Isp_vac'] -\
                Pa * Rocket_model['First_stage']['Nozzle_exit_area'] #thrust
        thrust_vac = Rocket_model['First_stage']['Mass_flow_rate']* Cst.g0*Rocket_model['First_stage']['Isp_vac'] 
        #Mass flow rate
        Mass_flow_rate=Rocket_model['First_stage']['Mass_flow_rate'] #engines mass flow rate
        
        #Aerodynamic forces
        CX = Parameters['Interp_drag'](Mach)  #drag coefficient as a function of Mach number
                
        # Reference area
        S_ref = Rocket_model['First_stage']['Reference_area']

        
        
    else:
        #Pitch angle
        theta= control_second_stage(t,x,Parameters,Rocket_model)
        #Angle of attack
        alpha = theta - gamma
        
        # Thrust
        thrust= Rocket_model['Second_stage']['Mass_flow_rate']* Cst.g0*Rocket_model['Second_stage']['Isp_vac'] -\
                Pa * Rocket_model['Second_stage']['Nozzle_exit_area']
        
        thrust_vac = Rocket_model['Second_stage']['Mass_flow_rate']* Cst.g0*Rocket_model['Second_stage']['Isp_vac']
        
        #Mass flow rat
        Mass_flow_rate=Rocket_model['Second_stage']['Mass_flow_rate'] #engines mass flow rate
         
        #Aerodynamic forces (assumption of vaccuum for the second stage)
        CX = 0.
        
        # Reference area
        S_ref = Rocket_model['Second_stage']['Reference_area']   
   
 
    #Equations of motion    
    r_dot = V* np.sin(gamma) #gradient of radius
    V_dot = -0.5*rho*S_ref*CX*V**2./m - g_current*np.sin(gamma) + thrust*np.cos(theta-gamma)/m #gradient of velocity
    gamma_dot = (V/r-g_current/V)*np.cos(gamma) + thrust*np.sin(theta-gamma)/(m*V)  #gradient of flight path angle
    longi_dot = V*np.cos(gamma)/r #gradient of longitude
    m_dot = - Mass_flow_rate #gradient of vehicle mass

    dx = np.zeros([1,5])
    dx[0][0] = r_dot
    dx[0][1] = V_dot
    dx[0][2] = gamma_dot
    dx[0][3] = longi_dot
    dx[0][4] = m_dot
    
    
    if alt<=0. : # Crash
        dx = np.zeros([1,5])
        
        
    if Integration == 1.: # Integration mode
        return dx[0]
    
    else :  # Simulation mode
        
        Pdyn = 0.5*rho*V**2
        flux= Pdyn*V
        
        CA = CX*np.cos(np.abs(alpha))  #aerodynamical force
        NX = (thrust-Pdyn*CA*S_ref)/(m*Cst.mu/((alt+Cst.RT)**2))  #axial load factor

        return (r,
                V,
                np.rad2deg(gamma),
                np.rad2deg(longi),
                m,
                NX,
                Pdyn,
                flux,
                alt,
                np.rad2deg(alpha),
                np.rad2deg(theta),
                rho,
                Pa,
                CX,
                thrust,
                thrust_vac,
                Mass_flow_rate,
                g_current)
    
    
    

def control_first_stage(t,x,Parameters,Rocket_model):
    """
    Function calculating the control law (pitch angle) for the first stage
        
    :t: time 
    
    :x: the state vector (radius, relative velocity, flight path angle ,longitude,mass)
    
    :Parameters: dictionary of parameters for the simulation (current number of stage, simulation mode)
    
    :Rocket_model: dictionary composed of rocket parameters
    
    """ 
    
    
    gamma = x[2]
    
    # Control law parameters
    Vertical_phase_duration = Parameters['Control']['First_stage']['Vertical_phase_duration']
    Pitch_over_duration =  Parameters['Control']['First_stage']['Pitch_over_duration']
    Delta_theta_pitch_over = Parameters['Control']['First_stage']['Delta_theta_pitch_over']
    Duration_decay = Parameters['Control']['First_stage']['Pitch_over_exp_decay_duration']
    
    
    # Vertical phase
    if t<Vertical_phase_duration : 
       alpha = 0.
       theta = np.pi/2.
       
    #Pitch over first phase  
    elif ((t>=Vertical_phase_duration) and (t<Vertical_phase_duration + Pitch_over_duration)):    ###pitch over manoeuver

        theta = (gamma*180/np.pi- Delta_theta_pitch_over * (t - Vertical_phase_duration) / Pitch_over_duration)*np.pi/180.
        alpha = theta- gamma
        
    # Exponential decay    
    elif ((t>=Vertical_phase_duration + Pitch_over_duration) and (t<Vertical_phase_duration + Pitch_over_duration+ Duration_decay)):    ###return to angle of attack = 0.
        theta = (gamma*180/np.pi - Delta_theta_pitch_over * np.exp(-(t - (Vertical_phase_duration + Pitch_over_duration)) / Duration_decay))*np.pi/180.
        alpha = theta- gamma
        
        
    # Gravity turn 
    else :
        #gravity turn phase
        theta = gamma
        alpha = theta - gamma 
        
    return theta
          
       
       
def control_second_stage(t,x,Parameters,Rocket_model):
    """
    Function calculating the control law (pitch angle) for the second stage
        
    :t: time 
    
    :x: the state vector (radius, relative velocity, flight path angle ,longitude,mass)
    
    :Parameters: dictionary of parameters for the simulation (current number of stage, simulation mode)
    
    :Rocket_model: dictionary composed of rocket parameters
    
    """
    
    # Control law parameters
    
    tf2 = Parameters['Control']['tf']
    theta_i = Parameters['Control']['Bilinear']['Initial_theta']
    theta_f = Parameters['Control']['Bilinear']['Final_theta']
    
    a = 100.
    ksi = Parameters['Control']['Bilinear']['ksi']
    
    # Bi linear tangent law
    theta = np.arctan2((pow(a,ksi)*np.tan(theta_i)+(np.tan(theta_f)-pow(a,ksi)*np.tan(theta_i))*(t)/tf2),(pow(a,ksi)+(1.-pow(a,ksi))*(t)/tf2))   
    return theta
       


def load_aerodynamics():
    """
    Function defining the drag coefficient as a function of the Mach number
    
    """
    Table_Mach = np.array([0.,0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.55, 0.6 , 0.65, 0.7 , 0.85, 0.95,
           1.2 , 1.7 , 2.7 , 3.7 , 4.7 , 5.7 , 6.7 , 7.7,40. ])
    
    Table_CX = np.array([0.086672,0.086672, 0.086123, 0.085465, 0.084549, 0.08352 , 0.082853,
           0.082218, 0.082354, 0.1385  , 0.25362 , 0.26857 , 0.26894 ,
           0.24551 , 0.20926 , 0.18949 , 0.17709 , 0.16715 , 0.16015 ,
           0.1548,0.1548  ])
    
    return Table_Mach, Table_CX



def trajectory_integration_first_stage(design_var, Parameters,Rocket_model):
    """
    Function computing the trajectory of the first stage
        
    :design_var: design variables of the first stage (payload mass and control law) 
        
    :Parameters: dictionary of parameters for the simulation
    
    :Rocket_model: dictionary composed of rocket parameters
    
    """
        
    # Calculation of payload mass
    Payload_mass = design_var[0]*1e3
    
    # Load aerodynamics
    Table_Mach, Table_CX = load_aerodynamics()
    Interp_drag = scipy.interpolate.interp1d(Table_Mach,Table_CX)    
    
    # Initial values of state variables for the trajectory integration
    mass_t0 = Rocket_model['First_stage']['Dry_mass'] + Rocket_model['First_stage']['Propellant_mass']\
              + Rocket_model['Second_stage']['Dry_mass']+ Rocket_model['Second_stage']['Propellant_mass']\
              + Rocket_model['Second_stage']['Fairing_mass'] + Payload_mass
                             
                             
    r_t0 = Cst.RT+1.
    v_t0 = 1.
    gamma_t0 = np.deg2rad(90.)
    longi_t0 = np.deg2rad(-52.65)
    
    # Definition of parameters for trajectory integration
    Parameters_integration = copy.deepcopy(Parameters)
    Parameters_integration['Mode_simu'] =1. # Integration : 1., simulation : 0.
    Parameters_integration['Number_stage'] = 1
    Parameters_integration['Interp_drag'] = Interp_drag
    Parameters_integration['Control']['First_stage']['Delta_theta_pitch_over'] = design_var[3]

    # ----- Integration of first stage trajectory
    
    #Integration parameters 
    atol_integration   = 1e-8
    rtol_integration   = 1e-8
    integration_method = 'RK45' #Implicit method based on backward-differentiation formulas see Scipy for more details
    step               = 2.
    
    
    Flight_duration =  Rocket_model['First_stage']['Propellant_mass'] / Rocket_model['First_stage']['Mass_flow_rate']
    
    initial_state = np.array([r_t0,\
                              v_t0,\
                              gamma_t0,\
                              longi_t0,\
                              mass_t0]) #initial state of launcher
        
    span_integration = (0.,Flight_duration)
    
    fonction_ode_integration = lambda t,x :flight_dynamics(t,x,Parameters_integration,Rocket_model)  #launcher simulation equations of motion
    
    Parameters_simu = copy.deepcopy(Parameters_integration)              
    Parameters_simu['Mode_simu']=0.
    
    fonction_ode_simu = lambda t,x :flight_dynamics(t,x,Parameters_simu,Rocket_model)  #launcher simulation equations of motion
    
    # Integration of trajectory
    solution = scipy.integrate.solve_ivp(fonction_ode_integration,
                                   span_integration,
                                   initial_state,
                                   atol         = atol_integration,
                                   rtol         = rtol_integration,
                                   dense_output = True,
                                   method       = 'RK45')
    
    
    # Simulation of trajectory to get history of variables 
    time = np.arange(solution.t[0],solution.t[-1],step)
    time = np.append(time, solution.t[-1])

    flight_history                      = {}
    flight_history['time']              = time
    flight_history['r']                 = np.zeros(len(time))
    flight_history['velocity']          = np.zeros(len(time))
    flight_history['flight_path_angle'] = np.zeros(len(time))
    flight_history['longitude']         = np.zeros(len(time))
    flight_history['mass']              = np.zeros(len(time))
    flight_history['axial_load_factor'] = np.zeros(len(time))
    flight_history['qdyn']              = np.zeros(len(time))
    flight_history['flux']              = np.zeros(len(time))
    flight_history['altitude']          = np.zeros(len(time))
    flight_history['AoA']               = np.zeros(len(time))
    flight_history['pitch_angle']       = np.zeros(len(time))
    flight_history['rho']               = np.zeros(len(time))
    flight_history['Pa']               = np.zeros(len(time))
    flight_history['Drag_coeff']        = np.zeros(len(time))
    flight_history['thrust']            = np.zeros(len(time))
    flight_history['thrust_vac']            = np.zeros(len(time))
    flight_history['mass_flow_rate']    = np.zeros(len(time))
    flight_history['g']                 = np.zeros(len(time))
    
    for i in range(len(time)):
        (flight_history['r'] [i],
         flight_history['velocity'][i],
         flight_history['flight_path_angle'][i],
         flight_history['longitude'][i],
         flight_history['mass'] [i],
         flight_history['axial_load_factor'][i],
         flight_history['qdyn'] [i],
         flight_history['flux'] [i],
         flight_history['altitude'][i],
         flight_history['AoA'][i],
         flight_history['pitch_angle'] [i],
         flight_history['rho'][i],
         flight_history['Pa'][i],
         flight_history['Drag_coeff'][i],
         flight_history['thrust'][i],
         flight_history['thrust_vac'][i],
         flight_history['mass_flow_rate'][i],
         flight_history['g'][i]) = fonction_ode_simu( time[i],
                                                                    solution.sol(time[i]))
    
    
    return flight_history
    
    


def trajectory_integration_second_stage(design_var,Parameters,Rocket_model,flight_history_first_stage):
    """
    Function computing the trajectory of the second stage
        
    :design_var: design variables of the first stage (payload mass and control law) 
        
    :Parameters: dictionary of parameters for the simulation
    
    :Rocket_model: dictionary composed of rocket parameters
    
    """
    
    
    # Calculation of payload mass
    
    Payload_mass = design_var[0]*1e3


    # Initial values of state variables for the trajectory integration

    # Jetissoning of fairing at the stage separation
    mass_t0 = Rocket_model['Second_stage']['Dry_mass'] + Rocket_model['Second_stage']['Propellant_mass']\
              + Payload_mass

    r_t0 = flight_history_first_stage['r'][-1]
    v_t0 = flight_history_first_stage['velocity'][-1]
    gamma_t0 = np.deg2rad(flight_history_first_stage['flight_path_angle'][-1])
    longi_t0 = np.deg2rad(flight_history_first_stage['longitude'][-1])
    
    # Definition of parameters for trajectory integration    
    Parameters_integration = copy.deepcopy(Parameters)
    Parameters_integration['Control']['Bilinear'] = {}
    Parameters_integration['Control']['Bilinear']['Initial_theta'] = np.deg2rad(design_var[1])
    Parameters_integration['Control']['Bilinear']['Final_theta']   =  0.
    Parameters_integration['Control']['Bilinear']['ksi']           = design_var[2]
    Parameters_integration['Mode_simu']    = 1. # Integration : 1., simulation : 0.
    Parameters_integration['Number_stage'] = 2

    # ----- Integration of first stage trajectory
    #Integration parameters 
    atol_integration   = 1e-8
    rtol_integration   = 1e-8
    integration_method = 'RK45' #Implicit method based on backward-differentiation formulas see Scipy for more details
    step               = 2.
    
    
    Flight_duration =  Rocket_model['Second_stage']['Propellant_mass'] / Rocket_model['Second_stage']['Mass_flow_rate']
    Parameters_integration['Control']['tf'] = Flight_duration

    initial_state = np.array([r_t0,\
                              v_t0,\
                              gamma_t0,\
                              longi_t0,\
                              mass_t0]) #initial state of launcher
        
    span_integration = (0.,Flight_duration)
    
    fonction_ode_integration = lambda t,x :flight_dynamics(t,x,Parameters_integration,Rocket_model)  #launcher simulation equations of motion
    
    Parameters_simu = copy.deepcopy(Parameters_integration)              
    Parameters_simu['Mode_simu']=0.
    
    fonction_ode_simu = lambda t,x :flight_dynamics(t,x,Parameters_simu,Rocket_model)  #launcher simulation equations of motion
    
    # Integration of trajectory
    solution = scipy.integrate.solve_ivp(fonction_ode_integration,
                                   span_integration,
                                   initial_state,
                                   atol         = atol_integration,
                                   rtol         = rtol_integration,
                                   dense_output = True,
                                   method       = 'RK45')
    
    
    # Simulation of trajectory to get history of variables 
    time = np.arange(solution.t[0],solution.t[-1],step)
    time = np.append(time, solution.t[-1])

    flight_history                      = {}
    flight_history['time']              = time
    flight_history['r']                 = np.zeros(len(time))
    flight_history['velocity']          = np.zeros(len(time))
    flight_history['flight_path_angle'] = np.zeros(len(time))
    flight_history['longitude']         = np.zeros(len(time))
    flight_history['mass']              = np.zeros(len(time))
    flight_history['axial_load_factor'] = np.zeros(len(time))
    flight_history['altitude']          = np.zeros(len(time))
    flight_history['qdyn']              = np.zeros(len(time))
    flight_history['flux']              = np.zeros(len(time))
    flight_history['AoA']               = np.zeros(len(time))
    flight_history['pitch_angle']       = np.zeros(len(time))
    flight_history['rho']               = np.zeros(len(time))
    flight_history['Pa']                = np.zeros(len(time))
    flight_history['Drag_coeff']        = np.zeros(len(time))
    flight_history['rho']               = np.zeros(len(time))
    flight_history['thrust']            = np.zeros(len(time))
    flight_history['thrust_vac']        = np.zeros(len(time))
    flight_history['mass_flow_rate']    = np.zeros(len(time))
    flight_history['g']                 = np.zeros(len(time))
    
    for i in range(len(time)):
        (flight_history['r'] [i],
         flight_history['velocity'][i],
         flight_history['flight_path_angle'][i],
         flight_history['longitude'][i],
         flight_history['mass'] [i],
         flight_history['axial_load_factor'][i],
         flight_history['qdyn'] [i],
         flight_history['flux'] [i],
         flight_history['altitude'][i],
         flight_history['AoA'][i],
         flight_history['pitch_angle'] [i],
         flight_history['rho'][i],
         flight_history['Pa'][i],
         flight_history['Drag_coeff'][i],
         flight_history['thrust'][i],
         flight_history['thrust_vac'][i],
         flight_history['mass_flow_rate'][i],
         flight_history['g'][i]) = fonction_ode_simu( time[i],
                                                      solution.sol(time[i]))
        
    return flight_history

def obj_function_trajectory(design_var):
    """
    Function defining the objective function for the optimization
    
    :design_var: design variables of the first stage (payload mass and control law) 
    
    """
    Payload_mass = design_var[0]*1e3
    return - Payload_mass


def constraint_function_trajectory(design_var, Parameters,Rocket_model):
    """
    Function defining the constraints function for the optimization
    
    :design_var: design variables of the first stage (payload mass and control law) 
    
    :Parameters: dictionary of parameters for the simulation
    
    :Rocket_model: dictionary composed of rocket parameters
    
    """
    
    # Integration of first stage trajectory
    flight_history_first_stage = trajectory_integration_first_stage(design_var, Parameters,Rocket_model)

    # Integration of second stage trajectory
    flight_history_second_stage = trajectory_integration_second_stage(design_var,Parameters,Rocket_model,flight_history_first_stage)
    
    
    # Calculation of discrepancy at the orbit

    discrepancy_altitude = np.abs(Parameters['Specifications']['Altitude']-flight_history_second_stage['altitude'][-1])
    discrepancy_velocity = np.abs(Parameters['Specifications']['Velocity']-flight_history_second_stage['velocity'][-1])
    discrepancy_flight_path_angle = np.abs(Parameters['Specifications']['Flight_path_angle']-flight_history_second_stage['flight_path_angle'][-1])
    
    #Definition of tolerances on state vector at the orbit injection
    tolerance_altitude = 1e3
    tolerance_velocity = 100.
    tolerance_flight_path_angle = 1.
    
    return np.array([(tolerance_altitude-discrepancy_altitude)/1e3, tolerance_velocity-discrepancy_velocity,
                     tolerance_flight_path_angle-discrepancy_flight_path_angle])




if __name__ == "__main__":
    
    
    Falcon_model={}
    Falcon_model['First_stage'] = {}
    Falcon_model['First_stage']['Propellant_mass'] = 395700.
    Falcon_model['First_stage']['Dry_mass'] = 25600.
    Falcon_model['First_stage']['Isp_vac'] = 312.
    Falcon_model['First_stage']['Mass_flow_rate'] = 9*298.8
    Falcon_model['First_stage']['Reference_area'] = np.pi*3.7 **2
    Falcon_model['First_stage']['Nozzle_exit_area'] = 9*(0.92/2)**2*np.pi

    Falcon_model['Second_stage'] = {}
    Falcon_model['Second_stage']['Propellant_mass'] = 92670.
    Falcon_model['Second_stage']['Dry_mass'] = 3900.
    Falcon_model['Second_stage']['Isp_vac'] = 348.
    Falcon_model['Second_stage']['Mass_flow_rate'] = 287.45
    Falcon_model['Second_stage']['Reference_area'] =  np.pi*3.7 **2
    Falcon_model['Second_stage']['Nozzle_exit_area'] = (3.3/2)**2*np.pi

    Falcon_model['Second_stage']['Fairing_mass'] = 1900.

    Parameters = {}
    Parameters['Control']={}
    Parameters['Control']['First_stage'] = {}
    Parameters['Control']['First_stage']['Vertical_phase_duration']=10.
    Parameters['Control']['First_stage']['Pitch_over_duration']=10.

    Parameters['Control']['First_stage']['Pitch_over_exp_decay_duration']=10.

    Parameters['Specifications'] = {}
    Parameters['Specifications']['Altitude']= 250*1e3
    Parameters['Specifications']['Flight_path_angle']= 0.
    Parameters['Specifications']['Velocity']= 7800.  

    Parameters['Control']['First_stage'] = {}
    Parameters['Control']['First_stage']['Vertical_phase_duration']=10.
    Parameters['Control']['First_stage']['Pitch_over_duration']=5.
    Parameters['Control']['First_stage']['Pitch_over_exp_decay_duration']=20.


    constraint_ = lambda x:  constraint_function_trajectory(x,Parameters,Falcon_model)

    design_var = np.array([3.,40.,5.,0.,1.])

    bounds_optim = ((0.,25.),(20.,50.),(-5.,10.),(-1,1.),(1.,2.))    


    t0 = time.time()
    sol_opt = scipy.optimize.fmin_slsqp(obj_function_trajectory,
                              design_var,
                              f_ieqcons = constraint_,
                              bounds = bounds_optim,disp=True)
    print(time.time()-t0)


    flight_history_first_stage = trajectory_integration_first_stage(sol_opt, Parameters,Falcon_model)
    flight_history_second_stage = trajectory_integration_second_stage(sol_opt,Parameters,Falcon_model,flight_history_first_stage)



    plt.subplot(2,3,1)
    plt.plot(np.concatenate((flight_history_first_stage['time'],flight_history_first_stage['time'][-1]+flight_history_second_stage['time'])),
             np.concatenate((flight_history_first_stage['altitude'],flight_history_second_stage['altitude'])))
    plt.title('Altitude')
    plt.grid()
    plt.subplot(2,3,2)

    plt.plot(np.concatenate((flight_history_first_stage['time'],flight_history_first_stage['time'][-1]+flight_history_second_stage['time'])),
             np.concatenate((flight_history_first_stage['velocity'],flight_history_second_stage['velocity'])))
    plt.title('Velocity')
    plt.grid()


    plt.subplot(2,3,3)
    plt.plot(np.concatenate((flight_history_first_stage['time'],flight_history_first_stage['time'][-1]+flight_history_second_stage['time'])),
             np.concatenate((flight_history_first_stage['flight_path_angle'],flight_history_second_stage['flight_path_angle'])))
    plt.title('Flight path angle')
    plt.grid()

    plt.subplot(2,3,4)
    plt.plot(np.concatenate((flight_history_first_stage['time'],flight_history_first_stage['time'][-1]+flight_history_second_stage['time'])),
             np.concatenate((flight_history_first_stage['pitch_angle'],flight_history_second_stage['pitch_angle'])))
    plt.title('Pitch angle')
    plt.grid()

    plt.subplot(2,3,5)
    plt.plot(np.concatenate((flight_history_first_stage['time'],flight_history_first_stage['time'][-1]+flight_history_second_stage['time'])),
             np.concatenate((flight_history_first_stage['mass'],flight_history_second_stage['mass'])))

    plt.title('Mass')
    plt.grid()

    plt.subplot(2,3,6)
    plt.plot(np.concatenate((flight_history_first_stage['time'],flight_history_first_stage['time'][-1]+flight_history_second_stage['time'])),
             np.concatenate((flight_history_first_stage['AoA'],flight_history_second_stage['AoA'])))

    plt.title('AoA')
    plt.grid()

    plt.figure()

    plt.subplot(2,3,1)
    plt.plot(np.concatenate((flight_history_first_stage['time'],flight_history_first_stage['time'][-1]+flight_history_second_stage['time'])),
             np.concatenate((flight_history_first_stage['qdyn'],flight_history_second_stage['qdyn'])))
    plt.title('qdyn')
    plt.grid()