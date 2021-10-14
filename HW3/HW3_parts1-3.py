"""
Created on Sat Sep 11 20:05:02 2021

@author: nikita
"""

import numpy as np
import matplotlib.pyplot as plt

#Set up the planet and physical laws
planet_pos = np.array([0.0, 0.0])
planet_mass = 20
planet_r = 2
G = 1

#Simulation parameters
max_steps = 200
dt = 0.05

#user inputs
launch_angle = 65
rocket_dry_mass = 0.001
rocket_fuel_mass = 0.001
thrust_vel = 2
thrust_dm_dt = 1e-4
vel0 = 3.4

#show output
print("Launch Angle: %g"%launch_angle)
#the arguments after the string get placed at the %g. There is an error if you 
#dont have the same number of %g's as arguments.
print("Rocket Mass: %g dry, %g fuel" %(rocket_dry_mass, rocket_fuel_mass))

#allocate arrays to store position and velocity
rocket_pos = np.zeros((max_steps, 2))
rocket_vel = np.zeros((max_steps, 2))

#set initial position to the "top" of the planet
rocket_pos[0,:] = [planet_pos[0], planet_pos[1] + planet_r]

##set initial velocity to be v0*[cos(angle), isn(angle)]
rocket_vel[0,:] = [vel0*np.sin(launch_angle*np.pi/180),
          vel0*np.cos(launch_angle*np.pi/180)]

#Simulation loop
for i in range(max_steps-1):
    #If the rocket has fuel
    if(rocket_fuel_mass > 0):
        dm = thrust_dm_dt*dt #The mass used is the mass flow times timestep
        
        #if that used mass is more than the rocket has, use remaining mass.
        if(dm > rocket_fuel_mass):
            dm = rocket_fuel_mass
            
        rocket_fuel_mass -= dm #subtract used mass
        F_thrust_mag = (dm/dt)*thrust_vel #calculate thrust force
        
    else:
        rocket_fuel_mass = 0
        F_thrust_mag = 0
        
    rocket_mass = rocket_dry_mass + rocket_fuel_mass
    
    vel_mag = np.sqrt(rocket_vel[i,0]**2+rocket_vel[i,1]**2)
    
    vel_dir = rocket_vel[i,:]/vel_mag
    
    F_thrust = F_thrust_mag*vel_dir
    
    r_vec = rocket_pos[i,:] - planet_pos
    
    r_mag = np.sqrt(np.dot(r_vec, r_vec))
    r_dir = r_vec/r_mag
    
    F_gravity = -G*planet_mass*rocket_mass*r_dir/r_mag**2
    
    F = F_thrust + F_gravity
    
    if(i == 0): rocket_vel[i,:] -= F/rocket_mass*0.5*dt
    
    rocket_vel[i+1,:] = rocket_vel[i,:] + F/rocket_mass*dt
    
    rocket_pos[i + 1,:] = rocket_pos[i,:] + rocket_vel[i+1,:]*dt

    #End of the for loop    

fig,ax = plt.subplots(figsize = (6,6))
plt.xlim([planet_pos[0] - 20, planet_pos[0] + 20])
plt.ylim([planet_pos[1] - 20, planet_pos[1] + 20])

cc = plt.Circle(planet_pos, planet_r, color=(0.4,0.3,0.1))
ax.add_artist(cc)

lp, = plt.plot(rocket_pos[:,0],rocket_pos[:,1])
rp, = plt.plot(rocket_pos[-1,0],rocket_pos[-1,1],'o', markersize = 2)

        
