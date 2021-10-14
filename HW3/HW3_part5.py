"""
Created on Sat Sep 11 20:05:02 2021

@author: nikita
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from pynput import keyboard
from pynput.keyboard import Listener as kl

#Set up the planet and physical laws
planet_pos = np.array([0.0, 0.0])
planet_mass = 20
planet_r = 2
G = 1

#Simulation parameters
max_steps = 5000
dt = 0.1

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

fig,ax = plt.subplots(figsize = (6,6))
plt.xlim([planet_pos[0] - 20, planet_pos[0] + 20])
plt.ylim([planet_pos[1] - 20, planet_pos[1] + 20])

cc = plt.Circle(planet_pos, planet_r, color=(0.4,0.3,0.1))
ax.add_artist(cc)

lp, = plt.plot([],[], markersize = 1)
rp, = plt.plot([],[],'o', markersize = 5)  

# Onboard thruster thrust initialization
F_mono = 0
F_ion = 0


def init() :
    pass

def advance(i):
    global rocket_fuel_mass, F_mono, F_ion
    
    lp.set_data(rocket_pos[0:i,0],rocket_pos[0:i,1])
    col = (0,1,0)
    lp.set_markerfacecolor(col)
    lp.set_color(col)
    
    #if the animation index reaches the end of the position data storage arrays
    if(i >= max_steps - 1):
        #move all the data back one index to make room.
        rocket_vel[0:max_steps - 1, :] = rocket_vel[1:max_steps, :]
        rocket_pos[0:max_steps - 1, :] = rocket_pos[1:max_steps, :]
        #also move the index back.
        i = max_steps - 2
    
    #If the rocket has fuel
    if(rocket_fuel_mass > 0):
        dm = thrust_dm_dt*dt #The mass used is the mass flow times timestep
        
        #if that used mass is more than the rocket has, use remaining mass.
        if(dm > rocket_fuel_mass): dm = rocket_fuel_mass
            
        rocket_fuel_mass -= dm #subtract used mass
        F_thrust_mag = (dm/dt)*thrust_vel #calculate thrust force
        
    else:
        rocket_fuel_mass = 0
        F_thrust_mag = 0
        
    rocket_mass = rocket_dry_mass + rocket_fuel_mass
    
    vel_mag = np.sqrt(rocket_vel[i,0]**2+rocket_vel[i,1]**2)
    
    vel_dir = rocket_vel[i,:]/vel_mag
    
    F_thrust_mag += F_mono + F_ion #add force from bonus thrusters
    
    F_thrust = F_thrust_mag*vel_dir
    
    r_vec = rocket_pos[i,:] - planet_pos
    
    r_mag = np.sqrt(np.dot(r_vec, r_vec))
    r_dir = r_vec/r_mag
    
    F_gravity = -G*planet_mass*rocket_mass*r_dir/r_mag**2
    
    F = F_thrust + F_gravity
    
    if(i == 0): rocket_vel[i,:] -= F/rocket_mass*0.5*dt
    
    rocket_vel[i+1,:] = rocket_vel[i,:] + F/rocket_mass*dt
    
    rocket_pos[i + 1,:] = rocket_pos[i,:] + rocket_vel[i+1,:]*dt

    rp.set_data(rocket_pos[i,0], rocket_pos[i,1])
    
    f = vel_mag/10;
    f = min(max(f,0),1)
    
    color = (0.5 + 0.5*f, 0.1 + 0.7*f, 0.1 + 0.7*f)
    
    if(F_mono > 0): #red for forward mono thrust
        color = (1, 0, 0)
    elif(F_mono < 0): #blue for backward mono thrust
        color = (0, 1, 0)
    elif(F_ion != 0): #purple for any ion thrust
        color = (1.0, 0, 0.7)
        
    #set mono thrust to zero after it is used, since it is an impulse
    F_mono = 0;
    
    rp.set_markerfacecolor(color)
    rp.set_color(color)

    plt.title("Distance: %5.2f, Speed: %5.2f"% (r_mag*planet_r, vel_mag))

run = True; #Thisi will stay true until the esc key is released
def gen(): #iterator that runs while the run variable is true
    i = 0 #initial loop counter value (local to the function?)
    while(run): #loop while run is true
        yield i #return to FuncAnimation, calls update (?)
        i += 1 #increment counter

ani = FuncAnimation(fig, advance, frames = gen, init_func = init, 
                    save_count = 0, interval = 0.1, repeat = False)

def on_press(key):
    global F_mono, F_ion
    if(key == keyboard.Key.up):
        F_mono = 0.01
    elif(key == keyboard.Key.down):
        F_mono = -0.01
    elif(key == keyboard.Key.insert): #the insert key switches ion on and off
        if(F_ion == 0):
            F_ion = 5e-7
        else:
            F_ion = 0
            
def on_release(key):
    global run, k_l
    
    if key == keyboard.Key.esc:
        run = False; #what does the semicolon do in python??
        k_l.stop()
        print("Donezo")
        return False

k_l = kl(on_press = on_press, on_release = on_release)
k_l.start()

print("Use 'insert' to toggle ion thruster, and up and down keys to fire "
      " forward and backward monoprop impulsive burns. Esc to quit.")

