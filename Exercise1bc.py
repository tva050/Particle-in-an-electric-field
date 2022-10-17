import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc

#Plot style and font

plt.style.use('ggplot')

""" Inital values """

q = 1.   # Charge
m = 1.   # Mass

dt = 1e-2 
k = 1000 # Time duration
t = 1000 # Time duration

""" Start inital """

# Creating a new array of given shape, filled with velocity and position. 
v = np.zeros((k+1, 2))  # Velocity with k rows and 2 colums 
v_1 = np.zeros((k+1, 2))
x = np.zeros((k+1, 2))  # Position   -----||-----
x_1 = np.zeros((k+1, 2))  
#Then we set the initial values for the velocity and position   
v[0] = np.array([0., 1.])
v_1[0] = np.array([0., 1.])
x[0] = np.array([1., 0.])
x_1[0] = np.array([1., 0.])

# Magnetic field set to zero and a array of electric field
B = 0.  # Magnetic field set to zero, defined as a constant
E = -np.zeros((k+1))  
E_1 = -np.zeros((k+1))

# electric field at particle initial position
E[0] = -x[0, 0] 
E_1[0] = -np.sin(x_1[0,0])

# Defining theta
theta = - ((q * B * dt) / m) 

""" Main """
# Making a function for the matrix that have rotation in it
def rotation(angel):
    return np.array([[np.cos(angel), -np.sin(angel)], [np.sin(angel), np.cos(angel)]])

# Then we take half a rotation bakwards , with the rotation (Thate_min = -theta/2).
V_min = np.matmul(rotation(- theta / 2), v[0])
V_min1 = np.matmul(rotation(- theta / 2), v_1[0])
# Half a rotation bakwards followed by half a decerlartion
v[0,0] = V_min[0] - (q * dt / 2 * m) * E[0]
v_1[0,0] = V_min1[0] - (q * dt / 2 * m) * E_1[0]
v[0,1] = V_min[1] 
v_1[0,1] = V_min1[1]

""" 
Boris mover algorithm, here the velocity is rotaded by half bakwards.
"""
for k in range(k):
    # Performing a halv acc and then store it in v_min
    v_min = np.array([v[k, 0] + ((q * dt) / (2 * m)) * E[k], v[k, 1]])
    
    #Performing a rotation with angle theta, then store it in v_pos 
    v_pos = np.matmul(rotation(theta), v_min)
    
    #Performing another half acc to get velocity in x and y direction
    v[k+1, 0] = (v_pos[0]) + (((q * dt) / (2 * m)) * E[k])
    v[k+1, 1] = v_pos[1]
    
    #Updating x to the next time step
    x[k+1, 0] = x[k, 0] + (dt * v[k+1, 0])
    x[k+1, 1] = x[k, 1] + (dt * v[k+1, 1])
    
    #Updating the electric field
    E[k+1] = -x[k+1,0]
 
# For loop for the second particle in E(x) = -sin(x)   
for t in range(t):
    # Performing a halv acc and then store it in v_min
    v_min1 = np.array([v_1[t, 0] + ((q * dt) / (2 * m)) * E_1[t], v_1[t, 1]])
    
    #Performing a rotation with angle theta, then store it in v_pos 
    v_pos1 = np.matmul(rotation(theta), v_min1)
    
    #Performing another half acc to get velocity in x and y direction
    v_1[t+1, 0] = (v_pos1[0]) + (((q * dt) / (2 * m)) * E_1[t])
    v_1[t+1, 1] = v_pos1[1]
    
    #Updating x to the next time step
    x_1[t+1, 0] = x_1[t, 0] + (dt * v_1[t+1, 0])
    x_1[t+1, 1] = x_1[t, 1] + (dt * v_1[t+1, 1])
    
    #Updating the electric field
    E_1[t+1] = -np.sin(x[t+1,0])
 
# Plotting the results   
plt.plot(x[:, 0], x[:,1], label = "$E(x)=-\epsilon x$" )
plt.plot(x_1[:, 0], x_1[:,1], label = "$E(x)=-\epsilon sin(x)$" )
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.show()
    