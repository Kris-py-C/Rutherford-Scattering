# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 14:10:11 2024

@author: kc671
"""

import numpy as np
import matplotlib.pyplot as plt

plt.figure(dpi=144)

q1 = 79
q2 = 2
m = 7294.3
#k = 1    ##atomic units. k will be omitted in all formulae
c = 137.035999
length = 1
time = 1
dt = 1e-5

r_y = -0.005
v_x = 0
v_y = 9.592520

inputs = input("Change initial conditions?(yes): ").lower()
if inputs == "yes":
    r_y = float(input("radius y: "))
    v_x = float(input("velo x: "))
    v_y = float(input("velo y: "))


#intial y = -0.005
#initial x = 0 ->> 0.001
#initial v_x = 0
#initial y_x = 9.592520 
def component(r,R):
    r_0 = r/R
    return r_0 
    
def trajectory(r_x,r_y,v_x,v_y):
    n = 0
    x_bound = r_x
    y_bound = r_y
    x_pos = []
    y_pos = []
    r_s = []
    deg = []
    bound = np.sqrt(x_bound**2 + y_bound**2)*1.1
    R = np.sqrt(r_x**2 + r_y**2)
    r_s.append(R)
    
    while R < bound:
        if n == 1000:
            return r_s
        else:
            F = (q1*q2) / (R**2)
            F_x = F * component(r_x,R)
            F_y = F * component(r_y,R)
        
            v_x = v_x + (F_x*dt)/m
            v_y = v_y + (F_y*dt)/m ## a = F / m. 
            deflect_deg = np.tan( v_y / v_x )
            deg.append(deflect_deg)
            r_s.append(R)       
            
            r_x = r_x + v_x*dt
            r_y = r_y + v_y*dt
            x_pos.append(r_x)
            y_pos.append(r_y)
            n+=1
            R = np.sqrt(r_x**2 + r_y**2)
    plt.plot(x_pos,y_pos)
    return r_s,deg

R_S,deg_RS = trajectory(0,-0.005,0,9.592520)


closest = np.min(R_S)
neg_close = np.negative(closest)
# print("theory = {}".format(closest))
closest_xs = np.arange(neg_close,closest+0.000001,0.000001)
closest_ys = np.sqrt(closest**2 - closest_xs**2)

closest_xs = closest * np.sin(np.linspace(0,2*np.pi,1000))
closest_ys = closest * np.cos(np.linspace(0,2*np.pi,1000))

# e_k_initial = e_v_closest

print(closest)

for value in np.arange(0, 0.001+0.0001 ,0.0001):
    rad,degs = trajectory(value,r_y,v_x,v_y)
    print(np.max(degs))


plt.scatter(0,0,s=10,color="y",label="Gold nucleus")
plt.gca().set_aspect('equal')
plt.plot(closest_xs,closest_ys,color="black",label="Closest approach")
#plt.plot(closest_xs,-closest_ys,color="y")
plt.legend()
plt.title("Model trajectory of incident alpha on gold nucleus")
plt.xlabel(r'$x$-position($r_0$)')
plt.ylabel(r'$y$-position($r_0$)')
plt.show()