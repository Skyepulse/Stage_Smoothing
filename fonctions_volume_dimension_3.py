%matplotlib inline
%matplotlib notebook

!pip install plotly

import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
import random as rand
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import colorsys
from matplotlib.tri import Triangulation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.special import betainc


def nSphereVolume(dim, iterations):
    count_in_sphere = 0

    for count_loops in range(iterations):
        point = np.random.uniform(-1.0, 1.0, dim)
        distance = np.linalg.norm(point)
        if distance < 1.0:
            count_in_sphere += 1

    return np.power(2.0, dim) * (count_in_sphere / iterations)

def nconevolume(iterations, c, theta):
    count_in_cone = 0
    points = []
    
    for count_loops in range(iterations):
        z = np.random.uniform(-1.0, 1.0)
        rho = np.random.uniform(0, 2.0)
        phi = np.random.uniform(0, np.pi)
        
        if (z > c) and (rho <= z*np.tan(theta)):
            count_in_cone += 1
            x1 = z
            x2 = rho*np.cos(phi)
            x3 = rho*np.sin(phi)
            tuple_ = (x1, x2, x3)
            points.append(tuple_)
           
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection='3d')
    for i in points:
        ax.scatter(i[0], i[1], i[2], color = 'green')
        
    return np.power(2.0, 3) * (count_in_cone / iterations)

def nconevolume_2(iterations, c, theta):
    count_in_cone = 0
    points = []
    
    for count_loops in range(iterations):
        x1 = np.random.uniform(-2.0, 2.0)
        x2 = np.random.uniform(-2.0, 2.0)
        x3 = np.random.uniform(-2.0, 2.0)
        
        z = x1
        rho = np.sqrt(x2**2 + x3**2) 
        if (z > c) and (rho <= z*np.tan(theta)):
            count_in_cone += 1
            tuple_ = (x1, x2, x3)
            points.append(tuple_)
           
    #fig = plt.figure(figsize=(6, 6))
    #ax = fig.add_subplot(111, projection='3d')
    #for i in points:
        #ax.scatter(i[0], i[1], i[2], color = 'green')
        
    return np.power(4.0, 3) * (count_in_cone / iterations)

#Testons avec une boule uniforme de rayon r
def ntest_ball_sample(rayon, c, theta, iterations):
    count_in_sphere = 0
    count_in_cone = 0
    points_sphere = []
    points_cone = []

    for count_loops in range(iterations):
        point = np.random.uniform(-rayon, rayon, 3)
        distance = np.linalg.norm(point)
        if distance < rayon:
            points_sphere.append(point)
            count_in_sphere += 1
    #fig = plt.figure(figsize=(6, 6))
    #ax = fig.add_subplot(111, projection='3d')
    #for i in points_sphere:
        #ax.scatter(i[0], i[1], i[2], color = 'green')
    
    vol_boule = np.power(2*rayon, 3) * (count_in_sphere / iterations)
    
    for point in points_sphere:
        x1 = point[0]
        x2 = point[1]
        x3 = point[2]
        
        z = x1
        rho = np.sqrt(x2**2 + x3**2) 
        if (z > c) and (rho <= z*np.tan(theta)):
            count_in_cone += 1
            tuple_ = (x1, x2, x3)
            points_cone.append(tuple_)
    vol_in_cone = vol_boule * (count_in_cone/count_in_sphere)
    
    #fig = plt.figure(figsize=(6, 6))
    #ax = fig.add_subplot(111, projection='3d')
    #for i in points_cone:
        #ax.scatter(i[0], i[1], i[2], color = 'orange')
    vol_out_cone = vol_boule - vol_in_cone
    
    return vol_in_cone, vol_out_cone, vol_boule

def single_noise_déplacement_epsilon_uniforme(rayon, epsilon, c, theta, iterations):
    count_in_sphere = 0
    count_in_cone = 0
    points_sphere = []
    points_cone = []

    for count_loops in range(iterations):
        point = np.random.uniform(-rayon, rayon, 3)
        distance = np.linalg.norm(point)
        if distance < rayon:
            points_sphere.append(point)
            count_in_sphere += 1
    #fig = plt.figure(figsize=(6, 6))
    #ax = fig.add_subplot(111, projection='3d')
    #for i in points_sphere:
        #ax.scatter(i[0], i[1], i[2], color = 'green')
    
    vol_boule = np.power(2*rayon, 3) * (count_in_sphere / iterations)
    
    for point in points_sphere:
        x1 = point[0]
        x2 = point[1]
        x3 = point[2]
        
        z = x1
        rho = np.sqrt(x2**2 + x3**2) 
        if (z > c) and (rho <= z*np.tan(theta)):
            count_in_cone += 1
            tuple_ = (x1, x2, x3)
            points_cone.append(tuple_)
    vol_in_cone = vol_boule * (count_in_cone/count_in_sphere)
    
    #fig = plt.figure(figsize=(6, 6))
    #ax = fig.add_subplot(111, projection='3d')
    #for i in points_cone:
        #ax.scatter(i[0], i[1], i[2], color = 'orange')
    vol_out_cone = vol_boule - vol_in_cone
    
    #probabilité que les points de la boule ne soient pas dans le cone:
    p1 = 1 - count_in_cone/count_in_sphere
    
    #on déplace de epsilon vers la droite de l'axe du cone
    a = (2*rayon + 1)/2
    b = 1/2
    x = 1 - (epsilon/(2*rayon))**2
    Vcap = (1/2) * vol_boule * betainc(a, b, x) #On utilise la fonction beta incomplete de la librairie scipy.
    
    Psn = p1 - 1 + 2*(Vcap/vol_boule)
    return Psn

def p_true_déplacement_epsilon_uniforme(rayon, epsilon, c, theta, iterations):
    #On parle de vol_deplace/vol_boule où vol_déplacé sont les points tels que après une translation epsilon soient dans la boule
    #mais pas dans le cone.
    count_in_sphere = 0
    count_in_transl = 0
    points_sphere = []
    points_transl = []
    points_not_in_cone = []
    points_in_cone = []
    count_not_in_cone = 0
    
    #Set de points dans la boule:
    for count_loops in range(iterations):
        point = np.random.uniform(c - rayon, c + rayon, 3)
        distance = np.linalg.norm(point)
        if distance < rayon:
            points_sphere.append(point)
            count_in_sphere += 1
    vol_boule = np.power(2*rayon, 3) * (count_in_sphere / iterations)
    
    #fig = plt.figure(figsize=(6, 6))
    #ax = fig.add_subplot(111, projection='3d')
    #for i in points_sphere:
        #ax.scatter(i[0], i[1], i[2], color = 'orange')
        
    #on fait la translation sur x, on teste si on est toujours dans la boule originelle
    for i in points_sphere:
        i[0] = i[0] + epsilon
        point = [i[0], i[1], i[2]]
        distance = np.linalg.norm(point)
        if distance < rayon:
            points_transl.append(point)
            count_in_transl +=1
    
    #on teste si on est en dehors ou non du cone:
    for point in points_transl:
        x1 = point[0]
        x2 = point[1]
        x3 = point[2]
        
        z = x1
        rho = np.sqrt(x2**2 + x3**2) 
        if (z > c) and (rho <= z*np.tan(theta)):
            points_in_cone.append(point)
        else:
            count_not_in_cone += 1
            points_not_in_cone.append(point)
    #on calcule alors le volume:
    vol_deplace = vol_boule*(count_not_in_cone/count_in_sphere)
        
    p_true = vol_deplace/vol_boule
        
    #for i in points_in_cone:
        #ax.scatter(i[0], i[1], i[2], color = 'green')
        
    return p_true

def visualise_cone(c, theta, iterations):
    points = []
    
    for count_loops in range(iterations):
        x1 = np.random.uniform(-2.0, 2.0)
        x2 = np.random.uniform(-2.0, 2.0)
        x3 = np.random.uniform(-2.0, 2.0)
        
        z = x1
        rho = np.sqrt(x2**2 + x3**2) 
        if (z > c) and (rho <= z*np.tan(theta)):
            tuple_ = (x1, x2, x3)
            points.append(tuple_)
           
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection='3d')
    for i in points:
        ax.scatter(i[0], i[1], i[2], color = 'green')
        
def diff_p_true_p_single_noise(rayon, epsilon, c, theta, iterations):
    return(single_noise_déplacement_epsilon_uniforme(rayon, epsilon, c, theta, iterations) - p_true_déplacement_epsilon_uniforme(rayon, epsilon, c, theta, iterations))

r = 1
epsilon = 0.5
c = 0
iterations = 100000
angles = np.linspace(0, np.pi/2, 100, endpoint = False)
list_results = []
for i in range(len(angles)):
    list_results.append(diff_p_true_p_single_noise(r, epsilon, c, angles[i], iterations))

plt.plot(angles, list_results)
plt.xlabel('Thetas (radians)')
plt.ylabel('P_single_noise - P_true')
plt.show()
