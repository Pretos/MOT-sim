import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import *
from scipy.special import erf
import ipywidgets as widgets
from scipy.interpolate import interp1d as interp
from matplotlib.animation import FuncAnimation, FFMpegWriter
from IPython import display
import matplotlib
from tqdm import tqdm

matplotlib.rcParams['animation.embed_limit'] = 2**128

amu = physical_constants['atomic mass constant'][0]
kb = physical_constants['Boltzmann constant'][0]
uB = physical_constants['Bohr magneton'][0]

m_rb = 87*amu
m_yb = 174*amu

k_yb = 2*pi/399e-9
k_ybt = 2*pi/556e-9
k_rb = 2*pi/780e-9

sat_rb = 1.7 # mW/cm^2
sat_yb = 63.1
sat_ybt = 0.139

# I/I_sat
i_yb = .043
i_ybt = 27
i_rb = 2

# linewidth 
gamma_yb = 28e6
gamma_ybt = 0.182e6
gamma_rb = 5.75e6

# detuning
delta_yb = -0.7*gamma_yb
delta_ybt = -5.5*gamma_ybt
delta_rb = -0.7*gamma_rb

# lande g factors
g_yb = 1
g_ybt = 1.5
g_rb = 1.5

def MB(v, T, m):
    '''
    Maxwell-Boltzmann distribution
    
    Parameters
    -----
    v: m/s
    T: k
    m: kg
    '''
    return (m/(2*pi*kb*T))**1.5 * 4*pi*v**2 * np.exp(-m*v**2/(2*kb*T))

def MB_CDF(v,m,T):
    """ Cumulative Distribution function of the Maxwell-Boltzmann speed distribution """
    a = np.sqrt(kb*T/m)
    return erf(v/(np.sqrt(2)*a)) - np.sqrt(2/np.pi)* v* np.exp(-v**2/(2*a**2))/a

def gaussian(r, d):
    '''
    d: beam width
    r: distance from beam
    '''
    return np.exp(-(r/d)**2/2)

def hist(arr, bins=200):
    '''
    returns (position, counts) for plotting histogram
    '''
    count, pos = np.histogram(arr, bins=bins)
    pos = pos[:-1] + (pos[1]-pos[0])/2
    return pos, count

'''
Variables
-----
dt: time interval of simulation
'''
T = 485+273
mass = 174*amu
Bp = 15*0.01 # G/cm to T/m

N = int(1e5)

dt = 1e-6

print(f'T: {T} K')
print(f'N: {N:,}')

# generate velocity magnitudes
vs = np.linspace(0,10,1000)
cdf = MB_CDF(vs, mass, T)
cdf /= np.max(cdf)
inv_cdf = interp(cdf,vs)
rand_nums = np.random.random(N)
speeds = inv_cdf(rand_nums)
    
# fig, ax = plt.subplots()
# ax.hist(speeds, bins=100)
# # ax.set_xlim(0,20)
# # ax.set_ylim(0,50)
# ax.set_title(f'Initial speeds N={N}')
# ax.set_xlabel('Speeds [m/s]')
# ax.set_ylabel('Counts')
# plt.savefig('fig.png', bbox_inches='tight')
# plt.show()

d1 = 10e-2
d3 = 5e-2
d2 = 3e-2
d4 = 10e-2
d5 = 10e-2

# MOT beam diameter
D = 1e-2

# initial velocity and position #

# RADIATE 
v = np.zeros((N,3))
theta = np.random.rand(N)*np.pi
phi = np.random.rand(N)*np.pi
v[:,0] = speeds*np.cos(theta)
v[:,1] = speeds*np.sin(theta)*np.cos(phi)
v[:,2] = speeds*np.sin(theta)*np.sin(phi)

x = np.zeros((N,3))
x[:,2] = -d1


steps_per_frame = 100

def MOTp(k, v, z, i, gamma, delta, g, Bp):
    '''
    MOT force (F+) in the direction of the laser beam
    
    Parameters
    -----
    g: g-factor
    Bp: B field gradient
    v: speed along direction of the beam
    z: displacement wrt to MOT center in the direction of the beam
    '''
    return hbar*k*gamma/2*i/(1+i+4*((delta-k*v-g*uB/h*Bp*z)/gamma)**2)

def MOTm(k, v, z, i, gamma, delta, g, Bp):
    '''
    MOT force (F-) in the direction of the laser beam
    
    Parameters
    -----
    g: g-factor
    Bp: B field gradient
    v: speed along direction of the beam
    z: displacement wrt to MOT center in the direction of the beam
    '''
    return -hbar*k*gamma/2*i/(1+i+4*((delta+k*v+g*uB/h*Bp*z)/gamma)**2)

def F(x, v): 
#     return 0*v
    F = np.zeros((len(x),3))
    
    #399
    r = np.sqrt(x[:,1]**2 + x[:,2]**2) # distance from beam center    
    i = i_yb*gaussian(r, D/2)
    F[:,0] += MOTp(k_yb, v[:,0], x[:,0], i, gamma_yb, delta_yb, g_yb, Bp)
    F[:,0] += MOTm(k_yb, v[:,0], x[:,0], i, gamma_yb, delta_yb, g_yb, Bp)
    
    r = np.sqrt(x[:,0]**2 + x[:,2]**2) # distance from beam center    
    i = i_yb*gaussian(r, D/2)
    F[:,1] += MOTp(k_yb, v[:,1], x[:,1], i, gamma_yb, delta_yb, g_yb, Bp)
    F[:,1] += MOTm(k_yb, v[:,1], x[:,1], i, gamma_yb, delta_yb, g_yb, Bp)
    
    r = np.sqrt(x[:,0]**2 + x[:,1]**2) # distance from beam center    
    i = i_yb*gaussian(r, D/2)
    F[:,2] += MOTp(k_yb, v[:,2], x[:,2], i, gamma_yb, delta_yb, g_yb, Bp)
    F[:,2] += MOTm(k_yb, v[:,2], x[:,2], i, gamma_yb, delta_yb, g_yb, Bp)
    
    #556
    r = np.sqrt(x[:,1]**2 + x[:,2]**2) # distance from beam center    
    i = i_ybt*gaussian(r, D/2)
    F[:,0] += MOTp(k_ybt, v[:,0], x[:,0], i, gamma_ybt, delta_ybt, g_ybt, Bp)
    F[:,0] += MOTm(k_ybt, v[:,0], x[:,0], i, gamma_ybt, delta_ybt, g_ybt, Bp)
    
    r = np.sqrt(x[:,0]**2 + x[:,2]**2) # distance from beam center    
    i = i_ybt*gaussian(r, D/2)
    F[:,1] += MOTp(k_ybt, v[:,1], x[:,1], i, gamma_ybt, delta_ybt, g_ybt, Bp)
    F[:,1] += MOTm(k_ybt, v[:,1], x[:,1], i, gamma_ybt, delta_ybt, g_ybt, Bp)
    
    r = np.sqrt(x[:,0]**2 + x[:,1]**2) # distance from beam center    
    i = i_ybt*gaussian(r, D/2)
    F[:,2] += MOTp(k_ybt, v[:,2], x[:,2], i, gamma_ybt, delta_ybt, g_ybt, Bp)
    F[:,2] += MOTm(k_ybt, v[:,2], x[:,2], i, gamma_ybt, delta_ybt, g_ybt, Bp)
    
    return F

# fig, ax = plt.subplots(2,3,figsize=(25,10), dpi=200)

# ax[0][0].set_xlabel('z [m]')
# ax[0][0].set_ylabel('x [m]')

# ax[0][0].set_xlim(-d1*1.1, (d2+d4)*1.1)
# ax[0][0].set_ylim(-d3*1.1, d3*1.1)

# ax[0][0].axvline(-d1, c='k')
# ax[0][0].axvline(d2+d4, c='k')
# ax[0][0].axhline(d3, c='k')
# ax[0][0].axhline(-d3, c='k')
# ax[0][0].plot(-d1, 0, '.k')

# ax[0][0].set_title('Time = 0')

# ax[0][1].set_xlabel('z [m]')
# ax[0][1].set_ylabel('y [m]')

# ax[0][1].set_xlim(-d1*1.1, (d2+d4)*1.1)
# ax[0][1].set_ylim(-d5*1.1, d5*1.1)

# ax[0][1].axvline(-d1, c='k')
# ax[0][1].axvline(d2+d4, c='k')
# ax[0][1].axhline(d5, c='k')
# ax[0][1].axhline(-d5, c='k')
# ax[0][1].plot(-d1, 0, '.k')

# ax[1][0].set_xlabel('Vx [m/s]')
# ax[1][0].set_ylabel('Counts')
# ax[1][1].set_xlabel('Vy [m/s]')
# ax[1][1].set_ylabel('Counts')
# ax[1][2].set_xlabel('Vz [m/s]')
# ax[1][2].set_ylabel('Counts')
# ax[0][2].set_xlabel('Speed [m/s]')
# ax[0][2].set_ylabel('Counts')

# 399
# ax[0][0].axvline(-D, c='b')
# ax[0][0].axvline(D, c='b')
# ax[0][1].axvline(-D, c='b')
# ax[0][1].axvline(D, c='b')
# ax[0][1].add_patch(plt.Circle((0, 0), D, color='b', fill=False))

# animated plots
# plot1, = ax[0][0].plot(x[:,2], x[:,0], '.r', alpha=0.15)
# plot2, = ax[0][1].plot(x[:,2], x[:,1], '.r', alpha=0.15)

# plot3, = ax[1][0].plot(*hist(v[:,0]))
# plot4, = ax[1][1].plot(*hist(v[:,1]))
# plot5, = ax[1][2].plot(*hist(v[:,2]))
# plot6, = ax[0][2].plot(*hist(np.sqrt(v[:,0]**2+v[:,1]**2+v[:,2]**2)))

# ax[1][0].set_xlim(-10,10)
# ax[1][0].set_ylim(-1e4,1e4)
# ax[1][1].set_xlim(-10,10)
# ax[1][1].set_ylim(-1e4,1e4)
# ax[1][2].set_xlim(-10,10)
# ax[1][2].set_ylim(-1e4,1e4)

# plt.tight_layout()

frames = 2
x_data = []
v_data = []

for i in range(frames):
    print(f'Frame: {i}', end='\r')
    
    # removing particles that go out of bounds
    oob = np.argwhere(np.abs(x[:,0])>d3)
    x = np.delete(x, oob, axis=0)
    v = np.delete(v, oob, axis=0)
    oob = np.argwhere(np.abs(x[:,1])>d5)
    x = np.delete(x, oob, axis=0)
    v = np.delete(v, oob, axis=0)
    oob = np.argwhere(x[:,2]>d2+d4)
    x = np.delete(x, oob, axis=0)
    v = np.delete(v, oob, axis=0)

    # advance steps_per_frame steps before plotting next frame
    for _ in range(steps_per_frame):
        x += v*dt
        v += (F(x,v)/mass - (0,9.8,0))*dt # cooling forces and gravity on Earth!

    x_data.append(x)
    v_data.append(v)

# # RUN SIMULATION
# def animate(i):
#     global x, v
#     print(f'Frame: {i}', end='\r')
    
#     # removing particles that go out of bounds
#     oob = np.argwhere(np.abs(x[:,0])>d3)
#     x = np.delete(x, oob, axis=0)
#     v = np.delete(v, oob, axis=0)
#     oob = np.argwhere(np.abs(x[:,1])>d5)
#     x = np.delete(x, oob, axis=0)
#     v = np.delete(v, oob, axis=0)
#     oob = np.argwhere(x[:,2]>d2+d4)
#     x = np.delete(x, oob, axis=0)
#     v = np.delete(v, oob, axis=0)
    
#     plot1.set_data(x[:,2], x[:,0])
#     plot2.set_data(x[:,2], x[:,1])
    
#     plot3.set_data(*hist(v[:,0]))
#     plot4.set_data(*hist(v[:,1]))
#     plot5.set_data(*hist(v[:,2]))
#     plot6.set_data(*hist(np.sqrt(v[:,0]**2+v[:,1]**2+v[:,2]**2)))

# #     acc = F(x,v)/mass - (0,9.8,0)
# #     plot3.set_data(v[:,0], acc[:,0])
# #     plot4.set_data(v[:,0], acc[:,0])
# #     plot5.set_data(v[:,0], acc[:,0])
# # #     plot6.set_data(*hist(np.sqrt(acc[:,0]**2+acc[:,1]**2+acc[:,2]**2)))
    
#     ax[0][2].relim()
#     ax[0][2].autoscale(True)
#     ax[1][0].relim()
#     ax[1][0].autoscale(True)
#     ax[1][1].relim()
#     ax[1][1].autoscale(True)
#     ax[1][2].relim()
#     ax[1][2].autoscale(True)

#     ax[0][0].set_title(f'Time = {i*steps_per_frame*dt*1e3:.1f} ms')

#     # advance steps_per_frame steps before plotting next frame
#     for _ in range(steps_per_frame):
#         x += v*dt
#         v += (F(x,v)/mass - (0,9.8,0))*dt # cooling forces and gravity on Earth!
    
#     return plot1, plot2,
    
# anim = FuncAnimation(fig,animate,frames=2,interval=50,blit=True)

# writervideo = FFMpegWriter(fps=20)
# anim.save('out.mp4', writer=writervideo)
# plt.close()    