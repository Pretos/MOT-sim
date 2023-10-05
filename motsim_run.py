from AtomEnsamble import *
from BoundingBox import *
from Beam import *
from MOTSim import *
import sys


args = {
	'mass': 174,
	'N': 1e5,
	'T': 485+273,
	'center': (0,0,0.5e-2),  # initial position
	'n_x': (0,0,1), 			  # initial position
	'radius': 2.5e-3,             # initial position
	'vmax': 10,			# initial velocity
	'sr': 0.3,          # initial velocity
	'n_v': (0,0,1),    # initial velocity
	'x1': -2e-2,			    # bounding box
	'x2': 2e-2,			        # bounding box
	'y1': -2e-2,			    # bounding box
	'y2': 2e-2,			        # bounding box
	'z2': 10e-2,			    # bounding box
	'dt': 1e-6,     # simulation
	'frames': 1500,  # simulation
}


name = 'sweep20230507_0/balanced_acute_6Gcm_10x'
args['save'] = name


# Yb 556
gamma = 182e3
S = 1600
Delta = -14*gamma
g = 1.5
lam = 556e-9
w0 = 0.5e-2

beams = []
# b = Beam(origin = (0,0,5e-2), n = (1,0,0), 
# 	     lam = lam, w0 = w0, gamma = gamma, S = S, Delta = Delta, g = g, Bp = 4*0.01)
# beams.append(b)
# b = Beam(origin = (0,0,5e-2), n = (-1,0,0), 
# 	     lam = lam, w0 = w0, gamma = gamma, S = S, Delta = Delta, g = g, Bp = 4*0.01)
# beams.append(b)
# b = Beam(origin = (0,0,5e-2), n = (0,1,0), 
# 	     lam = lam, w0 = w0/5, gamma = gamma, S = S, Delta = Delta, g = g, Bp = 8*0.01)
# beams.append(b)
# b = Beam(origin = (0,0,5e-2), n = (0,-1,0), 
# 	     lam = lam, w0 = w0/5, gamma = gamma, S = S, Delta = Delta, g = g, Bp = 8*0.01)
# beams.append(b)
# b = Beam(origin = (0,0,5e-2), n = (0,0,1), 
# 	     lam = lam, w0 = w0, gamma = gamma, S = S, Delta = Delta, g = g, Bp = 4*0.01)
# beams.append(b)
# b = Beam(origin = (0,0,5e-2), n = (0,0,-1), 
# 	     lam = lam, w0 = w0, gamma = gamma, S = S, Delta = Delta, g = g, Bp = 4*0.01)
# beams.append(b)

b = Beam(origin = (0,0,5e-2), n = (0,1,np.sqrt(3)), 
	     lam = lam, w0 = w0, gamma = gamma, S = S, Delta = Delta, g = g, Bp = 6*0.01)
beams.append(b)
b = Beam(origin = (0,0,5e-2), n = (0,-1,-np.sqrt(3)), 
	     lam = lam, w0 = w0, gamma = gamma, S = S, Delta = Delta, g = g, Bp = 6*0.01)
beams.append(b)

b = Beam(origin = (0,0,5e-2), n = (1,0,0), 
	     lam = lam, w0 = w0, gamma = gamma, S = S, Delta = Delta, g = g, Bp = 12*0.01)
beams.append(b)
b = Beam(origin = (0,0,5e-2), n = (-1,0,0), 
	     lam = lam, w0 = w0, gamma = gamma, S = S, Delta = Delta, g = g, Bp = 12*0.01)
beams.append(b)

b = Beam(origin = (0,0,5e-2), n = (0,1,-np.sqrt(3)), 
	     lam = lam, w0 = w0, gamma = gamma, S = S, Delta = Delta, g = g, Bp = 6*0.01)
beams.append(b)
b = Beam(origin = (0,0,5e-2), n = (0,-1,np.sqrt(3)), 
	     lam = lam, w0 = w0, gamma = gamma, S = S, Delta = Delta, g = g, Bp = 6*0.01)
beams.append(b)

args['beams'] = beams

box = MOTSim.simulate(args)

# box = BoundingBox.load(name)

box.animate_xyz_beams_sphereaperture(center=(0,0,5e-2), radius=0.5e-3, save=f'{name}')

# args = {
# 	'mass': 174,
# 	'N': 1e5,
# 	'T': 485+273,
# 	'center': (1e-2,-1e-2,5e-2),  # initial position
# 	'n_x': (-1,1,0), 			  # initial position
# 	'radius': 2.5e-3,             # initial position
# 	'vmax': 10,			# initial velocity
# 	'sr': 0.2,          # initial velocity
# 	'n_v': (-1,1,0),    # initial velocity
# 	'x1': -2e-2,			    # bounding box
# 	'x2': 2e-2,			        # bounding box
# 	'y1': -2e-2,			    # bounding box
# 	'y2': 2e-2,			        # bounding box
# 	'z2': 10e-2,			    # bounding box
# 	'dt': 1e-6,     # simulation
# 	'frames': 500,  # simulation
# }

# args = {
# 	'mass': 74,
# 	'N': 1e5,
# 	'T': 120+273,
# 	'center': (1e-2,-1e-2,5e-2),  # initial position
# 	'n_x': (-1,1,0), 			  # initial position
# 	'radius': 2.5e-3,             # initial position
# 	'vmax': 25,			# initial velocity
# 	'sr': 0.2,          # initial velocity
# 	'n_v': (-1,1,0),    # initial velocity
# 	'x1': -2e-2,			    # bounding box
# 	'x2': 2e-2,			        # bounding box
# 	'y1': -2e-2,			    # bounding box
# 	'y2': 2e-2,			        # bounding box
# 	'z2': 10e-2,			    # bounding box
# 	'dt': 1e-6,     # simulation
# 	'frames': 300,  # simulation
# }



# start = int(sys.argv[1])
# end = int(sys.argv[2])

# for Bp in np.linspace(0,75,25)[start:end]:
# 	Bp = Bp*0.01
# 	for Delta in np.linspace(0,3,20):
# 		# Yb 399
# 		# w0 = 0.5e-2
# 		# gamma = 28e6
# 		# S = 0.043
# 		# Delta = -Delta*gamma
# 		# g = 1
# 		# lam = 399e-9		

# 		# Rb 780
# 		w0 = 0.5e-2
# 		gamma = 36.1e6
# 		S = 2
# 		Delta = -Delta*gamma
# 		g = 1
# 		lam = 780e-9

# 		beams = []

# 		b = Beam(origin = (0,0,5e-2), n = (1,0,0), 
# 			     lam = lam, w0 = w0, gamma = gamma, S = S, Delta = Delta, g = g, Bp = Bp)
# 		beams.append(b)
# 		b = Beam(origin = (0,0,5e-2), n = (-1,0,0), 
# 			     lam = lam, w0 = w0, gamma = gamma, S = S, Delta = Delta, g = g, Bp = Bp)
# 		beams.append(b)
# 		b = Beam(origin = (0,0,5e-2), n = (0,1,0), 
# 			     lam = lam, w0 = w0, gamma = gamma, S = S, Delta = Delta, g = g, Bp = Bp)
# 		beams.append(b)
# 		b = Beam(origin = (0,0,5e-2), n = (0,-1,0), 
# 			     lam = lam, w0 = w0, gamma = gamma, S = S, Delta = Delta, g = g, Bp = Bp)
# 		beams.append(b)

# 		# push beam
# 		b = Beam(origin = (0,0,0), n = (0,0,1), 
# 			     lam = lam, w0 = 0.5e-2, gamma = gamma, S = S, Delta = Delta, g = g, Bp = 0)
# 		beams.append(b)

# 		args['beams'] = beams


# 		name = f'sweep20230419/780_{Bp*100:.2f}G_{Delta/gamma:.2f}_0.5w0_push'
# 		args['save'] = name

# 		print(name)

# 		box = MOTSim.simulate(args)

# simulate
# box.simulate(dt=2e-6, frames=250, cooling=True, save=f'out/{name}')

# Load previous simulation
# box = BoundingBox.load(f'out/{name}')


# name = '399_556_100G_0.5w0_push_small'
# box = BoundingBox.load(f'out/{name}')
# box.animate_xyz_beams_aperture(save=f'out/{name}', center=(0,0), radius=0.5e-2)

# animate
# box.animate_xyz_beams_aperture(save=f'out/{name}', center=(0,0), radius=0.5e-2)



# box = BoundingBox.load('simout')

# box.animate_xy()

# fig, ax = plt.subplots()

# x = np.array([[0,0,10e-2] for i in range(200)])
# speed = np.linspace(-10,10,200)
# v = np.array([[i*np.sqrt(2),0,i*np.sqrt(2)] for i in speed])

# ax.plot(speed, (box.beams[0].F(x,v)[:,2]+box.beams[1].F(x,v)[:,2])/(174*amu))

# plt.show()



# fig, ax = plt.subplots()

# v = np.array([[0,0,0] for i in range(200)])
# dist = np.linspace(-2e-2,2e-2,200)
# x = np.array([[i,0,10e-2+i] for i in dist])

# ax.plot(dist, box.beams[0].F(x,v)[:,2]+box.beams[1].F(x,v)[:,2])

# plt.show()
