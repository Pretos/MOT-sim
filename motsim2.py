from AtomEnsamble import *
from BoundingBox import *
from Beam import *

# Define an ensamble of atoms
atoms = AtomEnsamble(mass = 174, 
					 N = 1e5, 
					 T = 485 + 273)
atoms.init_pos(center=(3e-2,3e-2,5e-2),
			   n = (1,1,0),
			   radius = 2.5e-3)
atoms.init_vel(vmax = 10, 
			   sr = 0.1,
			   n = (-1,-1,0))

# Define a bounding box and associate an atom ensamble to the box
box = BoundingBox(-5e-2, 5e-2, -5e-2, 5e-2, 40e-2, atoms)

# Define MOT beams and associate with bounding box

# Yb 399
# w0 = 0.5e-2
# gamma = 28e6
# S = 0.043
# Delta = -0.7*28e6
# g = 1
# Bp = 50*0.01
# lam = 399e-9

# b = Beam(origin = (0,0,0), n = (0,0,1), 
# 	     lam = lam, w0 = 0.5e-2, gamma = gamma, S = S, Delta = gamma, g = g, Bp = 0)
# box.beams.append(b)


# fig, ax = plt.subplots()

# x = np.array([[0,0,10e-2] for i in range(200)])
# speed = np.linspace(-10,10,200)
# v = np.array([[0,0,i] for i in speed])

# ax.plot(speed, box.beams[0].F(x,v)[:,2]/(174*amu))

# plt.show()



# fig, ax = plt.subplots()

# v = np.array([[0,0,0] for i in range(200)])
# dist = np.linspace(-2e-2,2e-2,200)
# x = np.array([[i,0,10e-2+i] for i in dist])

# ax.plot(dist, box.beams[0].F(x,v)[:,2]+box.beams[1].F(x,v)[:,2])

# plt.show()




a = np.array([[1,2,3], [1.5,3,4]])


oob = np.argwhere(a[:,0]<1.5).flatten()

print(oob)


print(a[oob])