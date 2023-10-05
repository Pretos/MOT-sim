from AtomEnsamble import *
from BoundingBox import *
from Beam import *
from MOTSim import *
import sys

# Bps = np.linspace(0,75,25)
# Deltas = np.linspace(0,3,20)

# xx, yy = np.meshgrid(Bps, Deltas)
# zz = np.zeros(np.shape(xx))

# for xi in tqdm(range(len(Bps))):
# 	for yi in range(len(Deltas)):

# 		gamma = 28e6
# 		name = f'sweep20230419/780_{Bps[xi]:.2f}G_{-Deltas[yi]:.2f}_0.5w0_push'
# 		box = BoundingBox.load(name)

# 		center = (0,0)
# 		radius = 0.5e-2

# 		counts = []
# 		for i in range(len(box.xdata)):
# 			x = box.axdata[i]
# 			if len(x)>0:
# 				aper = np.argwhere(np.sqrt((x[:,0]-center[0])**2+ (x[:,1]-center[1])**2)<radius)
# 				counts.append(len(aper))

# 		counts = np.sum(counts)/box.atoms.N*box.atoms.frac
# 		zz[yi,xi] = counts

# with open(f'sweep20230419/counts.pkl', "wb") as f:
# 	pickle.dump((xx,yy,zz), f)

# with open(f'sweep20230419/counts.pkl', "rb") as f:
# 	load = pickle.load(f)
# xx, yy, zz = load[0], load[1], load[2]

# fig, ax = plt.subplots()
 
# p = ax.pcolormesh(xx, yy, zz)
# plt.colorbar(p)
# ax.set_xlabel("B' [G/cm]")
# ax.set_ylabel('$\Delta/\gamma$')
# ax.set_title('Fraction captured by aperture: Rb, 2D 780 nm, w0=0.5cm,\naperture radius=0.5, push beam (same $\Delta$), T=120C')

# plt.subplots_adjust(bottom=0.2)
# plt.text(0, -0.2, 'sweep20230419/counts.pkl', fontsize=9, c='lightgray', transform=ax.transAxes)

# plt.savefig('fig.png', bbox_inches='tight', dpi=300)

# plt.show()


########################################################
# velocity stuff

Bps = np.linspace(15,75,20)
Deltas = np.linspace(0,3,20)

xx, yy = np.meshgrid(Bps, Deltas)
zz = np.zeros(np.shape(xx))

# for xi in tqdm(range(len(Bps))):
# 	for yi in range(len(Deltas)):

# 		gamma = 28e6
# 		name = f'sweep20230418/399_{Bps[xi]:.2f}G_{-Deltas[yi]:.2f}_0.5w0_push'
# 		box = BoundingBox.load(name)

# 		center = (0,0)
# 		radius = 0.5e-2

# 		speeds = []
# 		for i in range(len(box.xdata)):
# 			x = box.axdata[i]
# 			v = box.avdata[i]
# 			if len(x)>0:
# 				aper = np.argwhere(np.sqrt((x[:,0]-center[0])**2+ (x[:,1]-center[1])**2)<radius)
# 				vaper = v[aper.flatten(),:]
# 				speeds.extend(np.sqrt(vaper[:,0]**2+vaper[:,1]**2+vaper[:,2]**2))

# 		zz[yi,xi] = np.mean(speeds)

# with open(f'sweep20230419/counts.pkl', "wb") as f:
# 	pickle.dump((xx,yy,zz), f)



with open(f'sweep20230418/counts.pkl', "rb") as f:
	load = pickle.load(f)
xx, yy, zz = load[0], load[1], load[2]

fig, ax = plt.subplots()
 
p = ax.pcolormesh(xx, yy, zz)
plt.colorbar(p)
ax.set_xlabel("B' [G/cm]")
ax.set_ylabel('$\Delta/\gamma$')
ax.set_title('Average velocity: Yb, 2D 399 nm, w0=0.5cm,\naperture radius=0.5, push beam (same $\Delta$), T=120C')

plt.subplots_adjust(bottom=0.2)
plt.text(0, -0.2, 'sweep20230419/counts.pkl', fontsize=9, c='lightgray', transform=ax.transAxes)

plt.savefig('fig.png', bbox_inches='tight', dpi=300)

plt.show()