from AtomEnsamble import *
from tqdm import tqdm
from matplotlib.animation import FuncAnimation, FFMpegWriter
import pickle
import time

class BoundingBox:

	def __init__(self, x1, x2, y1, y2, z2, atom_ensamble):
		'''
		Specify the bounding box in which the AtomEnsamble is restricted to be within. 
		There is no z1 because we define the oven at z=0

		Parameters
		-----
		x1: lowest allowed x value; m
		x2: highest allowed x value; m
		atoms: AtomEnsamble
		'''
		self.x1 = x1
		self.x2 = x2
		self.y1 = y1
		self.y2 = y2
		self.z1 = 0
		self.z2 = z2
		self.atoms = atom_ensamble

		self.beams = []
		# position, velocity, and time of the atoms over time
		self.xdata = []
		self.vdata = []
		self.tdata = []

	def simulate(self, dt, frames, steps_per_frame=100, save='simout', cooling=True):
		'''
		Simulates the evolution of the atomic enamble under the force of gravity and the MOT beams.

		Parameters
		-----
		dt: time step; s
		frames: the position/velocity of each particle is recorded at each frame and out-of-bounds atoms are removed
		steps_per_frame: the number of time steps before the next frame is recorded
		save: file to save simulation data to
		cooling: if False, the cooling forces are turned off
		'''
		# clear data
		self.xdata = []
		self.vdata = []
		self.tdata = np.arange(frames+1)*dt*steps_per_frame
		self.axdata = [[]] # position of atoms that pass through the aperture
		self.avdata = [[]] # velocity of atoms that pass through the aperture

		f = tqdm
		if not save:
			f = lambda x: x

		for i in f(range(frames)):
			# save a copy of the current state
			self.xdata.append(np.copy(self.atoms.x))
			self.vdata.append(np.copy(self.atoms.v))

			# evolve for steps_per_frame
			for _ in range(steps_per_frame):
				self.atoms.evolve(dt, self.beams, cooling)

			# records all atoms that pass through the right (+z) wall (i.e for counting atoms going through an aperture during animation)
			aper = np.argwhere(self.atoms.x[:,2]>self.z2)
			aperi = aper.flatten()
			self.axdata.append(self.atoms.x[aperi])
			self.avdata.append(self.atoms.v[aperi])
			self.atoms.x = np.delete(self.atoms.x, aper, axis=0)
			self.atoms.v = np.delete(self.atoms.v, aper, axis=0)

			# remove out-of-bounds atoms
			oob = np.argwhere((self.atoms.x[:,0]<self.x1) | (self.atoms.x[:,0]>self.x2) | (self.atoms.x[:,1]<self.y1) | (self.atoms.x[:,1]>self.y2) | (self.atoms.x[:,2]<self.z1))
			self.atoms.x = np.delete(self.atoms.x, oob, axis=0)
			self.atoms.v = np.delete(self.atoms.v, oob, axis=0)

		self.xdata.append(np.copy(self.atoms.x))
		self.vdata.append(np.copy(self.atoms.v))
		if save:
			with open(f'{save}.pkl', "wb") as f:
				pickle.dump(self, f)

	def load(file):
		'''
		Load previous simulation from pickle file.

		Parameters
		-----
		file: picke file name

		Returns
		-----
		BoundingBox object
		'''
		with open(f'{file}.pkl', "rb") as f:
			return pickle.load(f)

	def animate_xyz_beams(self, save='out'):
		'''
		Animates the currently loaded calculation in a 3 views.

		Parameters
		-----
		save: filename to save animation
		'''
		start_time = time.time()

		fig, ax = plt.subplots(1,3,figsize=(20,5),dpi=200)

		plot1, = ax[0].plot([],[], f'.{self.atoms.color}', alpha=0.25, mew=0, ms=2)
		plot2, = ax[1].plot([],[], f'.{self.atoms.color}', alpha=0.25, mew=0, ms=2)
		plot3, = ax[2].plot([],[], f'.{self.atoms.color}', alpha=0.25, mew=0, ms=2)

		p1 = []
		p2 = []
		p3 = []

		for beam in self.beams:
			plot, = ax[0].plot([],[], f'.{beam.color}', alpha=0.25, mew=0, ms=2)
			p1.append(plot)
			plot, = ax[1].plot([],[], f'.{beam.color}', alpha=0.25, mew=0, ms=2)
			p2.append(plot)
			plot, = ax[2].plot([],[], f'.{beam.color}', alpha=0.25, mew=0, ms=2)
			p3.append(plot)

		ax[0].set_ylim(self.x1*1.05, self.x2*1.05)
		ax[0].set_xlim(0, self.z2*1.05)
		ax[0].set_xlabel('z [m]')
		ax[0].set_ylabel('x [m]')
		ax[0].set_title('t=0')

		ax[1].set_xlim(self.x1*1.05, self.x2*1.05)
		ax[1].set_ylim(self.y1*1.05, self.y2*1.05)
		ax[1].set_xlabel('x [m]')
		ax[1].set_ylabel('y [m]')

		ax[2].set_xlim(self.z1*1.05, self.z2*1.05)
		ax[2].set_ylim(self.y1*1.05, self.y2*1.05)
		ax[2].set_xlabel('z [m]')
		ax[2].set_ylabel('y [m]')

		plt.tight_layout()

		def animate(i):
			xdata = self.xdata[i]
			for j in range(len(self.beams)):
				beam = self.beams[j]
				r = np.linalg.norm(np.cross(beam.n, xdata-beam.origin), axis=1)
				highlight = np.argwhere(r<beam.w0)

				p1[j].set_data(xdata[highlight,2], xdata[highlight,0])
				p2[j].set_data(xdata[highlight,0], xdata[highlight,1])
				p3[j].set_data(xdata[highlight,2], xdata[highlight,1])

				xdata = np.delete(xdata, highlight, axis=0)

			plot1.set_data(xdata[:,2], xdata[:,0])
			plot2.set_data(xdata[:,0], xdata[:,1])
			plot3.set_data(xdata[:,2], xdata[:,1])

			ax[0].set_title(f't={self.tdata[i]*1e3:.3f} ms')

			return plot1,plot2,plot3

		anim = FuncAnimation(fig,animate,frames=len(self.xdata),blit=True)

		writervideo = FFMpegWriter(fps=25)
		anim.save(f'{save}.mp4', writer=writervideo)
		plt.close()  
		print(f'Animation: {time.time()-start_time}')
	
	def animate_xyz_beams_aperture(self, center, radius, save='out'):
		'''
		Animates the currently loaded calculation in a 3 views.

		Parameters
		-----
		save: filename to save animation
		center: of the aperture located at the far right (+z); (x,y)
		radius: radius of the aperture cross-section 
		'''
		start_time = time.time()

		fig, ax = plt.subplots(2,3,figsize=(20,10),dpi=200)

		plot1, = ax[0][0].plot([],[], f'.{self.atoms.color}', alpha=0.25, mew=0, ms=2)
		plot2, = ax[0][1].plot([],[], f'.{self.atoms.color}', alpha=0.25, mew=0, ms=2)
		plot3, = ax[0][2].plot([],[], f'.{self.atoms.color}', alpha=0.25, mew=0, ms=2)

		p1 = []
		p2 = []
		p3 = []

		for beam in self.beams:
			plot, = ax[0][0].plot([],[], f'.{beam.color}', alpha=0.25, mew=0, ms=2)
			p1.append(plot)
			plot, = ax[0][1].plot([],[], f'.{beam.color}', alpha=0.25, mew=0, ms=2)
			p2.append(plot)
			plot, = ax[0][2].plot([],[], f'.{beam.color}', alpha=0.25, mew=0, ms=2)
			p3.append(plot)

		ax[0][0].set_ylim(self.x1*1.05, self.x2*1.05)
		ax[0][0].set_xlim(0, self.z2*1.05)
		ax[0][0].set_xlabel('z [m]')
		ax[0][0].set_ylabel('x [m]')
		ax[0][0].set_title('t=0')

		ax[0][1].set_xlim(self.x1*1.05, self.x2*1.05)
		ax[0][1].set_ylim(self.y1*1.05, self.y2*1.05)
		ax[0][1].set_xlabel('x [m]')
		ax[0][1].set_ylabel('y [m]')
		aperture = plt.Circle(center, radius ,fill = False, alpha=0.5)
		ax[0][1].add_artist(aperture)

		ax[0][2].set_xlim(self.z1*1.05, self.z2*1.05)
		ax[0][2].set_ylim(self.y1*1.05, self.y2*1.05)
		ax[0][2].set_xlabel('z [m]')
		ax[0][2].set_ylabel('y [m]')

		counts = []
		speeds = [[]]
		speeds_z = [[]]
		for i in range(len(self.xdata)):
			x = self.axdata[i]
			v = self.avdata[i]
			if len(x)>0:
				aper = np.argwhere(np.sqrt((x[:,0]-center[0])**2+ (x[:,1]-center[1])**2)<radius)
				counts.append(len(aper))
				vaper = v[aper.flatten(),:]
				speeds.append(np.append(speeds[-1], np.sqrt(vaper[:,0]**2+vaper[:,1]**2+vaper[:,2]**2)))
				speeds_z.append(np.append(speeds_z[-1], vaper[:,2]))
			else:
				counts.append(0)
				speeds.append(speeds[-1])
				speeds_z.append(speeds_z[-1])

		counts = np.cumsum(counts)/self.atoms.N*self.atoms.frac

		plot11, = ax[1][0].plot([],[], 'k')
		ax[1][0].set_xlabel('Time [s]')
		ax[1][0].set_ylabel('Cumulative count through aperture [fractional]')

		plot12, = ax[1][1].plot([],[], '.-k')
		ax[1][1].set_xlabel('Speed [m/s]')
		ax[1][1].set_ylabel('Cumulative fraction')

		plot13, = ax[1][2].plot([],[], '.-k')
		ax[1][2].set_xlabel('Z-velocity [m/s]')
		ax[1][2].set_ylabel('Cumulative fraction')
		plt.tight_layout()

		def animate(i):
			xdata = self.xdata[i]
			for j in range(len(self.beams)):
				beam = self.beams[j]
				r = np.linalg.norm(np.cross(beam.n, xdata-beam.origin), axis=1)
				highlight = np.argwhere(r<beam.w0)

				p1[j].set_data(xdata[highlight,2], xdata[highlight,0])
				p2[j].set_data(xdata[highlight,0], xdata[highlight,1])
				p3[j].set_data(xdata[highlight,2], xdata[highlight,1])

				xdata = np.delete(xdata, highlight, axis=0)

			plot1.set_data(xdata[:,2], xdata[:,0])
			plot2.set_data(xdata[:,0], xdata[:,1])
			plot3.set_data(xdata[:,2], xdata[:,1])

			ax[0][0].set_title(f't={self.tdata[i]*1e3:.3f} ms')

			plot11.set_data(self.tdata[:i], counts[:i])
			ax[1][0].relim()
			ax[1][0].autoscale(True)

			plot12.set_data(*self.hist(speeds[i]))
			ax[1][1].relim()
			ax[1][1].autoscale(True)

			plot13.set_data(*self.hist(speeds_z[i]))
			ax[1][2].relim()
			ax[1][2].autoscale(True)

		anim = FuncAnimation(fig,animate,frames=len(self.xdata),blit=False, init_func=lambda: None)

		writervideo = FFMpegWriter(fps=25)
		anim.save(f'{save}.mp4', writer=writervideo)
		plt.close()  
		print(f'Animation: {time.time()-start_time}')

	def animate_xyz_beams_sphereaperture(self, center, radius, save='out'):
		'''
		Animates the currently loaded calculation in a 3 views.

		Parameters
		-----
		save: filename to save animation
		'''
		start_time = time.time()

		fig, ax = plt.subplots(2,3,figsize=(20,10),dpi=200)

		plot1, = ax[0][0].plot([],[], f'.{self.atoms.color}', alpha=0.25, mew=0, ms=2)
		plot2, = ax[0][1].plot([],[], f'.{self.atoms.color}', alpha=0.25, mew=0, ms=2)
		plot3, = ax[0][2].plot([],[], f'.{self.atoms.color}', alpha=0.25, mew=0, ms=2)

		p1 = []
		p2 = []
		p3 = []

		for beam in self.beams:
			plot, = ax[0][0].plot([],[], f'.{beam.color}', alpha=0.25, mew=0, ms=2)
			p1.append(plot)
			plot, = ax[0][1].plot([],[], f'.{beam.color}', alpha=0.25, mew=0, ms=2)
			p2.append(plot)
			plot, = ax[0][2].plot([],[], f'.{beam.color}', alpha=0.25, mew=0, ms=2)
			p3.append(plot)

		ax[0][0].set_ylim(self.x1*1.05, self.x2*1.05)
		ax[0][0].set_xlim(0, self.z2*1.05)
		ax[0][0].set_xlabel('z [m]')
		ax[0][0].set_ylabel('x [m]')
		ax[0][0].set_title('t=0')

		ax[0][1].set_xlim(self.x1*1.05, self.x2*1.05)
		ax[0][1].set_ylim(self.y1*1.05, self.y2*1.05)
		ax[0][1].set_xlabel('x [m]')
		ax[0][1].set_ylabel('y [m]')
		aperture = plt.Circle(center, radius ,fill = False, alpha=0.5)
		ax[0][1].add_artist(aperture)

		ax[0][2].set_xlim(self.z1*1.05, self.z2*1.05)
		ax[0][2].set_ylim(self.y1*1.05, self.y2*1.05)
		ax[0][2].set_xlabel('z [m]')
		ax[0][2].set_ylabel('y [m]')

		counts = []
		speeds = [[]]
		# speeds_z = [[]]
		for i in range(len(self.xdata)):
			x = self.xdata[i]
			aper = np.argwhere(np.sqrt((x[:,0]-center[0])**2+(x[:,1]-center[1])**2+(x[:,2]-center[2])**2)<radius)
			counts.append(len(aper))
			# print(aper.flatten())
			# vaper = self.vdata[aper.flatten(),:]
			# speeds.append(np.append(speeds[-1], np.sqrt(vaper[:,0]**2+vaper[:,1]**2+vaper[:,2]**2)))

		counts = np.array(counts)/self.atoms.N*self.atoms.frac

		plot11, = ax[1][0].plot([],[], 'k')
		ax[1][0].set_xlabel('Time [s]')
		ax[1][0].set_ylabel('Counts in aperture region [fractional]')

		# plot12, = ax[1][1].plot([],[], '.-k')
		# ax[1][1].set_xlabel('Speed [m/s]')
		# ax[1][1].set_ylabel('Fraction')

		# plot13, = ax[1][2].plot([],[], '.-k')
		# ax[1][2].set_xlabel('Z-velocity [m/s]')
		# ax[1][2].set_ylabel('Cumulative fraction')
		plt.tight_layout()

		def animate(i):
			xdata = self.xdata[i]
			for j in range(len(self.beams)):
				beam = self.beams[j]
				r = np.linalg.norm(np.cross(beam.n, xdata-beam.origin), axis=1)
				highlight = np.argwhere(r<beam.w0)

				p1[j].set_data(xdata[highlight,2], xdata[highlight,0])
				p2[j].set_data(xdata[highlight,0], xdata[highlight,1])
				p3[j].set_data(xdata[highlight,2], xdata[highlight,1])

				xdata = np.delete(xdata, highlight, axis=0)

			plot1.set_data(xdata[:,2], xdata[:,0])
			plot2.set_data(xdata[:,0], xdata[:,1])
			plot3.set_data(xdata[:,2], xdata[:,1])

			ax[0][0].set_title(f't={self.tdata[i]*1e3:.3f} ms')

			plot11.set_data(self.tdata[:i], counts[:i])
			ax[1][0].relim()
			ax[1][0].autoscale(True)

			# plot12.set_data(*self.hist(speeds[i]))
			# ax[1][1].relim()
			# ax[1][1].autoscale(True)

			# plot13.set_data(*self.hist(speeds_z[i]))
			# ax[1][2].relim()
			# ax[1][2].autoscale(True)

		anim = FuncAnimation(fig,animate,frames=len(self.xdata),blit=False, init_func=lambda: None)

		writervideo = FFMpegWriter(fps=25)
		anim.save(f'{save}.mp4', writer=writervideo)
		plt.close()  
		print(f'Animation: {time.time()-start_time}')

	def hist(self, arr, bins=200):
	    '''
	    returns (position, counts) for plotting histogram
	    '''
	    count, pos = np.histogram(arr, bins=bins)
	    pos = pos[:-1] + (pos[1]-pos[0])/2
	    return pos, count/self.atoms.N*self.atoms.frac