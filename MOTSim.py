from BoundingBox import *
from multiprocessing import Process, Queue


class MOTSim:

	def __init__(self, args, queue=None):
		return

	def simulate(args):

		# Define an ensamble of atoms
		atoms = AtomEnsamble(mass = args['mass'], 
							 N = args['N'], 
							 T = args['T'])
		atoms.init_pos(center = args['center'],
					   n = args['n_x'],
					   radius = args['radius'])
		atoms.init_vel(vmax = args['vmax'], 
					   sr = args['sr'],
					   n = args['n_v'])

		# Define a bounding box and associate an atom ensamble to the box
		box = BoundingBox(args['x1'], args['x2'], args['y1'], args['y2'], args['z2'], atoms)
		box.beams = args['beams']

		box.simulate(dt=args['dt'], frames=args['frames'], cooling=True, save=f"{args['save']}")
		return box

# 	def simulate_multiprocess(args, num_process=5):
# 		'''
# 		Multiprocess simulation
# 		'''
# 		atoms = AtomEnsamble(mass = args['mass'], 
# 							 N = args['N'], 
# 							 T = args['T'])
# 		atoms.init_pos(center = args['center'],
# 					   n = args['n_x'],
# 					   radius = args['radius'])
# 		atoms.init_vel(vmax = args['vmax'], 
# 					   sr = args['sr'],
# 					   n = args['n_v'])

# 		Ns = np.ones(num_process)*args['N']//num_process
# 		Ns[-1] += args['N']%num_process

# 		queue = Queue()
# 		args['save'] = None
# 		processes = []
# 		for N in Ns:
# 			args['N'] = N
# 			p = Process(target=MOTSim.simulate, args=(args, queue))
# 			processes.append(p)

# 		for p in processes:
# 			p.start()
# 			p.join()

# if __name__ ==  '__main__':
# 	MOT