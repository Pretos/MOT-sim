from AtomEnsamble import *
from scipy.constants import pi, h, c, hbar

class Beam:

	def __init__(self, origin, n, lam, w0, gamma, S, Delta, g, Bp, color='b'):
		'''
		A single laser MOT beam (we do not assume retroreflection)

		Parameters
		-----
		origin: point along the beam where B=0. We assume a uniform field gradient Bp along the beam.
		n: the propagation direction of the beam 
		lam: wavelength in;
		w0: beam waist
		gamma: cooling transition linewidth; Hz
		S: saturation parameter I/I_sat
		Delta: detuning; Hz
		g: Lande g-factor
		Bp: magnetic field gradient; T/m
		color: color of beam when plotting
		'''
		self.origin = origin
		# normalize
		self.n = np.array(n)/np.linalg.norm(n)
		# wave vector
		self.k = 2*np.pi/lam * self.n

		self.w0 = w0
		self.gamma = gamma
		self.g = g
		self.Bp = Bp
		self.S = S
		self.Delta = Delta
		self.color = color

	def F(self, x, v):
		'''
		Cooling force for this beam
		Krzysztof, K., Xuan, K.D., Małgorzata, G., Huy, B.N. and Jerzy, S., 2010. Magneto-optical trap: fundamentals and realization. CMST, (2), pp.115-129.

		Parameters
		-----
		x: positions of atoms
		v: velocities of atoms
		'''

		# return np.zeros(np.shape(x))

		xp = x - self.origin

		# https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
		r = np.linalg.norm(np.cross(self.n, xp), axis=1)

		s = self.S * self.gaussian_S(r)

		F = hbar*self.gamma*s/(1 + s + 4*((self.Delta - np.dot(v, self.k) - self.g*uB/h*self.Bp*np.dot(xp, self.n))/self.gamma)**2)/2
		F = np.outer(F, self.k)

		return F

	# def F(self, x, v):
	# 	'''
	# 	Cooling force for this beam
	# 	Krzysztof, K., Xuan, K.D., Małgorzata, G., Huy, B.N. and Jerzy, S., 2010. Magneto-optical trap: fundamentals and realization. CMST, (2), pp.115-129.

	# 	Parameters
	# 	-----
	# 	x: positions of atoms
	# 	v: velocities of atoms
	# 	'''
	# 	x = np.array([[0.1e-2,0,0] for i in range(200)])
	# 	z = np.linspace(-4,4,200)
	# 	v = np.array([[0,0,i] for i in z])
	# 	# print(v)
	# 	r = self.dist_to_beam(x)
	# 	print(r)

	# 	s = self.gaussian_S(r)
	# 	# print(s)

	# 	F = hbar*self.gamma*s/(1 + s + 4*((self.Delta - np.dot(v, self.k))/self.gamma)**2)/2
	# 	F = np.outer(F, self.k)
	# 	# F = hbar*self.k*self.gamma/2*s/(1+s+4*((self.Delta - np.dot(v, self.k))/self.gamma)**2)

	# 	F2 = hbar*self.gamma*s/(1 + s + 4*((self.Delta - np.dot(v, -self.k))/self.gamma)**2)/2
	# 	F2 = np.outer(F2, -self.k)
	# 	F += F2

	# 	print(self.k)

	# 	fig, ax = plt.subplots()

	# 	ax.plot(z, F[:,2]/(174*amu))

	# 	plt.show()

	# def F2(self, x, v):
	# 	'''
	# 	Cooling force for this beam
	# 	Krzysztof, K., Xuan, K.D., Małgorzata, G., Huy, B.N. and Jerzy, S., 2010. Magneto-optical trap: fundamentals and realization. CMST, (2), pp.115-129.

	# 	Parameters
	# 	-----
	# 	x: positions of atoms
	# 	v: velocities of atoms
	# 	'''
	# 	xs = np.linspace(-1e-2, 1e-2, 200)
	# 	x = np.array([[0,0,i] for i in xs])
	# 	v = np.array([[0,0,0] for i in range(200)])
	# 	# print(v)
	# 	r = self.dist_to_beam(x)
	# 	# print(r)

	# 	s = self.gaussian_S(r)
	# 	print(s)

	# 	F = hbar*self.gamma*s/(1 + s + 4*((self.Delta - np.dot(v, self.k) - self.g*uB/h*self.Bp*np.dot(x, self.n))/self.gamma)**2)/2
	# 	F = np.outer(F, self.k)

	# 	# print(self.g*uB/h*self.Bp*np.dot(x, self.n))
	# 	# print(np.dot(x, self.n))
	# 	# F2 = hbar*self.gamma*s/(1 + s + 4*((self.Delta - np.dot(v, -self.k))/self.gamma)**2)/2
	# 	# F2 = np.outer(F2, -self.k)
	# 	# F += F2

	# 	# print(self.k)

	# 	# F = np.dot(x, self.n)

	# 	fig, ax = plt.subplots()

	# 	ax.plot(xs, F[:,2])

	# 	# ax.plot(z, F[:,2]/(174*amu))

	# 	ax.plot()

	# 	plt.show()

	# 	# I SATURATE
	# 	# print()

	def gaussian_S(self, r):
	    '''
	    Gaussian weighting of beam intensity

	    Parameters
	    -----
	    w0: beam waist
	    r: distance from beam
	    '''
	    return np.exp(-2*(r/self.w0)**2)