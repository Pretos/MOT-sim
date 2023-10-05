import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import physical_constants
from scipy.special import erf
from scipy.interpolate import interp1d as interp


amu = physical_constants['atomic mass constant'][0]
kb = physical_constants['Boltzmann constant'][0]
uB = physical_constants['Bohr magneton'][0]


class AtomEnsamble:
    '''
    Coordinates
    -----
    The center of the atom source is the origin. The optical axis is along +z and gravity is along -y (standard optical coordinates)
    '''

    def __init__(self, mass, N, T, color='r'):
        '''
        Parameters
        -----
        mass: amu
        N: number of atoms in ensamble
        T: temperature of the atoms out of the oven; K
        color: color of atoms when plotting
        '''
        self.mass = mass * amu
        self.N = int(N)
        self.T = T
        self.color = color

    def init_pos(self, center, radius, n):
        '''
        Initialize the positions of atoms. We assume the cross section of the atoms out of the oven
        is circular

        Parameters
        -----
        center: of the atomic beam 
        radius: of the cross section of atoms coming out of the oven in meters
        n: normal vector to the atomic cross section 
        '''
        # https://mathworld.wolfram.com/DiskPointPicking.html
        r = radius * np.sqrt(np.random.random(self.N))
        theta = np.random.random(self.N) * 2 * np.pi

        # generate positions assuming n=(0,0,1)
        self.x = np.zeros((self.N, 3))
        self.x[:, 0] = r * np.cos(theta)
        self.x[:, 1] = r * np.sin(theta)

        # set n
        R = self.R_matrix(n)
        for i in range(self.N):
            self.x[i] = R@self.x[i]

        # move to center position
        self.x += np.array(center)

    def init_vel(self, vmax, sr, n):
        '''
        Initialize the velocities of atoms. We assume a Maxwell-Boltzmann distribution

        Parameters
        -----
        vmax: maximum cutoff of the MB distribution to consider; m/s
        sr: steradians that the initial velocity radiates into. For example, if sr=0, the beam is collimated
        n: the (center) direction of the atoms
        '''
        # first draw speeds from the MD distribution
        vs = np.linspace(0, vmax, 10000)
        cdf = self._MB_CDF(vs, self.T)
        self.frac = cdf[-1]
        cdf /= cdf[-1]
        inv_cdf = interp(cdf, vs)
        rand_nums = np.random.random(self.N)
        speeds = inv_cdf(rand_nums)

        # assign a direction for each particle
        # https://mathworld.wolfram.com/SpherePointPicking.html
        theta = sr/2*np.random.random(self.N)
        phi = np.arccos(2*np.random.random(self.N)-1)

        # direction centered on (0,0,1)
        self.v = np.zeros((self.N, 3))
        self.v[:,2] = speeds*np.cos(theta)
        self.v[:,0] = speeds*np.sin(theta)*np.cos(phi)
        self.v[:,1] = speeds*np.sin(theta)*np.sin(phi)

        R = self.R_matrix(n)
        for i in range(self.N):
            self.v[i] = R@self.v[i]

    def R_matrix(self, n):
        '''
        Calculates rotation matrix from (0,0,1) to n
        '''
        # calculate rotation matrix to n direction
        # https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
        n /= np.linalg.norm(n)
        if n[2]==1:
            return np.eye(3)
        v = np.cross((0,0,1), n)
        s = np.linalg.norm(v)
        c = np.dot((0,0,1), n)

        v = np.array([[0,-v[2],v[1]],
                      [v[2],0,-v[0]],
                      [-v[1],v[0],0]])
        R = np.eye(3) + v + v@v*(1-c)/s**2
        return R

    def evolve(self, dt, beams, cooling):
        '''
        Evolve the atom ensamble for one time step WITH GRAVITY.

        Parameters
        -----
        dt: time step; s
        beams: mot beams
        cooling: if False, the cooling forces are turned off
        '''
        self.x += self.v*dt

        force = np.zeros(np.shape(self.x))
        if cooling:
            for beam in beams:
                force += beam.F(self.x, self.v)

        self.v += (force/self.mass + (0,-9.806,0))*dt

    def _MB_CDF(self, v, T):
        """ 
        Cumulative Distribution function of the Maxwell-Boltzmann speed distribution 

            Parameters
            -----
            T: temperature; K
        """
        a = np.sqrt(kb * T / self.mass)
        return erf(v / (np.sqrt(2) * a)) - np.sqrt(2 / np.pi) * v * np.exp(-v**2 / (2 * a**2)) / a
