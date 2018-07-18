## @package rwcp
#   Documentation for this module.
#
#   @author Abel Palafox.
#   @modified by Maria L. Daza-Torres.
#
#--------------------------------------------------------------------------------
#This module is a modification of the package rwcp.py of the library: lib_v1.0.zip. 
# You can find the original package in

# https://sites.google.com/a/cimat.mx/abel_palafox/codes

#--------------------------------------------------------------------------------

from alpha_shape import alpha_shape_connected_alpha,alpha_radius
from Polygon import Polygon, polygon_from_hull
import pymh

try :
	from numpy import random, array, linspace, pi,sqrt, log 
except :
	print 'Error. Module numpy cannot be imported.'
	exit(1)

try :
    import scipy.interpolate as si
except :
	print 'Error. Module scipy.interpolate cannot be imported.'
	exit(1)

try :
    from functools import partial
except :
	print 'Error. Module functools cannot be imported.'
	exit(1)



class rwcp(pymh.mh) :
    
    '''Class to implement a random walk Metropolis-Hastings over a set of points (Point Cloud)
       and a parameter that represent the refractive index.

       The alpha shape of the point cloud is interpolated with a B-spline to be later evaluated 
       within the forward map along with the refractive index.
       
       The refractive index is assumed constant on the scatterer support. Inside the code, 
       it is called elast. 
       
       This implementation uses the basic Metropolis-Hastings class implemented in 
       package pymh with affine invariant transition kernel for the point cloud 
       defined as in Palafox et, al. 2016, and a proposal with the same property 
       for the alpha shape parameter is added.
	'''	
    
    def __init__(self, energy, support, npts, status) :
        
        '''The constructor function.
        
			energy (function) : function to evaluate the 
				logarithm of the likelihood times the prior.
			support (function) : boolean function to determine whether
				any parameter is out of support.
			npts (int) : number of points in the point cloud
			status (int) : number of iterations to display partial 
				status of the chain. Default is zero for 
				not displaying.			
		'''
        #Transition kernel is a mix of three kernels.
        def q_eval(x) :
            
            '''Function to evaluate the transition kernel on the MH criterion for the
                proposals on the point cloud (Move Point and Translate).
			
				In this case is equals to 1 because the transition kernel is symmetric 
            '''
            
            return 1
        
        def q_elast(x):
            
            '''Function to evaluate the transition kernel on the MH criterion for the
                proposal on the refractive index.
                
                The transition kernel for the refractive index is its prior distribution.
			'''
            
            elast=x[3]
            k,THETA=1.66,30
            log_prior= -k*log(THETA)+(k-1)*log(elast)-(elast/THETA)

            return log_prior
    
        q_f = [q_eval,q_eval,q_eval,q_elast] #It is a mix kernel
        s_f = [self.point_proposal, self.translate_proposal, self.alpha_proposal,self.elast_proposal]
        
        prob=0.85
        w = array([prob*0.4,prob*0.4,prob*0.2,0.15]) # Weights to choose each kernel       
        
        self.q_spline = partial(q_spline,npts)
        
        #creates the pymh object
        pymh.mh.__init__(self, q_f, s_f, w, energy, support, npts, status)
        
    def point_proposal(self, params) :
        
        '''Kernel to move randomly a single point of the set.
            params(list[array]) : list of arrays with the state of the point cloud in t-1. 
			*params[0]* contains the point cloud.
			*params[1}* contains an interpolating B-spline.
			*params[2]* contains the alpha value for the alpha shape algorithm.
            *params[3]* contains an value for the refractive index, elast.
        
        Return: 
            
			[y,spline,alpha,elast] : sampled proposal for time *t*.

		'''
        
        x = params[0]
        alpha = params[2]
        elast= params[3]
        
        n = x.size / 2
        i = random.randint(n)
        y = x.copy() 
        #compute the mean of the pairwise distances
        dx = 0.0
        dy = 0.0
        counter = 0
        
        for j in range(n) :
            if j == i :
                continue
            for k in range(j+1,n) :
                if k != i :
                    dx += sqrt( (y[0,j]-y[0,k])**2.0)
                    dy += sqrt( (y[1,j]-y[1,k])**2.0)
                    counter += 1
        mean_d = [dx/counter,dy/counter]
        #update the point
        y[0,i] += random.uniform(-mean_d[0], mean_d[0])
        y[1,i] += random.uniform(-mean_d[1], mean_d[1])
        
        spline = self.q_spline(y,alpha)
            
        return [y,spline,alpha,elast]
    
    def translate_proposal(self, params) :
        
        '''Kernel to move randomly a the entire point cloud.
		
			params(list[array]) : list of arrays with the state of the 
				point cloud in t-1. 
				*params[0]* contains the point cloud.
				*params[1}* contains an interpolating B-spline.
				*params[2]* contains the alpha value for the alpha shape algorithm.
                *params[3]* contains an value for the refractive index, elast.
		
		Return: 
            
			[y,spline,alpha,elast] : sampled proposal for time *t*.

		'''
        
        x = params[0]
        alpha = params[2]
        elast= params[3]
        
        n = x.size / 2
        y = x.copy()
        
        max_dx = -1.0
        min_dx = 1e15
        max_dy = -1.0
        min_dy = 1e15
        #compute the pairwise distances, and the maximum and the minimum distances.
        dx = 0.0
        dy = 0.0
        counter = 0
        for j in range(n) :
            for k in range(j+1,n) :
                dx += sqrt( (y[0,j]-y[0,k])**2.0)
                dy += sqrt( (y[1,j]-y[1,k])**2.0)
                counter += 1
                dx_k = sqrt( (y[0,j]-y[0,k])**2.0)
                dy_k = sqrt( (y[1,j]-y[1,k])**2.0)
                if dx_k < min_dx :
                    min_dx = dx_k
                if dx_k > max_dx :
                    max_dx = dx_k
                if dy_k < min_dy :
                    min_dy = dy_k
                if dy_k > max_dy :
                    max_dy = dy_k
        u = random.uniform(0.5,2.0)
        
        mean_dx = u*max_dx
        mean_dy = u*max_dy
        
        vx = random.uniform(-mean_dx, mean_dx)
        vy = random.uniform(-mean_dy, mean_dy)         
        for i in range(n) :
            y[0,i] += vx
            y[1,i] += vy
        spline = self.q_spline(y,alpha)
        
        return [y,spline,alpha,elast]
        
    def alpha_proposal(self, params) :
        
        '''Kernel to move the alpha parameter.
		
			params(list[array]) : list of arrays with the state of the 
				point cloud in *t-1*. 
				*params[0] contains the point cloud.
				*params[1] contains an interpolating B-spline.
				*params[2] contains the alpha value for the alpha shape algorithm.
                *params[3] contains an value for the refractive index, elast.

		Return:

			[y,spline,alpha,elast] : sampled proposal for time *t*.
		
		'''
        
        x = params[0]
        y = x.copy()
        elast=params[3]
        
        alpha = params[2]
        alpha_R=alpha_radius(x)

        min_R=alpha_R[0]
        max_R=alpha_R[1]
        

        alpha = 0.5*alpha + 0.5*random.uniform(min_R, max_R) 

        spline = self.q_spline(y,alpha) 
        
        return [y,spline,alpha,elast]
        
        
    def elast_proposal(self, params) :
        '''Kernel to move the refractive index parameter.
        
            params(list[array]) : list of arrays with the state of the 
            point cloud in *t-1*. 
            *params[0] contains the point cloud.
            *params[1] contains an interpolating B-spline.
            *params[2] contains the alpha value for the alpha shape algorithm.
            *params[3] contains an value for the refractive index, elast.
        Return:
            [y,spline,alpha,elast] : sampled proposal for time *t*.
            '''
        x = params[0]
        spline=params[1]
        alpha=params[2]
        elast=random.gamma(shape=1.66,scale=30) 
        
        return [x,spline,alpha,elast]


def q_spline(n,p,alpha) :
	'''Function to compute an interpolating B-spline of 
		the alpha shape edges, for a given point cloud.
	
		@param n (int) : number of B-spline knots.
		@param p (array[]) : point cloud.
		@param alpha (float) : parameter for the alpha shape
			algorithm.
	
		@return [x,dx,ddx] (list[array]) : interpolating B-spline *x*,
			first derivative *dx* and second derivative *ddx*.
	
	'''	
	[hull,flag] = alpha_shape_connected_alpha(p,alpha)

	if flag == -1 :
		#The alpha shape cannot be computed with that alpha value
		return []
	
	#Provides the alpha shape as a polygon
	P = polygon_from_hull(hull,p)

	t = linspace(0.0,2.0*pi,n)
	
	#check for counter or clockwise orientation
	if P.area() < 0 :
		hull_ = P.x.T.tolist()
		hull_.reverse()
		Q = Polygon(array(hull_).T)
		x_sp = curve_spline(t,Q.x,0.0)


		return x_sp

	x_sp = curve_spline(t,P.x,0.0)
     

	return x_sp
  

def curve_spline(t_,p_,s) :
    
    	'''Computes an interpolating B-spline.
	
		t_ (array[]) : Knots of evaluation of the B-spline.
		p_ (array[]) : points to interpolate. 
		s (float) : smoothing parameter.
	
	Return :

	     x : interpolating spline.
	
	'''
	#repeat the last point
	p0 = p_.T.tolist()
	p0.append(p0[0])
	p = array(p0).T
	
	#evaluation points between 0 and 1
	t = t_/(2.0*pi)

	degree = 3
	
	tck,u = si.splprep(p,k=degree,per=True,quiet=2,s=s)
	
	#evaluate the B-spline
	x = array(si.splev(t,tck))


	return x







