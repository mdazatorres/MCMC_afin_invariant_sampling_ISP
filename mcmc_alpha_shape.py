
## Example 1 of the article:
##  A computational geometry method for the inverse scattering problem.
#   
#---------------------------------------------------------------------
#    Copyright (C) 2018 
#    @author Maria L. Daza-Torres

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Contact Information:
# Email: mdazatorres@cimat.mx

#---------------------------------------------------------------------

from numpy import loadtxt, save,random,log
from numpy.linalg import norm
from scipy.stats import norm as gaussian
from Polygon import Polygon, is_valid_spline_JA
from scipy.spatial import ConvexHull
import ScatteringWave_k1, ScatteringWave_k5
import rwcp
import time, os


def init(n) : #initial guess for the point cloud 
	x = random.uniform(0.5,0.6,size=(2,n))
	x[0,:] += random.uniform(-0.2,0.2)
	x[1,:] += random.uniform(-0.2,0.2)
	return x


def energy (p):
    ''' compute of energy function'''
    m=4
    [x,spline,alpha,elast] = p 
    U1 = ScatteringWave_k1.sol_ec_lippman_schwinger(spline,elast)
    U5 = ScatteringWave_k5.sol_ec_lippman_schwinger(spline,elast)  
    

    log_likehood_k1=0
    log_likehood_k5=0
    
    for i in range(m):
        log_likehood_k1+=norm(d_1[i]-U1[i])**2
    
    for i in range(m):
        log_likehood_k5+=norm(d_5[i]-U5[i])**2    
    
    k,THETA=1.66,30
    log_likehood=log_likehood_k1+log_likehood_k5
    log_prior= -k*log(THETA)+(k-1)*log(elast)-(elast/THETA)
    
    return log_likehood-log_prior

    
def support(p) : # false if not satisfies some the next conditions on the parameter spaces
    [x,spline,alpha,elast] = p
    
    if elast>60: 
        return False
        
    if len(spline) == 0 : 
        return False
        
    if (alpha <= 0) : 
        return False
    x_sp = spline
    
    if (x_sp > 0.78).any() | (x_sp < 0.0).any():
        return False
    Q = Polygon(x_sp)
    
    if Q.is_valid_polygon() == False :
        return False
    
    area_convex=ConvexHull(x.T).volume
    if (Q.area()/area_convex) < 0.8: 
        return False

    return is_valid_spline_JA(x_sp,h)     

start=time.clock()
npts=64 # number of interpolation points
it_max = 5 #maximum of mcmc iterations
n_p_cloud = 20 #points number on the cloud
sigma = 1e-2 #noise level



h=0.02 #size of step
n_data=1600 # number of data

#data for  k=1, k=5
d1=loadtxt('data_kite_ex1/near_field_data_kite_k1_40x40_dir1.txt',dtype=complex, converters={0: lambda s: complex(s.decode().replace('+-', '-'))})+gaussian.rvs(loc=0.0,scale=sigma,size=n_data)+1.0j*gaussian.rvs(loc=0.0,scale=sigma,size=n_data)
d2=loadtxt('data_kite_ex1/near_field_data_kite_k5_40x40_dir2.txt',dtype=complex, converters={0: lambda s: complex(s.decode().replace('+-', '-'))})+gaussian.rvs(loc=0.0,scale=sigma,size=n_data)+1.0j*gaussian.rvs(loc=0.0,scale=sigma,size=n_data)
d3=loadtxt('data_kite_ex1/near_field_data_kite_k1_40x40_dir3.txt',dtype=complex, converters={0: lambda s: complex(s.decode().replace('+-', '-'))})+gaussian.rvs(loc=0.0,scale=sigma,size=n_data)+1.0j*gaussian.rvs(loc=0.0,scale=sigma,size=n_data)
d4=loadtxt('data_kite_ex1/near_field_data_kite_k5_40x40_dir4.txt',dtype=complex, converters={0: lambda s: complex(s.decode().replace('+-', '-'))})+gaussian.rvs(loc=0.0,scale=sigma,size=n_data)+1.0j*gaussian.rvs(loc=0.0,scale=sigma,size=n_data)
d5=loadtxt('data_kite_ex1/near_field_data_kite_k1_40x40_dir5.txt',dtype=complex, converters={0: lambda s: complex(s.decode().replace('+-', '-'))})+gaussian.rvs(loc=0.0,scale=sigma,size=n_data)+1.0j*gaussian.rvs(loc=0.0,scale=sigma,size=n_data)
d6=loadtxt('data_kite_ex1/near_field_data_kite_k5_40x40_dir6.txt',dtype=complex, converters={0: lambda s: complex(s.decode().replace('+-', '-'))})+gaussian.rvs(loc=0.0,scale=sigma,size=n_data)+1.0j*gaussian.rvs(loc=0.0,scale=sigma,size=n_data)
d7=loadtxt('data_kite_ex1/near_field_data_kite_k1_40x40_dir7.txt',dtype=complex, converters={0: lambda s: complex(s.decode().replace('+-', '-'))})+gaussian.rvs(loc=0.0,scale=sigma,size=n_data)+1.0j*gaussian.rvs(loc=0.0,scale=sigma,size=n_data)
d8=loadtxt('data_kite_ex1/near_field_data_kite_k5_40x40_dir8.txt',dtype=complex, converters={0: lambda s: complex(s.decode().replace('+-', '-'))})+gaussian.rvs(loc=0.0,scale=sigma,size=n_data)+1.0j*gaussian.rvs(loc=0.0,scale=sigma,size=n_data)

d_1=[d1,d3,d5,d7]
d_5=[d2,d4,d6,d8]
n_d=len(d_1)+len(d_5)

##############################################################

path = 'result_kite'+'_'+'sim_'+str(it_max)+'_'+'num_dir='+str(n_d)+'_'+'h='+str(h)
try :
	os.mkdir(path)
except:
	#to avoid overwritte previous results
	print 'Existing folder. overwritte? (y/n)'
	pass
	op = input('op = ')
	if op == 'y' :
		pass
	else :
		raise
		
path += '/'

###############################################################
#initial guess for Metropolis-Hastings
init_guess_flag = False
while init_guess_flag == False :
    points = init(n_p_cloud)
    poly_temp = Polygon(points)
    alpha = 1e-6
    gamma = rwcp.q_spline(npts,points,1.0/alpha)
    elast= random.uniform(5,40)
    if len(gamma) == 0 :
        init_guess_flag = False
        continue
    if support([points,gamma,1.0/alpha,elast]) :
        init_guess_flag = True 
        init_guess = [points,gamma,1.0/alpha,elast]

#arguments that will be fixed along forward map evaluations

rwmh = rwcp.rwcp(energy, support,npts, 500)
rwmh.Run(init_guess, it_max)
burnin = it_max/2

U = rwmh.Walk(start=burnin)


##Output Data
save(path+'energy_total_ex1',rwmh.U)
save(path+'theta_ex1',rwmh.sims)


