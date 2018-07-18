## @package ScatteringWave_k1
#   Documentation for this module.
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


try :
	from numpy import array, log, pi, zeros, sqrt, exp, dot, identity, ceil,floor,eye
except :
	print 'Error. Module numpy cannot be imported.'
	exit(1)

try :
    from scipy import special
    from scipy import  linalg as lg
except :
	print 'Error. Module scipy.special or scipy.linalg cannot be imported.'
	exit(1)

try :
    from shapely.geometry import Polygon, Point
    
except :
	print 'Error. Module shapely.geometry cannot be imported.'
	exit(1)



def grid(N):
    '''Computes a square grid.
		N (int) : number of points in a side of the square domain. 
	Return :
		A(list(array([int,int]))) : NxN matrix with the grid points
	'''
    point_grid=[]
    for i in range (N):
        for j in range(N):
            point_grid.append(array([j,i]))
    return point_grid
   

   
def g(r):
    '''
    Computes the Hankel function of the first kind of order zero.
    
    '''
    return -(1j/4.0)*special.hankel1(0, k*r)  


def matrix_hankel(N,k,h):
    
    ''' Computes a matrix with the evaluation of the Hankel function in each
        point of the grid.
        N (int) : number of points in a side of the square domain. 
        k(int)  : numberwave
        h(float): mesh size
        
	Return :
        
		A(matrix(double,double)) :  NxN matrix 
        
    '''
    size=N*N
    beta1=-(1j/4.0)+(1/(2*pi))*(log(h*k/2.0)+ec+coef)
    A=zeros((size,size),dtype=complex)

    for i in range(size):
        A[i][i]=beta1
        for j in range(size-(i+1)):
            j=i+j+1
            a=point_grid[i]-point_grid[j]
            A[i][j]=g(lg.norm(a)*h)
            A[j][i]=A[i][j]
                            
    return A

m=4 #number of directions
def u_i(): ## incident wave
    ''' Computes the incident field for four incident directions.'''
    size=N*N
    
    d=[[sqrt(2)/2,sqrt(2)/2],[-sqrt(2)/2,sqrt(2)/2],[-sqrt(2)/2,-sqrt(2)/2],[sqrt(2)/2,-sqrt(2)/2]]
    Ui=zeros((size,m),dtype=complex)
    for j in range(m):    
        for i in range(size):
            vec_direction=array(d[j])
            Ui.T[j][i]=exp(1j*k*dot(vec_direction.T,h*point_grid[i]))
          
    return Ui #,d

N=40
point_grid=grid(N) 
size=N*N

k=5
ec=0.5772156649015328606    
h=0.02
h2,k2=h**2,k**2
k2h2=k2*h2
coef=-1.3105329259115095

A=matrix_hankel(N,k,h)
I=identity(N*N)
ui=u_i()

def sol_ec_lippman_schwinger(spline,elast):
    '''Evaluates the forward map using the corrected trapezoidal rule
       to approximate the Lippman-Schwinger equation.
       Based on J. Aguilar and Y.Chen 2004.

		spline (array([double,double])) : parameterized scatterer boundary
        elast(double): value of the refractive index
	
	Return :
	
		us(array) : the near-field evaluated at grid points in 
                     four incident directions.			
	'''
    x = spline
    max_x=max(x[0,:])
    min_x=min(x[0,:])
    max_y=max(x[1,:])
    min_y=min(x[1,:])
    poly=Polygon(x.T)
    alpha=k2h2*elast
    size=N*N
    sd=[]

    mx=ceil(min_x/h)
    Mx=floor(max_x/h)
    my=ceil(min_y/h)
    My=floor(max_y/h)
    
    k,s= int(Mx-mx+1), int(My-my+1)
    nd=range(0,int(my*N))
    for j in range(s):
        a1=range(int((my+j)*N),int((my+j)*N+mx))
        a2=range(int((my+j)*N+mx+k),int(N*((my+j)+1)))
        nd=nd+a1+a2
        for i in range(k):
            position=int((my+j)*N+mx+i)
            mesh_point=h*point_grid[position]
            point=Point(mesh_point)
            p1=poly.contains(point)
            if p1==True:
                sd.append(position)
            else:
                nd.append(position)
                
    nd1=range(int((My+1)*N),int(N*N))
    nd=nd+nd1
    
    C22=alpha*A[sd,:][:, sd]+ eye(len(sd),dtype=complex)
    C12=alpha*A[nd,:][:, sd]
    U2=lg.solve(C22,ui[sd,:])
    U1=ui[nd,:]-dot(C12,U2)
    U=zeros((size,m),dtype=complex)
    U[nd]=U1
    U[sd]=U2
             
    return U.T


