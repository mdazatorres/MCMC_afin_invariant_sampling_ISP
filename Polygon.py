## @package Polygon
#   Documentation for this module.
#
#   @author Abel Palafox.
#   @modified by Maria L. Daza-Torres.
#
#--------------------------------------------------------------------------------
#This module is a modification of the package Polygon.py of the library: lib_v1.0.zip. 
# You can find the original package in

# https://sites.google.com/a/cimat.mx/abel_palafox/codes

#--------------------------------------------------------------------------------
try:
    from numpy import zeros, ravel, unique, roll
    from numpy.linalg import norm
except :
    print 'Error. Module numpy cannot be imported.'
    exit(1)
    
##Class for handling two dimensional polygons.
class Polygon :
    def __init__(self, x) :
        
        '''The constructor function.
        
        x (array) : edges of the polygon (in counter or clockwise order).
        
        '''
        self.x = x
        self.xmax = max(x[0,:])
        self.ymax = max(x[1,:])
        self.xmin = min(x[0,:])
        self.ymin = min(x[1,:])
        self.n = x[0,:].size
    
    def copy(self) :
        
        '''Function to create a copy of the polygon.
        
            Return :
                
                P (Polygon) : a new polygon object.
        '''	
        #copy the data
        x = self.x.copy()
        return Polygon(x)
    

    def is_in_polygon(self, xp) :
        
        '''Function to determine if a point is inside of the polygon
            using crossing test method.
            
            xp (array) : test point.
            
            Return :
                
                flag (bool) : True whenever xp is inside of the polygon.
        '''     
        if (xp[0] > self.xmax) | (xp[0] < self.xmin) :
            return False
        if (xp[1] > self.ymax) | (xp[1] < self.ymin) :
            return False
        crossing_counter = 0
        
        for i in range(self.n) :
            x0 = self.x[:,i-1]-xp  
            x1 = self.x[:,i]-xp  
            
            if x0[1]*x1[1] < 0 :
                if (x0[0] > 0) & (x1[0] > 0) :
                    crossing_counter += 1
                    continue
                if (x0[0] < 0) & (x1[0] < 0) :
                    continue
                m = (x1[1] - x0[1]) / (x1[0] - x0[0])
                b = x1[1] - m*x1[0]
                x = -b / m
                
                if x > 0 :
                    crossing_counter += 1
        if crossing_counter % 2 == 0 :
            return False
        return True

        
    def is_valid_polygon(self) :
       
        '''Function to search for self-intersections.
        
            Return :
                
                flag (bool) : False whenever two sides intersect.
        '''
        def cross_product(u,v) :
            
            return u[0]*v[1] - u[1]*v[0]
        
        def intersects(line1, line2) :##To check for intersection of two sides.
            p = line1[0]
            r = [line1[1][0]-p[0],line1[1][1]-p[1]]
            q = line2[0]
            s = [line2[1][0]-q[0],line2[1][1]-q[1]]
            rs = cross_product(r,s)
            
            if rs !=0 :
                t = cross_product([q[0]-p[0],q[1]-p[1]],s) / rs
                u = cross_product([q[0]-p[0],q[1]-p[1]],r) / rs
                
                if (0.0 <= t <= 1.0) and (0.0 <= u <= 1.0) :
                    return True
            else :
                
                if cross_product([q[0]-p[0],q[1]-p[1]],r) == 0.0 :
                    return True
            return False
        
        #compare every pair of sides.
        for i in range(2,self.n-1) :
            
            line1 = [[self.x[0,i],self.x[1,i]],[self.x[0,i+1],self.x[1,i+1]]]
            for j in range(0,i-1) :				
                line2 = [[self.x[0,j],self.x[1,j]],[self.x[0,j+1],self.x[1,j+1]]]
                if intersects(line1,line2) :
                    return False
        line1 = [[self.x[0,-1],self.x[1,-1]],[self.x[0,0],self.x[1,0]]]
        
        for j in range(1,self.n-2) :
            line2 = [[self.x[0,j],self.x[1,j]],[self.x[0,j+1],self.x[1,j+1]]]
            if intersects(line1,line2) :
                return False
        return True
    
    
    def is_valid_polygon_hull(self) :
        
        '''Function to search for self-intersections for alpha-hull.
		
		Return :
            
			flag (bool) : False whenever two sides intersect.
		'''
        p=ravel(self.x, order='F').reshape(len(self.x[0,:]),2)
        uniq = unique(p.view(p.dtype.descr*p.shape[1]))
        if len(uniq)==len(p):
            return True
        else:
            return False
    
    def area(self):
        x_ = roll(self.x,-1,axis=1)
        return 0.5*(self.x[0,:]*x_[1,:] - x_[0,:]*self.x[1,:]).sum()


		
def polygon_from_hull(hull, points) :
    
	'''Provides a Polygon object from a unnordered sequence of segments.
		
		hull (array[[int,int],[int,int]]) : pair of segments
			with the indices of the points that form the new polygon.
		points (array[float]) : Point cloud.
	
	Return :

		P (Polygon) : a polygon.
	
	'''
	L = hull.tolist()
	
	n = len(L)
	x = zeros(shape=(2,n))	
	
	x[0,0] = points[0,L[0][0]]
	x[1,0] = points[1,L[0][0]]

	target = L[0][1]
	L.pop(0)
	
	counter = 1
	while len(L) > 0 :
		[i,j] = search_target(L,target)
		if i==-1 :
			raise ValueError('target not found.')
		x[0,counter] = points[0,L[i][j]]
		x[1,counter] = points[1,L[i][j]]

		target = L[i][1-j]
		
		L.pop(i)
		
		counter += 1
	
	return Polygon(x)


def search_target(L, target) :
    
	'''Function to locate the segment that has one vertex equals to target.
	
		L (array[int]) : array of segments.
		target (int) : index of the target vertex.
	
	Return :

		s ([int,]) : a segment where first entry is the 
			complement of the located segment [i,target]
	
	'''
	for i in range(len(L)) :
		if L[i][0] == target :
			return [i,0]
		if L[i][1] == target :
			return [i,1]

	return [-1,-1]


def is_valid_spline_JA(P,h):
    
    '''Function to determine if a closed curve has peaks width smaller than the mesh size.
    
        P (array[]) : Parametrization of a closed curve.
        h(float): size of the mesh.
        
        Return :
            
            flag (bool): True whenever there are no peaks in the closed curve.    
		'''
    
    P=P[:,0:-1]
    n=P.shape[1]
    for i in range(n):
        ll=0
        j=1
        while ll<=h:
            ll=norm(P[:,i]-P[:,(i+j)%n])
            j=j+1        
        inder=(i+j)%n
        ll=0
        j=1
        while ll<=h:
            ll=norm(P[:,i]-P[:,(i-j)%n])
            j=j+1
        inizq=(i-j)%n    
        if inizq > inder:
            for j in range(inder,inizq):
               if norm(P[:,i]-P[:,j])<h:
                   return False
        else: 
            for j in range(inder,n):
               if norm(P[:,i]-P[:,j])<h:
                   return False
              
    return True
