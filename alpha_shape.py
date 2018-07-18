## @package alpha_shape
#   Documentation for this module
#
# @author Abel Palafox 
# @modified by Maria L. Daza-Torres
#
#--------------------------------------------------------------------------------
#This module is a modification of the package alpha_shape.py of the library: lib_v1.0.zip. 
# You can find the original package in

# https://sites.google.com/a/cimat.mx/abel_palafox/codes

#--------------------------------------------------------------------------------
import math
import Polygon
try:
	from numpy import array
except :
	print 'Error. Module numpy cannot be imported.'
	exit(1)

try :
    from scipy.spatial import Delaunay
except :
	print 'Error. Module scipy.spatial cannot be imported.'
	exit(1)



#we calculate the mean of the circum radius
def alpha_radius(p):
    
    '''For a set of points in 2D, this function computes the maximum and minimum radius
       of the circumcircles associated to each triangle of the Delaunay triangulation.
       
        p (array) : point cloud.
        
        Return :
            
            hull (array[double,double]) : 2D-array with the minimum and maximum radius. 
    '''
    
    points = p.T
    tri = Delaunay(points)
    circum_R=[]

    for ia, ib, ic in tri.vertices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        # Lengths of sides of triangle
        
        a = math.sqrt((pa[0]-pb[0])**2 + (pa[1]-pb[1])**2)
        b = math.sqrt((pb[0]-pc[0])**2 + (pb[1]-pc[1])**2)
        c = math.sqrt((pc[0]-pa[0])**2 + (pc[1]-pa[1])**2)
        # Semiperimeter of triangle
        s = (a + b + c)/2.0
        # Area of triangle by Heron's formula
        area = math.sqrt(s*(s-a)*(s-b)*(s-c))
        #circum_r = a*b*c/(4.0*area) asi estaba
        if area == 0 :
            circum_r = 0
        else:
            circum_r = a*b*c/(4.0*area)
            
        circum_R.append(circum_r)
        
    max_circum_R=array(circum_R).max()
    min_circum_R=array(circum_R).min()

        
    return [min_circum_R, max_circum_R]

def alpha_shape(p,alpha) :
    
    '''Function to compute the alpha_shape of a set of points, for a
    given alpha value.
        
        p (array) : point cloud.
        alpha (float) : alpha value.	
    
    Return :
        
        hull (array[int]) : array with the indices of
        the points that corresponds to the edges of the alpha shape.	
    '''
    
    #computes the Delaunay triangulation
    points = p.T
    tri = Delaunay(points)
    edges_list = []
    for ia, ib, ic in tri.vertices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        
        # Lengths of sides of triangle
        a = math.sqrt((pa[0]-pb[0])**2 + (pa[1]-pb[1])**2)
        b = math.sqrt((pb[0]-pc[0])**2 + (pb[1]-pc[1])**2)
        c = math.sqrt((pc[0]-pa[0])**2 + (pc[1]-pa[1])**2)
        
        # Semiperimeter of triangle
        s = (a + b + c)/2.0
        
        # Area of triangle by Heron's formula
        area = math.sqrt(s*(s-a)*(s-b)*(s-c))

        if area == 0 :
            circum_r = 0
        else:
            circum_r = a*b*c/(4.0*area)
            
        # Here's the radius filter.
        if circum_r <= alpha :#ALERT
#        if circum_r < alpha :#ALERT
            edges_list.append([min(ia,ib),max(ia,ib)])
            edges_list.append([min(ib,ic),max(ib,ic)])
            edges_list.append([min(ic,ia),max(ic,ia)])

    #fill hull list with the indices of the edges without 
    #repeated elements.
    hull = []
    while len(edges_list) > 0 :
        item = edges_list.pop(0)
        
        i = 0
        counter = 0
        while i<len(edges_list) :
            if edges_list[i] == item :
                edges_list.pop(i)
                counter += 1
                continue
            i += 1
        if counter == 0 :
            hull.append(item)
                
    return array(hull)


def alpha_shape_connected_alpha(p,alpha) :
    
	'''Function to compute a valid alpha shape of 
		a point cloud for a value of alpha.
		That is, non-empty, non-self-intersected
		and non-disconnected alpha shape.
	
		p (array) : point cloud
		alpha (float) : alpha value
	
	Return :
        
		hull (array[int]) : a valid alpha shape	
	'''
    
	points = p
	val = alpha
	
	try :
		#try the alpha shape for a given alpha
		hull = alpha_shape(points,val)
		#check for non-empty hull
		if len(hull) == 0 :
			raise ValueError('Empty hull')

		#check for non-disconnedted hull
		P = Polygon.polygon_from_hull(hull,points)
		for i in range(points.shape[1]) :
			p = points[:,i]
			if P.is_in_polygon(p) == False :  
				edge = False
				for j in P.x.T :
					if (j == p).all() :
						edge = True
				if not edge :
					raise ValueError('Points outside')

		#check for non-self-intersected hull
            
		if P.is_valid_polygon_hull() == False : 
			raise ValueError('Invalid Polygon')

	except ValueError: #as err:
		return [[],-1]

	return [hull,0]



