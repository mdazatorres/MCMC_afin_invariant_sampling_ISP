## @package pymh
#   Documentation for this module
#
#   @author Abel Palafox 
#   @modified by Maria L. Daza-Torres
# 
#--------------------------------------------------------------------------------
#This module is a modification of the package pymh.py of the library: lib_v1.0.zip. 
# You can find the original package in

# https://sites.google.com/a/cimat.mx/abel_palafox/codes

#--------------------------------------------------------------------------------

try :
	from pylab import plot, figure, show
except :
	print 'Error. Module pylab cannot be imported.'
	exit(1)

try: 
	from numpy import ones, shape, mat, cov, ceil, matrix, sqrt, random, exp, zeros, array, floor, log 
except :
	print 'Error. Module numpy cannot be imported.'
	exit(1)

try :
	import copy
except :
	print 'Error. Module copy cannot be imported.'
	exit(1)

try: 
	from functools import partial 
except :
	print 'Error. Module functools cannot be imported.'
	exit(1)
    



class mh :
    
    '''
	Class to implement a standar Metropolis-Hasting algorithm.

	This implementation considers that the logarithm
	of the likelihood times the prior is evaluated as energy function. 

	The corresponding transition kernel is considered as a weighted mix of kernels:
		q(x,y) = \sum_i^n w_i q_i(x,y).
	That is, a proposal is obtained from simulating from kernel
	q_i with probability w_i.

	User must provide energy function, support for the parameters,
	kernels q_i, weights w_i and corresponding functions to simulate 
	from kernels.
	'''
    
    def __init__(self, q_f, s_f, w, energy, support, args, status=0) :
        
        '''The constructor function
			
			q_f (list[functions]) : list of functions to evaluate 
				transition kernel q(x,y).
			s_f (list[functions]) : list of functions to simulate
				a proposal from the transition kernel,
					y <- q(x,*).
			w (list[float]) : list of probability weights to 
				choose between transition kernels.
			energy (function) : function to evaluate the 
				logarithm of the likelihood times the prior.
			support (function) : boolean function to determine whether
				any parameter is out of support.
			args (list[params]) : list of the parameters that are 
				fixed within the energy function.
			status (int) : number of iterations to display partial 
				status of the chain. Default is zero for 
				not displaying.
			
		'''
        self.q = q_f
        self.sim_q = s_f
        self.w = w
        self.energy = partial(energy)
        self.support = support
        self.acc_rate = 0
        self.T = 0
        self.status = status
    
    def Run(self, x0, max_sims) : 
        
        '''The run function begins with the MH algorithm.

			x0 (list or array) : initial guess.
			max_sims (int) : maximum number of iterations of 
				the algorithm.

		This method fills the list *self.sims*, with the 
		simulated (accepted proposals) along the *max_sims* 
		iterations. Corresponding values of energy
		are stored in the list *self.U*.
		
		Acceptance rates for each kernel are also computed.	
		'''
        
        #accepted simulations
        self.Total_acc_rate=[]
        self.Acc_rate=[]
        self.sims = []
        #energy values 
        self.U = []
        #for reporting acceptance rates for every kernel
        self.acc_rate = zeros(len(self.sim_q))
        self.acc_rate_chain = []
        #to counting the number of rejected proposals 
        self.out_support = 0
        self.rejected = 0
        self.T = max_sims
        #energy value of the MAP
        self.umap = 0.0
        #store the initial guess
        #initial guess is assumed within the support
        self.sims.append(x0)
        u = self.energy(x0)
        self.U.append(u)
        # to consider list or numpy arrays
        if type(x0) == list :
            self.MAP = copy.deepcopy(x0)
        else :
            self.MAP = x0.copy()
        
        self.umap = u
        counter = 1
            
        while counter < max_sims :	
            #the simulation on t-1
            theta_0 = self.sims[-1]
            #choose the kernel to simulate from
            i = choose_q(self.w) 
            #simulate from i-eth kernel
            theta = self.sim_q[i](theta_0)

            #avoid out-supported proposals
            if not self.support(theta) :
                self.out_support += 1
                continue	

            #evaluate the proposal
            ut = self.energy(theta)
            #metropolis-hastings criterion
           
            a = self.U[-1] - ut + self.q[i](theta_0) - self.q[i](theta)
            u = random.uniform()	
            if u < exp(a) :
                #accept the proposal
                self.sims.append(theta)
                self.U.append(ut)
                self.acc_rate[i] += 1
                self.acc_rate_chain.append(i)
                #update the MAP
                if ut < self.umap :
                    if type(theta) == list :
                        self.MAP = copy.deepcopy(theta)
                    else :
                        self.MAP = theta.copy()
                    
                    self.umap = ut
            else :
                #reject
                self.sims.append(theta_0)
                self.U.append(self.U[-1])
                self.rejected += 1	
            
            counter += 1	
            if (self.status != 0) :
                if (counter % self.status == 0) :
                    self.Acc_rate.append(1.0*self.acc_rate/counter)
                    self.Total_acc_rate.append(1.0*self.acc_rate/(counter+self.rejected))
                    print 'it: ', counter, '  ,  energy: ', self.U[-1]
                    print '   acceptance rates: ', 1.0*self.acc_rate/counter
                    print '   total acceptance rates: ', 1.0*self.acc_rate/(counter+self.rejected)
                    print '   out of support: ', self.out_support, ' , rejected: ', self.rejected		
                    #self.acc_rate /= 1.0*max_sims
                                
        print 'Process terminated.'
        print '   MAP energy: ', self.umap, '  , P(MAP): ', exp(-self.umap)
        print '   acceptance rates: ', 1.0*self.acc_rate/counter
        print '   total acceptance rates: ', 1.0*self.acc_rate/(counter+self.rejected)
        print '   out of support: ', self.out_support, ' , rejected: ', self.rejected
        print '   total number of evaluations: ', (counter+self.rejected)
    
    def Acc(self) :
        
        '''Returns a list([int]) with the number 
		of times that every kernels was 
		accepted along the chain. Also, the 
		acceptance rate for each kernel with
		repect to the accepted proposal, and
		with respect to the total proposals
		(accepted and rejected) are reported.
		'''
        
        print '   acceptance rates: ', 1.0*self.acc_rate/self.T
        print '   total acceptance rates: ', 1.0*self.acc_rate/(self.T+self.rejected)
        return self.acc_rate_chain	

    def Out_Support(self) :
        
        '''Return the number of proposals
		rejected due to being out of support.
		'''
        
        return self.out_support

    def Walk(self, start=0, end=-1, step=1, plot_flag=True) :
        
        '''Returns the energy values for simulations
		along the chain from *start* to *end* with step 
		size *step*. 

			start (int) : initial iteration
			end (int) : final iteration
			step (int) : subsampling factor
			plot_flag (bool) : For plotting the chain. Default is True.
		'''
        
        if plot_flag :
            figure()
            plot(self.U[start:end:step])
            show()
        return self.U[start:end:step]	
    
    def IAT(self, par=-1, maxlag=0, start=0, end=0) :	
        
        '''Computes the Integrated Autocorrelation Time (IAT).
		This implementation is based on the code of 
		Christen, et al. 2010.  "A general purpose 
		sampling algorithm for continuous distributions 
		(the t-walk)".
		
			par (int) : parameter(s) number for which IAT is computed.
				Default is -1, that is all the parameters.
			maxlag (int) : Maximum lag allowed. Default zero is for 
				determining lag from the number of samples.
			start (int) : Initial simulation. Default is zero, the 
				first sample of the chain.
			end (int) : Final simulation. Default is zero for considering
				the entire chain.
			
		Returns:
			IAT (float) : IAT  factor.
			
		'''
        
        if (end == 0):
            end = self.T
            #To arrange the data in the format of Christen's code.
        Output = []
        counter = 0
        for s in self.sims :
            item = []
            for si in s :
                item.append(si)
            item.append(self.U[counter])
            Output.append(item)
            counter += 1
            #Running the Christen's code.
        iat = IAT( array(Output), cols=par, maxlag=maxlag, start=start, end=end)
        return iat
		
	
## Choose the proposal to simulate from
#   @param w list of weights
#   @return i the index of the selected proposal q_i
def choose_q(w) :
	
	n = w.size
	u = random.uniform()
	
	for i in range(n) :
		u -= w[i]
		if u <= 0 :
			return i
			
	return n-1


############################################################################################
#The following code corresponds to Christen's code for IAT computing. For 
#details see Christen et, al. 2010 and references therein

#### Auxiliary functions to calculate Integrated autocorrelation times of a time series 

####  Calculates an autocovariance 2x2 matrix at lag l in column c of matrix Ser with T rows
####  The variances of each series are in the diagonal and the (auto)covariance in the off diag.
def AutoCov( Ser, c, la, T=0):
	if (T == 0):
		T = shape(Ser)[0]  ### Number of rows in the matrix (sample size)

	return cov( Ser[0:(T-1-la), c], Ser[la:(T-1), c], bias=1)
	
	
	

## Calculates the autocorrelation from lag 0 to lag la of columns cols (list)
# for matrix Ser
def AutoCorr( Ser, cols=0, la=1):
	T = shape(Ser)[0]  ### Number of rows in the matrix (sample size)

	ncols = shape(mat(cols))[1] ## Number of columns to analyse (parameters)

	#if ncols == 1:
	#	cols = [cols]
		
	### Matrix to hold output
	Out = matrix(ones((la+1)*ncols)).reshape( la+1, ncols)
		
	for c in range(ncols):
		for l in range( 1, la+1):  
			Co = AutoCov( Ser, cols[c], l, T) 
			Out[l,c] = Co[0,1]/(sqrt(Co[0,0]*Co[1,1]))
	
	return Out
	

### Makes an upper band matrix of ones, to add the autocorrelation matrix
### gamma = auto[2*m+1,c]+auto[2*m+2,c] etc. 
### MakeSumMat(lag) * AutoCorr( Ser, cols=c, la=lag) to make the gamma matrix
def MakeSumMat(lag):
	rows = (lag)/2   ### Integer division!
	Out = mat(zeros([rows,lag], dtype=int))
	
	for i in range(rows): 
		Out[i,2*i] = 1
		Out[i,2*i+1] = 1
	
	return Out


## Finds the cutting time, when the gammas become negative
def Cutts(Gamma):
	cols = shape(Gamma)[1]
	rows = shape(Gamma)[0]
	Out = mat(zeros([1,cols], dtype=int))
	Stop = mat(zeros([1,cols], dtype=bool))
	
	if (rows == 1):
		return Out
		
	i = 0
	###while (not(all(Stop)) & (i < (rows-1))):
	for i in range(rows-1):
		for j in range(cols):  ### while Gamma stays positive and decreasing
			if (((Gamma[i+1,j] > 0.0) & (Gamma[i+1,j] < Gamma[i,j])) & (not Stop[0,j])):
				Out[0,j] = i+1 ## the cutting time for colomn j is i+i
			else:
				Stop[0,j] = True
		i += 1
	
	
	return Out


####  Automatically find a maxlag for IAT calculations
def AutoMaxlag( Ser, c, rholimit=0.05, maxmaxlag=20000):
	Co = AutoCov( Ser, c, la=1)
	rho = Co[0,1]/Co[0,0]  ### lag one autocorrelation
	
	### if autocorrelation is like exp(- lag/lam) then, for lag = 1
	lam = -1.0/log(abs(rho)) 
	
	### Our initial guess for maxlag is 1.5 times lam (eg. three times the mean life)
	maxlag = int(floor(3.0*lam))+1
	
	### We take 1% of lam to jump forward and look for the
	### rholimit threshold
	jmp = int(ceil(0.01*lam)) + 1
	
	T = shape(Ser)[0]  ### Number of rows in the matrix (sample size)

	while ((abs(rho) > rholimit) & (maxlag < min(T/2,maxmaxlag))):
		Co = AutoCov( Ser, c, la=maxlag)
		rho = Co[0,1]/Co[0,0]
		maxlag = maxlag + jmp
		###print("maxlag=", maxlag, "rho", abs(rho), "\n")
		
	maxlag = int(floor(1.3*maxlag));  #30% more
	
	if (maxlag >= min(T/2,maxmaxlag)): ###not enough data
		fixmaxlag = min(min( T/2, maxlag), maxmaxlag)
		print "AutoMaxlag: Warning: maxlag= %d > min(T/2,maxmaxlag=%d), fixing it to %d" % (maxlag, maxmaxlag, fixmaxlag)
		return fixmaxlag
	
	if (maxlag <= 1):
		fixmaxlag = 10
		print "AutoMaxlag: Warning: maxlag= %d ?!, fixing it to %d" % (maxlag, fixmaxlag)
		return fixmaxlag
		
	print "AutoMaxlag: maxlag= %d." % maxlag
	return maxlag
	
	
### Find the IAT
def IAT( Ser, cols=-1,  maxlag=0, start=0, end=0):

	ncols = shape(mat(cols))[1] ## Number of columns to analyse (parameters)
	if ncols == 1:
		if (cols == -1):
			cols = shape(Ser)[1]-1 ### default = last column
		cols = [cols]
	
	if (end == 0):
		end = shape(Ser)[0]

	if (maxlag == 0):
		for c in cols:
			maxlag = max(maxlag, AutoMaxlag( Ser[start:end,:], c))

	#print("IAT: Maxlag=", maxlag)

	#Ga = MakeSumMat(maxlag) * AutoCorr( Ser[start:end,:], cols=cols, la=maxlag)
	
	Ga = mat(zeros((maxlag/2,ncols)))
	auto = AutoCorr( Ser[start:end,:], cols=cols, la=maxlag)
	
	### Instead of producing the maxlag/2 X maxlag MakeSumMat matrix, we calculate the gammas like this
	for c in range(ncols):
		for i in range(maxlag/2):
			Ga[i,c] = auto[2*i,c]+auto[2*i+1,c]
	
	cut = Cutts(Ga)
	nrows = shape(Ga)[0]
		
	ncols = shape(cut)[1]
	Out = -1.0*mat(ones( [1,ncols] ))
	
	if any((cut+1) == nrows):
		print("IAT: Warning: Not enough lag to calculate IAT")
	
	for c in range(ncols):
		for i in range(cut[0,c]+1):
			Out[0,c] += 2*Ga[i,c]
	
	return Out



############################################################################################


