This module consists in the following files:

- pymh.py

	Implements a standar Metropolis-Hastings algorithm.
	Transition kernel is assumed as a mix of kernels.
	User must define each kernel by providing the 
	corresponding simulating functions, the kernel transition evaluation
	functions, and the weights on the mix.
	Also, user must provide the energy function as the 
	logarithm of likelihood times the prior, and the support 
	of the parameters.

- rwcp.py

	Implements a Random Walk Metropolis Hastings over a set 
	of points in 2D and a parameter that represent the refractive
	index. This code uses the pymh code, and it defines the 
	transition kernel with three affine invariant moves and
	an additional proposal to the refractive index. For a 
	sound obstacle, this implementation corresponds to the 
	method described in [Palafox et al. 2016]. 
	
	In this implementation, a proposal with the affine invariant
	property for the alpha shape parameter is added and the prior 
	distribution of the refractive index is used as its proposal.


- alpha_shape.py

	Implementation of the alpha-shape algorithm of 
	Edelsbrunner.

- Polygon.py

	For handling Polygonal objects suitablely for rcwp.

- ScatteringWave_k1.py
	
	This code implements the near-field computing with the corrected
	trapezoidal ruler to solve the Lippman-Schwinger equation. 

	This method is described in [Aguilar et al., 2014]. In this case, 
	is evaluated for numberwave equal to one and the incident field 
	is taken in the directions d_i, i=1,..,4.,(see article).


- ScatteringWave_k5.py

	In this code evaluates ScatteringWave_k1.py for a numberwave equal 
	to five, and the incident field is taken in the directions 
	d_2i, i=1,..,4., (see article).


For the example 1 of the article: A computational geometry method for the 
inverse scattering problem, see mcmc_alpha_shape.py

The data used for this example are computed using a method introduced 
in Vainikko 1997, and they are in the folder: data_kite_ex1.



