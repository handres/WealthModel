# WealthModel

This is a code archive of a simulation program built by Hunter Vallejos in the Summer of 2016 at Oak Ridge National Laboratory, while a HERE intern. Hunter was supervised by James Nutaro.

The code was used to simulate an agent-based model of the wealth distribution with two initializations (1) a uniform distribution of wealth and (2) a distribution of wealth similar to that of a real economy. The results were published in the following article:

https://link.springer.com/article/10.1007/s11403-017-0200-9

The program for the uniform distribution can be found in ```c_parallel```, and the program for the real economy initial distribution can be found in ```nonuniform```. In addition, an optimized model was constructed using a dynamic programming approach, which can be found in ```c_new_fast``` -- this model is for the uniform distribution initialization.

There is another portion of the code, written by Kalyan Perumalla, which has been omitted here. Kalyan's program simulates an approximate model in CUDA.

## Usage

### nonuniform:
	This folder contains the program that executes the single-packet model
	starting from a Boltzmann Gibbs exponential-Pareto distribution.
	
	To view the results that are presented in the paper, run the following code
	from inside the nonuniform directory:
	
	python run_model.py -c -rd -f 1000000_100years_diverge -p .9
	
	This will display a plot comparing the wealth share changes in the 
	single-packet model to actual changes in U.S. data. I zoomed in the figure
	to get the picture shown in the article.
	
	To get the model projection (not plotting real U.S. wealth share changes),
	do the same command without the "-c" flag:
	
	python run_model.py -rd -f 1000000_100years_diverge -p .9
	
	Description of flags:
		***RUNNING ARGS***
		-n: change the number of individuals [default: 100000]
		-gr: change the growth rate of the entire economy [default: 0.02588 
			(growth rate of US economy 1988-2012)]
		-dw: change the delta_omega (precision of the model) [default: 0.01]
		-t:	change the total time that the model runs [default: 25]
		-b:	change the exponent "beta" in wealth power function [default: 1.25]
		
		***INITIAL DISTRIBUTION ARGS***
		
		-ec: change where distribution begins transition from exponential to 
			pareto (called kappa_l in paper) [default: 0.75]
		-pc: change where distribution stops transitioning from exponential to 
			pareto (where the distribution is all pareto) (called kappa_u in 
			paper) [default: 0.9]
		-mp: the "minimum" percentile of the Pareto distribution; this should be
			very close to "ec" because if "mp" > "ec", then the distribution 
			will have an error; the code will generate a Pareto distribution 
			with w_m not set correctly. Essentially, it sets w_m to the wealth
			of the percentile "mp". It would make sense that this should be the
			beginning of the pareto distribution, which is at "ec", where the 
			initial distribution begins transition to pareto [default: 0.9]
		-g: changes the exponent (1/'wealth temperature') of the exponential 
			distribution (exp(-gt)); it is the inverse of the wealth temperature
			[default: 2.78]
			Note: I have not implemented a way to change the scaling constant 
			"c" that is mentioned in the paper because the pareto distribution
			is constructed relative to the exponential distribution, thus making
			it have no effect on the initial wealth shares. 
		-a: changes the initial distribution pareto exponent "alpha" [default: 
			1.55]
		
		***ANALYSIS ARGS***
		
		-r: tells the program to run a simulation; if the -r flag is not given,
			then the program defaults to an analysis of the data file specified 
			by the -f flag ("market" by default).
		
		-f: tells the program which file to analyze/save run data in; files are
			saved in "data/filename.data"
			[default: market]
			
		-p: tells the program what percentile interval to do power-law fit on
			[default: 0.9]
			
		-nc: tells the program to not compress the arranged individuals data;
			this means that cumulative distribution function will not be exactly
			a cumulative distribution function some of the time; this should
			only be used for debugging purposes.
			
		-rd: tells the program to do analysis/run the simulation in a periodic
			format; this doesn't actually run the program (-r is for that), but
			should be used when producing/analyzing data that involve multiple
			snapshots of the market
			
		-c: tells the program to compare simulation statistics with real U.S.
			statistics; it should only be used in conjunction with the -rd flag
	Example:
		
		python run_model.py -n 1000 -gr 0.03 -dw 0.01 -t 25 -b 1.3 -ec .75 -pc .9 -mp .75 -g 2.78 -a 1.5 -r -f my_run -p .9 -c -rd
		
		This command will run a new simulation over 25 years, 
		that provides periodic snapshots (every year) of the wealth shares and 
		pareto and gini indexes, with 1000 individuals in an economy of a growth
		rate of 3% annually, with a delta omega of 0.01. The distribution will
		be in the following form (Percentile Interval is in [] on left):
		
		"p" in the following 3 lines is the percentile:
		
		[0, .75] will be characterized by: log(p)/(-g)
		[.75, .9] will be characterized by: [(.9-p)/(.9-.75) * log(p)/(-g)] + [1 - (.9-p)/(.9-.75) * w_m * p^(-1/a)]
		[.9, 1] wll be characterized by: w_m * p^(-1/a)
		
		The run data will be saved into data/my_run.data and then analysis will
		be run on that file, comparing statistics with U.S. market share data,
		with the assumption that the pareto fit should be on the top 10% of the
		wealthy in the market.
		
		
	series_model.py:
		This program was built with the intention of validating the variability
		in the model's fit to U.S. market data. Run:
		
		python series_model.py
		
		to see statistics on the nonuniform initial start (i.e. variability)
		
		If you would like to run new variability checks, you must go into the
		source and modify some code, as I did not program in any arguments (We 
		only needed to run it once). If you would like to run a new validation
		of the model, change the "run" variable to "True". Comparable implies
		that the model should not calculate means squared error, but it is 
		to leave that variable untouched.
	
### c_parallel:
	This folder contains the program that executes the single-packet model
	starting from a uniform distribution.
	
	To view the results that are presented in the paper, run the following code
	from inside the c_parallel directory:
	
	python run_model.py -f great_run_10000_225_1.36 -p .9
	
	Description of flags:
		***RUNNING ARGS***
		-n: change the number of individuals [default: 1000]
		-g: change the growth rate of the entire economy [default: 0.03] 
		-i: initial wealth of each individual [default: 1]
		-w: change the delta_omega (precision of the model) [default: 0.01]
		-t:	change the total time that the model runs [default: 100]
		-e:	change the exponent "beta" in wealth power function [default: 1]
		
		***ANALYSIS ARGS***
		
		-r: tells the program to run a simulation; if the -r flag is not given,
			then the program defaults to an analysis of the data file specified 
			by the -f flag ("market" by default).
		
		-f: tells the program which file to analyze/save run data in; files are
			saved in "data/filename.data"
			[default: market]
			
		-p: tells the program what percentile interval to do power-law fit on
			[default: 0.0]
	
	Example:
	
		python run_model.py -n 100 -i 1 -w 0.01 -g 0.03 -t 100 -e 1 -r -f run -p .9
		
		This will run the simulation on 100 individuals each starting with 1
		for their wealth and wealth increments at 0.01. The market will grow 3% 
		annually for 100 years with an wealth power exponent of 1. The run will 
		be saved as data/run.data. The power-law fit will be on the percentile
		interval [.9, 1].
		
	find_model.py:
		This program was built to do a sort of "Monte Carlo" search for good
		parameters to run the model with. Look into the source code for more
		details.
