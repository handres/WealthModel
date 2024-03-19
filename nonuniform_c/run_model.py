import sys
import os
from matplotlib import pyplot as pp
from math import log, exp
from scipy import stats
from scipy.optimize import curve_fit
import numpy as np

def get_unsorted_worth(individual_worths):
	temp_arr = [0]*len(individual_worths)
	for i in range(len(individual_worths)):
		temp_arr[i] = individual_worths[i]
	return temp_arr

def get_arranged_worth(individual_worths, compress=True, percentile=0.0):
	s_temp_arr = sorted(get_unsorted_worth(individual_worths))
	a_temp_arr = [[], []]
	
	if compress:
		last_worth = 0
		for i in range(int(len(s_temp_arr)*percentile), len(s_temp_arr)):
			if last_worth < a_temp_arr[0]:
				a_temp_arr[0].append(s_temp_arr[i])
				a_temp_arr[1].append((len(s_temp_arr) - i) / float(len(s_temp_arr)))
			else:
				a_temp_arr[1][-1] += 1
	else:
		for i in range(len(s_temp_arr[int(percentile*len(s_temp_arr)):])):
			a_temp_arr[0].append(s_temp_arr[i])
			a_temp_arr[1].append((len(s_temp_arr) - i) / float(len(s_temp_arr)))
	return a_temp_arr

def get_log_worth(individual_worths, compress=True, percentile=0.0):
	a_temp_arr = get_arranged_worth(individual_worths, compress=compress, percentile=percentile)
	log_temp_arr = [[], []]*len(a_temp_arr)
	
	for i in range(len(a_temp_arr[0])):
		log_temp_arr[0].append(log(a_temp_arr[0][i], 10))
		log_temp_arr[1].append(log(a_temp_arr[1][i], 10))
	
	

	slope, intercept, r_value, p_value, std_err = stats.linregress(log_temp_arr[0], log_temp_arr[1])

	linear_fit = []	

	for entry in log_temp_arr[0]:
		linear_fit.append(entry*slope + intercept)
	
	return log_temp_arr, linear_fit, (slope, intercept, r_value, p_value, std_err)

def get_exp_worth(individual_worths, percentile=0.0):
	dat = []
	fit = []
	total = sum(individual_worths)
	dat.append(list(reversed(sorted(get_unsorted_worth(individual_worths))))[:int(percentile*len(individual_worths))])
	dat.append([])
	for i in range(len(dat[0])):
		# dat[0][i] = dat[0][i]/total
		dat[1].append(float(i)/len(dat[0]))
	
	#average = sum(dat[0])/len(dat[0])
	
	def func(x, a, b):
		return a*np.exp(-np.asarray(x)/b)
	
	popt, pcov = curve_fit(func, np.asarray(dat[0]), np.asarray(dat[1]), p0=(1, 5))
	
	print "Exponential normalization: {0}".format(popt)
	fit.append(dat[0])
	fit.append(func(dat[0], *popt))
	
	
	return dat, fit

def plot_individuals(individual_worths):
	pp.xlabel("Individuals")
	pp.ylabel("Wealth")
	pp.plot(get_unsorted_worth(individual_worths), linewidth=0.001, marker = 'o', markersize=2)

def plot_sorted_individuals(individual_worths):
	pp.plot(sorted(get_unsorted_worth(individual_worths))[:int(.9*len(individual_worths))], linewidth=0.001, marker = 'o', markersize=2)

def plot_arranged_individuals(individual_worths, points=True, percentile=0.0, compress=True, thickness=0.001):
	a_temp_arr = get_arranged_worth(individual_worths, percentile=percentile, compress=compress)
	if points:
		pp.plot(a_temp_arr[0], a_temp_arr[1], linewidth = thickness, marker = 'o', markersize=2)
	else:
		pp.plot(a_temp_arr[0], a_temp_arr[1], linewidth = thickness)
def plot_log_individuals(individual_worths, compress=True, percentile=0.0):

	log_temp_arr, linear_fit, stats = get_log_worth(individual_worths, percentile=percentile, compress=compress)	
		
	pp.plot(log_temp_arr[0], log_temp_arr[1], linewidth=0.001, marker='o', markersize=2, label="Data")
	pp.plot(log_temp_arr[0], linear_fit, linestyle='--', label="Fit")
	pp.xlabel("Log(Wealth)")
	pp.ylabel("Log(Percentile)")
	pp.legend()
	print "Linear Fit:\tSlope: {0}\tIntercept: {1}\tR^2 {2}\tP: {3}\tStandard_Error: {4}".format(stats[0], stats[1], stats[2]**2, stats[3], stats[4])
	print "Calculated Eta Value:\t {0}".format(-1.0/stats[0])

def plot_exp_individuals(indvidual_worths, percentile=0.0):

	act, fit = get_exp_worth(individual_worths, percentile)

	pp.plot(act[0], act[1], marker='o', linewidth=0.001, label="Data", markersize=2)
	pp.plot(fit[0], fit[1], linestyle='--', linewidth=1.5, label="Fit")
	pp.ylabel("Percentile")
	pp.xlabel("Wealth")
	pp.legend()


def plot_stats_compare(pstats, rstats, alphas, ginis):
	
	trans_pstats = []
	trans_rstats = []
	x_r = range(1988, 2013)
	x_p = range(1988, 1988 + len(pstats))
	for i in range(len(pstats[0])):
		trans_pstats.append([])
		for j in range(len(pstats)):
			trans_pstats[i].append(pstats[j][i])
	for i in range(len(rstats[0])):
		trans_rstats.append([])
		for j in range(len(rstats)):
			trans_rstats[i].append(rstats[j][i])
	for i, stat in enumerate(trans_pstats):
		pp.plot(x_p, stat, linestyle="--")
	for i, stat in enumerate(trans_rstats):
		if(i==0):
			l = "Data Bottom 90%"
		elif(i==1):
			l = "Data Top 10%"
		elif(i==2):
			l = "Data Top 5%"
		elif(i==3):
			l = "Data Top 1%"
		elif(i==4):
			l = "Data Top 0.5%"
		elif(i==5):
			l = "Data Top 0.1%"
		else:
			l = "Data Top 0.01%"
		pp.plot(x_r, stat, label=l)
		
	pp.plot(x_p, alphas, marker = 'o', linewidth = .001, label="Pareto index")
	pp.plot(x_p, ginis, marker = 's', linewidth = .001, label="Gini coefficient")
	pp.xlabel("Year")
	pp.ylabel("Wealth Share            Index Value")
	
def plot_stats_proj(stats, alphas, ginis):
	trans_pstats = []
	x = range(1988, 1988 + len(pstats))
	for i in range(len(pstats[0])):
		trans_pstats.append([])
		for j in range(len(pstats)):
			trans_pstats[i].append(pstats[j][i])
	

	for i, stat in enumerate(trans_pstats):
		if(i==0):
			l = "Model Bottom 90%"
		elif(i==1):
			l = "Model Top 10%"
		elif(i==2):
			l = "Model Top 5%"
		elif(i==3):
			l = "Model Top 1%"
		elif(i==4):
			l = "Model Top 0.5%"
		elif(i==5):
			l = "Model Top 0.1%"
		else:
			l = "Model Top 0.01%"
		pp.plot(x, stat, linestyle = "--", label = l)
		
	pp.plot(x, alphas, marker = 'o', linewidth = .001, label="Pareto index")
	pp.plot(x, ginis, marker = 's', linewidth = .001, label="Gini coefficient")
	pp.xlabel("Year")
	pp.ylabel("Wealth Share            Index Value")
def calc_gini(individual_worths):
	temp = 0.0
	
	s_individual_worths = sorted(individual_worths)
	n = len(individual_worths)
	for i in range(n):
		temp += i*s_individual_worths[i]
	
	temp *= 2
	
	return temp / (n*sum(individual_worths)) - (n+1)/n
	
	
def calc_stats(individual_worths):
	s_individual_worths = sorted(individual_worths)
	total_wealth = sum(individual_worths)
	
	bottom90 = sum(s_individual_worths[:int(.9*len(individual_worths))])/total_wealth
	top10 = sum(s_individual_worths[int(.9*len(individual_worths)):])/total_wealth
	top5 = sum(s_individual_worths[int(.95*len(individual_worths)):])/total_wealth
	top1 = sum(s_individual_worths[int(.99*len(individual_worths)):])/total_wealth
	topp5 = sum(s_individual_worths[int(.995*len(individual_worths)):])/total_wealth
	topp1 = sum(s_individual_worths[int(.999*len(individual_worths)):])/total_wealth
	toppp1 = sum(s_individual_worths[int(.9999*len(individual_worths)):])/total_wealth
	
	return (bottom90, top10, top5, top1, topp5, topp1, toppp1)

if __name__ == "__main__":
	
	# Default Paramters
	
	# Running Parameters
	num_individuals = 100000
	market_growth_rate = 0.02588
	
	market_growth_rate = 0.
	delta_omega = 0.01
	time_period = 25
	beta = 1.25 # Power function exponent
	
	
	# Initial Distribution Parameters
	exp_cutoff = .75 # When to start transforming into Pareto distribution
	pareto_cutoff = .999 # When to stop transforming into Pareto distribution
	min_percent = .75 # To calculate the "minimum" in the Pareto distribution
	gamma = 2.78 # Gamma in Exponential Distribution
	alpha = 1.55 # Alpha in Pareto distribution
	
	# Analysis Parameters
	run = False
	output_file = "market"
	percentile = 0.9
	compress = True
	record_all_data = 0
	compare_data = False


	for i, arg in enumerate(sys.argv):
		# Check for running parameter args
	
		if arg == "-n": # Change the number of individuals in the model
			num_individuals = int(sys.argv[i+1]) 
		elif arg == "-gr": # Change the growth rate of the entire economy
			market_growth_rate = float(sys.argv[i+1]) 
		elif arg == "-dw": # Change the delta_omega (precision of steady state model)
			delta_omega = float(sys.argv[i+1]) # 
		elif arg == "-t": # Change the time period over which to run the model
			time_period = float(sys.argv[i+1])
		elif arg == "-b":
			beta = float(sys.argv[i+1])
		
		# Check for initial distribution args	
		
		elif arg == "-ec":
			exp_cutoff = float(sys.argv[i+1])
		elif arg == "-pc":
			pareto_cutoff = float(sys.argv[i+1])
		elif arg == "-mp":
			min_percent = float(sys.argv[i+1])
		elif arg == "-g":
			gamma = float(sys.argv[i+1])
		elif arg == "-a":
			alpha = float(sys.argv[i+1])
		
		# Check for analysis args
		elif arg == "-r": # Run the c++ code (the model)
			run = True
		elif arg == "-f": # Set the filename to data/*name*.data
			output_file = sys.argv[i+1]
		elif arg == "-p": # Change which percentile to compute a linear regression on
			percentile = float(sys.argv[i+1])
		elif arg == "-nc": # Don't compress the arranged individuals (Not really useful to flag this)
			compress = False
		
		elif arg == "-rd":
			record_all_data = 1
		
		elif arg == "-c":
			compare_data = True
		

				
	if run:
		os.system("g++ nonuniform.cpp -std=c++11 -O3 -fopenmp -o nonuniform")
		os.system("./nonuniform {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}".format(num_individuals, market_growth_rate, delta_omega, time_period, beta, exp_cutoff, pareto_cutoff, min_percent, gamma, alpha, record_all_data, output_file))
		
	
	if record_all_data == 0:
	
		f = open("data/" + output_file + ".data")
	
		lines = f.readlines()
	
		individual_worths = []
	
		for line in lines:
			individual_worths.append(float(line))
	
		s_individual_worths = sorted(individual_worths)
		total_wealth = sum(individual_worths)
	
		bottom90 = sum(s_individual_worths[:int(.9*len(individual_worths))])
		top10 = sum(s_individual_worths[int(.9*len(individual_worths)):])
		top5 = sum(s_individual_worths[int(.95*len(individual_worths)):])
		top1 = sum(s_individual_worths[int(.99*len(individual_worths)):])
		topp5 = sum(s_individual_worths[int(.995*len(individual_worths)):])
		topp1 = sum(s_individual_worths[int(.999*len(individual_worths)):])
		toppp1 = sum(s_individual_worths[int(.9999*len(individual_worths)):])
	
		print "Bottom 90%: {0}\nTop 10% {1}\nTop 5%: {2}\nTop 1%: {3}\nTop 0.5%: {4}\nTop 0.1%: {5}\nTop 0.01%: {6}".format(bottom90/total_wealth, top10/total_wealth, top5/total_wealth, top1/total_wealth, topp5/total_wealth, topp1/total_wealth, toppp1 / total_wealth)
	
	
	
		pp.figure("Unsorted Individuals")
		pp.title("Unsorted Individuals")
		plot_individuals(individual_worths)


	#	pp.figure("Sorted Individuals")
	#	pp.title("Sorted Individuals")
	#	plot_sorted_individuals(individual_worths)

	#	pp.figure("Arranged Individuals")
	#	pp.title("Arranged Individuals")
	#	plot_arranged_individuals(individual_worths, compress=compress, percentile=percentile)

		pp.figure("Log-Log and Linear Fit Top {0}% of Market".format((1 - percentile)*100))
		pp.title("Log-Log and Linear Fit Top {0}% of Market".format((1 - percentile)*100))
		plot_log_individuals(individual_worths, compress=compress, percentile=percentile)
	
		pp.figure("Log-log of Cumulative Distribution of Entire Market")
		pp.title("Log-log of Cumulative Distribution of Entire Market")
		log_arr, temp1, temp2 = get_log_worth(individual_worths)
		pp.plot(log_arr[0], log_arr[1], marker='o', linewidth=0.001, markersize=2)
		pp.xlabel("Log(Wealth)")
		pp.ylabel("Log(Percentile)")
		
		pp.figure("Exponential and Exponential Fit of Market")
		pp.title("Exponential and Exponential Fit of Bottom 90% Market")
		plot_exp_individuals(individual_worths, percentile=.9)
		print "Gini Coefficient: " + str(calc_gini(individual_worths))
	else:
		alphas = []
		pstats = []
		ginis = []
		
		f = open("data/" + output_file + ".data")
		g = open("percentiles.src")
		
		
		# Read data file
		for line in f.readlines():
			individuals = []
			for entry in line.split(","):
				individuals.append(float(entry))
			pstats.append(calc_stats(individuals))
			alphas.append(-get_log_worth(individuals, percentile=percentile)[2][0])
			ginis.append(calc_gini(individuals))
		# Read src file
		
		rstats = []
		for line in g:
			if "bottom" in line:
				continue
			rstats.append([])
			for entry in line.split(","):
				rstats[-1].append(float(entry))
		pp.figure("Stats")
		#pp.title("Stats")
		if compare_data:
			plot_stats_compare(pstats, rstats, alphas, ginis)
		else:
			plot_stats_proj(pstats, alphas, ginis)
	pp.legend(fontsize=10, loc="upper center", bbox_to_anchor=(0.5,1.05), ncol=3, fancybox=True, shadow=True)
	pp.show()
