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
	dat.append(list(reversed(sorted(get_unsorted_worth(individual_worths))))[:int(percentile*len(individual_worths))])
	dat.append([])
	for i in range(len(dat[0])):
		dat[1].append(float(i)/len(dat[0]))
	
	average = sum(dat[0])/len(dat[0])
	
	def func(x, a, b):
		return a*np.exp(-np.asarray(x)*b)
	
	popt, pcov = curve_fit(func, np.asarray(dat[0]), np.asarray(dat[1]))
	
	print "Exponential normalization: {0}".format(popt)
	fit.append(dat[0])
	fit.append(func(dat[0], *popt))
	
	
	return dat, fit

def plot_individuals(individual_worths):
	pp.plot(get_unsorted_worth(individual_worths), linewidth=0.001, marker = 'o')

def plot_sorted_individuals(individual_worths):
	pp.plot(sorted(get_unsorted_worth(individual_worths))[:int(.9*len(individual_worths))], linewidth=0.001, marker = 'o')

def plot_arranged_individuals(individual_worths, points=True, percentile=0.0, compress=True, thickness=0.001):
	a_temp_arr = get_arranged_worth(individual_worths, percentile=percentile, compress=compress)
	if points:
		pp.plot(a_temp_arr[0], a_temp_arr[1], linewidth = thickness, marker = 'o')
	else:
		pp.plot(a_temp_arr[0], a_temp_arr[1], linewidth = thickness)
def plot_log_individuals(individual_worths, compress=True, percentile=0.0):

	log_temp_arr, linear_fit, stats = get_log_worth(individual_worths, percentile=percentile, compress=compress)	
		
	pp.plot(log_temp_arr[0], log_temp_arr[1], linewidth=0.001, marker='o', label="Data")
	pp.plot(log_temp_arr[0], linear_fit, linestyle='--', label="Fit")
	pp.xlabel("Log(Wealth)")
	pp.ylabel("Log(Percentile)")
	pp.legend()
	print "Linear Fit:\tSlope: {0}\tIntercept: {1}\tR^2 {2}\tP: {3}\tStandard_Error: {4}".format(stats[0], stats[1], stats[2]**2, stats[3], stats[4])
	print "Calculated Eta Value:\t {0}".format(-1.0/stats[0])

def plot_exp_individuals(indvidual_worths, percentile=0.0):

	act, fit = get_exp_worth(individual_worths, percentile)
	pp.plot(act[0], act[1], marker='o', linewidth=0.001, label="Data")
	pp.plot(fit[0], fit[1], linestyle='--', linewidth=1.5, label="Fit")
	pp.ylabel("Percentile")
	pp.xlabel("Wealth")
	pp.legend()

def calc_gini(individual_worths):
	temp = 0.0
	for i in individual_worths:
		for j in individual_worths:
			temp+= abs(i - j)
	return temp/(2*len(individual_worths)*sum(individual_worths))
	

if __name__ == "__main__":
	run = False
	output_file = "market"
	num_individuals = 1000
	init_worth_per_individual = 100
	market_growth_rate = 0.03
	delta_omega = 0.1
	time_period = 50
	market_type = 0 # 0 = linear, 1 = exponential, 2 = polynomial
	exponent_power = 1
	percentile = 0.0
	compress = True

	for i, arg in enumerate(sys.argv):
		if arg == "-r": # Run the c++ code (the model)
			run = True
		elif arg == "-f": # Set the filename to data/*name*.data
			output_file = sys.argv[i+1]
		elif arg == "-n": # Change the number of individuals in the model
			num_individuals = int(sys.argv[i+1]) 
		elif arg == "-g": # Change the growth rate of the entire economy
			growth_rate = float(sys.argv[i+1]) 
		elif arg == "-i": # Change the initial worth per individual
			init_worth_per_individual = float(sys.argv[i+1])
		elif arg == "-w": # Change the delta_omega (precision of steady state model)
			delta_omega = float(sys.argv[i+1]) # 
		elif arg == "-t": # Change the time period over which to run the model
			time_period = float(sys.argv[i+1])
		elif arg == "-p": # Change which percentile to compute a linear regression on
			percentile = float(sys.argv[i+1])
		elif arg == "-nc": # Don't compress the arranged individuals (Not really useful to flag this)
			compress = False
		elif arg == "-e":
			exponent_power = float(sys.argv[i+1])


		
	if run:
		os.system("g++ fast_model.cpp -std=c++11 -O3 -fopenmp -o market_model")
		os.system("./market_model {0} {1} {2} {3} {4} {5} {6} {7}".format(num_individuals, init_worth_per_individual, market_growth_rate, delta_omega, time_period, market_type, exponent_power, output_file))
		
	
	
	f = open("data/" + output_file + ".data")
	
	lines = f.readlines()
	
	individual_worths = []
	
	for line in lines:
		individual_worths.append(float(line))
	
	pp.figure("Unsorted Individuals")
	pp.title("Unsorted Individuals")
	plot_individuals(individual_worths)

	# pp.figure("Sorted Individuals")
	# pp.title("Sorted Individuals")
	# plot_sorted_individuals(individual_worths)

	pp.figure("Arranged Individuals")
	pp.title("Arranged Individuals")
	plot_arranged_individuals(individual_worths, compress=compress, percentile=percentile)

	pp.figure("Log-Log and Linear Fit Top {0}% of Market".format((1 - percentile)*100))
	pp.title("Log-Log and Linear Fit Top {0}% of Market".format((1 - percentile)*100))
	plot_log_individuals(individual_worths, compress=compress, percentile=percentile)
	
	pp.figure("Log-log of Cumulative Distribution of Entire Market")
	pp.title("Log-log of Cumulative Distribution of Entire Market")
	log_arr, temp1, temp2 = get_log_worth(individual_worths)
	pp.plot(log_arr[0], log_arr[1], marker='o', linewidth=0.001)
	pp.xlabel("Log(Wealth)")
	pp.ylabel("Log(Percentile)")
	
	pp.figure("Exponential and Exponential Fit of Market")
	pp.title("Exponential and Exponential Fit of Bottom 90% Market")
	plot_exp_individuals(individual_worths, percentile=.9)
	# print "Gini Coefficient: " + str(calc_gini(individual_worths))

	pp.show()
