import sys
import os
from matplotlib import pyplot as pp
from math import log, exp
from scipy import stats
from scipy.optimize import curve_fit, fsolve
import numpy as np

import multiprocessing

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
				a_temp_arr[1].append((len(s_temp_arr) - i))
			else:
				a_temp_arr[1][-1] += 1
	else:
		for i in range(len(s_temp_arr[int(percentile*len(s_temp_arr)):])):
			a_temp_arr[0].append(s_temp_arr[i])
			a_temp_arr[1].append(len(s_temp_arr) - i)
	return a_temp_arr

def get_log_worth(individual_worths, compress=True, percentile=0.0):
	a_temp_arr = get_arranged_worth(individual_worths, compress=compress, percentile=percentile)
	log_temp_arr = [[], []]*len(a_temp_arr)
	
	for i in range(len(a_temp_arr[0])):
		log_temp_arr[0].append(log(a_temp_arr[0][i]))
		log_temp_arr[1].append(log(a_temp_arr[1][i]))
	
	

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
	print average
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
		
	pp.plot(log_temp_arr[0], log_temp_arr[1], linewidth=0.001, marker='o')
	pp.plot(log_temp_arr[0], linear_fit)
	
	print "Linear Fit:\tSlope: {0}\tIntercept: {1}\tR^2 {2}\tP: {3}\tStandard_Error: {4}".format(stats[0], stats[1], stats[2]**2, stats[3], stats[4])
	print "Calculated Eta Value:\t {0}".format(-1.0/stats[0])

def plot_exp_individuals(indvidual_worths, percentile=0.0):

	act, fit = get_exp_worth(individual_worths, percentile)
	
	
	pp.plot(fit[0], fit[1])
	pp.plot(act[0], act[1], marker='o', linewidth=0.001)
	

def calc_gini(individual_worths):
	temp = 0.0
	for i in individual_worths:
		for j in individual_worths:
			temp+= abs(i - j)
	return temp/(2*len(individual_worths)*sum(individual_worths))
	
	
	
def run_model(args, p_id, returner):
	output_file = "market" + str(p_id)
	num_individuals = 1000
	init_worth_per_individual = 1
	market_growth_rate = 0.03
	time_period = 50
	market_type = 2 # 0 = linear, 1 = exponential, 2 = polynomial
	percentile = .9
	compress = True
	os.system("./market_model {0} {1} {2} {3} {4} {5} {6} {7}".format(num_individuals, init_worth_per_individual, market_growth_rate, args[0], time_period, market_type, args[1], output_file))
	f = open("data/" + output_file + ".data")
	lines = f.readlines()
	individual_worths = []
	for line in lines:
		individual_worths.append(float(line))
	stats = get_log_worth(individual_worths, compress=compress, percentile=percentile)[2]
	returner[p_id] = (stats[0], stats[2]**2)
	
	
	
	
def avg_metrics(seq):
	avg = []
	for i in range(len(seq[0])):
		total = 0
		for j in range(len(seq)):
			total += seq[j][i]
		avg.append(total/len(seq))
	return avg
		
if __name__ == "__main__":
	init_delta_omega = 0.01
	init_exp = 2
	end_exp = 2.2
	step_exp = .1
	n = int((end_exp - init_exp)/step_exp) + 1
	args = [init_delta_omega, []]
	proc_num = 4
	proc = [None]*proc_num
	os.system("g++ fast_model.cpp -o market_model -O3 -fopenmp -std=c++11")
	results = []
	manager = multiprocessing.Manager()
	returner = manager.dict()
	for i in range(n):
		args[1].append(init_exp + i*step_exp)
		
		for j in range(proc_num):
			proc[j] = multiprocessing.Process(target=run_model, args=((args[0], args[1][i]),j, returner))
			proc[j].start()
		for p in proc:
			p.join()
		results.append(returner.values())
			
	
	
	all_results = []
	for i, res in enumerate(results):
		avg = avg_metrics(res)
		all_results.append(avg)
		if (-3 < avg[0] < -1 and avg[1] > 0.9):
			print "Result: {0}\tdw={1}\texp={2}".format(avg, init_delta_omega, init_exp + i*step_exp)
	
	for res in all_results:
		print res
	
	
