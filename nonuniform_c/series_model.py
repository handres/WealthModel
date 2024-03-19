import sys
import os
from matplotlib import pyplot as pp
from math import log, exp
from scipy import stats
from scipy.optimize import curve_fit
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
		pp.plot(x_p, stat, linestyle = "--", label = l)
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

def run_model(p_id, returner):
	# Running Parameters
	num_individuals = 100000
	market_growth_rate = 0.02588
	delta_omega = 0.01
	time_period = 25
	beta = 1.25 # Power function exponent
	
	
	# Initial Distribution Parameters
	exp_cutoff = .75 # When to start transforming into Pareto distribution
	pareto_cutoff = .9 # When to stop transforming into Pareto distribution
	min_percent = .75 # To calculate the "minimum" in the Pareto distribution
	gamma = 2.78 # Gamma in Exponential Distribution
	alpha = 1.55 # Alpha in Pareto distribution
	
	# Analysis Parameters
	run = False
	output_file = "market" + str(p_id)
	percentile = 0.9
	compress = True
	record_all_data = 1
	compare_data = False
	
	
	
	alphas = []
	pstats = []
	ginis = []
	
	
	
	os.system("./nonuniform {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}".format(num_individuals, market_growth_rate, delta_omega, time_period, beta, exp_cutoff, pareto_cutoff, min_percent, gamma, alpha, record_all_data, output_file))
	
	# Read data file
	f = open("data/" + output_file + ".data")
	for line in f.readlines():
		individuals = []
		for entry in line.split(","):
			individuals.append(float(entry))
		pstats.append(calc_stats(individuals))
		alphas.append(-get_log_worth(individuals, percentile=percentile)[2][0])
		ginis.append(calc_gini(individuals))
	returner[p_id] = (pstats, alphas, ginis)


def read_comparison_data():
	f = open("percentiles.src")			
	
	rsrc = []
	for line in f:
		if "bottom" in line:
			continue
		rsrc.append([])
		for entry in line.split(","):
			rsrc[-1].append(float(entry))

	f.close()
	
	f = open("data/1.2_runs_compare_24/market.src")
	
	psrc = []
	
	palpha = []
	pgini = []
	
	for i, line in enumerate(f.readlines()):
		rem = i % 4
		if rem == 0:
			psrc.append([])
			for entry in line.split(")")[:-1]:
				psrc[-1].append([])
				for stat in entry.split(","):
					if "[" in stat:
						psrc[-1][-1].append(float(stat[2:]))
					elif "(" in stat:
						psrc[-1][-1].append(float(stat[2:]))
					elif stat != "":
						psrc[-1][-1].append(float(stat))
		elif rem == 1:
			palpha.append([])
			for entry in line.split(","):
				if "[" in entry:
					palpha[-1].append(float(entry[1:]))
				elif "]" in entry:
					palpha[-1].append(float(entry[:-2]))
				else:
					palpha[-1].append(float(entry))
		elif rem == 2:
			pgini.append([])
			for entry in line.split(","):
				if "[" in entry:
					pgini[-1].append(float(entry[1:]))
				elif "]" in entry:
					pgini[-1].append(float(entry[:-2]))
				else:
					pgini[-1].append(float(entry))	
		else:
			pass
			
	f.close()
	
	return (rsrc, psrc, palpha, pgini)

def calc_comparison_data():
	n = 25 # Number of times to run the model
	os.system("g++ nonuniform.cpp -std=c++11 -O3 -fopenmp -o nonuniform")

	proc_num = 4
	proc = [None]*proc_num
	results = []
	manager = multiprocessing.Manager()
	returner = manager.dict()

	for i in range(n):
		for j in range(proc_num):
			proc[j] = multiprocessing.Process(target=run_model, args=(j, returner))
			proc[j].start()
		for p in proc:
			p.join()
		results.append(returner.values())

	#print results
	# Write output to file

	#print results[0][0]
	f = open("data/1.2_runs_compare_24/market.src", 'w')
	for i in range(len(results)):
		for j in range(len(results[i])):
			for k in range(len(results[i][j])):
				f.write(str(results[i][j][k]) + "\n")
			f.write("--------------------------------\n")
	f.close()


def compare_data(comp, data):

	rsrc = zip(*data[0])
	psrc = data[1]
	tpsrc = zip(*data[1])
	alpha = zip(*data[2])
	gini = zip(*data[3])

	
	ms_error = []
	avg_mse = []
	src_var = []
	alpha_var = []
	gini_var = []
	
	
	if comp: # Mean Squared Error
		#print psrc
		for run in psrc:
			ms_error.append([])
			t_run = zip(*run)
			
			for pstat, rstat in zip(t_run[:len(rsrc)], rsrc):
				temp = 0
				for pentry, rentry in zip(pstat, rstat):
					temp += (pentry - rentry)**2
				
				ms_error[-1].append(temp/len(pstat))
	
		for mse in ms_error:
			avg_mse.append(sum(mse)/len(mse))
		
	for i in range(len(alpha)):
		alpha_var.append(calc_variance(alpha[i]))
		gini_var.append(calc_variance(gini[i]))
		ttpsrc = zip(*tpsrc[i])
		src_var.append([])
		for stat in ttpsrc:
			src_var[-1].append(calc_variance(stat))
			
	
	#raw_input()
		
			
	print "\nVariance of Model Wealth Shares:\n"
	print src_var
	print "\nVariance of Model Alpha Values:\n"
	print alpha_var
	print "\nVariance of Model Gini Values:\n"
	print gini_var
	print "\nMeans Squared error, [per-attribute][per-run]"
	print ms_error
	
	

def calc_variance(arr):
	#print arr
	m = float(sum(arr))/len(arr)
	temp = 0
	for x in arr:
		temp += (x-m)**2
	#print "-------"
	#print temp/len(arr)
	#print "-------"
	return temp/len(arr)
	

def run_single_check(output_file):
	pstats = []
	alphas = []
	ginis = []
	f = open("data/" + output_file + ".data")
	for line in f.readlines():
		individuals = []
		for entry in line.split(","):
			individuals.append(float(entry))
		pstats.append(calc_stats(individuals))
		alphas.append(-get_log_worth(individuals, percentile=.9)[2][0])
		ginis.append(calc_gini(individuals))
	return (pstats, alphas, ginis)
	

if __name__ == "__main__":
	
	
	run = False
	comparable = True

		
	
	if run:
		calc_comparison_data()
	
	#print read_dat_file("100000_100years_diverge")[1]
	#compare_data(True, [read_dat_file("100000_100years_diverge")])
	compare_data(comparable, read_comparison_data())
			
			
			
			
