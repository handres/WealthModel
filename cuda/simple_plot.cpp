#include <iostream>
#include <fstream>
#include <gsl/gsl_fit.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <string>
#include <cmath>

#define PERCENTILE 0.9
#define VERBOSE true
#define DEBUG false

#define rpt_err(_cond,_msg) if(_cond){cerr << "ERROR: " << _msg << endl;}
#define dbg(_msg) if(DEBUG){cout << _msg;}
typedef float WT;

using namespace std;


WT read_WT(FILE *f)
{
	WT v;
	fread((void*)(&v), sizeof(v), 1, f);
	return v;
}

vector<WT> get_raw_data(string filename)
{
	FILE* f = fopen(filename.c_str(), "rb");
	vector<WT> data(0);
	dbg("Reading in data from file: " << filename << "\n");
	
	while(!feof(f))
	{
		 data.push_back(read_WT(f));
	}
	return data;
}

void populate_data(vector<WT> &data, vector<long> &n_data, vector<WT> &w_data, vector<double> &log_n_data, vector<double> &log_w_data, vector<double> &linear_fit_x, vector<double> &linear_fit_y, vector<double> &exp_fit_x, vector<double> &exp_fit_y)
{
	dbg("Creating arranged worth and log array...");

	long idx = 0;
	long len = data.size();
	w_data.push_back(data[0]);
	n_data.push_back(len);
	for(long i=1;i<len;i++)
	{
		if(data[i] > w_data[idx])
		{
			n_data.push_back(len-i);
			w_data.push_back(data[i]);
			log_n_data.push_back(log(((double)(len-i))/len));
			log_w_data.push_back(log(data[i]));

			idx++;
		}
	}

	dbg("Done.\n");

	dbg("Creating array to linearly fit...");

	long ranks = log_n_data.size();
	long perc_cutoff = ranks*PERCENTILE;
	
	double* x = new double[(long)(len-perc_cutoff)];
	double* y = new double[(long)(len-perc_cutoff)];
	
	
	for(long i=perc_cutoff;i<ranks;i++)
	{
		x[i - perc_cutoff] = log_w_data[i];
		y[i - perc_cutoff] = log_n_data[i];
	}
	dbg("Done.\n");
	
	dbg("Fitting pareto distribution...");
	double c0, c1, cov00, cov01, cov11, sumsq;
	
	gsl_fit_linear(x, 1, y, 1, ranks-perc_cutoff, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

	dbg("Done.\n");
#if VERBOSE
	cout << "Linear fit: " << "\n\tSlope: " << c1 << "\n\tIntercept: " << c0 << "\n\tCov: [[" << cov00 << ", " << cov01 << "], [" << cov01 << ", " << cov11 << "]]" << "\n\tR^2: " << pow(cov01, 2)/(cov00*cov11) << endl;
#endif
	
	dbg("Creating linear fit data...");
	for(long i=0;i<ranks-perc_cutoff;i++)
	{
		linear_fit_x.push_back(x[i]);
		linear_fit_y.push_back(x[i]*c1 + c0);
	}
	dbg("Done.\n");
	
	delete x, y;
	
	dbg("Creating array to exponentially fit...");
	
	
	double ec0, ec1, ecov00, ecov01, ecov11, esumsq;
	x = new double[perc_cutoff];
	y = new double[perc_cutoff];

	for(long i=0;i<perc_cutoff;i++)
	{
		x[i] = log_n_data[i];
		y[i] = log_w_data[i];		
	}
	dbg("Done.\n");
	
	dbg("Fitting exponential...");
	
	
	gsl_fit_linear(x, 1, y, 1, ranks-perc_cutoff, &ec0, &ec1, &ecov00, &ecov01, &ecov11, &esumsq);
	
	dbg("Done.\n")
	
#if VERBOSE
	cout << "Exponential fit: " << "\n\tWealth Temperature: " << -1/ec1 << "\n\tScaling Factor: " << ec0 << "\n\tCov: [[" << ecov00 << ", " << ecov01 << "], [" << ecov01 << ", " << ecov11 << "]]" << "\n\tR^2: " << pow(ecov01, 2)/(ecov00*ecov11) << endl;
#endif
	
	dbg("Creating exponential fit data...");
	for(long i=0;i<perc_cutoff;i++)
	{
		exp_fit_x.push_back(w_data[i]);
		exp_fit_y.push_back(ec0 + ec1*w_data[i]);
	}
	dbg("Done.\n");
	delete x, y;
}

template <class x_type, class y_type>
void write_file(ofstream *out, const char* filename, vector<x_type> x, vector<y_type> y)
{

	dbg("Saving file: " << filename << "\n");
	out->open(filename);
	
	long len = x.size();
	
	rpt_err(x.size()!=y.size(), "xlen and ylen do not match for file: " << filename);	
	for(long i=0;i<len;i++)
	{
		*out << x[i] << "\t" << y[i] << endl;
	}
	
	out->close();
}

template<class T>
vector<T> grab_interval(long start, long end, vector<T> &arr)
{
	vector<T> new_arr(0);
	
	for(long i=start;i<end;i++)
	{
		new_arr.push_back(arr[i]);
	}
	return new_arr;
}

int main(int argc, char **argv)
{
	string filename = argv[1];
	long len;
	vector<WT> raw_data;
	
	vector<long> arr_n_data(0);
	vector<WT> arr_w_data(0);
	
	vector<double> log_n_data(0);
	vector<double> log_w_data(0);
	
	vector<double> lin_regress_x(0);
	vector<double> lin_regress_y(0);
	
	vector<double> exp_fit_x(0);
	vector<double> exp_fit_y(0);
	
	ofstream *out = new ofstream();

	
	raw_data = get_raw_data(filename);
	
	len = raw_data.size();
	WT total_wealth = 0;
	WT median_wealth = raw_data[len/2];
	
	for(long i=0;i<len;i++)
	{
		total_wealth += raw_data[i];
	}
#if VERBOSE
	cout << "Total Wealth: " << total_wealth << "\nAverage Wealth: " << total_wealth/len << "\nMedian Wealth: " << median_wealth << endl;
#endif 
	populate_data(raw_data, arr_n_data, arr_w_data, log_n_data, log_w_data, lin_regress_x, lin_regress_y, exp_fit_x, exp_fit_y);
	
	write_file<WT, long>(out, "plot/cum_arr.dat", arr_w_data, arr_n_data);
	
	
	
	write_file<double, double>(out, "plot/cum_log.dat", log_w_data, log_n_data);
	write_file<double, double>(out, "plot/lin_act.dat", 
		grab_interval<double>(PERCENTILE*log_w_data.size(), log_w_data.size(), log_w_data), 
		grab_interval<double>(PERCENTILE*log_n_data.size(), log_n_data.size(), log_n_data));

	
	write_file<double, double>(out, "plot/lin_fit.dat", lin_regress_x, lin_regress_y);
	write_file<double, double>(out, "plot/exp_fit.dat", exp_fit_x, exp_fit_y);


	
	
	return 0;
}




//int gsl_fit_linear(const double *x, const size_t xstride, const double *y, const size_t ystride, size_t n, double *c0, double *c1, double *cov00, double *cov01, double *cov11, double *sumsq)
