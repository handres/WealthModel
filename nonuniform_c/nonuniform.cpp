

#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <ctime>

using namespace std;

void record_data(double* ind, int len, string output_file)
{
	ofstream out("data/" + output_file + ".data");
	
	for(int i=0;i<len;i++)
	{
		out << ind[i] << endl;
	}
	
	out.close();
}

void record_data_all(double* data, unsigned int num_individuals, string output_file)
{
	ofstream out("data/" + output_file + ".data");
	
	unsigned int num_points = data[0]/num_individuals;
	
	
	
	for(int i=0;i<num_points-1;i++)
	{
		if(i==0)
		{
			for(int j=1;j<num_individuals;j++)
			{
				out << data[i*num_individuals + j] << ",";
			}
			out << data[i*num_points + num_individuals];
			out << endl;
		}
		else
		{
			for(int j=1;j<num_individuals;j++)
			{
				out << data[i*num_individuals + j] << ",";
			}
			out << data[(i+1)*num_individuals] << endl;
		}

	}
	
	out.close();
}

double* run_model_record(const int n, const double gr, const double dw, const double T, const double b, const double ec, const double pc, const double mp, const double g, const double a)
{
	/*
	 * Evolution Parameters
	 * n: the number of individuals in the market
	 * gr: the growth rate of the market (used to calculate how many distributions to do)
	 * dw: the size of the wealth increment each iteration
	 * T: the amount of time (in years) to run the model
	 * b: the exponent to use in the power function
	 *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * 
	 * Initial Distribution Parameters
	 *
	 * ec: cutoff for the exponential distribution
	 * pc: cutoff for the pareto distribution (given as a percentage)
	 * mp: minimum wage in pareto distribution (given as percentage)
	 * g: gamma parameter for exponential distribution
	 * a: alpha parameter for pareto distribution
	 */
	 
	double* ind = new double[n]; // Individual array
	
	const int exp_cutoff = (int) (ec*n);
	const int power_cutoff = (int) (pc*n);
	const int num_intermediate = power_cutoff - exp_cutoff;
	
	#pragma omp parallel for
	for(int i=1;i<exp_cutoff+1;i++) // Exponential Part of Distribution
	{
		ind[i-1] = - (1/g)*(log(1-(double)(i)/n));
	}
	#pragma omp parallel for
	for(int i=exp_cutoff;i<power_cutoff;i++) // Averaged out Part of Distribution
	{
		double percent = ((double)(i - exp_cutoff))/num_intermediate;
		
		
		ind[i] = (1 - percent) * (-(1/g)*(log(1-(double)(i)/n))) + percent * ind[(int)(mp*n)-1]*pow(1/(1 - ((double)i)/n), 1/a);
		
		//ind[i] = (1 - percent) * (-(1/g)*(log(1-(double)(i)/n))) + percent * pow(a*pow(ind[(int)(mp*n)], a)/(1 - ((double)i)/n), 1/(a+1));
	}
	#pragma omp parallel for
	for(int i=power_cutoff;i<n;i++) // Pareto Part of Distribution
	{
		
		ind[i] = ind[(int)(mp*n)-1]*pow(1/(1 - ((double)i)/n), 1/a);
		
		//ind[i] = pow(a*pow(ind[(int)(mp*n)], a)/(1 - ((double)i)/n), 1/(a+1));
	}
	
	double total_wealth = 0;
	
	for(int i=0;i<n;i++)
	{
		total_wealth +=	ind[i];
	}
	
	
	
	// Initialize the model
	
	const unsigned int years = (int)(T + 1);
	
	const int iter = (int) ((total_wealth)*(exp(gr*T) - 1)/dw + 1); // Calculate how many iterations (distributions) are needed
	
	double* data = new double[years*n];
	
	unsigned int* year_cutoffs = new unsigned int[years];
	
	
	for(int i=0;i<(int)T+1;i++)
	{
		year_cutoffs[i] = (int) ((total_wealth)*(exp(gr*i) - 1)/dw + 1);
	}
	
	#pragma omp parallel for
	for(int i=0;i<n;i++)
	{
		data[i] = ind[i];
	}
	
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> distribution(0, 1); // Create a uniform random number generator
	register double* ind_power = new double[n]; // Create individual power array
	clock_t start; // Use this to benchmark the program 
	double duration;
	start = clock();

	cout << "Iterations Needed: " << iter << endl;
	double mkt_power = 0.0;
	int count = 0;
	
	#pragma omp parallel for
	for(int j=0;j<n;j++) // Initialize the power arrray
	{
		ind_power[j] = pow(data[j], b);
	}

	for(int j=0;j<n;j++) // Calculate the market power
	{
		mkt_power += ind_power[j];
	}
	
	
	const int one_percent = (int) (iter/100); // How many iterations are in 1% of the entire needed amount
	
	
	// Run the model
	
	unsigned short int year_index = 1;
	
	for(int i=1;i<iter;i++) // Main loop
	{
		if(i % one_percent == 0)
		{
			cout << i << "\t\t (" << count << "%)" << endl;
			count++;

		}
		if(i == year_cutoffs[year_index])
		{
			// cout << "year_index: " << year_index << "\tmax: "<< (int)(T+1) << endl;
			for(int j=0;j<n;j++)
			{
				
				data[year_index*n + j] = ind[j];
				
			}
			year_index++;
		}
		const double rand = distribution(gen) * mkt_power;
		double pmf = ind_power[0];
	

		int j = 0;
		while(pmf < rand) // Distribute the wealth increment to an individual.
		{
			// This loop is what slows down the program the most now.
			j++;
			pmf += ind_power[j];	
		}
		ind[j] += dw;
		mkt_power -= ind_power[j];
		ind_power[j] = pow(ind[j], b);
		mkt_power += ind_power[j];
	}
	
	duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
	
	cout << "FINISHED. COMPUTED " << iter << " ITERATIONS IN " << duration << " SECONDS." << endl;
	cout << "IPS: " << iter/duration << endl;
	
	
	
	double* rdata = new double[years*n + 1];
	
	rdata[0] = (int)(T + 1)*n + 1;
	for(int i=0;i<(T+1)*n;i++)
	{
		rdata[i+1] = data[i];
	}
	
	delete ind;
	delete data;
	
	return rdata;

}



double* run_model(const unsigned int n, const double gr, const double dw, const double T, const double b, const double ec, const double pc, const double mp, const double g, const double a)
{
	/*
	 * Evolution Parameters
	 * n: the number of individuals in the market
	 * gr: the growth rate of the market (used to calculate how many distributions to do)
	 * dw: the size of the wealth increment each iteration
	 * T: the amount of time (in years) to run the model
	 * b: the exponent to use in the power function
	 *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * 
	 * Initial Distribution Parameters
	 *
	 * ec: cutoff for the exponential distribution
	 * pc: cutoff for the pareto distribution (given as a percentage)
	 * mp: minimum wage in pareto distribution (given as percentage)
	 * g: gamma parameter for exponential distribution
	 * a: alpha parameter for pareto distribution
	 */
	 
	 
	// Initialize the exponential-Pareto distribution
	
	double* ind = new double[n]; // Individual array
	
	const int exp_cutoff = (int) (ec*n);
	const int power_cutoff = (int) (pc*n);
	const int num_intermediate = power_cutoff - exp_cutoff;
	
	#pragma omp parallel for
	for(int i=1;i<exp_cutoff+1;i++) // Exponential Part of Distribution
	{
		ind[i-1] = - (1/g)*(log(1-(double)(i)/n));
	}
	#pragma omp parallel for
	for(int i=exp_cutoff;i<power_cutoff;i++) // Averaged out Part of Distribution
	{
		double percent = ((double)(i - exp_cutoff))/num_intermediate;
		
		
		ind[i] = (1 - percent) * (-(1/g)*(log(1-(double)(i)/n))) + percent * ind[(int)(mp*n)-1]*pow(1/(1 - ((double)i)/n), 1/a);
		
		//ind[i] = (1 - percent) * (-(1/g)*(log(1-(double)(i)/n))) + percent * pow(a*pow(ind[(int)(mp*n)], a)/(1 - ((double)i)/n), 1/(a+1));
	}
	#pragma omp parallel for
	for(int i=power_cutoff;i<n;i++) // Pareto Part of Distribution
	{
		
		ind[i] = ind[(int)(mp*n)-1]*pow(1/(1 - ((double)i)/n), 1/a);
		
		//ind[i] = pow(a*pow(ind[(int)(mp*n)], a)/(1 - ((double)i)/n), 1/(a+1));
	}
	
	double total_wealth = 0;
	
	for(int i=0;i<n;i++)
	{
		total_wealth +=	ind[i];
	}
	
	
	
	// Initialize the model
	
	const int iter = (int) ((total_wealth)*(exp(gr*T) - 1)/dw + 1); // Calculate how many iterations (distributions) are needed
	
	
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> distribution(0, 1); // Create a uniform random number generator
	register double* ind_power = new double[n]; // Create individual power array
	clock_t start; // Use this to benchmark the program 
	double duration;
	start = clock();

	cout << "Iterations Needed: " << iter << endl;
	double mkt_power = 0.0;
	int count = 0;
	
	#pragma omp parallel for
	for(int j=0;j<n;j++) // Initialize the power arrray
	{
		ind_power[j] = pow(ind[j], b);
	}

	for(int j=0;j<n;j++) // Calculate the market power
	{
		mkt_power += ind_power[j];
	}
	
	
	const int one_percent = (int) (iter/100); // How many iterations are in 1% of the entire needed amount
	
	
	// Run the model
	
	for(int i=0;i<iter;i++) // Main loop
	{
		if(i % one_percent == 0)
		{
			cout << i << "\t\t (" << count << "%)" << endl;
			count++;
		}
		const double rand = distribution(gen) * mkt_power;
		double pmf = ind_power[0];
		

		int j = 0;
		while(pmf < rand) // Distribute the wealth increment to an individual.
		{
			// This loop is what slows down the program the most now.
			j++;
			pmf += ind_power[j];	
		}
		ind[j] += dw;
		mkt_power -= ind_power[j];
		ind_power[j] = pow(ind[j], b);
		mkt_power += ind_power[j];
	}
	
	duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
	
	cout << "FINISHED. COMPUTED " << iter << " ITERATIONS IN " << duration << " SECONDS." << endl;
	cout << "IPS: " << iter/duration << endl;
	
	return ind;
}

int main(int argc, char* argv[])
{
	unsigned int num_individuals;
	double market_growth_rate;
	double delta_omega;
	double time_period;
	double beta; // Power function parameter

	
	double exp_cutoff; // When to start transforming into Pareto distribution
	double pareto_cutoff; // When to stop transforming into Pareto distribution
	double min_percent; // To calculate the "minimum" in the Pareto distribution
	double gamma; // Gamma in Exponential Distribution
	double alpha; // Alpha in Pareto distribution
	
	bool record_all_data;
	
	string output_file = "";
	
	
	if(argc == 13)
	{
		num_individuals = stoi(argv[1]);
		market_growth_rate = stod(argv[2]);
		delta_omega = stod(argv[3]);
		time_period = stod(argv[4]);
		beta = stod(argv[5]);
		
		exp_cutoff = stod(argv[6]);
		pareto_cutoff = stod(argv[7]);
		min_percent = stod(argv[8]);
		gamma = stod(argv[9]);
		alpha = stod(argv[10]);
		
		record_all_data = stoi(argv[11]);
		
		output_file = argv[12];
	}
	else
	{
		num_individuals = 100000;
		market_growth_rate = 0.03;
		delta_omega = 0.01;
		time_period = 0.001;
		beta = 1;
		
		exp_cutoff = .75;
		pareto_cutoff = .9;
		min_percent = 0.699;
		gamma = 1;
		alpha = 1.5;
		
		record_all_data = true;
		
		output_file = "market";
		
	}
	
	cout << "Running Parameters:" << endl << "n=" << num_individuals << endl << "gr=" << market_growth_rate << endl << "dw=" << delta_omega << endl << "T=" << time_period << endl << "b: " << beta << endl;
	
	cout << "Distribution Parameters:" << endl << "ec=" << exp_cutoff << endl << "pc=" << pareto_cutoff << endl << "mp=" << min_percent << endl << "g=" << gamma << endl << "a=" << alpha << endl;
	
	
	if(record_all_data==0)
	{
		cout << "Not recording all data..." << endl;
		double* individuals = run_model(num_individuals, market_growth_rate, delta_omega, time_period, beta, exp_cutoff, pareto_cutoff, min_percent, gamma, alpha);
	
		cout << "Output File: data/" << output_file << ".data" << endl;
	
		record_data(individuals, num_individuals, output_file);
	}
	else
	{
		cout << "Recording all data..." << endl;
		
		double* data = run_model_record(num_individuals, market_growth_rate, delta_omega, time_period, beta, exp_cutoff, pareto_cutoff, min_percent, gamma, alpha);
		cout << "Recording data..." << endl;
		record_data_all(data, num_individuals, output_file);
		delete data;
	}
	
	
	
}



