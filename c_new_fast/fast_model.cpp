

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


double* run_model(const int n, double w_i, const double a, const double dw, const double T, const double beta)
{
	/*
	 *
	 * n: the number of individuals in the market
	 * w_i: the initial wealth per individual in the market
	 * a: the growth rate of the market (used to calculate how many distributions to do)
	 * dw: the size of the wealth increment each iteration
	 * T: the amount of time (in years) to run the model
	 * beta: the exponent to use in the power function
	 *
	 */
	double* ind = new double[n]; // Individual array
	
	for(int i=0;i<n;i++) // Initialize the individual array
	{
		
	}
	
	const int iter = (int) ((n*w_i)*(exp(a*T) - 1)/dw + 1); // Calculate how many iterations (distributions) are needed
	
	
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> distribution(0, 1); // Create a uniform random number generator
	register double* ind_power = new double[n]; // Create individual power array
	clock_t start; // Use this to benchmark the program 
	double duration;
	start = clock();

	cout << "Iterations Needed: " << iter << endl;
	double mkt_worth = 0.0;
	double mkt_power = 0.0;
	int count = 0;
	
	for(int j=0;j<n;j++) // Initialize the power arrray
	{
		ind[j] = w_i;
		mkt_worth += w_i;
		
	}

	for(int j=0;j<n;j++) // Calculate the market power
	{
		ind_power[j] = pow(ind[j]/mkt_worth, beta);
		mkt_power += ind_power[j];
	}
	
	
	const int one_percent = (int) (iter/100); // How many iterations are in 1% of the entire needed amount
	
	for(int i=0;i<iter;i++) // Main loop
	{
		if(i % one_percent == 0)
		{
			cout << i << "\t\t (" << count << "%)" << endl;
			count++;
		}
		
		const double rand = distribution(gen) * mkt_power;
		double pmf = 0;
		mkt_worth += dw;
		mkt_power = 0;
		bool distributed = false;
		for(int j=0;j<n;j++) // Distribute wealth and update everyone's power
		{
			pmf += ind_power[j];
			if(pmf >= rand && !distributed)
			{
				ind[j] += dw;
				distributed = true;
			}
			
			ind_power[j] = pow(ind[j]/mkt_worth, beta);
			mkt_power += ind_power[j];
		}
	}
	
	duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
	
	cout << "FINISHED. COMPUTED " << iter << " ITERATIONS IN " << duration << " SECONDS." << endl;
	cout << "IPS: " << iter/duration << endl;
	
	return ind;
}

int main(int argc, char* argv[])
{
	int num_individuals;
	double init_worth_per_individual;
	double market_growth_rate;
	double delta_omega;
	double time_period;
	int market_type; // 0 = linear, 1 = exponential, 2 = linear-exponential
	string output_file = "";
	double exponent_power;

	if(argc == 9)
	{
		num_individuals = stoi(argv[1]);
		init_worth_per_individual = stod(argv[2]);
		market_growth_rate = stod(argv[3]);
		delta_omega = stod(argv[4]);
		time_period = stod(argv[5]);
		market_type = stoi(argv[6]);
		exponent_power = stod(argv[7]);
		output_file = argv[8];
	}
	else
	{
		num_individuals = 1000;
		init_worth_per_individual = 100;
		market_growth_rate = 0.03;
		delta_omega = 0.1;
		time_period = 5;
		market_type = 0; // 0 = linear, 1 = exponential, 2 = linear-exponential
		output_file = "data/market.data";
		exponent_power = 2;
	}
	cout << "Parameters:" << endl << "W_i=" << num_individuals*init_worth_per_individual << endl << "N=" << num_individuals << endl << "alpha=" << market_growth_rate << endl << "T=" << time_period << endl << "Exponent Power: " << exponent_power << endl;
	
	double* individuals = run_model(num_individuals, init_worth_per_individual, market_growth_rate, delta_omega, time_period, exponent_power);
	
	record_data(individuals, num_individuals, output_file);
}



