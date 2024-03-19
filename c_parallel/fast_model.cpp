

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
		ind[i] = w_i;
	}
	
	
	
	const unsigned long iter = (long) ((n*w_i)*(exp(a*T) - 1)/dw + 1); // Calculate how many iterations (distributions) are needed
	
	
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> distribution(0, 1); // Create a uniform random number generator
	register double* ind_power = new double[n]; // Create individual power array
	


	cout << "Iterations Needed: " << iter << endl;
	double mkt_power = 0.0;
	int count = 0;
	
	
	
	double* cum_power = new double[n];
	
	
	
	#pragma omp parallel for
	for(int j=0;j<n;j++) // Initialize the power arrray
	{
		ind_power[j] = pow(ind[j], beta);
	}

	for(int j=0;j<n;j++) // Calculate the market power
	{
		mkt_power += ind_power[j];
		cum_power[j] = mkt_power;
	}
	
	
	const int one_percent = (int) (iter/100); // How many iterations are in 1% of the entire needed amount
	
	// Benchmark Parameters
	clock_t start;
	double duration;
	start = clock();
	
	
	clock_t bin_search_clock;
	double bin_search_duration=0;
	
	clock_t update_clock;
	double update_duration=0;
	
	clock_t rand_gen_clock;
	double rand_gen_duration=0;
	
	for(long i=0;i<iter;i++) // Main loop
	{
		if(i % one_percent == 0)
		{
			cout << i << "\t\t (" << count << "%)" << endl;
			count++;
		}
		
		rand_gen_clock=clock();
		const double rand = distribution(gen) * mkt_power;
		rand_gen_duration += (clock() - rand_gen_clock) / (double) CLOCKS_PER_SEC;

		// Perform binary search
		bin_search_clock = clock();
		unsigned int ci;
		unsigned int high = n;
		unsigned int low = 0;
		if (cum_power[0] >= rand)
		{
			ci = 0;
		}
		else if (cum_power[n-1] <= rand)
		{
			ci = n-1;
		}
		else {

			while(high - low > 1)
			{
				if(cum_power[high - (high-low)/2]>=rand)
				{
					high -= (high-low)/2;
				}
			
				else if(cum_power[low + (high-low)/2 < rand])
				{
					low += (high-low)/2;
				}
			}
			ci = high;
		}
		
		bin_search_duration += (clock() - bin_search_clock) / (double) CLOCKS_PER_SEC;
		
		// Update variables
		ind[ci] += dw;
		double temp = ind_power[ci];
		
		ind_power[ci] = pow(ind[ci], beta);
		double change_in_power = ind_power[ci] - temp;
		mkt_power += change_in_power;
		
		// Update search array (this is the bottleneck)
		update_clock = clock();
		#pragma omp for
		for(int j=ci;j<n;j++)
		{
			cum_power[j] += change_in_power;
		}
		update_duration += (clock() - update_clock) / (double) CLOCKS_PER_SEC;

	}
	
	cout << "Time for random number generation: " << rand_gen_duration << endl;
	cout << "Time for binary search: " << bin_search_duration << endl;
	cout << "Time for cumulative power update: " << update_duration << endl;
	
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
	int market_type; // Deprecated variable - has no effect on program
	string output_file = "";
	double exponent_power;

	if(argc == 9)
	{
		num_individuals = stoi(argv[1]);
		init_worth_per_individual = stod(argv[2]);
		market_growth_rate = stod(argv[3]);
		delta_omega = stod(argv[4]);
		time_period = stod(argv[5]);
		market_type = stoi(argv[6]); // Deprecated variable - has no effect on program
		exponent_power = stod(argv[7]);
		output_file = argv[8];
	}
	else
	{
		num_individuals = 1000;
		init_worth_per_individual = 1;
		market_growth_rate = 0.03;
		delta_omega = 0.01;
		time_period = 5;
		market_type = 0; // 0 = linear, 1 = exponential, 2 = linear-exponential
		output_file = "data/market.data";
		exponent_power = 2;
	}
	cout << "Parameters:" << endl << "W_i=" << num_individuals*init_worth_per_individual << endl << "N=" << num_individuals << endl << "alpha=" << market_growth_rate << endl << "T=" << time_period << endl << "Exponent Power: " << exponent_power << endl;
	
	double* individuals = run_model(num_individuals, init_worth_per_individual, market_growth_rate, delta_omega, time_period, exponent_power);
	
	record_data(individuals, num_individuals, output_file);
}



