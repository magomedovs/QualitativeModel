#ifndef SOLUTION_STRUCT
#define SOLUTION_STRUCT

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

namespace TravelingWave
{

	struct Solution
	{
		Solution();
		Solution(const std::vector<double> &ix, const std::vector<double> &iu, 
							const std::vector<double> &iT, const std::vector<double> &ilambda, const double iD_0);
		Solution(const std::vector<double> &ix, const std::vector<double> &iu, 
							const std::vector<double> &iT, const std::vector<double> &ilambda);

		void reinit(const unsigned int number_of_elements);
		
		void save_to_file(std::string filename) const;

		std::vector<double> x;
		std::vector<double> u;
		std::vector<double> T;
		std::vector<double> lambda;

		double D_0;
	};

}

#endif