#include "SolutionStruct.h"

namespace TravelingWave
{
  
	Solution::Solution() {}
	Solution::Solution(const std::vector<double> &ix, const std::vector<double> &iu, 
											const std::vector<double> &iT, const std::vector<double> &ilambda, double iD_0)
	: x(ix)
	, u(iu)
	, T(iT)
	, lambda(ilambda)
	, D_0(iD_0)
	{}
	Solution::Solution(const std::vector<double> &ix, const std::vector<double> &iu, 
											const std::vector<double> &iT, const std::vector<double> &ilambda)
	: Solution(ix, iu, iT, ilambda, 0.)
	{}

	void Solution::reinit(const unsigned int number_of_elements)
	{
		D_0 = 0;
		x.clear();
		u.clear();
		T.clear();
		lambda.clear();

		x.resize(number_of_elements);
		u.resize(number_of_elements);
		T.resize(number_of_elements);
		lambda.resize(number_of_elements);
	}

	void Solution::save_to_file(std::string filename = "sol") const
	{
		const std::string file_for_solution = filename + ".txt";
		std::ofstream output(file_for_solution);

		output << std::scientific << std::setprecision(16);
		for (unsigned int i = 0; i < x.size(); ++i)
		{
			output << std::fixed << x[i]; 
			output << std::scientific << " " << u[i] << " " << T[i] << " " << lambda[i] << "\n";
		}
		output.close();

		// std::ofstream file_for_D_0_output("D_0_value-" + refinement_cycle_string + ".txt");
		std::ofstream file_for_D_0_output("D_0-" + file_for_solution);
		file_for_D_0_output << std::scientific << std::setprecision(16);
		file_for_D_0_output << D_0 << std::endl;
		file_for_D_0_output.close();
	}

}