#include "WaveConstructor.h"

void calculate_and_save_solution(TravelingWave::Parameters& parameters,
 																	const bool continuation_for_delta=false, 
																	const double delta_start=0.01, 
																	const unsigned int number_of_continuation_points=10)
{
	using namespace TravelingWave;

	Solution sol;

	if (parameters.problem.wave_type == 1) 				/* fast wave (Detonation, fast deflagration (root_sign=1) and supersonic compression wave (root_sign=-1).) */
	{
		compute_initial_guess_ode(parameters, sol);
		// compute_initial_guess_ode(parameters, sol, -1.);
	}
	else if (parameters.problem.wave_type == 0)		/* slow wave (Slow deflagration) */
	{
		compute_initial_guess_slow(parameters, sol);
	}
	
	if (continuation_for_delta == false)
	{
		WaveConstructor wave(parameters, sol);
		std::string filename = "solution_delta-" + Utilities::to_string(parameters.problem.delta) + "_eps-" 
																											+ Utilities::to_string(parameters.problem.epsilon);
		wave.run(parameters.mesh.adaptive, parameters.mesh.refinements_number, parameters.solver.tol, filename);
		wave.get_solution(sol);
	}
	else	/* Run with continuation_for_delta. */
	{
		double delta_target = parameters.problem.delta;
		parameters.problem.delta = delta_start;

		std::vector<double> delta_span(number_of_continuation_points);

		/* Generate a sequence of delta values being uniformly distributed in log10 scale. */
		{
			double delta_start_log10 = std::log10(delta_start);
			double delta_target_log10 = std::log10(delta_target);

			std::vector<double> delta_log_span(delta_span.size());
			linspace(delta_start_log10, delta_target_log10, delta_log_span);

			for (unsigned int i = 0; i < delta_span.size(); ++i)
			{
				delta_span[i] = std::pow(10, delta_log_span[i]);
			}
		}

		std::vector<double> wave_speed_array;
		wave_speed_array.reserve(delta_span.size());
		Triangulation<1> tria;
		bool first_iter_flag = true;

		for (double delta : delta_span)
		{
			parameters.problem.delta = delta;
			std::string filename = "solution_delta-" + Utilities::to_string(parameters.problem.delta) + "_eps-" 
																									+ Utilities::to_string(parameters.problem.epsilon);

			WaveConstructor wave(parameters, sol);

			if (first_iter_flag)
			{
				first_iter_flag = false;
			}
			else
			{
				wave.set_triangulation(tria);	
			}

			wave.run(parameters.mesh.adaptive, parameters.mesh.refinements_number, parameters.solver.tol, filename);
			wave.get_solution(sol);
			wave.get_triangulation(tria);
			wave_speed_array.push_back(sol.wave_speed);

		}

		{
			const std::string file_for_solution = "delta_wave_speed.txt";
			std::ofstream output(file_for_solution);

			output << std::scientific << std::setprecision(16);
			for (unsigned int i = 0; i < delta_span.size(); ++i)
			{
				output << std::scientific << delta_span[i] << " " << wave_speed_array[i] << "\n";
			}
			output.close();
		}
	}

	/* Boundary conditions check. */
	{
		unsigned int sol_length = sol.x.size();
		double u_r = sol.u[sol_length-1];
		double T_r = sol.T[sol_length-1];
		double u_l = sol.u[0];
		double T_l = sol.T[0];
		double wave_speed = sol.wave_speed;

		double residual_1 = ( u_r * (1 - wave_speed) + parameters.problem.epsilon / 2 * (u_r * u_r + T_r) ) - ( u_l * (1 - wave_speed) + parameters.problem.epsilon / 2 * (u_l * u_l + T_l) );
		double residual_2 = (T_r - u_r + parameters.problem.q) - (T_l - u_l);

		std::cout << "bc residual_1 = " << residual_1 << std::endl;
		std::cout << "bc residual_2 = " << residual_2 << std::endl;
	}

}


int main(int argc, char *argv[])
{

	try
	{
		using namespace TravelingWave;
	
		Parameters parameters;

		const std::string prm_file = argc > 1 ? argv[1] : "../ParametersList.prm";

		std::cout << "Reading parameters... " << std::flush;
		ParameterAcceptor::initialize(prm_file);
		std::cout << "done" << std::endl;
		
		// calculate_and_save_solution(parameters, true, 0.1, 13);	/* With continuation_for_delta */
		calculate_and_save_solution(parameters);
		

	}
	catch (std::exception &exc)
	{
		std::cerr << std::endl
							<< std::endl
							<< "----------------------------------------------------"
							<< std::endl;
		std::cerr << "Exception on processing: " << std::endl
							<< exc.what() << std::endl
							<< "Aborting!" << std::endl
							<< "----------------------------------------------------"
							<< std::endl;
		return 1;
	}
	catch (...)
	{
		std::cerr << std::endl
							<< std::endl
							<< "----------------------------------------------------"
							<< std::endl;
		std::cerr << "Unknown exception!" << std::endl
							<< "Aborting!" << std::endl
							<< "----------------------------------------------------"
							<< std::endl;
		return 1;
	}

	return 0;
}
