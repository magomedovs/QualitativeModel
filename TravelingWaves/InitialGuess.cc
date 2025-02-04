#include "InitialGuess.h"

namespace TravelingWave
{
	using namespace dealii;

	Interpolant::Interpolant(const std::vector<double> &ix_points, const std::vector<double> &iy_points)
	: interpolant(ix_points, iy_points)
	{}


	double Interpolant::value(const Point<1> &p, const unsigned int component) const
	{
		double x = p[0];
		double res = interpolant(x);

		return res;
	}


	template <typename FunctorType>
	InitGuessFunc<FunctorType>::InitGuessFunc(FunctorType iu_interpolant, FunctorType iT_interpolant, FunctorType ilambda_interpolant) 
	: Function<1>(3)
	, u_interpolant(iu_interpolant)
	, T_interpolant(iT_interpolant)
	, lambda_interpolant(ilambda_interpolant)
	{}


	template <typename FunctorType>
	double InitGuessFunc<FunctorType>::value(const Point<1> &p, const unsigned int component) const
	{
		double res = 0.;
		if (component == 0) 			{	res = u_interpolant.value(p); }
		else if (component == 1) 	{ res = T_interpolant.value(p);	}
		else if (component == 2)	{	res = lambda_interpolant.value(p);	}

		return res;
	}


	void compute_limit_sol_left_part(const Parameters &parameters, 
																		const double wave_speed, 
																		const double u_0, 
																		const double T_0, 
																		const double lambda_0, 
																		Solution &LimitSol, 
																		const double root_sign)
	{
		/* Computation of the limit case solution (corresponding to \delta = 0). */

		LimitSolution limit_sol(parameters, lambda_0, u_0, T_0, root_sign);
		limit_sol.set_wave_speed(wave_speed);
		
		{
			/* We take more integration points to better resolve the transition layer. */
			std::vector<double> t_span(static_cast<unsigned int>(std::abs( 0. - parameters.mesh.interval_left )));
			double finer_mesh_starting_value = -9.1;
			linspace(parameters.mesh.interval_left, finer_mesh_starting_value, t_span); 
			std::vector<double> fine_grid(10000);
			linspace(finer_mesh_starting_value + 1e-4, 0., fine_grid); 
			t_span.insert(t_span.end(), fine_grid.begin(), fine_grid.end());

			/* Reverse the order of the elements (because we need to perform back in time integration). */
			std::reverse(t_span.begin(), t_span.end());

			state_type lambda_val(1);
			lambda_val[0] = lambda_0; 		/* initial value */
			IntegrateSystemAtTimePoints(limit_sol.lambda_vec, limit_sol.t_vec, t_span,
				limit_sol, 
				lambda_val, //0., params.mesh.interval_left, 
				-1e-2, Integrator_Type::dopri5);
		}

		limit_sol.calculate_u_T_omega();

		std::reverse(limit_sol.t_vec.begin(), limit_sol.t_vec.end()); /* Reverse the order of elements */
		std::reverse(limit_sol.lambda_vec.begin(), limit_sol.lambda_vec.end()); /* Reverse the order of elements */
		std::reverse(limit_sol.u_vec.begin(), limit_sol.u_vec.end()); /* Reverse the order of elements */
		std::reverse(limit_sol.T_vec.begin(), limit_sol.T_vec.end()); /* Reverse the order of elements */
		std::reverse(limit_sol.omega_vec.begin(), limit_sol.omega_vec.end()); /* Reverse the order of elements */

		SaveSolutionIntoFile(limit_sol.lambda_vec, limit_sol.t_vec, "solution_lambda_limit.txt");
		SaveSolutionIntoFile(limit_sol.u_vec, limit_sol.t_vec, "solution_u_limit.txt");
		SaveSolutionIntoFile(limit_sol.T_vec, limit_sol.t_vec, "solution_T_limit.txt");
		SaveSolutionIntoFile(limit_sol.omega_vec, limit_sol.t_vec, "solution_omega_limit.txt");

		LimitSol.reinit(limit_sol.t_vec.size());
		LimitSol.wave_speed = wave_speed;
		for (unsigned int i=0; i < limit_sol.t_vec.size(); ++i)
		{
			LimitSol.x[i] = limit_sol.t_vec[i];
			LimitSol.u[i] = limit_sol.u_vec[i][0];
			LimitSol.T[i] = limit_sol.T_vec[i][0];
			LimitSol.lambda[i] = limit_sol.lambda_vec[i][0];
		}
	}


	/* Compute initial guess for Newton's iteration. The boundary values are also set appropriately inside the function. */
	void compute_initial_guess_ode(const Parameters &params, Solution &initial_guess, const double root_sign)
	{
		const Problem &problem = params.problem;
		double current_wave_speed(problem.wave_speed_init);

		double T_plus = 0.;
		if (problem.T_left > problem.T_ign + problem.q + std::sqrt(problem.q))  /* Detonation case */
		{
			T_plus = problem.T_right;
		}
		else 																																		/* Fast deflagration or supersonic compression wave case */
		{
			T_plus = problem.T_ign;
		}

		{
			double DeltaT = problem.T_left - T_plus;
			double qDT = problem.q - DeltaT;
			current_wave_speed = 1. + problem.epsilon * (problem.u_right - (qDT * qDT + DeltaT) / (2 * qDT));
		}

		double lambda_0 = 0.;
		double u_0 = problem.u_right;
		double T_0 = T_plus;

		compute_limit_sol_left_part(params, current_wave_speed, u_0, T_0, lambda_0, initial_guess, root_sign);

		initial_guess.wave_speed = current_wave_speed;

		for (int i = initial_guess.x.size() - 1; i > - 1; --i)
		{
			if (are_equal(initial_guess.x[i], 0.))
			{
				initial_guess.u[i] = problem.u_right;
				initial_guess.T[i] = problem.T_ign;
				initial_guess.lambda[i] = 0.;
				break;
			}
		}

		unsigned int number_of_additional_points = 5; 	/* Adding the points to the right part of the interval (w.r.t. \xi = 0). */
		for (unsigned int i = 0; i < number_of_additional_points; ++i)
		{
			initial_guess.x.push_back(params.mesh.interval_right / (std::pow(2., number_of_additional_points - 1 - i)));
			initial_guess.u.push_back(problem.u_right);
			initial_guess.T.push_back(problem.T_right);
			initial_guess.lambda.push_back(0.);
		}

	}


	void compute_initial_guess_slow(const Parameters &params, Solution &initial_guess)
	{
		const Problem &problem = params.problem;
		double current_wave_speed(problem.wave_speed_init);

		double del_Pr_eps = (problem.Pr * 4 * problem.delta / (3 * problem.epsilon));
		double del_Le = (problem.delta / problem.Le);

		auto u_init_guess_func = [&](double x) {
			if (x < 0.)
			{
				return problem.u_left;
			}
			else
			{
				return problem.u_right;
			}
		};

		auto T_init_guess_func = [&](double x) {
			if (x < 0.)
			{
				return problem.T_left;
			}
			else if (are_equal(x, 0.))
			{
				return problem.T_ign;
			}
			else
			{
				return problem.T_right;
			}
		};

		auto lambda_init_guess_func = [=](double x) {
			if (x < 0.)
			{
				return -std::exp(x * std::abs(1 - current_wave_speed) / del_Pr_eps) + 1;
			}
			else 
			{
				return 0.;
			}
		};

		unsigned int multiplier_for_number_of_points = 7;
		unsigned int number_of_points = multiplier_for_number_of_points * static_cast<unsigned int>(std::trunc(std::abs( params.mesh.interval_right - params.mesh.interval_left )));
		std::vector<double> x_span(number_of_points);
		linspace(params.mesh.interval_left, params.mesh.interval_right, x_span);

		std::vector<double> u_init_arr(number_of_points);
		std::vector<double> T_init_arr(number_of_points);
		std::vector<double> lambda_init_arr(number_of_points);

		for (unsigned int i = 0; i < number_of_points; ++i)
		{
			u_init_arr[i] = u_init_guess_func(x_span[i]);
			T_init_arr[i] = T_init_guess_func(x_span[i]);
			lambda_init_arr[i] = lambda_init_guess_func(x_span[i]);
		}

		initial_guess.x = x_span;
		initial_guess.u = u_init_arr;
		initial_guess.T = T_init_arr;
		initial_guess.lambda = lambda_init_arr;
		initial_guess.wave_speed = current_wave_speed;

	}


	template class InitGuessFunc<Interpolant>;

} // namespace TravelingWave