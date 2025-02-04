#include "TravelingWaveSolver.h"
#include "calculate_profile.h"

namespace TravelingWave
{
  using namespace dealii;


  // Construction of a piecewise constant initial guess for deflagration wave solution.
  void compute_initial_guess_deflagration(Parameters &params, SolutionStruct &initial_guess)
  {
    const Problem &problem = params.problem;
    double current_wave_speed(problem.wave_speed_init);

    double u_0l = 0.;

    if (problem.epsilon <= std::pow(10., -5))
    {
      // double kappa = 2. / problem.epsilon + 1. + 2. * problem.u_right;
      double kappa = (2. + problem.epsilon * (1. + 2. * problem.u_right)) / problem.epsilon;
      double u_0l_approx = problem.u_right - problem.q / kappa - std::pow(problem.q, 2) / std::pow(kappa, 3);   // Formula obtained using Taylor expansion
      u_0l = u_0l_approx;

      std::cout << "Computation with approximate formula for u_0l." << std::endl;
      std::cout << std::setprecision(16) << "u_0l_approx = " << u_0l_approx << std::endl;
    }
    else
    {
      u_0l = 1./2. * (-2./problem.epsilon - 1. + std::sqrt(std::pow(-2./problem.epsilon - 1. - 2. * problem.u_right, 2) - 4. * problem.q));

      std::cout << "Computation with exact formula for u_0l." << std::endl;
      std::cout << std::setprecision(16) << "u_0l = " << u_0l << std::endl;
    }
    
    params.problem.T_right = problem.T_left + problem.u_right - problem.q - u_0l;

    // std::cout << "u_0l = " << u_0l << std::endl;
    // std::cout << "T_right = " << problem.T_right << std::endl;

    auto T_init_guess_func = [&](double x) {
      if (x < 0.)
      {
        return (problem.T_left - problem.T_ign) * (1 - std::exp(current_wave_speed * x)) + problem.T_ign;
      }
      else if (isapprox(x, 0.))
      {
        return problem.T_ign;
      }
      else
      {
        return (problem.T_ign - problem.T_right) * std::exp(-current_wave_speed * x) + problem.T_right;
      }
    };

    auto lambda_init_guess_func = [=](double x) {
      double Lambda_0_at_zero_guess = 0.38;

      if (x < 0.)
      {
        // return 1 - std::exp(current_wave_speed * x);
        return (1 - Lambda_0_at_zero_guess) * (1 - std::exp(current_wave_speed * x)) + Lambda_0_at_zero_guess;
      }
      else 
      {
        // return 0.;
        return Lambda_0_at_zero_guess * std::exp(-current_wave_speed * x);
      }
    };

    unsigned int multiplier_for_number_of_points = 7;
    unsigned int number_of_points = multiplier_for_number_of_points * static_cast<unsigned int>(std::trunc(std::abs( params.mesh.interval_right - params.mesh.interval_left )));
    std::vector<double> x_span(number_of_points);
    linspace(params.mesh.interval_left, params.mesh.interval_right, x_span);

    std::vector<double> T_init_arr(number_of_points);
    std::vector<double> lambda_init_arr(number_of_points);

    for (unsigned int i = 0; i < number_of_points; ++i)
    {
      T_init_arr[i] = T_init_guess_func(x_span[i]);
      lambda_init_arr[i] = lambda_init_guess_func(x_span[i]);
    }

    initial_guess.x = x_span;
    initial_guess.T = T_init_arr;
    initial_guess.lambda = lambda_init_arr;
    initial_guess.wave_speed = current_wave_speed;

  }


  // Compute the traveling-wave profile. The continuation method can be switched on by setting the argument <code> continuation_for_delta </code> as <code> true </code>.
  void calculate_profile(Parameters& parameters)
  {
    SolutionStruct sol;

    compute_initial_guess_deflagration(parameters, sol);
    // sol.save_to_file("init_sol");

    TravelingWaveSolver wave(parameters, sol);
    std::string filename = "solution_defl_inner_eps-" + Utilities::to_string(parameters.problem.epsilon);
    wave.run(filename);
    wave.get_solution(sol);
  }

} // namespace TravelingWave
