#include "Parameters.h"

namespace TravelingWave
{
  using namespace dealii;
  
  Problem::Problem() 
    : ParameterAcceptor("Problem")
  {
	  add_parameter("delta", delta = 0.1); 
  	add_parameter("epsilon", epsilon = 0.1);            
  	add_parameter("Prandtl number", Pr = 0.75);        
  	add_parameter("Lewis number", Le = 1.0);      
  	add_parameter("Constant of reaction rate", k = 5.0);
  	add_parameter("Activation energy", theta = 1.65);
  	add_parameter("Heat release", q = 1.7);
    add_parameter("Ignition Temperature", T_ign = 1.0);
    add_parameter("Type of the wave", wave_type = 1);

    add_parameter("Type of Temperature right boundary condition", T_r_bc_type = 1);

    add_parameter("T_left", T_left = 4.13);
    add_parameter("T_right", T_right = 0.9);
    add_parameter("u_left", u_left = 0.);
    add_parameter("u_right", u_right = 0.);

    add_parameter("wave_speed initial value", wave_speed_init = 1.001);
  }

  // void Problem::parse_parameters(ParameterHandler &)
  // {
  //   del_Pr_eps 	= Pr * 4 * delta / (3 * epsilon);
  //   del_Le 		= delta / Le;

  //   // wave_speed_init = ( 1 + epsilon / 2 * (1 + 2 * std::sqrt(q)) ) * coef_for_wave_speed_init;
  // }

  FiniteElements::FiniteElements()
    : ParameterAcceptor("Finite elements")
  {
    add_parameter("Polynomial degree", poly_degree = 1);
    add_parameter("Number of quadrature points", quadrature_points_number = 2);
  }

  Mesh::Mesh()
	  : ParameterAcceptor("Mesh")
  {
	  add_parameter("Interval left boundary", interval_left = -7.0);  
  	add_parameter("Interval right boundary", interval_right = 3.0);
  	add_parameter<unsigned int>("Refinements number", refinements_number = 8);
   	add_parameter("Adaptive mesh refinement", adaptive = 1);
  }

  Solver::Solver()
    : ParameterAcceptor("Solver")
  {
    add_parameter("Tolerance", tol = 1e-8);
    add_parameter("Max iterations", max_iter = 1000);
  }

  Output::Output()
    : ParameterAcceptor("Output")
  {
    add_parameter("Log verbosity", verbosity = 0);
  }

} // TravelingWave
