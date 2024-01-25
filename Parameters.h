#ifndef PARAMETERS
#define PARAMETERS

#include <deal.II/base/parameter_acceptor.h>

namespace TravelingWave
{
  using namespace dealii;

  struct Problem : ParameterAcceptor
  {
    Problem();

    // void parse_parameters(ParameterHandler &prm) override;

    double delta, epsilon;
    double Pr, Le;
    double k, theta, q;
    double T_ign;
    int wave_type;
    int T_r_bc_type;
    int u_r_bc_type;
    double T_left, T_right;
    double u_left, u_right;

    // double del_Pr_eps;
    // double del_Le;
    
    double D_0_init;
  };

  struct FiniteElements : ParameterAcceptor
  {
	  FiniteElements();

    unsigned int poly_degree;
    unsigned int quadrature_points_number;
  };

  struct Mesh : ParameterAcceptor
  {
    Mesh();
    
    double interval_left;
    double interval_right;
    unsigned int refinements_number;
    int adaptive;
  };

  struct Solver : ParameterAcceptor
  {
    Solver();

    double       tol;
    unsigned int max_iter;
  };

  struct Output : ParameterAcceptor
  {
    Output();

    unsigned int verbosity;
  };

  struct Parameters
  {
    Problem           problem;
    FiniteElements 		fe;
    Mesh			 	      mesh;
    Solver   		      solver;
    Output         		output;
  };

} // TravelingWave

#endif
