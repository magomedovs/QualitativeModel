#ifndef WAVE_CONSTRUCTOR
#define WAVE_CONSTRUCTOR

#include <deal.II/base/timer.h>
#include <deal.II/base/function.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/sundials/kinsol.h>

#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>

#include "Parameters.h"
#include "SolutionStruct.h"
#include "InitialGuess.h"
#include "are_equal.h"

#include <cmath>
#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <set>

namespace TravelingWave
{
	using namespace dealii;

	class PhiFunction : public Function<1>
	{
	public:
		PhiFunction(double iinterval_start, double iinterval_end);
		virtual double value(const Point<1> &p, const unsigned int component = 0) const override;
		virtual Tensor<1, 1> gradient(const Point<1> &p, const unsigned int component = 0) const override;

	private:
		double interval_start;
		double interval_end;
		const std::vector<double> f{0., 0.5, 1.};
		boost::math::interpolators::cardinal_cubic_b_spline<double> spline;
	};


	class WaveConstructor
	{
	public:
		WaveConstructor(const Parameters &parameters, const Solution &initial_guess_input);
		void set_triangulation(const Triangulation<1> &itriangulation);
		WaveConstructor(const Parameters &parameters, const Solution &initial_guess_input, const Triangulation<1> &itriangulation);

		void run(const int mesh_refinement_type=1, 
							const unsigned int n_refinements=3, 
							const double tol=1e-5, 
							const std::string filename="solution", 
							const bool save_solution_to_file=true);
		void get_solution(Solution &solution) const;
		void get_triangulation(Triangulation<1> &otriangulation) const;

	private:
		void setup_system(const bool initial_step);
		void find_boundary_and_zero_dof_numbers();
		void set_boundary_and_zero_values();

		void set_initial_guess();

		const double steepness_const = 55.;
		double f_omega_mult(double x) const;
		double df_omega_mult(double x) const;
		
		void compute_and_factorize_jacobian(const Vector<double> &evaluation_point_extended);
		double compute_residual(const Vector<double> &evaluation_point_extended, Vector<double> &residual);
		void split_extended_solution_vector();

		void solve(const Vector<double> &rhs, Vector<double> &solution, const double /*tolerance*/);
		void refine_mesh();
		double run_newton_iterations(const double target_tolerance=1e-5);

		void output_results_txt(const Vector<double> &solution, const double wave_speed, const std::string filename="solution");

		unsigned int extended_solution_dim;
		std::map<std::string, unsigned int> boundary_and_zero_dof_numbers;

		const Parameters &params;
		const Problem		 &problem;

		PhiFunction phiFunction;

		unsigned int number_of_quadrature_points;

		Triangulation<1> triangulation;
		bool 						 triangulation_uploaded;
		FESystem<1>      fe;
		DoFHandler<1>    dof_handler;

		AffineConstraints<double> zero_boundary_constraints;	/* Constraints for boundary conditions. Non-homogeneous, if we start from zero initial guess; homogeneous, if we set appropriate boundary conditions in the initial guess. */

		SparsityPattern		  									sparsity_pattern_extended;
		SparseMatrix<double>  								jacobian_matrix_extended;
	  std::unique_ptr<SparseDirectUMFPACK> 	jacobian_matrix_extended_factorization;


		Vector<double> 	current_solution;
		double 					current_wave_speed;
		Vector<double> 	current_solution_extended;		/* Solution with an additional term, corresponding to the velocity variable wave_speed. */

		Solution initial_guess;

		TimerOutput computing_timer;

	};

} // namespace TravelingWave

#endif