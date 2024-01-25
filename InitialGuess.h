#ifndef INITIAL_GUESS
#define INITIAL_GUESS

#include <deal.II/base/function.h>

#include "linspace.h"
#include "LinearInterpolator.h"
#include "Parameters.h"
#include "SolutionStruct.h"
#include "LimitSolution.h"
#include "IntegrateSystem.h"
#include "are_equal.h"

namespace TravelingWave
{
	using namespace dealii;

	/* Linear interpolation class */
	class Interpolant : public Function<1>
	{
	public:
		Interpolant(const std::vector<double> &ix_points, const std::vector<double> &iy_points);
		virtual double value(const Point<1> &p, const unsigned int component = 0) const override;

	private:
		LinearInterpolator<double> interpolant;
	};


	template <typename FunctorType>
	class InitGuessFunc : public Function<1>
	{
	public:
		InitGuessFunc(FunctorType iu_interpolant, FunctorType iT_interpolant, FunctorType ilambda_interpolant);
		virtual double value(const Point<1> &p, const unsigned int component = 0) const override;

	private:
		FunctorType u_interpolant;
		FunctorType T_interpolant;
		FunctorType lambda_interpolant;
	};


  void compute_limit_sol_left_part(const Parameters &parameters, 
                                    const double D_0, 
                                    const double u_0, 
                                    const double T_0, 
                                    const double lambda_0, 
                                    Solution &LimitSol, 
                                    const double root_sign = 1.);

	void compute_initial_guess_ode(const Parameters &params, Solution &initial_guess, const double root_sign = 1.);
	void compute_initial_guess_slow(const Parameters &params, Solution &initial_guess);


} // namespace TravelingWave

#endif