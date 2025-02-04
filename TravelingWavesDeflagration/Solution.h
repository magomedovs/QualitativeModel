#ifndef SOLUTION
#define SOLUTION

#include <deal.II/base/function.h>

#include "LinearInterpolator.h"

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>


namespace TravelingWave
{
  using namespace dealii;

  // The structure for keeping the solution: arrays of coordinates $\xi$, solution $u$, $T$, $\lambda$, and the wave speed $c$. 
  struct SolutionStruct
  {
    SolutionStruct();
    SolutionStruct(const std::vector<double> &ix, 
                    const std::vector<double> &iT, const std::vector<double> &ilambda, const double iwave_speed);
    SolutionStruct(const std::vector<double> &ix, 
                    const std::vector<double> &iT, const std::vector<double> &ilambda);

    void reinit(const unsigned int number_of_elements);
    
    void save_to_file(std::string filename) const;

    std::vector<double> x;        // mesh coordinates (must be an increasing sequence)
    std::vector<double> T;        // array of T components
    std::vector<double> lambda;   // array of lambda components

    double wave_speed;            // speed of the wave
  };

  // Interpolation class
  class Interpolant : public Function<1>
  {
  public:
    Interpolant(const std::vector<double> &ix_points, const std::vector<double> &iy_points);
    virtual double value(const Point<1> &p, const unsigned int component = 0) const override;

  private:
    LinearInterpolator<double> interpolant;
  };

  // Vector function $(u(p), T(p), \lambda(p))$
  template <typename InterpolantType>
  class SolutionVectorFunction : public Function<1>
  {
  public:
    SolutionVectorFunction(InterpolantType iT_interpolant, InterpolantType ilambda_interpolant);
    virtual double value(const Point<1> &p, const unsigned int component = 0) const override;

  private:
    InterpolantType T_interpolant;
    InterpolantType lambda_interpolant;
  };

} // namespace TravelingWave

#endif
