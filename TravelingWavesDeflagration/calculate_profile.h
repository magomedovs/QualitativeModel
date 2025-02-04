#ifndef CALCULATE_PROFILE
#define CALCULATE_PROFILE

#include "Parameters.h"
#include "Solution.h"
#include "AuxiliaryFunctions.h"

namespace TravelingWave
{
	void compute_initial_guess_deflagration(Parameters &params, SolutionStruct &initial_guess);

  void calculate_profile(Parameters& parameters);

} // namespace TravelingWave

#endif
