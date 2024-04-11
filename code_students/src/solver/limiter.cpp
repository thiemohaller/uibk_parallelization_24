#include "solver/limiter.hpp"

#include <algorithm>
#include <iostream>
#include <cmath>

limiter_base::limiter_base() {}

limiter_minmod::limiter_minmod(double theta) { this->theta = theta; }

double limiter_minmod::compute(double first, double second, double third) {
	double minmod = 0.0;
	bool sign_is_same = first * second > 0.0 && second * third > 0.0;

	if (sign_is_same) {
		double abs_first = std::abs(first);
		double abs_second = std::abs(second);
		double abs_third = std::abs(third);
		
		// take minimum value
		double min_value = std::min({abs_first, abs_second, abs_third});
		minmod = std::copysign(min_value, first);
	}

	return minmod;
}