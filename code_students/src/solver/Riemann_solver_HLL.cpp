#include "solver/Riemann_solvers.hpp"

HLL_solver::HLL_solver(std::size_t num_fields_) : Riemann_solver(num_fields_) { epsilon = 1.e-50; }

// * HLL solver approximates the Riemann solution by using the slowest and fastest characteristic speeds (slowest and fastest waves)
// * The slowest and fastest characteristic speeds are the minimum and maximum eigenvalues of the Jacobian matrix
void HLL_solver::get_num_flux(fluid_cell &fluid_left_cell, fluid_cell &fluid_right_cell, const std::vector<double> &phys_flux_left_cell,
                              const std::vector<double> &phys_flux_right_cell, std::vector<double> &num_flux, double v_char_slowest, double v_char_fastest) {
	// Apply HLL fluxes
	// * S_L > 0 => this is the slowest? wave, if greater than 0, all waves flow to the right -> left cell is the fastest? and we copy the flux from the left cell
	if (v_char_slowest > 0.0) {
		// 
		for (std::size_t i_field = 0; i_field < num_fields; i_field++) {
			num_flux[i_field] = phys_flux_left_cell[i_field];
		}
	// * S_R < 0 => this is the fastest? wave, if less than 0, all waves flow to the left -> right cell is the fastest? and we copy the flux from the right cell
	} else if (v_char_fastest < 0.0) {
		for (std::size_t i_field = 0; i_field < num_fields; i_field++) {
			num_flux[i_field] = phys_flux_right_cell[i_field];
		}
	// * some waves are moving to the left, some are moving to the right -> we need to compute the weighted average
	// weights are determined by the characteristic velocities (lambdas?) and the difference in the fluid data between the right and left cells.
	// chatgpt: three characteristic velocities:
	// the fluid velocity minus the sound speed, the fluid velocity, and the fluid velocity plus the sound speed.
	} else { // slowest < 0 < fastest => S_L <= 0 <= S_R
		// * readable format
		// double sr_minus_sl = v_char_fastest - v_char_slowest;

		// for (std::size_t i_field = 0; i_field < num_fields; i_field++) {
		// 	double sr_times_fl = v_char_fastest * phys_flux_left_cell[i_field];
		// 	double sr_times_fr = v_char_slowest * phys_flux_right_cell[i_field];
		// 	double sl_times_sr = v_char_slowest * v_char_fastest;
		// 	double qr_minus_ql = fluid_right_cell.fluid_data[i_field] - fluid_left_cell.fluid_data[i_field];

		// 	num_flux[i_field] = (sr_times_fl - sr_times_fr + sl_times_sr * qr_minus_ql) / sr_minus_sl;
		// }

		double sr_minus_sl = 1.0 / (v_char_fastest - v_char_slowest);
		for (std::size_t i_field = 0; i_field < num_fields; i_field++) {
			num_flux[i_field] = (v_char_fastest * phys_flux_left_cell[i_field] - v_char_slowest * phys_flux_right_cell[i_field] +
			                     v_char_slowest * v_char_fastest * (fluid_right_cell.fluid_data[i_field] - fluid_left_cell.fluid_data[i_field])) *
			                    sr_minus_sl;
		}
	}
}
