#include "core/config.hpp"
#include "setup/fluid.hpp"
#include "setup/grid.hpp"
#include "solver/finite_volume_solver.hpp"
#include "solver/time_integrator.hpp"
#include "util/matrix.hpp"
#include "util/utility_functions.hpp"
#include "setup/mpi_handler.hpp"

#include <cmath>
#include <functional>
#include <iostream>
#include <mpi.h>

double Sedov_volume;

void printDebug(const grid_3D& global_grid, const grid_3D& local_grid, const int rank) {
	std::cout << "---------------------- DEBUG GRID RANK " << rank << "----------------------\n";
	std::cout << "Global Grid Information:\n";
	std::cout << "Number of Cells: " << global_grid.get_num_cells(0) << " " << global_grid.get_num_cells(1) << " " << global_grid.get_num_cells(2) << "\n";
	std::cout << "Lower Bound: " << global_grid.x_grid.get_left(0) << " " << global_grid.y_grid.get_left(0) << " " << global_grid.z_grid.get_left(0) << "\n";
	std::cout << "Upper Bound: " << global_grid.x_grid.get_left(global_grid.get_num_cells(0)) << " " << global_grid.y_grid.get_left(global_grid.get_num_cells(1)) << " " << global_grid.z_grid.get_left(global_grid.get_num_cells(2)) << "\n";
	std::cout << "Cell Size: " << global_grid.x_grid.get_dx() << " " << global_grid.y_grid.get_dx() << " " << global_grid.z_grid.get_dx() << "\n";

	std::cout << "----\n";
	std::cout << "Local Grid Information:\n";
	std::cout << "Number of Cells: " << local_grid.get_num_cells(0) << " " << local_grid.get_num_cells(1) << " " << local_grid.get_num_cells(2) << "\n";
	std::cout << "Lower Bound: " << local_grid.x_grid.get_left(0) << " " << local_grid.y_grid.get_left(0) << " " << local_grid.z_grid.get_left(0) << "\n";
	std::cout << "Upper Bound: " << local_grid.x_grid.get_left(local_grid.get_num_cells(0)) << " " << local_grid.y_grid.get_left(local_grid.get_num_cells(1)) << " " << local_grid.z_grid.get_left(local_grid.get_num_cells(2)) << "\n";
	std::cout << "Cell Size: " << local_grid.x_grid.get_dx() << " " << local_grid.y_grid.get_dx() << " " << local_grid.z_grid.get_dx() << "\n";
	std::cout << "---------------------- end of DEBUG GRID RANK " << rank << "----------------------\n";
}

void init_Sedov(fluid_cell &fluid, double x_position, double y_position, double z_position) {
	double radius = sqrt(sim_util::square(x_position) + sim_util::square(y_position) + sim_util::square(z_position));
	double radius_init = 0.05;
	double volume_init = 4.0 / 3.0 * M_PI * radius_init * radius_init * radius_init;
	volume_init = Sedov_volume;
	double E_init = 1.0;
	double e_dens_init = E_init / volume_init;
	if (radius < 0.1) {
		fluid.fluid_data[fluid.get_index_energy()] = e_dens_init;
		fluid.fluid_data[fluid.get_index_tracer()] = 1.0;
	} else {
		fluid.fluid_data[fluid.get_index_energy()] = 1.e-5 / 0.4;
		fluid.fluid_data[fluid.get_index_tracer()] = 0.0;
	}
	fluid.fluid_data[fluid.get_index_density()] = 1.0;
	fluid.fluid_data[fluid.get_index_v_x()] = 0.0;
	fluid.fluid_data[fluid.get_index_v_y()] = 0.0;
	fluid.fluid_data[fluid.get_index_v_z()] = 0.0;
}

int main(int argc, char **argv) {

	MPI_Init(&argc, &argv);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	if(world_rank == 0) {
		std::cout << "World size: " << world_size << std::endl;
	}
	
	std::vector<double> bound_low(3), bound_up(3);
	bound_low[0] = -0.5;
	bound_low[1] = -0.5;
	bound_low[2] = -0.5;

	bound_up[0] = 0.5;
	bound_up[1] = 0.5;
	bound_up[2] = 0.5;

	std::vector<int> num_cells(3);
	num_cells[0] = 32;
	num_cells[1] = 32;
	num_cells[2] = 32;

	std::vector<int> num_tasks(3);
	num_tasks[0] = 2;
	num_tasks[1] = 2;
	num_tasks[2] = 2;
	mpi_handler mpi(num_tasks);

	grid_3D global_grid(bound_low, bound_up, num_cells, 2);
	grid_3D local_grid = mpi.make_local_grid(global_grid);
	printDebug(global_grid, local_grid, world_rank);

	// Get neighboring processors
	// int left_rank, right_rank;
	// MPI_Cart_shift(mpi.comm3D, 0, 1, &left_rank, &right_rank);
	// std::cout << "Left neighbor: " << left_rank << std::endl;
	// std::cout << "Right neighbor: " << right_rank << std::endl;

	// Get number of Sedov cells
	Sedov_volume = 0.0;
	int num_Sedov_cells_local = 0;
	double volume_cell = local_grid.x_grid.get_dx() * local_grid.y_grid.get_dx() * local_grid.z_grid.get_dx();

	for (int ix = 0; ix < local_grid.get_num_cells(0); ++ix) {
		double x_position = local_grid.x_grid.get_center(ix);
		for (int iy = 0; iy < local_grid.get_num_cells(1); ++iy) {
			double y_position = local_grid.y_grid.get_center(iy);
			for (int iz = 0; iz < local_grid.get_num_cells(2); ++iz) {
				double z_position = local_grid.z_grid.get_center(iz);
				double dist = sqrt(sim_util::square(x_position) + sim_util::square(y_position) + sim_util::square(z_position));
				if (dist < 0.1) {
					Sedov_volume += volume_cell;
					num_Sedov_cells_local++;
				}
			}
		}
	}

	// Gather the number of Sedov cells and the total volume
	int total_Sedov_cells = 0;
	double total_Sedov_volume = 0.0;
	MPI_Allreduce(&num_Sedov_cells_local, &total_Sedov_cells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&Sedov_volume, &total_Sedov_volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	// Print the volume of Sedov region on all ranks
	std::cout << " Volume of Sedov region: " << total_Sedov_volume << " in " << total_Sedov_cells << " cells for rank " << world_rank << "\n";

	MPI_Finalize();
	return 0;
	
	// Now, I will create a HD fluid only on ranks with Sedov cells
	if (num_Sedov_cells_local > 0) {
		fluid hd_fluid(parallelisation::FluidType::adiabatic);
		hd_fluid.setup(local_grid);

		std::function<void(fluid_cell &, double, double, double)> function_init = init_Sedov;

		finite_volume_solver solver(hd_fluid);
		solver.set_init_function(function_init);

		double t_final = 0.1;
		double dt_out = 0.005;

		solver.run(local_grid, hd_fluid, t_final, dt_out);
	}

	MPI_Finalize();
	return 0;
}