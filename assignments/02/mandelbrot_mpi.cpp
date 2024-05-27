#include <mpi.h>
#include <bits/chrono.h>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <sys/time.h>
#include <tuple>
#include <vector>

// Include that allows to print result as an image
// Also, ignore some warnings that pop up when compiling this as C++ mode
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#pragma GCC diagnostic pop

// 8K by default
constexpr int default_size_x = 7680;
constexpr int default_size_y = 4320;

// RGB image will hold 3 color channels
constexpr int num_channels = 3;
// max iterations cutoff
constexpr int max_iterations = 10000;

#define IND(Y, X, SIZE_Y, SIZE_X, CHANNEL) (Y * SIZE_X * num_channels + X * num_channels + CHANNEL)

size_t index(int y, int x, int /*size_y*/, int size_x, int channel) {
	return y * size_x * num_channels + x * num_channels + channel;
}

using Image = std::vector<uint8_t>;

auto HSVToRGB(double H, const double S, double V) {
	if (H >= 1.0) {
		V = 0.0;
		H = 0.0;
	}

	const double step = 1.0 / 6.0;
	const double vh = H / step;

	const int i = (int)floor(vh);

	const double f = vh - i;
	const double p = V * (1.0 - S);
	const double q = V * (1.0 - (S * f));
	const double t = V * (1.0 - (S * (1.0 - f)));
	double R = 0.0;
	double G = 0.0;
	double B = 0.0;

	// clang-format off
	switch (i) {
	case 0: { R = V; G = t; B = p; break; }
	case 1: { R = q; G = V; B = p; break; }
	case 2: { R = p; G = V; B = t; break; }
	case 3: { R = p; G = q; B = V; break; }
	case 4: { R = t; G = p; B = V; break; }
	case 5: { R = V; G = p; B = q; break; }
	}
	// clang-format on

	return std::make_tuple(R, G, B);
}

void calcMandelbrot(Image &image, int size_x, int size_y, int start_y, int num_rows, int rank) {
	std::cout << "Rank " << rank << " is calculating " << num_rows << " rows starting from " << start_y << " until " << start_y + num_rows << std::endl;

	const float left = -2.5, right = 1;
	const float bottom = -1, top = 1;

	// TODOs for MPI parallelization
	// 1) domain decomposition
	//   - decide how to split the image into multiple parts
	// -> use rows for splits
	//   - ensure every rank is computing its own part only
	// -> use subimages based on rank
	// 2) result aggregation
	//   - aggregate the individual parts of the ranks into a single, complete image on the root rank (rank 0)
	// -> use gather to collect all subimages on root

	auto time_start = std::chrono::high_resolution_clock::now();

	for (int pixel_y = start_y; pixel_y < start_y + num_rows; pixel_y++) {
		// scale y pixel into mandelbrot coordinate system
		const float cy = (pixel_y / (float)size_y) * (top - bottom) + bottom;
		for (int pixel_x = 0; pixel_x < size_x; pixel_x++) {
			// scale x pixel into mandelbrot coordinate system
			const float cx = (pixel_x / (float)size_x) * (right - left) + left;
			float x = 0;
			float y = 0;
			int num_iterations = 0;

			// Check if the distance from the origin becomes
			// greater than 2 within the max number of iterations.
			while ((x * x + y * y <= 2 * 2) && (num_iterations < max_iterations)) {
				float x_tmp = x * x - y * y + cx;
				y = 2 * x * y + cy;
				x = x_tmp;
				num_iterations += 1;
			}

			// Normalize iteration and write it to pixel position
			double value = fabs((num_iterations / (float)max_iterations)) * 200;
			auto [red, green, blue] = HSVToRGB(value, 1.0, 1.0);

			int channel = 0;
			image[index(pixel_y - start_y, pixel_x, size_y, size_x, channel++)] = (uint8_t)(red * UINT8_MAX);
			image[index(pixel_y - start_y, pixel_x, size_y, size_x, channel++)] = (uint8_t)(green * UINT8_MAX);
			image[index(pixel_y - start_y, pixel_x, size_y, size_x, channel++)] = (uint8_t)(blue * UINT8_MAX);
		}
	}

	auto time_end = std::chrono::high_resolution_clock::now(); // Stop time for each node
	auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
	std::cout << "Mandelbrot set calculation for Rank " << rank << " took: " << time_elapsed << " ms." << std::endl;
}

int main(int argc, char **argv) {
	// mpi stuff
	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	std::cout << "Ranks " << rank << " of " << size << " started." << std::endl;

	int size_x = default_size_x;
	int size_y = default_size_y;

	if (argc == 3) {
		size_x = atoi(argv[1]);
		size_y = atoi(argv[2]);
		if (rank == 0) {
			std::cout << "Using size " << size_x << "x" << size_y << std::endl;
		}
	} else {
		if (rank == 0) {
			std::cout << "No arguments given, using default size " << size_x << "x" << size_y << std::endl;
		}
	}

	// -----
	int num_rows_per_rank = size_y / size;
	int start_y = rank * num_rows_per_rank; // slice based on rank
	int num_rows = (rank == size - 1) ? (size_y - start_y) : num_rows_per_rank; // handle remainder if y is not divisible by size

	int sub_image_size = num_channels * size_x * num_rows;
	Image sub_image(sub_image_size);
	auto time_start = std::chrono::high_resolution_clock::now();

	// compute mandelbrot locally, based on subimage
	calcMandelbrot(sub_image, size_x, size_y, start_y, num_rows, rank);

	// originally I had this logic in the calcMandelbrot function, however, chatGPT suggested to move it here
	// also used https://www.mpi-forum.org/docs/mpi-1.1/mpi-11-html/node70.html for examples of mpi_gather
	// on root, gather all subimages and write to file
	if (rank == 0) {
		Image final_image(num_channels * size_x * size_y);
		MPI_Gather(sub_image.data(), sub_image_size, MPI_UINT8_T,
				   final_image.data(), sub_image_size, MPI_UINT8_T, 0, MPI_COMM_WORLD);

		// stop time
		auto time_end = std::chrono::high_resolution_clock::now();
		auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
		std::cout << "In total Mandelbrot set calculation for image size" << size_x << "x" << size_y << " took: " << time_elapsed << " ms." << std::endl;

		constexpr int stride_bytes = 0;
		stbi_write_png("mandelbrot_mpi.png", size_x, size_y, num_channels, final_image.data(), stride_bytes);
	} else {
		// if not root, send data to root (use nullptr to receive as seen in the lecture)
		MPI_Gather(sub_image.data(), sub_image_size, MPI_UINT8_T,
				   nullptr, 0, MPI_UINT8_T, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return EXIT_SUCCESS;
}
