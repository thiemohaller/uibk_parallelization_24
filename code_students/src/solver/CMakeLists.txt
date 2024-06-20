add_library(solver
    singlestep.cpp
    Riemann_solvers.cpp
    Riemann_solver_HLL.cpp
    finite_volume_solver.cpp
    reconstruction.cpp
    limiter.cpp
    time_integrator.cpp
)

target_include_directories(solver
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    PUBLIC ${MPI_CXX_INCLUDE_PATH}
)

target_link_libraries(solver
    PUBLIC utilities
    PUBLIC sim_setup
    PUBLIC IO
)

target_compile_options(mpi_stuff PRIVATE ${MPI_CXX_COMPILE_FLAGS})
