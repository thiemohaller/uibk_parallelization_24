cmake_minimum_required(VERSION 3.10)


add_executable (run_full_code
    main_full_code.cpp
    )

target_include_directories(run_full_code
    PUBLIC ${CODE_INCLUDE_DIR}
)

target_link_libraries(run_full_code
    PUBLIC utilities
    PUBLIC sim_setup
    PUBLIC solver
)

install(TARGETS run_full_code
    RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin
)

# #################### PARALLEL ####################

find_package(HDF5 REQUIRED)
MESSAGE(STATUS "HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}")
find_package(MPI REQUIRED)
MESSAGE(STATUS "MPI_INCLUDE_PATH: ${MPI_INCLUDE_PATH}")


add_executable (run_full_code_parallel
    main_full_code_parallel.cpp
    )

# target_compile_definitions(run_full_code.parallel
#     PRIVATE PARALLEL_VERSION
#     PUBLIC OMPI_SKIP_MPICXX
# )

# target_compile_options(run_full_code.parallel PRIVATE 
#     ${MPI_CXX_COMPILE_FLAGS}
# )

target_compile_options(
    run_full_code_parallel
    PUBLIC ${MPI_C_COMPILE_FLAGS}
)

target_include_directories(run_full_code_parallel
    PUBLIC ${CODE_INCLUDE_DIR}
    PUBLIC ${MPI_CXX_INCLUDE_PATH}
)

target_link_libraries(run_full_code_parallel
    PUBLIC utilities
    PUBLIC sim_setup
    PUBLIC mpi_stuff
    PUBLIC solver
    PUBLIC ${MPI_CXX_LIBRARIES}
)

install(TARGETS run_full_code_parallel
    RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin
)

# #################### PARALLEL ####################

add_executable (run_singlestep
    main_test_singlesteps.cpp
    )

target_include_directories(run_singlestep
    PUBLIC ${CODE_INCLUDE_DIR}
)

target_link_libraries(run_singlestep
    PUBLIC utilities
    PUBLIC sim_setup
    PUBLIC solver
    PUBLIC IO
)

install(TARGETS run_singlestep
    RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin
)