add_library(IO 
    data_storage.cpp
    data_storage_parallel.cpp
)

target_include_directories(IO
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    PUBLIC ${MPI_CXX_INCLUDE_PATH}
)

target_link_libraries(IO
    PUBLIC SerialHdf5
)

target_compile_options(mpi_stuff PRIVATE ${MPI_CXX_COMPILE_FLAGS})
