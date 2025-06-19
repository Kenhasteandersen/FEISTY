add_library(fabm_models_feisty OBJECT
            feisty_model_library.F90
	    input.F90
            FEISTY.F90
            setup.F90
            fish.F90
            spectrum.F90
            globals.F90
            FEISTY_FABM.F90
            setupbasic2.F90
           )

#target_link_libraries(fabm_models_feisty PRIVATE fabm_base)
add_dependencies(fabm_models_feisty fabm_base)

# Enable preprocessing for different compilers and define _FABM_ for the OBJECT library
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(fabm_models_feisty PRIVATE /fpp)
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(fabm_models_feisty PRIVATE -cpp)
endif()

target_compile_definitions(fabm_models_feisty PRIVATE _FABM_)
