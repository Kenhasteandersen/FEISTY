add_library(fabm_models_feisty OBJECT
            feisty_model_library.F90
            FEISTY.F90
            setup.F90
            fish.F90
            spectrum.F90
            globals.F90
            setupbasic.F90
            setupbasic2.F90
           )

#target_link_libraries(fabm_models_feisty PRIVATE fabm_base)
add_dependencies(fabm_models_feisty fabm_base)
