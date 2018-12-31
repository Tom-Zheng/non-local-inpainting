cmake_minimum_required(VERSION 2.6)

# DEFAULT BUILD: RELEASE
if (NOT CMAKE_BUILD_TYPE)
   message(STATUS "No build type selected, default to Release")
   set(CMAKE_BUILD_TYPE "Release")
endif()



# IIO (soon to be used)
if( NOT TARGET IIOLIB )
   INCLUDE_DIRECTORIES("3rdparty/iio")
   ADD_SUBDIRECTORY(3rdparty/iio)
ENDIF() 
SET( LIBS ${LIBS} IIOLIB )



## FFTW3 (not used)
#SET( LIBS "${LIBS}" fftw3 m )



# OpenMP (not used yet)
FIND_PACKAGE( OpenMP )
IF(OPENMP_FOUND)
   MESSAGE("OPENMP FOUND")
   SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}" )
ENDIF()

# IF THE COMPILER IS NOT GCC FAIL!
if(CMAKE_COMPILER_IS_GNUCC)
   SET( EXTRA_COMPILER_FLAGS "${EXTRA_COMPILER_FLAGS} -Wall" )

   # compiler flags for 3rdparty/simpois
   set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
else()
   MESSAGE("This project needs c99")
endif()

# UNCOMMENT TO ENABLE DEBUG OUTPUT
#SET( EXTRA_COMPILER_FLAGS "${EXTRA_COMPILER_FLAGS} -DDBG_OUTPUT" )

# SET THE FLAG FOR THE IPOL DEMO 
if (IPOL_DEMO)
   message(STATUS "Generating the exacutable for the IPOL demo")
   SET( EXTRA_COMPILER_FLAGS "${EXTRA_COMPILER_FLAGS} -DIPOL_DEMO" )
endif()

#
add_executable (Inpainting
                a_image_updating.cpp        
                image_inpainting.cpp        
                l1_norm_patch_distance.cpp          
                l2_norm_patch_distance.cpp          
                l2_combined_patch_distance.cpp
                a_patch_distance.cpp
                point.cpp
                distance_transform.cpp      
                patch_match.cpp             
                gaussian_weights.cpp        
                main.cpp                    
                patch_non_local_means.cpp   
                patch_non_local_poisson.cpp
                sampling.cpp
                io_utility.cpp              
                mask_iterator.cpp           
                patch_non_local_medians.cpp 
                shape.cpp
                gradient.cpp
                mask.cpp
                image.hpp
                image.h
                a_image_updating.h        
                image_inpainting.h        
                l1_norm_patch_distance.h          
                l2_norm_patch_distance.h          
                l2_combined_patch_distance.h
                a_patch_distance.h
                point.h
                distance_transform.h      
                patch_match.h             
                gaussian_weights.h        
                patch_non_local_means.h   
                patch_non_local_poisson.h
                sampling.h
                io_utility.h              
                mask_iterator.h           
                patch_non_local_medians.h 
                shape.h
                gradient.h
                mask.h
                3rdparty/simpois/simpois.c
               )
target_link_libraries(Inpainting ${LIBS})
SET_TARGET_PROPERTIES(Inpainting PROPERTIES
      COMPILE_FLAGS "${EXTRA_COMPILER_FLAGS} ${OpenMP_CXX_FLAGS}")

