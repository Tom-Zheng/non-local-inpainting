Variational Framework for Non-Local Inpainting
==============================================

Build instructions
==================
To compile this code cmake, make and a C99 compiler (ie. gcc) are required.
The libraries, libpng, libjpeg and libtiff are optional but WE RECOMMEND THEM.

Compile the binaries by running the following line in the current directory: 

    $ mkdir build; cd build; cmake ..; make

the process will produce the binary: Inpainting



Usage
=====
Running the program without parameters prints its usage instructions:

    $ build/Inpainting

      USAGE:
      ./build/Inpainting input mask output [OPTIONS]
      
      Available options are:
       -depth     Specify depth input
       -method    method name (nlmeans)
       -patch     patch side (9)
       -iters     inpainting iterations (200)
       -scales    scales amount (11)
       -coarse    coarsest rate (0.1)
       -conft     confidence decay time (5)
       -confa     confidence asymptotic value (0.1)
       -lambda    lambda (0.05)
       -lambdar   lambda for RGB channel (only for nldeph method)
       -lambdad   lambda for depth channel (only for nldeph method)
       -init      initialization type [poisson/black/avg/none] (poisson)
       -psigma    Gaussian patch weights (10000)
       -showpyr   PREFIX write intermediate pyramid results
       -shownnf   FILENAME write illustration of the final NNF


## Example call:
    
    $ build/Inpainting data/img.png data/mask.png output.png -depth data/depth.png -patch 9 -method nldepth -lambdar 1.0 -lambdad 0.0
    
    Note: lambda ranges in [0..1], specifying the proportion of intensity information over gradient information.

Files and main functions 
========================

    image.hpp           : FixedImage and Image classes. Image with some defined
                          pixel type (e.g. float, Point, etc.). 
                          The word 'fixed' in the name should be considered as 
                          'immutable', because this class does not provide any 
                          capabilities to modify its data
    mask.cpp            : FixedMask and Mask classes. A binary mask container 
    mask_iterator.cpp   : STL-like bidirectional input (constant) iterator

    point.cpp           : container for 2D coordinates (x, y)
    shape.cpp           : container for a 2D rectangular shape (width, height)

    a_image_updating.cpp         : abstract class and insances of the image 
    patch_non_local_medians.cpp    update step of the algorithm: patch-NLmeans,
    patch_non_local_means.cpp      patch-NLmedians, and patch-NLpoisson
    patch_non_local_poisson.cpp

    a_patch_distance.cpp         : abstract class and instance of the patch
    l1_norm_patch_distance.cpp     distances: l1, l2, and gradient-based l2
    l2_norm_patch_distance.cpp
    l2_combined_patch_distance.cpp

    image_inpainting.cpp         : ImageInpainting algorithm 

    patch_match.cpp              : PatchMatch algorithm

    distance_transform.cpp       : compute the distance function to a set
    gaussian_weights.cpp         : compute gaussian weighted patches
    gradient.cpp                 : compute image gradients
    sampling.cpp                 : up/down-sampling utils for multiscale
    io_utility.cpp               : read and write image files
    main.cpp                     : main program (see next section)


## Main functions

Setup the components of the algorithm and launch the inpainting.
   
   
     // setup patchmatch
     Shape  patch_size =  Shape(patch_side, patch_side);
     APatchDistance *patch_distance = new  L2CombinedPatchDistance(lambda, 
                                            patch_size, patch_sigma);
     PatchMatch *patch_match        = new  PatchMatch( patch_distance, 
                                            patch_match_iterations, 
                                            random_shots_limit, -1);
                                   
     // setup image update                                
     AImageUpdating *image_updating = new PatchNonLocalPoisson( patch_size, 
                                           patch_sigma, lambda, conj_grad_tol, 
                                           conj_grad_iter);
   
     // create inpainting object
     ImageInpainting image_inpainting = ImageInpainting(inpainting_iterations,
                                          tolerance,
                                          scales_amount,
                                          subsampling_rate,
                                          confidence_decay_time,
                                          confidence_asymptotic_value,                                       
                                          initialization_type);     // enum
   
     image_inpainting.set_weights_updating(patch_match);
     image_inpainting.set_image_updating(image_updating);


     // image of type Image, and mask of type Mask
     Image<float> output = image_inpainting.process(input, mask);



Third party components
======================

This software uses the following 3rd party components:

* 3rdparty/iio/ contains Enric Meinhardt's iio library for image input/output: 
    https://github.com/mnhrdt/iio

* 3rdparty/simpois/ contains simpois.c (and other files used by simpois.c) 
  from Enric Meinhardt's image script collection:
    https://github.com/mnhrdt/imscript

  simpois.c provides a simpler Poisson solver which is used to compute the
  initial condition in the coarsest scale.

Both 3rdparty components are distributed under the terms of the BSD lincense,
see 3rdparty/iio/LICENSE and 3rdparty/simpois/LICENSE. 
