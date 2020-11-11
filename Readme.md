FROG : Fast Registration Of image Groups 
========================================
<p align="center">
  <img src="https://www.creatis.insa-lyon.fr/~valette/public/project/frog/featured_huf3b80675463e3a95253905cb1f0a6681_392849_720x0_resize_lanczos_2.png">
</p>

### Info ###
This code is the implementation deriving from those papers:

[1] Rémi Agier, Sébastien Valette, Laurent Fanton, Pierre Croisille and Rémy Prost, Hubless 3D Medical Image Bundle Registration, In Proceedings of the 11th Joint Conference on Computer Vision, Imaging and Computer Graphics Theory and Applications - Volume 3: VISAPP, 265-272, 2016, Rome, Italy. https://hal.archives-ouvertes.fr/hal-01284240

[2] Rémi Agier, Sébastien Valette, Razmig Kéchichian, Laurent Fanton, Rémy Prost, Hubless keypoint-based 3D deformable groupwise registration, Medical Image Analysis, Volume 59, January 2020


Authors:
* Rémi Agier : match code
* Sébastien Valette : registration and tools code

### Licence ###
This code is distributed under the CeCILL-B license (BSD-compatible)
(copyright CNRS, INSA-Lyon, UCBL, INSERM.)

###  Dependencies ###
* Boost www.boost.org
* CMAKE www.cmake.org
* OpenCV www.opencv.org
* VTK www.vtk.org

###  Simple compilation HowTo under Linux ###
	git clone https://github.com/valette/FROG.git
	cd FROG
	git submodule init
	git submodule update
	cmake . -DCMAKE_BUILD_TYPE=Release
	make

Note for FEDORA linux : in case of crashes, please use a self-compiled version of VTK.

### Usage with global script run.sh ###

Groupwise registration is computed via the run.sh script, in three steps:
* Keypoint extraction from input images
* Keypoint matching
* groupwise registration of keypoints

As parameter, run.sh should be fed with a params.sh script containing all input parameters. The run.sh script should not be modified whereas one can copy, rename and edit the params.sh script at will.

To launch registration, execute

	./run.sh ./params.sh

After registration is computed, input images can be transformed into the common space via the transform.sh script, which should also be fed with the parameters file as well as the desired spacing for the transformed images. As an example:

	./transform.sh ./params.sh 2

will transform input images and compute their average with a spacing equal to 2

To check that the output transforms are diffeomorphic:

	./checkDiffeomorphism.sh ./params.sh


### Usage with custom commandline execution ###

Three steps are needed to compute groupwise registration:
* Keypoint extraction from input images

use the surf3d executable to extract keypoints from each 3D image:

	bin/surf3d inputImage [options]

Documentation for this program is available here: https://github.com/valette/vtkOpenSURF3D

* Keypoint matching

Edit a text file (e.g. myPointFiles.txt) containing the list of keypoint files to match.
use the match executable. Then:

	bin/match myPointFiles.txt [options]

Options:

	* -d distance : set max distance between matches. Default : 0.22
	* -d2 distance2secondRatio : maximum closestDistance / distance2second ratio. Default : 1

This program outputs a binary pairs.bin file describing matches between keypoints. Reading the binary file can easily be done, as shown in the frog code: https://github.com/valette/FROG/blob/master/registration/imageGroup.cxx#L953

* Groupwise registration of keypoints

Finally, groupwise registration can occur, with the frog program:

	bin/frog pairs.bin [options]

Options : 

	 * -da <value>  : set alpha for deformable registration. Default : 0.02
	 * -dlinear 0/1 : display linear parameters during registration. default : 0
	 * -dstats 0/1  : display stats during registration. Default : 0
	 * -di number   : number of iterations for each deformable level. Default : 200
	 * -dl number   : number of deformable levels. Default : 3
	 * -emi number  : max number of iterations for EM weighting. Default : 10000
	 * -fi number   : number of fixed images. Default : 0
	 * -fd path     : fixed images transforms directory.
	 * -g spacing   : initial grid spacing for deformable. Default : 100
	 * -gd 0/1      : guaranteed diffeomorphism. Default : 1
	 * -gm ratio    : maximal displacement ratio to guarantee diffeomorphism. Default : 0.4
	 * -il 0/1      : invert landmarks x and y coordinates. Default : 1
	 * -la <value>  : set alpha for linear registration. Default : 0.5
	 * -li number   : number of iterations for linear registration. Default : 50
	 * -nt <number> : set number of threads. Default : number of cores
	 * -s 0/1       : use scale for linear registration. Default : 1
	 * -se number   : stats epsilon. Default : 1e-06
	 * -si number   : interval update for statistics. Default : 10
	 * -ss number   : stats maximal sample size. Default : 10000
	 * -t threshold : inlier probability threshold. Default : 0.5
	 * -l path      : path containing reference landmarks.

This results in a set of transform files (transformXX.json) which can be used to transform images using the VolumeTransform executable

comments, suggestions : https://github.com/valette/FROG/issues
