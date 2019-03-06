FROG : Fast Registration Of image Groups 
========================================

### Info ###
This code is the implementation deriving from those papers:

[1] Rémi Agier, Sébastien Valette, Laurent Fanton, Pierre Croisille and Rémy Prost, Hubless 3D Medical Image Bundle Registration, In Proceedings of the 11th Joint Conference on Computer Vision, Imaging and Computer Graphics Theory and Applications - Volume 3: VISAPP, 265-272, 2016, Rome, Italy. https://hal.archives-ouvertes.fr/hal-01284240

[2] Rémi Agier, Sébastien Valette, Razmig Kéchichian, Laurent Fanton, Rémy Prost, Hubless keypoint-based 3D deformable groupwise registration, arXiv:1809.03951 https://arxiv.org/abs/1809.03951

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

### Usage ###

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

comments, suggestions : https://github.com/valette/FROG/issues
