FROG : Fast Registration Of image Groups 
========================================

### Info ###
This code is the implementation deriving from those papers:

[1] Rémi Agier, Sébastien Valette, Laurent Fanton, Pierre Croisille and Rémy Prost, Hubless 3D Medical Image Bundle Registration, In Proceedings of the 11th Joint Conference on Computer Vision, Imaging and Computer Graphics Theory and Applications - Volume 3: VISAPP, 265-272, 2016, Rome, Italy.

[2] Rémi Agier, Sébastien Valette, Razmig Kéchichian, Laurent Fanton, Rémy Prost, Hubless keypoint-based 3D deformable groupwise registration, arXiv:1809.03951 https://arxiv.org/abs/1809.03951

Authors:
* Rémi Agier : match code
* Sébastien Valette : registration code

### Licence ###
This code is distributed under the CeCILL-B license (BSD-compatible)
(copyright CNRS, INSA-Lyon, UCBL, INSERM.)

###  Dependencies ###
* VTK www.vtk.org
* CMAKE www.cmake.org

###  Simple compilation HowTo under Linux ###
	git clone https://github.com/valette/FROG.git
	cd FROG
	git submodule init
	git submodule update
	cmake . -DCMAKE_BUILD_TYPE=Release
	make

the executables (ACVD, ACVDQ, AnisotropicRemeshingQ and others should be found under the "bin" subdirectory)

### Usage ###

comments, suggestions : https://github.com/valette/FROG/issues