#Raphaël LEBRUN
# 22/12/22

The two scripts here can be used to sample cloud fraction profiles from cloudy atmospheric files with the NETCDF format. Both ERO.py and degrad.py include an help that can be accessed with the -h option :

	python ERO.py -h
	python degrad.py -h

********************************************
*** Exponential-Random Overlap Sampling ***

ERO.py

This script takes as input a netcdf file consisting of an atmospheric cloudy scene containing liquid water. It produces an output netcdf file consisting of a sample a cloudy subcolumns generated with the wanter spatial resolution and using Exp-Rand Overlap. It can be used this way :

	(base) home:~$ python ERO.py -h
	usage: ERO.py [-h] [-p PDF] [-v VOL] [-l] [-r RANK] [-ra]
		      GCMfile Dx Dy Dz outputfile dx dy dz

	ERO generation from a single GCM column

	positional arguments:
	  GCMfile               GCM Netcdf file path
	  Dx                    Dx resolution of the GCM netcdf
	  Dy                    Dy resolution of the GCM netcdf
	  Dz                    Dz resolution of the GCM netcdf
	  outputfile            output Netcdf file path
	  dx                    dx resolution after ERO
	  dy                    dy resolution after ERO
	  dz                    dz resolution after ERO

	optional arguments:
	  -h, --help            show this help message and exit
	  -p PDF, --pdf PDF     use a pdf of liquid water at each level ('-p 1',
				default), homogeneous version('-p 0'), or single pdf
				version('-p 2')(tests only)
	  -v VOL, --vol VOL     use CF_vol for the ERO generation ('-v 1',default),
				else use CF_surf ('-v 0')
	  -l, --lwp             use of LWP gamma distribtion instead of LWC pdf.
				(intended to be used without subgriding.)
	  -r RANK, --rank RANK  Liquid Water Content rank correlation ('-r 0.9' for
				example) (unused by default)
	  -ra, --rank_overlap   To have the lwc rank correlation parameter equal to
				the overlap parameter (False by default) 

Two single-column files are given (one for the ARM-SGP cloud case, the other for the BOMEX cloud-case) to test the ERO.py script. Their vertical resolution is 100m and they were obtained from high resolution simulations (25m*25m*25m) of 6.4km*6.4km*4km cloud scenes.
Similar single-column files can be produced by the coarsening of the 3D high-resolution files with the degrad.py script.

********************************************
*** Coarsening a high-resolution LES file to a 'GCM-like' single column ***


degrad.py

This script, given a high resolution netcdf fie, transforms it into a single column with the wanted vertical resolution and the same global properties.


	(base) home:~$ python degrad.py -h
	Help to use degradation.py :
	Example :
		python degradation.py input.nc 25 25 25 100 (to coarsen the input.nc file with a  25m*25m*25m resolution,
		into a single column of vertical resolution 100m, with the file name input_small_GCM.nc 
	Usage: python degrad.py input.nc dx dy dz Dz,   output name will be input_small_GCM.nc



