#22/12/22
#Raphael LEBRUN

import netCDF4 as ncdf
import numpy as np
import scipy as sp
import time, datetime
import warnings
from netCDF4 import Dataset


def write_netcdf(filename,T_mat,P_mat,rl_mat,rv_mat,*args, **kwargs):
	print('Writing .nc file...')
	#creating the file with general attributes
	rootgrp = Dataset(filename, "w", format="NETCDF4")
	rootgrp.title = "ZER"

	#optionnal arguments
	alpha_mean_LES_scal = kwargs.get('alpha_mean_LES',None)
	alpha_mean_LES_ERO_scal = kwargs.get('alpha_mean_LES_ERO',None)
	alpha_opti_scal = kwargs.get('alpha_opti',None)
	dP_mat = kwargs.get('dP', None)    #option argument to fill in dP variable
	z_lvls_vec=kwargs.get('zlvls',None)
	S_N_vec=kwargs.get('S_N',None)
	W_E_vec=kwargs.get('W_E',None)
	print('arg len(W_E)=%i'%len(W_E_vec))
	CF_mat = kwargs.get('CF',None)
	CF_surf_mat = kwargs.get('CF_surf',None)   #to memorize the total cloud cover of the LES being coarsened
	CF_surf_scal = kwargs.get('CF_surf_tot',None)
	param_gamma_alpha= kwargs.get('gamma_alpha',None)
	param_gamma_loc= kwargs.get('gamma_loc',None)
	param_gamma_beta= kwargs.get('gamma_beta',None)
	param_gamma_alpha_total= kwargs.get('gamma_alpha_total',None)
	param_gamma_loc_total= kwargs.get('gamma_loc_total',None)
	param_gamma_beta_total= kwargs.get('gamma_beta_total',None)
	rct_max_loc_vec = kwargs.get('rct_max_loc',None)
	rct_max_tot_scal = kwargs.get('rct_max_tot',None)

	param_lwp_gamma_alpha= kwargs.get('lwp_gamma_alpha',None)
	param_lwp_gamma_loc= kwargs.get('lwp_gamma_loc',None)
	param_lwp_gamma_beta= kwargs.get('lwp_gamma_beta',None)

	qut_mat=kwargs.get('quantiles',None)
	#get the dimensions right
	if len(np.shape(T_mat))==1:
		Nx=Ny=1
		Nz=np.size(T_mat)
	else :
		(Nz,Ny,Nx) = np.shape(T_mat)

	print("SHAPE (Nz,Ny,Nx)",Nz,Ny,Nx)

	#creating the dimensions
	vertical_levels = rootgrp.createDimension("vertical_levels", Nz)
	time = rootgrp.createDimension("time", 1)    #unlimited mais en fait à 1
	S_N_direction = rootgrp.createDimension("S_N_direction", Ny)
	W_E_direction = rootgrp.createDimension("W_E_direction", Nx)
	print('Nx=%i'%Nx)


	#creating the variables
	time = rootgrp.createVariable("time","i4",("time",))
	time.units = "seconds since 1999-1-1 00:00:00"
	W_E_direction = rootgrp.createVariable("W_E_direction","f4",("W_E_direction",))
	W_E_direction.units = "km"
	S_N_direction = rootgrp.createVariable("S_N_direction","f4",("S_N_direction",))
	S_N_direction.units = "km"
	vertical_levels = rootgrp.createVariable("vertical_levels","f4",("vertical_levels",))
	vertical_levels.units = "km"
	LAT = rootgrp.createVariable("LAT","f4",("S_N_direction",))
	LAT.units = "degrees_north"
	LAT.longname = "latitudes"
	LON = rootgrp.createVariable("LON","f4",("W_E_direction",))
	LON.units = "degrees_east"
	LON.longname = "longitudes"
	THT = rootgrp.createVariable("THT","f4",("time","vertical_levels","S_N_direction","W_E_direction",))
	THT.units = "K"
	PABST = rootgrp.createVariable("PABST","f4",("time","vertical_levels","S_N_direction","W_E_direction",))
	PABST.units = "Pa"
	dP = rootgrp.createVariable("dP","f4",("time","vertical_levels","S_N_direction","W_E_direction",))
	dP.units = "Pa"
	RCT = rootgrp.createVariable("RCT","f4",("time","vertical_levels","S_N_direction","W_E_direction",))
	RCT.units = "kg/kg"
	RCT.longname = "liquid water mixing ratio"
	RVT = rootgrp.createVariable("RVT","f4",("time","vertical_levels","S_N_direction","W_E_direction",))
	RVT.units = "kg/kg"
	RVT.longname = "vapor water mixing ratio"
	CF = rootgrp.createVariable("CF","f4",("time","vertical_levels","S_N_direction","W_E_direction",))
	CF.units = "m^-3"
	CF.longname = "Cloud fraction (per volume)"
	CF_surf = rootgrp.createVariable("CF_surf","f4",("time","vertical_levels",))
	CF_surf.units = "m^-2"
	CF_surf.longname = "Surface cloud fraction (per area)"
	CF_surf_total = rootgrp.createVariable("CF_surf_total","f4",("time",))
	CF_surf_total.units = "m^-2"
	CF_surf_total.longname = "Total surface cloud fraction"

	gamma_alpha=rootgrp.createVariable("gamma_alpha","f4",("time","vertical_levels",))
	gamma_loc=rootgrp.createVariable("gamma_loc","f4",("time","vertical_levels",))
	gamma_beta=rootgrp.createVariable("gamma_beta","f4",("time","vertical_levels",))

	lwp_gamma_alpha=rootgrp.createVariable("lwp_gamma_alpha","f4",("time","vertical_levels",))
	lwp_gamma_loc=rootgrp.createVariable("lwp_gamma_loc","f4",("time","vertical_levels",))
	lwp_gamma_beta=rootgrp.createVariable("lwp_gamma_beta","f4",("time","vertical_levels",))

	gamma_alpha_total=rootgrp.createVariable("gamma_alpha_total","f4",("time",))
	gamma_loc_total=rootgrp.createVariable("gamma_loc_total","f4",("time",))
	gamma_beta_total=rootgrp.createVariable("gamma_beta_total","f4",("time",))

	alpha_mean_LES = rootgrp.createVariable("alpha_mean_LES","f4",("time",))
	alpha_mean_LES_ERO = rootgrp.createVariable("alpha_mean_LES_ERO","f4",("time",))
	alpha_opti = rootgrp.createVariable("alpha_opti","f4",("time",))

	QUT = rootgrp.createVariable("QUT","f4",("time","vertical_levels","S_N_direction","W_E_direction",))
	QUT.longname = "quantile of liquid water content"

	rct_max_tot=rootgrp.createVariable("rct_max_tot","f4",("time",))
	rct_max_loc=rootgrp.createVariable("rct_max_loc","f4",("time","vertical_levels",))

	#writting each data
	THT[0,:,:,:]=T_mat
	PABST[0,:,:,:]=P_mat
	dP[0,:,:,:]=dP_mat
	RCT[0,:,:,:]=rl_mat
	RVT[0,:,:,:]=rv_mat
	vertical_levels[:]=z_lvls_vec
	S_N_direction[:]=S_N_vec
	W_E_direction[:]=W_E_vec
	CF[:]=CF_mat
	CF_surf[0,:]=CF_surf_mat
	CF_surf_total[:] = CF_surf_scal
	gamma_alpha[0,:]=param_gamma_alpha
	gamma_loc[0,:]=param_gamma_loc
	gamma_beta[0,:]=param_gamma_beta
	alpha_mean_LES[:]=alpha_mean_LES_scal
	alpha_mean_LES_ERO[:]=alpha_mean_LES_ERO_scal
	alpha_opti[:]=alpha_opti_scal
	gamma_alpha_total[:]=param_gamma_alpha_total
	gamma_loc_total[:]=param_gamma_loc_total
	gamma_beta_total[:]=param_gamma_beta_total

	lwp_gamma_alpha[0,:]=param_lwp_gamma_alpha
	lwp_gamma_loc[0,:]=param_lwp_gamma_loc
	lwp_gamma_beta[0,:]=param_lwp_gamma_beta

	QUT[0,:,:,:]=qut_mat

	rct_max_loc[0,:] = rct_max_loc_vec
	rct_max_tot[:] = rct_max_tot_scal

	rootgrp.close()
