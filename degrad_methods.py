#22/12/22
#Raphael LEBRUN

import os,sys
import netCDF4 as ncdf
import numpy as np
import scipy as sp
import time, datetime
import warnings

###########################################################################################################################################################################################################
def prop_sub(T,P,dx,dy,dz,rl,rv):

	g=9.81
	Ra=287.0  # [J.K^-1.kg^-1]
	Rv=462.0	# [J.K^-1.kg^-1]
	V = dx*dy*dz   #[m^3]
	S=dx*dy	       #[m^2]
	T_P=T/P

	inv_rho_air= Ra*T_P       #inv[kg*m^3]
	inv_rho_liq=0.001
	inv_rho_vap= Rv*T_P

	m_a=V/(inv_rho_air+rl*inv_rho_liq+rv*inv_rho_vap)
	m_l=rl*m_a
	m_v=rv*m_a
	dP=g*(m_a+m_l+m_v)/S

	return [m_a,m_l,m_v,dP]

###########################################################################################################################################################################################################
def temp_moy(vec_T,vec_ma,vec_ml,vec_mv):

	cp_a=1004.0   # [J.kg^-1.K^-1]
	cp_e=2010.0	# [J.kg^-1.K^-1]

	num = np.sum(vec_T*( vec_ma*cp_a + (vec_ml+vec_mv)*cp_e))
	den = np.sum(vec_ma)*cp_a + np.sum(vec_ml+vec_mv)*cp_e

	return num/den

###########################################################################################################################################################################################################
def fusion(mat_T,mat_P,mat_rl,mat_rv,dx,dy,dz):   #with 3D matrices
	g=9.81

	[mat_ma,mat_ml,mat_mv] = prop_sub(mat_T,mat_P,dx,dy,dz,mat_rl,mat_rv)[0:3]
	M_a =np.sum(mat_ma)
	M_l =np.sum(mat_ml)
	M_v =np.sum(mat_mv)

	new_T = temp_moy(mat_T.flatten(),mat_ma.flatten(),mat_ml.flatten(),mat_mv.flatten())

	Rl = M_l/M_a
	Rv = M_v/M_a
	S=dx*np.shape(mat_T)[2]*dy*np.shape(mat_T)[1]

	dP=g*(M_a+M_l+M_v)/S

	return [new_T,dP,Rl,Rv,M_a,M_l,M_v]

###########################################################################################################################################################################################################
def P_from_dP(dP_vec,Pstd):
	# from dP_vec compute P_vec

	P_vec=np.zeros(np.size(dP_vec))
	P_vec[0]=Pstd
	for i in range(np.size(dP_vec)-1):
		P_vec[i+1]=P_vec[i]-dP_vec[i]/2.0-dP_vec[i+1]/2.0

	return P_vec
###########################################################################################################################################################################################################
def P_from_dP_3D(dP_mat,Pstd):
	# from dP_vec compute P_vec

	P_mat=dP_mat*0.0
	P_mat[0,:,:]=Pstd
	for i in range(np.size(dP_mat[:,0,0])-1):
		P_mat[i+1,:,:]=P_mat[i,:,:]-dP_mat[i,:,:]/2.0-dP_mat[i+1,:,:]/2.0

	return P_mat
