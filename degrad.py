#22/12/22
#Raphael LEBRUN

import os,sys
import plotly.graph_objects as go
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4
from scipy import stats
from degrad_methods import *
from ERO_methods import *
from write_netcdf import *


if len(sys.argv)==2 and sys.argv[1]=="-h":
	print("""Help to use degradation.py :
Example :
	python degradation.py input.nc 25 25 25 100 (to coarsen the input.nc file with a  25m*25m*25m resolution,
	into a single column of vertical resolution 100m, with the file name input_small_GCM.nc """)

if (len(sys.argv)<6) | (len(sys.argv)>6) :
	print("Usage: python %s input.nc dx dy dz Dz,   output name will be input_small_GCM.nc "%sys.argv[0]); exit(1)

LES_filename=sys.argv[1]
LES=netCDF4.Dataset(LES_filename)

rct=LES.variables['RCT'][0,:,:,:].data
rvt=LES.variables['RVT'][0,:,:,:].data
T_pot=LES.variables['THT'][0,:,:,:].data
P=LES.variables['PABST'][0,:,:,:].data

CF_surf_total = surf_cloud_fraction(rct)
print('surface cloud fraction is :',CF_surf_total)

alpha_mean_LES = function_alpha_mean(rct)
print('alpha mean of the LES  = ',alpha_mean_LES)

Pstd=P.max()
print('P(z=0)_min=',np.min(P[0,:,:]))
print('Pstd=',Pstd)
print('Pmax=',P.max())
T=T_pot*((Pstd/P)**(-0.2857))
LES.close()

rho_air = 1.225

#resolution of the input file
dx=int(sys.argv[2])
dy=int(sys.argv[3])
dz=int(sys.argv[4])
#wanted coarsened vertical resolution
Dz=int(sys.argv[5])

rct_vec=rct[:,:,:].flatten()
out_alpha_gama_total,out_loc_gamma_total,out_beta_gamma_total = stats.gamma.fit(rct_vec[rct_vec>0.0])
print('total params gamma : alpha,loc,beta : '
,out_alpha_gama_total,out_loc_gamma_total,out_beta_gamma_total)

out_rct_max_tot=np.max(rct_vec)

Nz=np.size(rct[:,0,0])
f=int(Dz//dz)
A=Nz%f   # si A>0 on distribue sur ceux du haut
B=Nz//f-A    # B blocs de taille f*dz en bas, A blocs de taille (f+1)*dz en haut

out_T=np.zeros(B+A)
out_dP=np.zeros(B+A)
out_rl=np.zeros(B+A)
out_rv=np.zeros(B+A)
out_CF=np.zeros(B+A)
out_CF_surf=np.zeros(B+A)
out_alpha_gama=np.zeros(B+A)
out_loc_gamma=np.zeros(B+A)
out_beta_gamma=np.zeros(B+A)
out_rct_max_loc=np.zeros(B+A)
out_lwp_alpha_gamma=np.zeros(A+B)
out_lwp_beta_gamma=np.zeros(A+B)
out_lwp_loc_gamma=np.zeros(A+B)


for i in range(B):
	print('processing ',i,'/',B,' remains ',A)
	print('fusion du bloc de coordonnees [',i*f , ' , ',(i+1)*f,']')
	[out_T[i],out_dP[i],out_rl[i],out_rv[i],M_a,M_l,M_v]=fusion(T[i*f:(i+1)*f,:,:],
														P[i*f:(i+1)*f,:,:],
														rct[i*f:(i+1)*f,:,:],
														rvt[i*f:(i+1)*f,:,:],dx,dy,dz)
	out_CF[i]=cloud_fraction(rct[i*f:(i+1)*f,:,:])
	out_CF_surf[i]=surf_cloud_fraction(rct[i*f:(i+1)*f,:,:])
	rct_vec=rct[i*f:(i+1)*f,:,:].flatten()
	if np.max(rct_vec)>0.0:
		out_alpha_gama[i],out_loc_gamma[i],out_beta_gamma[i] = stats.gamma.fit(rct_vec[rct_vec>0.0])
		print('param gamma : ',out_alpha_gama[i],out_loc_gamma[i],out_beta_gamma[i])
		out_rct_max_loc[i]=np.max(rct_vec)
		lwp_mat=np.sum(rct[i*f:(i+1)*f,:,:],axis=0)*dz*rho_air
		lwp_vec=lwp_mat.flatten()
		out_lwp_alpha_gamma[i],out_lwp_loc_gamma[i],out_lwp_beta_gamma[i]=stats.gamma.fit(lwp_vec[lwp_vec>0.0])


for i in range(A):
	print('processing ',i,'/',A)
	[out_T[B+i],out_dP[B+i],out_rl[B+i],out_rv[B+i],M_a,M_l,M_v]=fusion(T[B+i*(f+1):B+(i+1)*(f+1),:,:],P[B+i*(f+1):B+(i+1)*(f+1),:,:],rct[B+i*(f+1):B+(i+1)*(f+1),:,:],rvt[B+i*(f+1):B+(i+1)*(f+1),:,:],dx,dy,dz)
	out_CF[i]=cloud_fraction(rct[B+i*(f+1):B+(i+1)*(f+1),:,:])
	out_CF_surf[i]=surf_cloud_fraction(rct[B+i*(f+1):B+(i+1)*(f+1),:,:])
	rct_vec=rct[B+i*(f+1):B+(i+1)*(f+1),:,:].flatten()
	if np.max(rct_vec)>0:
		out_alpha_gama[i],out_loc_gamma[i],out_beta_gamma[i] = stats.gamma.fit(rct_vec[rct_vec>0.0])
		print('param gamma : ',out_alpha_gama[i],out_loc_gamma[i],out_beta_gamma[i])
		out_rct_max_loc[i]=np.max(rct_vec)
		lwp_mat=np.sum(rct[i*f:(i+1)*f,:,:],axis=0)*dz*rho_air
		lwp_vec=lwp_mat.flatten()
		out_lwp_alpha_gamma[i],out_lwp_loc_gamma[i],out_lwp_beta_gamma[i]=stats.gamma.fit(lwp_vec[lwp_vec>0.0])

new_P=P_from_dP(out_dP,Pstd)
new_T_pot=out_T*((Pstd/new_P)**0.2857)
x=y=[50*64/1000]
Lz = (A+B)*Dz
z_vec=np.linspace(Dz/2.0,Lz-Dz/2.0,B+A)/1000



write_netcdf(LES_filename[:-3]+'_small_GCM_'+str(Dz)+'.nc',new_T_pot,new_P,out_rl,out_rv,
			dP=out_dP,CF=out_CF,CF_surf=out_CF_surf,CF_surf_tot = CF_surf_total,
			zlvls=z_vec,S_N=y,W_E=x,
			gamma_alpha=out_alpha_gama,gamma_loc=out_loc_gamma,gamma_beta=out_beta_gamma,
			alpha_mean_LES=alpha_mean_LES,
			gamma_alpha_total=out_alpha_gama_total,gamma_loc_total=out_loc_gamma_total,gamma_beta_total=out_beta_gamma_total,
			rct_max_tot=out_rct_max_tot,rct_max_loc=out_rct_max_loc,
			lwp_gamma_alpha=out_lwp_alpha_gamma,lwp_gamma_loc=out_lwp_loc_gamma,lwp_gamma_beta=out_lwp_beta_gamma)
