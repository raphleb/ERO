#22/12/22
#Raphael LEBRUN

import argparse
import sys, os
import netCDF4
import numpy as np
from write_netcdf import *
from ERO_methods import *

######################### Main ######################################################
if __name__ == "__main__" :

################################## PARSER #############################################
  parser = argparse.ArgumentParser(description = "ERO generation from a single GCM column")

  parser.add_argument("GCMfile", type=str,
  help="GCM Netcdf file path")
  parser.add_argument('Dx', metavar='Dx', type=int,
  help='Dx resolution of the GCM netcdf')
  parser.add_argument('Dy', metavar='Dy', type=int,
  help='Dy resolution of the GCM netcdf')
  parser.add_argument('Dz', metavar='Dz', type=int,
  help='Dz resolution of the GCM netcdf')

  parser.add_argument("outputfile", type=str,
  help="output Netcdf file path")
  parser.add_argument('dx', metavar='dx', type=int,
  help='dx resolution after ERO')
  parser.add_argument('dy', metavar='dy', type=int,
  help='dy resolution after ERO')
  parser.add_argument('dz', metavar='dz', type=int,
  help='dz resolution after ERO')

  parser.add_argument("-p","--pdf",type=int,
  default=1,
  help="use a pdf of liquid water at each level ('-p 1', default), homogeneous version('-p 0'), or single pdf version('-p 2')(tests only)")

  parser.add_argument("-v","--vol",type=int,
  default=1,
  help="use CF_vol for the ERO generation ('-v 1',default), else use CF_surf ('-v 0')")

  parser.add_argument("-l","--lwp",action='store_true', default=False,
  help= "use of LWP gamma distribtion instead of LWC pdf. (intended to be used without subgriding.)")

  parser.add_argument("-r","--rank",type=float, default=-1,
  help= "Liquid Water Content rank correlation ('-r 0.9' for example) (unused by default)")

  parser.add_argument("-ra","--rank_overlap",action="store_true", default=False,
  help= "To have the lwc rank correlation parameter equal to the overlap parameter (False by default)")

  args = parser.parse_args()

  input_name=args.GCMfile
  Dx=args.Dx
  Dy=args.Dy
  Dz=args.Dz
  output_name=args.outputfile
  dx=args.dx
  dy=args.dy
  dz=args.dz


  #Default parameters
  mode_vol=True
  mode_single_pdf=False
  mode_homog=False
  mode_pdf=True
  mode_lwp=args.lwp
  rank_a=args.rank_overlap
  rank=args.rank
  rank_correl=(rank>=0.0)

  print('\n \n LWC rank correlation used : r=',rank,'\n \n')

  print("\n\n Running '{}'".format(__file__))
  print("netcdf used " + args.GCMfile + " with resolution " + str((Dx,Dy,Dz))+
  " to generate with ERO " + args.outputfile + " with resolution " + str((dx,dy,dz)))
  if args.pdf==1:
      print(" Using a liquid water pdf at each vertical level\n "+str(args.pdf))
  else :
      if args.pdf==0:
          print(" Homogeneous liquid water distribution at each vertical level\n ")
          mode_homog=True
          mode_pdf=False
      else: #args.pdf==2
          print("Single liquid water pdf\n ")
          mode_single_pdf=True
          mode_pdf=False

  if args.vol==1:
      print('vol==',args.vol)
      print("ERO used with CF_vol (default)\n")

  else:
      mode_vol=False
      print('vol==',args.vol)
      print("ERO used with CF_surf\n")
######################################################################################
  print('Reading the input file...\n')
  GCM=netCDF4.Dataset(input_name)
  rct=GCM.variables['RCT'][0,:,:,:].data
  rvt=GCM.variables['RVT'][0,:,:,:].data
  T=GCM.variables['THT'][0,:,:,:].data
  P=GCM.variables['PABST'][0,:,:,:].data

###########################################################################################################################################################################################################
  #Computing data that might be missing from the input file
  LES_file = 'dP' not in GCM.variables.keys()
  #computing dP and CF (not included in the LES) :
  if LES_file:
      print('Computing dP et and the gamma-pdf of LWC if the input is directly LES files...\n')
      dP=np.zeros(np.shape(rct))
      param_alpha_vec=np.zeros(np.shape(rct[:,0,0]))
      param_loc_vec=np.zeros(np.shape(rct[:,0,0]))
      param_beta_vec=np.zeros(np.shape(rct[:,0,0]))

      #global gamma pdf for the LWC
      rct_vec=rct.flatten()
      param_alpha_tot,param_beta_tot,param_loc_tot= stats.gamma.fit(rct_vec[rct_vec>0.0])
      rct_max_tot=np.max(rct)
      rct_max_vec=np.zeros(np.shape(rct[:,0,0]))

      dP[0,:,:] = P[0,:,:]-P[1,:,:]
      for i in range(np.shape(rct)[0]-1):
          dP[i+1,:,:] = 2.0*(P[i,:,:]-P[i+1,:,:])-dP[i]
      print('Computing the cloud fractions for LES files...\n')
      CF = np.zeros(np.shape(rct)[0]) #colonne !
      for i in range(np.shape(rct)[0]):
          CF[i]=cloud_fraction(rct[i,:,:])
          rct_vec=rct[i,:,:].flatten()
          if np.max(rct_vec)>0:
              param_alpha_vec[i],param_loc_vec[i],param_beta_vec[i]=stats.gamma.fit(rct_vec[rct_vec>0.0])
              rct_max_vec[i]=np.max(rct[i,:i,:])

  else : #file 'GCM_like'
      print('GCM-like FILE \n \n \n')
      CF=GCM.variables['CF'][0,:,:,:].data
      dP=GCM.variables['dP'][0,:,:,:].data
      param_alpha_vec=GCM.variables['gamma_alpha'][0,:].data
      param_loc_vec=GCM.variables['gamma_loc'][0,:].data
      param_beta_vec=GCM.variables['gamma_beta'][0,:].data

      param_lwp_alpha_vec=GCM.variables['lwp_gamma_alpha'][0,:].data
      param_lwp_loc_vec=GCM.variables['lwp_gamma_loc'][0,:].data
      param_lwp_beta_vec=GCM.variables['lwp_gamma_beta'][0,:].data

      param_alpha_tot=GCM.variables['gamma_alpha_total'][:].data
      param_beta_tot=GCM.variables['gamma_beta_total'][:].data
      param_loc_tot=GCM.variables['gamma_beta_total'][:].data

      rct_max_tot=GCM.variables['rct_max_tot'][:].data
      rct_max_vec=GCM.variables['rct_max_loc'][0,:].data

  file_CF_surf = 'CF_surf' in GCM.variables.keys()
  if file_CF_surf:
      CF_surf = GCM.variables['CF_surf'][0,:].data
      CF_surf_total = GCM.variables['CF_surf_total'][0].data
      print("CF_surf already computed\n")
  else:
      CF_surf=np.zeros(np.shape(rct[:,0,0]))
      for i in range(np.shape(rct)[0]):
          CF_surf[i]=surf_cloud_fraction(rct[i,:,:])
      CF_surf_total=surf_cloud_fraction(rct)

  if 'alpha_mean_LES' in GCM.variables.keys():
      alpha_mean_LES_via_GCM=GCM.variables['alpha_mean_LES'][0].data
  else:
      alpha_mean_LES_via_GCM=-1   # indicates the value has not been computed
  GCM.close()

 ###########################################################################################################################################################################################################
  # Multiplying factors (for the numbers of grid cells of the output scene)
  Nx=int(Dx/dx)
  Ny=int(Dy/dy)
  Nz=int(Dz/dz)
  # Shape (in m and number of cells) of the initial input domain
  (nz,ny,nx) = np.shape(rct)
  Lx = nx*Dx
  Ly = ny*Dy
  Lz = nz*Dz
  # Coordinates of the grid cells' centers (KM)
  x_vec=np.linspace(dx/2.0,Lx-dx/2.0,Nx*nx)/1000
  y_vec=np.linspace(dy/2.0,Ly-dy/2.0,Ny*ny)/1000
  z_vec=np.linspace(dz/2.0,Lz-dz/2.0,Nz*nz)/1000
  # changing the resolution of the temperature, liquide water and vapor matrices
  new_T=duplication(T,Nz,Ny,Nx)
  new_dP=np.zeros((Nz*nz,Ny*ny,Nx*nx))
  new_rvt=duplication(rvt,Nz,Ny,Nx)
  new_rct = np.zeros((Nz*nz,Ny*ny,Nx*nx))

#################################################################################################################################
  # Computing the indices of the different cloud blocks of the input column
  (strt_ind,stp_ind) = clouds_blocks(rct)
  print('indices of the base of each cloud block:',strt_ind)
  #print('CF start =',CF[strt_ind[0]])
  print('indices of the top of each cloud block',stp_ind)
  #print('CF sotp=',CF[stp_ind[0]])
  print('%i cloud blocks were identified'%(len(strt_ind)))

  new_dP[:strt_ind[0]*Nz,:,:] = duplication(dP[:strt_ind[0],:,:],Nz,Ny,Nx)/Nz #masse répartie de façon homogène

#############################################################################
  # We compute the overlap parameter for the given cloud cover
  if mode_vol:
      print("VOLUMIC MODE\n\n")
      alpha=overlap_param_opt(CF,CF_surf_total,Nz)
  else:
      print("SURFACIC MODE\n\n")
      alpha=overlap_param_opt(CF_surf,CF_surf_total,Nz)
  if mode_lwp:
      print("LWP being used \n\n")
  else:
      print("LWC being used \n\n")

  #When using Maximum Overlap
  #alpha=1.0

  ### When using the subgrid reconstruction as well as interlayer overlap ###
  #Computes the subgrid overlap parameter
  #alpha_sg=

  #Computes the subgrid surface cloud fraction
  #CF_surf=

  #Computes the interlayer overlap for the given cloud cover
  #alpha_inter=overlap_param_opt(CF_surf,CF_surf_total,1)
  #alpha=alpha_inter

  print('\n \n \n CF_surf_total=',CF_surf_total,'\n \n \n')
  print('overlap parameter =',alpha)
  prob=prob_empty_cols(CF,alpha,Nz)
  print(' \n Probability to generate a clear-sky column=',prob,'prob + CF_surf_total=',prob+CF_surf_total)

  if rank_a: # same rank correlation parameter than the overlap parameter
      rank=alpha

#################################################################################################
# PRINCIPAL ERO LOOP (the new scene is generated from top to bottom, the cloud blocks are treated one after the other )
  for i in range(len(strt_ind)-1):   # cloud bloc i: lvl in [strt_ind[i],stp_ind[i]]
      CF_ERO_prev=0.0 #the layer just above the current cloud block is always clear-sky
      cf_top=np.zeros((Ny*ny,Nx*nx))
      print('Bloc nuageux %i/%i, [%i,%i]'%(i+1,len(strt_ind),strt_ind[i],stp_ind[i]))
      x2D=np.random.rand(Ny*ny,Nx*nx)  # layer above
      y2D=np.random.rand(Ny*ny,Nx*nx)  #layer above
      for lvl in range(stp_ind[i],strt_ind[i]-1,-1): #looping over the layers of the current cloud block
          print('\n \n ***********LVL= *****************\n \n',lvl)
          print('Layer %i/%i,'%(lvl,stp_ind[i]))
          cf_lvl = CF[lvl]
          print('CF_vol = %f'%cf_lvl)
          cf_surf_lvl=CF_surf[lvl]
          print('CF_surf = %f'%cf_surf_lvl)
          if mode_vol:
              CF_ERO=cf_lvl
          else:
              CF_ERO=cf_surf_lvl

          if mode_single_pdf:  # LWC total pdf
              param_alpha_ERO = param_alpha_tot
              param_beta_ERO = param_beta_tot
              param_loc_ERO = param_loc_tot
              rct_max_ERO = rct_max_tot
          else: # LWC pdf local every level
              param_alpha_ERO = param_alpha_vec[lvl]
              param_beta_ERO = param_beta_vec[lvl]
              param_loc_ERO = param_loc_vec[lvl]
              rct_max_ERO = rct_max_vec[lvl]

          param_lwp_alpha_ERO = max(0.0,param_lwp_alpha_vec[lvl])
          param_lwp_beta_ERO = max(0.0,param_lwp_beta_vec[lvl])
          param_lwp_loc_ERO = param_lwp_loc_vec[lvl]

          #version mode_single _pdf option
          new_rct[lvl*Nz:(lvl+1)*Nz,:,:],x_base,y_base,cf_base = ERO(Nz,Ny,Nx,CF_ERO,CF_ERO_prev,
                                                                rct[lvl,:,:],rvt[lvl,:,:],x2D,y2D,
                                                                cf_top,alpha,param_alpha_ERO,param_loc_ERO,
                                                                param_beta_ERO,rct_max_ERO,mode_homog,
                                                                mode_lwp,rank_correl,rank,dz,param_lwp_alpha_ERO,
                                                                param_lwp_loc_ERO,param_lwp_beta_ERO)

          CF_ERO_prev=CF_ERO  # memory of the previous cloud fraction (above)
          cf_top=cf_base

          new_dP[lvl*Nz:(lvl+1)*Nz,:,:] = duplication(dP[lvl,:,:],Nz,Ny,Nx)/Nz
          x2D=x_base
          y2D=y_base

      #between too cloud blocks everything is clear-sky
      new_dP[(stp_ind[i]+1)*Nz:strt_ind[i+1]*Nz,:,:] = duplication(dP[stp_ind[i]+1:strt_ind[i+1],:,:],Nz,Ny,Nx)/Nz

  #Loop for the lowest cloud block (if the scene has one cloud block only this loop will be used)
  x2D=np.random.rand(Ny*ny,Nx*nx)
  y2D=np.random.rand(Ny*ny,Nx*nx)
  CF_ERO_prev=0.0  #the layer above is clear-sky
  cf_top=np.zeros((Ny*ny,Nx*nx))
  for lvl in range(stp_ind[len(strt_ind)-1],strt_ind[len(strt_ind)-1]-1,-1):
      print('\n \n ***** LVL= ******',lvl)
      cf_lvl = CF[lvl]
      print('CF_vol = %f'%cf_lvl)
      cf_surf_lvl=CF_surf[lvl]
      print('CF_surf = %f'%cf_surf_lvl)
      if mode_vol:
          CF_ERO=cf_lvl
      else:
          CF_ERO=cf_surf_lvl

      print('Cloud block (%i/%i): [%i,%i], layer %i/%i' %(len(strt_ind),len(strt_ind),strt_ind[len(strt_ind)-1],stp_ind[len(strt_ind)-1],lvl - strt_ind[len(strt_ind)-1]+1,stp_ind[len(strt_ind)-1]-strt_ind[len(strt_ind)-1]+1))

      if mode_single_pdf:
          param_alpha_ERO = param_alpha_tot
          param_beta_ERO = param_beta_tot
          param_loc_ERO = param_loc_tot
          rct_max_ERO= rct_max_tot
      else: # LWC pdf local every level
          param_alpha_ERO = param_alpha_vec[lvl]
          param_beta_ERO = param_beta_vec[lvl]
          param_loc_ERO = param_loc_vec[lvl]
          rct_max_ERO = rct_max_vec[lvl]

      param_lwp_alpha_ERO = max(0.0,param_lwp_alpha_vec[lvl])
      param_lwp_beta_ERO = max(0.0,param_lwp_beta_vec[lvl])
      param_lwp_loc_ERO = param_lwp_loc_vec[lvl]
      #version mode_single _pdf option
      new_rct[lvl*Nz:(lvl+1)*Nz,:,:],x_base,y_base,cf_base = ERO(Nz,Ny,Nx,CF_ERO,CF_ERO_prev,
                                                            rct[lvl,:,:],rvt[lvl,:,:],x2D,y2D,cf_top,
                                                            alpha,param_alpha_ERO,param_loc_ERO,
                                                            param_beta_ERO,rct_max_ERO,mode_homog,
                                                            mode_lwp,rank_correl,rank,dz,param_lwp_alpha_ERO,
                                                            param_lwp_loc_ERO,param_lwp_beta_ERO)

      CF_ERO_prev=CF_ERO # memory of the previous cloud fraction (above)
      cf_top=cf_base

      print('Computing the new pressure levels')
      new_dP[lvl*Nz:(lvl+1)*Nz,:,:] = duplication(dP[lvl,:,:],Nz,Ny,Nx)/Nz

      x2D=x_base
      y2D=y_base

  #dealing with the lowest clear-sky block (if it exists)
  if stp_ind[len(strt_ind)-1] < nz-1:
      print('Lowest clear-sky block:',stp_ind[len(stp_ind)-1])
      new_dP[(stp_ind[len(strt_ind)-1]+1)*Nz:,:,:] = duplication(dP[stp_ind[len(strt_ind)-1]+1:,:,:],Nz,Ny,Nx)/Nz #masse répartie de façon homogène (air et eau vap)

###########################################################################################################################################################################################################
  # The new pressure levels are compunted with the matrix dP with the function P_from_dP
  P_mat = P_from_dP(new_dP,duplication(P[0,:,:],1,Ny,Nx))
  alpha_mean_LES_ERO=function_alpha_mean(new_rct)

  # Writting the NETCDF file
  write_netcdf(output_name,new_T,P_mat,new_rct,new_rvt,zlvls=z_vec,S_N=y_vec,W_E=x_vec,alpha_opti=alpha,alpha_mean_LES=alpha_mean_LES_via_GCM,alpha_mean_LES_ERO=alpha_mean_LES_ERO)

  print('\n \n LWC rank correlation used : r=',rank,'\n \n')

  print("\n\n Running '{}'".format(__file__))
  print("netcdf used " + args.GCMfile + " with resolution " + str((Dx,Dy,Dz))+
  " to generate with ERO " + args.outputfile + " with resolution " + str((dx,dy,dz)))
  if args.pdf==1:
      print(" \n Using a liquid water pdf at each vertical level "+str(args.pdf))
  else :
      if args.pdf==0:
          print("\n  Homogeneous liquid water distribution at each vertical level ")
          mode_homog=True
          mode_pdf=False
      else: #args.pdf==2
          print("\n Single liquid water pdf ")
          mode_single_pdf=True
          mode_pdf=False

  if args.vol==1:
      print('\n vol==',args.vol)
      print("\n ERO en mode volumique (defaut)")

  else:
      mode_vol=False
      print('\n vol==',args.vol)
      print("\n ERO en mode surfacique (pour tests)")
  print('\n ALPHA OPTI =',alpha)
  print('\n rank = ',rank)
  print('\n rank_correl:',rank_correl)
  print('\n rank_overlap:',rank_a,'(when rank=overlap)')

  print('\n \n \n CF_surf_total=',CF_surf_total,'\n \n \n')
  print('Overlap parameter used =',alpha)
  prob=prob_empty_cols(CF,alpha,Nz)
  print(' \n Probability to generate a clear-sky subcolumn=',prob,'prob + CF_surf_total=',prob+CF_surf_total)
  print('\n Computed cloud cover :',1-prob)

  print(" \n  Done!")
