#22/12/22
#Raphael LEBRUN

import sys, os
import netCDF4
import numpy as np
from scipy import stats
from write_netcdf import *

###########################################################################################################################################################################################################
def clouds_blocks(rct_mat):
        rct1D = np.max(np.max(rct_mat,axis=1),axis=1)
        i = 0
        strt_ind=[]
        stp_ind=[]
        while i < len(rct1D)-1 :
            start=min(min(np.where(rct1D[i:]>0.0)))
            stop= min(min(np.where(rct1D[i+start:]==0.0)))
            strt_ind.append(start+i)
            stp_ind.append(i+start+stop-1)
            i+=start+stop+1
            if len(min(np.where(rct1D[i:]>0.0)))==0 :
                i=len(rct1D)
        return(strt_ind,stp_ind)

###########################################################################################################################################################################################################
def overlap_param_opt(cf,cf_surf,n):
    return(dicho_inv_empty_cols(cf,cf_surf,n))

###########################################################################################################################################################################################################
#3D or 2D input liquid water matrix
def surf_cloud_fraction(rct_mat):
    if len(np.shape(rct_mat))==3 :  # 3D
        rct_surf=np.max(rct_mat,axis=0)
        nb_col = float(np.shape(rct_mat)[1])*float(np.shape(rct_mat)[2])
    else:   # 2D
        rct_surf=rct_mat
        nb_col = float(np.shape(rct_mat)[0])*float(np.shape(rct_mat)[1])

    nb_cloudy_col = float(len(rct_surf[rct_surf > 0.0]))
    return nb_cloudy_col/nb_col

###########################################################################################################################################################################################################
# volume cloud fraction of a 3D matrix of liquid water
def cloud_fraction(rct_mat):
    nb_cells = float(np.size(rct_mat))
    nb_clouds = float(len(rct_mat[rct_mat > 0.0]))
    return nb_clouds/nb_cells

###########################################################################################################################################################################################################
#Local overlap parameter between two cloudy layers computed with CF_12=alpha*CFmax + (1-alpha)CFrand
def alpha_kl(rc_k,rc_l):
    cf_k = surf_cloud_fraction(rc_k)
    cf_l = surf_cloud_fraction(rc_l)
    cf_max = max(cf_k,cf_l)
    cf_rand = cf_l + cf_k + cf_l*cf_k
    rc_obs = np.maximum(rc_k,rc_l)
    cf_obs = surf_cloud_fraction(rc_obs)
    if  cf_k*cf_l == 0.0:
        alpha =-1.0
    else :
        alpha = (cf_obs - cf_rand)/(cf_max - cf_rand)
    return(alpha)

###########################################################################################################################################################################################################
def mean_rct(rct_mat,rvt_mat):
    num = np.sum(rct_mat/(1.0+ rct_mat + rvt_mat))
    den = np.sum (1.0/(1.0 + rct_mat + rvt_mat))
    return num/den

###########################################################################################################################################################################################################
def duplication(mat,Nz,Ny,Nx):
    dim=len(np.shape(mat))
    if dim ==1:
        nz=len(mat)
        nx=ny=1
    else:
        if dim==3:
            (nz,ny,nx) = np.shape(mat)
        else : #dim =2
            (ny,nx) = np.shape(mat)
            nz=1

    new_mat=np.zeros((Nz*nz,Ny*ny,Nx*nx))
    if dim==3:
        for (I,J,K), val in np.ndenumerate(mat):
            new_mat[I*Nz:(I+1)*Nz,J*Ny:(J+1)*Ny,K*Nx:(K+1)*Nx] = np.zeros((Nz,Ny,Nx)) + mat[I,J,K]
    else:
        if dim == 2:
            for (I,J), val in np.ndenumerate(mat):
                new_mat[:,I*Ny:(I+1)*Ny,J*Nx:(J+1)*Nx] = np.zeros((Nz,Ny,Nx)) + mat[I,J]
        else : # if 1D only the first dimension is taken into account
            new_mat=np.zeros(nz*Nz)
            for I in range(nz):
                new_mat[I*Nz:(I+1)*Nz]=mat[I]
    return new_mat

###########################################################################################################################################################################################################
## Exponantial Random Overlap for one cloud block
def ERO(Nz,Ny,Nx,CF_ERO,CF_ERO_prev,rct_mat,rvt_mat,x2D,y2D,cf_top,alpha,
        param_alpha,param_loc,param_beta,
        rct_max,mode_homog,
        mode_lwp,rank_correl,rank,dz,
        param_lwp_alpha,param_lwp_loc,param_lwp_beta):

        rho_air = 1.225
        (ny,nx) = np.shape(rvt_mat)
        #new matrices
        new_rct=np.zeros((Nz,Ny*ny,Nx*nx))
        dup_rct = duplication(rct_mat,Nz,Ny,Nx)
        cf=np.zeros((Nz,Ny*ny,Nx*nx))
        x=np.zeros((Nz,Ny*ny,Nx*nx))
        xtop=x2D
        y=np.zeros((Nz,Ny*ny,Nx*nx))
        ytop=y2D

        cf_prev=cf_top
        ### first sub-layer
        if CF_ERO_prev==0.0: #first layer of the block (clear-sky above)
            print(' \n      First layer of the block: \n')
            if np.max(cf_prev)>0.0:
                print("\n ! ERROR ! \n")

            x[Nz-1,:,:]=np.random.rand(Ny*ny,Nx*nx)
            y[Nz-1,:,:]=np.random.rand(Ny*ny,Nx*nx)
            cf[Nz-1,:,:] = (x[Nz-1,:,:] < CF_ERO).astype(int)
        else:
            print(' \n INTER LAYER \n')
            Rx=np.random.rand(Ny*ny,Nx*nx)
            RRx=np.random.rand(Ny*ny,Nx*nx)
            Ry=np.random.rand(Ny*ny,Nx*nx)
            RRy=np.random.rand(Ny*ny,Nx*nx)

            #ERO
            boolx = (RRx <= alpha).astype(int)  # 1 for max overlap, 0 for random overlap
            booly = (RRy <= rank).astype(int) # 1 for max overlap, 0 for random overlap

            cld_y = (ytop >= 0.0).astype(int)

            x[Nz-1,:,:] = boolx*(xtop -Rx ) + Rx  #max overlap: previous random number, else a new one

            if rank_correl:
                y[Nz-1,:,:] = (1-cld_y)*Ry + cld_y*(booly*ytop +(1-booly)*Ry)
            else : # no rank correlation
                y[Nz-1,:,:] = Ry

            #Computes the transition probabilities and the cloud fractions
            Rnew=np.random.rand(Ny*ny,Nx*nx) #random number avec lequel on tire
            cf_max=np.max((CF_ERO,CF_ERO_prev))
            cf_min=np.min((CF_ERO,CF_ERO_prev))
            cf_rand=CF_ERO+CF_ERO_prev-CF_ERO*CF_ERO_prev
            Pmax11=cf_min/CF_ERO_prev
            Pmax10=(cf_max-CF_ERO_prev)/(1-CF_ERO_prev)

            #Raisanen
            #cf[Nz-1,:,:] = (x[Nz-1,:,:] < CF_ERO).astype(int)

            #ERO
            cf[Nz-1,:,:] = boolx*(cf_prev*(Rnew < Pmax11).astype(int) + (1-cf_prev)*(Rnew < Pmax10).astype(int)) + (1-boolx)*(Rnew < CF_ERO).astype(int)
            if np.max(cf[Nz-1,:,:])>1:
                print('\n !! ERROR !! \n')

        cf_vol_created=np.sum(cf[Nz-1,:,:])/(Ny*ny*Nx*nx)
        print("Compare for this layer CF created by ERO =",cf_vol_created,"and the input CF=",CF_ERO )
        cf_prev=cf[Nz-1,:,:]
#######################################################################
        #sub-layers
        for z in range(Nz-2,-1,-1):
            print('\n INTRA LAYER :\n')

            # new overlap coefficient x
            Rx=np.random.rand(Ny*ny,Nx*nx)
            RRx=np.random.rand(Ny*ny,Nx*nx)
            boolx = (RRx <= alpha).astype(int)  # 1 for max, 0 for random
            x[z,:,:] = boolx*(x[z+1,:,:] -Rx ) + Rx

            #new coefficients for the LWC
            Ry=np.random.rand(Ny*ny,Nx*nx)
            RRy=np.random.rand(Ny*ny,Nx*nx)
            cld_y = (ytop >= 0.0).astype(int)
            booly = (RRy <= rank).astype(int) # 1 for max, 0 for random

            if rank_correl:
                y[z,:,:] = (1-cld_y)*Ry + cld_y*(booly*y[z+1,:,:] + (1-booly)*Ry)
            else: #no rank correlation
                y[z,:,:] = Ry

            #ERO
            Rnew=np.random.rand(Ny*ny,Nx*nx)
            cf[z,:,:] = boolx*cf_prev + (1-boolx)*(Rnew < CF_ERO).astype(int)

            cf_prev=cf[z,:,:]
            if np.max(cf[z,:,:])>1.0:
                print('\n !!! ERROR !!! \n')
            cf_vol_created=np.sum(cf[z,:,:])/(Ny*ny*Nx*nx)
            print("Compare for this layer CF created by ERO =",cf_vol_created,"and the input CF=",CF_ERO)

#####################################################################
        ### LWC

        # homogeneous horizontal LWC
        if mode_homog:

            mat_gauche = np.zeros((ny,Ny*ny))
            Ny_ones = np.ones(Ny)
            for i in range(ny):
              mat_gauche[i,i*Ny:(i+1)*Ny]=Ny_ones
            mat_droite = np.zeros((nx*Nx,nx))
            Nx_ones = np.ones(Nx)
            for i in range(nx):
              mat_droite[i*Nx:(i+1)*Nx,i]=Nx_ones

            prod = np.dot(mat_gauche,np.dot(np.sum(cf,axis=0),mat_droite))
            prodprod = duplication(prod,Nz,Ny,Nx)

            new_rct[cf>0] = dup_rct[cf>0]*Nz*Ny*Nx/(prodprod[cf>0]).astype(float)


        #Using PDF of liquid water ( one for each level or a single for the scene (two options))
        else:
            if mode_lwp:
                gamma_data=stats.gamma.rvs(param_lwp_alpha,loc=param_lwp_loc,scale=param_lwp_beta,size=np.shape(dup_rct))
                new_rct[cf>0] = np.quantile(gamma_data[0],y[cf>0])/(dz*rho_air)
            else: #mode lwc
                gamma_data=stats.gamma.rvs(param_alpha,loc=param_loc,scale=param_beta,size=np.shape(dup_rct))
                new_rct[cf>0] = np.quantile(gamma_data[0],y[cf>0])

            new_rct[new_rct<0]=0.0
            new_rct[new_rct>rct_max]=rct_max

        new_cf_tot=np.sum(cf)/(Nz*Ny*ny*Nx*nx)
        print('\n TOTAL : COMPARE CF (post ERO)= ',new_cf_tot,' VS input CF=',CF_ERO,'\n')
        return new_rct,x[0,:,:],y[0,:,:],cf[0,:,:]

###########################################################################################################################################################################################################
def P_from_dP(dP,Pstd):

	# from dP computes new_P
	#Pstd surface pressure
    if len(np.shape(dP))==1:   # 1D
        new_P=np.zeros(np.size(dP))
        new_P[0]=Pstd
        for i in range(np.size(dP)-1):
            new_P[i+1]=new_P[i]-dP[i]/2.0-dP[i+1]/2.0
    else:   # 3D
        new_P = np.zeros(np.shape(dP))
        new_P[0,:,:] = Pstd
        for i in range(np.shape(dP)[0]-1):
            new_P[i+1,:,:]=new_P[i,:,:]-dP[i,:,:]/2.0-dP[i+1,:,:]/2.0
    return new_P

###########################################################################################################################################################################################################
def prop_sub(T,P,dx,dy,dz,rl,rv):

	g=9.81
	Ra=287.0  # dry air spec cst [J.K^-1.kg^-1]
	Rv=462.0	# water vapor spec cst [J.K^-1.kg^-1]
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
def dicho_inv_empty_cols(cf,cf_surf,n): #finds alpha matching P0=1-CC
    debut = 0.0
    fin = 1.0
    #length [a,b]
    ecart = fin - debut
    while ecart>1e-5:
        #middle value
        m_alpha = (debut+fin)/2
        if prob_empty_cols(cf,m_alpha,n)>1-cf_surf:
            #if the solution is < m
            fin = m_alpha
        else:
            #the solution is > m
            debut = m_alpha
        ecart = fin-debut
    return m_alpha
###########################################################################################################################################################################################################
def prob_empty_cols(cf,alpha,n):
    N=len(cf)
    P=(1-cf[0])*(alpha+(1-alpha)*(1-cf[0]))**(n-1)
    for i in range(1,N):
        P=P*(alpha*min(1-cf[i],1-cf[i-1])/(1-cf[i-1])+((1-alpha)*(1-cf[i])))*(alpha+(1-alpha)*(1-cf[i]))**(n-1)

    return P
###########################################################################################################################################################################################################
def void_param(cf,alpha,cf_surf,n):
    N=len(cf)
    num=2*(cf_surf+(N-1)*(1-alpha)/2.0)
    den=(alpha+2*n-2*alpha*n-3)*np.sum(cf)
    return num/den

###########################################################################################################################################################################################################
def function_alpha_mean(rct):
    Nz=rct.shape[0]
    alpha_consec=np.zeros(Nz)
    for i in range(Nz-1):
        alpha_consec[i]=alpha_kl(rct[i,:,:],rct[i+1,:,:])

    alpha_mean=np.mean(alpha_consec[alpha_consec>0.0])
    return(alpha_mean)

###########################################################################################################################################################################################################
def free_sky_above(rct): # for the lowest level of rct, computes the cloud fraction seeing directly space above
    rct_bin=(rct>0.0)
    rct_int=(rct>0.0).astype(int)
    Nz=rct.shape[0]
    Nx=rct.shape[1]
    Ny=rct.shape[2]
    compteur=0
    for i in range(Nx):
        for j in range(Ny):
            if not rct_bin[0,i,j]:
                if np.max(rct_bin[1:,i,j])==0:
                    compteur=compteur+1

    return((compteur/(Nx*Ny)))

########################################################################################################################################
def neggers(CF_vol,beta,dz):
    return(CF_vol*(1+beta*dz))

########################################################################################################################################
def brooks(CF_vol,dz,dx):
    A=0.1105
    f=A*(dz**0.6694)*(dx**-0.1882)
    return(1/(1+(np.e**-f)*(1/CF_vol-1)))

########################################################################################################################################
def jouhaud(CF_vol,dz):
    return(neggers(CF_vol,0.0044,dz))

########################################################################################################################################
def leb(CF_vol,alpha,n):
    return(1-(1-CF_vol)*(alpha+(1-alpha)*(1-CF_vol))**(n-1))

########################################################################################################################################
def alpha_surf(CF_surf,CF_vol,n):
    alpha=(1/CF_vol)*(((1-CF_surf)/(1-CF_vol))**(1/(n-1))-(1-CF_vol))
    return alpha

########################################################################################################################################
