#import matplotlib   # import matplotlib for next statement
#matplotlib.use('Agg')   # this agg backend for plotting supports pdf, pngimport numpy as np

import os
import sys
import math
import readsnap as rs
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import LogNorm
#import wrapper2


#from pylab import * #import pylab for ioff()
#ioff()  # set interactive to off so no plotting to x-window


def drawMassFromIMF(x):
    A = 0.126512
    y = 0.0
    if( x < 0.5):
        y = 2. * A * x**(-1.3)
    elif(x >= 0.5):
        y = A * x**(-2.3)

    return y


#filename_path="/mnt/ceph/users/chu/snapshots/cooling_test/wss_cie_cool"
#filename_path="/mnt/ceph/users/chu/snapshots/cooling_test/wss_cie_cool_const_rho"




#filename_path="/mnt/ceph/users/chu/snapshots/dwarf_chem/ng1e7HI4e7gsl2_cutDM_PE_PI_SN_localSh_eps1e20_soft05pc_sIMF2myr4pc_wsscie"
#filename_path="/mnt/ceph/users/chu/snapshots/dwarf_chem/ng1e7HI4e7gsl2_cutDM_PE_PI_SN_localSh_eps1e20_soft05pc_sIMF1myr10pc"
#filename_path="/mnt/ceph/users/chu/snapshots/dwarf_chem/ng1e7HI4e7gsl2_cutDM_PE_PI_SN_localSh_eps1e20_soft05pc_sIMF1myr4pc"
#filename_path="/mnt/ceph/users/chu/snapshots/dwarf_chem/ng1e7HI4e7gsl2_cutDM_PE_PI_SN_localSh_eps1_soft05pc_sIMF2myr4pc_wsscie_adpsoft"

#filename_path="/mnt/ceph/users/chu/snapshots/dwarf_chem/ng1e7HI4e7gsl2_cutDM_PE_PI_SN_localSh_eps1e20"
#filename_path="/mnt/ceph/users/chu/snapshots/dwarf_chem/ng1e7HI4e7gsl2_cutDM_PE_PI_SN_localSh_eps0p02"
#filename_path="/mnt/ceph/users/chu/snapshots/dwarf_chem/ng1e7HI4e7gsl2_cutDM_PE_PI_SN_localSh_eps0p02_soft05pc_sIMF1myr500pc"
#filename_path="/mnt/ceph/users/chu/snapshots/dwarf_chem/ng1e7HI1e7gsl1_PE_PI_SN_localSh_cutDM"
#filename_path="/mnt/ceph/users/chu/snapshots/isoTdisk/ng8e7HI1e7gsl1_isoT_cutDM"
#filename_path="/mnt/ceph/users/chu/snapshots/sampleIMF/sampleIMFtest_noAGB_N64_rng"
filename_path = "/mnt/ceph/users/chu/snapshots/dwarf_chem/ng1e7HI4e7gsl2_cutDM_PE_PI_SN_localSh_eps1_soft05pc_sIMF2myr4pc_wsscie_adpsoft"





filename_base = filename_path + "/snap_"
#filename_base = filename_path + "/snap_MW_mres_nort_"

Hydrogen_massfrac=0.76
XH=Hydrogen_massfrac
yhelium=(1-XH)/(4*XH)
GAMMA=5.0/3.0
GAMMA_MINUS1=GAMMA-1
BOLTZMANN=1.3806e-16
PROTONMASS=1.6726e-24
HUBBLE=0.65
GRAVCON=6.67e-8
UnitMass_in_g = 1.989e43
UnitDensity_in_cgs = 6.76991e-22
UnitDensity_in_pccm = UnitDensity_in_cgs/PROTONMASS 
UnitLength_in_cm    = 3.085678e21
UnitTime_in_s = 3.08568e+16
Year_in_s = 31556926.

#rs.list_format2_blocks(filename)


jump=10 #What is this
N_snap = 1000 / jump

sfr_tot = np.zeros(len(range(0,1000/jump)))

thisSnap = 182

for k in range(0,1000/jump):
#for k in range(thisSnap, thisSnap+1):
    kk = k*jump
    if (kk < 10):
        num = '00' + str(kk)
    elif (kk < 100):
        num = '0' + str(kk)
    elif (kk < 10000):
        num = str(kk)
        

    
    filename = filename_base + num


    plt.clf()

    print ''
    print 'Read snapshot ', num
    head= rs.snapshot_header(filename)

    Ngas   = head.npart[0]
    Ndm    = head.npart[1]
    Ndisk  = head.npart[2]
    Nbulge = head.npart[3]
    Nstar  = head.npart[4]
    


    #-------- Dark matter --------
    if (Ndm > 0):
        pos_dm = rs.read_block(filename, "POS ", parttype=1) 
        x_dm = pos_dm[:,0]
        y_dm = pos_dm[:,1]
        z_dm = pos_dm[:,2]

        vel_dm = rs.read_block(filename, "VEL ", parttype=1) 
        vx_dm = vel_dm[:,0]
        vy_dm = vel_dm[:,1]
        vz_dm = vel_dm[:,2]

        m_dm = rs.read_block(filename, "MASS", parttype=1) 
    
    #-------- Disk --------
    if (Ndisk > 0):
        pos_disk = rs.read_block(filename, "POS ", parttype=2) 
        x_disk = pos_disk[:,0]
        y_disk = pos_disk[:,1]
        z_disk = pos_disk[:,2]
    
        vel_disk = rs.read_block(filename, "VEL ", parttype=2) 
        vx_disk = vel_disk[:,0]
        vy_disk = vel_disk[:,1]
        vz_disk = vel_disk[:,2]

        m_disk   = rs.read_block(filename, "MASS", parttype=2) 

    #-------- Bulge --------
    if (Nbulge > 0):
        pos_bulge = rs.read_block(filename, "POS ", parttype=3) 
        x_bulge = pos_bulge[:,0]
        y_bulge = pos_bulge[:,1]
        z_bulge = pos_bulge[:,2]
        
        vel_bulge = rs.read_block(filename, "VEL ", parttype=3) 
        vx_bulge = vel_bulge[:,0]
        vy_bulge = vel_bulge[:,1]
        vz_bulge = vel_bulge[:,2]

    #-------- Star --------
    if (Nstar > 0):
        pos_star = rs.read_block(filename, "POS ", parttype=4) 
        x_star = pos_star[:,0]
        y_star = pos_star[:,1]
        z_star = pos_star[:,2]
    
        vel_star = rs.read_block(filename, "VEL ", parttype=4) 
        vx_star = vel_star[:,0]
        vy_star = vel_star[:,1]
        vz_star = vel_star[:,2]

        id_star  = rs.read_block(filename, "ID  ", parttype=4) 
        m_star   = rs.read_block(filename, "MASS", parttype=4) 
        Zm_star = rs.read_block(filename, "Z   ", parttype=4, csformat = 1) 
        Zm_star_total = (Zm_star[:,1] + Zm_star[:,2] + Zm_star[:,3] + Zm_star[:,4] + Zm_star[:,5] + Zm_star[:,7] + Zm_star[:,8] + Zm_star[:,9] + Zm_star[:,10] + Zm_star[:,11]) / m_star[:]
	age   = head.time - rs.read_block(filename, "AGE ", parttype=4)
        m_imf = rs.read_block(filename, "MIMF", parttype=4)


    #-------- Gas --------
    if (Ngas > 0):
        pos_gas = rs.read_block(filename, "POS ", parttype=0) 
        x_gas = pos_gas[:,0]
        y_gas = pos_gas[:,1]
        z_gas = pos_gas[:,2]
    
        vel_gas = rs.read_block(filename, "VEL ", parttype=0) 
        vx_gas = vel_gas[:,0]
        vy_gas = vel_gas[:,1]
        vz_gas = vel_gas[:,2]
    
        id_gas  = rs.read_block(filename, "ID  ", parttype=0) 
        m_gas   = rs.read_block(filename, "MASS", parttype=0) 
        
        u       = rs.read_block(filename, "U   ", parttype=0) 
        rho     = rs.read_block(filename, "RHO ", parttype=0) 
#        ne      = rs.read_block(filename, "NE  ", parttype=0) 
#        nh      = rs.read_block(filename, "NH  ", parttype=0) 
        hsml    = rs.read_block(filename, "HSML", parttype=0) 
        sfr_pp  = rs.read_block(filename, "SFR ", parttype=0)   #star formation rate per particle
        Zm_gas  = rs.read_block(filename, "Z   ", parttype=0, csformat = 1) 
        #cs_temp = rs.read_block(filename, "CSTE", parttype=0) 
        chemT   = rs.read_block(filename, "CHET", parttype=0)   #####gas temperature
#        v_disp  = rs.read_block(filename, "CSSI", parttype=0)   
#        G0      = rs.read_block(filename, "G0  ", parttype=0)
	#coln  = rs.read_block(filename, "COLN", parttype=0)   


        print sum(sfr_pp)
        sfr_tot[k]=sum(sfr_pp)


	#col = np.empty(Ngas)
	#col[:] = rho[:] * hsml[:] + coln[:,0]+coln[:,1]+coln[:,2]+coln[:,3]+coln[:,4]+coln[:,5]+coln[:,6]+coln[:,7]+coln[:,8]+coln[:,9]+coln[:,10]+coln[:,11]
	#col_local = rho * hsml
        
	#col       *= UnitDensity_in_pccm * UnitLength_in_cm / (1.0 + 4.0 * 0.1)
        #col_local *= UnitDensity_in_pccm * UnitLength_in_cm / (1.0 + 4.0 * 0.1)


        Zm_gas_total  = (Zm_gas[:,1]  + Zm_gas[:,2]  + Zm_gas[:,3]  + Zm_gas[:,4]  + Zm_gas[:,5]  + Zm_gas[:,7]  + Zm_gas[:,8]  + Zm_gas[:,9]  + Zm_gas[:,10]  + Zm_gas[:,11] ) / m_gas[:]
        





###########################################
    """
    cmx = np.sum(x_disk*m_disk) / np.sum(m_disk)  
    cmy = np.sum(y_disk*m_disk) / np.sum(m_disk)  
    cmz = np.sum(z_disk*m_disk) / np.sum(m_disk)  
 
    vcmx = np.sum(vx_disk*m_disk) / np.sum(m_disk)
    vcmy = np.sum(vy_disk*m_disk) / np.sum(m_disk)
    vcmz = np.sum(vz_disk*m_disk) / np.sum(m_disk) 
n

    print 'disk+gas cm position = ', cmx, cmy, cmz
    print 'disk+gas cm velocity = ', vcmx,vcmy,vcmz

    
    x_gas -= cmx
    y_gas -= cmy
    z_gas -= cmz
    vx_gas -= vcmx
    vy_gas -= vcmy
    vz_gas -= vcmz

    x_dm -= cmx
    y_dm -= cmy
    z_dm -= cmz
    vx_dm -= vcmx
    vy_dm -= vcmy
    vz_dm -= vcmz
    
    x_disk -= cmx
    y_disk -= cmy
    z_disk -= cmz
    vx_disk -= vcmx
    vy_disk -= vcmy
    vz_disk -= vcmz

    if(Nstar > 0):
       x_star -= cmx
       y_star -= cmy
       z_star -= cmz
       vx_star -= vcmx
       vy_star -= vcmy
       vz_star -= vcmz
    """   
###########################################



    r2d_gas = np.sqrt(x_gas**2+y_gas**2)

#    f_dust=rs.read_block(filename, "SHDU", parttype=0)


#    tracAbundOut  = rs.read_block(filename, "CHEM", parttype=0) 
    """
    x_h2 = tracAbundOut[:,0]
    x_hp = tracAbundOut[:,1]
    x_co = tracAbundOut[:,2]

    x_HI = 1. - x_hp - 2.*x_h2
    x_cp = Zm_gas[:,1] / 12. / Zm_gas[:,6] - x_co
    x_si = Zm_gas[:,5] / 28. / Zm_gas[:,6]
    x_o  = Zm_gas[:,3] / 16. / Zm_gas[:,6]

    x_e = x_hp + x_cp + x_si
    """
    dust_to_gas_ratio = 0.1
    Ns = 1e5
    idplot =  np.random.random(Ns) * (Ngas-1) 
    idplot = idplot.astype(int)

    n_pccm = rho * UnitDensity_in_pccm *XH  #####gas density
    """
    plt.subplot(211)
    plt.plot(np.log10(n_pccm[idplot]), np.log10(chemT[idplot]), '.', markersize=0.1, c='black')
    """

    velocity = np.sqrt( vx_gas**2 + vy_gas**2 + vz_gas**2 ) #gas velocity
    """
    if(Nstar > 0):
        M_min = 1
        M_max = 50.
        Nbin=70
        binsize = np.float(np.log10(M_max) - np.log10(M_min)) / Nbin

        #M_clus = np.sum(np.abs(m_imf)) * (4.3441 / 2.3025) #account for the M<1M_sun stars that were discarded          
        M_clus = np.sum(np.abs(m_imf)) / 0.385224 #account for the M<1M_sun stars that were discarded      
        N_clus = M_clus * 1.819
    
        normal=False
        plt.hist(np.log10(np.abs(m_imf[m_imf!=0])), bins=Nbin,normed=normal, alpha=0.5, color='blue' , range=[np.log10(1), np.log10(50)], histtype='step', label='initial mass', linewidth=2)
        plt.hist(np.log10(np.abs(m_imf[m_imf>0])),  bins=Nbin,normed=normal, alpha=0.5, color='red'  , range=[np.log10(1), np.log10(50)], histtype='step', label='current mass', linewidth=2)
        plt.yscale('log')
        plt.axis([0, np.log10(50), 10**-(0.2), 1e4])

        aaa=np.array([M_min, 0.5, M_max])
        bbb=N_clus*np.log(10.)*binsize * np.array([aaa[0] * drawMassFromIMF(aaa[0]), aaa[1] * drawMassFromIMF(aaa[1]), aaa[2] * drawMassFromIMF(aaa[2])])
        plt.plot(np.log10(aaa), bbb, '-', markersize=2.0, linewidth=2, color='black')
    """




    """
    n_cp = x_cp * n_pccm       #####number density of C+
    n_e  = x_e  * n_pccm       #####number density of electrons
    n_HI = x_HI * n_pccm       #####number density of atomic hydrogen
    n_h2 = x_h2 * n_pccm       #####number density of molecular hydrogen



    gamma = G0 * f_dust * chemT**0.5 / n_e
    eps = 0.049 / (1.0 + 0.004*gamma**0.73) + (0.037*(chemT/10000)**0.7) / (1.0 + 2e-4*gamma)
    PEheat = 1.3e-24*eps*dust_to_gas_ratio*G0*f_dust
    """
#    lam        = rs.read_block(filename, "CHC1", parttype=0)
#    lam_chem   = rs.read_block(filename, "CHC2", parttype=0)
#    dustT   = rs.read_block(filename, "DUST", parttype=0)


#    photoelec_heat = -n_pccm*lam[:,11]
#    CII_cool       =  n_pccm*lam[:,15]
#    OI_cool        =  n_pccm*lam[:,12]


#    CII    = np.zeros(Ngas)
#    OI_63  = np.zeros(Ngas)
#    OI_145 = np.zeros(Ngas)
#    for i in range(Ngas):
#        CII[i]    = wrapper.get_cii_spec_emiss(   chemT[i], rho[i]*UnitDensity_in_cgs, x_cp[i], x_o[i], x_HI[i], x_h2[i], x_e[i], x_hp[i], G0[i])
#        OI_63[i]  = wrapper.get_oi_spec_emiss_63( chemT[i], rho[i]*UnitDensity_in_cgs, x_cp[i], x_o[i], x_HI[i], x_h2[i], x_e[i], x_hp[i], G0[i])
#        OI_145[i] = wrapper.get_oi_spec_emiss_145(chemT[i], rho[i]*UnitDensity_in_cgs, x_cp[i], x_o[i], x_HI[i], x_h2[i], x_e[i], x_hp[i], G0[i])
    """
    print 'calculating CII...'
    CII    = wrapper2.get_cii_spec_emiss_array(chemT, rho*XH*UnitDensity_in_cgs, x_cp, x_o, x_HI, x_h2, x_e, x_hp, G0, Ngas)
    print 'calculating OI 63...'
    OI_63  = wrapper2.get_oi_spec_emiss_63_array(chemT, rho*XH*UnitDensity_in_cgs, x_cp, x_o, x_HI, x_h2, x_e, x_hp, G0, Ngas)
    print 'calculating OI 145...'
    OI_145 = wrapper2.get_oi_spec_emiss_145_array(chemT, rho*XH*UnitDensity_in_cgs, x_cp, x_o, x_HI, x_h2, x_e, x_hp, G0, Ngas)
    """

#    dust_cool = 4.68e-31 * dustT**6 * n_pccm
#    l_dust = dust_to_gas_ratio * dust_cool / n_pccm
#    l_dust = dust_to_gas_ratio * dust_cool / n_pccm * (m_gas * UnitMass_in_g / PROTONMASS )

#    L_dust[counter] = np.sum( l_dust[idx_ism] ) / SolarLuminosity

#    L_CII[counter]   = np.sum(CII[idx_ism] * (m_gas[idx_ism] * UnitMass_in_g) )    / SolarLuminosity
#    L_OI63[counter]  = np.sum(OI_63[idx_ism] * (m_gas[idx_ism] * UnitMass_in_g) )  / SolarLuminosity
#    L_OI145[counter] = np.sum(OI_145[idx_ism] * (m_gas[idx_ism] * UnitMass_in_g) ) / SolarLuminosity
    



    """
    plt.subplot(212)
    plt.plot(n_pccm[idplot], n_pccm[idplot], '.', markersize=0.1, c='black')
    plt.plot(n_pccm[idplot], n_cp[idplot], '.', markersize=0.1, c='green')
    plt.plot(n_pccm[idplot], n_e[idplot] , '.', markersize=0.1, c='red')
    plt.plot(n_pccm[idplot], n_HI[idplot], '.', markersize=0.1, c='blue')
    plt.plot(n_pccm[idplot], n_h2[idplot], '.', markersize=0.1, c='cyan')
    plt.axis([1e-10,1e5,1e-10,1e5])

    plt.xscale('log')  
    plt.yscale('log')
    
    plt.show()
    """

"""
G_code = 43022.
N = 200    
dr = 0.005
pot_gas = np.zeros(N)
pot_disk = np.zeros(N)
pot_dm = np.zeros(N)
zzz = np.zeros(N)

for i in range(N):
    zzz[i] = i*dr
    dist_gas = np.sqrt( (x_gas - 0)**2 + (y_gas - 0)**2 + (z_gas - zzz[i])**2 )
    pot_gas[i] = -np.sum( m_gas / dist_gas )
    dist_disk = np.sqrt( (x_disk - 0)**2 + (y_disk - 0)**2 + (z_disk - zzz[i])**2 )
    pot_disk[i] = -np.sum( m_disk / dist_disk )
    dist_dm = np.sqrt( (x_dm - 0)**2 + (y_dm - 0)**2 + (z_dm - zzz[i])**2 )
    pot_dm[i] = -np.sum( m_dm / dist_dm )


pot_gas  *= G_code
pot_disk *= G_code
pot_dm   *= G_code


plt.plot(zzz, -pot_dm, label='dark matter')
plt.plot(zzz, -pot_gas, label='gas')        
plt.plot(zzz, -pot_disk, label='disk')      
plt.yscale('log')
plt.legend(loc='best')
plt.xlabel('z (kpc)', fontsize=16)
plt.ylabel('-potential [1e10 M_sun / kpc]', fontsize=16)



acc_gas = np.zeros(N-1)
acc_disk = np.zeros(N-1)
acc_dm = np.zeros(N-1)
for i in range(N-1):
    acc_gas[i]  =  (pot_gas[i+1] - pot_gas[i])  / dr    
    acc_disk[i] = (pot_disk[i+1] - pot_disk[i]) / dr    
    acc_dm[i]   =   (pot_dm[i+1] - pot_dm[i])   / dr

plt.plot(zzz[:-1]+0.5*dr, acc_dm, label='dark matter')
plt.plot(zzz[:-1]+0.5*dr, acc_gas, label='gas')
plt.plot(zzz[:-1]+0.5*dr, acc_disk, label='disk')
plt.yscale('log')
plt.legend(loc='best')
plt.xlabel('z (kpc)', fontsize=16)
plt.ylabel('acc [1e10 M_sun / (kpc**2)]', fontsize=16)
"""
print sfr_tot

time = jump*np.arange(1000/jump)

plt.plot(time, sfr_tot)
