import matplotlib   # import matplotlib for next statement
matplotlib.use('Agg')   # this agg backend for plotting supports pdf, pngimport numpy as np

import os
import sys
import math
import readsnap as rs
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import LogNorm


from pylab import * #import pylab for ioff()
ioff()  # set interactive to off so no plotting to x-window


#%reset


filename_path="/mnt/ceph/users/chu/snapshots/dwarf_chem/ng1e7HI4e7gsl2_cutDM_PE_PI_SN_localSh_eps1_soft05pc_sIMF2myr4pc_wsscie_adpsoft"                                                                                                      
#filename_path="/mnt/ceph/users/chu/snapshots/dwarf_chem/ng1e7HI4e7gsl2_cutDM_PE_PI_SN_localSh_eps1e20_soft05pc_sIMF1myr10pc" 



filename_base = filename_path + "/snap_"

varG0 = True
#varG0 = False
print 'varG0 = ', varG0

if(varG0 == False):
    G0 = 0.
    print 'G0 = ', G0

save_data = True
#save_data = False



dust_to_gas_ratio = 0.1
print ' dust_to_gas_ratio = ', dust_to_gas_ratio

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

SolarLuminosity = 3.839e33

#rs.list_format2_blocks(filename)


jump=10

N_snap = 1000 / jump

sfr_total = np.zeros(N_snap)
H2_frac   = np.zeros(N_snap)
H2_frac_eq   = np.zeros(N_snap)
CO_frac   = np.zeros(N_snap)
time      = np.zeros(N_snap)

ofr_3     = np.zeros(N_snap)
ifr_3     = np.zeros(N_snap)
ofr_2     = np.zeros(N_snap)
ifr_2     = np.zeros(N_snap)

ofr_vir = np.zeros(N_snap)
ifr_vir = np.zeros(N_snap)

vel_in_2    = np.zeros(N_snap)
vel_out_2   = np.zeros(N_snap)
vel_out_vir = np.zeros(N_snap)

Zm_out_vir = np.zeros(N_snap)

mass_disc   = np.zeros(N_snap)
mass_halo   = np.zeros(N_snap)
mass_escape = np.zeros(N_snap)

Zm_disc   = np.zeros(N_snap)
Zm_halo   = np.zeros(N_snap)
Zm_escape = np.zeros(N_snap)

sfr_total_sc = np.zeros(N_snap)

M_frac_cold = np.zeros(N_snap)
M_frac_hot  = np.zeros(N_snap)
M_frac_warm = np.zeros(N_snap)

V_frac_cold = np.zeros(N_snap)
V_frac_hot  = np.zeros(N_snap)
V_frac_warm = np.zeros(N_snap)
V_frac_hot_r2z0p5 = np.zeros(N_snap)
V_frac_hot_r2z1 = np.zeros(N_snap)


F_H2_diff = np.zeros(N_snap)
F_H2_SF = np.zeros(N_snap)

f_PDR = np.zeros(N_snap)
f_dense_PDR = np.zeros(N_snap)
mean_G0 = np.zeros(N_snap)

temp_of = np.zeros(N_snap)


m_gas_SF = np.zeros(N_snap)

L_dust = np.zeros(N_snap)
L_CII = np.zeros(N_snap)
L_OI = np.zeros(N_snap)

F_gas_SF = np.zeros(N_snap)

thisSnap = 182

for k in range(0, N_snap):
#for k in range(thisSnap, thisSnap+1):
    kk = k*jump
    if (kk < 10):
        num = '00' + str(kk)
    elif (kk < 100):
        num = '0' + str(kk)
    elif (kk < 1000):
        num = str(kk)

    filename = filename_base + num
    plt.clf()

    print ''
    print 'Read snaphot ', num
    head= rs.snapshot_header(filename)

    Ngas   = head.npart[0]
    Ndm    = head.npart[1]
    Ndisk  = head.npart[2]
    Nbulge = head.npart[3]
    Nstar  = head.npart[4]
    


    #-------- Dark matter --------
    """
    if (Ndm > 0):
        pos_dm = rs.read_block(filename, "POS ", parttype=1) 
        x_dm = pos_dm[:,0]
        y_dm = pos_dm[:,1]
        z_dm = pos_dm[:,2]

        vel_dm = rs.read_block(filename, "VEL ", parttype=1) 
        vx_dm = vel_dm[:,0]
        vy_dm = vel_dm[:,1]
        vz_dm = vel_dm[:,2]
    """
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

#        id_star  = rs.read_block(filename, "ID  ", parttype=4) 
        m_star   = rs.read_block(filename, "MASS", parttype=4) 
#        Zm_star = rs.read_block(filename, "Z   ", parttype=4, csformat = 1) 
#        Zm_star_total = (Zm_star[:,1] + Zm_star[:,2] + Zm_star[:,3] + Zm_star[:,4] + Zm_star[:,5] + Zm_star[:,7] + Zm_star[:,8] + Zm_star[:,9] + Zm_star[:,10] + Zm_star[:,11]) / m_star[:]
	age   = head.time - rs.read_block(filename, "AGE ", parttype=4)


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
    
#        id_gas  = rs.read_block(filename, "ID  ", parttype=0) 
        m_gas   = rs.read_block(filename, "MASS", parttype=0) 
        
#        u       = rs.read_block(filename, "U   ", parttype=0) 
        rho     = rs.read_block(filename, "RHO ", parttype=0) 
#        ne      = rs.read_block(filename, "NE  ", parttype=0) 
#        nh      = rs.read_block(filename, "NH  ", parttype=0) 
#        hsml    = rs.read_block(filename, "HSML", parttype=0) 
        sfr_pp  = rs.read_block(filename, "SFR ", parttype=0)   #star formation rate per particle
        Zm_gas  = rs.read_block(filename, "Z   ", parttype=0, csformat = 1) 
        #cs_temp = rs.read_block(filename, "CSTE", parttype=0) 
        chemT   = rs.read_block(filename, "CHET", parttype=0) 
        #sfr_ff  = rs.read_block(filename, "SFFF", parttype=0)   #star formation rate per particle
        v_disp  = rs.read_block(filename, "CSSI", parttype=0)   
	#coln  = rs.read_block(filename, "COLN", parttype=0)   
        if(varG0 == True):
            G0      = rs.read_block(filename, "G0  ", parttype=0)

	n_pccm = rho * UnitDensity_in_pccm * XH
        n_H = n_pccm

	#col = np.zeros(Ngas)
	#col[:] = rho[:] * hsml[:] + coln[:,0]+coln[:,1]+coln[:,2]+coln[:,3]+coln[:,4]+coln[:,5]+coln[:,6]+coln[:,7]+coln[:,8]+coln[:,9]+coln[:,10]+coln[:,11]
	#col_local = rho * hsml
        
	#col       *= UnitDensity_in_pccm * UnitLength_in_cm / (1.0 + 4.0 * 0.1)
        #col_local *= UnitDensity_in_pccm * UnitLength_in_cm / (1.0 + 4.0 * 0.1)


        Zm_gas_total  = (Zm_gas[:,1]  + Zm_gas[:,2]  + Zm_gas[:,3]  + Zm_gas[:,4]  + Zm_gas[:,5]  + Zm_gas[:,7]  + Zm_gas[:,8]  + Zm_gas[:,9]  + Zm_gas[:,10]  + Zm_gas[:,11] ) / m_gas[:]
        





###########################################
#    cmx = ( np.sum(x_disk*m_disk) + np.sum(x_gas*m_gas) ) / ( np.sum(m_disk) + np.sum(m_gas) )
    if(Ndisk > 0):
        
        cmx = np.sum(x_disk*m_disk) / np.sum(m_disk)  
        cmy = np.sum(y_disk*m_disk) / np.sum(m_disk)  
        cmz = np.sum(z_disk*m_disk) / np.sum(m_disk)  
        
        vcmx = np.sum(vx_disk*m_disk) / np.sum(m_disk)
        vcmy = np.sum(vy_disk*m_disk) / np.sum(m_disk)
        vcmz = np.sum(vz_disk*m_disk) / np.sum(m_disk) 
        
        
        print 'disk+gas cm position = ', cmx, cmy, cmz
        print 'disk+gas cm velocity = ', vcmx,vcmy,vcmz
        
        
        x_gas -= cmx
        y_gas -= cmy
        z_gas -= cmz
        vx_gas -= vcmx
        vy_gas -= vcmy
        vz_gas -= vcmz
        
        
        
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
            
###########################################


    i_plot = np.random.randint(Ngas, size=1000 )

    """
    if(Nstar > 0):
        flux = np.zeros(Ngas)
        G0_direct = np.zeros(Ngas)
        for ii in range(Nstar):
            flux = m_star[ii] / ((x_star[ii]-x_gas)**2 + (y_star[ii]-y_gas)**2 + (z_star[ii]-z_gas)**2)
            G0_direct += flux
            
    """
    


    velocity = np.sqrt(vx_gas**2 + vy_gas**2 + vz_gas**2)

    


    tracAbundOut  = rs.read_block(filename, "CHEM", parttype=0) 

    x_h2 = tracAbundOut[:,0] * 2.0
    x_hp = tracAbundOut[:,1]
    x_co = tracAbundOut[:,2]

#    x_HI = 1. - x_hp - 2.*x_h2
    x_HI = 1. - x_hp - x_h2
#    x_cp = Zm_gas[:,1] / 12. / Zm_gas[:,6] - x_co
#    x_si = Zm_gas[:,5] / 28. / Zm_gas[:,6]
#    x_o  = Zm_gas[:,3] / 16. / Zm_gas[:,6]

#    x_e = x_hp + x_cp + x_si

    print 'total gas mass = ', np.sum(m_gas)

    inplane = np.abs(z_gas) < 1.
    print 'gas farction within +- 1kpc = ', np.sum(m_gas[inplane]) / np.sum(m_gas)
    
    #print 'total star mass = ', np.sum(m_star)

    #print 'total cold gas mass = ', np.sum(m_gas[id_cold])
    

    time[k] = head.time
    sfr_total[k] = np.sum(sfr_pp)   #total star formation rate
    

    tlife_ys = 5e-3

    if (Nstar > 0):
        sfr_total_sc[k] = np.sum(m_star[age < tlife_ys])*1e10 / (tlife_ys*1e9)


#    N_grid = 200  #number of grids

#    flux_in  = np.zeros([N_grid, N_grid])  #inflow flux
#    flux_out = np.zeros([N_grid, N_grid])  #outflow flux
    
 
    r2d_max = 2.    
    r2d_gas = np.sqrt( x_gas**2 + y_gas**2 )

    cylin = r2d_gas < r2d_max
    


    z_cr = 3.0
    dz = 0.1

    upper_in = (np.abs(z_gas - z_cr) < dz) & (vz_gas*z_cr < 0.0)
    lower_in = (np.abs(z_gas + z_cr) < dz) & (vz_gas*z_cr > 0.0)    
    idx_in   = (upper_in | lower_in) & cylin

    upper_out = (np.abs(z_gas - z_cr) < dz) & (vz_gas*z_cr > 0.0)
    lower_out = (np.abs(z_gas + z_cr) < dz) & (vz_gas*z_cr < 0.0)
    idx_out   = (upper_out | lower_out) & cylin

    ifr_3[k] = np.sum( m_gas[idx_in]  * np.abs(vz_gas[idx_in])  ) / dz 
    ofr_3[k] = np.sum( m_gas[idx_out] * np.abs(vz_gas[idx_out]) ) / dz 

    ifr_3[k] *= 1e10 * 1e5 / 3.08e21 * 31556926
    ofr_3[k] *= 1e10 * 1e5 / 3.08e21 * 31556926




    z_cr = 2.0
    dz = 0.1

    upper_in = (np.abs(z_gas - z_cr) < dz) & (vz_gas*z_cr < 0.0)
    lower_in = (np.abs(z_gas + z_cr) < dz) & (vz_gas*z_cr > 0.0)    
    idx_in   = (upper_in | lower_in) & cylin

    upper_out = (np.abs(z_gas - z_cr) < dz) & (vz_gas*z_cr > 0.0)
    lower_out = (np.abs(z_gas + z_cr) < dz) & (vz_gas*z_cr < 0.0)
    idx_out   = (upper_out | lower_out) & cylin

    ifr_2[k] = np.sum( m_gas[idx_in]  * np.abs(vz_gas[idx_in])  ) / dz 
    ofr_2[k] = np.sum( m_gas[idx_out] * np.abs(vz_gas[idx_out]) ) / dz 

    ifr_2[k] *= 1e10 * 1e5 / 3.08e21 * 31556926
    ofr_2[k] *= 1e10 * 1e5 / 3.08e21 * 31556926








#    temp_of[k] = np.sum( m_gas[idx_out] * chemT[idx_out] ) / np.sum( m_gas[idx_out] )

    if(ofr_2[k] > 0):
        temp_of[k] = 10**np.mean( np.log10(chemT[idx_out]) )

 
    if (len(m_gas[idx_in]) > 0):
        vel_in_2[k]  = np.sum( m_gas[idx_in]  * np.abs(vz_gas[idx_in])  ) / np.sum( m_gas[idx_in])

    if (len(m_gas[idx_out]) > 0):
        vel_out_2[k] = np.sum( m_gas[idx_out] * np.abs(vz_gas[idx_out]) ) / np.sum( m_gas[idx_out])




    #----------- R_vir ----------------------
    radius   = np.sqrt( x_gas**2 +  y_gas**2 +  z_gas**2)


    dr = 1.0
    r_vir = 44.0

    v_r = (vx_gas*x_gas + vy_gas*y_gas + vz_gas*z_gas) / radius


    idx_out_vir = (np.abs(radius - r_vir) < dr) & (v_r > 0.0)
    idx_in_vir  = (np.abs(radius - r_vir) < dr) & (v_r < 0.0)

    ifr_vir[k] = np.sum( m_gas[idx_in_vir]  * np.abs(v_r[idx_in_vir])  ) / dr 
    ofr_vir[k] = np.sum( m_gas[idx_out_vir] * np.abs(v_r[idx_out_vir]) ) / dr 

    ifr_vir[k] *= 1e10 * 1e5 / 3.08e21 * 31556926
    ofr_vir[k] *= 1e10 * 1e5 / 3.08e21 * 31556926


    if (len(m_gas[idx_out_vir]) > 0):
        vel_out_vir[k] = np.sum( m_gas[idx_out_vir] * np.abs(v_r[idx_out_vir]) ) / np.sum( m_gas[idx_out_vir])

    if (len(Zm_gas_total[idx_out_vir]) > 0):
        Zm_out_vir[k] = np.mean(Zm_gas_total[idx_out_vir])


    
#    r2d_discBC = 2.0 * r2d_max
    r2d_discBC = 6.0

    idx_disc = (r2d_gas < r2d_discBC) & (np.abs(z_gas) < z_cr)
    idx_halo = ( (r2d_gas > r2d_discBC) | (np.abs(z_gas) > z_cr) ) & (radius < r_vir)
    idx_escape = (radius > r_vir)

    mass_disc[k]   = np.sum(m_gas[idx_disc])
    mass_halo[k]   = np.sum(m_gas[idx_halo])
    mass_escape[k] = np.sum(m_gas[idx_escape])

    if (len(Zm_gas_total[idx_disc]) > 0):
        Zm_disc[k] = np.mean(Zm_gas_total[idx_disc])

    if (len(Zm_gas_total[idx_halo]) > 0):
        Zm_halo[k] = np.mean(Zm_gas_total[idx_halo])

    if (len(Zm_gas_total[idx_escape]) > 0):
        Zm_escape[k] = np.mean(Zm_gas_total[idx_escape])



    #------- Mass fraction -----------
    T_cold = 100
    T_hot  = 3e4
    idx_ism = (r2d_gas < 1.5) & (np.abs(z_gas) < 0.2)


    M_frac_cold[k] = np.sum(m_gas[idx_ism&(chemT<T_cold)]) / np.sum(m_gas[idx_ism])
    M_frac_hot[k]  = np.sum(m_gas[idx_ism&(chemT>T_hot)]) / np.sum(m_gas[idx_ism])
    M_frac_warm[k] = np.sum(m_gas[idx_ism&(chemT<=T_hot)&(chemT>=T_cold)]) / np.sum(m_gas[idx_ism])
    
    print 'cold/hot/warm mass fraction = ', M_frac_cold[k], M_frac_hot[k], M_frac_warm[k]

    #------- Volume-weighted Temperature PDF -----------
    T_PDF = plt.hist(np.log10(chemT[idx_ism]), bins=140, histtype='step', normed=True, range=[1,8], weights=1./rho[idx_ism])
#    T_PDF = plt.hist(np.log10(chemT[idx_disc]), bins=140, histtype='step', normed=True, range=[1,8], weights=1./rho[idx_disc])
    Tbin = T_PDF[1][:-1]
    Tcount = T_PDF[0]

    binsize = (Tbin[1]-Tbin[0])



    V_frac_cold[k] = np.sum(Tcount[(10**Tbin)<T_cold]) * binsize
    V_frac_hot[k]  = np.sum(Tcount[(10**Tbin)>T_hot]) * binsize
    V_frac_warm[k] = np.sum(Tcount[((10**Tbin)<=T_hot)&((10**Tbin)>=T_cold)]) * binsize
    print 'cold/hot/warm volume fraction = ', V_frac_cold[k], V_frac_hot[k], V_frac_warm[k]


    V_gas = m_gas / rho
    idx_ism = (r2d_gas < 2.0) & (np.abs(z_gas) < 0.5)
    V_frac_hot_r2z0p5[k] = np.sum( V_gas[idx_ism&(chemT>T_hot)] ) / np.sum(V_gas[idx_ism])
    print '!!!hot volume fraction r2z0p5 = ', V_frac_hot_r2z0p5[k]

    idx_ism = (r2d_gas < 2.0) & (np.abs(z_gas) < 1.0)
    V_frac_hot_r2z1[k] = np.sum( V_gas[idx_ism&(chemT>T_hot)] ) / np.sum(V_gas[idx_ism])
    print '???hot volume fraction r2z1 = ', V_frac_hot_r2z1[k]


    idx_ism = (r2d_gas < 2.0) | (np.abs(z_gas) < 1.0)
    plt.hist2d(np.log10(n_H[idx_ism]), np.log10(chemT[idx_ism]), range = [(-4, 2), (0, 8)], bins = 100, norm = LogNorm()); plt.xlabel('Density (cm^-3)'); plt.ylabel('Temperature (K)'); plt.title('Phase Diagram')
    
    if k < 10:
        plt.savefig('phase_diagram_cropped00{k}.png'.format(k=k))
    elif 10 <= k < 100:
        plt.savefig('phase_diagram_cropped0{k}.png'.format(k=k))

    else:
        plt.savefig('phase_diagram_cropped{k}.png'.format(k=k))

    """
    dz = 0.1
    plt.hist2d(np.log10(n_H[idx_out]), np.log10(chemT[idx_out]), range = [(-4, 2), (0, 8)], bins = 100, norm = LogNorm()); plt.xlabel('log_10 Density (cm^-3)'); plt.ylabel('Temperature (K)'); plt.title('Phase Diagram')
    
    if k < 10:
        plt.savefig('of_phase_diagram00{k}.png'.format(k=k))
    elif 10 < k < 100:
        plt.savefig('of_phase_diagram0{k}.png'.format(k=k))

    else:
        plt.savefig('of_phase_diagram{k}.png'.format(k=k))
    
    
    """


    plt.clf()

    print 'ifr_vir = ', ifr_vir[k]
    print 'ofr_vir = ', ofr_vir[k]
    print 'vel_out_vir = ', vel_out_vir[k]

    print 'ifr (2 kpc) = ', ifr_2[k]
    print 'ofr (2 kpc) = ', ofr_2[k]
    print 'mass loading (2 kpc) = ', ofr_2[k] / sfr_total[k]
    print 'outflow temperature = ', temp_of[k]

    print 'vel_in_2 = ', vel_in_2[k]
    print 'vel_out_2 = ', vel_out_2[k]


    print 'mass fraction (disc) = ', mass_disc[k] / np.sum(m_gas)
    print 'mass fraction (halo) = ', mass_halo[k] / np.sum(m_gas)
    print 'mass fraction (escape) = ', mass_escape[k] / np.sum(m_gas)
    print 'total fraction add up to :', ( mass_disc[k] + mass_halo[k] + mass_escape[k] ) / np.sum(m_gas)

    print 'Star formation rate = ', sfr_total[k], ' M_sol/yr'
    print 'Star formation rate (star count) = ', sfr_total_sc[k], ' M_sol/yr'
    
##########################################################################

    """
    if(save_data == True):
        print "Saving outputs..."
        np.savetxt('./soifr_time_jump1.txt', (time, sfr_total, ofr_vir, ifr_vir, ofr_2, ifr_2, H2_frac, H2_frac_eq, Zm_out_vir, vel_out_vir, vel_out_2, vel_in_2, mass_disc, mass_halo, mass_escape, Zm_disc, Zm_halo, Zm_escape, sfr_total_sc, M_frac_cold, M_frac_hot, M_frac_warm, V_frac_cold, V_frac_hot, V_frac_warm, F_H2_diff, F_H2_SF, f_PDR, f_dense_PDR, mean_G0, ofr_3, ifr_3, V_frac_hot_r2z0p5, V_frac_hot_r2z1, F_gas_SF) )
        print "done."
    """


time_axis = jump*np.arange(1000/jump)
plt.plot(time, ofr_2)
