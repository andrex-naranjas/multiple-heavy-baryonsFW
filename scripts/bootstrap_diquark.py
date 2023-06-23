#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
---------------------------------------------------------------------------
 Script to obtain uncertainties of heavy baryon mass spectrum via bootstrap
 Authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch) and
          H. Garcia-Tecocoatzi
---------------------------------------------------------------------------
"""
import sys
import os
from iminuit import Minuit
import numpy as np
import datetime
import pandas as pd
# framework modules
import bottomfw.baryons.data_preparation as dp
from bottomfw.baryons.bottom_diquark import BottomDiquark


if len(sys.argv) <= 1:
    sys.exit('Provide bottom states group name. Try again!')

#states = 'omega' # All, omega, cascades, sigmaLamb
# states = sys.argv[1]
run_baryons = sys.argv[1]

# for running batch jobs with htcondor
batch_number = None
if len(sys.argv) == 4:
    batch_number = sys.argv[2]
    workpath = sys.argv[3]
else:
    workpath = os.getcwd()

fit_type = "trad" # sext or trip, All, trad

# input parameters
param_v,param_w,param_x,param_y,param_z,\
    param_is_omega,param_is_cascade,param_is_sigma = dp.fetch_data_diquark(fit_type=fit_type)


def model(is_omega, is_cascade, is_sigma, v, w, x, y, z, md1, md2, md3, mb, k, a, b, e, g):
    return is_omega*md1 + is_cascade*md2 + is_sigma*md3 + mb + \
        v*k*np.sqrt(1./( (is_omega*md1 + is_cascade*md2 + is_sigma*md3)*mb / ( is_omega*md1 + is_cascade*md2 + is_sigma*md3 + mb ) ) ) +\
        w*a + x*b + y*e + z*g

def least_squares(md1, md2, md3, mb, k, a, b, e, g):
    # y_var_0 = sigma_0 # best sigma_0= 13.63
    # yvar_0 = y_var_0*np.ones(17)
    # yvar = y_errors_exp
    # yvar_2 = np.power(yvar_0, 2) + np.power(yvar, 2)
    yvar_2 = 0.001
    pred_m = model(param_is_omega, param_is_cascade, param_is_sigma,
                   param_v, param_w, param_x, param_y, param_z,
                   md1, md2, md3, mb, k, a, b, e, g)
    yval_2 = np.power( (pred_m - exp_m), 2)
    return np.sum( np.divide(yval_2, yvar_2) )
    # return np.sum((pred_m - exp_m)**2 / (yvar_2**2)) #**2

def fit(least_squares):
    m = Minuit(least_squares, md1=900, md2=300, md3=250, mb=1400, k=0, a=0, b=0, e=0, g=0)
    # m.limits['mb'] = (4000,6000)

    m.limits['md1'] =  (900, 1200)  # (850, 950) #omega (500, 1500)
    m.limits['md2'] =  (700, 950)   # (650, 850) #cascade prime
    m.limits['md3'] =  (500, 700)   # (500, 700) #sigma
    # m.limits['md1'] = (850,950)
    # m.limits['md2'] = (650,850)
    # m.limits['md3'] = (500,700)
    
    
    m.errordef = Minuit.LEAST_SQUARES
    m.migrad()
    return m

def sample_gauss(mu, sigma):
    return np.random.normal(mu, sigma, 10000)

def random(sample):
    return np.random.choice(sample, size=None)

# arrays to store the sampled parameters
sampled_k,sampled_a,sampled_b,sampled_e,sampled_g = ([]),([]),([]),([]),([])
sampled_md1,sampled_md2,sampled_md3,sampled_mb = ([]),([]),([]),([])

# arrays to store sampled correlation coeficients
rho_md2md1,rho_md3md1,rho_mbmd1,rho_kmd1,rho_amd1,rho_bmd1 = ([]),([]),([]),([]),([]),([])
rho_emd1,rho_gmd1,rho_md3md2,rho_mbmd2,rho_kmd2,rho_amd2   = ([]),([]),([]),([]),([]),([])
rho_bmd2,rho_emd2,rho_gmd2,rho_mbmd3,rho_kmd3,rho_amd3     = ([]),([]),([]),([]),([]),([])
rho_bmd3,rho_emd3,rho_gmd3,rho_kmb,rho_amb,rho_bmb         = ([]),([]),([]),([]),([]),([])
rho_emb,rho_gmb,rho_ak,rho_bk,rho_ek,rho_gk                = ([]),([]),([]),([]),([]),([])
rho_ba,rho_ea,rho_ga,rho_eb,rho_gb,rho_ge                  = ([]),([]),([]),([]),([]),([])


# start bootstrap
start = datetime.datetime.now()
sigma_model = 10**2 # to be obtained with optimization (Li.Jin)
# gaussian pdf with the measured value and with experimental and model(sigma_model) uncertainties
# Omega states
gauss_6061 = sample_gauss(6045.2, np.power((0.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2702), PDG
gauss_6316 = sample_gauss(6315.6, np.power((0.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2767), PDG
gauss_6330 = sample_gauss(6330.3, np.power((0.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 3015), PDG
gauss_6340 = sample_gauss(6339.7, np.power((0.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 3044), PDG
gauss_6350 = sample_gauss(6349.8, np.power((0.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 3051), PDG
# Cascade b sextet
gauss_5935 = sample_gauss(5935.0, np.power((0.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2461), PDG Average
gauss_5953 = sample_gauss(5953.8, np.power((0.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2796), PDG Average
gauss_6328 = sample_gauss(6227.9, np.power((0.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2832), PDG Average
# Sigma b
gauss_5813 = sample_gauss(5813.1, np.power((0.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2453), PDG Average
gauss_5837 = sample_gauss(5837.0, np.power((0.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2517), PDG Average
gauss_6097 = sample_gauss(6096.9, np.power((0.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2819), PDG Average
# Lambda b
gauss_5617 = sample_gauss(5619.6, np.power((0.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2283), PDG
gauss_5912 = sample_gauss(5912.2, np.power((0.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2649), PDG
gauss_5920 = sample_gauss(5920.1, np.power((0.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2685), PDG
gauss_6146 = sample_gauss(6146.2, np.power((0.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2283), PDG
gauss_6152 = sample_gauss(6152.5, np.power((0.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2649), PDG
gauss_6070 = sample_gauss(6072.3, np.power((0.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2685), PDG
# Cascades anti-3-plet
gauss_5794 = sample_gauss(5794.5, np.power((0.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2570), PDG Average
gauss_6100 = sample_gauss(6100.0, np.power((0.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2635), PDG Average
gauss_6327 = sample_gauss(6327.0, np.power((0.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2934), LHCb PhysRevLett
gauss_6333 = sample_gauss(6333.0, np.power((0.00**2 + sigma_model), 0.5 ))  # all OK (corresponds to predicted 2941), LHCb PhysRevLett

# plug here the sigma_0 optimization lines from data_utils.py

# construct the simulated sampling distribution (bootstrap technique)
# for _ in range(10000): # max 10000 with decays included, computationally expensive

count = 0
while count < 10000: # max 10000 with decays included, computationally expensive

    if fit_type=="All" or fit_type=="trad":
        exp_m = np.array([ # measured baryon masses        
            # omegas
            random(gauss_6061),
            random(gauss_6316),
            random(gauss_6330),
            random(gauss_6340),
            random(gauss_6350),
            # Cascade 
            random(gauss_5935),
            random(gauss_5953),
            random(gauss_6328),
            # Sigma b         
            random(gauss_5813),
            random(gauss_5837),
            random(gauss_6097),        
            # Lambda b
            random(gauss_5617),
            random(gauss_5912),
            random(gauss_5920),
            # random(gauss_6146),
            # random(gauss_6152),
            # Cascades
            random(gauss_5794),
            random(gauss_6100),
            # random(gauss_6327),
            # random(gauss_6333)        
        ])
    elif fit_type=="sext":
        exp_m = np.array([ # measured baryon masses        
            # omegas
            random(gauss_6061),
            random(gauss_6316),
            random(gauss_6330),
            random(gauss_6340),
            random(gauss_6350),
            # Cascade 
            random(gauss_5935),
            random(gauss_5953),
            random(gauss_6328),
            # Sigma b         
            random(gauss_5813),
            random(gauss_5837),
            random(gauss_6097),
        ])
    elif fit_type=="trip":
        exp_m = np.array([ # measured baryon masses        
            # Lambda b
            random(gauss_5617),
            random(gauss_5912),
            random(gauss_5920),
            # random(gauss_6146),
            # random(gauss_6152),
            # Cascades antitriplet
            random(gauss_5794),
            random(gauss_6100),
            # random(gauss_6327),
            # random(gauss_6333)        
        ])
        
    # print(len(exp_m))
    # input()
    # perform the parameter fitting (via minimizing squared distance)
    m = fit(least_squares)

    if type(m.covariance) != type(None):
        count += 1
    else:
        continue

    sampled_md1 = np.append(sampled_md1, m.values['md1'])
    sampled_md2 = np.append(sampled_md2, m.values['md2'])
    sampled_md3 = np.append(sampled_md3, m.values['md3'])
    sampled_mb  = np.append(sampled_mb,  m.values['mb'])

    sampled_k = np.append(sampled_k, m.values['k'])
    sampled_a = np.append(sampled_a, m.values['a'])
    sampled_b = np.append(sampled_b, m.values['b'])
    sampled_e = np.append(sampled_e, m.values['e'])
    sampled_g = np.append(sampled_g, m.values['g'])

    # correlation matrix
    corr = m.covariance.correlation()

    rho_md2md1 = np.append(rho_md2md1, corr['md2','md1'])
    rho_md3md1 = np.append(rho_md3md1, corr['md3','md1'])
    rho_mbmd1  = np.append(rho_mbmd1, corr['mb','md1'])
    rho_kmd1   = np.append(rho_kmd1,  corr['k','md1'])
    rho_amd1   = np.append(rho_amd1,  corr['a','md1'])
    rho_bmd1   = np.append(rho_bmd1,  corr['b','md1'])
    rho_emd1   = np.append(rho_emd1,  corr['e','md1'])
    rho_gmd1   = np.append(rho_gmd1,  corr['g','md1'])

    rho_md3md2 = np.append(rho_md3md2, corr['md3','md2'])
    rho_mbmd2  = np.append(rho_mbmd2,  corr['mb','md2'])
    rho_kmd2   = np.append(rho_kmd2 , corr['k','md2'])
    rho_amd2   = np.append(rho_amd2 , corr['a','md2'])
    rho_bmd2   = np.append(rho_bmd2 , corr['b','md2'])
    rho_emd2   = np.append(rho_emd2 , corr['e','md2'])
    rho_gmd2   = np.append(rho_gmd2 , corr['g','md2'])

    rho_mbmd3  = np.append(rho_mbmd3,corr['mb','md3'])
    rho_kmd3   = np.append(rho_kmd3, corr['k','md3'])
    rho_amd3   = np.append(rho_amd3, corr['a','md3'])
    rho_bmd3   = np.append(rho_bmd3, corr['b','md3'])
    rho_emd3   = np.append(rho_emd3, corr['e','md3'])
    rho_gmd3   = np.append(rho_gmd3, corr['g','md3'])

    rho_kmb    = np.append(rho_kmb, corr['k','mb'])
    rho_amb    = np.append(rho_amb, corr['a','mb'])
    rho_bmb    = np.append(rho_bmb, corr['b','mb'])
    rho_emb    = np.append(rho_emb, corr['e','mb'])
    rho_gmb    = np.append(rho_gmb, corr['g','mb'])
    
    rho_ak     = np.append(rho_ak, corr['a','k'])
    rho_bk     = np.append(rho_bk, corr['b','k'])
    rho_ek     = np.append(rho_ek, corr['e','k'])
    rho_gk     = np.append(rho_gk, corr['g','k'])

    rho_ba     = np.append(rho_ba, corr['b','a'])
    rho_ea     = np.append(rho_ea, corr['e','a'])
    rho_ga     = np.append(rho_ga, corr['g','a'])

    rho_eb     = np.append(rho_eb, corr['e','b'])
    rho_gb     = np.append(rho_gb, corr['g','b'])
    
    rho_ge     = np.append(rho_ge, corr['g','e'])
    
# save bootstrap results
df = pd.DataFrame({"Md1" : sampled_md1,"Md2" : sampled_md2,"Md3" : sampled_md3,
                   "MB"  : sampled_mb,"K" : sampled_k, "A" : sampled_a,
                   "B": sampled_b,    "E" : sampled_e, "G" : sampled_g})
if batch_number is None:
    if not os.path.exists(workpath+"/tables/"):
        os.mkdir(workpath+"/tables/")
    df.to_csv(workpath+"/tables/bootstrap_param_diquark_"+run_baryons+".csv", index=False)
else:
    if not os.path.exists(workpath+"/batch_results_diquark/"+run_baryons+"/parameters/"):
        os.mkdir(workpath+"/batch_results_diquark/"+run_baryons+"/parameters/")
    df.to_csv(workpath+"/batch_results_diquark/"+run_baryons+"/parameters/"+str(batch_number)+".csv", index=False)

# create dictionaries
param   = {'is_omega':param_is_omega,
           'is_cascade':param_is_cascade, 'is_sigma':param_is_sigma,'V':param_v, 'W':param_w, 'X':param_x, 'Y':param_y, 'Z':param_z}

sampled = {'sampled_md1':sampled_md1,'sampled_md2':sampled_md2,'sampled_md3':sampled_md3,'sampled_mb':sampled_mb,
           'sampled_k':sampled_k,'sampled_a':sampled_a, 'sampled_b':sampled_b, 'sampled_e':sampled_e, 'sampled_g':sampled_g}
corr_mat_diquark ={
    'rho_md2md1':rho_md2md1,'rho_md3md1':rho_md3md1,'rho_mbmd1':rho_mbmd1, 'rho_kmd1':rho_kmd1, 'rho_amd1':rho_amd1, 'rho_bmd1':rho_bmd1,
    'rho_emd1':rho_emd1, 'rho_gmd1':rho_gmd1, 'rho_md3md2':rho_md3md2, 'rho_mbmd2':rho_mbmd2, 'rho_kmd2':rho_kmd2, 'rho_amd2':rho_amd2,
    'rho_bmd2':rho_bmd2, 'rho_emd2':rho_emd2, 'rho_gmd2':rho_gmd2, 'rho_mbmd3':rho_mbmd3, 'rho_kmd3':rho_kmd3, 'rho_amd3':rho_amd3,
    'rho_bmd3':rho_bmd3, 'rho_emd3':rho_emd3, 'rho_gmd3':rho_gmd3, 'rho_kmb':rho_kmb, 'rho_amb':rho_amb,'rho_bmb':rho_bmb,
    'rho_emb':rho_emb,'rho_gmb':rho_gmb,'rho_ak':rho_ak,'rho_bk':rho_bk,'rho_ek':rho_ek,'rho_gk':rho_gk,
    'rho_ba':rho_ba, 'rho_ea':rho_ea, 'rho_ga':rho_ga, 'rho_eb':rho_eb, 'rho_gb':rho_gb, 'rho_ge':rho_ge}

# save bootstrap correlation parameters
df = pd.DataFrame(corr_mat_diquark)
if batch_number is None:
    if not os.path.exists(workpath+"/tables/"):
        os.mkdir(workpath+"/tables/")
    df.to_csv(workpath+"/tables/bootstrap_correlation_diquark_"+run_baryons+".csv", index=False)
else:
    if not os.path.exists(workpath+"/batch_results_diquark/"+run_baryons+"/correlation/"):
        os.mkdir(workpath+"/batch_results_diquark/"+run_baryons+"/correlation/")
    df.to_csv(workpath+"/batch_results_diquark/"+run_baryons+"/correlation/"+str(batch_number)+".csv", index=False)

# calculate the results using bootstrap simulation above
results = BottomDiquark(baryons=run_baryons, params=param, sampled=sampled, corr_mat=corr_mat_diquark, asymmetric=True, batch_number=batch_number, workpath=workpath)
results.fetch_values()

print('Getting paper results for:', run_baryons)
# # omegas,cascades,sigmas,lambdas,cascades_anti3
results.paper_results_predictions(bootstrap=True, bootstrap_width=False, prev_params=False, decay_width=False)

end = datetime.datetime.now()
elapsed_time = end - start
print("Elapsed total time = " + str(elapsed_time))
print(round(sampled_md1.mean()), "md1 omega")
print(round(sampled_md2.mean()), "md2 cascade")
print(round(sampled_md3.mean()), "md3 sigma, lambda")
print(round(sampled_mb.mean()), "mb")

print("K", pow(sampled_k.mean(), 2)/(1000**3))
print("A", sampled_a.mean())
print("B", sampled_b.mean())
print("E", sampled_e.mean())
print("G", sampled_g.mean())
