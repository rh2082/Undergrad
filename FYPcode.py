#import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from lmfit import Parameters, Model
from astropy.io import fits

#define plancks constant and c as universal
h = 4.135e-18 # units of keV/Hz
c = 3e8 # units of m/s

# Dictionary with Pei (1992) parameters for each extinction curve
extcurves = dict(MW = dict(BKG = [165., 4.7e-2, 2., 90.], 
			   FUV = [14., 0.08e0, 6.5, 4.], 
                           two = [0.045, 0.22e0, 2., -1.95], 
                           nine = [0.002, 9.7e0, 2., -1.95], 
                           eighteen = [0.002, 18e0, 2., -1.8], 
                           FIR = [0.012, 25e0, 2., 0.]),
		 LMC = dict(BKG = [175., 4.6e-2, 2., 90.], 
                	    FUV = [19., 0.08e0, 4.5, 5.5], 
                            two = [0.023, 0.22e0, 2., -1.95], 
                            nine = [0.005, 9.7e0, 2., -1.95], 
                            eighteen = [0.006, 18e0, 2., -1.8], 
                            FIR = [0.03, 25e0, 2., 0.]),
        	 SMC = dict(BKG = [185., 4.2e-2, 2., 90.], 
                	    FUV = [27., 0.08e0, 4., 5.5], 
                            two = [0.005, 0.22e0, 2., -1.95], 
                            nine = [0.01, 9.7e0, 2., -1.95], 
                            eighteen = [0.012, 18e0, 2., -1.8], 
                            FIR = [0.03, 25e0, 2., 0.]))


#evaluate totsigma for attenuation in the xray region according to modelling by Wein
def evaluate_sigma():
    # Read in abundances
    ISMabund = {'H':12.0,'He':10.99,'C':8.38,'N':7.88,'O':8.69,'Fe':7.43,'Ne':7.94,'Mg':7.40,'Si':7.27,'S':7.09,'Ar':6.41,'Ca':6.20}
    dustdplt = {'H':1.0,'He':1.0,'C':0.5,'N':1.0,'O':0.6,'Fe':0.3,'Ne':1.0,'Mg':0.2,'Si':0.1,'S':0.6,'Ar':1.0,'Ca':0.003}
    #read fits file
    xsections = fits.getdata('IsmabsAtomicData_reduced.fits',1)
    # Read in required cross-sections
    H = xsections['H']*10.**(ISMabund['H']-12)*dustdplt['H']
    He = xsections['HeI']*10.**(ISMabund['He']-12)*dustdplt['He']
    C = xsections['CI']*10.**(ISMabund['C']-12)*dustdplt['C']
    N = xsections['NI']*10.**(ISMabund['N']-12)*dustdplt['N']
    O = xsections['OI']*10.**(ISMabund['O']-12)*dustdplt['O']
    Fe = xsections['FeI']*10.**(ISMabund['Fe']-12)*dustdplt['Fe']
    Ne = xsections['NeI']*10.**(ISMabund['Ne']-12)*dustdplt['Ne']
    Mg = xsections['MgI']*10.**(ISMabund['Mg']-12)*dustdplt['Mg']
    Si = xsections['SiI']*10.**(ISMabund['Si']-12)*dustdplt['Si']
    S = xsections['SI']*10.**(ISMabund['S']-12)*dustdplt['S']
    Ar = xsections['ArI']*10.**(ISMabund['Ar']-12)*dustdplt['Ar']
    Ca = xsections['CaI']*10.**(ISMabund['Ca']-12)*dustdplt['Ca']
    keV = xsections['Energy']/1000.
    #sum
    totsigma = (H+He+C+O+Fe+Ne+Mg+Si+S+Ar+Ca)
    return totsigma, keV

#define extinction curve as described by Pei
def extinction_curve(curve, wavelength):
    #generates the selected extinction curve
    BKG = extinction_fit(extcurves[curve]['BKG'], wavelength)
    FUV = extinction_fit(extcurves[curve]['FUV'], wavelength)
    two = extinction_fit(extcurves[curve]['two'], wavelength)
    nine = extinction_fit(extcurves[curve]['nine'], wavelength)
    eighteen = extinction_fit(extcurves[curve]['eighteen'], wavelength)
    FIR = extinction_fit(extcurves[curve]['FIR'], wavelength)
    return BKG + FUV + two + nine + eighteen + FIR

def extinction_fit(pars, wl):
    #generates individual component of the extinction curve (ie FUV, BKG) given wavelength range and parameters
    wl_lambda = wl / pars[1]
    lambda_wl = pars[1] / wl
    extinction = pars[0] / ((wl_lambda**pars[2])+(lambda_wl**pars[2])+pars[3])
    return extinction

def dust_tau(wl_rest,Rv,curve):
    #return tau (optical depth)
    Ev = extinction_curve(curve, 0.55)
    extcurve = extinction_curve(curve, wl_rest)
    return (Rv*extcurve)/(1.086*Ev)

# X-ray absorption functions
def compute_xraysig(xray_keV):
    #create a loop to obtain sigma for each value of xray energy
    dE = 0.15
    xray_sigma = np.zeros(len(xray_keV))
    for i,E in enumerate(xray_keV):
        if E < 10 and E > 0.1:
            Eindex = np.where((WilmskeV<E+dE) & (WilmskeV>E))
            sigma = totsigma[Eindex[0]]
            xray_sigma[i] = sigma[0]
        else:
            xray_sigma[i] = 0
    return xray_sigma

#xray attenuation function
def xray_attenuation(I0, sigma, NH):
    return I0 * np.exp(-sigma * NH)

def xray_atten(NH, xray_wl_rest, xray_flux):
    xray_keV_rest = (h*c)/(xray_wl_rest/1e6) #convert units
    xray_sigma = compute_xraysig(xray_keV_rest)
    attenuated_xray = xray_attenuation(xray_flux, xray_sigma, NH)   #determine attenuation
    return attenuated_xray

#define extinction curve specific attenuation functions
def dust_atten(wavelength_rest, I0, EBV, Rv, extcurve): #rest wavelength, power law, EBV, Rv, usercurve
    tau_opt = dust_tau(wavelength_rest, Rv, extcurve) #obtain the optical depth
    Iatt = I0 * np.exp(-1*EBV*tau_opt) #attenuate the power law
    Iatt[(1/wavelength_rest)>10.9]=0.0 #drop values where curve has blueshifted outside validity of Pei's model
    return Iatt

#find tau in optical region and xray region
def gettau(wl_rest, Rv, NHAv, curve):
    tau_opt = dust_tau(wl_rest[wl_rest>0.06],Rv,curve)
    sigma = compute_xraysig((h*c)/(wl_rest[wl_rest<0.06]/1e6))
    tau_xray = sigma*NHAv*Rv
    tau = np.concatenate((tau_opt,tau_xray))
    return tau

#fitting for MW type extinction curve
def MW_fitting(x, norm, beta, Ebv, redshift):
    Rv=3.08
    NHAv = 1.8e22
    wl_rest = x/(1+redshift)
    tau = gettau(wl_rest,Rv,NHAv,'MW')
    # Set flux to zero where wl_rest<912.5A
    x[(wl_rest>0.06)*(wl_rest<0.09)] = 0.0
    return (norm * x**beta) * np.exp(-1 * Ebv*tau)

#fitting for LMC type extinction curve
def LMC_fitting(x, norm, beta, Ebv, redshift):
    Rv = 3.16
    NHAv = 5.9e22
    wl_rest = x/(1+redshift)
    tau = gettau(wl_rest,Rv,NHAv,'LMC')
    # Set flux to zero where wl_rest<912.5A
    x[(wl_rest>0.06)*(wl_rest<0.09)] = 0.0
    return (norm * x**beta) * np.exp(-1 * Ebv*tau)

#fitting for SMC type extinction curve
def SMC_fitting(x, norm, beta, Ebv, redshift):
    Rv = 2.93
    NHAv = 2.8e22
    wl_rest = x/(1+redshift)
    tau = gettau(wl_rest,Rv,NHAv,'SMC')
    # Set flux to zero where wl_rest<912.5A
    x[(wl_rest>0.06)*(wl_rest<0.09)] = 0.0
    return (norm * x**beta) * np.exp(-1 * Ebv*tau)

#model setup using lmfit library
def model_setup(function,norm,beta,Ebv,redshift):
    atten_model = Model(function)
    params = atten_model.make_params(norm=norm,beta=beta,Ebv=Ebv,redshift=redshift)
    params['redshift'].vary=0
    params['beta'].vary=0
    return atten_model,params

def read_param(result,param):
    par = result.params[param].value
    if result.params[param].stderr:
        err = result.params[param].stderr
    else:
        err = 0.0        
    return par,err

#simulates noise
def addnoise(intensity):
    sigma = 3*np.sqrt(intensity)
    noise = np.random.normal(0, sigma)
    noisy_intensity = intensity + noise
    return noisy_intensity, noise

#calculates reduced chi squared
def chisquared(fitted_val, exp_intensity, err):
    each_chisq = (exp_intensity - fitted_val)**2 / err**2
    chisqs = np.nansum(each_chisq)
    dof = len(each_chisq) - 3
    red_chisq = chisqs/dof
    return red_chisq


#arrays with relevant wavelengths for each region
optical_wavelength = np.array([2.2,1.08,1.03,0.9,0.8,0.65,0.55,0.48,0.45,0.35,0.28,0.24,0.2]) 
xray_keV = np.arange(0.3, 10, 0.4)
xray_wavelength = ((h*c)/(xray_keV/1e6))

#takes userinput on curve parameters Rv, beta, EBV and the extinction curve (MW, LMC or SMC)
userbeta = float(input('What is beta? '))
userhostext = float(input('What is the host galaxy extinction? '))
userinterext = float(input('What is the intervening galaxy extinction? '))
hostz = float(input('What is the host galaxy redshift? '))
intz= float(input('What is the intervening galaxy redshift? '))
hostcurve = str(input('What type of extinction curve does the host galaxy have? '))
intercurve = str(input('What type of extinction curve does the intervening galaxy have? '))

#set up normalisation for determining noise
usernorm = 1e5

#sets up totsigma 
totsigma, WilmskeV = evaluate_sigma()

#simulate observed power law
flux_opt = usernorm*optical_wavelength**userbeta
flux_xray = usernorm*xray_wavelength**userbeta

#NH and Rv values correspond with each curve type
if hostcurve == 'MW':
    hostRv = 3.08
    hostNHAv = 1.8e22
    function = MW_fitting
elif hostcurve == 'LMC':
    hostRv = 3.16
    hostNHAv = 5.9e22
    function = LMC_fitting
elif hostcurve == 'SMC':
    hostRv = 2.93
    hostNHAv = 2.8e22
    function = SMC_fitting

#AG attenuated by host ISM
flux_hostext_opt = dust_atten(optical_wavelength/(1+hostz), flux_opt, userhostext, hostRv, hostcurve)
flux_hostext_xray = xray_atten(hostNHAv*hostRv*userhostext, xray_wavelength/(1+hostz), flux_xray)

#NH and Rv values correspond with each curve type
if intercurve == 'MW':
    intRv = 3.08
    intNHAv = 1.8e22
elif intercurve == 'LMC':
    intRv = 3.16
    intNHAv = 5.9e22
elif intercurve == 'SMC':
    intRv = 2.93
    intNHAv = 2.8e22

# AG attenuated by intervening galaxy ISM
flux_dblext_opt = dust_atten(optical_wavelength/(1+intz), flux_hostext_opt, userinterext, intRv, intercurve)
flux_dblext_xray = xray_atten(intNHAv*intRv*userinterext, xray_wavelength/(1+hostz), flux_hostext_xray)

#combine optical and xray region
wl_rest = np.concatenate((optical_wavelength/(1+hostz), xray_wavelength/(1+hostz)))
wl_obs = np.concatenate((optical_wavelength, xray_wavelength))
keV_obs = np.concatenate(((h*c)/(optical_wavelength/1e6),xray_keV))
flux_atten = np.concatenate((flux_hostext_opt, flux_hostext_xray))
flux_dbl_atten = np.concatenate((flux_dblext_opt, flux_dblext_xray))

#add noise to attenuation
noisy_flux_dbl_atten, dbl_atten_noise = addnoise(flux_dbl_atten)
noisy_flux_atten, atten_noise = addnoise(flux_atten)

#plot simulated curves for comparison
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1)

ax.errorbar(keV_obs,           
              noisy_flux_atten, 
              yerr = 3*np.sqrt(flux_atten),
              marker = '.',
              linestyle = ':',
              markersize = 6,
              color = 'blue',
              label = ('Host attenuation'),
              zorder = 4
              )

ax.errorbar(keV_obs,
            noisy_flux_dbl_atten,
            yerr = 3*np.sqrt(flux_dbl_atten),
            marker = '.',
            linestyle = '--',
            markersize = 6,
            color = 'black',
            label = ('Double attenuation'),
            zorder = 4
            )

ax.set_xlabel('Energy [keV]', fontsize = 20)
ax.set_ylabel('Attenuated Flux', fontsize = 20)
ax.set_xscale("log", base=10)
ax.set_yscale("log", base=10)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.tick_params(axis='both', which='minor', labelsize=15)
#ax.set_xlim(10e-20, 10e-15)
ax.set_ylim(-1000, 1e5)

#attempts to fit simulated data to a singly attenuated curve and prints & plots the result
try:
    hostmod, hostpars = model_setup(function, norm=usernorm, beta=userbeta,Ebv=userhostext,redshift=hostz)
    result = hostmod.fit(flux_dbl_atten, x=wl_obs, norm=usernorm, beta=userbeta, Ebv=userhostext, redshift=hostz)
    fit_beta, fit_beta_err = read_param(result,'beta')
    fit_ext, fit_ext_err = read_param(result,'Ebv')
    fit_norm, fit_norm_err = read_param(result,'norm')
    opt_fitted_curve = dust_atten(optical_wavelength/(1+hostz), fit_norm*optical_wavelength**fit_beta, fit_ext, hostRv, hostcurve)
    xray_fitted_curve = xray_atten(hostNHAv*hostRv*fit_ext, xray_wavelength/(1+hostz), fit_norm*xray_wavelength**fit_beta)
    fitted_curve = np.concatenate((opt_fitted_curve, xray_fitted_curve))
    chisq = chisquared(fitted_curve, noisy_flux_dbl_atten, 3*np.sqrt(flux_dbl_atten))
    
    print("Best-fit beta = %.1f +/- %.1e" % (fit_beta, fit_beta_err))
    print("Best-fit E(B-V) = %.1f +/- %.1e" % (fit_ext, fit_ext_err))
    print("Chi Squared = ", chisq)
    
    ax.errorbar(keV_obs,           
    fitted_curve, 
    marker = '.',
    linestyle = '-',
    markersize = 6,
    color = 'red',
    label = 'Fitted curve',
    zorder = 4
    )
 
#ValueError is thrown when lmfit fails so exception is included to highlight when this happens
except ValueError:
    print('fit failed')
    
ax.legend()