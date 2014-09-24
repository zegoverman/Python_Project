'''Module to Find the Best-Fit Model for the Geometry of the Universe using Measurements of Supernovae'''

from __future__ import division
import numpy
import scipy
import matplotlib.pyplot as pyplot
from scipy import integrate

USER = 'Elliot Byrne'
USER_ID = 'jtdl85'

def data_from_file(filename):
    '''Function to obtain the Supernovae data'''
    data = numpy.loadtxt(filename,usecols=(1, 2, 3) )
    return data

def observed_flux(magnitude):
    '''Function to convert observed magnitudes to fluxes'''
    observed_flux = []
    for i in magnitude:
        f = 10**((m0-i)/2.5) # converted magnitudes to fluxes
        observed_flux.append(f)
    flux =  numpy.array(observed_flux) #converted fluxes to an array
    return flux

def error_flux(mag_error, observed_magnitude):
    '''Function to find the error on flux'''
    observed_error = []
    for i in range(0, len(mag_error)):
        divided = (mag_error[i])/2.5
        flux_error = 10**((m0-observed_magnitude[i])/2.5)*numpy.log(10)*divided
        observed_error.append(flux_error)
    observed_error = numpy.array(observed_error)
    return observed_error 

def find_Sk_rc(redshift):
    '''Function to calculate Sk_rc'''
    Sk_rc = (c*redshift)/H0
    return Sk_rc

def find_Lpeak(redshifts, obs_flux):
    '''Function to calculate Lpeak'''
    Lpeak = []
    for i in range(0, len(redshifts)):
        Sk_rc = find_Sk_rc(redshifts[i])
        LP =  4*numpy.pi*(Sk_rc**2)*((1 + redshifts[i])**2)*obs_flux[i] 
        Lpeak.append(LP)
    Lpeak = numpy.array(Lpeak)
    return Lpeak

def find_Lpeak_trials(max, min):
    '''Function to calculate the trial values of Lpeak to be used'''
    trial_values = []
    chosen_range = 50 # arbitary spacing
    spacing = ((max_Lpeak - min_Lpeak)/chosen_range)
    for i in range(0, chosen_range+1):
        number = min_Lpeak + (i*spacing)
        trial_values.append(number)
    trial_values = numpy.array(trial_values)
    return trial_values

def find_model_flux(Lpeak, redshift):
    '''Function to calculate the model flux using Lpeak'''
    Sk_rc = find_Sk_rc(redshift)
    flux = Lpeak/((4*numpy.pi*Sk_rc**2)*(1 + redshift)**2)
    return flux

def find_chi_squared(fobserved, fmodel, error):
    '''Function to calculate chi squared'''
    chi = (fobserved - fmodel)**2 / (error**2)
    return chi

def choose_best_fit_value(list_of_Lpeaks, chi_squared_list):
    '''Function to choose the value with the lowest chi squared'''
    best_fit = chi_squared_list.min()
    index_lowest_chi = numpy.where(chi_squared_list == best_fit)
    value = list_of_Lpeaks[index_lowest_chi]
    return value

def fluxes_of_a_variable(big_list, list_of_Omegas_or_Lpeaks ,value_to_match):
    '''Function to choose the fluxes associated with the given OmegaLambda or Lpeak'''
    index = numpy.where(list_of_Omegas_or_Lpeaks == value_to_match)
    fluxes = big_list[index]
    return fluxes

def find_OmegaLambda_trial_values(trial_number):
    '''Function to find the trial values for Omega_Lambda'''
    OmegaLambda_trial_spacing = 1/number_of_trials
    OmegaLambda_trials = []
    for index in range(0, number_of_trials+1):
        value = index*OmegaLambda_trial_spacing
        OmegaLambda_trials.append(value)
    OmegaLambda_trials = numpy.array(OmegaLambda_trials)
    return OmegaLambda_trials

def find_H(OmegaLambda, OmegaM0, redshifts):
    '''Function to calculate H'''
    a = 1 / (1 + redshifts)
    H = (H0**2*OmegaM0*(1/a)**3 + H0**2*OmegaLambda)**0.5
    return H

def integral(a, H):
    '''Function to integrate in order to find rc'''
    integrand = lambda a: (c/(H*a**2))
    rH0 = (integrate.quad(integrand, a, 1))
    return rH0

def find_model_flux_using_rc(rc_value, redshift, Lpeak):
    '''Function to calculate flux using rc, redshift and Lpeak'''
    flux = Lpeak / (4*numpy.pi*rc_value**2*(1+redshift)**2)
    #print type(flux) # now a number
    return flux

def find_fluxes(OmegaLambda, OmegaM0, redshifts, Lpeak):
    '''Function to calculate the model fluxes values given Omega_Lambdas, redshifts and Lpeak'''
    flux_values = []
    for index in range(0, len(OmegaLambda)):
        small_list = []
        for p in range(0, len(redshifts)):
            a = 1 / (1 + redshifts[p]) 
            H = find_H(OmegaLambda[index], OmegaM0[index], redshifts[p]) # now have H ,find rc for each 
            #print H #numbers 
            rc = integral(a, H) # value for rc 
            rc = rc[0] 
            mod_flux = find_model_flux_using_rc(rc, redshifts[p], Lpeak[0]) # flux values for all redshifts at different OmegaLambda
            small_list.append(mod_flux)
        flux_values.append(small_list)
    return flux_values

def fluxes_to_mags(fluxes):
    '''Function to convert an array of fluxes to magnitudes'''
    mags = []
    for element in range(0, len(fluxes)):
        mag = m0 - 2.5*numpy.log10(fluxes[element])
        mags.append(mag)
    return numpy.array(mags)

def plot_graph():
    '''Function to plot a graph of the observed and model magnitudes against redshift'''
    pyplot.figure(1)
    pyplot.title('magnitude against Redshift for Observed and Model Data')
    pyplot.xlabel('Redshift')
    pyplot.ylabel('magnitude')
    # plot the observed and model magnitudes against redshift for both models and observed data, with a legend
    pyplot.plot(high_red_shift, OmegaLambda_model_magnitudes, 'bo', label = ' Omega_Lambda Model Data')
    pyplot.plot(low_red_shift, Lpeak_model_magnitudes,'yo', label = 'Lpeak model data')
    pyplot.plot(high_red_shift, high_red_shift_mag,'ro', label = 'Observed Data high')
    pyplot.plot(low_red_shift, low_red_shift_mag,'go', label = 'Observed data low')
    pyplot.legend(loc='center right')
    pyplot.show()

data = data_from_file('sn_data.txt')
high_red_shift_mag = data[:42,1] 
low_red_shift_mag = data[42:,1]
size_data = len(data)
high_red_shift_mag_error = data[:42,2]
low_red_shift_mag_error = data[42:,2]
low_red_shift = data[42:,0]
high_red_shift = data[:42,0]
H0 = 75000 # ms-1/Mpc
c = 3*10**8 # ms-1
m0 = -20.45

#mag to flux
observed_fluxes_high = observed_flux(high_red_shift_mag)
observed_fluxes_low = observed_flux(low_red_shift_mag)

#error mag to error flux
mag_error = data[:,2]
observed_fluxes_high_error = error_flux(high_red_shift_mag_error, high_red_shift_mag)
observed_fluxes_low_error = error_flux(low_red_shift_mag_error, low_red_shift_mag)

#find a value for Lpeak using each low read shift data
Lpeak_low = find_Lpeak(low_red_shift,observed_fluxes_low) #correct

#choose max and min Lpeak and make trial values in this range
max_Lpeak = Lpeak_low.max()
min_Lpeak = Lpeak_low.min()
Lpeak_trials = find_Lpeak_trials(max_Lpeak, min_Lpeak) # Lpeak values I will try 

huge_list = [] 
for i in range (0, len(Lpeak_trials)):
    temp_list = []
    for j in range(0, len(low_red_shift)):
        model_flux_low_red_shift = find_model_flux(Lpeak_trials[i],low_red_shift[j])
        temp_list.append(model_flux_low_red_shift) 
    huge_list.append(temp_list)
huge_list = numpy.array(huge_list) 

#calc chi squared 
chi_for_Lpeak = []
for index in huge_list: #gets a list of fluxes at each Lpeak
    each_Lpeak = []
    for j in range(0, len(index)):
        d = find_chi_squared(observed_fluxes_low[j], index[j], observed_fluxes_low_error[j]) 
        each_Lpeak.append(d)
    summation_chi = sum(each_Lpeak)
    chi_for_Lpeak.append(summation_chi)
chi_for_Lpeak = numpy.array(chi_for_Lpeak)

# choose your Lpeak value
final_Lpeak = choose_best_fit_value(Lpeak_trials, chi_for_Lpeak)
print final_Lpeak, 'The value for Lpeak' # This is my final Lpeak value 

# now find OmegaLambda
number_of_trials = 50 # It's actually 51 as 0 and 50 are included
OmegaLambda_trial_values = find_OmegaLambda_trial_values(number_of_trials)
OmegaM0_trial_values = 1 - OmegaLambda_trial_values

model_fluxes_high_red_shift = find_fluxes(OmegaLambda_trial_values, OmegaM0_trial_values, high_red_shift, final_Lpeak) # now have values for flux at each OmegaLambda
model_fluxes_high_red_shift = numpy.array(model_fluxes_high_red_shift)

#find chi squared for each OmegaLambda
chi_for_OmegaLambda = []
for element in model_fluxes_high_red_shift:
    new_list = []
    for k in range(0, len(element)):
        dd = find_chi_squared(observed_fluxes_high[k], element[k], observed_fluxes_high_error[k])
        new_list.append(dd)
    sum_chi = sum(new_list)
    chi_for_OmegaLambda.append(sum_chi)
chi_for_OmegaLambda = numpy.array(chi_for_OmegaLambda)

#choose my OmegaLambda
final_OmegaLambda = choose_best_fit_value(OmegaLambda_trial_values, chi_for_OmegaLambda)
print final_OmegaLambda, 'The value for OmegaLambda'

# make lists of fluxes using the decided OmegaLambda and Lpeak
Lpeak_model_fluxes = fluxes_of_a_variable(huge_list, Lpeak_trials, final_Lpeak)[0]
Lpeak_model_magnitudes = fluxes_to_mags(Lpeak_model_fluxes)
OmegaLambda_model_fluxes = fluxes_of_a_variable(model_fluxes_high_red_shift, OmegaLambda_trial_values, final_OmegaLambda)[0]
OmegaLambda_model_magnitudes = fluxes_to_mags(OmegaLambda_model_fluxes) 

graph = plot_graph()
