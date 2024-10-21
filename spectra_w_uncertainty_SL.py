import streamlit as st
import numpy as np
import pandas as pd
#import csv
#import time
import io
import matplotlib.pyplot as plt
from scipy.stats import norm, skew
from voigtfwhm import voigtfwhm

st.title("Absorbance spectra with uncertainty")
# Sidebar to take user inputs
st.sidebar.header("Simulation Controls")
k_logo = "images/kaust_2.png"
#a_logo = "images/aramco_logo.png"
#k_icon = "images/kaust_2.png"
st.logo(k_logo,icon_image=k_logo,size='large')
plt.style.use('dark_background')

wn_validation_flag = 1
wn_change_flag = 1
def wn_validation():
    if (wnstart >= wnend) or ((wnend - wnstart) > 10):
        st.error('Wavelength upper limit should be (up to 10 cm-1) larger than lower limit', icon="🚨")
        wn_validation_flag = 0
        wn_change_flag = 0
    else:
        wn_validation_flag = 1
        wn_change_flag = 1

    return wn_validation_flag, wn_change_flag

species_options = ['CH4 - HITRAN', 'H2O - HITRAN', 'CO2 - HITRAN', 'CO - HITRAN']

with st.sidebar:
    selected_species = st.selectbox("Species", species_options, 0)
    temperature = st.number_input("Temperature (K)", min_value=300, max_value=1000, value=300)
    pressure = st.number_input("Pressure (bar)", min_value=0.01, max_value=10.00, value=1.00)
    molefraction = st.number_input("Mole Fraction", min_value=0.00, max_value=1.00, value=0.01)
    pathlength = st.number_input('Pathlength (cm)', min_value=1, max_value=1000, step=1, value=10)
    N_simulations = st.number_input('Number of simulations', min_value=1, max_value=2000, step=1, value=1000)
    wnstart = st.number_input('Wavelength start (cm-1)', min_value=500.00, max_value=5000.00, step=0.01, value=1331.00)
    wnend = st.number_input('Wavelength end (cm-1)', min_value=500.00, max_value=5000.00, step=1.00, value=1334.00)
    s0_min_input = st.number_input("Line strength threshold (cm-1/(molec.cm-2))", min_value=1E-23, max_value=1E-19, value=1E-21, format="%f")


# Constants
h = 6.626070E-34  # Planck's constant (J.s)
kb = 1.380649E-23  # Boltzmann constant (J/K)
c = 2.99792458E+10  # Speed of light (cm/s)
R = 1.36259479E-22  # Gas constant (cm^3.atm.K^-1.molecule^-1)

# %% Initial parameters (replace MATLAB's "clear")
n_simulations = N_simulations
T = temperature  # Temperature in K
M = 16  # Molar mass of CH4 (g/mol)
mole_fraction = molefraction
P = pressure  # Pressure in atm
L = pathlength  # Path length in cm
s0_min = s0_min_input # minimum line strength

# Range
x = np.arange(wnstart, wnend, 0.001)  # Similar to MATLAB's 1331:0.001:1334

def import_data():
    # Load data (replace readmatrix and readtable)
    CH4lines = pd.read_csv('CH4_lines_formatted.csv').values
    tips = pd.read_csv('q32.csv', delim_whitespace=True).values
    return CH4lines, tips

# Define interpolation function for tips data (equivalent to tips1 = @(z) in MATLAB)
def tips1(z):
    return np.interp(z, tips[:, 0], tips[:, 1])

np.set_printoptions(legacy='1.25')

def find_range():
    #t = time.time()
    start_x = 0
    for i in range(0, len(CH4lines), 100):
        if (CH4lines[i, 0] < min(x)):
            start_x = start_x + 100
        else:
            break
    
    start_x = start_x - 100
    for i in range((start_x), len(CH4lines), 1):
        if (CH4lines[i, 0] < min(x)):
            start_x = start_x + 1
        else:
            break

    #print(time.time() - t)
    #t = time.time()

    end_x = start_x + 1
    for i in range(start_x,len(CH4lines)):
        if (CH4lines[i, 0] < max(x)):
            end_x = end_x + 1
        else:
            break
    #print(time.time() - t)
    return start_x, end_x

def extract_lines():
    # %% Get parameters and their uncertainties
    lines = []
    j = 0
    for i in range(start_x,end_x): #range(len(CH4lines)):
        if CH4lines[i,1] > s0_min: #(CH4lines[i, 0] > min(x)) and (CH4lines[i, 0] < max(x)):
            j += 1
            line = [
                CH4lines[i, 0],  # line position
                CH4lines[i, 1],  # line strength
                CH4lines[i, 2],  # gamma air
                CH4lines[i, 3],  # gamma self
                CH4lines[i, 4],  # LES
                CH4lines[i, 5],  # n_air
                CH4lines[i, 6]   # delta_air
            ]

            # Uncertainty handling (equivalent to switch cases in MATLAB)
            for k in [1]:
                uncertainty = CH4lines[i, 6 + k]
                if uncertainty == 0 or uncertainty == 1:
                    line.append(1)
                elif uncertainty == 2:
                    line.append(1E-1)
                elif uncertainty == 3:
                    line.append(1E-2)
                elif uncertainty == 4:
                    line.append(1E-3)
                elif uncertainty == 5:
                    line.append(1E-4)
                elif uncertainty == 6:
                    line.append(1E-5)
                elif uncertainty == 7:
                    line.append(1E-6)

            for k in [2, 3, 4, 5]:
                uncertainty = CH4lines[i, 6 + k]
                if uncertainty in [0, 1, 2]:
                    line.append(1)
                elif uncertainty == 3:
                    line.append(2.0E-1)
                elif uncertainty == 4:
                    line.append(2.0E-1)
                elif uncertainty == 5:
                    line.append(1.0E-1)
                elif uncertainty == 6:
                    line.append(5.0E-2)
                elif uncertainty == 7:
                    line.append(2.0E-2)
                elif uncertainty == 8:
                    line.append(1.0E-2)

            for k in [6]:
                uncertainty = CH4lines[i, 6 + k]
                if uncertainty == 0 or uncertainty == 1:
                    line.append(1)
                elif uncertainty == 2:
                    line.append(1E-1)
                elif uncertainty == 3:
                    line.append(1E-2)
                elif uncertainty == 4:
                    line.append(1E-3)
                elif uncertainty == 5:
                    line.append(1E-4)
                elif uncertainty == 6:
                    line.append(1E-5)
                elif uncertainty == 7:
                    line.append(1E-6)

            lines.append(line)
    #st.text('Number of simulated lines: '+str(len(lines)))
    return lines

# %% Simulate spectra with uncertainty
def rand_distribution(mu, sigma):
    """Equivalent to rand_distrition in MATLAB."""
    range_vals = np.linspace(mu - 3 * sigma, mu + 3 * sigma, 100)
    cdf_vals = norm.cdf(range_vals, loc=mu, scale=sigma)
    return cdf_vals, range_vals

def random_value(cdf, range_vals):
    """Equivalent to random_value in MATLAB."""
    return np.interp(np.random.rand(), cdf, range_vals)

def extract_parameters():
    x0_rand_cdf = np.zeros((len(lines),100))
    x0_rand_range = np.zeros((len(lines),100))
    S0_rand_cdf = np.zeros((len(lines),100))
    S0_rand_range = np.zeros((len(lines),100))
    gamma_air_rand_cdf = np.zeros((len(lines),100))
    gamma_air_rand_range = np.zeros((len(lines),100))
    gamma_self_rand_cdf = np.zeros((len(lines),100))
    gamma_self_rand_range = np.zeros((len(lines),100))
    n_air_rand_cdf = np.zeros((len(lines),100))
    n_air_rand_range = np.zeros((len(lines),100))
    delta_air_rand_cdf = np.zeros((len(lines),100))
    delta_air_rand_range = np.zeros((len(lines),100))

    x0 = np.zeros(len(lines))
    s0 = np.zeros(len(lines))
    gamma_air_0 = np.zeros(len(lines))
    gamma_self_0 = np.zeros(len(lines))
    n_air = np.zeros(len(lines))
    delta_air = np.zeros(len(lines))


    #fig, ax = plt.subplots()

    j=0
    for line in lines:
        x0_rand_cdf[j], x0_rand_range[j] = rand_distribution(line[0], line[7] / 3)
        S0_rand_cdf[j], S0_rand_range[j] = rand_distribution(line[1], line[1] * line[8] / 3)
        gamma_air_rand_cdf[j], gamma_air_rand_range[j] = rand_distribution(line[2], line[2] * line[9] / 3)
        gamma_self_rand_cdf[j], gamma_self_rand_range[j] = rand_distribution(line[3], line[3] * line[10] / 3)
        n_air_rand_cdf[j], n_air_rand_range[j] = rand_distribution(line[5], line[5] * line[11] / 3)
        delta_air_rand_cdf[j], delta_air_rand_range[j] = rand_distribution(line[6], line[12] / 3)

        x0[j] = line[0]
        s0[j] = line[1]
        gamma_air_0[j] = line[2]
        gamma_self_0[j] = line[3]
        n_air[j] = line[5]
        delta_air[j] = line[6]

        j = j + 1

    return  x0_rand_cdf, x0_rand_range,S0_rand_cdf, S0_rand_range,gamma_air_rand_cdf, gamma_air_rand_range\
    ,gamma_self_rand_cdf, gamma_self_rand_range,n_air_rand_cdf, n_air_rand_range,\
    delta_air_rand_cdf, delta_air_rand_range, x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air

def MC_simulation():
    # Run the simulations
    for i in range(n_simulations):
        j=0
        for line in lines:
            # Line position

            x0_rand = random_value(x0_rand_cdf[j], x0_rand_range[j])

            delta_air_rand = random_value(delta_air_rand_cdf[j], delta_air_rand_range[j])
            x0_shifted_rand = x0_rand + P * (1 - mole_fraction) * delta_air_rand

            # Line strength

            S0_rand = random_value(S0_rand_cdf[j], S0_rand_range[j])
            S_rand = S0_rand * (tips1(296) / tips1(T)) * np.exp(-(h * c * line[4] / kb) * (1 / T - 1 / 296)) \
                     * (1 - np.exp(-h * c * x0_rand / (kb * T))) / (1 - np.exp(-h * c * x0_rand / (kb * 296)))

            A_rand = S_rand * L * mole_fraction * (P / (R * T))  # cm-1

            # Doppler broadening
            wG_rand = x0_shifted_rand * (7.1623E-7) * np.sqrt(T / M)

            # Pressure broadening
            n_air_rand = random_value(n_air_rand_cdf[j], n_air_rand_range[j])


            gamma_self_rand = random_value(gamma_self_rand_cdf[j], gamma_self_rand_range[j])

            gamma_air_rand = random_value(gamma_air_rand_cdf[j], gamma_air_rand_range[j])

            gamma_self_rand = gamma_self_rand * (296 / T) ** n_air_rand
            gamma_air_rand = gamma_air_rand * (296 / T) ** n_air_rand

            wL_rand = P * (mole_fraction * 2 * gamma_self_rand + (1 - mole_fraction) * 2 * gamma_air_rand)


            #spectra[:, i] += (x, [A_rand, x0_shifted_rand, wG_rand, wL_rand])

            spectra[:, i] += voigtfwhm(x, [A_rand, x0_shifted_rand, wG_rand, wL_rand])

            j=j+1

        if np.isnan(np.mean(spectra[:, i])):
            spectra[:, i] = spectra[:, i - 1]
        else:
            ax.plot(x, spectra[:, i], '.', markersize=0.1, color='red')
    return spectra

def mean_spectrum_simulation():
    spectrum_mean_parameters = np.zeros(len(x))
    j = 0
    for line in lines:
        # Line position

        x0_shifted = x0[j] + P * (1 - mole_fraction) * delta_air[j]
        # Line strength

        S = s0[j] * (tips1(296) / tips1(T)) * np.exp(-(h * c * line[4] / kb) * (1 / T - 1 / 296)) \
                 * (1 - np.exp(-h * c * x0[j] / (kb * T))) / (1 - np.exp(-h * c * x0[j] / (kb * 296)))
        A = S * L * mole_fraction * (P / (R * T))  # cm-1
        # Doppler broadening
        wG = x0_shifted * (7.1623E-7) * np.sqrt(T / M)
        # Pressure broadening

        gamma_self = gamma_self_0[j] * (296 / T) ** n_air[j]
        gamma_air = gamma_air_0[j] * (296 / T) ** n_air[j]
        wL = P * (mole_fraction * 2 * gamma_self + (1 - mole_fraction) * 2 * gamma_air)

        #spectra[:, i] += (x, [A_rand, x0_shifted_rand, wG_rand, wL_rand])
        spectrum_mean_parameters +=  np.transpose(voigtfwhm(x, [A, x0_shifted, wG, wL]))

        j=j+1
    return spectrum_mean_parameters

def calc_error_bars():
    error_bars = np.zeros(len(x))
    relative_uncertainty = np.zeros(len(x))
    skewness = np.zeros(len(x))
    std_residuals = np.zeros(N_simulations)
    for i in range(len(x)):
        relative_uncertainty[i] = 100*3*np.std(spectra[i][:])/np.mean(spectrum_mean_parameters[i])
        skewness[i] = skew(spectra[i][:])
        error_bars[i] = 3*np.std(spectra[i][:])
    
    max_std_index = np.argmax(error_bars)
    print(max_std_index)
    
    for i in range(2,n_simulations):
        std_residuals[i] = np.std(spectra[max_std_index][range(i)])#/np.std(spectra[i][range(n_simulations)])
    return relative_uncertainty, error_bars, skewness, std_residuals

def plotting_commands():
    # Plot the mean spectrum
    #mean_spectrum = np.mean(spectra, axis=1)
    #ax.plot(x, mean_spectrum, '--', linewidth=2, color='blue')
    #ax.plot(x, spectrum_mean_parameters, '-', linewidth=2, color='black')
    #st.pyplot(fig)
    fig, ax = plt.subplots()
    ax.set_title('Absorbance based on mean parameters and standard deviation bars')
    ax.set_xlabel('Wavenumbers (cm-1)')
    ax.set_ylabel('Absorbance')
    ax.errorbar(x, spectrum_mean_parameters, yerr=error_bars, fmt='-', color='black', ecolor='red')
    #ax.plot(x, error_bars, '-', linewidth=2, color='black')
    st.pyplot(fig)

def plot_uncertainty():
    fig, ax1 = plt.subplots()

    color = (1,1,1,1)
    ax1.set_xlabel('Wavenumbers (cm-1)')
    ax1.set_ylabel('Relative uncertainty (%)', color=color)
    ax1.plot(x, relative_uncertainty, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    if max(relative_uncertainty) < 100:
        ax1.set_ylim(0,100)

    ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

    color = (1,1,0,0.8)
    ax2.set_ylabel('Skewness', color=color)  # we already handled the x-label with ax1
    ax2.plot(x, skewness, color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    if abs(max(skewness)) < 1:
        ax2.set_ylim(-1,1)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    st.text('Uncertainty spectrum (3 x std.) and skewness spectrum:')
    st.pyplot(fig)

    fig, ax = plt.subplots()
    #print(std_residuals)
    ax.plot(range(n_simulations), std_residuals)
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Standard deviation')
    max_std_index = np.argmax(error_bars)
    st.text('Standard deviation with iterations at ('+str(round(x[max_std_index],2))+' cm-1):')
    st.pyplot(fig)


wn_validation_flag, wn_change_flag = wn_validation()
print(wn_validation_flag)
if wn_validation_flag == 1:
    
    CH4lines, tips = import_data()   
    start_x, end_x = find_range()
    lines = extract_lines()
    number_of_lines = len(lines)
    x0_rand_cdf, x0_rand_range,S0_rand_cdf, S0_rand_range,gamma_air_rand_cdf, gamma_air_rand_range\
    ,gamma_self_rand_cdf, gamma_self_rand_range,n_air_rand_cdf, n_air_rand_range,\
    delta_air_rand_cdf, delta_air_rand_range, x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air = extract_parameters()
    spectra = np.zeros((len(x), n_simulations))
    
    st.text('Absorbance based on mean line-parameters of '+str(len(lines)) + ' lines,\nand '+str(n_simulations)+' spectra based on randomly sampled line-parameters:')
    fig, ax = plt.subplots()
    #ax.set_title('Absorbance based on mean parameters and '+str(n_simulations)+' simulated spectra')
    ax.set_xlabel('Wavenumbers (cm-1)')
    ax.set_ylabel('Absorbance')
    
    spectra = MC_simulation()
    spectrum_mean_parameters = mean_spectrum_simulation()

    ax.plot(x, spectrum_mean_parameters, '-', color='white')
    st.pyplot(fig)

    relative_uncertainty, error_bars, skewness, std_residuals = calc_error_bars()
    #plotting_commands()
    
    plot_uncertainty()    

    arr = np.asarray(np.transpose([x,spectrum_mean_parameters,error_bars,relative_uncertainty,skewness]))
    # Create an in-memory buffer
    with io.BytesIO() as buffer:
        # Write array to buffer
        np.savetxt(buffer, arr, delimiter=",")
        st.download_button(
            label="Download result as CSV",
            data = buffer, # Download buffer
            file_name = 'results.csv',
            mime='text/csv'
        ) 

