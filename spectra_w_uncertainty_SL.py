import streamlit as st
import numpy as np
import pandas as pd
#import csv
#import time
import matplotlib.pyplot as plt
from scipy.stats import norm, skew
from voigtfwhm import voigtfwhm
import io
import base64
import uuid
import datetime
import csv

st.set_page_config(
    page_title="MonteSpectra",
    page_icon="🗻",
    layout="centered",
    initial_sidebar_state="expanded",
    menu_items={
        'Get Help': 'https://www.kaust.edu.sa',
        'Report a bug': "mailto:ihsan.farouki@kaust.edu.sa",
        'About': "# Absorbance spectra simulations with uncertainty quantification"
    }
)

# Sidebar to take user inputs
k_logo = "images/kaust_2.png"
MS_logo = "images/montespectra_logo.png"
MS_logo_2 = "images/montespectra_logo_2.png"
empty_image = "images/Empty.png"
st.title("Simulation Results")
#a_logo = "images/aramco_logo.png"
#k_icon = "images/kaust_2.png"
#st.logo(k_logo,icon_image=k_logo,size='large')
st.logo(empty_image,icon_image=MS_logo_2,size='large')
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

if 'dek' not in st.session_state:
    st.session_state.dek = str(uuid.uuid4())

def change_wn_range():
    if st.session_state.selected_species == '(12)CH4 - HITRAN':
        st.session_state.wn_start = 1331
        st.session_state.wn_end = 1334
    elif st.session_state.selected_species == 'H2(16)O - HITRAN':
        st.session_state.wn_start = 3742
        st.session_state.wn_end = 3747
    elif st.session_state.selected_species == '(12)CO2 - HITRAN':
        st.session_state.wn_start = 2300
        st.session_state.wn_end = 2305
    elif st.session_state.selected_species == '(14)N2O - HITRAN':
        st.session_state.wn_start = 1285
        st.session_state.wn_end = 1290

def molar_mass():
    if st.session_state.selected_species == '(12)CH4 - HITRAN':
        M = 16 # Molar mass of CH4 (g/mol)
    elif st.session_state.selected_species == 'H2(16)O - HITRAN':
        M = 18 # Molar mass of CH4 (g/mol)
    elif st.session_state.selected_species == '(12)CO2 - HITRAN':
        M = 44 # Molar mass of CH4 (g/mol)
    elif st.session_state.selected_species == '(14)N2O - HITRAN':
        M = 44 # Molar mass of CH4 (g/mol)
    
    return M
        
species_options = ['(12)CH4 - HITRAN', 'H2(16)O - HITRAN', '(12)CO2 - HITRAN', '(14)N2O - HITRAN']

with st.sidebar:
    st.image(MS_logo, width=250)
    st.divider() 
    st.header("Simulation Controls")
    selected_species = st.selectbox("Species", species_options, 0, on_change=change_wn_range,key='selected_species')
    temperature = st.number_input("Temperature (K)", min_value=300, max_value=1000, value=300, step=100)
    pressure = st.number_input("Pressure (bar)", min_value=0.01, max_value=10.00, value=1.00, step=0.2)
    molefraction = st.number_input("Mole Fraction", min_value=0.00, max_value=1.00, value=0.01, step=0.01)
    pathlength = st.number_input('Pathlength (cm)', min_value=1, max_value=1000, step=1, value=10)
    wnstart = st.number_input('Wavelength start (cm-1)', min_value=500.00, max_value=5000.00, step=0.01, value=1331.00, key='wn_start')
    wnend = st.number_input('Wavelength end (cm-1)', min_value=500.00, max_value=5000.00, step=1.00, value=1334.00, key='wn_end')
    wnres = st.number_input('Resolution (cm-1)', min_value=0.001, max_value=0.1, step=0.001, value=0.005, key='wn_res', format="%.3f")
    st.divider() 
    st.subheader('Advanced simulation controls')
    manual_control = st.toggle("Enable manual control of line parameters",key='manual_control')
    N_simulations = st.number_input('Number of simulations', min_value=1, max_value=2000, step=100, value=1000)
    s0_min_input = st.number_input("Line strength threshold (cm-1/(molec.cm-2))", min_value=1E-23, max_value=1E-19, value=1E-21, format="%.3e")
    manual_control = st.toggle("Plot standard deviation vs iterations to test convergence",key='conv_test')
    if st.session_state.conv_test:
        convergence_frequency = st.number_input('Convergence test position (cm-1)', min_value=st.session_state.wn_start, max_value=st.session_state.wn_end, value=st.session_state.wn_start, key='wn_conv')
    st.divider() 
    st.subheader('Warnings')

# Constants
h = 6.626070E-34  # Planck's constant (J.s)
kb = 1.380649E-23  # Boltzmann constant (J/K)
c = 2.99792458E+10  # Speed of light (cm/s)
R = 1.36259479E-22  # Gas constant (cm^3.atm.K^-1.molecule^-1)

# %% Initial parameters (replace MATLAB's "clear")
n_simulations = N_simulations
T = temperature  # Temperature in K
mole_fraction = molefraction
P = pressure  # Pressure in atm
L = pathlength  # Path length in cm
s0_min = s0_min_input # minimum line strength
M = molar_mass()  # Molar mass of selected species

# Range
x = np.arange(wnstart, wnend, wnres)  # Similar to MATLAB's 1331:0.001:1334

@st.cache_data
def import_data(selected_species):
    print('importing data')
    # Load data (replace readmatrix and readtable)
    if selected_species == '(12)CH4 - HITRAN':
        CH4lines = pd.read_csv('CH4_lines_formatted.csv').values
        tips = pd.read_csv('q32_12CH4.csv', sep='\s+').values
    elif selected_species == 'H2(16)O - HITRAN':
        CH4lines = pd.read_csv('H2O_lines_formatted.csv').values
        tips = pd.read_csv('q1_H2O16.csv', sep='\s+').values
    elif selected_species == '(12)CO2 - HITRAN':
        CH4lines = pd.read_csv('CO2_lines_formatted.csv').values
        tips = pd.read_csv('q7_12CO2.csv', sep='\s+').values
    elif selected_species == '(14)N2O - HITRAN':
        CH4lines = pd.read_csv('N2O_lines_formatted.csv').values
        tips = pd.read_csv('q21_14N2O.csv', sep='\s+').values
    return CH4lines, tips

# Define interpolation function for tips data (equivalent to tips1 = @(z) in MATLAB)
def tips1(z):
    return np.interp(z, tips[:, 0], tips[:, 1])

np.set_printoptions(legacy='1.25')

@st.cache_data
def find_range(x,CH4lines):
    print('finding start and end indices')
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

@st.cache_data
def extract_lines(start_x,end_x,CH4lines,s0_min):
    print('extracting lines')
    st.session_state.dek = str(uuid.uuid4())
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
                if uncertainty == 0:
                    with st.sidebar:
                        st.sidebar.warning('Uncertainty in parameter '+str(k)+' is unavailable for '+str(CH4lines[i, 0]), icon="⚠️")
                    line.append(0)
                elif uncertainty == 1:
                    line.append(1.0)
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
                    with st.sidebar:
                        st.sidebar.warning('Uncertainty in parameter '+str(k)+' is unavailable for '+str(np.round(CH4lines[i, 0],2)), icon="⚠️")
                    line.append((100)*0)
                elif uncertainty == 3:
                    line.append((100)*2.0E-1)
                elif uncertainty == 4:
                    line.append((100)*2.0E-1)
                elif uncertainty == 5:
                    line.append((100)*1.0E-1)
                elif uncertainty == 6:
                    line.append((100)*5.0E-2)
                elif uncertainty == 7:
                    line.append((100)*2.0E-2)
                elif uncertainty == 8:
                    line.append((100)*1.0E-2)

            for k in [6]:
                uncertainty = CH4lines[i, 6 + k]
                if uncertainty == 0:
                    with st.sidebar:
                        st.sidebar.warning('Uncertainty in parameter '+str(k)+' is unavailable for '+str(CH4lines[i, 0]), icon="⚠️")
                    line.append(0)
                elif uncertainty == 1:
                    line.append(1.0)
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

@st.cache_data
def extract_parameters(lines):
    print('extracting parameters and generating distributions')
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
        S0_rand_cdf[j], S0_rand_range[j] = rand_distribution(line[1], line[1] * (1/100)*line[8] / 3)
        gamma_air_rand_cdf[j], gamma_air_rand_range[j] = rand_distribution(line[2], line[2] * (1/100)*line[9] / 3)
        gamma_self_rand_cdf[j], gamma_self_rand_range[j] = rand_distribution(line[3], line[3] * (1/100)*line[10] / 3)
        n_air_rand_cdf[j], n_air_rand_range[j] = rand_distribution(line[5], line[5] * (1/100)*line[11] / 3)
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

@st.cache_data
def MC_simulation(lines,n_simulations,T,P,mole_fraction,L,x):
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
        #else:
            #ax.plot(x, spectra[:, i], '.', markersize=0.1, color="#A87BF9")
    return spectra

@st.cache_data
def mean_spectrum_simulation(lines,T,P,mole_fraction,L,x):
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

@st.cache_data
def calc_error_bars(spectra,spectrum_mean_parameters):
    error_bars = np.zeros(len(x))
    relative_uncertainty = np.zeros(len(x))
    skewness = np.zeros(len(x))
    for i in range(len(x)):
        relative_uncertainty[i] = 100*3*np.std(spectra[i][:])/np.mean(spectrum_mean_parameters[i])
        skewness[i] = skew(spectra[i][:])
        error_bars[i] = 3*np.std(spectra[i][:])
    
    return relative_uncertainty, error_bars, skewness
  
def std_deviation_with_iterations():
    std_residuals = np.zeros(N_simulations)
    #max_std_index = np.argmax(error_bars)
    std_index = np.where(np.round(x,3) == st.session_state.wn_conv)[0][0]
    #print(std_index)
    #print('at ('+str(round(x[std_index],2))+' cm-1')
    #print(max_std_index)
    for i in range(2,n_simulations):
        std_residuals[i] = np.std(spectra[std_index][range(i)])#/np.std(spectra[i][range(n_simulations)])
    
    fig, ax = plt.subplots()
    #print(std_residuals)
    ax.plot(range(n_simulations), std_residuals,color="#ECBC7A")
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Standard deviation')
    
    st.divider() 
    st.write('_Standard deviation with iterations at ('+str(round(x[std_index],2))+' cm-1):_')
    st.pyplot(fig)
    
    #return std_residuals

@st.cache_data
def plot_MC_spectra(spectra, spectrum_mean_parameters):
    
    ax.plot(x, spectrum_mean_parameters, '-', color='white',zorder=n_simulations+1)

    for i in range(n_simulations):
        ax.plot(x, spectra[:, i], '-', lw=1, color="#A87BF9")
    
    ax.legend(['Mean parameters', 'Unceratinty envelope'])

    textstr = (str(100*mole_fraction)+'% ' + selected_species + '\n' + str(T) + ' K\n' + str(P) + ' atm\n'+ str(L) + ' cm')
    props = dict(boxstyle='round', facecolor="#A87BF9", alpha=0.1)
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)

    secax = ax.secondary_xaxis('top', functions=(wn2wl, wl2wn))
    secax.set_xlabel('Wavelength (µm)')
    st.pyplot(fig)

def wn2wl(x):
    return 10000 / x

def wl2wn(x):
    return 10000 / x

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

@st.cache_data
def plot_uncertainty(relative_uncertainty,skewness):
    fig, ax1 = plt.subplots()

    color = (1,1,1,1)
    ax1.set_xlabel('Wavenumbers (cm-1)')
    ax1.set_ylabel('Relative uncertainty (%)', color=color)
    ax1.plot(x, relative_uncertainty, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    if max(relative_uncertainty) < 100:
        ax1.set_ylim(0,100)

    secax = ax1.secondary_xaxis('top', functions=(wn2wl, wl2wn))
    secax.set_xlabel('Wavelength (µm)')

    ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

    color = ("#57D2E9")
    ax2.set_ylabel('Skewness', color=color)  # we already handled the x-label with ax1
    ax2.plot(x, skewness, color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    if abs(max(skewness)) < 1:
        ax2.set_ylim(-1,1)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    st.divider() 
    st.write('_Uncertainty spectrum (3 x std.) and skewness spectrum:_')
    
    st.pyplot(fig)

wn_validation_flag, wn_change_flag = wn_validation()
#print(wn_validation_flag)
if wn_validation_flag == 1:
    
    CH4lines, tips = import_data(selected_species)   
    start_x, end_x = find_range(x,CH4lines)
    lines = extract_lines(start_x,end_x,CH4lines,s0_min)
    number_of_lines = len(lines)
    
    if st.session_state.manual_control:
        #print(lines)
        #print('data editor key: '+str(st.session_state.dek))
        st.divider() 
        st.write('_Editable list of parameters for lines within the selected range:_')
        lines_column_config = {
        1: st.column_config.NumberColumn(
            "v0",
            help="Line position (cm-1)",
            min_value=500,
            max_value=5000,
            step=0.0001,
            format="%.4f",
        ),
        2: st.column_config.NumberColumn(
            "S0",
            help="Line strength at 296 K (cm-1/molec.cm-2)",
            min_value=1.000e-25,
            max_value=1.000e+1,
            #step=1.000e-25,
            format="%.4e",
        ),
        3: st.column_config.NumberColumn(
            "γ_bath",
            help="Bath gas broadning coeffecient (cm-1/atm) at 296 K. Value imported from HITRAN is for air.",
            min_value=0,
            #max_value=10,
            step=0.0001,
            format="%.5f",
        ),
        4: st.column_config.NumberColumn(
            "γ_self",
            help="Self gas broadning coeffecient (cm-1/atm) at 296 K.",
            min_value=0,
            #max_value=10,
            step=0.0001,
            format="%.5f",
        ),
        5: st.column_config.NumberColumn(
            "E''",
            help="Lower state energy (cm-1)",
            min_value=0,
            #max_value=end_x,
            step=0.001,
            format="%.3f",
        ),
        6: st.column_config.NumberColumn(
            "n",
            help="Temperature exponent for broadning coeffecients. Value imported from HITRAN is for air.",
            #min_value=0,
            #max_value=end_x,
            step=0.0001,
            format="%.4f",
        ),
        7: st.column_config.NumberColumn(
            "δ",
            help="Line position pressure shift (cm-1.atm-1). Value imported from HITRAN is for air.",
            #min_value=start_x,
            #max_value=end_x,
            step=0.0001,
            format="%.4f",
        ),
        8: st.column_config.NumberColumn(
            "± v0",
            help="Unceratinty in line position (cm-1)",
            min_value=1E-6,
            max_value=1,
            step=1E-6,
            format="%.1e",
        ),
        9: st.column_config.NumberColumn(
            "± S0",
            help="Relative unceratinty in line strength",
            min_value=0,
            max_value=100,
            step=0.1,
            format="%.1f %%",
        ),
        10: st.column_config.NumberColumn(
            "± γ_bath",
            help="Relative unceratinty in bath gas broadning coeffecient",
            min_value=0,
            max_value=100,
            step=0.1,
            format="%.1f %%",
        ),
        11: st.column_config.NumberColumn(
            "± γ_self",
            help="Relative unceratinty in self broadning coeffecient",
            min_value=0,
            max_value=100,
            step=0.1,
            format="%.1f %%",
        ),
        12: st.column_config.NumberColumn(
            "± n",
            help="Relative unceratinty in temperature exponents",
            min_value=0,
            max_value=100,
            step=0.1,
            format="%.1f %%",
        ),
        13: st.column_config.NumberColumn(
            "± δ",
            help="Unceratinty in pressure shift parameter (cm-1.atm-1)",
            min_value=1E-6,
            max_value=1,
            step=1E-6,
            format="%.1e",
        )
    }
        edited_lines = st.data_editor(lines,column_config=lines_column_config,key=st.session_state.dek)
    else:
        edited_lines = lines
    
    x0_rand_cdf, x0_rand_range,S0_rand_cdf, S0_rand_range,gamma_air_rand_cdf, gamma_air_rand_range\
    ,gamma_self_rand_cdf, gamma_self_rand_range,n_air_rand_cdf, n_air_rand_range,\
    delta_air_rand_cdf, delta_air_rand_range, x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air = extract_parameters(edited_lines)
    spectra = np.zeros((len(x), n_simulations))
    
    st.divider() 
    st.write('_Absorbance based on mean line-parameters of '+str(len(lines)) + ' lines,\nand '+str(n_simulations)+' spectra based on randomly sampled line-parameters:_')
    fig, ax = plt.subplots()
    #ax.set_title('Absorbance based on mean parameters and '+str(n_simulations)+' simulated spectra')
    ax.set_xlabel('Wavenumbers (cm-1)')
    ax.set_ylabel('Absorbance')
    
    spectra = MC_simulation(edited_lines,n_simulations,T,P,mole_fraction,L,x)
    spectrum_mean_parameters = mean_spectrum_simulation(edited_lines,T,P,mole_fraction,L,x)

    plot_MC_spectra(spectra, spectrum_mean_parameters)

    relative_uncertainty, error_bars, skewness = calc_error_bars(spectra,spectrum_mean_parameters)
    #plotting_commands()

    plot_uncertainty(relative_uncertainty,skewness)    

    if st.session_state.conv_test:
        std_deviation_with_iterations()

    simulation_info = [datetime.datetime.now(),selected_species,T,P,mole_fraction,L,wnstart,wnend,wnres,n_simulations,s0_min,st.session_state.manual_control,st.session_state.conv_test]
    #print(simulation_info)
    with open('simulation_history.csv','a') as fd:
        #fd.write(np.array2string(simulation_info))
        writer = csv.writer(fd)
        writer.writerow(simulation_info)

    arr = np.array(np.transpose([x,spectrum_mean_parameters,error_bars,relative_uncertainty,skewness]))
    arr_df = pd.DataFrame(arr,columns=['wavenumbers (cm-1)','absorbance - mean parameters','3 x standard deviation','relative uncertainty','skewness'])
    # Create an in-memory buffer
    with io.BytesIO() as buffer:
        # Write array to buffer
        #np.savetxt(buffer, arr_df, delimiter=",")
        textstr = ('MonteSpectra: '+selected_species + '_' + str(100*mole_fraction)+'% _' + str(T) + ' K_' + str(P) + ' atm_'+ str(L) + ' cm')
        buffer = pd.DataFrame.to_csv(arr_df, sep=',', index=False, encoding='utf-8')
        st.download_button(
            label="Download results as CSV",
            data = buffer, # Download buffer
            file_name = textstr+'.csv',
            mime='text/csv'
        ) 

st.divider() 
st.subheader('Sponsored by:')
#st.image(k_logo,width=50)
st.markdown(
    """<a href="https://www.kaust.edu.sa/">
    <img src="data:image/png;base64,{}" width="50">
    </a>""".format(
        base64.b64encode(open("images/kaust_2.png", "rb").read()).decode()
    ),
    unsafe_allow_html=True,
)