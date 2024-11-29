import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm, skew
from voigtfwhm import voigtfwhm
from voigtfwhm_fast import voigtfwhm_fast
import io
import base64
import uuid
import datetime
import csv
import time
#import mpld3
#import streamlit.components.v1 as components
#from streamlit_js_eval import streamlit_js_eval

st.set_page_config(
    page_title="MCSpectra",
    page_icon="🗻",
    layout="centered",
    initial_sidebar_state="expanded",
    menu_items={
        'Get Help': 'https://www.kaust.edu.sa',
        'Report a bug': "mailto:ihsan.farouki@kaust.edu.sa",
        'About': "# Absorbance spectra simulations with uncertainty quantification"
    }
)

#print(f"Screen width is {streamlit_js_eval(js_expressions='window.innerWidth')}")

# Sidebar to take user inputs
k_logo = "images/kaust_2.png"
MS_logo = "images/mcspectra_logo.png"
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

# change wavelength range depending on selected species
def change_wn_range():
    if st.session_state.selected_species == '(12)CH4 - HITRAN':
        st.session_state.wn_start = 1331
        st.session_state.wn_end = 1334
        #st.session_state.self_shift_toggle = 0
    elif st.session_state.selected_species == 'H2(16)O - HITRAN':
        st.session_state.wn_start = 3742
        st.session_state.wn_end = 3747
        #st.session_state.self_shift_toggle = 0
    elif st.session_state.selected_species == '(12)CO2 - HITRAN':
        st.session_state.wn_start = 2300
        st.session_state.wn_end = 2305
        #st.session_state.self_shift_toggle = 1
    elif st.session_state.selected_species == '(13)CO2 - HITRAN':
        st.session_state.wn_start = 2300
        st.session_state.wn_end = 2305
        #st.session_state.self_shift_toggle = 1
    elif st.session_state.selected_species == '(14)N2O - HITRAN':
        st.session_state.wn_start = 1285
        st.session_state.wn_end = 1290
        #st.session_state.self_shift_toggle = 1
    elif st.session_state.selected_species == '(12)CO - HITRAN':
        st.session_state.wn_start = 2172
        st.session_state.wn_end = 2182
        #st.session_state.self_shift_toggle = 1
    elif st.session_state.selected_species == '(14)NH3 - HITRAN':
        st.session_state.wn_start = 850
        st.session_state.wn_end = 856
        #st.session_state.self_shift_toggle = 0
    elif st.session_state.selected_species == '(12)C2H6 - HITRAN':
        st.session_state.wn_start = 3009
        st.session_state.wn_end = 3012
        #st.session_state.self_shift_toggle = 0

# molar mass of selected species
def molar_mass():
    if st.session_state.selected_species == '(12)CH4 - HITRAN':
        M = 16 # Molar mass of CH4 (g/mol)        
    elif st.session_state.selected_species == 'H2(16)O - HITRAN':
        M = 18 # Molar mass of CH4 (g/mol)       
    elif st.session_state.selected_species == '(12)CO2 - HITRAN':
        M = 44 # Molar mass of CH4 (g/mol)
    elif st.session_state.selected_species == '(13)CO2 - HITRAN':
        M = 45 # Molar mass of CH4 (g/mol)        
    elif st.session_state.selected_species == '(14)N2O - HITRAN':
        M = 44 # Molar mass of CH4 (g/mol)       
    elif st.session_state.selected_species == '(12)CO - HITRAN':
        M = 28 # Molar mass of CH4 (g/mol)
    elif st.session_state.selected_species == '(14)NH3 - HITRAN':
        M = 17 # Molar mass of CH4 (g/mol)
    elif st.session_state.selected_species == '(12)C2H6 - HITRAN':
        M = 30 # Molar mass of CH4 (g/mol)
        
    
    return M
        
species_options = ['(12)CH4 - HITRAN', 'H2(16)O - HITRAN', '(12)CO2 - HITRAN', '(13)CO2 - HITRAN', '(14)N2O - HITRAN','(12)CO - HITRAN','(14)NH3 - HITRAN','(12)C2H6 - HITRAN']

# preprogrammed list of broadeners for different species
if not(hasattr(st.session_state,'selected_species')):
    broadener_options = ['Air','H2O']
    self_shift_available = False
elif st.session_state.selected_species == '(12)CH4 - HITRAN':
    broadener_options = ['Air','H2O']
    self_shift_available = False
elif st.session_state.selected_species == 'H2(16)O - HITRAN':
    broadener_options = ['Air']
    self_shift_available = False
elif st.session_state.selected_species == '(12)CO2 - HITRAN':
    broadener_options = ['Air','H2','He','H2O']
    self_shift_available = True
elif st.session_state.selected_species == '(13)CO2 - HITRAN':
    broadener_options = ['Air','H2','He','H2O']
    self_shift_available = True
elif st.session_state.selected_species == '(14)N2O - HITRAN':
    broadener_options = ['Air','He','H2O']
    self_shift_available = True
elif st.session_state.selected_species == '(12)CO - HITRAN':
    broadener_options = ['Air','H2','He','CO2','H2O']
    self_shift_available = True
elif st.session_state.selected_species == '(14)NH3 - HITRAN':
    broadener_options = ['Air','H2','He','CO2','H2O']
    self_shift_available = False
elif st.session_state.selected_species == '(12)C2H6 - HITRAN':
    broadener_options = ['Air']
    self_shift_available = False

with st.sidebar:
    st.image(MS_logo, width=250)
    #st.divider() 
    #st.header("Simulation Controls")
    with st.expander('Basic simulation controls',True):
        simulation_type = st.selectbox("Spectrum type", ['Absorbance','Emission'], 0,key='simulation_type')
        selected_species = st.selectbox("Species", species_options, 0, on_change=change_wn_range,key='selected_species')
        selected_broadener = st.selectbox("Broadener", broadener_options, 0,key='selected_broadener')
        temperature = st.number_input("Temperature (K)", min_value=300, max_value=3000, value=300, step=100)
        pressure = st.number_input("Pressure (bar)", min_value=0.001, max_value=50.00, value=1.00, step=0.2)
        molefraction = st.number_input("Mole Fraction", min_value=0.00, max_value=1.00, value=0.01, step=0.001, format="%.3e")
        pathlength = st.number_input('Pathlength (cm)', min_value=1, max_value=50000, step=1, value=10)
        wnstart = st.number_input('Wavelength start (cm-1)', min_value=500.00, max_value=5000.00, step=0.01, value=1331.00, key='wn_start')
        wnend = st.number_input('Wavelength end (cm-1)', min_value=500.00, max_value=5000.00, step=1.00, value=1334.00, key='wn_end')
        wnres = st.number_input('Resolution (cm-1)', min_value=0.001, max_value=0.1, step=0.001, value=0.005, key='wn_res', format="%.3f")
    #st.divider() 
    #st.subheader('Advanced simulation controls')
    with st.expander('Advanced simulation controls',True):
        N_simulations = st.number_input('Number of simulations', min_value=1, max_value=2000, step=100, value=1000)
        N_PDF_points = st.number_input('PDF resolution (Number of points)', min_value=50, max_value=100, step=1, value=100)
        s0_min_input = st.number_input("Line strength threshold (cm-1/(molec.cm-2))", min_value=1E-25, max_value=1E-19, value=1E-21, format="%.3e", disabled=True,key='s0_min')
        max_residual = st.number_input("Max. allowed risidual due to frequency cut-off", min_value=1E-3, max_value=1E-2, value=1E-2, format="%.3e")
        #conv_test = st.toggle("Plot standard deviation vs iterations to test convergence",key='conv_test',value=1, disabled=1)
        conv_test = 1
        if conv_test == 1:
            convergence_frequency = st.number_input('Cursor position for convergence test and PDF (cm-1)', min_value=min(st.session_state.wn_start,st.session_state.wn_end), max_value=max(st.session_state.wn_start,st.session_state.wn_end), value=min(st.session_state.wn_start,st.session_state.wn_end), key='wn_conv')
        manual_control = st.toggle("Enable manual control of line parameters",key='manual_control')
        calc_method = st.toggle("Less accurate evaluation of the Voigt function",key='calc_method')
        exp_unc = st.toggle("Account for uncertainty in experimental conditions",key='exp_unc')
        #self_shift_toggle = st.toggle("Account for pressure self-shift",key='self_shift_toggle') #,disabled=self_shift_available
        if st.session_state.exp_unc:
            molefraction_unc = st.number_input('Uncertainty in mole fraction (%)', min_value=0, max_value=100, value=0, key='molefraction_unc')
            pathlength_unc = st.number_input('Uncertainty in pathlength (%)', min_value=0, max_value=100, value=0, key='pathlength_unc')
            pressure_unc = st.number_input('Uncertainty in pressure (%)', min_value=0, max_value=100, value=0, key='pressure_unc')
            temperature_unc = st.number_input('Uncertainty in temperature (%)', min_value=0, max_value=100, value=0, key='temperature_unc')

    #st.divider() 
    st.subheader('Warnings')

tab1, tab2, tab3, tab4 = st.tabs(["Simulated spectra","Global statistics","Convergence","PDF"])

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
exp_unc_values = [0,0,0,0]
if st.session_state.exp_unc:
    exp_unc_values = [molefraction_unc,pathlength_unc,pressure_unc,temperature_unc]
    #print(exp_unc_values)
np.set_printoptions(precision=16)


# import lines and tips data for the selected lines
@st.cache_resource(show_spinner=False,max_entries=3)
def import_data(selected_species):
    print('importing data')
    # Load data (replace readmatrix and readtable)
    if selected_species == '(12)CH4 - HITRAN':
        CH4lines = pd.read_csv('12CH4_lines_formatted.csv').values
        tips = pd.read_csv('q32_12CH4.csv', sep='\s+').values
    elif selected_species == 'H2(16)O - HITRAN':
        CH4lines = pd.read_csv('H216O_lines_formatted.csv').values
        tips = pd.read_csv('q1_H2O16.csv', sep='\s+').values
    elif selected_species == '(12)CO2 - HITRAN':
        CH4lines = pd.read_csv('12CO2_lines_formatted.csv').values
        tips = pd.read_csv('q7_12CO2.csv', sep='\s+').values
    elif selected_species == '(13)CO2 - HITRAN':
        CH4lines = pd.read_csv('13CO2_lines_formatted.csv').values
        tips = pd.read_csv('q8_13CO2.csv', sep='\s+').values
    elif selected_species == '(14)N2O - HITRAN':
        CH4lines = pd.read_csv('14N2O_lines_formatted.csv').values
        tips = pd.read_csv('q21_14N2O.csv', sep='\s+').values
    elif selected_species == '(12)CO - HITRAN':
        CH4lines = pd.read_csv('12CO_lines_formatted.csv').values
        tips = pd.read_csv('q26_12CO.csv', sep='\s+').values
    elif selected_species == '(14)NH3 - HITRAN':
        CH4lines = pd.read_csv('14NH3_lines_formatted.csv').values
        tips = pd.read_csv('q45_14NH3.csv', sep='\s+').values
    elif selected_species == '(12)C2H6 - HITRAN':
        CH4lines = pd.read_csv('12C2H6_lines_formatted.csv').values
        tips = pd.read_csv('q78_12C2H6.csv', sep='\s+').values
    return CH4lines, tips

# Define interpolation function for tips data (equivalent to tips1 = @(z) in MATLAB)
def tips1(z):
    return np.interp(z, tips[:, 0], tips[:, 1])

def plank_emission(x,T):
    Ibb = 2*h*(c**2)*np.divide(np.power(x,3),(np.exp(c*h*x/(kb*T))-1))
    return Ibb

np.set_printoptions(legacy='1.25')

# start and end indices for the lines going into the calculation
@st.cache_resource(show_spinner=False,max_entries=3)
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

# Get parameters and their uncertainties for lines within the range
@st.cache_resource(show_spinner=False,max_entries=3)
def extract_lines(start_x,end_x,CH4lines,s0_min,selected_broadener, testing_range):
    print('extracting lines')
    flags_array = np.zeros(23)
    # %% Get parameters and their uncertainties
    lines = []
    j = 0
    for i in range(start_x,end_x): #range(len(CH4lines)):
        if CH4lines[i,1] > s0_min: #(CH4lines[i, 0] > min(x)) and (CH4lines[i, 0] < max(x)):
            #print(CH4lines[i])
            j += 1
            #print(flags_array[11] == 0)
            if selected_broadener == 'Air':
                line = [
                    CH4lines[i, 0],  # line position
                    CH4lines[i, 1],  # line strength
                    CH4lines[i, 2],  # gamma air
                    CH4lines[i, 3],  # gamma self
                    CH4lines[i, 4],  # LES
                    CH4lines[i, 5],  # n_air
                    CH4lines[i, 6]   # delta_air
                ]
            elif selected_broadener == 'H2':
                line = [
                    CH4lines[i, 0],  # line position
                    CH4lines[i, 1],  # line strength
                    CH4lines[i, 11],  # gamma H2
                    CH4lines[i, 3],  # gamma self
                    CH4lines[i, 4],  # LES
                    CH4lines[i, 12],  # n_H2
                    CH4lines[i, 13]   # delta_H2
                ]
                if str(line[5]) == 'nan':
                    line[5] = 0
                    if (not(testing_range)) & (flags_array[1] == 0):
                        with st.sidebar:
                            flags_array[1] = 1
                            st.sidebar.warning('Broadener temperature exponent is unavailable for '+str(CH4lines[i, 0])+'. Assumed to be equal to zero.', icon="⚠️")
                if str(line[6]) == 'nan':
                    line[6] = 0
                    if (not(testing_range)) & (flags_array[2] == 0):
                        with st.sidebar:
                            flags_array[2] =1
                            st.sidebar.warning('Broadener pressure shift parameter is unavailable for '+str(CH4lines[i, 0])+'. Assumed to be equal to zero.', icon="⚠️")
            elif selected_broadener == 'He':
                line = [
                    CH4lines[i, 0],  # line position
                    CH4lines[i, 1],  # line strength
                    CH4lines[i, 15],  # gamma He
                    CH4lines[i, 3],  # gamma self
                    CH4lines[i, 4],  # LES
                    CH4lines[i, 16],  # n_He
                    CH4lines[i, 17]   # delta_He
                ]
                if str(line[5]) == 'nan':
                    line[5] = 0
                    if (not(testing_range)) & (flags_array[3] == 0):
                        with st.sidebar:
                            flags_array[3] =1
                            st.sidebar.warning('Broadener temperature exponent is unavailable for '+str(CH4lines[i, 0])+'. Assumed to be equal to zero.', icon="⚠️")
                if str(line[6]) == 'nan':
                    line[6] = 0
                    if (not(testing_range)) & (flags_array[4] == 0):
                        with st.sidebar:
                            flags_array[4] =1
                            st.sidebar.warning('Broadener pressure shift parameter is unavailable for '+str(CH4lines[i, 0])+'. Assumed to be equal to zero.', icon="⚠️")
            elif selected_broadener == 'CO2':
                line = [
                    CH4lines[i, 0],  # line position
                    CH4lines[i, 1],  # line strength
                    CH4lines[i, 18],  # gamma CO2
                    CH4lines[i, 3],  # gamma self
                    CH4lines[i, 4],  # LES
                    CH4lines[i, 19],  # n_CO2
                    CH4lines[i, 20]   # delta_CO2
                ]
                if str(line[5]) == 'nan':
                    line[5] = 0
                    if (not(testing_range)) & (flags_array[5] == 0):
                        with st.sidebar:
                            flags_array[5] =1
                            st.sidebar.warning('Broadener temperature exponent is unavailable for '+str(CH4lines[i, 0])+'. Assumed to be equal to zero.', icon="⚠️")
                if str(line[6]) == 'nan':
                    line[6] = 0
                    if (not(testing_range)) & (flags_array[6] == 0):
                        with st.sidebar:
                            flags_array[6] =1
                            st.sidebar.warning('Broadener pressure shift parameter is unavailable for '+str(CH4lines[i, 0])+'. Assumed to be equal to zero.', icon="⚠️")
            elif selected_broadener == 'H2O':
                line = [
                    CH4lines[i, 0],  # line position
                    CH4lines[i, 1],  # line strength
                    CH4lines[i, 21],  # gamma H2O
                    CH4lines[i, 3],  # gamma self
                    CH4lines[i, 4],  # LES
                    CH4lines[i, 22],  # n_H2O
                    0   # delta_He
                ]
                if str(line[5]) == 'nan':
                    line[5] = 0
                    if (not(testing_range)) & (flags_array[7] == 0):
                        with st.sidebar:
                            flags_array[7] =1
                            st.sidebar.warning('Broadener temperature exponent is unavailable for '+str(CH4lines[i, 0])+'. Assumed to be equal to zero.', icon="⚠️")
                
                if (not(testing_range)):
                    with st.sidebar:
                        st.sidebar.warning('Broadener pressure shift parameter is unavailable for '+str(CH4lines[i, 0])+'. Assumed to be equal to zero.', icon="⚠️")

            # Uncertainty handling (equivalent to switch cases in MATLAB)
            for k in [1]:
                uncertainty = CH4lines[i, 22 + k]
                if uncertainty == 0:
                    if (not(testing_range)) & (flags_array[8] == 0):
                        with st.sidebar:
                            flags_array[8] =1
                            st.sidebar.warning('Uncertainty in line position is unavailable for '+str(CH4lines[i, 0]), icon="⚠️")
                    line.append(0)
                elif uncertainty == 1:
                    if (not(testing_range)) & (flags_array[9] == 0):
                        with st.sidebar:
                            flags_array[9] =1
                            st.sidebar.warning('Uncertainty in line position (parameter '+str(k)+') for '+str(np.round(CH4lines[i, 0],2)) + ' might be larger than visible in the simulation result', icon="⚠️")
                    line.append(1E-1)
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
                    line.append(1E-5)

            for k in [2]:
                uncertainty = CH4lines[i, 22 + k]
                if uncertainty in [0, 1, 2]:
                    if (not(testing_range))  & (flags_array[10] == 0):
                        with st.sidebar:
                            flags_array[10] =1
                            st.sidebar.warning('Uncertainty in line strength is unavailable for '+str(np.round(CH4lines[i, 0],2)), icon="⚠️")
                    line.append((100)*0)
                elif uncertainty == 3:
                    if (not(testing_range)) & ((flags_array[11] == 0)):
                        with st.sidebar:
                            #print('not testing range')
                            #print(not(testing_range))
                            #print('flag down')
                            #print(flags_array[11] == 0)
                            #print('overall')
                            #print(not(testing_range) & (flags_array[11] == 0))
                            flags_array[11] =1
                            st.sidebar.warning('Uncertainty in line strength for '+str(np.round(CH4lines[i, 0],2)) + ' might be larger than visible in the simulation result', icon="⚠️")
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
            
            if selected_broadener == 'Air':
                temp_index = 3
            elif selected_broadener == 'H2':
                temp_index = 11
            elif selected_broadener == 'He':
                temp_index = 15
            elif selected_broadener == 'CO2':
                temp_index = 18
            elif selected_broadener == 'H2O':
                temp_index = 21

            for k in [temp_index]:
                uncertainty = CH4lines[i, 22 + k]
                if uncertainty in [0, 1, 2]:
                    if (not(testing_range)) & (flags_array[12] == 0):
                        with st.sidebar:
                            flags_array[12] =1
                            st.sidebar.warning('Uncertainty in bath gas broadening coeffecient is unavailable for '+str(np.round(CH4lines[i, 0],2)), icon="⚠️")
                    line.append((100)*0)
                elif uncertainty == 3:
                    if (not(testing_range)) & (flags_array[13] == 0):
                        with st.sidebar:
                            flags_array[13] =1
                            st.sidebar.warning('Uncertainty in bath gas broadening coeffecient for '+str(np.round(CH4lines[i, 0],2)) + ' might be larger than visible in the simulation result', icon="⚠️")
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

            for k in [4]:
                uncertainty = CH4lines[i, 22 + k]
                if uncertainty in [0, 1, 2]:
                    if (not(testing_range)) & (flags_array[14] == 0):
                        with st.sidebar:
                            flags_array[14] =1
                            st.sidebar.warning('Uncertainty in self broadening coeffecient is unavailable for '+str(np.round(CH4lines[i, 0],2)), icon="⚠️")
                    line.append((100)*0)
                elif uncertainty == 3:
                    if (not(testing_range)) & (flags_array[15] == 0):
                        with st.sidebar:
                            flags_array[15] =1
                            st.sidebar.warning('Uncertainty in self broadening coeffecient is unavailable for '+str(np.round(CH4lines[i, 0],2)) + ' might be larger than visible in the simulation result', icon="⚠️")
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

            if selected_broadener == 'Air':
                temp_index = 5
            elif selected_broadener == 'H2':
                temp_index = 12
            elif selected_broadener == 'He':
                temp_index = 16
            elif selected_broadener == 'CO2':
                temp_index = 19
            elif selected_broadener == 'H2O':
                temp_index = 22

            for k in [temp_index]:
                uncertainty = CH4lines[i, 22 + k]
                if uncertainty in [0, 1, 2]:
                    if (not(testing_range)) & (flags_array[16] == 0):
                        with st.sidebar:
                            flags_array[16] =1
                            st.sidebar.warning('Uncertainty in bath gas broadening coeffecient temperature exponent is unavailable for '+str(np.round(CH4lines[i, 0],2)), icon="⚠️")
                    line.append((100)*0)
                elif uncertainty == 3:
                    if (not(testing_range)) & (flags_array[17] == 0):
                        with st.sidebar:
                            flags_array[17] =1
                            st.sidebar.warning('Uncertainty in bath gas broadening coeffecient temperature exponent for '+str(np.round(CH4lines[i, 0],2)) + ' might be larger than visible in the simulation result', icon="⚠️")
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
            
            if selected_broadener == 'Air':
                temp_index = 6
            elif selected_broadener == 'H2':
                temp_index = 13
            elif selected_broadener == 'He':
                temp_index = 17
            elif selected_broadener == 'CO2':
                temp_index = 20
            elif selected_broadener == 'H2O':
                temp_index = 22
            
            for k in [temp_index]:
                uncertainty = CH4lines[i, 22 + k]
                if uncertainty == 0:
                    if (not(testing_range)) & (flags_array[18] == 0):
                        with st.sidebar:
                            flags_array[18] = 1
                            st.sidebar.warning('Uncertainty in bath gas pressure shift is unavailable for '+str(CH4lines[i, 0]), icon="⚠️")
                    line.append(0)
                elif uncertainty == 1:
                    if (not(testing_range)) & (flags_array[19] == 0):
                        with st.sidebar:
                            flags_array[19] = 1
                            st.sidebar.warning('Uncertainty in bath gas pressure shift for '+str(np.round(CH4lines[i, 0],2)) + ' might be larger than visible in the simulation result', icon="⚠️")
                    line.append(1E-1)
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
                    line.append(1E-5)
            

            if not(self_shift_available):
                if (not(testing_range)) & (flags_array[20] == 0):
                    with st.sidebar:
                        flags_array[20] =1
                        st.sidebar.warning('Self pressure shift data is not available for this species. Assumed equal to zero, but can be modified using manual control mode.', icon="⚠️")
                line.append(0)
                line.append(0)
            else:
                if (not(testing_range)) & (flags_array[21] == 0):
                    with st.sidebar:
                        flags_array[21] =1
                        st.sidebar.info('Self pressure shift data is available for this species.', icon="ℹ️")
                line.append(CH4lines[i, 8])
                for k in [8]:
                    uncertainty = CH4lines[i, 22 + k]
                    if uncertainty == 0:
                        if (not(testing_range)) & (flags_array[22] == 0):
                            with st.sidebar:
                                flags_array[22] = 1
                                st.sidebar.warning('Uncertainty in self pressure shift is unavailable for '+str(CH4lines[i, 0]), icon="⚠️")
                        line.append(0)
                    elif uncertainty == 1:
                        if (not(testing_range)) & (flags_array[23] == 0):
                            flags_array[23] = 1
                            with st.sidebar:
                                st.sidebar.warning('Uncertainty in self pressure shift for '+str(np.round(CH4lines[i, 0],2)) + ' might be larger than visible in the simulation result', icon="⚠️")
                        line.append(1E-1)
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
                        line.append(1E-5)

            
            #print(line)
            lines.append(line)
    #st.text('Number of simulated lines: '+str(len(lines)))
    #my_spinner.empty()
    return lines

# %% Simulate spectra with uncertainty
def rand_distribution(mu, sigma):
    """Equivalent to rand_distrition in MATLAB."""
    range_vals = np.linspace(mu - 3 * sigma, mu + 3 * sigma, num_of_PDF_points)
    cdf_vals = norm.cdf(range_vals, loc=mu, scale=sigma)
    return cdf_vals, range_vals

def random_value(cdf, range_vals):
    """Equivalent to random_value in MATLAB."""
    return np.interp(np.random.rand(), cdf, range_vals)

# Generate distributions for line parameters
@st.cache_resource(show_spinner=False,max_entries=3)
def extract_parameters(lines):
    print('extracting parameters and generating distributions')
    x0_rand_cdf = np.zeros((len(lines),num_of_PDF_points))
    x0_rand_range = np.zeros((len(lines),num_of_PDF_points))
    S0_rand_cdf = np.zeros((len(lines),num_of_PDF_points))
    S0_rand_range = np.zeros((len(lines),num_of_PDF_points))
    gamma_air_rand_cdf = np.zeros((len(lines),num_of_PDF_points))
    gamma_air_rand_range = np.zeros((len(lines),num_of_PDF_points))
    gamma_self_rand_cdf = np.zeros((len(lines),num_of_PDF_points))
    gamma_self_rand_range = np.zeros((len(lines),num_of_PDF_points))
    n_air_rand_cdf = np.zeros((len(lines),num_of_PDF_points))
    n_air_rand_range = np.zeros((len(lines),num_of_PDF_points))
    delta_air_rand_cdf = np.zeros((len(lines),num_of_PDF_points))
    delta_air_rand_range = np.zeros((len(lines),num_of_PDF_points))
    delta_self_rand_cdf = np.zeros((len(lines),num_of_PDF_points))
    delta_self_rand_range = np.zeros((len(lines),num_of_PDF_points))

    x0 = np.zeros(len(lines))
    s0 = np.zeros(len(lines))
    gamma_air_0 = np.zeros(len(lines))
    gamma_self_0 = np.zeros(len(lines))
    n_air = np.zeros(len(lines))
    delta_air = np.zeros(len(lines))
    delta_self = np.zeros(len(lines))


    #fig, ax = plt.subplots()

    j=0
    for line in lines:
        #print(line)

        x0_rand_cdf[j], x0_rand_range[j] = rand_distribution(line[0], line[7] / 3)
        S0_rand_cdf[j], S0_rand_range[j] = rand_distribution(line[1], line[1] * (1/100)*line[8] / 3)
        gamma_air_rand_cdf[j], gamma_air_rand_range[j] = rand_distribution(line[2], line[2] * (1/100)*line[9] / 3)
        gamma_self_rand_cdf[j], gamma_self_rand_range[j] = rand_distribution(line[3], line[3] * (1/100)*line[10] / 3)
        n_air_rand_cdf[j], n_air_rand_range[j] = rand_distribution(line[5], line[5] * (1/100)*line[11] / 3)
        delta_air_rand_cdf[j], delta_air_rand_range[j] = rand_distribution(line[6], line[12] / 3)
        delta_self_rand_cdf[j], delta_self_rand_range[j] = rand_distribution(line[13], line[14] / 3)

        x0[j] = line[0]
        s0[j] = line[1]
        gamma_air_0[j] = line[2]
        gamma_self_0[j] = line[3]
        n_air[j] = line[5]
        delta_air[j] = line[6]
        delta_self[j] = line[13]

        j = j + 1

    return  x0_rand_cdf, x0_rand_range,S0_rand_cdf, S0_rand_range,gamma_air_rand_cdf, gamma_air_rand_range\
    ,gamma_self_rand_cdf, gamma_self_rand_range,n_air_rand_cdf, n_air_rand_range, delta_self_rand_cdf, delta_self_rand_range,\
    delta_air_rand_cdf, delta_air_rand_range, x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self

@st.cache_resource(show_spinner=False,max_entries=3)
def extract_mean_parameters(lines):
    print('extracting mean parameters')

    x0 = np.zeros(len(lines))
    s0 = np.zeros(len(lines))
    gamma_air_0 = np.zeros(len(lines))
    gamma_self_0 = np.zeros(len(lines))
    n_air = np.zeros(len(lines))
    delta_air = np.zeros(len(lines))
    delta_self = np.zeros(len(lines))

    #fig, ax = plt.subplots()

    j=0
    temp_flag_start = False
    temp_flag_end = False
    visible_start = 0
    visible_end = 0
    for line in lines:
        #print(line)
        if (line[0] >= wnstart):
            if temp_flag_start == False:
                visible_start = j
                temp_flag_start = True
        if (line[0] >= wnend):
            if temp_flag_end == False:
                visible_end = j
                temp_flag_end = True
        
        x0[j] = line[0]
        s0[j] = line[1]
        gamma_air_0[j] = line[2]
        gamma_self_0[j] = line[3]
        n_air[j] = line[5]
        delta_air[j] = line[6]
        delta_self[j] = line[13]

        j = j + 1
    
    if temp_flag_end == False:
        visible_end = j

    return  x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self, visible_start, visible_end

# Generate distributions for experimental conditions
def exp_unc_distributions(mole_fraction, pathlength, pressure, temperature,exp_unc_values):
    print('generating distributions for experimental uncertainties')
    molefraction_cdf, molefraction_range = rand_distribution(mole_fraction, mole_fraction * (1/100)*exp_unc_values[0] / 3)
    pathlength_cdf, pathlength_range = rand_distribution(pathlength, pathlength * (1/100)*exp_unc_values[1] / 3)
    pressure_cdf, pressure_range = rand_distribution(pressure, pressure * (1/100)*exp_unc_values[2] / 3)
    temperature_cdf, temperature_range = rand_distribution(temperature, temperature * (1/100)*exp_unc_values[3] / 3)

    return molefraction_cdf, molefraction_range, pathlength_cdf, pathlength_range, pressure_cdf, pressure_range, temperature_cdf, temperature_range

# Run the simulations
@st.cache_resource(show_spinner=False)
def MC_simulation(lines,n_simulations,T,P,mole_fraction,L,x,exp_unc_values,calc_method,simulation_type):
    # Run the simulations
    with tab1:
        my_bar = st.progress(0, text='Monte Carlo simulation progress:')
    spectra = np.zeros((len(x), n_simulations))
    spectra_1 = np.zeros((len(x_limited), n_simulations))

    if not(st.session_state.exp_unc & (exp_unc_values != [0,0,0,0])):
        mole_fraction_1 = mole_fraction
        L_1 = L
        P_1 = P
        T_1 = T

    for i in range(n_simulations):
        j=0
        if (st.session_state.exp_unc & (exp_unc_values != [0,0,0,0])):
            mole_fraction_1 = random_value(molefraction_cdf, molefraction_range)
            L_1 = random_value(pathlength_cdf, pathlength_range)
            P_1 = random_value(pressure_cdf, pressure_range)
            T_1 = random_value(temperature_cdf, temperature_range)
            
        for line in lines:
            # Line position
            #t = time.time()
            x0_rand = random_value(x0_rand_cdf[j], x0_rand_range[j])

            delta_air_rand = random_value(delta_air_rand_cdf[j], delta_air_rand_range[j])
            delta_self_rand = random_value(delta_self_rand_cdf[j], delta_self_rand_range[j])
            x0_shifted_rand = x0_rand + P_1 * ((1 - mole_fraction_1) * delta_air_rand + mole_fraction_1*delta_self_rand)

            # Line strength

            S0_rand = random_value(S0_rand_cdf[j], S0_rand_range[j])
            S_rand = S0_rand * (tips1(296) / tips1(T)) * np.exp(-(h * c * line[4] / kb) * (1 / T_1 - 1 / 296)) \
                     * (1 - np.exp(-h * c * x0_rand / (kb * T_1))) / (1 - np.exp(-h * c * x0_rand / (kb * 296)))

            A_rand = S_rand * L_1 * mole_fraction_1 * (P_1 / (R * T_1))  # cm-1

            # Doppler broadening
            wG_rand = x0_shifted_rand * (7.1623E-7) * np.sqrt(T_1 / M)

            # Pressure broadening
            n_air_rand = random_value(n_air_rand_cdf[j], n_air_rand_range[j])

            gamma_self_rand = random_value(gamma_self_rand_cdf[j], gamma_self_rand_range[j])

            gamma_air_rand = random_value(gamma_air_rand_cdf[j], gamma_air_rand_range[j])

            gamma_self_rand = gamma_self_rand * (296 / T_1) ** n_air_rand
            gamma_air_rand = gamma_air_rand * (296 / T_1) ** n_air_rand

            wL_rand = P_1 * (mole_fraction_1 * 2 * gamma_self_rand + (1 - mole_fraction_1) * 2 * gamma_air_rand)

            #elapsed = time.time() - t
            #print('time to sample parameters for a single line')
            #print(elapsed)

            #t = time.time()

            #spectra[:, i] += (x, [A_rand, x0_shifted_rand, wG_rand, wL_rand])
            if not (st.session_state.calc_method):
                spectra[:, i] += voigtfwhm(x, [A_rand, x0_shifted_rand, wG_rand, wL_rand])
            else:
                X = np.sqrt(np.log(2))*(x-x0_shifted_rand)/(0.5*wG_rand)
                Y = np.sqrt(np.log(2))*((0.5*wL_rand)/(0.5*wG_rand))
                spectra[:, i] += voigtfwhm_fast(x, [A_rand, X, Y])/ ((0.5*wG_rand) / np.sqrt(np.log(2)/np.pi))

            #elapsed = time.time() - t
            #print('time to calculate and add spectrum for a single line based on sampled parameters')
            #print(elapsed)

            j=j+1

        if np.isnan(np.mean(spectra[:, i])):
            spectra[:, i] = spectra[:, i - 1]
        
        if simulation_type == 'Emission':
            spectra[:, i] = np.multiply(plank_emission(x,T),(1-np.exp(-spectra[:, i])))
        #else:
            #ax.plot(x, spectra[:, i], '.', markersize=0.1, color="#A87BF9")

        spectra_1[:, i] = np.interp(x_limited, x, spectra[:, i])
        my_bar.progress(i/n_simulations, text='Monte Carlo simulation progress')
        
    my_bar.empty()
    return spectra_1

# Calculate mean spectrum
@st.cache_resource(show_spinner=False,max_entries=3)
def mean_spectrum_simulation(lines,T,P,mole_fraction,L,x,calc_method,simulation_type):
    spectrum_mean_parameters = np.zeros(len(x))
    j = 0
    for line in lines:
        # Line position

        x0_shifted = x0[j] + P * ((1 - mole_fraction) * delta_air[j] + mole_fraction * delta_self[j])
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
        if not (st.session_state.calc_method):
            spectrum_mean_parameters +=  np.transpose(voigtfwhm(x, [A, x0_shifted, wG, wL]))
        else:
            X = np.sqrt(np.log(2))*(x-x0_shifted)/(0.5*wG)
            Y = np.sqrt(np.log(2))*((0.5*wL)/(0.5*wG))
            spectrum_mean_parameters += voigtfwhm_fast(x, [A, X, Y])/ ((0.5*wG) / np.sqrt(np.log(2)/np.pi))        
        
        j=j+1

    if simulation_type == 'Emission':      
        spectrum_mean_parameters = np.multiply(plank_emission(x,T),(1-np.exp(-spectrum_mean_parameters)))
        
    return spectrum_mean_parameters

@st.cache_resource(show_spinner=False,max_entries=3)
def calc_error_bars(spectra,spectrum_mean_parameters):
    error_bars = np.zeros(len(x_limited))
    relative_uncertainty = np.zeros(len(x_limited))
    skewness = np.zeros(len(x_limited))
    for i in range(len(x_limited)):
        relative_uncertainty[i] = 100*3*np.std(spectra[i][:])/np.mean(spectrum_mean_parameters[i])
        skewness[i] = skew(spectra[i][:])
        error_bars[i] = 3*np.std(spectra[i][:])
    
    return relative_uncertainty, error_bars, skewness
  
def std_deviation_with_iterations(spectra,spectrum_mean_parameters):
    std_residuals = np.zeros(N_simulations)
    #max_std_index = np.argmax(error_bars)
    std_index = np.where(np.round(x_limited,3) == st.session_state.wn_conv)[0][0]
    #print(std_index)
    #print('at ('+str(round(x[std_index],2))+' cm-1')
    #print(max_std_index)
    for i in range(2,n_simulations):
        std_residuals[i] = 100*3*np.std(spectra[std_index][range(i)])/spectrum_mean_parameters[std_index]#/np.std(spectra[i][range(n_simulations)])
    
    #np.savetxt('sample_results/convergence_for_analysis.csv', std_residuals, delimiter=',')

    fig_3, ax = plt.subplots()
    #print(std_residuals)
    ax.plot(range(n_simulations), std_residuals,color="#ECBC7A")
    ax.set_xlabel('Iteration')
    ax.set_ylabel('3 x std. dev. (%)')
    ax.grid(visible=True, linestyle='--', linewidth=0.5)
    ax.set_xlim(0,n_simulations)

    #tab3 = st.tabs(['Covergence'])
    #st.divider() 
    return fig_3, std_index
    
    #return std_residuals

@st.cache_resource(show_spinner=False,max_entries=3)
def uncertainty_PDF(spectra):
    #std_residuals = np.zeros(N_simulations)
    #max_std_index = np.argmax(error_bars)
    std_index = np.where(np.round(x_limited,3) == st.session_state.wn_conv)[0][0]
    #print(std_index)
    #print('at ('+str(round(x[std_index],2))+' cm-1')
    #print(max_std_index)
    #for i in range(2,n_simulations):
    #    std_residuals[i] = np.std(spectra[std_index])#/np.std(spectra[i][range(n_simulations)])
    N_points = len(spectra[std_index])
    n_bins = 50
    fig_4, ax = plt.subplots()
    #print(std_residuals)
    if simulation_type == 'Absorbance':
        ax.set_xlabel('Absorbance')
        type_scale = 1
    else:
        ax.set_xlabel('Emission (µW/(cm-1-cm2-sr)')
        type_scale = 1E+6
    ax.set_ylabel('Frequency')
    ax.grid(visible=True, linestyle='--', linewidth=0.5)
    ax.hist(type_scale*spectra[std_index], bins=n_bins, color="#A87BF9")
    ax.axvline(type_scale*spectrum_mean_parameters[std_index], color='white', linestyle='dashed', linewidth=2)
    #textstr = (str(100*mole_fraction)+'% ' + selected_species + '\n' + str(T) + ' K\n' + str(P) + ' atm\n'+ str(L) + ' cm\n' + 'broadener: '+ selected_broadener +'\n' + '# of simulations: ' + str(n_simulations))
    textstr = (str(100*mole_fraction)+'% ' + selected_species + '\n' + str(T) + ' K\n' + str(P) + ' atm\n'+ str(L) + ' cm\n' + 'broadener: '+ selected_broadener +'\n' + '# of simulations: ' + str(n_simulations) +'\n' + 'Skewness: ' + str(np.round(skewness[std_index],4))+'\n' + 'Std. Dev.: ' + str(np.round(0.333*error_bars[std_index],4)))
    props = dict(boxstyle='round', facecolor="#A87BF9", alpha=0)
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props) 
    #ax.set_xlim(0,n_simulations)

    #tab3 = st.tabs(['Covergence'])
    #st.divider() 
    return fig_4, std_index
    
    #return std_residuals

@st.cache_resource(show_spinner=False,max_entries=3)
def plot_MC_spectra(spectra, spectrum_mean_parameters):
    fig_1, ax = plt.subplots()
    #ax.set_title('Absorbance based on mean parameters and '+str(n_simulations)+' simulated spectra')
    ax.set_xlabel('Wavenumbers (cm-1)')
    if simulation_type == 'Absorbance':
        ax.set_ylabel('Absorbance')
        type_scale = 1
    else:
        ax.set_ylabel('Emission (µW/(cm-1-cm2-sr))')
        type_scale = 1E+6

    ax.plot(x_limited, type_scale*spectrum_mean_parameters, '-', color='white',zorder=n_simulations+1)
    #ax.plot(x, type_scale*spectrum_mean_parameters, '-', color='black', zorder=n_simulations+1)


    for i in range(n_simulations):
        ax.plot(x_limited, type_scale*spectra[:, i], '-', lw=1, color="#A87BF9")
    
    ax.set_xlim(wnstart,wnend)
    #print(np.ceil(10*spectra.max())/10)
    #print(spectra.max())
    #scaling_factor = 1*1/(spectra.max())
    ax.set_ylim(0,type_scale*spectra.max())#np.round(scaling_factor*spectra.max())/scaling_factor)
    ax.legend(['Mean parameters', 'Unceratinty envelope'])

    textstr = (str(100*mole_fraction)+'% ' + selected_species + '\n' + str(T) + ' K\n' + str(P) + ' atm\n'+ str(L) + ' cm\n' + 'broadener: '+ selected_broadener +'\n' + '# of simulations: ' + str(n_simulations))
    props = dict(boxstyle='round', facecolor="#A87BF9", alpha=0)
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)

    ax.grid(visible=True, linestyle='--', linewidth=0.5)

    secax = ax.secondary_xaxis('top', functions=(wn2wl, wl2wn))
    secax.set_xlabel('Wavelength (µm)')
    
    return fig_1

    #fig_html = mpld3.fig_to_html(fig)
    #components.html(fig_html, height=600)

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

@st.cache_resource(show_spinner=False,max_entries=3)
def plot_uncertainty(relative_uncertainty,skewness):
    fig_2, ax1 = plt.subplots()

    color = (1,1,1,1)
    ax1.set_xlabel('Wavenumbers (cm-1)')
    ax1.set_ylabel('3 x std. dev. (%)', color=color)
    ax1.plot(x_limited, relative_uncertainty, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    if max(relative_uncertainty) < 100:
        ax1.set_ylim(0,100)
    ax1.set_xlim(wnstart,wnend)
    ax1.grid(visible=True, linestyle='--', linewidth=0.5)

    secax = ax1.secondary_xaxis('top', functions=(wn2wl, wl2wn))
    secax.set_xlabel('Wavelength (µm)')

    ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

    color = ("#57D2E9")
    ax2.set_ylabel('Skewness', color=color)  # we already handled the x-label with ax1
    ax2.plot(x_limited, skewness, color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    if abs(max(skewness)) < 1:
        ax2.set_ylim(-1,1)

    textstr = (str(100*mole_fraction)+'% ' + selected_species + '\n' + str(T) + ' K\n' + str(P) + ' atm\n'+ str(L) + ' cm\n' + 'broadener: '+ selected_broadener +'\n' + '# of simulations: ' + str(n_simulations))
    props = dict(boxstyle='round', facecolor="#A87BF9", alpha=0.1)
    #ax1.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10, verticalalignment='top', bbox=props)
    ax1.annotate(textstr, xy=(0, 1), xytext=(12, -12), va='top', xycoords='axes fraction', textcoords='offset points')

    fig_2.tight_layout()  # otherwise the right y-label is slightly clipped
    #st.divider() 
    
    #tab2 = st.tabs(['Global statistics'])
    return fig_2

@st.cache_resource(show_spinner=False)
def update_manual_control(lines):
    st.session_state.dek = str(uuid.uuid4()) # refresh key to reset lines


wn_validation_flag, wn_change_flag = wn_validation()
#print(wn_validation_flag)
if wn_validation_flag == 1:
    t = time.time()    
    CH4lines, tips = import_data(selected_species)
    num_of_PDF_points = N_PDF_points
    
    testing_range = True
    wn_cutoff = 0
    x_limited = np.arange(wnstart, wnend, wnres)
    x = np.arange(wnstart - wn_cutoff, wnend + wn_cutoff, wnres)
    start_x, end_x = find_range(x,CH4lines)
    start_x_limited, end_x_limited = start_x, end_x
    lines = extract_lines(start_x,end_x,CH4lines,s0_min, selected_broadener, testing_range)
    number_of_lines_limited = len(lines)
    #edited_lines = lines
    x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self, visible_start, visible_end = extract_mean_parameters(lines)
    # loop to satisfy the risidual requirement
    spectrum_mean_parameters = mean_spectrum_simulation(lines,T,P,mole_fraction,L,x,calc_method,simulation_type)
    
    residual = 1
    with tab1:
        with st.spinner('Computing lines strength threshold ...'):
            while residual > max_residual:
                s0_min = 0.1*s0_min
                lines = extract_lines(start_x,end_x,CH4lines,s0_min, selected_broadener, testing_range)
                #print(len(lines))
                #edited_lines = lines
                x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self, visible_start, visible_end = extract_mean_parameters(lines)
                extended_spectrum_mean_parameters = mean_spectrum_simulation(lines,T,P,mole_fraction,L,x,calc_method,simulation_type)
                residual = max(extended_spectrum_mean_parameters - spectrum_mean_parameters)/max(spectrum_mean_parameters)
                #('threshold residual')
                #print(residual)
                spectrum_mean_parameters = extended_spectrum_mean_parameters

    s0_min = 10*s0_min
    with st.sidebar:
        st.sidebar.info(f'Line strength threshold set at {s0_min:.1e} (cm-1/(molec.cm-2))', icon="ℹ️")
    #st.session_state.s0_min = s0_min
    #print('line strength threshold')
    #print(s0_min)
    number_of_lines_limited = len(lines)
    if number_of_lines_limited == 0:
        max_residual = 1
        st.warning('No lines within the selected frequency range.', icon="⚠️")
    else:
        with tab1:
            with st.spinner('Computing cut-off frequency ...'):
                residual = 1
                while residual > max_residual:

                    wn_cutoff = wn_cutoff + 10
                    x = np.arange(wnstart - wn_cutoff, wnend + wn_cutoff, wnres)
                    extended_start_x, extended_end_x = find_range(x,CH4lines)
                    lines = extract_lines(extended_start_x,extended_end_x,CH4lines,s0_min, selected_broadener, testing_range)
                    #edited_lines = lines
                    x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self, visible_start, visible_end = extract_mean_parameters(lines)
                    extended_spectrum_mean_parameters = mean_spectrum_simulation(lines,T,P,mole_fraction,L,x,calc_method,simulation_type)

                    extended_spectrum_mean_parameters = np.interp(x_limited, x, extended_spectrum_mean_parameters)

                    # reevaluate extended_spectrum_mean_parameters at original x_range points and store to extended_spectrum_mean_parameters
                    #residual = max(np.divide(extended_spectrum_mean_parameters - spectrum_mean_parameters,extended_spectrum_mean_parameters))
                    residual = max(extended_spectrum_mean_parameters - spectrum_mean_parameters)/max(spectrum_mean_parameters)
                    start_x = extended_start_x
                    end_x = extended_end_x
                    spectrum_mean_parameters = extended_spectrum_mean_parameters
                    #print(wn_cutoff)
                    #print(residual)

        testing_range = False
        wn_cutoff = wn_cutoff - 10

        with st.sidebar:
            st.success('Lines within '+str(wnstart - wn_cutoff)+' - '+str(wnend + wn_cutoff)+' cm-1 will be included in the simulation for enhanced accuracy within the selected wavenumber range.', icon="✅")

        x = np.arange(wnstart - wn_cutoff, wnend + wn_cutoff, wnres)
        extended_start_x, extended_end_x = find_range(x,CH4lines)
        lines = extract_lines(extended_start_x,extended_end_x,CH4lines,s0_min, selected_broadener, testing_range)
        x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self, visible_start, visible_end = extract_mean_parameters(lines)

        # Range
        #x = np.arange(wnstart - wn_cutoff, wnend + wn_cutoff, wnres)  # Similar to MATLAB's 1331:0.001:1334

        #start_x, end_x = find_range(x,CH4lines)
        #lines = extract_lines(start_x,end_x,CH4lines,s0_min, selected_broadener)

        #print(lines)
        if st.session_state.manual_control:
            update_manual_control(lines)
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
            help="Bath gas broadning coeffecient (cm-1/atm) at 296 K.",
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
            help="Temperature exponent for broadning coeffecients.",
            #min_value=0,
            #max_value=end_x,
            step=0.0001,
            format="%.4f",
        ),
        7: st.column_config.NumberColumn(
            "δ",
            help="Line position pressure shift (cm-1.atm-1).",
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
        ),
        14: st.column_config.NumberColumn(
            "δ_s",
            help="Line position pressure self shift (cm-1.atm-1).",
            #min_value=start_x,
            #max_value=end_x,
            step=0.0001,
            format="%.4f",
        ),
        15: st.column_config.NumberColumn(
            "± δ_s",
            help="Unceratinty in self pressure shift parameter (cm-1.atm-1)",
            min_value=1E-6,
            max_value=1,
            step=1E-6,
            format="%.1e",
        )
    }
            #print(lines[:])
            #print(lines[(start_x_limited-start_x):(end_x-start_x_limited)])
            #print(visible_start)
            #print(visible_end)
            edited_lines_limited = st.data_editor(lines[visible_start:visible_end],column_config=lines_column_config,key=st.session_state.dek)
            st.divider()
            #print(type(edited_lines_limited))
            edited_lines = lines[0:(visible_start)] + edited_lines_limited + lines[(visible_end):len(lines)]
        else:
            edited_lines = lines

        #print(edited_lines == lines)
        with tab1:
            with st.spinner('Extracting line parameters ...'):
                number_of_lines = len(edited_lines)
                x0_rand_cdf, x0_rand_range,S0_rand_cdf, S0_rand_range,gamma_air_rand_cdf, gamma_air_rand_range\
                ,gamma_self_rand_cdf, gamma_self_rand_range,n_air_rand_cdf, n_air_rand_range,delta_self_rand_cdf, delta_self_rand_range,\
                delta_air_rand_cdf, delta_air_rand_range, x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self = extract_parameters(edited_lines)

        #print(st.session_state.exp_unc)
        #print(exp_unc_values != [0,0,0,0])
        if (st.session_state.exp_unc & (exp_unc_values != [0,0,0,0])):
            molefraction_cdf, molefraction_range, pathlength_cdf, pathlength_range, pressure_cdf, pressure_range, temperature_cdf, temperature_range = exp_unc_distributions(mole_fraction, pathlength, pressure, temperature,exp_unc_values)

        #st.divider()

        spectra_limited = np.zeros((len(x_limited), n_simulations))    
        spectra_limited = MC_simulation(edited_lines,n_simulations,T,P,mole_fraction,L,x,exp_unc_values, calc_method,simulation_type)

        #np.savetxt('sample_results/spectra_for_analysis.csv', spectra_limited, delimiter=',')

        with tab1:
            with st.spinner('Computing spectrum based on mean parameters ...'):
                extended_spectrum_mean_parameters = mean_spectrum_simulation(edited_lines,T,P,mole_fraction,L,x,calc_method,simulation_type)
                extended_spectrum_mean_parameters = np.interp(x_limited, x, extended_spectrum_mean_parameters)
                spectrum_mean_parameters = extended_spectrum_mean_parameters

        with tab1:
            with st.spinner('Plotting simulated spectra ...'):
                fig_1 = plot_MC_spectra(spectra_limited, spectrum_mean_parameters)
        with tab1:
            st.write('_'+simulation_type+' spectrum based on mean line-parameters of '+str(number_of_lines) + ' lines,\nand '+str(n_simulations)+' spectra based on randomly sampled line-parameters:_')
            st.pyplot(fig_1)

        with tab1:
            with st.spinner('Calculating and plotting global statistics ...'):
                relative_uncertainty, error_bars, skewness = calc_error_bars(spectra_limited,spectrum_mean_parameters)
                #plotting_commands()
                fig_2 = plot_uncertainty(relative_uncertainty,skewness)
        with tab2:
            st.write('_Uncertainty spectrum (3 x std.) and skewness spectrum:_')
            st.pyplot(fig_2)    

        if conv_test == 1:
            fig_3, std_index = std_deviation_with_iterations(spectra_limited,spectrum_mean_parameters)
            with tab3:
                st.write('_Standard deviation with iterations at ('+str(round(x_limited[std_index],2))+' cm-1):_')
                st.pyplot(fig_3)

        with tab1:
            with st.spinner('Calculating and plotting PDF at ('+str(round(x_limited[std_index],2))+' cm-1)...'):
                fig_4, std_index = uncertainty_PDF(spectra_limited)
        with tab4:
            st.write('_Histogram of predicted absorbance at ('+str(round(x_limited[std_index],2))+' cm-1). Dashed line indicates predicted absorbance based on mean parameters:_')
            st.pyplot(fig_4)   

        print('Elapsed time')
        print(time.time() - t)

        simulation_info = [datetime.datetime.now(),selected_species,T,P,mole_fraction,L,wnstart,wnend,wnres,n_simulations,s0_min,st.session_state.manual_control,conv_test]
        #print(simulation_info)
        with open('simulation_history.csv','a') as fd:
            #fd.write(np.array2string(simulation_info))
            writer = csv.writer(fd)
            writer.writerow(simulation_info)

        arr = np.array(np.transpose([x_limited,spectrum_mean_parameters,error_bars,relative_uncertainty,skewness]))
        arr_df = pd.DataFrame(arr,columns=['wavenumbers (cm-1)','absorbance - mean parameters','3 x standard deviation','relative uncertainty','skewness'])
        # Create an in-memory buffer

        #st.divider()
        with tab1:
            with st.spinner('Preparing downloadable data ...'):
                with io.BytesIO() as buffer:
                    # Write array to buffer
                    #np.savetxt(buffer, arr_df, delimiter=",")
                    textstr = ('MCSpectra_'+selected_species + '_' + str(100*mole_fraction)+'_' + '_'+ str(simulation_type) + '_'+ str(T) + '_K_' + str(P) + '_atm_'+ str(L) + '_cm_' + str(wnstart) + '_'+ str(wnend)+ '_broadener_' + str(selected_broadener))
                    buffer = pd.DataFrame.to_csv(arr_df, sep=',', index=False, encoding='utf-8')

                    st.download_button(
                        label="Download results as CSV",
                        data = buffer, # Download buffer
                        file_name = textstr+'.csv',
                        mime='text/csv'
                    )
elif (wnend - wnstart) == 1990 :
    with open("simulation_history.csv", "rb") as file:
        btn = st.download_button(
            label="Download log",
            data=file,
            file_name="simulation_history.csv",
            mime="text/csv",
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