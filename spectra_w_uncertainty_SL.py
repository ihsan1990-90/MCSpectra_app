import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm, skew, median_abs_deviation
from voigtfwhm import voigtfwhm
from voigtfwhm_fast import voigtfwhm_fast
import io
import base64
import uuid
import datetime
import csv
import time
import math

# General page configuration
st.set_page_config(
    page_title="MCSpectra",
    page_icon="🗻",
    layout="centered",
    initial_sidebar_state="expanded",
    menu_items={
        'Get Help': 'https://www.kaust.edu.sa',
        'Report a bug': "mailto:ihsan.farouki@kaust.edu.sa",
        'About': "# Absorbance and emission spectra simulations with uncertainty quantification"
    }
)


# paths to logos
k_logo = "images/kaust_2.png"
MS_logo = "images/mcspectra_logo.png"
MS_logo_2 = "images/montespectra_logo_2.png"
empty_image = "images/Empty.png"

# Page setup
st.title("Simulation Results")
st.logo(empty_image,icon_image=MS_logo_2,size='large')
plt.style.use('dark_background')

# 'wn_validation()' validates the wavnumber range selected by the user
# (start of the range should be smaller than end of the range)
# (for simulations with uncertainty quantification the range is limited to 10 cm-1)
# (ub survey mode, no limit on the wavenumber range exists)
# the function returns a flag which determines whether the simulation run proceeds or not
wn_validation_flag = 1
wn_change_flag = 1
def wn_validation():
    # st.session_state.survey_mode returns whether the user has selected survey mode
    if not(st.session_state.survey_mode):
        # wnstart and wnend are the start and end wavenumber range as selected by the user
        if (wnstart >= wnend) or ((wnend - wnstart) > 10):
            st.error('Wavelength upper limit should be (up to 10 cm-1) larger than lower limit', icon="🚨")
            wn_validation_flag = 0
            wn_change_flag = 0
        else:
            wn_validation_flag = 1
            wn_change_flag = 1
    else:
        if (wnstart >= wnend):
            st.error('Wavelength upper limit should be larger than lower limit', icon="🚨")
            wn_validation_flag = 0
            wn_change_flag = 0
        else:
            wn_validation_flag = 1
            wn_change_flag = 1

    return wn_validation_flag, wn_change_flag

# 'xaxis_validation()' performs validation checks on the x-axis boundaries selected by the user
# The boundaries need to be within the range covered by the simulation, and the start of the range needs to be lower than the end of the range
def xaxis_validation():
    if ((st.session_state.xaxis_start > st.session_state.xaxis_end) or (st.session_state.xaxis_start < wnstart)  or (st.session_state.xaxis_end > wnend) or (st.session_state.xaxis_start > wnend)  or (st.session_state.xaxis_end < wnstart)) and (wn_validation_flag == 1):
        st.session_state.xaxis_start = st.session_state.wn_start
        st.session_state.xaxis_end = st.session_state.wn_end

# 'dek' is a session state variable linked to the manual line parameter control table.
# this key helps in managing the line list within the (manual control) line parameter table
# the following conditional initiates the session state variable
# later in the code, a function is provisioned to repopulate the line list when needed
if 'dek' not in st.session_state:
    st.session_state.dek = str(uuid.uuid4())

# 'change_wn_range()' adjusts the simulation wavelength range when the selected species is changed
# the simulation range is adjusted to a region which bears interesting features
def change_wn_range():
    if st.session_state.selected_species == 'CH4':
        st.session_state.wn_start = 1331
        st.session_state.wn_end = 1334
    elif st.session_state.selected_species == '(12)CH4':
        st.session_state.wn_start = 1331
        st.session_state.wn_end = 1334
    elif st.session_state.selected_species == 'H2(16)O':
        st.session_state.wn_start = 3742
        st.session_state.wn_end = 3747
    elif st.session_state.selected_species == 'CO2':
        st.session_state.wn_start = 2300
        st.session_state.wn_end = 2305
    elif st.session_state.selected_species == '(12)CO2':
        st.session_state.wn_start = 2300
        st.session_state.wn_end = 2305
    elif st.session_state.selected_species == '(13)CO2':
        st.session_state.wn_start = 2300
        st.session_state.wn_end = 2305
    elif st.session_state.selected_species == '(14)N2O':
        st.session_state.wn_start = 1285
        st.session_state.wn_end = 1290
    elif st.session_state.selected_species == 'NO':
        st.session_state.wn_start = 1810
        st.session_state.wn_end = 1820
    elif st.session_state.selected_species == '(12)CO':
        st.session_state.wn_start = 2000
        st.session_state.wn_end = 2010
    elif st.session_state.selected_species == 'CO':
        st.session_state.wn_start = 2000
        st.session_state.wn_end = 2010
    elif st.session_state.selected_species == '(14)NH3':
        st.session_state.wn_start = 850
        st.session_state.wn_end = 856
    elif st.session_state.selected_species == '(12)C2H6':
        st.session_state.wn_start = 3009
        st.session_state.wn_end = 3012
    elif st.session_state.selected_species == 'C2H6':
        st.session_state.wn_start = 3009
        st.session_state.wn_end = 3012
    elif st.session_state.selected_species == 'CH3OH':
        st.session_state.wn_start = 1030
        st.session_state.wn_end = 1035
    elif st.session_state.selected_species == 'O3':
        st.session_state.wn_start = 1010    
        st.session_state.wn_end = 1015
    elif st.session_state.selected_species == 'HF':
        st.session_state.wn_start = 3876
        st.session_state.wn_end = 3879
    elif st.session_state.selected_species == 'SO2':
        st.session_state.wn_start = 1135
        st.session_state.wn_end = 1140
    elif st.session_state.selected_species == 'NO2':
        st.session_state.wn_start = 1600
        st.session_state.wn_end = 1605
    # adjust x-axis range along with the simulation range
    st.session_state.xaxis_start = st.session_state.wn_start
    st.session_state.xaxis_end = st.session_state.wn_end

# 'change_xaxis_range()' is linked to the wavenumber range inputs
# it is called whenever the selected wavenumber range is modified
# the function adjust the x-axis range to match the selected wavenumber range
def change_xaxis_range():
    st.session_state.xaxis_start = st.session_state.wn_start
    st.session_state.xaxis_end = st.session_state.wn_end

# 'molar_mass()' assigns the molar mass of selected species to 'M'
def molar_mass():
    if st.session_state.selected_species == '(12)CH4':
        M = 16 # Molar mass of CH4 (g/mol)
    elif st.session_state.selected_species == 'CH4':
        M = 16.04 # Molar mass of CH4 (g/mol)        
    elif st.session_state.selected_species == 'H2(16)O':
        M = 18 # Molar mass of CH4 (g/mol)  
    elif st.session_state.selected_species == 'CO2':
        M = 44.01 # Molar mass of CH4 (g/mol)     
    elif st.session_state.selected_species == '(12)CO2':
        M = 44 # Molar mass of CH4 (g/mol)
    elif st.session_state.selected_species == '(13)CO2':
        M = 45 # Molar mass of CH4 (g/mol)        
    elif st.session_state.selected_species == '(14)N2O':
        M = 44 # Molar mass of CH4 (g/mol)
    elif st.session_state.selected_species == 'NO':
        M = 30.01 # Molar mass of CH4 (g/mol)       
    elif st.session_state.selected_species == '(12)CO':
        M = 28 # Molar mass of CH4 (g/mol)
    elif st.session_state.selected_species == 'CO':
        M = 28.01 # Molar mass of CH4 (g/mol)
    elif st.session_state.selected_species == '(14)NH3':
        M = 17 # Molar mass of CH4 (g/mol)
    elif st.session_state.selected_species == '(12)C2H6':
        M = 30 # Molar mass of CH4 (g/mol)
    elif st.session_state.selected_species == 'C2H6':
        M = 30 # Molar mass of CH4 (g/mol)
    elif st.session_state.selected_species == 'O3':
        M = 48 # Molar mass of CH4 (g/mol)
    elif st.session_state.selected_species == 'HF':
        M = 20 # Molar mass of CH4 (g/mol)
    elif st.session_state.selected_species == 'CH3OH':
        M = 32.04 # Molar mass of CH4 (g/mol)
    elif st.session_state.selected_species == 'SO2':
        M = 64.066 # Molar mass of CH4 (g/mol)
    elif st.session_state.selected_species == 'NO2':
        M = 46.0055 # Molar mass of CH4 (g/mol)     
    
    return M

# list of species for which species are available in the /HITRAN_data        
species_options = ['CH4', '(12)CH4', 'H2(16)O', 'CO2', '(12)CO2', '(13)CO2', '(14)N2O', 'NO','(12)CO','CO','(14)NH3','(12)C2H6','C2H6','O3','HF','CH3OH','SO2','NO2']

# pre-programmed list of broadeners for different species
# also indicates whether self-shift parameter data is available 
# based on the availablility of data in the HITRAN database
if not(hasattr(st.session_state,'selected_species')):
    broadener_options = ['Air','H2O']
    self_shift_available = False
elif st.session_state.selected_species == 'CH4':
    broadener_options = ['Air','H2O']
    self_shift_available = False
elif st.session_state.selected_species == '(12)CH4':
    broadener_options = ['Air','H2O']
    self_shift_available = False
elif st.session_state.selected_species == 'H2(16)O':
    broadener_options = ['Air']
    self_shift_available = False
elif st.session_state.selected_species == 'CO2':
    broadener_options = ['Air','H2','He','H2O']
    self_shift_available = True
elif st.session_state.selected_species == '(12)CO2':
    broadener_options = ['Air','H2','He','H2O']
    self_shift_available = True
elif st.session_state.selected_species == '(13)CO2':
    broadener_options = ['Air','H2','He','H2O']
    self_shift_available = True
elif st.session_state.selected_species == '(14)N2O':
    broadener_options = ['Air','He','H2O']
    self_shift_available = True
elif st.session_state.selected_species == 'NO':
    broadener_options = ['Air']
    self_shift_available = False
elif st.session_state.selected_species == '(12)CO':
    broadener_options = ['Air','H2','He','CO2','H2O']
    self_shift_available = True
elif st.session_state.selected_species == 'CO':
    broadener_options = ['Air','H2','He','CO2','H2O']
    self_shift_available = True
elif st.session_state.selected_species == '(14)NH3':
    broadener_options = ['Air','H2','He','CO2','H2O']
    self_shift_available = False
elif st.session_state.selected_species == '(12)C2H6':
    broadener_options = ['Air']
    self_shift_available = False
elif st.session_state.selected_species == 'C2H6':
    broadener_options = ['Air']
    self_shift_available = False
elif st.session_state.selected_species == 'O3':
    broadener_options = ['Air']
    self_shift_available = False
elif st.session_state.selected_species == 'HF':
    broadener_options = ['Air','H2','He']
    self_shift_available = False
elif st.session_state.selected_species == 'CH3OH':
    broadener_options = ['Air']
    self_shift_available = False
elif st.session_state.selected_species == 'SO2':
    broadener_options = ['Air','H2','He']
    self_shift_available = False
elif st.session_state.selected_species == 'NO2':
    broadener_options = ['Air']
    self_shift_available = False

with open('simulation_history.csv') as f:
    row_count = sum(1 for line in f)

#history = pd.read_csv('simulation_history.csv', on_bad_lines='skip')
# the following section of the code contains the components listed in the sidebar
# mainly consist of number inputs and toggles to control the simulation
with st.sidebar:
    st.image(MS_logo, width=250)
    # 'Total number of simulations: '+str(row_count/2)
    st.text(f'Total # of simulations: {(row_count/2):.2e}')

    # basic simulation controls within a collapsable container
    with st.expander('Basic simulation controls',True):
        survey_mode = st.toggle("Survey mode", help='Simulate spectra over a broad wavenumber range without uncertainty quantification.', key='survey_mode')
        simulation_type = st.selectbox("Spectrum type", ['Absorbance', 'Transmittance','Emission'], 0,key='simulation_type')
        selected_species = st.selectbox("Species", species_options, 0, on_change=change_wn_range,key='selected_species')
        selected_broadener = st.selectbox("Bath-gas", broadener_options, 0,key='selected_broadener')
        temperature = st.number_input("Temperature (K)", min_value=300, max_value=3000, value=300, step=100)
        pressure = st.number_input("Pressure (atm)", min_value=0.001, max_value=100.00, value=1.00, step=0.2)
        molefraction = st.number_input("Mole Fraction", min_value=0.00, max_value=1.00, value=0.01, step=0.001, format="%.3e")
        pathlength = st.number_input('Pathlength (cm)', min_value=1, max_value=50000, step=1, value=10)
        wnstart = st.number_input('Wavelength start (cm-1)', min_value=500.00, max_value=5000.00, step=1.00, value=1331.00, on_change=change_xaxis_range, key='wn_start')
        wnend = st.number_input('Wavelength end (cm-1)', min_value=500.00, max_value=5000.00, step=1.00, value=1334.00, on_change=change_xaxis_range, key='wn_end')
        wnres = st.number_input('Resolution (cm-1)', min_value=0.001, max_value=0.1, step=0.001, value=0.005, key='wn_res', format="%.3f")
    # plotting controls
    # useful for zooming in and out on the x-axis
    # also to select the wavenumber at which the standard deviation convergence is tested
    with st.expander('Plotting controls',True):
        if 'xaxis_end' not in st.session_state:
            st.session_state.xaxis_end = wnend
        xaxis_start = st.number_input('x-axis min (cm-1)', help="Use x-axis min and max to zoom into portions of simulated spectra without having to rerun the simulations." , step=1.0, value=1331.00, on_change=xaxis_validation, key='xaxis_start')
        xaxis_end = st.number_input('x-axis max (cm-1)', step=1.0, value=1334.00, on_change=xaxis_validation, key='xaxis_end')
        conv_test = 1
        if conv_test == 1:
            convergence_frequency = st.number_input('Cursor position for convergence test and PDF (cm-1)', step=1.0, min_value=float(min(st.session_state.wn_start,st.session_state.wn_end)), max_value=float(max(st.session_state.wn_start,st.session_state.wn_end)), value=float(min(st.session_state.wn_start,st.session_state.wn_end)), key='wn_conv')
    
    # Advanced simulation controls
    with st.expander('Advanced simulation controls',True):
        # number of MC spectra instances adjustable only when not in survey mode 
        if not(st.session_state.survey_mode):
            N_simulations = st.number_input('Number of simulations', min_value=1, max_value=2000, step=100, value=1000)
        max_residual = st.number_input("Max. allowed risidual due to frequency cut-off or line-strength threshold", help='Risidual calculated as a fraction of the maximum absorbance/emission within the selected wavenumber range.', min_value=1E-4, max_value=1E-2, value=1E-2, format="%.3e")
        manual_control = st.toggle("Enable manual control of line parameters",key='manual_control')
        calc_method_wofz = st.toggle("More accurate evaluation of the Voigt function",key='calc_method_wofz')
        # manually adjusting line parameters is possible only when not in survey mode 
        if not(st.session_state.survey_mode):
            exp_unc = st.toggle("Account for uncertainty in experimental conditions",key='exp_unc')
            if st.session_state.exp_unc:
                molefraction_unc = st.number_input('Uncertainty in mole fraction (%)', min_value=0, max_value=100, value=0, key='molefraction_unc')
                pathlength_unc = st.number_input('Uncertainty in pathlength (%)', min_value=0, max_value=100, value=0, key='pathlength_unc')
                pressure_unc = st.number_input('Uncertainty in pressure (%)', min_value=0, max_value=100, value=0, key='pressure_unc')
                temperature_unc = st.number_input('Uncertainty in temperature (%)', min_value=0, max_value=100, value=0, key='temperature_unc')
        else:
            st.session_state.exp_unc = 0
    
    # warning and information messages are displayed below this subheader in the sidebar
    st.subheader('Warnings')

tab1, tab2, tab3, tab4 = st.tabs(["Simulated spectra","Global statistics","Convergence","PDF"])

# Constants
h = 6.626070E-34  # Planck's constant (J.s)
kb = 1.380649E-23  # Boltzmann constant (J/K)
c = 2.99792458E+10  # Speed of light (cm/s)
R = 1.36259479E-22  # Gas constant (cm^3.atm.K^-1.molecule^-1)

# %% Initial parameters (replace MATLAB's "clear")
if st.session_state.survey_mode:
    N_simulations = 1
    n_simulations = N_simulations
else:
    n_simulations = N_simulations


T = temperature  # Temperature in K
mole_fraction = molefraction
P = pressure  # Pressure in atm
L = pathlength  # Path length in cm
# S0_min value is automatically adjusted depending on the selected wavelength range
s0_min = 1E-21 # initial minimum line strength (cm-1/(molecule.cm-2)), 
M = molar_mass()  # Molar mass of selected species g/mol

# set uncertainty in experimental physical variables
exp_unc_values = [0,0,0,0]
if st.session_state.exp_unc:
    exp_unc_values = [molefraction_unc,pathlength_unc,pressure_unc,temperature_unc]

#np.set_printoptions(precision=16)
#np.set_printoptions(legacy='1.25')


# import HITRAN and TIPS (total internal partition sums) data for the selected species
@st.cache_resource(show_spinner=True,max_entries=3)
def import_data(selected_species):
    # print('importing data')

    if selected_species == '(12)CH4':
        selected_species_lines = pd.read_csv('HITRAN_data/12CH4_lines_formatted.csv').values
        tips = pd.read_csv('HITRAN_data/q32_12CH4.csv', sep='\s+').values
        num_of_isotopologues = 1
        first_isotopologue = 0
        isotopologue_abundance = 0.988274
        rotational_constant = 5.24 #cm-1
    elif selected_species == 'CH4':
        selected_species_lines = pd.read_csv('HITRAN_data/CH4_natural_lines_formatted.csv').values
        tips = np.genfromtxt('HITRAN_data/q_CH4_natural.csv', delimiter=',')
        num_of_isotopologues = 2
        first_isotopologue = 32
        isotopologue_abundance = 1
        rotational_constant = 5.24 #cm-1
    elif selected_species == 'H2(16)O':
        selected_species_lines = pd.read_csv('HITRAN_data/H216O_lines_formatted.csv').values
        tips = pd.read_csv('HITRAN_data/q1_H2O16.csv', sep='\s+').values
        num_of_isotopologues = 1
        first_isotopologue = 0
        isotopologue_abundance = 0.997317
        rotational_constant = 9.28 #cm-1
    elif selected_species == 'CO2':
        selected_species_lines = pd.read_csv('HITRAN_data/CO2_natural_lines_formatted.csv').values
        #tips = pd.read_csv('HITRAN_data/q_CO2_natural.csv', sep='\s+').values
        tips = np.genfromtxt('HITRAN_data/q_CO2_natural.csv', delimiter=',')
        num_of_isotopologues = 3
        first_isotopologue = 7
        isotopologue_abundance = 1
        rotational_constant = 0.39 #cm-1
    elif selected_species == '(12)CO2':
        selected_species_lines = pd.read_csv('HITRAN_data/12CO2_lines_formatted.csv').values
        tips = pd.read_csv('HITRAN_data/q7_12CO2.csv', sep='\s+').values
        num_of_isotopologues = 1
        first_isotopologue = 0
        isotopologue_abundance = 0.984204
        rotational_constant = 0.39 #cm-1
    elif selected_species == '(13)CO2':
        selected_species_lines = pd.read_csv('HITRAN_data/13CO2_lines_formatted.csv').values
        tips = pd.read_csv('HITRAN_data/q8_13CO2.csv', sep='\s+').values
        num_of_isotopologues = 1
        first_isotopologue = 0
        isotopologue_abundance = 0.0110574
        rotational_constant = 0.39 #cm-1
    elif selected_species == '(14)N2O':
        selected_species_lines = pd.read_csv('HITRAN_data/14N2O_lines_formatted.csv').values
        tips = pd.read_csv('HITRAN_data/q21_14N2O.csv', sep='\s+').values
        num_of_isotopologues = 1
        first_isotopologue = 0
        isotopologue_abundance = 0.990333
        rotational_constant = 0.42 #cm-1
    elif selected_species == 'NO':
        selected_species_lines = pd.read_csv('HITRAN_data/NO_natural_lines_formatted.csv').values
        #tips = pd.read_csv('HITRAN_data/q_CO2_natural.csv', sep='\s+').values
        tips = np.genfromtxt('HITRAN_data/q_NO_natural.csv', delimiter=',')
        num_of_isotopologues = 3
        first_isotopologue = 39
        isotopologue_abundance = 1
        rotational_constant = 1.7 #cm-1
    elif selected_species == '(12)CO':
        selected_species_lines = pd.read_csv('HITRAN_data/12CO_lines_formatted.csv').values
        tips = pd.read_csv('HITRAN_data/q26_12CO.csv', sep='\s+').values
        num_of_isotopologues = 1
        first_isotopologue = 0
        isotopologue_abundance = 0.986544
        rotational_constant = 1.93 #cm-1
    elif selected_species == 'CO':
        selected_species_lines = pd.read_csv('HITRAN_data/CO_natural_lines_formatted.csv').values
        tips = np.genfromtxt('HITRAN_data/q_CO_natural.csv', delimiter=',')
        num_of_isotopologues = 3
        first_isotopologue = 26
        isotopologue_abundance = 1
        rotational_constant = 1.93 #cm-1
    elif selected_species == '(14)NH3':
        selected_species_lines = pd.read_csv('HITRAN_data/14NH3_lines_formatted.csv').values
        tips = pd.read_csv('HITRAN_data/q45_14NH3.csv', sep='\s+').values
        num_of_isotopologues = 1
        first_isotopologue = 0
        isotopologue_abundance = 0.995872
        rotational_constant = 9.93 #cm-1
    elif selected_species == '(12)C2H6':
        selected_species_lines = pd.read_csv('HITRAN_data/12C2H6_lines_formatted.csv').values
        tips = pd.read_csv('HITRAN_data/q78_12C2H6.csv', sep='\s+').values
        num_of_isotopologues = 1
        first_isotopologue = 0
        isotopologue_abundance = 0.976990
        rotational_constant = 2.86 #cm-1
    elif selected_species == 'C2H6':
        selected_species_lines = pd.read_csv('HITRAN_data/12C2H6_lines_formatted.csv').values
        tips = pd.read_csv('HITRAN_data/q78_12C2H6.csv', sep='\s+').values
        num_of_isotopologues = 1
        first_isotopologue = 0
        isotopologue_abundance = 1
        rotational_constant = 2.86 #cm-1
    elif selected_species == 'O3':
        selected_species_lines = pd.read_csv('HITRAN_data/O3_natural_lines_formatted.csv').values
        #tips = pd.read_csv('HITRAN_data/q_CO2_natural.csv', sep='\s+').values
        tips = np.genfromtxt('HITRAN_data/q_O3_natural.csv', delimiter=',')
        num_of_isotopologues = 2
        first_isotopologue = 16
        isotopologue_abundance = 1
        rotational_constant = 0.45 #cm-1
    elif selected_species == 'HF':
        selected_species_lines = pd.read_csv('HITRAN_data/HF_natural_lines_formatted.csv').values
        tips = np.genfromtxt('HITRAN_data/q_HF_natural.csv', delimiter=',')
        num_of_isotopologues = 1
        first_isotopologue = 0
        isotopologue_abundance = 1
        rotational_constant = 20.96 #cm-1
    elif selected_species == 'CH3OH':
        selected_species_lines = pd.read_csv('HITRAN_data/CH3OH_natural_lines_formatted.csv').values
        tips = np.genfromtxt('HITRAN_data/q_CH3OH_natural.csv', delimiter=',')
        num_of_isotopologues = 1
        first_isotopologue = 0
        isotopologue_abundance = 1
        rotational_constant = 4.26 #cm-1
    elif selected_species == 'SO2':
        selected_species_lines = pd.read_csv('HITRAN_data/SO2_natural_lines_formatted.csv').values
        tips = np.genfromtxt('HITRAN_data/q_SO2_natural.csv', delimiter=',')
        num_of_isotopologues = 3
        first_isotopologue = 42
        isotopologue_abundance = 1
        rotational_constant = 0.001451 #cm-1
    elif selected_species == 'NO2':
        selected_species_lines = pd.read_csv('HITRAN_data/NO2_natural_lines_formatted.csv').values
        tips = np.genfromtxt('HITRAN_data/q_NO2_natural.csv', delimiter=',')
        num_of_isotopologues = 2
        first_isotopologue = 44
        isotopologue_abundance = 1
        rotational_constant = 8.0012 #cm-1
    
    return selected_species_lines, tips, num_of_isotopologues, first_isotopologue, isotopologue_abundance, rotational_constant


# Define interpolation function for TIPS data
# allows calculating the partetion function ratio at any temperature, not just the points listed in the table
def tips1(z,tips_index,tips):
    return np.interp(z, tips[:, 0], tips[:, int(tips_index+1)])

# Blackbody radiance formula, function of wavenumber and temperature
def plank_emission(x,T):
    Ibb = 2*h*(c**2)*np.divide(np.power(x,3),(np.exp(c*h*x/(kb*T))-1))
    return Ibb



# start and end indices for the lines going into the calculation
# find the index range over which the lines of interest exist in the imported data
@st.cache_resource(show_spinner=False,max_entries=3)
def find_range(x,selected_species_lines):
    #print('finding start and end indices')
    #start_x = np.where(np.round(selected_species_lines[:, 0],0) == np.round(min(x),0))[0][0]
    
    #print('where test')
    #print(np.where(np.abs(min(x) - selected_species_lines[:, 0]) == min(np.abs(min(x) - selected_species_lines[:, 0])))[0])
    #start_x = np.where(((min(x) - selected_species_lines[:, 0]) < 10) & ((min(x) - selected_species_lines[:, 0]) >= 0))[0][-1]
    start_x = np.where(np.abs(min(x) - selected_species_lines[:, 0]) == min(np.abs(min(x) - selected_species_lines[:, 0])))[0][0]

    #print(np.where(np.abs(selected_species_lines[:, 0] - max(x)) == min(np.abs(selected_species_lines[:, 0] - max(x))))[0] )
    #end_x = np.where(((selected_species_lines[:, 0] - max(x)) < 10) & ((selected_species_lines[:, 0] - max(x)) > 0) )[0][0]
    end_x = np.where(np.abs(selected_species_lines[:, 0] - max(x)) == min(np.abs(selected_species_lines[:, 0] - max(x))))[0][0]
    


    return start_x, end_x

# Extract the line parameters and their uncertainties for lines within the range
@st.cache_resource(show_spinner=False,max_entries=3)
def extract_lines(start_x,end_x,selected_species_lines,s0_min,selected_broadener, testing_range,isotopologue_abundance):
    #print('extracting lines')
    flags_array = np.zeros(23)
    lines = []
    j = 0
    
    # loops over all the lines within the range of interest
    # extracts parameters and arranges them into 'lines' which is returned
    for i in range(start_x,end_x): #range(len(selected_species_lines)):
        if selected_species_lines[i,1] > s0_min: #(selected_species_lines[i, 0] > min(x)) and (selected_species_lines[i, 0] < max(x)):
            #print(selected_species_lines[i])
            j += 1
            
            # the following conditional statement extracts seven parameters
            # the extraction indices are different depending on the type of bath-gas
            if selected_broadener == 'Air':
                line = [
                    selected_species_lines[i, 0],  # line position
                    selected_species_lines[i, 1]/isotopologue_abundance,  # line strength
                    selected_species_lines[i, 2],  # gamma air
                    selected_species_lines[i, 3],  # gamma self
                    selected_species_lines[i, 4],  # LES
                    selected_species_lines[i, 5],  # n_air
                    selected_species_lines[i, 6]   # delta_air
                ]
            elif selected_broadener == 'H2':
                line = [
                    selected_species_lines[i, 0],  # line position
                    selected_species_lines[i, 1]/isotopologue_abundance,  # line strength
                    selected_species_lines[i, 11],  # gamma H2
                    selected_species_lines[i, 3],  # gamma self
                    selected_species_lines[i, 4],  # LES
                    selected_species_lines[i, 12],  # n_H2
                    selected_species_lines[i, 13]   # delta_H2
                ]
                # following conditional detects and warns whether temperature dependent exponent or bath-gas pressure shift parameters are missing 
                # 'testing_range' and 'flags_array' are intended to prevent displaying redundunt warning messages in the sidebar
                if str(line[5]) == 'nan':
                    line[5] = 0
                    if (not(testing_range)) & (flags_array[1] == 0):
                        with st.sidebar:
                            flags_array[1] = 1
                            st.sidebar.warning('Bath-gas temperature exponent is unavailable for '+str(selected_species_lines[i, 0])+'. Assumed to be equal to zero.', icon="⚠️")
                if str(line[6]) == 'nan':
                    line[6] = 0
                    if (not(testing_range)) & (flags_array[2] == 0):
                        with st.sidebar:
                            flags_array[2] =1
                            st.sidebar.warning('Bath-gas pressure shift parameter is unavailable for '+str(selected_species_lines[i, 0])+'. Assumed to be equal to zero.', icon="⚠️")
            
            elif selected_broadener == 'He':
                line = [
                    selected_species_lines[i, 0],  # line position
                    selected_species_lines[i, 1]/isotopologue_abundance,  # line strength
                    selected_species_lines[i, 15],  # gamma He
                    selected_species_lines[i, 3],  # gamma self
                    selected_species_lines[i, 4],  # LES
                    selected_species_lines[i, 16],  # n_He
                    selected_species_lines[i, 17]   # delta_He
                ]
                # following conditional detects and warns whether temperature dependent exponent or bath-gas pressure shift parameters are missing
                # 'testing_range' and 'flags_array' are intended to prevent displaying redundunt warning messages in the sidebar 
                if str(line[5]) == 'nan':
                    line[5] = 0
                    if (not(testing_range)) & (flags_array[3] == 0):
                        with st.sidebar:
                            flags_array[3] =1
                            st.sidebar.warning('Bath-gas temperature exponent is unavailable for '+str(selected_species_lines[i, 0])+'. Assumed to be equal to zero.', icon="⚠️")
                if str(line[6]) == 'nan':
                    line[6] = 0
                    if (not(testing_range)) & (flags_array[4] == 0):
                        with st.sidebar:
                            flags_array[4] =1
                            st.sidebar.warning('Bath-gas pressure shift parameter is unavailable for '+str(selected_species_lines[i, 0])+'. Assumed to be equal to zero.', icon="⚠️")
            
            elif selected_broadener == 'CO2':
                line = [
                    selected_species_lines[i, 0],  # line position
                    selected_species_lines[i, 1]/isotopologue_abundance,  # line strength
                    selected_species_lines[i, 18],  # gamma CO2
                    selected_species_lines[i, 3],  # gamma self
                    selected_species_lines[i, 4],  # LES
                    selected_species_lines[i, 19],  # n_CO2
                    selected_species_lines[i, 20]   # delta_CO2
                ]
                if str(line[5]) == 'nan':
                    line[5] = 0
                    if (not(testing_range)) & (flags_array[5] == 0):
                        with st.sidebar:
                            flags_array[5] =1
                            st.sidebar.warning('Bath-gas temperature exponent is unavailable for '+str(selected_species_lines[i, 0])+'. Assumed to be equal to zero.', icon="⚠️")
                if str(line[6]) == 'nan':
                    line[6] = 0
                    if (not(testing_range)) & (flags_array[6] == 0):
                        with st.sidebar:
                            flags_array[6] =1
                            st.sidebar.warning('Bath-gas pressure shift parameter is unavailable for '+str(selected_species_lines[i, 0])+'. Assumed to be equal to zero.', icon="⚠️")
            elif selected_broadener == 'H2O':
                line = [
                    selected_species_lines[i, 0],  # line position
                    selected_species_lines[i, 1]/isotopologue_abundance,  # line strength
                    selected_species_lines[i, 21],  # gamma H2O
                    selected_species_lines[i, 3],  # gamma self
                    selected_species_lines[i, 4],  # LES
                    selected_species_lines[i, 22],  # n_H2O
                    0   # delta_He
                ]
                # following conditional detects and warns whether temperature dependent exponent or bath-gas pressure shift parameters are missing
                # 'testing_range' and 'flags_array' are intended to prevent displaying redundunt warning messages in the sidebar 
                if str(line[5]) == 'nan':
                    line[5] = 0
                    if (not(testing_range)) & (flags_array[7] == 0):
                        with st.sidebar:
                            flags_array[7] =1
                            st.sidebar.warning('Bath-gas temperature exponent is unavailable for '+str(selected_species_lines[i, 0])+'. Assumed to be equal to zero.', icon="⚠️")
                
                if (not(testing_range)):
                    with st.sidebar:
                        st.sidebar.warning('Bath-gas pressure shift parameter is unavailable for '+str(selected_species_lines[i, 0])+'. Assumed to be equal to zero.', icon="⚠️")

            # The following conditionals read uncertainty parameters and convert them to relative/absolute uncertainties to be used later on
            # extraction indices are different depending on the type of bath-gas
            # the following 'for' statement extract uncertainty in the line position
            for k in [1]:
                uncertainty = selected_species_lines[i, 22 + k]
                if uncertainty == 0:
                    if (not(testing_range)) & (flags_array[8] == 0):
                        with st.sidebar:
                            flags_array[8] =1
                            st.sidebar.warning('Uncertainty in line position is unavailable for '+str(selected_species_lines[i, 0]), icon="⚠️")
                    line.append(0)
                elif uncertainty == 1:
                    if (not(testing_range)) & (flags_array[9] == 0):
                        with st.sidebar:
                            flags_array[9] =1
                            st.sidebar.warning('Uncertainty in line position (parameter '+str(k)+') for '+str(np.round(selected_species_lines[i, 0],2)) + ' might be larger than visible in the simulation result', icon="⚠️")
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
                    line.append(1E-6)
                elif uncertainty == 8:
                    line.append(1E-7)
            # the following 'for' statement extract uncertainty in line strength
            for k in [2]:
                uncertainty = selected_species_lines[i, 22 + k]
                if uncertainty in [0, 1, 2]:
                    if (not(testing_range))  & (flags_array[10] == 0):
                        with st.sidebar:
                            flags_array[10] =1
                            st.sidebar.warning('Uncertainty in line strength is unavailable for '+str(np.round(selected_species_lines[i, 0],2)), icon="⚠️")
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
                            st.sidebar.warning('Uncertainty in line strength for '+str(np.round(selected_species_lines[i, 0],2)) + ' might be larger than visible in the simulation result', icon="⚠️")
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
            
            # the following 'for' statement extracts uncertainty in bath-gas broadening coeffecient
            for k in [temp_index]:
                uncertainty = selected_species_lines[i, 22 + k]
                if uncertainty in [0, 1, 2]:
                    if (not(testing_range)) & (flags_array[12] == 0):
                        with st.sidebar:
                            flags_array[12] =1
                            st.sidebar.warning('Uncertainty in bath-gas broadening coeffecient is unavailable for '+str(np.round(selected_species_lines[i, 0],2)), icon="⚠️")
                    line.append((100)*0)
                elif uncertainty == 3:
                    if (not(testing_range)) & (flags_array[13] == 0):
                        with st.sidebar:
                            flags_array[13] =1
                            st.sidebar.warning('Uncertainty in bath-gas broadening coeffecient for '+str(np.round(selected_species_lines[i, 0],2)) + ' might be larger than visible in the simulation result', icon="⚠️")
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

            # the following 'for' statement extract uncertainty in the line position
            for k in [4]:
                uncertainty = selected_species_lines[i, 22 + k]
                if uncertainty in [0, 1, 2]:
                    if (not(testing_range)) & (flags_array[14] == 0):
                        with st.sidebar:
                            flags_array[14] =1
                            st.sidebar.warning('Uncertainty in self-broadening coeffecient is unavailable for '+str(np.round(selected_species_lines[i, 0],2)), icon="⚠️")
                    line.append((100)*0)
                elif uncertainty == 3:
                    if (not(testing_range)) & (flags_array[15] == 0):
                        with st.sidebar:
                            flags_array[15] =1
                            st.sidebar.warning('Uncertainty in self-broadening coeffecient is unavailable for '+str(np.round(selected_species_lines[i, 0],2)) + ' might be larger than visible in the simulation result', icon="⚠️")
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

            # the following 'for' statement extract uncertainty in the temperature exponent of bath-gas broadening coefficient
            for k in [temp_index]:
                uncertainty = selected_species_lines[i, 22 + k]
                if uncertainty in [0, 1, 2]:
                    if (not(testing_range)) & (flags_array[16] == 0):
                        with st.sidebar:
                            flags_array[16] =1
                            st.sidebar.warning('Uncertainty in temperature exponent of bath-gas broadening coeffecient is unavailable for '+str(np.round(selected_species_lines[i, 0],2)), icon="⚠️")
                    line.append((100)*0)
                elif uncertainty == 3:
                    if (not(testing_range)) & (flags_array[17] == 0):
                        with st.sidebar:
                            flags_array[17] =1
                            st.sidebar.warning('Uncertainty in temperature exponent of bath-gas broadening coeffecient for '+str(np.round(selected_species_lines[i, 0],2)) + ' might be larger than visible in the simulation result', icon="⚠️")
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
            
            # the following 'for' statement extract uncertainty in pressure shift due to bath-gas
            for k in [temp_index]:
                uncertainty = selected_species_lines[i, 22 + k]
                if uncertainty == 0:
                    if (not(testing_range)) & (flags_array[18] == 0):
                        with st.sidebar:
                            flags_array[18] = 1
                            st.sidebar.warning('Uncertainty in bath-gas pressure shift is unavailable for '+str(selected_species_lines[i, 0]), icon="⚠️")
                    line.append(0)
                elif uncertainty == 1:
                    if (not(testing_range)) & (flags_array[19] == 0):
                        with st.sidebar:
                            flags_array[19] = 1
                            st.sidebar.warning('Uncertainty in bath-gas pressure shift for '+str(np.round(selected_species_lines[i, 0],2)) + ' might be larger than visible in the simulation result', icon="⚠️")
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
            
            # the following 'for' statement extract uncertainty in pressure shift due absorbing/emitting species
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
                line.append(selected_species_lines[i, 8])
                for k in [8]:
                    uncertainty = selected_species_lines[i, 22 + k]
                    if uncertainty == 0:
                        if (not(testing_range)) & (flags_array[22] == 0):
                            with st.sidebar:
                                flags_array[22] = 1
                                st.sidebar.warning('Uncertainty in self pressure shift is unavailable for '+str(selected_species_lines[i, 0]), icon="⚠️")
                        line.append(0)
                    elif uncertainty == 1:
                        if (not(testing_range)) & (flags_array[23] == 0):
                            flags_array[23] = 1
                            with st.sidebar:
                                st.sidebar.warning('Uncertainty in self pressure shift for '+str(np.round(selected_species_lines[i, 0],2)) + ' might be larger than visible in the simulation result', icon="⚠️")
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

            #isotopologue ID
            line.append(selected_species_lines[i, 45])
            #print(line)
            lines.append(line)
    #st.text('Number of simulated lines: '+str(len(lines)))
    #my_spinner.empty()
    return lines

# the implementation which utilized the following two functions have been replaced by np.random.Generator.normal
# this implementation is more general and maybe restored if parameters' possible values do not follow a normal distribution

#def rand_distribution(mu, sigma):
#    range_vals = np.linspace(mu - 3 * sigma, mu + 3 * sigma, num_of_PDF_points)
#    cdf_vals = norm.cdf(range_vals, loc=mu, scale=sigma)
#    return cdf_vals, range_vals
#
#def random_value(cdf, range_vals):
#    return np.interp(np.random.rand(), cdf, range_vals)

# extract parameter values and uncertainties from the 'lines' array generated in 'extract_lines()' in preparation for spectra calculations
# The function also used to generate distributions for each line parameter in an earlier version
# lines which are not needed in the new implementation are commented out
@st.cache_resource(show_spinner=False,max_entries=3)
def extract_parameters(lines):
    print('extracting parameters and generating distributions')
    #x0_rand_cdf = np.zeros((len(lines),num_of_PDF_points))
    #x0_rand_range = np.zeros((len(lines),num_of_PDF_points))
    #S0_rand_cdf = np.zeros((len(lines),num_of_PDF_points))
    #S0_rand_range = np.zeros((len(lines),num_of_PDF_points))
    #gamma_air_rand_cdf = np.zeros((len(lines),num_of_PDF_points))
    #gamma_air_rand_range = np.zeros((len(lines),num_of_PDF_points))
    #gamma_self_rand_cdf = np.zeros((len(lines),num_of_PDF_points))
    #gamma_self_rand_range = np.zeros((len(lines),num_of_PDF_points))
    #n_air_rand_cdf = np.zeros((len(lines),num_of_PDF_points))
    #n_air_rand_range = np.zeros((len(lines),num_of_PDF_points))
    #delta_air_rand_cdf = np.zeros((len(lines),num_of_PDF_points))
    #delta_air_rand_range = np.zeros((len(lines),num_of_PDF_points))
    #delta_self_rand_cdf = np.zeros((len(lines),num_of_PDF_points))
    #delta_self_rand_range = np.zeros((len(lines),num_of_PDF_points))

    # initialize uncertainties arrays
    x0_sigma = np.zeros(len(lines))
    s0_sigma = np.zeros(len(lines))
    gamma_air_sigma = np.zeros(len(lines))
    gamma_self_sigma = np.zeros(len(lines))
    n_air_sigma = np.zeros(len(lines))
    delta_air_sigma = np.zeros(len(lines))
    delta_self_sigma = np.zeros(len(lines))
    # initialize parameters arrays
    x0 = np.zeros(len(lines))
    s0 = np.zeros(len(lines))
    gamma_air_0 = np.zeros(len(lines))
    gamma_self_0 = np.zeros(len(lines))
    n_air = np.zeros(len(lines))
    delta_air = np.zeros(len(lines))
    delta_self = np.zeros(len(lines))
    isotopologue_ID = np.zeros(len(lines))

    j=0
    for line in lines:
        #print(line)

        #x0_rand_cdf[j], x0_rand_range[j] = rand_distribution(line[0], line[7] / 3)
        #S0_rand_cdf[j], S0_rand_range[j] = rand_distribution(line[1], line[1] * (1/100)*line[8] / 3)
        #gamma_air_rand_cdf[j], gamma_air_rand_range[j] = rand_distribution(line[2], line[2] * (1/100)*line[9] / 3)
        #gamma_self_rand_cdf[j], gamma_self_rand_range[j] = rand_distribution(line[3], line[3] * (1/100)*line[10] / 3)
        #n_air_rand_cdf[j], n_air_rand_range[j] = rand_distribution(line[5], line[5] * (1/100)*line[11] / 3)
        #delta_air_rand_cdf[j], delta_air_rand_range[j] = rand_distribution(line[6], line[12] / 3)
        #delta_self_rand_cdf[j], delta_self_rand_range[j] = rand_distribution(line[13], line[14] / 3)
        
        # uncertainties
        x0_sigma[j] = line[7] / 3
        s0_sigma[j] = line[1] * (1/100)*line[8] / 3
        gamma_air_sigma[j] = line[2] * (1/100)*line[9] / 3
        gamma_self_sigma[j] = line[3] * (1/100)*line[10] / 3
        n_air_sigma[j] = line[5] * (1/100)*line[11] / 3
        delta_air_sigma[j] = line[12] / 3
        delta_self_sigma[j] = line[14] / 3

        # line parameter values
        x0[j] = line[0]
        s0[j] = line[1]
        gamma_air_0[j] = line[2]
        gamma_self_0[j] = line[3]
        n_air[j] = line[5]
        delta_air[j] = line[6]

        if  np.isnan(line[13]):
            delta_self[j] = 0
        else:
            delta_self[j] = line[13]
        
        isotopologue_ID[j] = line[15]

        j = j + 1

    #return  x0_rand_cdf, x0_rand_range,S0_rand_cdf, S0_rand_range,gamma_air_rand_cdf, gamma_air_rand_range\
    #,gamma_self_rand_cdf, gamma_self_rand_range,n_air_rand_cdf, n_air_rand_range, delta_self_rand_cdf, delta_self_rand_range,\
    #delta_air_rand_cdf, delta_air_rand_range, x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self, isotopologue_ID

    return  x0_sigma, s0_sigma, gamma_air_sigma,gamma_self_sigma, n_air_sigma, delta_self_sigma,delta_air_sigma, x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self, isotopologue_ID

# extract parameter values only without uncertainties from the 'lines' array generated in 'extract_lines()'
# has additional measures needed for the cut-off frequency determination (visible start and visible end)
@st.cache_resource(show_spinner=False,max_entries=3)
def extract_mean_parameters(lines):
    print('extracting mean parameters')

    # initialize parameters arrays
    x0 = np.zeros(len(lines))
    s0 = np.zeros(len(lines))
    gamma_air_0 = np.zeros(len(lines))
    gamma_self_0 = np.zeros(len(lines))
    n_air = np.zeros(len(lines))
    delta_air = np.zeros(len(lines))
    delta_self = np.zeros(len(lines))
    isotopologue_ID = np.zeros(len(lines))


    j=0
    temp_flag_start = False
    temp_flag_end = False
    visible_start = 0
    visible_end = 0

    # extract parameters
    for line in lines:
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


        if  np.isnan(line[13]):
            delta_self[j] = 0
        else:
            delta_self[j] = line[13]

        isotopologue_ID[j] = line[15]
        j = j + 1
    
    if temp_flag_end == False:
        visible_end = j

    return  x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self, visible_start, visible_end, isotopologue_ID

# Fetches uncertainties in physical parameters from user input
# Used to generate probablity distributions in a previous implementation
def exp_unc_uncertainties(mole_fraction, pathlength, pressure, temperature,exp_unc_values):
    print('generating distributions for experimental uncertainties')
    #molefraction_cdf, molefraction_range = rand_distribution(mole_fraction, mole_fraction * (1/100)*exp_unc_values[0] / 3)
    #pathlength_cdf, pathlength_range = rand_distribution(pathlength, pathlength * (1/100)*exp_unc_values[1] / 3)
    #pressure_cdf, pressure_range = rand_distribution(pressure, pressure * (1/100)*exp_unc_values[2] / 3)
    #temperature_cdf, temperature_range = rand_distribution(temperature, temperature * (1/100)*exp_unc_values[3] / 3)

    mole_fraction_sigma = mole_fraction * (1/100)*exp_unc_values[0] / 3
    pathlength_sigma = pathlength * (1/100)*exp_unc_values[1] / 3
    pressure_sigma = pressure * (1/100)*exp_unc_values[2] / 3
    temperature_sigma = temperature * (1/100)*exp_unc_values[3] / 3

    return mole_fraction_sigma, pathlength_sigma,pressure_sigma,temperature_sigma

# 'MC_simulation()' is the core function which performs the random sampling and calls the line shape calculation function
@st.cache_resource(show_spinner=False)
def MC_simulation(lines,n_simulations,T,P,mole_fraction,L,x,exp_unc_values,calc_method_wofz,simulation_type,num_of_isotopologues, first_isotopologue,x_limited,tips,mole_fraction_sigma, pathlength_sigma,pressure_sigma,temperature_sigma,x0_sigma, s0_sigma, gamma_air_sigma,gamma_self_sigma, n_air_sigma, delta_self_sigma,delta_air_sigma, x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self):
    # Run the simulations
    #t_1 = time.time()
    with tab1:
        my_bar = st.progress(0, text='Monte Carlo simulation progress')
    spectra = np.zeros((len(x), n_simulations))
    spectra_1 = np.zeros((len(x_limited), n_simulations))

    if not(st.session_state.exp_unc & (exp_unc_values != [0,0,0,0])):
        mole_fraction_1 = mole_fraction
        L_1 = L
        P_1 = P
        T_1 = T
    #t_sampling = 0
    #t_lineshape = 0
    # Outer loop, repeats as many times as requested by the user (n_simulations)
    t_1 = time.time()
    for i in range(n_simulations):
        
        j=0
        # if uncertainty in physical conditions is of interest, the lines within the conditional are executed
        if (st.session_state.exp_unc & (exp_unc_values != [0,0,0,0])):
            mole_fraction_1 = rng.normal(mole_fraction,mole_fraction_sigma)
            #random_value(molefraction_cdf, molefraction_range)
            L_1 = rng.normal(pathlength,pathlength_sigma)
            #random_value(pathlength_cdf, pathlength_range)
            P_1 = rng.normal(pressure,pressure_sigma)
            #random_value(pressure_cdf, pressure_range)
            T_1 = rng.normal(temperature,temperature_sigma)
            #random_value(temperature_cdf, temperature_range)

        # inner loop, calculates the absorbance/emission associated with each line contributing to the selected range
        for line in lines:

            # additional commands for cases with multiple isotopologues
            if num_of_isotopologues > 1:
                tips_index = line[15]-first_isotopologue
            else:
                tips_index = 0
            
            #t_2 = time.time()

            # == line position ==
            #x0_rand = random_value(x0_rand_cdf[j], x0_rand_range[j])
            x0_rand = rng.normal(x0[j],x0_sigma[j])

            # == bath-gas and self pressure shift ==
            #delta_air_rand = random_value(delta_air_rand_cdf[j], delta_air_rand_range[j])
            delta_air_rand = rng.normal(delta_air[j],delta_air_sigma[j])
            delta_self_rand = rng.normal(delta_self[j],delta_self_sigma[j])

            x0_shifted_rand = x0_rand + P_1 * ((1 - mole_fraction_1) * delta_air_rand + mole_fraction_1*delta_self_rand)

            # == Line strength ==
            #S0_rand = random_value(S0_rand_cdf[j], S0_rand_range[j])
            S0_rand = rng.normal(s0[j],s0_sigma[j])
            S_rand = S0_rand * (tips1(296,tips_index,tips) / tips1(T,tips_index,tips)) * np.exp(-(h * c * line[4] / kb) * (1 / T_1 - 1 / 296)) \
                     * (1 - np.exp(-h * c * x0_rand / (kb * T_1))) / (1 - np.exp(-h * c * x0_rand / (kb * 296)))

            A_rand = S_rand * L_1 * mole_fraction_1 * (P_1 / (R * T_1))  # cm-1

            # == Doppler broadening HWHM ==
            wG_rand = x0_shifted_rand * (7.1623E-7) * np.sqrt(T_1 / M)

            # == Pressure broadening ==
            # == temperature exponent for pressure broadening ==
            #n_air_rand = random_value(n_air_rand_cdf[j], n_air_rand_range[j])
            n_air_rand = rng.normal(n_air[j],n_air_sigma[j])

            # == collisional HWHM ==
            #gamma_self_rand = random_value(gamma_self_rand_cdf[j], gamma_self_rand_range[j])
            gamma_self_rand = rng.normal(gamma_self_0[j],gamma_self_sigma[j])
            #gamma_air_rand = random_value(gamma_air_rand_cdf[j], gamma_air_rand_range[j])
            gamma_air_rand = rng.normal(gamma_air_0[j],gamma_air_sigma[j])

            # == tempearture adjustment - collisional HWHM ==
            gamma_self_rand = gamma_self_rand * (296 / T_1) ** n_air_rand
            gamma_air_rand = gamma_air_rand * (296 / T_1) ** n_air_rand

            # == summation of collisional broadening coeffecients ==
            wL_rand = P_1 * (mole_fraction_1 * 2 * gamma_self_rand + (1 - mole_fraction_1) * 2 * gamma_air_rand)
            #t_sampling = t_sampling + (time.time() - t_2)
            #spectra[:, i] += (x, [A_rand, x0_shifted_rand, wG_rand, wL_rand])
            #t_3 = time.time()
            # detect of lower state energy is not listed correctly
            if not(np.isnan(line[4])):
                # lower state energy is not nan for this line
                # if it is, the spectrum will not be added
                # line shap function calculated using one of two methods depending on user choice
                if st.session_state.calc_method_wofz:
                    spectra[:, i] += voigtfwhm(x, [A_rand, x0_shifted_rand, wG_rand, wL_rand])
                else:
                    X = np.sqrt(np.log(2))*(x-x0_shifted_rand)/(0.5*wG_rand)
                    Y = np.sqrt(np.log(2))*((0.5*wL_rand)/(0.5*wG_rand))
                    spectra[:, i] += voigtfwhm_fast(x, [A_rand, X, Y])/ ((0.5*wG_rand) / np.sqrt(np.log(2)/np.pi))

                #elapsed = time.time() - t
                #print('time to calculate and add spectrum for a single line based on sampled parameters')
                #print(elapsed)
            #t_lineshape = t_lineshape + (time.time() - t_3)
            j=j+1

        # detect if any point in the calculated spectrum is corrupted, if so the spectrum is discarded 
        # / replaced by the spectrum from the previous iteration
        if np.isnan(np.mean(spectra[:, i])):
            #print('spectrum has nan values')
            spectra[:, i] = spectra[:, i - 1]
        
        # perform calculations to convert calculated absorbance spectrum to emission or transmittance, if requested
        if simulation_type == 'Emission':
            spectra[:, i] = np.multiply(plank_emission(x,T),(1-np.exp(-spectra[:, i])))
        elif simulation_type == 'Transmittance':
            spectra[:, i] = np.exp(-spectra[:, i])
        #else:
            #ax.plot(x, spectra[:, i], '.', markersize=0.1, color="#A87BF9")
        
        # add the spectrum instance to the array of spectra
        spectra_1[:, i] = np.interp(x_limited, x, spectra[:, i])
        # update progress bar for user
        
        if (i/10)%1 == 0:
            #print(i)
            #estimated_time = (time.time() - t_1)*(n_simulations-i)
            #my_bar.progress(i/n_simulations, text='Monte Carlo simulation progress (estimated remaining time: '+str(np.round(estimated_time/10,0))+' sec.)')
            my_bar.progress(i/n_simulations, text='Monte Carlo simulation progress ('+str(i)+' spectra instances simulated)')
            #t_1 = time.time()

    my_bar.empty()
    #print('t_sampling')
    #print(t_sampling)
    #print('t_lineshape')
    #print(t_lineshape)
    #print('t_total')
    #print(time.time()-t_1)
    return spectra_1

# Calculate spectrum based on mean parameters
@st.cache_resource(show_spinner=False,max_entries=3)
def mean_spectrum_simulation(lines,T,P,mole_fraction,L,x,calc_method_wofz,simulation_type,num_of_isotopologues,first_isotopologue,tips,x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self):
    spectrum_mean_parameters = np.zeros(len(x))
    j = 0
    for line in lines:
        # Line position
        if num_of_isotopologues > 1:
                tips_index = line[15]-first_isotopologue
        else:
            tips_index = 0
        
        # line position / including pressure shift
        x0_shifted = x0[j] + P * ((1 - mole_fraction) * delta_air[j] + mole_fraction * delta_self[j])
        # Line strength
        S = s0[j] * (tips1(296,tips_index,tips) / tips1(T,tips_index,tips)) * np.exp(-(h * c * line[4] / kb) * (1 / T - 1 / 296)) \
                 * (1 - np.exp(-h * c * x0[j] / (kb * T))) / (1 - np.exp(-h * c * x0[j] / (kb * 296)))
        A = S * L * mole_fraction * (P / (R * T))  # cm-1
        # Doppler broadening
        wG = x0_shifted * (7.1623E-7) * np.sqrt(T / M)
        # Pressure/collisional broadening
        gamma_self = gamma_self_0[j] * (296 / T) ** n_air[j]
        gamma_air = gamma_air_0[j] * (296 / T) ** n_air[j]
        wL = P * (mole_fraction * 2 * gamma_self + (1 - mole_fraction) * 2 * gamma_air)

        #spectra[:, i] += (x, [A_rand, x0_shifted_rand, wG_rand, wL_rand])
        # calculate line shape function depending on the function choice by user
        if not(np.isnan(line[4])):
            if st.session_state.calc_method_wofz:
                spectrum_mean_parameters +=  np.transpose(voigtfwhm(x, [A, x0_shifted, wG, wL]))
            else:
                X = np.sqrt(np.log(2))*(x-x0_shifted)/(0.5*wG)
                Y = np.sqrt(np.log(2))*((0.5*wL)/(0.5*wG))
                spectrum_mean_parameters += voigtfwhm_fast(x, [A, X, Y])/ ((0.5*wG) / np.sqrt(np.log(2)/np.pi))        
        
        j=j+1
    
    # perform calculations to convert calculated absorbance spectrum to emission or transmittance, if requested
    if simulation_type == 'Emission':      
        spectrum_mean_parameters = np.multiply(plank_emission(x,T),(1-np.exp(-spectrum_mean_parameters)))
    elif simulation_type == 'Transmittance':
        spectrum_mean_parameters = np.exp(-spectrum_mean_parameters)
        
    return spectrum_mean_parameters

# Calculate relative standard deviation, skewness and percentile based coeffecient of variation
# Percentiles are not presented in the interface of MCSpectra
@st.cache_resource(show_spinner=False,max_entries=3)
def calc_error_bars(spectra,spectrum_mean_parameters,x_limited):
    error_bars = np.zeros(len(x_limited))
    relative_uncertainty = np.zeros(len(x_limited))
    skewness = np.zeros(len(x_limited))
    
    P90 = np.zeros(len(x_limited))
    P10 = np.zeros(len(x_limited))
    PCV = np.zeros(len(x_limited))
    RMAD = np.zeros(len(x_limited))
    for i in range(len(x_limited)):
        relative_uncertainty[i] = 100*3*np.std(spectra[i][:])/np.mean(spectrum_mean_parameters[i])
        skewness[i] = skew(spectra[i][:])
        error_bars[i] = 3*np.std(spectra[i][:])
        P90[i] = np.percentile(spectra[i][:], 90)
        P10[i] = np.percentile(spectra[i][:], 10)
        # Devide RIQR by 1.35 to make it comparable to standard deviation
        PCV[i] = (3/1.35)*100*(P90[i] - P10[i])/np.median(spectra[i][:])
        RMAD[i] = (3/0.6745)*100*median_abs_deviation(spectra[i][:])/np.median(spectra[i][:])
    
    return relative_uncertainty, error_bars, skewness, P90, P10, PCV, RMAD

# build array with standard deviation value vs the number of iterations
# useful for testing convergence of the MC simulation  
def std_deviation_with_iterations(spectra,spectrum_mean_parameters,x_limited):
    std_residuals = np.zeros(N_simulations)
    if (st.session_state.wn_conv == wnend):
        std_index = np.where(np.round(x_limited,3) == (st.session_state.wn_conv - wnres))[0][0]
    else:
        std_index = np.where(np.round(x_limited,3) == (st.session_state.wn_conv))[0][0]
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
    
# calculate and plot distribution of results at wn_conv
@st.cache_resource(show_spinner=False,max_entries=3)
def uncertainty_PDF(convergence_frequency,spectra,spectrum_mean_parameters,skewness,error_bars,x_limited):
    
    # find the index at which the selected wavenumber position exists
    if (st.session_state.wn_conv == wnend):
        std_index = np.where(np.round(x_limited,3) == (st.session_state.wn_conv - wnres))[0][0]
    else:
        std_index = np.where(np.round(x_limited,3) == (st.session_state.wn_conv))[0][0]
    
    N_points = len(spectra[std_index])
    n_bins = 50
    fig_4, ax = plt.subplots()
    #print(std_residuals)
    if simulation_type == 'Absorbance':
        ax.set_xlabel('Absorbance')
        type_scale = 1
    elif simulation_type == 'Transmittance':
        ax.set_xlabel('Transmittance')
        type_scale = 1
    else:
        ax.set_xlabel('Emission (µW/(cm-1-cm2-sr)')
        type_scale = 1E+6
    ax.set_ylabel('Frequency')
    ax.grid(visible=True, linestyle='--', linewidth=0.5)
    
    # plottin command:
    ax.hist(type_scale*spectra[std_index], bins=n_bins, color="#A87BF9")

    ax.axvline(type_scale*spectrum_mean_parameters[std_index], color='white', linestyle='dashed', linewidth=2)
    textstr = (str(100*mole_fraction)+'% ' + selected_species + '\n' + str(T) + ' K\n' + str(P) + ' atm\n'+ str(L) + ' cm\n' + 'Bath-gas: '+ selected_broadener +'\n' + '# of simulations: ' + str(n_simulations) +'\n' + 'Skewness: ' + str(np.round(skewness[std_index],4))+'\n' + 'Std. Dev.: ' + str(np.round(0.333*error_bars[std_index],4)))
    props = dict(boxstyle='round', facecolor="#A87BF9", alpha=0)
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props) 

    return fig_4, std_index

# function which plots the figure which shows the absorbance/emission spectrum along with confidence interval
@st.cache_resource(show_spinner=False,max_entries=3)
def plot_MC_spectra(spectra, spectrum_mean_parameters, P90, P10, xaxis_start,xaxis_end,x_limited):
    fig_1, ax = plt.subplots()
    #ax.set_title('Absorbance based on mean parameters and '+str(n_simulations)+' simulated spectra')
    ax.set_xlabel('Wavenumbers (cm-1)')
    if simulation_type == 'Absorbance':
        ax.set_ylabel('Absorbance')
        type_scale = 1
    elif simulation_type == 'Transmittance':
        ax.set_ylabel('Transmittance')
        type_scale = 1
    else:
        ax.set_ylabel('Emission (µW/(cm-1-cm2-sr))')
        type_scale = 1E+6

    # plot the spectrum based on mean parameters, zorder makes sure it shows on top of the CI
    ax.plot(x_limited, type_scale*spectrum_mean_parameters, '-', color='white',zorder=n_simulations+1)
    #ax.plot(x_limited, type_scale*P90, '--', color="#ECBC7A",zorder=n_simulations+2)
    #ax.plot(x_limited, type_scale*P10, '--', color="#ECBC7A",zorder=n_simulations+3)
    #ax.plot(x, type_scale*spectrum_mean_parameters, '-', color='black', zorder=n_simulations+1)

    # for loop to plot all random spectra
    for i in range(n_simulations):
        ax.plot(x_limited, type_scale*spectra[:, i], '-', lw=1, color="#A87BF9")
    
    # adjust a-axis zoom based on the choice of the user
    ax.set_xlim(xaxis_start,xaxis_end)
    #round(time.time() - t,2)
    temp_index_start = np.where(np.round(x_limited,4) == xaxis_start)[0]
    temp_index_end = np.where(np.round(x_limited,4) == xaxis_end-wnres)[0]

    # scale -y-axis automatically
    if simulation_type == 'Transmittance':
        ax.set_ylim(0,1)
    else:
        ax.set_ylim(0,1.1*type_scale*spectra[temp_index_start[0]:temp_index_end[0]].max())#np.round(scaling_factor*spectra.max())/scaling_factor)

    ax.legend(['Mean parameters', 'Unceratinty envelope'])

    textstr = (str(100*mole_fraction)+'% ' + selected_species + '\n' + str(T) + ' K\n' + str(P) + ' atm\n'+ str(L) + ' cm\n' + 'Bath-gas: '+ selected_broadener +'\n' + '# of simulations: ' + str(n_simulations))
    props = dict(boxstyle='round', facecolor="#A87BF9", alpha=0)
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)

    ax.grid(visible=True, linestyle='--', linewidth=0.5)

    secax = ax.secondary_xaxis('top', functions=(wn2wl, wl2wn))
    secax.set_xlabel('Wavelength (µm)')
    
    return fig_1

# the following functions convert from wavenumber to wavelength, for the secondary x-axis
def wn2wl(x):
    return 10000 / x

def wl2wn(x):
    return 10000 / x

# plots mean spectrum for survey mode
def plot_mean_spectrum(spectrum_mean_parameters, xaxis_start,xaxis_end,x_limited):
    
    fig_1, ax = plt.subplots()
    #ax.set_title('Absorbance based on mean parameters and '+str(n_simulations)+' simulated spectra')
    ax.set_xlabel('Wavenumbers (cm-1)')
    if simulation_type == 'Absorbance':
        ax.set_ylabel('Absorbance')
        type_scale = 1
    elif simulation_type == 'Transmittance':
        ax.set_ylabel('Tranmittance')
        type_scale = 1
    else:
        ax.set_ylabel('Emission (µW/(cm-1-cm2-sr))')
        type_scale = 1E+6

    ax.plot(x_limited, type_scale*spectrum_mean_parameters, '-', color="#57D2E9",zorder=n_simulations+1)
    
    ax.set_xlim(xaxis_start,xaxis_end)

    temp_index_start = np.where(np.round(x_limited,4) == xaxis_start)[0]
    temp_index_end = np.where(np.round(x_limited,4) == xaxis_end-wnres)[0]
    
    if simulation_type == 'Transmittance':
        ax.set_ylim(0,1)
    else:
        ax.set_ylim(0,1.1*type_scale*spectrum_mean_parameters[temp_index_start[0]:temp_index_end[0]].max())#np.round(scaling_factor*spectra.max())/scaling_factor)

    ax.legend(['Mean parameters'])

    textstr = (str(100*mole_fraction)+'% ' + selected_species + '\n' + str(T) + ' K\n' + str(P) + ' atm\n'+ str(L) + ' cm\n' + 'Bath-gas: '+ selected_broadener +'\n')
    props = dict(boxstyle='round', facecolor="#A87BF9", alpha=0)
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)
    ax.grid(visible=True, linestyle='--', linewidth=0.5)

    secax = ax.secondary_xaxis('top', functions=(wn2wl, wl2wn))
    secax.set_xlabel('Wavelength (µm)')
    
    return fig_1

# relative uncertainty and skewness plots
@st.cache_resource(show_spinner=False,max_entries=3)
def plot_uncertainty(relative_uncertainty,skewness, PCV, RMAD, xaxis_start,xaxis_end,x_limited):
    fig_2, ax1 = plt.subplots()

    color = (1,1,1,1)
    ax1.set_xlabel('Wavenumbers (cm-1)')
    ax1.set_ylabel('Relative uncertainty - Δ(v) (%)', color=color)
    ax1.plot(x_limited, relative_uncertainty, color=color)
    #ax1.plot(x_limited, PCV, color="#A87BF9")
    #ax1.plot(x_limited, RMAD, linestyle='--', color="#A87BF9")
    ax1.tick_params(axis='y', labelcolor=color)
    if max(relative_uncertainty) < 100:
        ax1.set_ylim(0,100)
    ax1.set_xlim(xaxis_start,xaxis_end)
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

    textstr = (str(100*mole_fraction)+'% ' + selected_species + '\n' + str(T) + ' K\n' + str(P) + ' atm\n'+ str(L) + ' cm\n' + 'Bath-gas: '+ selected_broadener +'\n' + '# of simulations: ' + str(n_simulations))
    props = dict(boxstyle='round', facecolor="#A87BF9", alpha=0.1)
    #ax1.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10, verticalalignment='top', bbox=props)
    ax1.annotate(textstr, xy=(0, 1), xytext=(12, -12), va='top', xycoords='axes fraction', textcoords='offset points')

    fig_2.tight_layout()  # otherwise the right y-label is slightly clipped
    #st.divider() 
    
    #tab2 = st.tabs(['Global statistics'])
    return fig_2

# The following function updates 'dek' (key of manual parameter control table) when the (species, wavenumber range, temperature ... etc) are changed by the user 
# updating the key results in repopulating the table using HITRAN data
@st.cache_resource(show_spinner=False)
def update_manual_control(lines):
    st.session_state.dek = str(uuid.uuid4()) # refresh key to reset lines

@st.cache_resource(show_spinner=False)
def find_line_strength_thresh(s0_min, selected_broadener, T,P,mole_fraction, L, x, selected_species_lines,isotopologue_abundance,num_of_isotopologues,calc_method_wofz,simulation_type,first_isotopologue,tips,extended_start_x,extended_end_x,lines, x_limited, max_residual):
    # extract lines of interest from imported data
    testing_range = True
        
    # import line parameters from extracted lines
    x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self, visible_start, visible_end, isotopologue_ID = extract_mean_parameters(lines)
    
    # calculate a spectrum based on mean parameter to initiate the residual calculation loop 
    spectrum_mean_parameters = mean_spectrum_simulation(lines,T,P,mole_fraction,L,x,calc_method_wofz,simulation_type,num_of_isotopologues,first_isotopologue,tips,x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self)
    spectrum_mean_parameters = np.interp(x_limited, x, spectrum_mean_parameters)
    residual = 1
    with tab1:
        with st.spinner('Computing line strength threshold ...'):
            # loop to satisfy the risidual requirement
            while residual > max_residual:
                # decrease line strength threshold by a factor of 10
                s0_min = s0_min/2#/10
                # repeat line extraction
                lines = extract_lines(extended_start_x,extended_end_x,selected_species_lines,s0_min, selected_broadener, testing_range,isotopologue_abundance)
                # repeat parameter extraction
                x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self, visible_start, visible_end, isotopologue_ID = extract_mean_parameters(lines)
                # repeat spectrum calculation
                extended_spectrum_mean_parameters = mean_spectrum_simulation(lines,T,P,mole_fraction,L,x,calc_method_wofz,simulation_type,num_of_isotopologues,first_isotopologue,tips,x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self)
                extended_spectrum_mean_parameters = np.interp(x_limited, x, extended_spectrum_mean_parameters)
                # calculate maximum residual
                residual = max(extended_spectrum_mean_parameters - spectrum_mean_parameters)/max(extended_spectrum_mean_parameters)
                
                #residual = max(np.divide(extended_spectrum_mean_parameters - spectrum_mean_parameters,spectrum_mean_parameters))
                
                # update spectrum
                spectrum_mean_parameters = extended_spectrum_mean_parameters

    s0_min = 2*s0_min
    with st.sidebar:
        st.sidebar.info(f'Line strength threshold set at {s0_min:.1e} (cm-1/(molec.cm-2)). Max residual: {residual:.2%}.', icon="ℹ️")
    #st.session_state.s0_min = s0_min
    #print('line strength threshold')
    #print(s0_min)

    return s0_min, lines, spectrum_mean_parameters

@st.cache_resource(show_spinner=False)
def find_cut_off_wn(rotational_constant,max_residual,wnstart,wnend,wnres,selected_species_lines,s0_min, selected_broadener,isotopologue_abundance,T,P,mole_fraction,L,calc_method_wofz,simulation_type,num_of_isotopologues,first_isotopologue,tips,x_limited):
    # evaluate line spacing for frequency cut-off calculations
    testing_range = True
    line_spacing = 2*rotational_constant
    if line_spacing > 10:
        wn_step = line_spacing
    else:
        wn_step = 10
    # gradual expansion of covered wavelength range
    wn_cutoff = 0
    x = np.arange(wnstart - wn_cutoff, wnend + wn_cutoff, wnres)
    extended_start_x, extended_end_x = find_range(x,selected_species_lines)
    lines = extract_lines(extended_start_x,extended_end_x,selected_species_lines,s0_min, selected_broadener, testing_range,isotopologue_abundance)
    x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self, visible_start, visible_end, isotopologue_ID = extract_mean_parameters(lines)
    spectrum_mean_parameters = mean_spectrum_simulation(lines,T,P,mole_fraction,L,x,calc_method_wofz,simulation_type,num_of_isotopologues,first_isotopologue,tips,x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self)
    with tab1:
        with st.spinner('Computing cut-off frequency ...'):
            residual = 1
            while residual > max_residual:
                start_x = extended_start_x
                end_x = extended_end_x
                wn_cutoff = wn_cutoff + wn_step
                x = np.arange(wnstart - wn_cutoff, wnend + wn_cutoff, wnres)
                extended_start_x, extended_end_x = find_range(x,selected_species_lines)
                lines = extract_lines(extended_start_x,extended_end_x,selected_species_lines,s0_min, selected_broadener, testing_range,isotopologue_abundance)
                #edited_lines = lines
                x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self, visible_start, visible_end, isotopologue_ID = extract_mean_parameters(lines)
                extended_spectrum_mean_parameters = mean_spectrum_simulation(lines,T,P,mole_fraction,L,x,calc_method_wofz,simulation_type,num_of_isotopologues,first_isotopologue,tips,x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self)
                extended_spectrum_mean_parameters = np.interp(x_limited, x, extended_spectrum_mean_parameters)
                # reevaluate extended_spectrum_mean_parameters at original x_range points and store to extended_spectrum_mean_parameters
                #residual = max(np.divide(extended_spectrum_mean_parameters - spectrum_mean_parameters,extended_spectrum_mean_parameters))
                residual = max(extended_spectrum_mean_parameters - spectrum_mean_parameters)/max(extended_spectrum_mean_parameters)
                
                spectrum_mean_parameters = extended_spectrum_mean_parameters
                #print(wn_cutoff)
                #print(residual)
    
    wn_cutoff = wn_cutoff - wn_step
    # inform user about the selected cut-off frequency 
    with st.sidebar:
        st.sidebar.success('Lines within '+str(wnstart - wn_cutoff)+' - '+str(wnend + wn_cutoff)+' cm-1 will be included in the simulation for enhanced accuracy within the selected wavenumber range.', icon="✅")

    return wn_cutoff, start_x, end_x


#xaxis_validation()
#print(wn_validation_flag)
def main(s0_min,max_residual,selected_species,wnstart, wnend, wnres, selected_broadener,T,P,mole_fraction,L,calc_method_wofz,simulation_type):
    t = time.time()
    # import HITRAN data from csv file    
    selected_species_lines, tips, num_of_isotopologues, first_isotopologue, isotopologue_abundance, rotational_constant = import_data(selected_species)

    # 'x_limited' is an array representing the wavenumber range selected by the user
    # used below for interpolating spectra calculated over an extended range
    x_limited = np.arange(wnstart, wnend, wnres)

    # expand wavenumber range before testing for threshold range
    # this would help in case a range is chosen where no lines exist within range, but do exist just outside of the range
    # start_x and end_x are indices which define all lines from the database that are within the selected range (before extending the range)
    # extended_start_x and extended end_x are indices which define all lines from the database that are relavant to the selected range (after extending the range)
    # s0_min set at '0' to include all lines during frequency cut-off determination
    if not(st.session_state.survey_mode):
        wn_cutoff, extended_start_x, extended_end_x = find_cut_off_wn(rotational_constant,max_residual,wnstart,wnend,wnres,selected_species_lines,0, selected_broadener,isotopologue_abundance,T,P,mole_fraction,L,calc_method_wofz,simulation_type,num_of_isotopologues,first_isotopologue,tips,x_limited)
    else:
        wn_cutoff = 20*P
        x = np.arange(wnstart - wn_cutoff, wnend + wn_cutoff, wnres)
        extended_start_x, extended_end_x = find_range(x,selected_species_lines)
        with st.sidebar:
            st.sidebar.success('Lines within '+str(wnstart - wn_cutoff)+' - '+str(wnend + wn_cutoff)+' cm-1 will be included in the simulation for enhanced accuracy within the selected wavenumber range.', icon="✅")
    # 'lines' is the list of lines that will go into the spectra simulations after excluding lines below the line-strength threshold
    testing_range = True
    lines = extract_lines(extended_start_x, extended_end_x,selected_species_lines,s0_min, selected_broadener, testing_range,isotopologue_abundance)
    number_of_lines_limited = len(lines)
    
    # while loop to determine the starting line strength threshold
    while (number_of_lines_limited == 0) and (s0_min > 1E-40):
        s0_min = 0.1*s0_min
        lines = extract_lines(extended_start_x, extended_end_x,selected_species_lines,s0_min, selected_broadener, testing_range,isotopologue_abundance)
        number_of_lines_limited = len(lines)


    number_of_lines_limited = len(lines)
    if number_of_lines_limited == 0:
        max_residual = 1
        st.warning('No lines within the selected frequency range.', icon="⚠️")
    else:
        x = np.arange(wnstart - wn_cutoff, wnend + wn_cutoff, wnres)
        s0_min, lines, spectrum_mean_parameters = find_line_strength_thresh(s0_min, selected_broadener, T,P,mole_fraction, L, x, selected_species_lines,isotopologue_abundance,num_of_isotopologues,calc_method_wofz,simulation_type,first_isotopologue,tips,extended_start_x,extended_end_x,lines, x_limited, max_residual)
    
        testing_range = False
        # final extraction of lines and parameters after theshold and cut-off have been determined
        x = np.arange(wnstart - wn_cutoff, wnend + wn_cutoff, wnres)
        extended_start_x, extended_end_x = find_range(x,selected_species_lines)
        lines = extract_lines(extended_start_x,extended_end_x,selected_species_lines,s0_min, selected_broadener, testing_range,isotopologue_abundance)
        x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self, visible_start, visible_end, isotopologue_ID = extract_mean_parameters(lines)

        # setup for manual control table
        if st.session_state.manual_control:
            # reset values if list of simulation parameters have changed
            update_manual_control(lines)
                       
            st.divider() 
            st.write('_Editable list of parameters for lines within the selected range:_')
            # table configuration
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
        ),
        16: st.column_config.NumberColumn(
            "ID.",
            help="Isotopologue ID",
            disabled=True,
        )
    }
            
            # instantiate the table
            # the table only shows (and allows modifying) lines within the selected/viewable wavelength range
            edited_lines_limited = st.data_editor(lines[visible_start:visible_end],column_config=lines_column_config,key=st.session_state.dek)
            st.divider()
            # bring back the lines which are outside of the selected range
            edited_lines = lines[0:(visible_start)] + edited_lines_limited + lines[(visible_end):len(lines)]
        else:
            edited_lines = lines

        #print(edited_lines == lines)
        with tab1:
            with st.spinner('Extracting line parameters ...'):
                # redo parameter extraction (needed for cases where manual control is avtive)
                number_of_lines = len(edited_lines)
                x0_sigma, s0_sigma, gamma_air_sigma\
                ,gamma_self_sigma, n_air_sigma, delta_self_sigma,\
                delta_air_sigma, x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self, isotopologue_ID = extract_parameters(edited_lines)

        # Prepare physical parameter uncertainties
        if (st.session_state.exp_unc & (exp_unc_values != [0,0,0,0]) & (st.session_state.survey_mode == 0)):
            mole_fraction_sigma, pathlength_sigma,pressure_sigma,temperature_sigma = exp_unc_uncertainties(mole_fraction, pathlength, pressure, temperature,exp_unc_values)
        else:
            mole_fraction_sigma, pathlength_sigma,pressure_sigma,temperature_sigma = [0,0,0,0]
        
        # call monte carlo simulations function
        if not(st.session_state.survey_mode):
            spectra_limited = np.zeros((len(x_limited), n_simulations))    
            spectra_limited = MC_simulation(edited_lines,n_simulations,T,P,mole_fraction,L,x,exp_unc_values, calc_method_wofz,simulation_type,num_of_isotopologues, first_isotopologue, x_limited,tips,mole_fraction_sigma, pathlength_sigma,pressure_sigma,temperature_sigma,x0_sigma, s0_sigma, gamma_air_sigma,gamma_self_sigma, n_air_sigma, delta_self_sigma,delta_air_sigma, x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self)

        # the following line may be used to store the randomly generated spectra locally
        #np.savetxt('sample_results/spectra_for_analysis.csv', spectra_limited, delimiter=',')

        # calculate spectrum based on mean parameters
        with tab1:
            with st.spinner('Computing spectrum based on mean parameters ...'):
                extended_spectrum_mean_parameters = mean_spectrum_simulation(edited_lines,T,P,mole_fraction,L,x,calc_method_wofz,simulation_type,num_of_isotopologues,first_isotopologue,tips,x0, s0, gamma_air_0, gamma_self_0, n_air, delta_air, delta_self)
                extended_spectrum_mean_parameters = np.interp(x_limited, x, extended_spectrum_mean_parameters)
                spectrum_mean_parameters = extended_spectrum_mean_parameters

        # plot and show absorbance/emission spectra
        with tab1:
            with st.spinner('Plotting simulated spectra ...'):
                if not(st.session_state.survey_mode):
                    relative_uncertainty, error_bars, skewness, P90, P10, PCV, RMAD = calc_error_bars(spectra_limited,spectrum_mean_parameters,x_limited)
                    fig_1 = plot_MC_spectra(spectra_limited, spectrum_mean_parameters, P90, P10, xaxis_start,xaxis_end,x_limited)
                else:
                    fig_1 = plot_mean_spectrum(spectrum_mean_parameters, xaxis_start,xaxis_end,x_limited)
        with tab1:
            if not(st.session_state.survey_mode):
                st.write('_'+simulation_type+' spectrum based on mean line-parameters of '+str(number_of_lines) + ' lines,\nand '+str(n_simulations)+' spectra based on randomly sampled line-parameters:_')
            else:
                st.write('_'+simulation_type+' spectrum based on mean line-parameters of '+str(number_of_lines) + ' lines\n:_')
                        
            st.pyplot(fig_1)


        if not(st.session_state.survey_mode):
            with tab1:
                with st.spinner('Calculating and plotting global statistics ...'):
                    fig_2 = plot_uncertainty(relative_uncertainty,skewness, PCV, RMAD, xaxis_start,xaxis_end,x_limited)
            with tab2:
                st.write('_Relative uncertainty and skewness spectra:_')
                st.pyplot(fig_2)    

            if conv_test == 1:
                fig_3, std_index = std_deviation_with_iterations(spectra_limited,spectrum_mean_parameters,x_limited)
                with tab3:
                    st.write('_Standard deviation with iterations at ('+str(round(x_limited[std_index],2))+' cm-1):_')
                    st.pyplot(fig_3)

            with tab1:
                with st.spinner('Calculating and plotting PDF at ('+str(round(x_limited[std_index],2))+' cm-1)...'):
                    fig_4, std_index = uncertainty_PDF(convergence_frequency,spectra_limited,spectrum_mean_parameters,skewness,error_bars,x_limited)
            with tab4:
                st.write('_Histogram of predicted absorbance at ('+str(round(x_limited[std_index],2))+' cm-1). Dashed line indicates predicted absorbance based on mean parameters:_')
                st.pyplot(fig_4)
        else:
            with tab2:
                st.write('_Disable survey mode for uncertainty quantification._')
            with tab3:
                st.write('_Disable survey mode for uncertainty quantification._')
            with tab4:
                st.write('_Disable survey mode for uncertainty quantification._')

        with st.sidebar:
            st.sidebar.info('Total computation time: '+str(round(time.time() - t,2))+' seconds.', icon="ℹ️")


        simulation_info = [datetime.datetime.now(),selected_species,T,P,mole_fraction,L,wnstart,wnend,wnres,n_simulations,s0_min,st.session_state.manual_control,conv_test,st.session_state.survey_mode]
        #print(simulation_info)
        with open('simulation_history.csv','a') as fd:
            #fd.write(np.array2string(simulation_info))
            writer = csv.writer(fd)
            writer.writerow(simulation_info)

        # prepare data in csv format for download
        if not(st.session_state.survey_mode):
            arr = np.array(np.transpose([x_limited,spectrum_mean_parameters,error_bars,relative_uncertainty,skewness]))
            arr_df = pd.DataFrame(arr,columns=['wavenumbers (cm-1)','absorbance - mean parameters','3 x standard deviation','relative uncertainty','skewness'])
        else:
            arr = np.array(np.transpose([x_limited,spectrum_mean_parameters]))
            arr_df = pd.DataFrame(arr,columns=['wavenumbers (cm-1)','absorbance / emission - mean parameters'])
        # Create an in-memory buffer

        # download button
        with tab1:
            with st.spinner('Preparing downloadable data ...'):
                with io.BytesIO() as buffer:
                    # Write array to buffer
                    #np.savetxt(buffer, arr_df, delimiter=",")
                    textstr = ('MCSpectra_'+selected_species + '_' + str(100*mole_fraction)+'_' + '_'+ str(simulation_type) + '_'+ str(T) + '_K_' + str(P) + '_atm_'+ str(L) + '_cm_' + str(wnstart) + '_'+ str(wnend)+ '_Bath-gas_' + str(selected_broadener))
                    buffer = pd.DataFrame.to_csv(arr_df, sep=',', index=False, encoding='utf-8')

                    st.download_button(
                        label="Download results as CSV",
                        data = buffer, # Download buffer
                        file_name = textstr+'.csv',
                        mime='text/csv'
                    )

wn_validation_flag, wn_change_flag = wn_validation()
# run main() function if simulation parameters are valid
if wn_validation_flag == 1:
    if not(st.session_state.survey_mode):
        rng = np.random.default_rng()
    main(s0_min,max_residual,selected_species,wnstart, wnend, wnres, selected_broadener,T,P,mole_fraction,L,calc_method_wofz,simulation_type)
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