//The Warming, Acidification and Sea-level Projector: WASP (Septembe 2020)
//WASP is an 8-box representation of the Earth system for efficient very-large ensemble simulations.
//Coded by Philip Goodwin, This version 01-September-2020

//Citations for use of code:

//This version of the WASP model:

//Goodwin, P. and Cael, B.B. (2020 submitted) Bayesian estimation of Earth’s climate sensitivity and transient climate response from observational warming and heat content datasets, submitted to Earth System Dynamics discussions.


//Previous version (including time-varying climate feedbacks):
//Goodwin (2018) On the time Evolution of Climate Sensitivity and Future Warming, Earth's Future

////Goodwin, P., A. Katavouta,, V. M. Roussenov, G. L. Foster, E. J. Rohling and R. G. Williams (2018), Pathways to 1.5 and 2 °C warming based on observational and geological constraints, Nature Geoscience

//WASP Model initial description of basic equations:
//Goodwin, P. (2016) How historic simulation-observation discrepancy affects future warming projections in a very large model ensemble, Climate Dynamics, CLDY-D-15-00368R2, doi: 10.1007/s00382-015-2960-z.

//Utilising the input distribution for climate sensitivity, S=1/lambda, from geological data requires citation of:
//Rohling, E. J. et al - Palaeosens Project Members (2012), Making sense of Palaeoclimate sensitivity, Nature 491, 683-691, doi:10.1038/nature11574.

//The theory behind how WASP calculated atmospheric CO2 and surface temperature anomaly is given here:
//Goodwin, P., R.G. Williams and A. Ridgwell (2015), Sensitivity of climate to cumulative carbon emissions due to compensation of ocean heat and carbon uptake, Nature Geoscience, Vol. 8, p29-34. doi:10.1038/ngeo2304.

//If you use the sea level component of WASP (commented out in the code), also cite:
//Goodwin, P., I. D. Haigh, E. J. Rohling, and A. Slangen (2017), A new approach to projecting 21st century sea-level changes and extremes, Earth’s Future, 5, doi:10.1002/2016EF000508.



////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Code compiled on a Mac Darwin environment, operating system OS X version 10.15.6, using command line compilation on the Terminal.
//To compile: place the files WAS2_ESM_main.cpp and WASP2_ESM_functions.cpp together in a directory. Add a folder called 'RESULTS' for the results files to be written into.
//Using compiler g++ code can be compiled with command line: g++ -O3 -o filename Main_WASP_ESM.cpp

//This code has *not* been checked across a range of C++ compilers.
//Also, the sequence of pseudo-random numbers may be compiler and platform dependent





//Main programme - include the following...
#include <cmath>
#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include "WASPv3_ESM_functions.cpp"
#include "WASPv3_ESM_stored_inputs.cpp"


using namespace std;

int main()
{
    //Use the clock to generate a random number seed...
    std::time_t SEED = std::time(nullptr);
    //Use the following line to submit multiple jobs at to the queue at once and know that none will have the same SEED - by varying the integer of days
    SEED = SEED - 1*(24*60*60); //integer before bracket = number of days ago to vary the SEED
    
    cout << SEED << endl;
    
    InitialiseGlobalVectors();
    set_radiative_forcing_coeffs();
    
    setEm_RCP85();
    
    //includes one of the random number generators.
    std::minstd_rand0 generator (SEED);
 
    //seeds rand() random number generator related to time.
    srand(SEED);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //Choose the RCP forcing scenario//////////////////////////////////////////////////////////////////////////////
    int Future_scenario;
    int Use_stored_params = 0; //Set equal to 0 for Monte Carlo run, set equal to 1 for usinbg stored values
    
    //Set Future_scenario = 1 for RCP2.6, =2 for RCP4.5, =3 for RCP6.0 and =4 for RCP8.5, 5=AMP, 6=ssp119, 7=ssp126, 8=ssp245, 9=ssp370, 10=ssp370-lowNTCF, 11=ssp434, 12=ssp460, 13=ssp534-over, 14=ssp585, 15=1pcentCO2 for 151 years, 16=1pcent CO2 4xCO2 extension, 17 = abrupt 2xCO2, 18=abrupt 4xCO2, 19=abrupt 0.5xCO2
    Future_scenario = 4;
    
    int T_data =  0; //set T_data = 0 for HadCRUT4 and T_data = 1 for Cowtan&Way v2.0.0 data constraints
    int OHC_data = 1; //set OHC_data = 0 for Cheng et al., and OHC_data = 1 for NODC  data constraints. NODC uses IPCC AR5 (2015) whole system heat content change since NODC only goes to 2000m.
    
    //Choose the number of simulations in the initial ensemble/////////////////////////////////////////////////////
    int SCENARIOS; //This is the number of simulations in the initial ensemble - from which a smaller number will be extracted by the observational consistency tests
    
    SCENARIOS = 3250000; //Perform this number of Monte Carlo simulations in the initial ensemble
    if(Future_scenario >= 15 || Use_stored_params == 1)
    {
        InitialiseInputVectors();
        SCENARIOS = a_CO2_obs.size();
    }
    if(Future_scenario < 15 && Use_stored_params == 0)
    {
        //Stores observation-consistent parameter values
        a_CO2_obs.clear();
        dNPPdT_obs.clear();
        gamma_K_obs.clear();
        dtaudT_obs.clear();
        ratio_obs.clear();
        ratio2_obs.clear();
        lambda_Planck_obs.clear();
        lambda_WVLR_obs.clear();
        lambda_Cloud_Fast_obs.clear();
        lambda_Cloud_SST_obs.clear();
        lambda_albedo_obs.clear();
        f_heat_ocean_obs.clear();
        tau_C_mixed_obs.clear();
        tau_C_upper_obs.clear();
        tau_C_inter_obs.clear();
        tau_C_deep_obs.clear();
        tau_C_bottom_obs.clear();
        Ib_obs.clear();
        tau_WVLR_obs.clear();
        tau_Cloud_Fast_obs.clear();
        tau_albedo_obs.clear();
        tau_Cloud_SST_obs.clear();
        R_aerosol_Uncert_obs.clear();
        R_WMGHG_nonCO2_Uncert_obs.clear();
        R_volcanic_coeff_obs.clear();
        gamma_aero_SOx_obs.clear();
        gamma_aero_BC_obs.clear();
        gamma_aero_OC_obs.clear();
        gamma_aero_ENMVOC_obs.clear();
        gamma_aero_NOx_obs.clear();
        gamma_aero_NH3_obs.clear();
        AeroRF_ind_2011_obs.clear();
        Uncert_CH4_obs.clear();
        Uncert_N2O_obs.clear();
        Uncert_Halocarbons_obs.clear();
    }
    
    double OHC_Cheng_etal;
    double OHC700m_Cheng_etal;
    double DT1_GISTEMP;
    double DT2_GISTEMP;
    double DT1_HadCRUT4;
    double DT2_HadCRUT4;
    double DT1_Cowtan_Way;
    double DT2_Cowtan_Way;
    
    
    
    //This is target warming if on AMP scenario
    double Target1 = 2.0;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //To weight Rnon-CO2 from WMGHG in between two scenarios
    double R_WMGHG_nonCO2_b[tmax];
    double RnonCO2_weight = 1.0;
    

    
    
    //print to screen which RCP scenario is chosen and how many ensemble members are being generated
    if(Future_scenario == 1)
        cout << "Generating an initial ensemble of " <<  SCENARIOS << " simulations following RCP2.6" << endl;
    if(Future_scenario == 2)
        cout << "Generating an initial ensemble of " <<  SCENARIOS << " simulations following RCP4.5" << endl;
    if(Future_scenario == 3)
        cout << "Generating an initial ensemble of " << SCENARIOS <<  " simulations following RCP6.0" << endl;
    if(Future_scenario == 4)
        cout << "Generating an initial ensemble of " << SCENARIOS <<  " simulations following RCP8.5" << endl;
    if(Future_scenario == 5)
        cout << "Generating an initial ensemble of " << SCENARIOS <<  " simulations stabilising at " << Target1 << "degrees C warming" << endl;
    if(Future_scenario == 6)
        cout << "Generating an initial ensemble of " << SCENARIOS <<  " simulations following ssp119" << endl;
    if(Future_scenario == 7)
        cout << "Generating an initial ensemble of " << SCENARIOS <<  " simulations following ssp126" << endl;
    if(Future_scenario == 8)
        cout << "Generating an initial ensemble of " << SCENARIOS <<  " simulations following ssp245" << endl;
    if(Future_scenario == 9)
        cout << "Generating an initial ensemble of " << SCENARIOS <<  " simulations following ssp370" << endl;
    if(Future_scenario == 10)
        cout << "Generating an initial ensemble of " << SCENARIOS <<  " simulations following ssp370-lowNTCF" << endl;
    if(Future_scenario == 11)
        cout << "Generating an initial ensemble of " << SCENARIOS <<  " simulations following ssp434" << endl;
    if(Future_scenario == 12)
        cout << "Generating an initial ensemble of " << SCENARIOS <<  " simulations following ssp460" << endl;
    if(Future_scenario == 13)
        cout << "Generating an initial ensemble of " << SCENARIOS <<  " simulations following ssp534-over" << endl;
    if(Future_scenario == 14)
        cout << "Generating an initial ensemble of " << SCENARIOS <<  " simulations following ssp585" << endl;
    if(Future_scenario == 15)
        cout << "Using and observation-consistent ensemble of " << SCENARIOS <<  " simulations following 1pctCO2" << endl;
    if(Future_scenario == 16)
        cout << "Using and observation-consistent ensemble of " << SCENARIOS <<  " simulations following 1pctCO2 4xCO2 extension" << endl;
    if(Future_scenario == 17)
        cout << "Using and observation-consistent ensemble of " << SCENARIOS <<  " simulations following abrupt 2xCO2" << endl;
    if(Future_scenario == 18)
        cout << "Using and observation-consistent ensemble of " << SCENARIOS <<  " simulations following abrupt 4xCO2" << endl;
    if(Future_scenario == 19)
        cout << "Using and observation-consistent ensemble of " << SCENARIOS <<  " simulations following abrupt 0.5xCO2" << endl;
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    
    
    //Sets up the output files///////////////////////////////////////////////////////////////////
    
    //The model must be run from a location containing a folder called 'RESULTS' for the output files below to be written into//////////////
    
    //Outputs results of the model at the end of the year 2016, when the observational consistency tests are applied////////////////////
    ofstream results ("./RESULTS/Outputs.txt", std::ios::app);
    
    //Outputs the model input values at 2016 in the observation-consistent ensemble members/////////////////////////////////////////////
    
    //Format:
    ofstream inputs ("./RESULTS/Inputs.txt", std::ios::app);
    
    //Outputs specific projections of temperature anonaly and carbon emissions for the 21st century/////////////////////////////////////
    
    //Format: ensemble-member number, mean temperature anonaly relative to preindustrial at 2081-2100, change in temperature anomaly 1986-2005 to 2081-2100, change in temperature anomaly 1986-2005 to 2046-2065, cumulative emissions 2017 to 2100, cumulative emissions 2012 to 2100
    ofstream projection ("./RESULTS/Projection.txt", std::ios::app);
    
  
   
    //The following files output annual properties by year from 1765 to 2100 for all observation-consistent ensemble members////////////
    //All anomnalies are relative to the preindustrial//////////////////////////////////////////////////////////////////////////////////
    
    //Formats: Rows are values by year from 1765 to 2100, values are separated by a tab/////////////////////////////////////////////////
    
    //Effective climate sensitivity in K using air temperature
    ofstream annual_S_effective ("./RESULTS/EffCS.txt", std::ios::app);
    //Effective lambda in Wm/2/K over time for all sources of radiative forcing using air temperature
    ofstream annual_lambda ("./RESULTS/lambda.txt", std::ios::app);
    //Instantaneous TCRE - using air-temperature (annual average?)
    ofstream Instantaneous_TCRE ("./RESULTS/Instantaneous_TCRE.txt", std::ios::app);
    
    //Surface air temperature anomaly relative to 1850-1900 average///////
    ofstream annual_warming ("./RESULTS/Warming.txt", std::ios::app);
    //Surface blended temperature anomaly relative to 1850-1900 average///////
    ofstream annual_warming_blended ("./RESULTS/Warming_blended.txt", std::ios::app);
    //Sea Surface temperature anomaly relative to 1850-1900 average///////
    ofstream annual_warming_SST ("./RESULTS/Warming_SST.txt", std::ios::app);
    
    
    //Radiative forcing from all sources
    ofstream annual_R_tot ("./RESULTS/R_tot.txt", std::ios::app);
    //Radiative forcing from all anthropogenic sources
    ofstream annual_R_anth ("./RESULTS/R_anth.txt", std::ios::app);
    //Radiative forcing from solar activity
    ofstream annual_R_solar ("./RESULTS/R_solar.txt", std::ios::app);
    //Radiative forcing from volcanic eruptions
    ofstream annual_R_volcanic("./RESULTS/R_volcanic.txt", std::ios::app);
    //Radiative forcing from CO2
    ofstream annual_RCO2 ("./RESULTS/RCO2.txt", std::ios::app);
    //efficacy-weighted radiative forcing from aerosols and other non-Kyoto agents
    ofstream annual_R_aerosol ("./RESULTS/R_aerosol.txt", std::ios::app);
    //Radiative forcing from WMGHG other than CO2)
    ofstream annual_R_WMGHG_nonCO2 ("./RESULTS/R_WMGHG_nonCO2.txt", std::ios::app);
    //Net TOA energy imbalance
    ofstream annual_N_tot ("./RESULTS/N_tot.txt", std::ios::app);
    
    //Annual heat content anonaly of the upper 700m of the ocean
    ofstream annual_heat700m ("./RESULTS/OHC_700m.txt", std::ios::app);
    //Annual heat content anonaly of the whole ocean
    ofstream annual_heatOcean ("./RESULTS/OHC_Ocean.txt", std::ios::app);
    //Annual heat content anonaly of the Earth system
    ofstream annual_heatEarth ("./RESULTS/OHC_700_2000m.txt", std::ios::app);
    
    //Residual terrestrial carbon storage anonaly relative to previous year
    ofstream annual_landC ("./RESULTS/LandC.txt", std::ios::app);
    //Residual terrestrial carbon storage anonaly relative to previous year
    ofstream annual_oceanC ("./RESULTS/OceanC.txt", std::ios::app);
    //Compatible cumulative emissions
    ofstream annual_Em ("./RESULTS/CumulativeEm.txt", std::ios::app);
    //Equivalent carbon emissions from ocean heating
    ofstream annual_I_eq_Heat ("./RESULTS/I_eq_heat.txt", std::ios::app);
    //Amount of carbon in the soil carbon pool
    ofstream landC_soil ("./RESULTS/SoilC.txt", std::ios::app);
    //Amount of carbon in the plant carbon pool
    ofstream landC_plant ("./RESULTS/PlantC.txt", std::ios::app);
    //Amount of carbon in the atmospheric CO2 pool (MtC)
    ofstream atmosCO2 ("./RESULTS/CO2.txt", std::ios::app);
    //NPP carbon flux (MtC/yr)
    ofstream NPP_flux ("./RESULTS/NPP_flux.txt", std::ios::app);
    //Airborne fraction of total emitted carbon
    ofstream Airborne_fraction ("./RESULTS/Airborne_fraction.txt", std::ios::app);
    
    
    //File for storing the observation-consistent input parameter vectors, which can then be used to run the idealised experiments
    ofstream stored_inputs ("./RESULTS/WASPv3_ESM_stored_inputs.cpp", std::ios::app);
    
    stored_inputs << "//Stored inputs for the Warming Acidification and Sea level Projector, version 3 or higher" << endl;
    stored_inputs << "#include <cmath>" << endl << "#include <iostream>" << endl << "#include <fstream>" << endl << "#include <random>" << endl << "#include <vector>" << endl << "using namespace std;" << endl << "std::vector<double> a_CO2_obs;" << endl << "std::vector<double> dNPPdT_obs;" << endl << "std::vector<double> gamma_K_obs;" << endl << "std::vector<double> dtaudT_obs;" << endl << "std::vector<double> ratio_obs;" << endl << "std::vector<double> ratio2_obs;" << endl << "std::vector<double> lambda_Planck_obs;" << endl << "std::vector<double> lambda_WVLR_obs;" << endl << "std::vector<double> lambda_Cloud_Fast_obs;" << endl << "std::vector<double> lambda_Cloud_SST_obs;" << endl << "std::vector<double> lambda_albedo_obs;" << endl << "std::vector<double> f_heat_ocean_obs;" << endl << "std::vector<double> tau_C_mixed_obs;" << endl << "std::vector<double> tau_C_upper_obs;" << endl << "std::vector<double> tau_C_inter_obs;" << endl << "std::vector<double> tau_C_deep_obs;" << endl << "std::vector<double> tau_C_bottom_obs;" << endl << "std::vector<double> Ib_obs;" << endl << "std::vector<double> tau_WVLR_obs;" << endl << "std::vector<double> tau_Cloud_Fast_obs;" << endl << "std::vector<double> tau_albedo_obs;" << endl << "std::vector<double> tau_Cloud_SST_obs;" << endl << "std::vector<double> R_aerosol_Uncert_obs;" << endl << "std::vector<double> R_WMGHG_nonCO2_Uncert_obs;" << endl << "std::vector<double> R_volcanic_coeff_obs;" << endl <<  "std::vector<double> gamma_aero_SOx_obs;" << endl << "std::vector<double> gamma_aero_BC_obs;" << endl << "std::vector<double> gamma_aero_OC_obs;"  << endl << "std::vector<double> gamma_aero_ENMVOC_obs;"  << endl << "std::vector<double> gamma_aero_NOx_obs;"  << endl << "std::vector<double> gamma_aero_NH3_obs;" << endl << "std::vector<double> AeroRF_ind_2011_obs;"  << endl << "std::vector<double> Uncert_CH4_obs;" << endl <<  "std::vector<double> Uncert_N2O_obs;" << endl <<  "std::vector<double> Uncert_Halocarbons_obs;" << endl << "void InitialiseInputVectors()" << endl << "{" << endl;
   
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    //The following files output properties every 5PgC emitted relative to the carbon emitted at the start of 2017/////////////////////
    
    //Format: values are given every 5Pg cumulative emissions relative to start of 2017, starting at -400PgC. Values are tab-separated//
    
    //Surface temperature anomnaly relative to the 1850-1900 average
    ofstream em_warming_T_2020 ("./RESULTS/Em_Warm_T_2020.txt", std::ios::app);
    
    //Cumulative emissions (this output of cumulative emissions is to check that the file above is working correctly, and clarify what the cumulative emissions relative to 2017 are in the file above)
    ofstream em_warming_I_2017 ("./RESULTS/Em_Warm_I_2020.txt", std::ios::app);
    
    
    //The sensitivities of temperature, radiative forcing and carbon emissions/////////////////////////////
    
    //Format: rows are values given by year from 1765 to 2100. Values are tab separated////////////////////
    
    //Note in above files the un-weighted radiative forcing is used (i.e. no efficacy weighting is applied)
    
    //ofstream Adjustment_TCRE ("./RESULTS/AdjustmentTCRE_INDC_SimTSL_AP2_coeff0.txt", std::ios::app);
    
    //ofstream Adjustment_Emrate ("./RESULTS/AdjustmentEmrate_INDC_SimTSL_AP2_coeff0.txt", std::ios::app);
    
    //ofstream Adjustment_DT ("./RESULTS/AdjustmentDT_INDC_SimTSL_AP2_coeff0.txt", std::ios::app);
    
    //ofstream Adjustment_Emtime ("./RESULTS/AdjustmentEmtime_INDC_SimTSL_AP2_coeff0.txt", std::ios::app);
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
   
    
    //This counts the number of observation-consistent simulations
    int count_obs = 0;

    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //////////////////////////////////////////////////////////////////////////////
    //Parameters for assessing the TCRE and making a responsive emissions decision
    double Assess_year1 = 2030.0;
    double Assess_year2 = 2050.0;
    double Assess_year3 = 2070.0;
    
    double TCRE_estimate = 0.0;
    double Warming_in_period = 0.0;
    double Emissions_in_period = 0.0;
    double EmRate_Assessment=0.0;
    
    double I_remaining = 0.0;
    double t_remaining = 0.0;
    
    double Coeff_Target = 0.0;
    double TCRE_safety = 0.0;//0.5/1.0e3;
    
    double Restore_year = 2018.0;
    int Start_year = 1765;
    if(Future_scenario >= 6)
    {
        Start_year = 1700;
    }
    int Restore_year_int = 2018;
    
    
    int overshoot = 0;
    //////////////////////////////////////////////////////////////////////////////

    //Definition of variables - while initial values are given, many parameters are re-assigned different values later
    const int n_forcings = 6;
    
    int count_annual = 0;   //this is for writing output
    int year_count = 0;     //this is for requesting forcing values at sub-annual time steps from forcing vectors at annual resolution
    int t_step_count = 0;
    
    const int t_step_per_year = 48; //Number of timsteps per year
    const int t_step_per_month = t_step_per_year/12; //Number of timesteps per month
    
    int data_fit = 0; //Integer that counts every time a simulation passes an observation consistency test
    int data_fit2 = 0; //Integer that counts every time a simulation misses an observational consistency test
    double error_term = 0.0;//The total error outside the 90%/95% ranges of the observational consistency tests
    double cost_function = 1.0;
    double cost_function_max = 0.0;
    
    const int RESTORE = 1; //if RESTORE = 1, then emissions are set to restore CO2 to RCP scenarios. If RESTORE=0 then emissions are set regardless of CO_2 values.
    
    //These provide the ratios of lambda values for diffferent feedback types (relative to CO2) at equilibrium, for different forcing types [n]
    double ratioCO2_WVLR[n_forcings];
    double ratioCO2_CloudFast[n_forcings];
    double ratioCO2_albedo[n_forcings];
    double ratioCO2_CloudSST[n_forcings];
    double ratioCO2_1000yr[n_forcings];
    
    
    
    //Set all ratios to 1.0 initially, then can vary below in the programme if desired.
    for(int n=0; n<n_forcings; n++)
    {
        ratioCO2_WVLR[n]=1.0;
        ratioCO2_CloudFast[n]=1.0;
        ratioCO2_albedo[n]=1.0;
        ratioCO2_CloudSST[n]=1.0;
        ratioCO2_1000yr[n]=1.0;
    }
    
    //Equilibrium climate variables/////////////
    double a_CO2 = 5.35; //W m-2 from Myhre et al [1998] 
    double lambda = 1.18;// 1.25; // W m-2 K-1 - this is equilibrium climate parameter, equal to 1/climate sensitivity, S
    double Ib = 3500.0;  //PgC - the buffered carbon inventory of the air-sea system
    
    double lambda_time[n_forcings][tmax];
    //double lambda_volcanic_time[tmax];
    
    //The Planck climate sensitivity feedback parameter acts instantaneously
    double lambda_Planck = 3.15;// W m-2 K-1 - the Plank sensitivity
    
    //These are the climate sensitivity feedback parameter values that act over longer timescales, at equilibrium for each feedback
    double lambda_WVLR = -1.15;         //Combined Water Vapour - Lapse Rate
    double lambda_Cloud_Fast = -0.44;   //Clouds with original SST patterns
    double lambda_Cloud_SST = -0.47;    //Clouds with final SST patterns
    double lambda_albedo = -0.37;       //Fast snow and ice-albedo
    double lambda_1000yr = -0.5;        //Longer timescale processes, vegetation, dust, ice-sheet etc
    
    
    //These are the values that the feedbacks take over time //Could make these arrays with n number of radiative forcings
    double lambda_WVLR_time[n_forcings];
    double lambda_Cloud_Fast_time[n_forcings];
    double lambda_Cloud_SST_time[n_forcings];
    double lambda_albedo_time[n_forcings];
    double lambda_1000yr_time[n_forcings];
    
    
    
    //e-folding timescales for feedbacks to operate in years
    double tau_WVLR = 10.0/365.0;
    double tau_Cloud_Fast = 10.0/365.0;
    double tau_albedo = 5.0;
    double tau_Cloud_SST = 30.0;
    
    double NoverR = 1.0;
    double RF_fraction_slow;
    double RF_fraction_total;
    
    
    double ratio = 0.8; //Ratio of global mean ocean warming to atmopsheric surface warming at equilibrium
    double ratio2 = 0.6; //Ratio of SST change to sub-surface ocean temperature change at equilibrium
    double S_ocean = (1.0/lambda)*ratio; //effective SST ocean climate sensitivity (K [Wm-2]-1)
    
    //Required constants////////////////////////////////
    double CO2init = 590; //initial atmoispheric CO2 (PgC)
    double PgCtoppm = 278.0/CO2init; //constant converting CO2 in PgC to ppm
    
    double dt = 1.0/double (t_step_per_year); //time step in years
    double area = 5.1e8*1.0e6; //surface area of planet earth over which radiative forcing acts.
    double area_ocean = 3.5e8*1.0e6; //surface area of the ocean to allow conversion of N to effective ocean heat flux (in W m-2)
    
    
    //Parameters linked to radiative foricng
    double R_aerosol_Uncert=0.0;
    double R_aerosol_2011 = -0.6506;
    double R_volcanic_coeff = 0.0;
    double Uncert_CH4 = 1.0;
    double Uncert_N2O = 1.0;
    double Uncert_Halocarbons = 1.0;
    
    double R_WMGHG_nonCO2_Uncert=0.0;
    double R_WMGHG_nonCO2_2011 = 0.695;
    
    //Separate aerosol radiative forcing factors
    double gamma_aero_SOx=0.0;
    double gamma_aero_BC=0.0;
    double gamma_aero_OC=0.0;
    double gamma_aero_ENMVOC=0.0;
    double gamma_aero_NOx=0.0;
    double gamma_aero_NH3=0.0;
    
    //Indirect aerosol radiative forcing in 2011
    double AeroRF_ind_2011;
    
    
    std::vector<double> R_WMGHG_nonCO2;//The non-CO2 from the Kyoto protocol gasses
    R_WMGHG_nonCO2.assign(tmax, 0.0);
    
    //Tunable parameters/////////////////////////////////
    double tau_C_mixed = 0.5; //years. e-folding timescale for carbon equilibration of surface mixed layer
    
    //years. e-folding timescales for tracer equilibration of upper, intermediate deep and bottom ocean boxes with surface mixed layer.
    double tau_C_upper = 25.0;
    double tau_C_inter = 25.0;
    double tau_C_deep = 400;
    double tau_C_bottom = 1500.0;
    
    
    double vol_mixed = 0.03;   //fraction of the ocean volume in the surface mixed layer
                        
    double vol_bottom = 0.40;  //fraction of ocean volume in bottom water section
    
    double vol_deep =  0.21;   //fraction of the ocean volume originating as North Atlantic Deep Water
    
    double vol_inter = 0.20;   //fraction of the ocean volume originating as intermediate waters
    
    double vol_upper =  (1.0 - vol_bottom - vol_mixed - vol_deep - vol_inter); //Fraction of ocean volume in upper (e.g. ventilatted thermocline) waters
    
    double vol_total = 1.3e18; //m3 vol of total ocean
    
    //Global mean heat capacity of seawater
    double c_p = 3.91e3*1026.0;// 3.91e3*1026; //J K-1 kg-1 * kg m-3  - values from Williams et al 2012
    
    
    double H_mixed_equil; //The equilibrium heat content anomaly for the surface mixed layer with respect to the instantaneous radiative forcing
    double H_mixed_equil_comp[n_forcings];
    
    //Fraction of total ocean heat uptake going into the ocean
    double f_heat_ocean=0.93;
    
    double year[tmax]; //year
    
    
    
    
    
    //Fluxes of carbon undersaturation between boxes///////////////////////
    double I_Usat_MixToDeep;
    double I_Usat_MixToUpper;
    double I_Usat_MixToInter;
    double I_Usat_MixToBottom;
    
    
    double N_mixed[tmax];      //Mixed layer heat uptake (W m-2 averaged across *whole planet's* surface area)
    
    double CO2_restore[tmax];  //CO2 concentration used to restore towards (e.g. to set to RCP pathway)
    
    //The T, Heat and Carbon undersaturation can be declared as vectors (slow, but uses less stack memory) ot arrays (faster but can only use if there is enough stack memory)
    
    //double DeltaT[tmax];       //Mean surface temperature anomaly relative to preindustrial
    
    /*
    std::vector<double> DeltaT;
    DeltaT.assign(tmax, 0.0);
    */
    double DeltaT[tmax];
    /*
    std::vector<double> Heat_mixed; //Mixed layer cumulative anthropogenic heat content
    Heat_mixed.assign(tmax,0.0);
    std::vector<double> Heat_upper; //Upper ocean cumulative anthropogenic heat content increase
    Heat_upper.assign(tmax,0.0);
    std::vector<double> Heat_inter; //Intermediate box heat content anomaly
    Heat_inter.assign(tmax,0.0);
    std::vector<double> Heat_deep;  //deep ocean cumulative heat content increase
    Heat_deep.assign(tmax,0.0);
    std::vector<double> Heat_bottom;//bottom ocean cumulative anthrpogenic heat uptake
    Heat_bottom.assign(tmax,0.0);
    */
    double Heat_mixed[tmax];
    double Heat_upper[tmax];
    double Heat_inter[tmax];
    double Heat_deep[tmax];
    double Heat_bottom[tmax];
    
    //Carbon undersaturation inventories of the ocean boxes////////////////
    /*
    std::vector<double> I_Usat_mixed;
    I_Usat_mixed.assign(tmax,0.0);
    std::vector<double> I_Usat_upper;
    I_Usat_upper.assign(tmax,0.0);
    std::vector<double> I_Usat_inter;
    I_Usat_inter.assign(tmax,0.0);
    std::vector<double> I_Usat_deep;
    I_Usat_deep.assign(tmax,0.0);
    std::vector<double> I_Usat_bottom;
    I_Usat_bottom.assign(tmax,0.0);
    */
    double I_Usat_mixed[tmax];
    double I_Usat_upper[tmax];
    double I_Usat_inter[tmax];
    double I_Usat_deep[tmax];
    double I_Usat_bottom[tmax];
    
    
    double Em_heat = 0.0; //equivalent emission during time-step from temperature-CO2 solubility feedback (after Goodwin and Lenton, 2009 in GRL)
    
    double dCsatdT = -0.1/1.0e15; //PgC per K-1 m-3
    
    
    std::vector<double> I_eq_heat;     //Cumulative equivalent emissions from heating-CO2 solubility feedback
    I_eq_heat.assign(tmax, 0.0);
    
    
    std::vector<double> Steric_rise;    //To estimate steric contribution to sea level rise (m)
    Steric_rise.assign(tmax, 0.0);
    std::vector<double> ice_melt_rise; //To estimate ice_melt contribution to sea level rise (m)
    ice_melt_rise.assign(tmax, 0.0);
    std::vector<double> Sea_surf_acid; //to estimate sea surface acidification (Delta pH units)
    Sea_surf_acid.assign(tmax, 0.0);
    
    
    std::vector<double> I_veg;          //terrestrial cerbon stored in vegetation (PgC)
    I_veg.assign(tmax, 0.0);
    std::vector<double> I_soil;         //terrestrial carbon stored in soil (PgC)
    I_soil.assign(tmax, 0.0);
    std::vector<double> I_ter;  //total terrestrial carbon (PgC)
    I_ter.assign(tmax, 0.0);

    cout << "Hello" << endl;
    
    //Initial values at preindustrial/////
    double I_veg_init = 450.0;   //PgC
    double I_soil_init = 1500.0; //PgC
    
    //Net Primary Productivity in terrestrial ecosystem, PgC yr-1
    double NPP_init = 60; //Initial net primary productivity PgC yr-1
    
    double dIvegdt, dIsoildt; //rates of change of soil and vegetation carbon pools with time.
    
    //Keeling formula growth factor
    double gamma_K = 0.36; //growth factor 
    
    //Partial sensitivity of NPP to T
    double dNPPdT = -1.5; //in PgC yr-1 K-1
    
    double dtaudT = -0.50; //(years K-1) sensitivity of global mean soil carbon residence time to global mean temperature
    
    
    double ice_coeff = 0.0020; //To estimate ice_melt sea level rise, coefficient in m (K.yr)-1. Model after Rahmstorf (2007) and Vermeer and Rahmstorf (2009) 
    double steric_coeff=1.236; //Steric sea level rise coefficient from Williams et al (2012). 1.236 [mm yr-1 (W m-2)-1] is the sea level rise per year per unit ocean heat flux.

    double land_water=0.0;
    
    //Surface warming noise parameters, AR(2) noise
    ///////////////////////////////////////////////
    
    //if applying external noise to T
    std::vector<double> T_noise;            //external noise applied
    T_noise.assign(tmax, 0.0);
    
    double gamma1_noise = 0.3;    //0.8//0.65
    double gamma2_noise = 0.4;
    double alpha_noise = 0.062;   //0.0075//0.01
    double z_noise = 0.0;
    ///////////////////////////////////////////////
    
    //if applying external noise to RF
    R_noise.assign(tmax, 0.0);              //external noise on radiative forcing
    double gamma_Rnoise = 0.75; //fractional relaxation in internal RF noise per timestep
    double alpha_Rnoise = 0.78; //maximum growth of internal RF noise in W per m2 per timestep
    
    //////////////////////////////////////////////////////////////////////////////
    
    int T_restore_check = 0;
    
    for(int i = 0; i<tmax; i++) //Defines the year in terms of timesteps
    {
        year[i] = double(Start_year) + double(i)*dt;
    }
    
    //Set first year conditions - initialise variables to zero/////////////////////
    I_Usat_mixed[0] = 0.0;
    I_Usat_upper[0]=0.0; 
    I_Usat_inter[0] = 0.0;
    I_Usat_deep[0]=0.0;  
    I_Usat_bottom[0]=0.0;
    
    N_mixed[0]=0.0;
   
    CO2[0] = CO2init;
    DeltaT[0] = 0.0;
    
    Steric_rise[0] = 0.0;   
    Sea_surf_acid[0] = 0.0;
    
    
    I_veg[0] = I_veg_init;
    I_soil[0] = I_soil_init;
    I_ter[0] = I_veg[0] + I_soil[0];
    
    //////////////////////////////////////////////////////////////////////////////
    
    setRnonCO2_zero(); //initially zero the non-CO2 radiative forcing
    
    
    
    int output = 0;
    
    //Make volcanoes have larger feedback (smaller sensitivity) after Gregory et al (2016)
    double ratio_volcanic_CO2=1.0;
    ratioCO2_WVLR[0] = 1.0;//ratio_volcanic_CO2;
    ratioCO2_CloudFast[0] = 1.0;//ratio_volcanic_CO2;
    ratioCO2_albedo[0] = ratio_volcanic_CO2;//ratio_volcanic_CO2;
    lambda_Cloud_SST_time[0] = 1.0;//ratio_volcanic_CO2;
    
    
    for (int j = 0; j<SCENARIOS; j++) //This loop runs 'SCENARIOS' number of ensemble members
    {
        output = 0;
        T_restore_check = 0;
        year_count = 0; //set to first year of simulation
        t_step_count = 0;
       
    
        
        
        //Monte-carlo approach to selecting the internal model parameters such as climate sensitivity///////
        if(RESTORE == 1 )
        {
            
            //Set up the RCP scenario - note need to choose CO2, radiative frocing from other Kyoto agents and radiative forcing from non-Kyoto agents individually
            //Note that RCP3PD = RCP2.6
            
            //Select a scneario by commenting out all but one scnearios from the setCO2RCPXX, setRnonCO2_RCPXXK and setRnonCO2_RCPXX functions/////////////////////
            //These must all match!!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
            if(Future_scenario==1) //Sets up RCP2.6 (RCP3PD) as the forcing
            {
                setCO2RCP3PD(PgCtoppm);
                convertCO2(t_step_per_year); //Converts the CO2 scenario to the sub-annual timesteps;
                
                setR_WMGHG_nonCO2_RCP3PD();
                //setRnonCO2_RCP85_K();
                convertRnonCO2(t_step_per_year);
                //Set R_WMGHG_nonCO2 gasses array, refers to gasses inside the kyoto protocol
                for(int n=0; n<tmax; n++)
                {
                    R_WMGHG_nonCO2[n] = RnonCO2[n]; //Radiative frocing from CH4, N2O and Halogens
                }
                
                
                //Set Aerosol Optical Depth for stratospheric (volcanic) aerosols (Must have 1 month timestep to use this function!!!!!)
                
                setAOD(Start_year, t_step_per_month);
                convertAOD(t_step_per_month);
                
                setR_volcanicIPCC();
                convertRnonCO2(t_step_per_year);
                for(int n=0; n<tmax; n++)
                {
                    R_volcanic[n] = RnonCO2[n]; //Radiative forcing from volcanic aerosols
                }
                
                
                //Set tropospheric (largely ntrhopogenic) aerosol RF///////////////////
                setR_aerosol_RCP3PD();
                convertRnonCO2(t_step_per_year);
                //Set R_aerosols
                for(int n=0; n<tmax; n++)
                {
                    R_aerosol[n] = RnonCO2[n];//Radiative frocing from CH4, N2O and Halogens
                }
                
                setR_solar();
                convertRnonCO2(t_step_per_year);
                //Set R_solar array
                for(int n=0; n<tmax; n++)
                {
                    R_solar[n] = RnonCO2[n]; //Radiative frocing from solar changes
                }
                

                
            }
            
            if(Future_scenario == 2)
            {
                setCO2RCP45(PgCtoppm);
                convertCO2(t_step_per_year); //Converts the CO2 scenario to the sub-annual timesteps;
                
                setR_WMGHG_nonCO2_RCP45();
                //setRnonCO2_RCP85_K();
                convertRnonCO2(t_step_per_year);
                //Set R_WMGHG_nonCO2 gasses array, refers to gasses inside the kyoto protocol
                for(int n=0; n<tmax; n++)
                {
                    R_WMGHG_nonCO2[n] = RnonCO2[n]; //Radiative frocing from CH4, N2O and Halogens
                }
                
                
                //Set Aerosol Optical Depth for stratospheric (volcanic) aerosols (Must have 1 month timestep to use this function!!!!!)
                
                setAOD(Start_year, t_step_per_month);
                convertAOD(t_step_per_month);
                
                
                //Below uses RF definition for volcanic aerosols
                setR_volcanicIPCC();
                convertRnonCO2(t_step_per_year);
                for(int n=0; n<tmax; n++)
                {
                    R_volcanic[n] = RnonCO2[n]; //Radiative forcing from volcanic aerosols
                }
                
                
                //Set tropospheric (largely ntrhopogenic) aerosol RF///////////////////
                setR_aerosol_RCP45();
                convertRnonCO2(t_step_per_year);
                //Set R_aerosols
                for(int n=0; n<tmax; n++)
                {
                    R_aerosol[n] = RnonCO2[n];//Radiative frocing from CH4, N2O and Halogens
                }
                
                setR_solar();
                convertRnonCO2(t_step_per_year);
                //Set R_solar array
                for(int n=0; n<tmax; n++)
                {
                    R_solar[n] = RnonCO2[n]; //Radiative frocing from solar changes
                }
                

                
            }
            
            if(Future_scenario == 3)
            {
                setCO2RCP6(PgCtoppm);
                convertCO2(t_step_per_year); //Converts the CO2 scenario to the sub-annual timesteps;
                
                setR_WMGHG_nonCO2_RCP6();
                //setRnonCO2_RCP85_K();
                convertRnonCO2(t_step_per_year);
                //Set R_WMGHG_nonCO2 gasses array, refers to gasses inside the kyoto protocol
                for(int n=0; n<tmax; n++)
                {
                    R_WMGHG_nonCO2[n] = RnonCO2[n]; //Radiative frocing from CH4, N2O and Halogens
                }
                
                
                //Set Aerosol Optical Depth for stratospheric (volcanic) aerosols (Must have 1 month timestep to use this function!!!!!)
                
                setAOD(Start_year, t_step_per_month);
                convertAOD(t_step_per_month);
                
                //Set alternative radiative forcing from aerosols
                setR_volcanicIPCC();
                convertRnonCO2(t_step_per_year);
                for(int n=0; n<tmax; n++)
                {
                    R_volcanic[n] = RnonCO2[n]; //Radiative forcing from volcanic aerosols
                }
                
                
                //Set tropospheric (largely ntrhopogenic) aerosol RF///////////////////
                setR_aerosol_RCP6();
                convertRnonCO2(t_step_per_year);
                //Set R_aerosols
                for(int n=0; n<tmax; n++)
                {
                    R_aerosol[n] = RnonCO2[n];//Radiative frocing from CH4, N2O and Halogens
                }
                
                setR_solar();
                convertRnonCO2(t_step_per_year);
                //Set R_solar array
                for(int n=0; n<tmax; n++)
                {
                    R_solar[n] = RnonCO2[n]; //Radiative frocing from solar changes
                }
                

                
            }
            
            if(Future_scenario == 4)
            {
                //Set emissions profiles
                setEm_RCP85();
                
                setCO2RCP85(PgCtoppm);
                convertCO2(t_step_per_year); //Converts the CO2 scenario to the sub-annual timesteps;
                
                setR_WMGHG_nonCO2_RCP85();
                //setRnonCO2_RCP85_K();
                convertRnonCO2(t_step_per_year);
                //Set R_WMGHG_nonCO2 gasses array, refers to gasses inside the kyoto protocol
                for(int n=0; n<tmax; n++)
                {
                    R_WMGHG_nonCO2[n] = RnonCO2[n]; //Radiative frocing from CH4, N2O and Halogens
                }
                
               
                //Set Aerosol Optical Depth for stratospheric (volcanic) aerosols (Must have 1 month timestep to use this function!!!!!)
                
                setAOD(Start_year, t_step_per_month);
                convertAOD(t_step_per_month);
                
                //Set alternative radiative forcing from aerosols
                setR_volcanicIPCC();
                convertRnonCO2(t_step_per_year);
                for(int n=0; n<tmax; n++)
                {
                    R_volcanic[n] = RnonCO2[n]; //Radiative forcing from volcanic aerosols
                }
                
                
                //Set tropospheric (largely ntrhopogenic) aerosol RF///////////////////
                setR_aerosol_RCP85();
                //setRnonCO2_RCP85();
                convertRnonCO2(t_step_per_year);
                //Set R_aerosols
                for(int n=0; n<tmax; n++)
                {
                    R_aerosol[n] = RnonCO2[n] ;//- R_WMGHG_nonCO2[n];//Radiative frocing from CH4, N2O and Halogens
                }
                
                setR_solar();
                convertRnonCO2(t_step_per_year);
                //Set R_solar array
                for(int n=0; n<tmax; n++)
                {
                    R_solar[n] = RnonCO2[n]; //Radiative frocing from solar changes
                }
                
                
            }
            
            
            if(Future_scenario==5) //Sets up RCP2.6 (RCP3PD) as the forcing
            {
                setCO2RCP45_plusGCB(PgCtoppm);
                convertCO2(t_step_per_year); //Converts the CO2 scenario to the sub-annual timesteps;
                
                
                setR_WMGHG_nonCO2_RCP45();
                convertRnonCO2(t_step_per_year);
                
                for (int n = 0; n<tmax; n++)
                {
                    R_WMGHG_nonCO2_b[n] = RnonCO2[n];
                }
                
                
                setR_WMGHG_nonCO2_RCP3PD();
                convertRnonCO2(t_step_per_year);
                //Set R_WMGHG_nonCO2 gasses array, refers to gasses inside the kyoto protocol
                for(int n=0; n<tmax; n++)
                {
                    if(year[n] < 2018.0 + dt/2.0)
                    {
                        R_WMGHG_nonCO2[n] = RnonCO2[n];
                    }
                    if(year[n] > 2018.0 + dt/2.0)
                    {
                        R_WMGHG_nonCO2[n] = (1.0-RnonCO2_weight)* RnonCO2[n] + (RnonCO2_weight)*R_WMGHG_nonCO2_b[n]; //Radiative frocing from CH4, N2O and Halogens
                    }
                    
                }
                
                
                //Set Aerosol Optical Depth for stratospheric (volcanic) aerosols (Must have 1 month timestep to use this function!!!!!)
                
                setAOD(Start_year, t_step_per_month);
                convertAOD(t_step_per_month);
                
                //Set alternative radiative forcing from aerosols
                setR_volcanicIPCC();
                convertRnonCO2(t_step_per_year);
                for(int n=0; n<tmax; n++)
                {
                    R_volcanic[n] = RnonCO2[n]; //Radiative forcing from volcanic aerosols
                }
                
                
                //Set tropospheric (largely anthropogenic) aerosol RF///////////////////
                setR_aerosol_RCP3PD();
                convertRnonCO2(t_step_per_year);
                //Set R_aerosols
                for(int n=0; n<tmax; n++)
                {
                    R_aerosol[n] = RnonCO2[n];//Radiative frocing from CH4, N2O and Halogens
                }
                
                setR_solar();
                convertRnonCO2(t_step_per_year);
                //Set R_solar array
                for(int n=0; n<tmax; n++)
                {
                    R_solar[n] = RnonCO2[n]; //Radiative frocing from solar changes
                }
                
                
                
            }
            
            if(Future_scenario == 6)
            {
                setCO2_ssp119(PgCtoppm);
                set_concs_ssp119();
                convertCO2(t_step_per_year);
                for (int n = 0; n<tmax; n++)
                {
                    CO2[n] = CO2[n]/PgCtoppm;
                }
            }
            if(Future_scenario == 7)
            {
                setCO2_ssp126(PgCtoppm);
                set_concs_ssp126();
                convertCO2(t_step_per_year);
                for (int n = 0; n<tmax; n++)
                {
                    CO2[n] = CO2[n]/PgCtoppm;
                }
            }
            if(Future_scenario == 8)
            {
                setCO2_ssp245(PgCtoppm);
                set_concs_ssp245();
                convertCO2(t_step_per_year);
                for (int n = 0; n<tmax; n++)
                {
                    CO2[n] = CO2[n]/PgCtoppm;
                }
            }
            if(Future_scenario == 9)
            {
                setCO2_ssp370(PgCtoppm);
                set_concs_ssp370();
                convertCO2(t_step_per_year);
                for (int n = 0; n<tmax; n++)
                {
                    CO2[n] = CO2[n]/PgCtoppm;
                }
            }
            if(Future_scenario == 10)
            {
                setCO2_ssp370lowNTCF(PgCtoppm);
                set_concs_ssp370_lowNTCF();
                convertCO2(t_step_per_year);
                for (int n = 0; n<tmax; n++)
                {
                    CO2[n] = CO2[n]/PgCtoppm;
                }
            }
            if(Future_scenario == 11)
            {
                setCO2_ssp434(PgCtoppm);
                set_concs_ssp434();
                convertCO2(t_step_per_year);
                for (int n = 0; n<tmax; n++)
                {
                    CO2[n] = CO2[n]/PgCtoppm;
                }
            }
            if(Future_scenario == 12)
            {
                setCO2_ssp460(PgCtoppm);
                set_concs_ssp460();
                convertCO2(t_step_per_year);
                for (int n = 0; n<tmax; n++)
                {
                    CO2[n] = CO2[n]/PgCtoppm;
                }
            }
            if(Future_scenario == 13)
            {
                setCO2_ssp534(PgCtoppm);
                set_concs_ssp534_over();
                convertCO2(t_step_per_year);
                for (int n = 0; n<tmax; n++)
                {
                    CO2[n] = CO2[n]/PgCtoppm;
                }
            }
            if(Future_scenario == 14)
            {
                //setCO2_ssp585(PgCtoppm);
                set_concs_ssp585();
                convertCO2(t_step_per_year);
                for (int n = 0; n<tmax; n++)
                {
                    CO2[n] = CO2[n]/PgCtoppm;
                }
                
                
                
                
            }
                
            if(Future_scenario >=6 && Future_scenario <= 14)
            {
                
                //Sets RF from solar and volcanic sources, automatically linearly interpolates adjusts to time step
                //Using RCMIP prior to 1850, and NASA AOD after 1850
                setR_solar_volcanic_1700(t_step_per_year);
                setAOD(Start_year, t_step_per_month);
                convertAOD(t_step_per_month);
                
            }
            
            
            
            if(Future_scenario == 15)
            {
                setCO2_1pcent(PgCtoppm);
                
                convertCO2(t_step_per_year); //Converts the CO2 scenario to the sub-annual timesteps;
                
                for(int n=0; n<tmax; n++)
                {
                    AOD[n] = 0.0;
                    R_solar[n] = 0.0;
                    R_WMGHG_nonCO2[n] = 0.0;
                    RnonCO2[n] = 0.0;
                    R_volcanic[n] = 0.0;
                    R_aerosol[n] = 0.0;
                    
                }
            }
            
            if(Future_scenario == 16)
            {
                setCO2_1pcent4Xext(PgCtoppm);
                convertCO2(t_step_per_year); //Converts the CO2 scenario to the sub-annual timesteps;
                
                for(int n=0; n<tmax; n++)
                {
                    AOD[n] = 0.0;
                    R_solar[n] = 0.0;
                    R_WMGHG_nonCO2[n] = 0.0;
                    RnonCO2[n] = 0.0;
                    R_volcanic[n] = 0.0;
                    R_aerosol[n] = 0.0;
                    
                }
            }
            
            if(Future_scenario == 17)
            {
                setCO2_abrupt2X(PgCtoppm);
                convertCO2(t_step_per_year); //Converts the CO2 scenario to the sub-annual timesteps;
                
                for(int n=0; n<tmax; n++)
                {
                    AOD[n] = 0.0;
                    R_solar[n] = 0.0;
                    R_WMGHG_nonCO2[n] = 0.0;
                    RnonCO2[n] = 0.0;
                    R_volcanic[n] = 0.0;
                    R_aerosol[n] = 0.0;
                    
                }
            }
            
            if(Future_scenario == 18)
            {
                setCO2_abrupt4X(PgCtoppm);
                convertCO2(t_step_per_year); //Converts the CO2 scenario to the sub-annual timesteps;
                
                for(int n=0; n<tmax; n++)
                {
                    AOD[n] = 0.0;
                    R_solar[n] = 0.0;
                    R_WMGHG_nonCO2[n] = 0.0;
                    RnonCO2[n] = 0.0;
                    R_volcanic[n] = 0.0;
                    R_aerosol[n] = 0.0;
                    
                }
            }
            
            if(Future_scenario == 19)
            {
                setCO2_abrupt0_5X(PgCtoppm);
                convertCO2(t_step_per_year); //Converts the CO2 scenario to the sub-annual timesteps;
                
                for(int n=0; n<tmax; n++)
                {
                    AOD[n] = 0.0;
                    R_solar[n] = 0.0;
                    R_WMGHG_nonCO2[n] = 0.0;
                    RnonCO2[n] = 0.0;
                    R_volcanic[n] = 0.0;
                    R_aerosol[n] = 0.0;
                    
                }
            }
            
            
            
            //Calculate 2011 values for the standard case with no uncertainty and scaling to match observations///
            R_WMGHG_nonCO2_2011 = R_WMGHG_nonCO2[(2011-1765)*t_step_per_year];
            R_aerosol_2011 = R_aerosol[(2011-1765)*t_step_per_year];
            
            //////////////////////////////////////////////////////////////////////////////////////////////////////
            
            ////////////////////////////////////////////////////////////
            
            //The section chooses pseudo-random values for internal model parameters, such that each of the 'SCENARIOS' simulations is different
            if(Future_scenario<15 && Use_stored_params == 0)
            {
                
                //Random normal value for radiative forcing coefficient
                a_CO2 = getRandomNormal2(5.35,0.27,generator() );
                
                //lambda is the equilibrium climate parameter in W m-2 K-1
                //Randomly select values for lambda from the palaeo prior distribution of lambda from log-normal S distribution from Palaeosens
                //double random_palaeosens1 = 0.0;
                //random_palaeosens1 = getRandomLinear2(0.0, 1.0027623323, generator());
                //lambda = get_lambda_palaeosens1(random_palaeosens1);
                
                //Use this function for palaeo prior distribution from standard error S distribution from Palaeosens
                //random_palaeosens1 = getRandomLinear2(0.0, 9.99022973240979E-01, generator());
                //lambda = get_lambda_palaeosens2(random_palaeosens1);
                ///////////////////////////////////////////////////////////////////////////////////
                
                
                ///////////////////////////////////////////////////////////////////////////////////
                
                
                //dNPPdT =  getRandomLinear2(-5.0, 1.0, generator());
                //gamma_K = getRandomLinear2(0, 1.0, generator());
                //dtaudT =  getRandomNormal2(-1.36, 0.45, generator());// mean +/- stdev of DGVMs in Pugh et al (2018)//getRandomLinear2(-2.0, 1.0, generator());
                ratio =    getRandomLinear2(0.2, 1.5, generator());//(0.25,1.1);
                ratio2 =  getRandomLinear2(0.1, 1.0, generator());
                //epsilon = 1.0;//getRandomNormal(1.28,0.25);/////getRandomLinear2(0.83, 1.82, generator());///
                //epsilon_aero = 1.0;//getRandomLinear2(0.33, 3.0, generator());
                
                lambda_Planck = getRandomNormal2(3.3,0.1, generator());// Planck feedback from 4f_greenhouse_sigmaT^3// getRandomNormal(3.15,0.04);        //Planck feedback from CMIP5 models
                
                lambda_WVLR = getRandomLinear2(-3.0, 1.0, generator());// getRandomNormal(-1.15,0.09); //getRandomLinear2(-3.15,2.0, generator());//         //Combined Water Vapour - Lapse Rate
                lambda_Cloud_Fast = 0.0;// Use 0.0 when doing just one getRandomNormal(-0.43,0.33);   //Clouds with original SST patterns
                lambda_Cloud_SST = getRandomLinear2(-3.0, 2.0, generator()); //getRandomNormal(0.0,3.3);//getRandomLinear2(-2.0, 2.0, generator()); //getRandomNormal(-0.47,0.30);    //Clouds with final SST patterns
                lambda_albedo = 0.0;//getRandomNormal(0.0,3.3);//getRandomLinear2(-2.0, 2.0, generator()); //getRandomNormal(-0.37,0.10);      //getRandomLinear2(-3.15,2.0, generator());// Fast snow and ice-albedo
                lambda_1000yr = 0.0;//lambda - (lambda_WVLR + lambda_Cloud_Fast + lambda_Cloud_SST + lambda_albedo + lambda_Planck);
                
                lambda = lambda_Planck + lambda_WVLR + lambda_Cloud_Fast + lambda_albedo + lambda_Cloud_SST + lambda_1000yr;
                
                f_heat_ocean = 0.93;//getRandomLinear2(0.9, 0.96, generator());
                
                tau_C_mixed = getRandomLinear2(0.5, 1.0, generator()); //getRandomLinear2(0.1, 0.5, generator());
                tau_C_upper = getRandomLinear2(5.0, 40.0, generator());
                tau_C_inter = getRandomLinear2(15.0, 60.0, generator());
                tau_C_deep =  getRandomLinear2(100.0, 500.0, generator());
                tau_C_bottom = getRandomLinear2(400.0, 1500.0, generator());
                Ib = getRandomLinear2(3100, 3900, generator());
                
                //Set timescales from random distribution
                tau_WVLR = getRandomNormal2(0.024383562, 0.00109589, generator()); //8.9±0.4 days from residnece time of water vapour in Ent and Tuinenberg (2017)
                
                //If negative value then set to small positive value of ~ 0.4 days/////////
                if(tau_WVLR < 0.0)
                {
                    tau_WVLR = 0.001;
                }
                
                tau_Cloud_Fast = tau_WVLR; //Same as WVLR
                tau_albedo = getRandomLinear2(0.5,1.5, generator());//respopnses on approx. annula timescales //getRandomLinear2(0.5,5.0, generator()); //half a year to 5 years: fast snow feedback to slow multi-year sea-ice feedback. Colman 2013 find multiple timescales from seasonal to decadal.
                tau_Cloud_SST = getRandomLinear2(20.0, 45.0, generator()); //Andres uses first 20 years for fast - giving lower limit. Fine et al find 45 year timescale for tracers going through ventilated thermocline - giving upper limit.
                
                
                if(Future_scenario <=5)
                {
                    //For random normal aerosol etc distribution to get total from aerosol etc
                    //R_aerosol_Uncert = getRandomNormal(-0.14,0.6079);///getRandomNormal(0.074,0.6079);///
                    //Get to -0.35 +-0.45 Wm-2 with standard dev
                    R_aerosol_Uncert = getRandomNormal2(-0.057,0.4298, generator());//
                    //R_WMGHG_nonCO2_Uncert = getRandomNormal(0.0,0.06079);////getRandomNormal(0.32,0.06079); //
                    //Use 0.13 to get to Etminen et al.'s updated CH4 radiative forcing (increase by 0.13 Wm-2 for 2011
                    R_WMGHG_nonCO2_Uncert = getRandomNormal2(0.13,0.06079, generator());////getRandomNormal(0.32,0.06079); //
                }
                /*
                if(Future_scenario >=6 && Future_scenario <15)
                {
                    //For random normal aerosol etc distribution to get total from aerosol etc
                   // R_aerosol_Uncert = getRandomNormal(0.43026,0.6079);///getRandomNormal(0.074,0.6079);///
                    //Get to -0.35 + -0.45 Wm-2 with standard dev
                    R_aerosol_Uncert = getRandomNormal2(0.43026,0.6079, generator());///getRandomNormal(0.074,0.6079);///
                    R_WMGHG_nonCO2_Uncert = getRandomNormal2(-0.03,0.06079, generator());////getRandomNormal(0.32,0.06079); //
                    
                    
                }
                */
                //For uninformed distribution (-0.14 +/- 3 standard deviations)
                //R_aerosol_Uncert = getRandomLinear2(-1.9637,1.6837, generator());
                
                //R_volcanic_coeff = getRandomNormal(-19.0, 0.5); //(-17.0, 1.0);   ////Refers to AOD to RF conversion with error at 1-sigma //AR5 is -25 Wm-2, but Gregory et al (2016) suggest -19 to -17 Wm-2
                //R_volcanic_coeff = getRandomNormal(-19.0, 0.5); //5% stdev uncertainty, or ~+/- 1W/m2 in 20W/m2
                R_volcanic_coeff = getRandomNormal2(-19.0, 0.5, generator());
                //17+/-1.0 is the AMIP period in Gregory et al., 19+/-0.5 is the HadCM3 in Gregory et al.
                
                //Use top line below IF using AOD to calculate RF from volcanoes. Otherwise use second line to introduce similar uncertianty in volcanic RF effect from volcanic RF scenarios
                for(int n=0; n<tmax; n++)
                {
                    if(year[n] >= 1850.0)
                        R_volcanic[n] = R_volcanic_coeff * AOD[n];
                    if(year[n] < 1850.0)
                        R_volcanic[n] = -(R_volcanic_coeff/19.0) * R_volcanic[n] ;
                    //the number 19.0 or 17.0 above needs to fit the value in R_volcanic_coeff
                }
                
                //Using 14% is 2stdev uncertainty in CH4 and 10% is 2stdev uncertainty in N2O after Entiman
                Uncert_CH4 = getRandomNormal2(1.0,0.07, generator());
                Uncert_N2O = getRandomNormal2(1.0,0.05, generator());
                //Uncert Halocarbons is from +/- 10% uncertainty in IPCC 2013 - where 10% is taken as 2stdevs
                Uncert_Halocarbons = getRandomNormal2(1.0,0.05, generator());
                
                
                //Caclulate gamma for each separate aerosol radiative forcing
                //gamma's are in Wm-2 per unit emission - emissions are in same units as Aero_XX_em are defined
                gamma_aero_SOx = getRandomNormal2(-0.31,0.11, generator())/Aero_SOx_em[(2010-Start_year)];
                gamma_aero_BC = getRandomNormal2(0.18,0.07, generator())/Aero_BC_em[(2010-Start_year)];
                gamma_aero_OC = getRandomNormal2(-0.03,0.01, generator())/Aero_OC_em[(2010-Start_year)];
                gamma_aero_ENMVOC = getRandomNormal2(-0.06,0.09, generator())/Aero_ENMVOC_em[(2010-Start_year)];
                gamma_aero_NOx = 0.4*getRandomNormal2(-0.08,0.04, generator())/Aero_NOx_em[(2010-Start_year)];
                gamma_aero_NH3 = 0.6*getRandomNormal2(-0.08,0.04, generator())/Aero_NH3_em[(2010-Start_year)];
                
                AeroRF_ind_2011 = getRandomSkewNormnal2(-0.55, 0.37, -2.0, generator());
               
                ratio_volcanic_CO2=1.0;//getRandomLinear2(0.25,1.25,generator());
                ratioCO2_WVLR[0] = ratio_volcanic_CO2;//ratio_volcanic_CO2;
                ratioCO2_CloudFast[0] = ratio_volcanic_CO2;//ratio_volcanic_CO2;
                ratioCO2_albedo[0] = ratio_volcanic_CO2;//ratio_volcanic_CO2;
                lambda_Cloud_SST_time[0] = ratio_volcanic_CO2;//ratio_volcanic_CO2;
                
            }
            
            
            cout << a_CO2_obs.size() << endl;
            if(Future_scenario>=15 || Use_stored_params == 1)
            {
                
                //Radnom normal value for radiative forcing coefficient
                a_CO2 = a_CO2_obs[j];
                
                                
                dNPPdT =  dNPPdT_obs[j];
                gamma_K = gamma_K_obs[j];
                dtaudT =  dtaudT_obs[j];// mean +/- stdev of DGVMs in Pugh et al (2018)//getRandomLinear2(-2.0, 1.0, generator());
                ratio =   ratio_obs[j];//(0.25,1.1);
                ratio2 =  ratio2_obs[j];
                //epsilon = 1.0;//getRandomNormal(1.28,0.25);/////getRandomLinear2(0.83, 1.82, generator());///
                //epsilon_aero = 1.0;//getRandomLinear2(0.33, 3.0, generator());
                
                
                lambda_Planck =  lambda_Planck_obs[j];        //Planck feedback
                lambda_WVLR =  lambda_WVLR_obs[j];         //Combined Water Vapour - Lapse Rate
                lambda_Cloud_Fast = lambda_Cloud_Fast_obs[j];   //Clouds with original SST patterns
                lambda_Cloud_SST = lambda_Cloud_SST_obs[j];    //Clouds with final SST patterns
                lambda_albedo = lambda_albedo_obs[j];      //Fast snow and ice-albedo
                lambda_1000yr = 0.0;//lambda - (lambda_WVLR + lambda_Cloud_Fast + lambda_Cloud_SST + lambda_albedo + lambda_Planck);
                lambda = lambda_Planck + lambda_WVLR + lambda_Cloud_Fast + lambda_albedo + lambda_Cloud_SST + lambda_1000yr;
                
                f_heat_ocean = f_heat_ocean_obs[j];
                
                tau_C_mixed = tau_C_mixed_obs[j]; //getRandomLinear2(0.1, 0.5, generator());
                tau_C_upper = tau_C_upper_obs[j];
                tau_C_inter = tau_C_inter_obs[j];
                tau_C_deep =  tau_C_deep_obs[j];
                tau_C_bottom = tau_C_bottom_obs[j];
                Ib = Ib_obs[j];
                
                //Set timescales from random distribution
                tau_WVLR = tau_WVLR_obs[j];; //8.9±0.4 days from residnece time of water vapour in Ent and Tuinenberg (2017)
                
                //If negative value then set to small positive value of ~ 0.4 days/////////
                if(tau_WVLR < 0.0)
                {
                    tau_WVLR = 0.001;
                }
                
                tau_Cloud_Fast = tau_Cloud_Fast_obs[j]; //Same as WVLR
                tau_albedo = tau_albedo_obs[j]; //half a year to 5 years: fast snow feedback to slow multi-year sea-ice feedback. Colman 2013 find multiple timescales from seasonal to decadal.
                tau_Cloud_SST = tau_Cloud_SST_obs[j]; //Andres uses first 20 years for fast - giving lower limit. Fine et al find 45 year timescale for tracers going through ventilated thermocline - giving upper limit.
                
                R_aerosol_Uncert = R_aerosol_Uncert_obs[j];
                R_WMGHG_nonCO2_Uncert = R_WMGHG_nonCO2_Uncert_obs[j];
                R_volcanic_coeff = R_volcanic_coeff_obs[j];
                
                //Using 14% is 2stdev uncertainty in CH4 and 10% is 2stdev uncertainty in N2O after Entiman
                Uncert_CH4 = Uncert_CH4_obs[j];
                Uncert_N2O = Uncert_N2O_obs[j];
                //Uncert Halocarbons is from +/- 10% uncertainty in IPCC 2013 - where 10% is taken as 2stdevs
                Uncert_Halocarbons = Uncert_Halocarbons_obs[j];
                
                
                //Caclulate gamma for each separate aerosol radiative forcing
                //gamma's are in Wm-2 per unit emission - emissions are in same units as Aero_XX_em are defined
                gamma_aero_SOx = gamma_aero_SOx_obs[j];
                gamma_aero_BC = gamma_aero_BC_obs[j];
                gamma_aero_OC = gamma_aero_OC_obs[j];
                gamma_aero_ENMVOC = gamma_aero_ENMVOC_obs[j];
                gamma_aero_NOx = gamma_aero_NOx_obs[j];
                gamma_aero_NH3 = gamma_aero_NH3_obs[j];
                
                AeroRF_ind_2011 = AeroRF_ind_2011_obs[j];
                
            }
            
            
            
            
            
            
            
            
            
            
            
            
            
            //If use sea level component then initalise random sea level coefficients here//////////////////
            //Sea level component described in Goodwin et al, Earth's Future 2017//////////////////////////
            
            //steric_coeff = getRandomNormal(1.5,0.3);//Random normal inputs for Williams et al (2012) in GRL approach to steric sea level rise (see Goodwin et al, Earth's Future 2017
            
            //ice_coeff = getRandomLinear(0.0,0.005); //For T, Rahmstorf (2007) coefficient in units mm K-1 yr-1
            
            //land_water = getRandomLinear(-0.01,0.09);
            
            
            //Now calculate the effective climate sensitivity of sea surface temperatures (in K [Wm-2]-1)
            S_ocean = (1.0/lambda)*ratio; //ocean climate sensitivity (K [Wm-2]-1)
            
            
            ////////////////////////////////////////////////////////////////
            
        }
        
        
        
        
        for (int i = 0; i<tmax; i++) //Initialise heat content changes in upper and deep ocean
        {
            
            
            Heat_upper[i] = 0.0;
            Heat_inter[i] = 0.0;
            Heat_deep[i] = 0.0;
            Heat_mixed[i] = 0.0;
            Heat_bottom[i] = 0.0;
            
            I_eq_heat[0] = 0.0; //Initialise equivalent emissions from heat content change
            
        }
        
        if(RESTORE == 1) //Need to initialise emission profile to 0.0 if using a restoring CO2 pathway
        {
            for(int i = 0; i<tmax; i++)
            {
                I_em[i] = 0.0;
                CO2_restore[i] = CO2[i];
                CO2[i] = 0.0;
            }
            CO2[0] = CO2init;
        }
        
        //Initialise lambda feedbacks from different timescales
        for(int n=0; n<n_forcings; n++)
        {
            lambda_WVLR_time[n] = 0.0;
            lambda_Cloud_Fast_time[n] = 0.0;
            lambda_Cloud_SST_time[n] = 0.0;
            lambda_albedo_time[n] = 0.0;
            lambda_1000yr_time[n] = 0.0;
            
        }
        
        overshoot = 0;
        int i_last_check_year = 0;
        
        for (int i = 1; i<tmax; i++) 
             
        {
            
            //Define the internal radiative budget noise
            R_noise[i] = R_noise[i-1]*gamma_Rnoise + alpha_Rnoise*getRandomLinear(-1.0, 1.0);
            /*
            if(Future_scenario>= 15)
            {
                R_noise[i] = 0.0;
            }
            */
            //This counts on the time step and year counter
            t_step_count +=1;
            if (t_step_count == t_step_per_year)
            {
                t_step_count = 0;
                year_count += 1;
            }
            
            
            //check if lambda is physically plausible - must be positive for all timescales and give ECS less than 20
            if( i == 1 && (lambda_Planck + lambda_WVLR + lambda_Cloud_Fast  < (a_CO2 * log(2.0) / 20.0) || lambda_Planck +lambda_WVLR + lambda_Cloud_Fast + lambda_Cloud_SST < (a_CO2 * log(2.0) / 20.0) || lambda <  (a_CO2 * log(2.0) / 20.0) ) )
            {
                i = tmax - 1;
            }
            
            //Scaling for non-CO2 forcing agents, scaled by uncertainties in 2011
            if(Future_scenario < 15)
            {
                
                //First non-Kyoto and nonCO2 (e.g. aerosols)
                //R_aerosol[i] = R_aerosol[i]*(1.0+(R_aerosol_Uncert/R_aerosol_2011));
                
                //cout << year[i] << '\t' << R_aerosol[i]  << '\t' << getR_aerosol(year_count, t_step_count, gamma_aero_BC, gamma_aero_OC, gamma_aero_SOx, gamma_aero_NOx, gamma_aero_NH3, gamma_aero_ENMVOC, Start_year, AeroRF_ind_2011) << endl;
           
                //Aerosol radiative forcing from separate aerosol emission components
                R_aerosol[i] =getR_aerosol(year_count, t_step_count, gamma_aero_BC, gamma_aero_OC, gamma_aero_SOx, gamma_aero_NOx, gamma_aero_NH3, gamma_aero_ENMVOC, Start_year, AeroRF_ind_2011);
                //Now Kyoto but nonCO2 (e.g. CH4, N2O and halogens)
                //R_WMGHG_nonCO2[i] = R_WMGHG_nonCO2[i]*(1.0 + (R_WMGHG_nonCO2_Uncert/R_WMGHG_nonCO2_2011) );
                
                R_WMGHG_nonCO2[i] = getRnonCO2(year_count, i, Uncert_N2O, Uncert_CH4, Uncert_Halocarbons, Start_year, PgCtoppm);
                //cout << CH4.size() << '\t' << R_WMGHG_nonCO2[i] << '\t' << year_count << '\t' << i << '\t' << Uncert_N2O << '\t' << Uncert_CH4 << '\t' << Start_year << '\t' << PgCtoppm <<  endl;
                
                //Alternative getR_aerosol and getR_WMGHG_nonCO2 functions here
                //R_aerosol[i] = getR_aerosol(year_count, t_step_count, coefficients);
                //R_WMHGH_nonCO2[i] = getR_WMGHG_nonCO2(year_count, t_step_count, coefficients);
            }
            
            
            
            //NOTE: if have internal noise, need to re-tune since volcanic eruptions are now contributing to total noise.
            //Noise in surface warming///////////////////////////////////////////////////
            z_noise = getRandomLinear(-1.0, 1.0);//0.0;
            
            if(i<=2*t_step_per_year)
            {
                T_noise[i] = (gamma1_noise * T_noise[i-1]) + (alpha_noise * z_noise);
                
            }
            if(i>2*t_step_per_year)
            {
                T_noise[i] = (gamma1_noise * T_noise[i-1]) + (gamma2_noise*T_noise[i-2]) + (alpha_noise*z_noise);
            }
            T_noise[i] = 0.0;
            /////////////////////////////////////////////////////////////////////////////
            //No internal noise for CO2 only scenarios
            if(Future_scenario >= 15)
            {
                T_noise[i] = 0.0;
            }
        
            //Carbon component////////////////////////////////////////////////////
            
            //Terrestrial Carbon Component////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            
            dIvegdt = getdCvdt(CO2[i-1]*PgCtoppm, CO2init*PgCtoppm, DeltaT[i-1], I_veg[i-1], I_veg_init, NPP_init, gamma_K, dNPPdT);
            
            dIsoildt = getdCsdt(I_veg[i-1], I_veg_init, I_soil[i-1], I_soil_init, NPP_init, DeltaT[i-1], dtaudT);
            
            
             
            I_veg[i] = I_veg[i-1] + dIvegdt*dt;
            
            I_soil[i] = I_soil[i-1] + dIsoildt*dt;
            
            //Set minimum I_veg and I_soil carbon inventories of 10 PgC
            if(I_veg[i] < 10.0)
                I_veg[i] = 10.0;
            if(I_soil[i] < 10.0)
                I_soil[i] = 10.0;
            
            I_ter[i] = I_veg[i] + I_soil[i];
            
            
            
            //Ocean carbon component /////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            
            
            //The temperature - CO2 solubility feedback (after Goodwin and Lenton 2009 in GRL)/////////////////////////////////
            dCsatdT = -0.1 + 0.025*(log(CO2[i-1]/CO2init));
            Em_heat = (-dCsatdT/1.0e15)*f_heat_ocean*N_mixed[i-1]*area*(60.0*60.0*24.0*365.0)*dt/c_p;
            I_eq_heat[i] = I_eq_heat[i-1] + Em_heat;
            
            ////////////////////////////////////////////////////////////////////////////
            
            
           //Calculate emission to restore towards a given CO2 curve (e.g. from the RCP scenarios)
            I_em[i] = I_em[i-1] + 0.9*( CO2_restore[i] - CO2[i-1] ) + I_Usat_mixed[i-1]*(1.0-exp(-(dt/tau_C_mixed))) ;
                
                //Note: 0.9 coefficient in above equation is needed for stability - if you do not have it then you can get I_em in an unstable oscillation about CO2_restore[i]
            
            
            
            
            
            //Adjusting Mitigation Pathways code////////////////////////////
            /////////////////////////////////////////////////////////////////
            
            //Start of AMP scenarios, linear reductions to stabilise climate at given level of surface warming.
            
            if(Future_scenario == 5 && year[i] > Restore_year - (dt/2.0) )
            {
                
                //New scheme////////////////////////////////////////////
                
                double EmRate =0.0; //Annual carbon emission rate in PgC yr-1
                double EmRate_assess = 0.0; //Annual carbon emissions rate at the Assess year
                
                //Total annual emissions from land use change and fossil fuels = (I_em[i] + I_ter[i]) - (I_em[i-1] + I_ter[i-1])
                
                
                
                EmRate = getRateINDCs2017(year[i]);
                
               
                
                if( year[i] > (Restore_year-(dt/2.0)) ) //Estimate emissions needed to
                {
                    
                    
                    
                    

                    
                    //
                    if( i == (2030-Start_year)*t_step_per_year || i == (2040-Start_year)*t_step_per_year || i == (2050-Start_year)*t_step_per_year || i == (2060-Start_year)*t_step_per_year || i == (2070-Start_year)*t_step_per_year || i == (2080-Start_year)*t_step_per_year || i == (2090-Start_year)*t_step_per_year || i == (2100-Start_year)*t_step_per_year || i == (2125-Start_year)*t_step_per_year )
                    {
                        i_last_check_year = i;
                        
                        //Calculate average temperatures
                        double DT_ave_1961_1990 = 0.0;
                        double DT_ave_2003_2012 = 0.0;
                        double DT_ave2=0.0;
                        double Emissions_1975_2025= 460.6; //Estimate from Carbon Budget plus INDC
                        double Emissions_2008_minus5years=0.0;

                        
                        if(i == (2030-Start_year)*t_step_per_year)
                        {
                            EmRate_Assessment = getRateINDCs(2030.0);
                        }
                        
                        if(i > (2030-Start_year)*t_step_per_year)
                        {
                            //Rate of emission over last 12 months
                            EmRate_Assessment = ( (I_em[i-1]+I_ter[i-1])-(I_em[i-t_step_per_year-1]+I_ter[i-t_step_per_year-1]) ) ;// ( (I_em[i-1]+I_ter[i-1])-(I_em[i-2]+I_ter[i-2]) )* double(t_step_per_year);
                        }
                        
                        
                            
                            
                        for(int n = (1961-Start_year)*t_step_per_year; n< (1991-Start_year)*t_step_per_year; n++)
                        {
                            DT_ave_1961_1990 += DeltaT[n];
                        }
                        DT_ave_1961_1990 = DT_ave_1961_1990/(30.0*double(t_step_per_year));
                            
                            
                            
                        for(int n = (2003-Start_year)*t_step_per_year; n< (2013-Start_year)*t_step_per_year; n++)
                        {
                            DT_ave_2003_2012 += DeltaT[n];
                        }
                        DT_ave_2003_2012 = DT_ave_2003_2012/(10.0*double(t_step_per_year));
                            
                            
                            
                            
                        for(int n = 0; n< 10*t_step_per_year; n++)
                        {
                            DT_ave2 += DeltaT[i-11*t_step_per_year+n];
                        }
                        DT_ave2 = DT_ave2/(10.0*double(t_step_per_year));
                            
                        Emissions_2008_minus5years = (I_em[i-5*t_step_per_year] + I_ter[i-5*t_step_per_year]) - (I_em[(2008-Start_year)*t_step_per_year] + I_ter[(2008-Start_year)*t_step_per_year]);
                            
                        //Estimate TCRE from 2003_2012 and ten year period prior to Assess_year1
                        TCRE_estimate = (DT_ave2 - DT_ave_2003_2012)/Emissions_2008_minus5years + TCRE_safety;
                            
                            
                        //Estimate remaining allowable emissions////////////////////////////////
                        if(Target1 > DT_ave2 - DT_ave_2003_2012 + 0.78 + Coeff_Target)
                        {
                            I_remaining = ((Target1-0.78-Coeff_Target) - (DT_ave2-DT_ave_2003_2012) )/TCRE_estimate;
                        }
                        if(Target1 <= DT_ave2 - DT_ave_2003_2012 + 0.78 + Coeff_Target)
                        {
                            I_remaining = 0.0;
                            overshoot = 1;
                        }
                            
                        //Estimate phase-out time
                        t_remaining = 2.0*I_remaining/EmRate_Assessment;
                            
                            
                        cout << year[i] << '\t' << TCRE_estimate << '\t' << DT_ave2 - DT_ave_2003_2012 + 0.78 << '\t' << t_remaining << '\t' << I_remaining << '\t' << DT_ave_2003_2012 << '\t' << DT_ave2 << '\t' << Emissions_2008_minus5years << '\t' << (I_em[i-5*t_step_per_year] + I_ter[i-5*t_step_per_year]) << '\t' << I_em[(2008-Start_year)*t_step_per_year] + I_ter[(2008-Start_year)*t_step_per_year] << endl;
                            
                        /*if(data_fit>=16)
                        {
                            Adjustment_TCRE << j << '\t' << TCRE_estimate;
                            Adjustment_DT << j << '\t' << DT_ave2 - DT_ave_2003_2012 + 0.78;
                            if(I_remaining > 0.0)
                            {
                                Adjustment_Emtime << j << '\t' << year[i]+t_remaining ;
                                    
                            }
                            Adjustment_Emrate << j << '\t' << EmRate_Assessment;
                        }*/
                    }
                        
                    //Calculate emissions between assessment points 1 and 2//////////////////////
                    /////////////////////////////////////////////////////////////////////////////
                    if( year[i] > 2030.0 && year[i] < year[i_last_check_year]+t_remaining && year[i]<=2150.0)
                    {
                            
                        EmRate = EmRate_Assessment*(1.0 - (year[i] - year[i_last_check_year])/(t_remaining));
                            //cout << EmRate_assess << '\t' << EmRate << endl;
                    }
                    ////////////////////////////////////////////////////////
                        
                    if(year[i] > 2030.0 && year[i] >= year[i_last_check_year]+t_remaining && year[i] <= 2150.0)
                    {
                        EmRate = 0.0;
                    }
                        
                
                
                //Below ammends if over warming target////////
                
                //From 2150, add in carbon if below warming target, and remove carbon if above target
                if( ( (year[i] >= 2030.0 && overshoot == 1) || year[i] >= 2150.0 ) && year[i] < 2300.0)
                {
                    double DT_ave_2003_2012 = 0.0;
                    double DT_ave2 = 0.0;
                    
                    for(int n = 0; n< 10*t_step_per_year; n++)
                    {
                        DT_ave_2003_2012 += DeltaT[(2003-Start_year)*t_step_per_year+n];
                    }
                    DT_ave_2003_2012 = DT_ave_2003_2012/(10.0*double(t_step_per_year));
                    
                    //Last 10 years (11 years ago to 1 year ago///////////////////
                    for(int n = 0; n< 10*t_step_per_year; n++)
                    {
                        DT_ave2 += DeltaT[i-11*t_step_per_year+n];
                    }
                    DT_ave2 = DT_ave2/(10.0*double(t_step_per_year));
                    
                    
                    if( DT_ave2 - DT_ave_2003_2012 + 0.78 > Target1+0.00)
                    {
                        EmRate = -1.0;
                        //cout << "negative  emissions" << endl;
                        overshoot=1;
                    }
                    
                    if( DT_ave2 - DT_ave_2003_2012 + 0.78 > Target1+0.05)
                    {
                        EmRate = -2.0;
                        //cout << "negative  emissions" << endl;
                        overshoot=1;
                    }
                    
                    if( DT_ave2 - DT_ave_2003_2012 + 0.78 > Target1+0.1)
                    {
                        EmRate = -3.0;
                        //cout << "negative  emissions" << endl;
                        overshoot=1;
                    }
                    
                    if(EmRate < 0.01 && DT_ave2 - DT_ave_2003_2012 + 0.78 < Target1 - 0.2)
                    {
                        EmRate = 1.0;
                        //cout << "negative  emissions" << endl;
                        
                    }
                }
                
                    
                    I_em[i] = (I_em[i-1] + I_ter[i-1] - I_ter[i]) + EmRate*dt;
                    
                }
            
            }
            
                
            
                
                
                
        
        
            ///End of AMP scenario (linear reductions)
            ///////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////
            
            
            
            
            
                //Now calculate the undersaturation of DIC in the ocean with respect to instantaneous CO2 concentration (after Goodwin, 2016).
            
                //See Goodwin et al (2015) in Nature Geoscience for description of I_Usat.
            
            
                //Following adds in new I_Usat if CO2 rises, advects I_Usat between ocean boxes and reduces I_Usat due to air-sea carbon fluxes. see Goodwin (2016) in Climate Dynamics for description of how this works.
                
                I_Usat_MixToUpper =    (I_Usat_upper[i-1] - (I_Usat_mixed[i-1]*vol_upper/vol_mixed)  )*exp(-(dt/(tau_C_upper )))+ (I_Usat_mixed[i-1]*vol_upper/vol_mixed) - I_Usat_upper[i-1]   ;
                
                I_Usat_MixToInter =    (I_Usat_inter[i-1] - (I_Usat_mixed[i-1]*vol_inter/vol_mixed)  )*exp(-(dt/(tau_C_inter )))+ (I_Usat_mixed[i-1]*vol_inter/vol_mixed) - I_Usat_inter[i-1]   ;
                
                I_Usat_MixToDeep =     (I_Usat_deep[i-1] - (I_Usat_mixed[i-1]*vol_deep /vol_mixed)  )*exp(-(dt/tau_C_deep))   +  (I_Usat_mixed[i-1]*vol_deep/vol_mixed) - I_Usat_deep[i-1]   ; 
                
                I_Usat_MixToBottom =   (I_Usat_bottom[i-1]-(I_Usat_mixed[i-1]*vol_bottom/vol_mixed) )*exp(-(dt/tau_C_bottom)) +  (I_Usat_mixed[i-1]*vol_bottom/vol_mixed) - I_Usat_bottom[i-1]  ;
                
                ////////////////////////////////////////////////////////////////////////////
                
                //I_Usat_mixed tends towards zero with an e-folding timescale of tau_C_mixed, and there are fluxes of I_Usat into the sub-surface regions.
                
                I_Usat_mixed[i] = I_Usat_mixed[i-1]*exp(-(dt/tau_C_mixed))     //The exponential decay towards zero with air-sea exchange
                - I_Usat_MixToUpper                               //exchange of I_Usat with upper
                - I_Usat_MixToInter                               //exchange of I_Usat with intermediate
                - I_Usat_MixToDeep                                //exchange of I_Usat with deep
                - I_Usat_MixToBottom                              //exchange of I_Usat with bottom
                + vol_mixed* (Ib* log( (CO2[i-1] + (I_em[i]-I_em[i-1] + Em_heat)  )/CO2[i-1]) - (I_em[i]-I_em[i-1]+Em_heat)) ;   //The addition of carbon into the atmosphere increasing global I_Usat, partitioned by volume
                
                
                I_Usat_upper[i] = I_Usat_upper[i-1] + I_Usat_MixToUpper     //Exchange of I_Usat with mixed layer
                          +  vol_upper*(Ib* log( (CO2[i-1] + (I_em[i]-I_em[i-1] + Em_heat) )/CO2[i-1]) - (I_em[i]-I_em[i-1]+Em_heat) - I_Usat_mixed[i-1]*(1.0-exp(-(dt/tau_C_mixed)) )) ;                          //volume weighted contribution of carbon emissions increasing global I_Usat
                
                I_Usat_inter[i] = I_Usat_inter[i-1] + I_Usat_MixToInter     //Exchange of I_Usat with mixed layer
                          +  vol_inter*(Ib* log( (CO2[i-1] + (I_em[i]-I_em[i-1] + Em_heat) )/CO2[i-1]) - (I_em[i]-I_em[i-1]+Em_heat) - I_Usat_mixed[i-1]*(1.0-exp(-(dt/tau_C_mixed)) )) ;     
            
                I_Usat_deep[i] = I_Usat_deep[i-1]  + I_Usat_MixToDeep       //Exchange of I_Usat with mixed layer
                                        +  vol_deep*(Ib*  log( (CO2[i-1] + (I_em[i]-I_em[i-1] + Em_heat) )/CO2[i-1]) - (I_em[i]-I_em[i-1]+Em_heat) - I_Usat_mixed[i-1]*(1.0-exp(-(dt/tau_C_mixed)) )) ;                          //volume weighted contribution of carbon emissions increasing global I_Usat
                
                I_Usat_bottom[i] = I_Usat_bottom[i-1] + I_Usat_MixToBottom  //Exchange of I_Usat with mixed layer
                                                     +  vol_bottom*(Ib*  log( (CO2[i-1] + (I_em[i]-I_em[i-1] + Em_heat) )/CO2[i-1]) - (I_em[i]-I_em[i-1]+Em_heat) - I_Usat_mixed[i-1]*(1.0-exp(-(dt/tau_C_mixed)) )) ;                        //volume weighted contribution of carbon emissions increasing global I_Usat
            
            
            
                //Caclulation of atmospheric CO2 using I_Usat, I_em and Ib after Goodwin et al (2015) in Nature Geoscience, also see Goodwin (2016) in climate dynamics
                
                CO2[i] = CO2init*exp( (I_em[i] + I_Usat_mixed[i] + I_Usat_upper[i] + I_Usat_inter[i] + I_Usat_deep[i] + I_Usat_bottom[i] + I_eq_heat[i])/Ib);
                
                
                //////////////////////////////////////////////////////////////////////
            
            //set lambda equal to lambda_Planck for first few timesteps
            if(i <= 2)
            {
                for(int n =0; n<n_forcings; n++)
                {
                    lambda_time[n][i] = lambda_Planck;
                    //S_ocean = (1.0/lambda_time[i])*ratio; //ocean climate sensitivity (K [Wm-2]-1)
                    
                    
                }
                
            }
            //Independently time-varying lambda's for each source of radiative forcing
            if(i > 2)
            {
                
                //Now modify lambda from different timescale-processes by exponential decay from current value towards equilibrium value
                
                
                for(int n = 0; n<n_forcings; n++)
                {
                    lambda_WVLR_time[n] = lambda_WVLR_time[n]     + (1.0-exp(-dt/tau_WVLR))*( (ratioCO2_WVLR[n]*lambda_WVLR)-lambda_WVLR_time[n]);
                    lambda_Cloud_Fast_time[n] = lambda_Cloud_Fast_time[n]     + (1.0-exp(-dt/tau_Cloud_Fast))*( (ratioCO2_CloudFast[n]*lambda_Cloud_Fast)-lambda_Cloud_Fast_time[n]);
                    lambda_Cloud_SST_time[n] = lambda_Cloud_SST_time[n]   + (1.0-exp(-dt/tau_Cloud_SST))*((ratioCO2_CloudSST[n]*lambda_Cloud_SST)-lambda_Cloud_SST_time[n]);
                    lambda_albedo_time[n] = lambda_albedo_time[n]   + (1.0-exp(-dt/tau_albedo))*((ratioCO2_albedo[n]*lambda_albedo)-lambda_albedo_time[n]);
                    lambda_1000yr_time[n] = lambda_1000yr_time[n] + (1.0-exp(-dt/1000.0))*((ratioCO2_1000yr[n]*lambda_1000yr)-lambda_1000yr_time[n]);
                }
                
                //New radiative forcing further from R=zero contributes at lambda_Planck only, so modify down the lambda_process_time terms
                if(abs(R_volcanic[i]) - abs(R_volcanic[i-1]) > 0.0)
                {
                    lambda_WVLR_time[0] = lambda_WVLR_time[0]*(abs( R_volcanic[i-1] /R_volcanic[i]  ));
                    lambda_Cloud_Fast_time[0] = lambda_Cloud_Fast_time[0]*(abs( R_volcanic[i-1] /R_volcanic[i]  ) );
                    lambda_Cloud_SST_time[0] = lambda_Cloud_SST_time[0]*(abs( R_volcanic[i-1] /R_volcanic[i]  ) );
                    lambda_albedo_time[0] = lambda_albedo_time[0]*(abs( R_volcanic[i-1] /R_volcanic[i]  ) );
                    lambda_1000yr_time[0] = lambda_1000yr_time[0]*(abs( R_volcanic[i-1] /R_volcanic[i]  ) );
                }
                if(abs(R_solar[i]) - abs(R_solar[i-1]) > 0.0)
                {
                    lambda_WVLR_time[1] = lambda_WVLR_time[1]*(abs( R_solar[i-1] /R_solar[i]  ));
                    lambda_Cloud_Fast_time[1] = lambda_Cloud_Fast_time[1]*(abs( R_solar[i-1] /R_solar[i]  ) );
                    lambda_Cloud_SST_time[1] = lambda_Cloud_SST_time[1]*(abs( R_solar[i-1] /R_solar[i]  ) );
                    lambda_albedo_time[1] = lambda_albedo_time[1]*(abs( R_solar[i-1] /R_solar[i]  ) );
                    lambda_1000yr_time[1] = lambda_1000yr_time[1]*(abs( R_solar[i-1] /R_solar[i]  ) );
                }
                if( log(CO2[i]) - log(CO2[i-1]) > 0.0)
                {
                    lambda_WVLR_time[2] = lambda_WVLR_time[2]*(abs( log(CO2[i-1]/CO2init) /log(CO2[i]/CO2init)  ));
                    lambda_Cloud_Fast_time[2] = lambda_Cloud_Fast_time[2]*(abs( log(CO2[i-1]/CO2init) /log(CO2[i]/CO2init)  ) );
                    lambda_Cloud_SST_time[2] = lambda_Cloud_SST_time[2]*(abs( log(CO2[i-1]/CO2init) /log(CO2[i]/CO2init)  ) );
                    lambda_albedo_time[2] = lambda_albedo_time[2]*(abs( log(CO2[i-1]/CO2init) /log(CO2[i]/CO2init)  ) );
                    lambda_1000yr_time[2] = lambda_1000yr_time[2]*(abs( log(CO2[i-1]/CO2init) /log(CO2[i]/CO2init)  ) );
                }
                if(abs(R_aerosol[i]) - abs(R_aerosol[i-1]) > 0.0)
                {
                    lambda_WVLR_time[3] = lambda_WVLR_time[3]*(abs( R_aerosol[i-1] /R_aerosol[i]  ));
                    lambda_Cloud_Fast_time[3] = lambda_Cloud_Fast_time[3]*(abs( R_aerosol[i-1] /R_aerosol[i]  ) );
                    lambda_Cloud_SST_time[3] = lambda_Cloud_SST_time[3]*(abs( R_aerosol[i-1] /R_aerosol[i]  ) );
                    lambda_albedo_time[3] = lambda_albedo_time[3]*(abs( R_aerosol[i-1] /R_aerosol[i]  ) );
                    lambda_1000yr_time[3] = lambda_1000yr_time[3]*(abs( R_aerosol[i-1] /R_aerosol[i]  ) );
                }
                if(abs(R_WMGHG_nonCO2[i]) - abs(R_WMGHG_nonCO2[i-1]) > 0.0)
                {
                    lambda_WVLR_time[4] = lambda_WVLR_time[4]*(abs( R_WMGHG_nonCO2[i-1] /R_WMGHG_nonCO2[i]  ));
                    lambda_Cloud_Fast_time[4] = lambda_Cloud_Fast_time[4]*(abs( R_WMGHG_nonCO2[i-1] /R_WMGHG_nonCO2[i]  ) );
                    lambda_Cloud_SST_time[4] = lambda_Cloud_SST_time[4]*(abs( R_WMGHG_nonCO2[i-1] /R_WMGHG_nonCO2[i]  ) );
                    lambda_albedo_time[4] = lambda_albedo_time[4]*(abs( R_WMGHG_nonCO2[i-1] /R_WMGHG_nonCO2[i]  ) );
                    lambda_1000yr_time[4] = lambda_1000yr_time[4]*(abs( R_WMGHG_nonCO2[i-1] /R_WMGHG_nonCO2[i]  ) );
                }
                if(abs(R_noise[i]) - abs(R_noise[i-1]) > 0.0)
                {
                    lambda_WVLR_time[5] = lambda_WVLR_time[5]*(abs( R_noise[i-1] /R_noise[i]  ));
                    lambda_Cloud_Fast_time[5] = lambda_Cloud_Fast_time[5]*(abs( R_noise[i-1] /R_noise[i]  ) );
                    lambda_Cloud_SST_time[5] = lambda_Cloud_SST_time[5]*(abs( R_noise[i-1] /R_noise[i]  ) );
                    lambda_albedo_time[5] = lambda_albedo_time[5]*(abs( R_noise[i-1] /R_noise[i]  ) );
                    lambda_1000yr_time[4] = lambda_1000yr_time[5]*(abs( R_noise[i-1] /R_noise[i]  ) );
                }
                
                
                for(int n = 0; n<n_forcings; n++)
                {
                    //Sum up feedback contributions to calculate total lambda_time
                    lambda_time[n][i] = lambda_Planck + lambda_WVLR_time[n] + lambda_Cloud_Fast_time[n] + lambda_Cloud_SST_time[n] + lambda_albedo_time[n] + lambda_1000yr_time[n];
                }
                
                
                //S_ocean = (1.0/lambda_time[i])*ratio; //ocean climate sensitivity (K [Wm-2]-1)
                
                
                
                
            }
            
            ////////////////////////////////////////////////////////////////////////////////////////////////////
        
            
            
            //Heat component//////////////////////////////////////////////////////
            
            //The following code advects heat content anomaly round the ocean, and heats the surface ocean due to radiative forcing
            
            //Calculates the equilibrium sea surface heat content change for the instantaneous efficacy-weighted radiative forcing
            //New, time dependent
            H_mixed_equil = ratio*(vol_total * vol_mixed) * c_p *( ( (a_CO2*(log(CO2[i]/CO2init)))/lambda_time[2][i] ) + (R_aerosol[i]/lambda_time[3][i]) + (R_WMGHG_nonCO2[i]/lambda_time[4][i]) + (R_volcanic[i]/lambda_time[0][i]) + (R_solar[i]/lambda_time[1][i]) + (R_noise[i]/lambda_time[1][i])  );
            
            
            
            double Heat_mixed_to_upper =  + (Heat_upper[i-1] - (ratio2*Heat_mixed[i-1]*vol_upper/vol_mixed) )*exp(-(dt/((tau_C_upper)))) + (ratio2*Heat_mixed[i-1]*vol_upper/vol_mixed) -  Heat_upper[i-1] ;
            
            double Heat_mixed_to_inter =  + (Heat_inter[i-1] - (ratio2*Heat_mixed[i-1]*vol_inter/vol_mixed) )*exp(-(dt/((tau_C_inter)))) + (ratio2*Heat_mixed[i-1]*vol_inter/vol_mixed) -  Heat_inter[i-1] ;
            
            double Heat_mixed_to_deep =  + (Heat_deep[i-1] - (ratio2*Heat_mixed[i-1]*vol_deep/vol_mixed) )*exp(-(dt/(tau_C_deep))) +  (ratio2*Heat_mixed[i-1]*vol_deep/vol_mixed) - Heat_deep[i-1] ;
            
            double Heat_mixed_to_bottom = + (Heat_bottom[i-1] - (ratio2*Heat_mixed[i-1]*vol_bottom/vol_mixed) )*exp(-(dt/(tau_C_bottom))) +  (ratio2*Heat_mixed[i-1]*vol_bottom/vol_mixed) - Heat_bottom[i-1]  ; 
            
            
            //Caclulate heat imbalance from distance from equilibrium SST due to efficacy weighted radiative forcing
            //Heat imbalance = total(non-efficacy weighted) radiative forcing times fraction of sea surface heat content anomaly towards equilibrium (See Goodwin, 2016 in Climate Dynamics)
            
            double surface_heat_ratio;
            
            surface_heat_ratio = ((H_mixed_equil - Heat_mixed[i-1])/(H_mixed_equil));
            
            
            if(surface_heat_ratio > 15.0)
            {
                
                surface_heat_ratio = 15.0; //prevents from going too negative after eruption
            }
            if(surface_heat_ratio < -15.0)
            {
                
                surface_heat_ratio = -  15.0; //prevents from going too positive after eruption
            }
            N_mixed[i] =  ( a_CO2*(log(CO2[i]/CO2init)) + R_aerosol[i] + R_WMGHG_nonCO2[i] + R_solar[i] + R_volcanic[i] + R_noise[i] ) * surface_heat_ratio ;
            
            

            //Check to see if N_mixed is outside +/- 20W/m2
            if(N_mixed[i] > 20.0)
            {
                N_mixed[i] = 20.0;
            }
            if(N_mixed[i] < -20.0)
            {
                N_mixed[i] = -20.0;
            }
            
        
            
            

            
            //New, remove divide by H_mixed_equil as when total RF goes near zero this gives very large N
            //Now, N_mixed = DeltaH_mix / (area * dt), where DeltaH_mix is solved via an equilibrium restoring timescale and does not include the fluxes below
           // cout << N_mixed[i] << endl;
            
            //Insert f_heat_ocean fraction of energy imbalance into surface ocean mixed layer//////////////////////////
            Heat_mixed[i] = Heat_mixed[i-1] + (f_heat_ocean*N_mixed[i] * area*60.0*60.0*24.0*365.0*dt) + ((T_noise[i] - T_noise[i-1])*vol_total*vol_mixed*c_p)
            - Heat_mixed_to_upper - Heat_mixed_to_inter - Heat_mixed_to_deep - Heat_mixed_to_bottom ;
            
            
            
            Heat_upper[i] = Heat_upper[i-1] + Heat_mixed_to_upper ;   //Exchange of Heat with mixed layer
            
            Heat_inter[i] = Heat_inter[i-1] + Heat_mixed_to_inter ; 
            
            Heat_deep[i] = Heat_deep[i-1] + Heat_mixed_to_deep;  //Exchange of Heat with mixed layer
            
            
            Heat_bottom[i] = Heat_bottom[i-1] + Heat_mixed_to_bottom;          //Exchange of Heat with mixed layer
                        
        
            
            
            //Calculate warming/////////////////////////////////////////////////////////////////////////////////
           
            //old warming no independent time-varying lambda's
           //DeltaT[i] = (1.0/lambda_time[2][i])*(1.0 - (epsilon*N_mixed[i])/(a_CO2*log(CO2[i]/CO2init) + epsilon_aero*R_aerosol[i] + R_WMGHG_nonCO2[i] + R_solar[i] + epsilon_volcanic*R_volcanic[i] ) ) * ( (a_CO2/Ib)* (I_em[i] + I_Usat_mixed[i] + I_Usat_deep[i] + I_Usat_upper[i] + I_Usat_inter[i] + I_Usat_bottom[i] + I_eq_heat[i]) + epsilon_aero*R_aerosol[i] + R_WMGHG_nonCO2[i] + R_solar[i] + R_volcanic[i] ) + T_noise[i];
            
            //New warming, independent time varying lambda's
            DeltaT[i] = (1.0/lambda_time[2][i]) * ((1.0 - N_mixed[i]/(a_CO2*log(CO2[i]/CO2init) + R_aerosol[i] + R_WMGHG_nonCO2[i] + R_solar[i] + R_volcanic[i] + R_noise[i])))*(a_CO2/Ib)*(I_em[i] + I_Usat_mixed[i] + I_Usat_deep[i] + I_Usat_upper[i] + I_Usat_inter[i] + I_Usat_bottom[i] + I_eq_heat[i])
                        + (1.0/lambda_time[0][i])*((1.0 - N_mixed[i]/(a_CO2*log(CO2[i]/CO2init) + R_aerosol[i] + R_WMGHG_nonCO2[i] + R_solar[i] + R_volcanic[i] + R_noise[i] )))*R_volcanic[i]
                        + (1.0/lambda_time[1][i])*((1.0 - N_mixed[i]/(a_CO2*log(CO2[i]/CO2init) + R_aerosol[i] + R_WMGHG_nonCO2[i] + R_solar[i] + R_volcanic[i] + R_noise[i] )))*R_solar[i]
                        + (1.0/lambda_time[3][i])*((1.0 - N_mixed[i]/(a_CO2*log(CO2[i]/CO2init) + R_aerosol[i] + R_WMGHG_nonCO2[i] + R_solar[i] + R_volcanic[i] + R_noise[i] )))*R_aerosol[i]
                        + (1.0/lambda_time[4][i])*((1.0 - N_mixed[i]/(a_CO2*log(CO2[i]/CO2init) + R_aerosol[i] + R_WMGHG_nonCO2[i] + R_solar[i] + R_volcanic[i] + R_noise[i] )))*R_WMGHG_nonCO2[i]
                        + (1.0/lambda_time[5][i])*((1.0 - N_mixed[i]/(a_CO2*log(CO2[i]/CO2init) + R_aerosol[i] + R_WMGHG_nonCO2[i] + R_solar[i] + R_volcanic[i] + R_noise[i] )))*R_noise[i]
                        + T_noise[i];
            
            
            if(DeltaT[i] - DeltaT[i-1] > 2.0)
            {
                DeltaT[i] = DeltaT[i-1] +2.0;
            }
            if(DeltaT[i] - DeltaT[i-1] < -2.0)
            {
                DeltaT[i] = DeltaT[i-1] -2.0;
            }
            
            
            //Calculate sea level rise//////////////////////////////////////////////////////////////////////////
            
            //Williams et al coefficient
            //Steric_rise[i] = Steric_rise[i-1] + (steric_coeff/1000.0)*(area/area_ocean)*(f_heat_ocean*N_mixed[i])*dt; //steric_coeff [mm yr-1 (W m-2)-1] is the sea level rise per year per unit ocean heat flux from Williams et al (2012) in GRL
            
            
            
            //Rahmstorf (2007) T related ice_melt sea level rise
            //ice_melt_rise[i] = ice_melt_rise[i-1] + ice_coeff*DeltaT[i]*dt; //Modelled after Rahmstorf (2007) and Vermeer and Rahmstorf (2009) with ice_coeff in m (K.yr)-1
            
            
            //Calculate surface ocean acidification/////////////////////////////////////////////////////////////
            
            Sea_surf_acid[i] = -0.387*log(CO2[i]/CO2init) - 2.41*(-(I_Usat_mixed[i]*1.0e15/12.0)/(vol_mixed*1.3e18));       //Using DpH = -0.387Dln (CO2) is at DIC saturation, then DpH = -2.41Delta DIC (DIC in moles C per m3) to account for distance from DIC saturation. Coefficients obtained from numerical carbonate chemistry solver for characteristic seasurface values.
            
            
            ///////////////////////////////////////////////////////////////////////
            
            
            
            
            
            
            
              
            //Data check at start of 2018////////////////////////////////////////////////////////////////////
            if ( year[i] >= 2020.0-dt/2.0 && year[i] <= 2020.0+dt/2.0  )
            {
                
                
                //Start with data_fit =0, add +2 to data_fit when within best estimate of observational constraint. Add +1 to data_fit when just outside best estimate range (up to 50% either way)
                //Accept simulation as observationally consistent based on value of data_fit once all observational consistency tests are conducted
                cost_function = 1.0; //This is multiplied by
                data_fit = 0;
                data_fit2 = 0;
                double random_test;
                
                //Set the error term outside the 90/95% ranges equal to zero initially
                error_term = 0.0;
                double obs_mid = 0.0;
                double obs_range = 0.0;
                
                double DT_1961_1990 = 0.0;
                for (int n= (1961-Start_year)*t_step_per_year; n<(1991-Start_year)*t_step_per_year; n++)
                {
                    DT_1961_1990 += DeltaT[n];
                }
                DT_1961_1990 = DT_1961_1990/double((1991-1961)*t_step_per_year);
                
                
                double DT_1850_1899 = 0.0;
                for (int n= (1850-Start_year)*t_step_per_year; n<(1900-Start_year)*t_step_per_year; n++)//Add up contributions
                {
                    DT_1850_1899 += DeltaT[n];
                    
                }
                DT_1850_1899  =  DT_1850_1899/double((1900-1850)*t_step_per_year);
                
                
                //2010 to 2019 inclusive
                double DT_2010_2019 = 0.0;
                
                for (int n= (2010-Start_year)*t_step_per_year; n<(2020-Start_year)*t_step_per_year; n++)//Add up contributions
                {
                    DT_2010_2019 += DeltaT[n];
                }
                DT_2010_2019  =  DT_2010_2019/double((2010-2000)*t_step_per_year);
                
                //2000 to 2009 inclusive
                double DT_2000_2019 = 0.0;
                
                for (int n= (2000-Start_year)*t_step_per_year; n<(2010-Start_year)*t_step_per_year; n++)//Add up contributions
                {
                    DT_2000_2019 += DeltaT[n];
                }
                DT_2000_2019  =  DT_2000_2019/double((2010-2000)*t_step_per_year);
                
                //1990 to 1999 inclusive
                double DT_1990_1999 = 0.0;
                
                for (int n= (1990-Start_year)*t_step_per_year; n<(2000-Start_year)*t_step_per_year; n++)//Add up contributions
                {
                    DT_1990_1999 += DeltaT[n];
                    
                }
                DT_1990_1999  =  DT_1990_1999/double((2000-1990)*t_step_per_year);
                
                //1980 to 1990 inclusive
                double DT_1980_1999 = 0.0;
                
                for (int n= (1980-Start_year)*t_step_per_year; n<(1990-Start_year)*t_step_per_year; n++)//Add up contributions
                {
                    DT_1980_1999 += DeltaT[n];
                    
                }
                DT_1980_1999  =  DT_1980_1999/double((1990-1980)*t_step_per_year);
                
                //1960 to 1979 inclusive
                double DT_1960_1979=0.0;
                for (int n= (1960-Start_year)*t_step_per_year; n<(1980-Start_year)*t_step_per_year; n++)
                {
                    DT_1960_1979 += DeltaT[n];
                }
                DT_1960_1979 = DT_1960_1979/double((1980-1960)*t_step_per_year);
                
                //1940 to 1959 inclusive
                double DT_1940_1959=0.0;
                for (int n= (1940-Start_year)*t_step_per_year; n<(1960-Start_year)*t_step_per_year; n++)
                {
                    DT_1940_1959 += DeltaT[n];
                }
                DT_1940_1959 = DT_1940_1959/double((1960-1940)*t_step_per_year);
                
                //1920 to 1939 inclusive
                double DT_1920_1939=0.0;
                for (int n= (1920-Start_year)*t_step_per_year; n<(1940-Start_year)*t_step_per_year; n++)
                {
                    DT_1920_1939 += DeltaT[n];
                }
                DT_1920_1939 = DT_1920_1939/double((1940-1920)*t_step_per_year);
                
                //1900 to 1919 inclusive
                double DT_1900_1919=0.0;
                for (int n= (1900-Start_year)*t_step_per_year; n<(1920-Start_year)*t_step_per_year; n++)
                {
                    DT_1900_1919 += DeltaT[n];
                }
                DT_1900_1919 = DT_1900_1919/double((1920-1900)*t_step_per_year);
                
                //Now rescale all to be change relative to the 1961 to 1990 period (for HadCRUT4 analysis)
                DT_2010_2019 = DT_2010_2019 - DT_1961_1990;
                DT_2000_2019 = DT_2000_2019 - DT_1961_1990;
                DT_1990_1999 = DT_1990_1999 - DT_1961_1990;
                DT_1980_1999 = DT_1980_1999 - DT_1961_1990;
                DT_1960_1979 = DT_1960_1979 - DT_1961_1990;
                DT_1940_1959 = DT_1940_1959 - DT_1961_1990;
                DT_1920_1939 = DT_1920_1939 - DT_1961_1990;
                DT_1900_1919 = DT_1900_1919 - DT_1961_1990;
                DT_1850_1899 = DT_1850_1899 - DT_1961_1990;
                
                double SST_1850_1899=0.0;//Average SS temperature 1880 to 1900
                for (int n= (1850-Start_year)*t_step_per_year; n<(1900-Start_year)*t_step_per_year; n++)
                {
                    SST_1850_1899 += Heat_mixed[n]/(vol_mixed*vol_total*c_p);
                }
                SST_1850_1899 = SST_1850_1899/double((1900-1850)*t_step_per_year);
                
                double SST_2000_2019=0.0;//Average SS temperature 1880 to 1900
                for (int n= (2000-Start_year)*t_step_per_year; n<(2020-Start_year)*t_step_per_year; n++)
                {
                    SST_2000_2019 += Heat_mixed[n]/(vol_mixed*vol_total*c_p);
                }
                SST_2000_2019 = SST_2000_2019/double((2020-2000)*t_step_per_year);
                
                
                double SST_1961_1990 = 0.0;//average SS temperature
                
                for (int n= (1961-Start_year)*t_step_per_year; n<(1990-Start_year)*t_step_per_year; n++)//Add up contributions
                {
                    SST_1961_1990 += Heat_mixed[n]/(vol_mixed*vol_total*c_p);
                    
                }
                SST_1961_1990  =  SST_1961_1990/double((1991-1961)*t_step_per_year);
                
                //Now re-calibreate to 1961-1990 average
                SST_1850_1899 = SST_1850_1899 - SST_1961_1990;
                SST_2000_2019 = SST_2000_2019 - SST_1961_1990;
                
                
                //Ocean heat content increase 700 to 2000m
                double Ocean_Heat_700_2000m_06_15 = 0.0;
                
                for (int n= (2006-Start_year)*t_step_per_year; n<(2016-Start_year)*t_step_per_year; n++)//Add up contributions
                {
                    Ocean_Heat_700_2000m_06_15 += Heat_inter[n] + (0.15/0.21)*Heat_deep[n];
                    
                }
                Ocean_Heat_700_2000m_06_15  =  Ocean_Heat_700_2000m_06_15/double((2016-2006)*t_step_per_year);
                
                double Ocean_Heat_700_2000m_1960_69= 0.0;
                for (int n= (1960-Start_year)*t_step_per_year; n<(1970-Start_year)*t_step_per_year; n++)//Add up contributions
                {
                    Ocean_Heat_700_2000m_1960_69 += Heat_inter[n] + (0.15/0.21)*Heat_deep[n];
                }
                Ocean_Heat_700_2000m_1960_69 =  Ocean_Heat_700_2000m_1960_69/double((1970-1960)*t_step_per_year);
                
                
                //Ocean heat content increase 700 to 2000m
                double Ocean_Heat_700m_06_15 = 0.0;
                
                for (int n= (2006-Start_year)*t_step_per_year; n<(2016-Start_year)*t_step_per_year; n++)//Add up contributions
                {
                    Ocean_Heat_700m_06_15 += Heat_mixed[n] + Heat_upper[n];
                }
                Ocean_Heat_700m_06_15  =  Ocean_Heat_700m_06_15/double((2016-2006)*t_step_per_year);
                
                
                double Ocean_Heat_700m_1960_69 = 0.0;
                for (int n= (1960-Start_year)*t_step_per_year; n<(1970-Start_year)*t_step_per_year; n++)//Add up contributions
                {
                    Ocean_Heat_700m_1960_69 += Heat_mixed[n] + Heat_upper[n];
                }
                Ocean_Heat_700m_1960_69 = Ocean_Heat_700m_1960_69/double((1970-1960)*t_step_per_year);
                
                Ocean_Heat_700m_1960_69 = Ocean_Heat_700m_06_15 - Ocean_Heat_700m_1960_69;
                Ocean_Heat_700_2000m_1960_69  = Ocean_Heat_700_2000m_06_15 - Ocean_Heat_700_2000m_1960_69;
                
                double OHC_1960_2016 =  (Heat_mixed[(2016-Start_year)*t_step_per_year-1] + Heat_upper[(2016-Start_year)*t_step_per_year-1] + Heat_inter[(2016-Start_year)*t_step_per_year-1] + Heat_deep[(2016-Start_year)*t_step_per_year-1] + Heat_bottom[(2016-Start_year)*t_step_per_year-1])
                -  (Heat_mixed[(1960-Start_year)*t_step_per_year] + Heat_upper[(1960-Start_year)*t_step_per_year] + Heat_inter[(1960-Start_year)*t_step_per_year] + Heat_deep[(1960-Start_year)*t_step_per_year] + Heat_bottom[(1960-Start_year)*t_step_per_year]);
                
                
                //The linear slope in ocean heat content from 1955 to 2017
                double OHC700m_55_17_slope = 0.0;
                double SumX = 0.0;
                double SumY= 0.0;
                double SumXY = 0.0;
                double SumX2 = 0.0;
                
                for(int n=(1955-Start_year)*t_step_per_year;n<(2018-Start_year)*t_step_per_year; n++)
                {
                    SumX += year[n];
                    SumX2 += year[n]*year[n];
                    SumY += Heat_mixed[n]+Heat_upper[n];
                    SumXY += (Heat_mixed[n]+Heat_upper[n])*year[n];
                }
                OHC700m_55_17_slope = (double((2018-1955)*t_step_per_year)*SumXY - (SumX*SumY))/(double((2018-1955)*t_step_per_year)*SumX2 - SumX*SumX);
                OHC700m_55_17_slope = OHC700m_55_17_slope/(area*24.0*60.0*60.0*365.0);
                
                double OHC2000m_55_17_slope = 0.0;
                SumX = 0.0;
                SumY= 0.0;
                SumXY = 0.0;
                SumX2 = 0.0;
                
                for(int n=(1955-Start_year)*t_step_per_year;n<(2018-Start_year)*t_step_per_year; n++)
                {
                    SumX += year[n];
                    SumX2 += year[n]*year[n];
                    SumY += Heat_mixed[n]+Heat_upper[n]+Heat_inter[n]+(0.15/0.21)*Heat_deep[n];
                    SumXY += (Heat_mixed[n]+Heat_upper[n]+Heat_inter[n]+(0.15/0.21)*Heat_deep[n])*year[n];
                }
                OHC2000m_55_17_slope = (double((2018-1955)*t_step_per_year)*SumXY - (SumX*SumY))/(double((2018-1955)*t_step_per_year)*SumX2 - SumX*SumX);
                OHC2000m_55_17_slope = OHC2000m_55_17_slope/(area*24.0*60.0*60.0*365.0);
                
                
                
                //Ocean carbon uptake 1982 to 2018:
                double Ocean_C_1982_2018 = 0.0;
                Ocean_C_1982_2018 = (I_em[(2018-Start_year)*t_step_per_year] -I_em[(1982-Start_year)*t_step_per_year])  - (CO2[(2018-Start_year)*t_step_per_year] - CO2[(1982-Start_year)*t_step_per_year]);
                
                
                
                
                //Now caclulate cost function for warming
                cost_function = 1.0;
                
                //Choose HadCRUT4 or Cowtan and Way///
                
                if (T_data == 0)
                {
                    //Global Mean Surface Temperature, HadCRUT4
                    
                     cost_function = cost_function * exp( -0.5*(((DT_2000_2019 - 0.537)*(DT_2000_2019 - 0.537))/(0.045*0.045)));
                    cost_function = cost_function * exp( -0.5*(((DT_1980_1999 - 0.185)*(DT_1980_1999 - 0.185))/(0.043*0.043)));
                    cost_function = cost_function * exp( -0.5*(((DT_1960_1979 - (-0.066))*(DT_1960_1979 - (-0.066)))/(0.042*0.042)));
                    cost_function = cost_function * exp( -0.5*(((DT_1940_1959 - (-0.034))*(DT_1940_1959 - (-0.034)))/(0.053*0.053)));
                    cost_function = cost_function * exp( -0.5*(((DT_1920_1939 - (-0.178))*(DT_1920_1939 - (-0.178)))/(0.057*0.057)));
                    cost_function = cost_function * exp( -0.5*(((DT_1900_1919 - (-0.388))*(DT_1900_1919 - (-0.388)))/(0.066*0.066)));
                    cost_function = cost_function * exp( -0.5*(((DT_1850_1899 - (-0.314))*(DT_1850_1899 - (-0.314)))/(0.084*0.084)));
                    
                }
                
                if (T_data == 1)
                {
                    //Global Mean Surface Temperature, Cowtan and Way Had4krig annual v2.0.0
                    cost_function = cost_function * exp( -0.5*(((DT_2000_2019 - 0.577)*(DT_2000_2019 - 0.577))/(0.03*0.03)));
                    cost_function = cost_function * exp( -0.5*(((DT_1980_1999 - 0.190)*(DT_1980_1999 - 0.190))/(0.021*0.021)));
                    cost_function = cost_function * exp( -0.5*(((DT_1960_1979 - (-0.061))*(DT_1960_1979 - (-0.061)))/(0.016*0.016)));
                    cost_function = cost_function * exp( -0.5*(((DT_1940_1959 - (-0.021))*(DT_1940_1959 - (-0.021)))/(0.028*0.028)));
                    cost_function = cost_function * exp( -0.5*(((DT_1920_1939 - (-0.168))*(DT_1920_1939 - (-0.168)))/(0.031*0.031)));
                    cost_function = cost_function * exp( -0.5*(((DT_1900_1919 - (-0.378))*(DT_1900_1919 - (-0.378)))/(0.043*0.043)));
                    cost_function = cost_function * exp( -0.5*(((DT_1850_1899 - (-0.351))*(DT_1850_1899 - (-0.351)))/(0.058*0.058)));
                    
                }
                
                
                
                //Global Mean SST, HadSST3
                cost_function = cost_function*exp( -0.5*(((SST_1850_1899 - (-0.276))*(SST_1850_1899 - (-0.276)))/(0.076*0.076)));
                
                //Choose Cheng et al. or NODC
                
                if (OHC_data == 0)
                {
                    //Cheng et al. ocean heat content
                    
                    //Ocean heat content from *start* of Jan 1960 t0 2006-2015 average
                    cost_function = cost_function*exp( -0.5*(((Ocean_Heat_700_2000m_1960_69 - 7.56e22)*(Ocean_Heat_700_2000m_1960_69 - 7.56e22))/(1.23e22*1.23e22)));
                    cost_function = cost_function*exp( -0.5*(((Ocean_Heat_700m_1960_69  - 17.78e22)*(Ocean_Heat_700m_1960_69  - 17.78e22))/(1.38e22*1.38e22)));
                    
                    //Whole ocean heat content 1960 to 2015 from Cheng et al
                    cost_function = cost_function*exp( -0.5*(((OHC_1960_2016- 360.0e21)*(OHC_1960_2016- 360.0e21))/(35.0e21*35.0e21)));
                    
                }
                
                if (OHC_data == 1)
                {
                    //This is NODC and IPCC (2013)/////////////////
                    cost_function = cost_function*exp( -0.5*(((Ocean_Heat_0_700m_2015_2019 - 2.119e23)*(Ocean_Heat_0_700m_2015_2019 - 2.119e23))/(1.85e22*1.85e22)));
                    cost_function = cost_function*exp( -0.5*(((Ocean_Heat_0_2000m_2015_2019  - 3.145e23)*(Ocean_Heat_0_2000m_2015_2019  - 3.145e23))/(2.03e22*2.03e22)));
                    
                    //Total Earth system heat content 1971 to 2010 from IPCC (2013) AR5
                    cost_function = cost_function*exp( -0.5*(((Total_ES_HC_1971_2010- 274.0e21)*(Total_ES_HC_1971_2010 - 274.0e21))/(48.0e21*48.0e21)));
                    
                    ////////////////////////////////////////////////
                    
                }
                
                //rate of change in OHC in top 700m and top 2000m from Zanna et al. 2019
                //cost_function = cost_function*exp( -0.5*(((OHC700m_55_17_slope- 0.22)*(OHC700m_55_17_slope - 0.22))/(0.05*0.05)));
                //cost_function = cost_function*exp( -0.5*(((OHC2000m_55_17_slope  - 0.30)*(OHC2000m_55_17_slope  - 0.30))/(0.06*0.06)));
                
                //Whole ocean heat content 1960 to 2017 from Zanna et al. 2019
                //cost_function = cost_function*exp( -0.5*(((OHC_1960_2016- 323.0e21)*(OHC_1960_2016- 323.0e21))/(66.0e21*66.0e21)));
                
                //Ocean carbon uptake
                cost_function = cost_function*exp( -0.5*(((Ocean_C_1982_2018  - 71.2)*(Ocean_C_1982_2018  - 71.2))/(24.3*24.3)));
                
                //cout << cost_function << '\t' << cost_function_max << endl;
                if(cost_function > cost_function_max)
                {
                    cost_function_max = cost_function;
                    cout << j << '\t' << cost_function_max << '\t' << exp( -0.5*(((Ocean_Heat_700_2000m_1960_69 - 7.56e22)*(Ocean_Heat_700_2000m_1960_69 - 7.56e22))/(1.23e22*1.23e22))) << '\t' << exp( -0.5*(((Ocean_Heat_700m_1960_69  - 17.78e22)*(Ocean_Heat_700m_1960_69  - 17.78e22))/(1.38e22*1.38e22))) << '\t' << exp( -0.5*(((OHC_1960_2016- 360.0e21)*(OHC_1960_2016- 360.0e21))/(35.0e21*35.0e21))) << endl;
                }
                random_test = getRandomLinear(0.0,0.001);
                //random_test = 0.01;
                
                if(cost_function > random_test)
                {
                    data_fit = 12;
                }
                
                
                
                
                
                
                
                //Cumulative carbon emissions 1750 to 2011: data is 555 [470 to 640] PgC
                double Emissions_1750_2011 = 0.0;
                Emissions_1750_2011 =  I_em[i-7*t_step_per_year] + (I_ter[i-7*t_step_per_year] - I_ter[0]);
                
                
                //OLD constraint:
                //Terrestrial carbon uptake 1750 to 2011: data is 160 [70 to 250] PgC
                
                //New constraint:
                //Terrestrial carbon budget in GCB 2018's 16 DGVMs have mean and +/- 2 sigma of 213 [96 to 331] PgC from 1750 to start of 2017
                double Terrestrial_C_1959_2018 = 0.0;
                Terrestrial_C_1959_2018 = (I_ter[i] - I_ter[i-59*t_step_per_year]);
                /*
                if(Terrestrial_C_1959_2018 >= 61.8 && Terrestrial_C_1959_2018 <= 192.3)
                {
                    data_fit += 2;
                    
                }
                obs_mid = 0.5*(61.8+192.3);
                obs_range = 192.3-61.8;
                if(Terrestrial_C_1959_2018 < 61.8 )
                {
                    data_fit2 += 2;
                    error_term += ( abs(obs_mid - Terrestrial_C_1959_2018)  / (0.5*obs_range) ) - 1.0;
                }
                if(Terrestrial_C_1959_2018 > 192.3 )
                {
                    data_fit2 += 2;
                    error_term += ( abs(obs_mid - Terrestrial_C_1959_2018)  / (0.5*obs_range) ) - 1.0;
                }
               */
                
                //Terrestrial carbon uptake
                double Terrestrial_C_2000_2009 = 0.0;
                Terrestrial_C_2000_2009 = (I_ter[i-8*t_step_per_year] - I_ter[i-18*t_step_per_year]);
                /*
                if(Terrestrial_C_2000_2009 >= 10.0*1.4 && Terrestrial_C_2000_2009 <= 10.0*3.8)
                {
                    data_fit += 2;
                }
                obs_mid = 0.5*(10.0*1.4+10.0*3.8);
                obs_range = (10.0*3.8)-(10.0*1.4);
                if( Terrestrial_C_2000_2009 < 1.4*10.0)
                {
                    data_fit2 +=2;
                    error_term += ( abs(obs_mid - Terrestrial_C_2000_2009)  / (0.5*obs_range) ) - 1.0;
                }
                if(Terrestrial_C_2000_2009 > 3.8*10.0)
                {
                    data_fit2 +=2;
                    error_term += ( abs(obs_mid - Terrestrial_C_2000_2009)  / (0.5*obs_range) ) - 1.0;
                }
                */
            
                
                
                
                
                
               //If use sea level rise component, add in additional sea level observational constraints here
                
                
                
                
                
                
                //Global mean SST increase 1870-1900 average to 1985-2005 average estimated from fig. 5 in Jha et al 2014 (10 CMIP5 models range from ~0.27 increase to ~0.66 increase, say 0.2 to 0.7 is consistent)
                double SST_ave_1870_1900=0.0;
                double SST_ave_1885_1915=0.0;
                double SST_ave_1976_2005=0.0;
                double SST_increase_1885_1915_to_1976_2005=0.0;
                
                for (int n=0; n<30; n++)
                {
                    SST_ave_1870_1900+= (Heat_mixed[i-(147-n)*t_step_per_year]/(vol_mixed*vol_total*c_p) );
                    SST_ave_1976_2005+= (Heat_mixed[i-(41-n)*t_step_per_year]/(vol_mixed*vol_total*c_p) ) ;
                    SST_ave_1885_1915+= (Heat_mixed[i-(117-n)*t_step_per_year]/(vol_mixed*vol_total*c_p) ) ;
                    if(n==29)//divide through at the end to get the average anthropogenic SST contribution for the 30-year periods
                    {
                        SST_ave_1870_1900 = SST_ave_1870_1900/30.0;
                        SST_ave_1976_2005 = SST_ave_1976_2005/30.0;
                        SST_ave_1885_1915 = SST_ave_1885_1915/30.0;
                    }
                }
                
                
    
                
                //make it pass data test if on idealised (non-historic) scenario
                if(Future_scenario >=15 || Use_stored_params == 1)
                {
                    data_fit = 16;
                }
                
                //Only accept the ensemble member into the data-consistent ensemble if data_fit is 17 or more
                //datafit <= m, where m should be tuned to '2*n - 1', where n is the number of observational constraints applied - see Goodwin (2016).
                if(data_fit >= 10)
                {
                    //prints to screen when observationally consistent ensemble member is found, and what number simulation is observationally consistent
                    cout << "Simulation" << '\t' <<  j << '\t' << "is observation consistent" << '\t' << cost_function  << '\t' << cost_function_max << endl;
                    
                    //Writes the historic simulated climate properties out to file
                    results << j << '\t' << cost_function << '\t' << DT_1850_1899 << '\t' <<  DT_1900_1919 << '\t' <<  DT_1920_1939 << '\t' << DT_1940_1959 << '\t' << DT_1960_1979 <<'\t' << DT_1980_1999 << '\t' << DT_2000_2019 << '\t' << SST_1850_1899 << '\t' << SST_2000_2019 << '\t' << Ocean_Heat_700m_1960_69 << '\t'  << Ocean_Heat_700_2000m_1960_69 << '\t' << Ocean_C_1982_2018 << '\t' << OHC_1960_2016<< '\t' << Terrestrial_C_1959_2018 << '\t' << Terrestrial_C_2000_2009/10.0 << endl;
                    
                    //Writes model input parameters to file/////////////////////////////////////////
                    inputs << j << '\t' << lambda << '\t' <<  lambda_Planck << '\t' <<  lambda_WVLR << '\t' << lambda_Cloud_Fast << '\t' << lambda_albedo << '\t' <<  lambda_Cloud_SST << '\t' << lambda_1000yr << '\t' << tau_WVLR << '\t' << tau_Cloud_Fast << '\t' << tau_albedo << '\t' <<  tau_Cloud_SST << '\t' << ratio << '\t' << ratio2 << '\t' << tau_C_mixed << '\t' << tau_C_upper << '\t' << tau_C_inter <<'\t' << tau_C_deep << '\t' << tau_C_bottom << '\t' << Ib << '\t' << a_CO2 << '\t' << R_aerosol_Uncert << '\t' << R_WMGHG_nonCO2_Uncert << '\t' << R_volcanic_coeff << '\t' << (a_CO2*log(CO2[i-6*t_step_per_year]/CO2init) + R_aerosol[i-6*t_step_per_year] + R_WMGHG_nonCO2[i-6*t_step_per_year] ) << '\t' << a_CO2*log(CO2[i-6*t_step_per_year]/CO2init) << '\t' << R_WMGHG_nonCO2[i-6*t_step_per_year] << '\t' << R_aerosol[i-6*t_step_per_year] << '\t' << f_heat_ocean << '\t' << dNPPdT << '\t' << dtaudT << '\t' << gamma_K << '\t' << ratio_volcanic_CO2 << '\t' << gamma_aero_SOx*Aero_SOx_em[(2010-Start_year)] << '\t' << gamma_aero_NOx*Aero_NOx_em[(2010-Start_year)] << '\t' << gamma_aero_OC*Aero_OC_em[(2010-Start_year)] << '\t' << gamma_aero_BC*Aero_BC_em[(2010-Start_year)] << '\t' << gamma_aero_NH3*Aero_NH3_em[(2010-Start_year)] << '\t' << gamma_aero_ENMVOC*Aero_ENMVOC_em[(2010-Start_year)] << '\t' << AeroRF_ind_2011  << '\t' << Uncert_CH4 << '\t' << Uncert_N2O << '\t' << Uncert_Halocarbons << endl;
                    
                    //For scenarios with historic constraints, store observation-consistent parameters
                    if(Future_scenario <15 && Use_stored_params == 0)
                    {
                        //Stores observation-consistent parameter values
                        a_CO2_obs.push_back(a_CO2);
                        dNPPdT_obs.push_back(dNPPdT);
                        gamma_K_obs.push_back(gamma_K);
                        dtaudT_obs.push_back(dtaudT);
                        ratio_obs.push_back(ratio);
                        ratio2_obs.push_back(ratio2);
                        lambda_Planck_obs.push_back(lambda_Planck);
                        lambda_WVLR_obs.push_back(lambda_WVLR);
                        lambda_Cloud_Fast_obs.push_back(lambda_Cloud_Fast);
                        lambda_Cloud_SST_obs.push_back(lambda_Cloud_SST);
                        lambda_albedo_obs.push_back(lambda_albedo);
                        f_heat_ocean_obs.push_back(f_heat_ocean);
                        tau_C_mixed_obs.push_back(tau_C_mixed);
                        tau_C_upper_obs.push_back(tau_C_upper);
                        tau_C_inter_obs.push_back(tau_C_inter);
                        tau_C_deep_obs.push_back(tau_C_deep);
                        tau_C_bottom_obs.push_back(tau_C_bottom);
                        Ib_obs.push_back(Ib);
                        tau_WVLR_obs.push_back(tau_WVLR);
                        tau_Cloud_Fast_obs.push_back(tau_Cloud_Fast);
                        tau_albedo_obs.push_back(tau_albedo);
                        tau_Cloud_SST_obs.push_back(tau_Cloud_SST);
                        R_aerosol_Uncert_obs.push_back(R_aerosol_Uncert);
                        R_WMGHG_nonCO2_Uncert_obs.push_back(R_WMGHG_nonCO2_Uncert);
                        R_volcanic_coeff_obs.push_back(R_volcanic_coeff);
                        gamma_aero_SOx_obs.push_back(gamma_aero_SOx);
                        gamma_aero_BC_obs.push_back(gamma_aero_BC);
                        gamma_aero_OC_obs.push_back(gamma_aero_OC);
                        gamma_aero_ENMVOC_obs.push_back(gamma_aero_ENMVOC);
                        gamma_aero_NOx_obs.push_back(gamma_aero_NOx);
                        gamma_aero_NH3_obs.push_back(gamma_aero_NH3);
                        AeroRF_ind_2011_obs.push_back(AeroRF_ind_2011);
                        Uncert_CH4_obs.push_back(Uncert_CH4);
                        Uncert_N2O_obs.push_back(Uncert_N2O);
                        Uncert_Halocarbons_obs.push_back(Uncert_Halocarbons);
                    }
                        
                    
                    
                        count_obs++;
                    
                    }
                
                
                
                
                
                //If not data-consistent then end simulation and move on to the next - no need to run past 2017 if not data-consistent
               if(data_fit <10)
               {
                   
                   i = tmax - 1;
                   
                   
               }
                     
               
                                         
            
                                     
            }
            
                
            
            
            //Output at end of year 2300/start of 2301.
            if (/*output == t_step_per_year-1 &&*/ year[i] >= 2301.0-dt/2.0 && year[i] <= 2301+dt/2.0 )
            {
                
                
                //Calculate 20-year averages for warming and sea level 
                double DT_2000_2020=0.0;
                double DT_46_65=0.0;
                double DT_81_00=0.0;
                double DS_86_05=0.0;
                double DS_46_65=0.0;
                double DS_81_00=0.0;
                double DS_S_86_05=0.0;
                double DS_S_46_65=0.0;
                double DS_S_81_00=0.0;
                
                double DTaverage_1850_1900 = 0.0;
                
                double DTaverage_2280_2300 =0.0;
                
                
                
                double DT_ave_2003_2012=0.0;
                
                for (int n=0; n<51*t_step_per_year; n++)
                {
                    DTaverage_1850_1900 += DeltaT[i-(451*t_step_per_year-n)]/(51.0*double(t_step_per_year));
                    
                }
                
                for (int n=0; n<20*t_step_per_year; n++)
                {
                    DTaverage_2280_2300 += DeltaT[i-(20*t_step_per_year-n)]/(20.0*double(t_step_per_year));
                    
                }
                
                
                
                for (int n=0; n<20*t_step_per_year; n++)//Add up contributions
                {
                    DT_2000_2020 += DeltaT[i-(301*t_step_per_year - n)]/(20.0*double(t_step_per_year));
                    DS_86_05 += (ice_melt_rise[i-(315*t_step_per_year-n)] + Steric_rise[i-(315*t_step_per_year-n)])/(20.0*double(t_step_per_year));
                    DS_S_86_05 += Steric_rise[i-(315*t_step_per_year-n)]/(20.0*double(t_step_per_year));
                    
                    DT_46_65 += DeltaT[i-(256*t_step_per_year-n)]/(20.0*double(t_step_per_year));
                    DS_46_65 += (ice_melt_rise[i-(256*t_step_per_year-n)] + Steric_rise[i-(256*t_step_per_year-n)])/(20.0*double(t_step_per_year));
                    DS_S_46_65 += Steric_rise[i-(256*t_step_per_year-n)]/(20.0*double(t_step_per_year));

                    DT_81_00 += DeltaT[i-(220*t_step_per_year-n)]/(20.0*double(t_step_per_year));
                    DS_81_00 += (ice_melt_rise[i-(220*t_step_per_year-n)] + Steric_rise[i-(220*t_step_per_year-n)])/(20.0*double(t_step_per_year));
                    DS_S_81_00 += Steric_rise[i-(220*t_step_per_year-n)]/(20.0*double(t_step_per_year));
                    

                    
                }
                
                for(int n = 0; n< 10*t_step_per_year; n++)
                {
                    DT_ave_2003_2012 += DeltaT[i-298*t_step_per_year+n];
                }
                DT_ave_2003_2012 = DT_ave_2003_2012/(10.0*double(t_step_per_year));
                
                double count_em = 0;
                if(data_fit>=10)
                {
                    //writes: ensemble-member number, mean temperature anonaly relative to preindustrial at 2081-2100, change in temperature anomaly 1986-2005 to 2081-2100, change in temperature anomaly 1986-2005 to 2046-2065, cumulative emissions 2017 to 2100, cumulative emissions 2012 to 2100
                    //projection << j << '\t' <<  DT_81_00 - DTaverage_1850_1900  << '\t' <<  DTaverage_2280_2300 - DTaverage_1850_1900  << '\t' << DT_81_00 - DT_2000_2020 << '\t' << DT_46_65 - DT_2000_2020 <<  '\t' << I_em[i] + I_ter[i] - (I_em[i-281*t_step_per_year] + I_ter[i-281*t_step_per_year]) << '\t' << I_em[i-200*t_step_per_year] + I_ter[i-200*t_step_per_year] - (I_em[i-281*t_step_per_year] + I_ter[i-281*t_step_per_year]) << '\t' << I_em[i] + I_ter[i] - (I_em[i-291*t_step_per_year] + I_ter[i-291*t_step_per_year]) << '\t' << I_em[i-200*t_step_per_year] + I_ter[i-200*t_step_per_year] - (I_em[i-291*t_step_per_year] + I_ter[i-291*t_step_per_year]) << '\t' << 1000.0*(DT_81_00)/(I_em[i-210*t_step_per_year]+I_ter[i-210*t_step_per_year]-I_ter[0]) << '\t' << 1000.0*(DTaverage_2280_2300)/(I_em[i-10*t_step_per_year]+I_ter[i-10*t_step_per_year]-I_ter[0]) << endl;
                    
                    int count_month = 0;
                    double T_ave = 0.0;
                    double N_ave = 0.0;
                    for (int i2 = 0; i2<=535*t_step_per_year; i2++)
                    {
                        
                        
                        if(count_month == 0)
                        {
                            T_ave = 0.0;
                            N_ave = 0.0;
                        }
                        count_month++;
                        T_ave = T_ave + DeltaT[i2]/double(t_step_per_month);
                        N_ave = N_ave + N_mixed[i2]/double(t_step_per_month);
                        /*
                        if(year[i2] >1969.999 && year[i2] < 2016.999 && count_month == t_step_per_month)
                        {
                            annual << year[i2] << '\t' << T_ave << endl;
                            
                        }
                        
                        
                        if(year[i2] > 1960-dt/2.0 && year[i2] < 2000.0+dt/2.0 && count_month == t_step_per_month)
                        {
                            
                            
                            Pinatubo_T << T_ave << '\t';
                            Pinatubo_N << N_ave << '\t' ;
                            Pinatubo_lambda_CO2 << lambda_time[2][i2] << '\t' ;
                            Pinatubo_lambda_volcanic << lambda_time[0][i2] << '\t' ;
                            Pinatubo_R_CO2 << (a_CO2*log(CO2[i2]/CO2init) ) << '\t' ;
                            Pinatubo_R << (a_CO2*log(CO2[i2]/CO2init) + R_aerosol[i2] + R_WMGHG_nonCO2[i2] + R_volcanic[i2] + R_solar[i2]) << '\t' ;
                            Pinatubo_R_volcanic << R_volcanic[i2] << '\t';
                            
                        }
                        */
                        if (count_month == t_step_per_month)
                        {
                            count_month = 0;
                        }
                       
                        
                        
                    }
                    /*Pinatubo_T << endl;
                    Pinatubo_N << endl;
                    Pinatubo_lambda_CO2 << endl;
                    Pinatubo_lambda_volcanic << endl;
                    Pinatubo_R_CO2 << endl;
                    Pinatubo_R << endl;
                    Pinatubo_R_volcanic << endl;*/
                    
                    
                    int years_back = 535;
                    if(Future_scenario >= 15)
                    {
                        //cout << year[i] << '\t' << Start_year << endl;
                        years_back = year[i] - Start_year - 1;
                    }
                    //Output stored properties by year
                        for (int i2 = -years_back; i2 <=0; i2 ++)
                        {
                            double T_ave=0.0;
                            double SAT_ave = 0.0;
                            double SST_ave=0.0;
                            double N_ave=0.0;
                            double Rtot_ave =0.0;
                            double RCO2_ave =0.0;
                            double Raerosol_ave =0.0;
                            double Rvolcanic_ave =0.0;
                            double Rsolar_ave =0.0;
                            double RWMGHG_nonCO2_ave = 0.0;
                            double NPP_ave = 0.0;
                            double year_ave = 0.0;
                            
                            for(int n=0; n<t_step_per_year; n++)
                            {
                                year_ave += year[i+i2*t_step_per_year - n]/double(t_step_per_year);
                                T_ave += DeltaT[i+i2*t_step_per_year - n]/double(t_step_per_year);
                            
                                SST_ave += (Heat_mixed[i+i2*t_step_per_year - n]/(vol_mixed*vol_total*c_p))/double(t_step_per_year);
                                N_ave += N_mixed[i+i2*t_step_per_year - n]/double(t_step_per_year);
                                Rtot_ave += (a_CO2*log(CO2[i+i2*t_step_per_year - n]/CO2init) + R_aerosol[i+i2*t_step_per_year - n] + R_WMGHG_nonCO2[i+i2*t_step_per_year - n] + R_volcanic[i+i2*t_step_per_year - n] + R_solar[i+i2*t_step_per_year - n])/double(t_step_per_year);
                                
                                RCO2_ave += (a_CO2*log(CO2[i+i2*t_step_per_year - n]/CO2init) )/double(t_step_per_year);
                                
                                Raerosol_ave += ( R_aerosol[i+i2*t_step_per_year - n] )/double(t_step_per_year);
                                
                                Rvolcanic_ave += ( R_volcanic[i+i2*t_step_per_year - n] )/double(t_step_per_year);
                                Rsolar_ave += ( R_solar[i+i2*t_step_per_year - n] )/double(t_step_per_year);
                                RWMGHG_nonCO2_ave += ( R_WMGHG_nonCO2[i+i2*t_step_per_year - n] )/double(t_step_per_year);
                                
                                NPP_ave += (getNPP(CO2[i-1+i2*t_step_per_year - n]*PgCtoppm, CO2init*PgCtoppm, DeltaT[i-1+i2*t_step_per_year - n], I_veg[i-1+i2*t_step_per_year - n], I_veg_init, NPP_init, gamma_K, dNPPdT) ) / double(t_step_per_year);
                                
                                
                            }
                            
                            SAT_ave = T_ave/0.95;
                            
                            annual_N_tot << N_ave << '\t' ;
                            
                            annual_warming << SAT_ave << '\t' ;
                            annual_warming_blended << T_ave << '\t' ;
                            annual_warming_SST << SST_ave << '\t';
                            
                            annual_lambda << (Rtot_ave - N_ave)/(T_ave) << '\t' ;
                            
                            annual_S_effective << (T_ave*a_CO2*log(2) ) / (Rtot_ave - N_ave) << '\t' ;
                            
                            annual_RCO2 << RCO2_ave << '\t' ;
                            annual_R_tot << Rtot_ave << '\t';
                            annual_R_anth << Rtot_ave - Rvolcanic_ave - Rsolar_ave << '\t';
                            annual_R_aerosol << Raerosol_ave << '\t' ;
                            annual_R_WMGHG_nonCO2 << RWMGHG_nonCO2_ave << '\t' ;
                            annual_R_volcanic << Rvolcanic_ave << '\t' ;
                            annual_R_solar << Rsolar_ave << '\t';
                            
                            annual_heatOcean << (Heat_mixed[i+i2*t_step_per_year] + Heat_upper[i+i2*t_step_per_year] + Heat_inter[i+i2*t_step_per_year] + Heat_deep[i+i2*t_step_per_year] + Heat_bottom[i+i2*t_step_per_year]) << '\t';
                            
                            annual_heat700m << (Heat_mixed[i+i2*t_step_per_year] + Heat_upper[i+i2*t_step_per_year] ) << '\t';
                            
                            annual_heatEarth <<   Heat_inter[i+i2*t_step_per_year] + (0.15/0.21)*Heat_deep[i+i2*t_step_per_year] << '\t' ;
                            
                            annual_landC << I_ter[i+i2*t_step_per_year] - I_ter[i+(i2-1)*t_step_per_year] << '\t' ;
                            annual_oceanC << I_em[i+i2*t_step_per_year] - (CO2[i+i2*t_step_per_year] - CO2[0]) - (I_em[i+(i2-1)*t_step_per_year] - (CO2[i+(i2-1)*t_step_per_year] - CO2[0]) ) << '\t' ;
                            //annual_Em << (I_em[i+i2*t_step_per_year] + I_ter[i+i2*t_step_per_year]) - I_ter[0]  << '\t' ;
                            annual_I_eq_Heat << I_eq_heat[i+i2*t_step_per_year] << '\t' ;
                            
                            //Instantaneous_TCRE << SAT_ave/(I_em[i+i2*t_step_per_year] + I_ter[i+i2*t_step_per_year] - I_ter[0]) << '\t';
                            
                            atmosCO2 << CO2[i+i2*t_step_per_year] << '\t' ;
                            //Airborne_fraction << (CO2[i+i2*t_step_per_year] - CO2[0]) / (I_em[i+i2*t_step_per_year] + I_ter[i+i2*t_step_per_year] - I_ter[0]) << '\t' ;
                            
                            //landC_plant << I_veg[i+i2*t_step_per_year] << '\t' ;
                            //landC_soil  << I_soil[i+i2*t_step_per_year] << '\t' ;
                            
                            //NPP_flux << NPP_ave << '\t' ;
                            
                            
                            
                            if(i2 == 0)
                            {
                                annual_warming << endl ;
                                annual_warming_blended << endl ;
                                annual_warming_SST << endl ;
                                
                                annual_lambda << endl ;
                                
                                annual_S_effective << endl ;
                                
                                annual_RCO2 << endl ;
                                annual_R_tot << endl ;
                                annual_R_anth << endl ;
                                annual_R_aerosol << endl ;
                                annual_R_WMGHG_nonCO2 << endl ;
                                annual_R_volcanic << endl ;
                                annual_R_solar << endl ;
                                
                                annual_N_tot << endl;
                                
                                annual_heatOcean << endl ;
                                
                                annual_heat700m << endl ;
                                
                                annual_heatEarth << endl ;
                                
                                //annual_landC << endl ;
                                annual_oceanC << endl ;
                                //annual_Em << endl ;
                                
                                //Instantaneous_TCRE << endl ;
                                
                                atmosCO2 << endl ;
                                //Airborne_fraction << endl ;
                                
                                //landC_plant << endl ;
                                //landC_soil  << endl ;
                                
                                //NPP_flux << endl ;
                                
                            }
                        }
                    //}
                    count_annual ++;
                    
                    
                    
                    count_em = 0;
                    int count_budget = 0;
                    for(int i2 = -450*t_step_per_year; i2<=0; i2++)
                    {
                        //Calculate emissions as of the START of 2020 (i.e. emissions as of the last time step in 2016).
                        double Em_2020 = I_em[i-281*t_step_per_year] + I_ter[i-281*t_step_per_year] - I_ter[0];
                        double T_ave = 0.0;
                        double DT_1850_1900_analysis=0.0;
                        
                        //Now do each 5 PgC intervals
                        if( (I_em[i+i2] + I_ter[i+i2] - I_ter[0]) - Em_2020 > 5.0*double(count_em) - 400.0 )
                        {
                            
                            
                            for (int n=0; n<51*t_step_per_year; n++)//Add up contributions
                            {
                                DT_1850_1900_analysis += DeltaT[i-451*t_step_per_year+n];
                                
                            }
                            DT_1850_1900_analysis  =  DT_1850_1900_analysis/(double(51*t_step_per_year));
                            
                            //Calculates year-average temperatures centred on January 1st each year
                            //This is so when plot year-average temperatures, plotted points are on Jan 1st.
                            for(int n=0; n<t_step_per_year; n++)
                            {
                                T_ave += DeltaT[i+i2 - n + 6 ]/double(t_step_per_year);
                            }
                            //em_warming_T_2020 << T_ave - DT_1850_1900_analysis << '\t';
                            //em_warming_I_2017 << I_em[i+i2] + I_ter[i+i2] - I_ter[0] - Em_2020 << '\t';
                            count_em ++;
                        }
                        if(i2 == 0)
                        {
                            //em_warming_T_2020 << endl;
                            //em_warming_I_2017 << endl;
                        }
                        
                        
                        
                        
                        
                    
                    
                        
                        
                        
                        
                    }
                    
                    
                    
                    
                    
                    
                }
                
                
                
                
                //use if don;t want all the other outputs
               i = tmax - 1;
                
                
                
            }
            
            
            
            output++;
            
            
            //Set output to zero at start of each year for first 1000 years of simulation
            if (i<=t_step_per_year*1000 && output == t_step_per_year)
            {
                output = 0;
            }
            
            
            
        }
    
        
    }
    
    //If there have been any observation-consistent simulations then output the results
    if(a_CO2_obs.size()>=1)
    {
        
        //Write out stored inputs
        stored_inputs << "a_CO2_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << a_CO2_obs[i] << "," ;
        }
        stored_inputs << a_CO2_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "dNPPdT_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << dNPPdT_obs[i] << "," ;
        }
        stored_inputs << dNPPdT_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "gamma_K_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << gamma_K_obs[i] << "," ;
        }
        stored_inputs << gamma_K_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "dtaudT_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << dtaudT_obs[i] << "," ;
        }
        stored_inputs << dtaudT_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "ratio_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << ratio_obs[i] << "," ;
        }
        stored_inputs << ratio_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "ratio2_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << ratio2_obs[i] << "," ;
        }
        stored_inputs << ratio2_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "lambda_Planck_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << lambda_Planck_obs[i] << "," ;
        }
        stored_inputs << lambda_Planck_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "lambda_WVLR_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << lambda_WVLR_obs[i] << "," ;
        }
        stored_inputs << lambda_WVLR_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "lambda_Cloud_Fast_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << lambda_Cloud_Fast_obs[i] << "," ;
        }
        stored_inputs << lambda_Cloud_Fast_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "lambda_Cloud_SST_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << lambda_Cloud_SST_obs[i] << "," ;
        }
        stored_inputs << lambda_Cloud_SST_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "lambda_albedo_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << lambda_albedo_obs[i] << "," ;
        }
        stored_inputs << lambda_albedo_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "f_heat_ocean_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << f_heat_ocean_obs[i] << "," ;
        }
        stored_inputs << f_heat_ocean_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "tau_C_mixed_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << tau_C_mixed_obs[i] << "," ;
        }
        stored_inputs << tau_C_mixed_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "tau_C_upper_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << tau_C_upper_obs[i] << "," ;
        }
        stored_inputs << tau_C_upper_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "tau_C_inter_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << tau_C_inter_obs[i] << "," ;
        }
        stored_inputs << tau_C_inter_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "tau_C_deep_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << tau_C_deep_obs[i] << "," ;
        }
        stored_inputs << tau_C_deep_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "tau_C_bottom_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << tau_C_bottom_obs[i] << "," ;
        }
        stored_inputs << tau_C_bottom_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "Ib_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << Ib_obs[i] << "," ;
        }
        stored_inputs << Ib_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "tau_WVLR_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << tau_WVLR_obs[i] << "," ;
        }
        stored_inputs << tau_WVLR_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "tau_Cloud_Fast_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << tau_Cloud_Fast_obs[i] << "," ;
        }
        stored_inputs << tau_Cloud_Fast_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "tau_albedo_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << tau_albedo_obs[i] << "," ;
        }
        stored_inputs << tau_albedo_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "tau_Cloud_SST_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << tau_Cloud_SST_obs[i] << "," ;
        }
        stored_inputs << tau_Cloud_SST_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "R_aerosol_Uncert_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << R_aerosol_Uncert_obs[i] << "," ;
        }
        stored_inputs << R_aerosol_Uncert_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "R_WMGHG_nonCO2_Uncert_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << R_WMGHG_nonCO2_Uncert_obs[i] << "," ;
        }
        stored_inputs << R_WMGHG_nonCO2_Uncert_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "R_volcanic_coeff_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << R_volcanic_coeff_obs[i] << "," ;
        }
        stored_inputs << R_volcanic_coeff_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "gamma_aero_SOx_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << gamma_aero_SOx_obs[i] << "," ;
        }
        stored_inputs << gamma_aero_SOx_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "gamma_aero_BC_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << gamma_aero_BC_obs[i] << "," ;
        }
        stored_inputs << gamma_aero_BC_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "gamma_aero_OC_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << gamma_aero_OC_obs[i] << "," ;
        }
        stored_inputs << gamma_aero_OC_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "gamma_aero_ENMVOC_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << gamma_aero_ENMVOC_obs[i] << "," ;
        }
        stored_inputs << gamma_aero_ENMVOC_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "gamma_aero_NOx_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << gamma_aero_NOx_obs[i] << "," ;
        }
        stored_inputs << gamma_aero_NOx_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "gamma_aero_NH3_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << gamma_aero_NH3_obs[i] << "," ;
        }
        stored_inputs << gamma_aero_NH3_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "AeroRF_ind_2011_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << AeroRF_ind_2011_obs[i] << "," ;
        }
        stored_inputs << AeroRF_ind_2011_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "Uncert_N2O_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << Uncert_N2O_obs[i] << "," ;
        }
        stored_inputs << Uncert_N2O_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "Uncert_CH4_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << Uncert_CH4_obs[i] << "," ;
        }
        stored_inputs << Uncert_CH4_obs[count_obs-1] << "};" << endl;
        
        stored_inputs << "Uncert_Halocarbons_obs={";
        for (int i = 0; i<count_obs-1; i++)
        {
            stored_inputs << Uncert_Halocarbons_obs[i] << "," ;
        }
        stored_inputs << Uncert_Halocarbons_obs[count_obs-1] << "};" << endl;
        
        
        stored_inputs << "return;" << endl << "}" << endl;
        
    }
        
   

    
    
    
    //Close open files
    
    results.close();
    inputs.close();
    projection.close();
    
    
    annual_warming.close();
    annual_lambda.close();
    
    annual_S_effective.close();
    
    annual_R_solar.close();
    annual_RCO2.close();
    annual_R_tot.close();
    annual_R_aerosol.close();
    annual_R_WMGHG_nonCO2.close();
    annual_R_volcanic.close();
    //annual_SST.close();
    //annual_heat.close();
    annual_heat700m.close();
    annual_oceanC.close();
    annual_landC.close();
    annual_Em.close();
    //em_warming_T_2020.close();
    //em_warming_I_2017.close();
    stored_inputs.close();
    

    
    
    
    
    return 0;

}

