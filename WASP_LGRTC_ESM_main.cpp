// Main code for the WASP_LGRTC model, described in A computationally efficient model for probabilistic local warming projections constrained by history matching and pattern scaling

// Manuscript for Geoscientific Model Development Discussions


// Philip Goodwin1,*, Martin Leduc2, Antti-Ilari Partanen3, H. Damon Matthews4 and Alex Rogers5


//1School of Ocean and Earth Science, University of Southampton, Southampton, SO14 3ZH, UK
//2Ouranos, Montreal, Canada
//3Climate System Research, Finnish Meteorological Institute, Helsinki, Finland
//4Department of Geography, Planning and Environment, Concordia University, Montreal, Canada
//5Department of Computer Science, University of Oxford, Oxford, UK


// * Corresponding author, email: p.a.goodwin@soton.ac.uk



//The Warming, Acidification and Sea-level Projector with Local to Global Temperautre Change: WASP/LGRT (September 2019)

//WASP is an 8-box representation of the Earth system for efficient very-large ensemble simulations.
//Coded by Philip Goodwin, This version 16-September-2019


//Citations for use of code:

//This combined WASP/LGRTC model version:
//Goodwin, P., M. Leduc, A.-I. Partanen, H.D. Matthews and A. Rogers (2019, submitted), A computationally efficient model for probabilistic local warming projections constrained by history matching and pattern scaling, submitted to Geoscientific Model Development Discussions.

//Previous version of the WASP model:
//Goodwin, P., V. M. Roussenov, A. Katavouta, G. L. Foster, E. J. Rohling and R. G. Williams (2018), Pathways to 1.5 and 2 °C warming based on observational and geological constraints, Nature Geoscience, doi:10.1038/s41561-017-0054-8.

//WASP Model description:
//Goodwin, P. (2016) How historic simulation-observation discrepancy affects future warming projections in a very large model ensemble, Climate Dynamics, CLDY-D-15-00368R2, doi: 10.1007/s00382-015-2960-z.

//Utilising the input distribution for climate sensitivity, S=1/lambda, from geological data requires citation of:
//Rohling, E. J. et al - Palaeosens Project Members (2012), Making sense of Palaeoclimate sensitivity, Nature 491, 683-691, doi:10.1038/nature11574.

//The theory behind how WASP calculated atmospheric CO2 and surface temperature anomaly is given here:
//Goodwin, P., R.G. Williams and A. Ridgwell (2015), Sensitivity of climate to cumulative carbon emissions due to compensation of ocean heat and carbon uptake, Nature Geoscience, Vol. 8, p29-34. doi:10.1038/ngeo2304.

//If you use the sea level component of WASP (commented out in the code), also cite:
//Goodwin, P., I. D. Haigh, E. J. Rohling, and A. Slangen (2017), A new approach to projecting 21st century sea-level changes and extremes, Earth’s Future, 5, doi:10.1002/2016EF000508.



////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Code compiled on a Mac Darwin environment, operating system OS X version 10.11.6, using command line compilation on the Terminal.
//To compile: place the files WASP_LGRTC_ESM_main.cpp and WASP_LGRTC_ESM_functions.cpp together in a directory. Add a folder called 'RESULTS' for the results files to be written into.
//Using compiler g++ code can be compiled with command line: g++ -O3 -o filename WASP_LGRTC_ESM_main.cpp

//This code has *not* been checked across a range of C++ compilers.
//Also, the sequence of pseudo-random numbers may be compiler and platform dependent





//Main programme - include the following...
#include <cmath>
#include <iostream>
#include <fstream>
#include <random>
#include "WASP_LGRTC_ESM_functions.cpp"

//#include "scenarios.h"
//#include "terrestrial.h"
//#include "Random_normal.h"
//#include "ensemble_parameters.h"



using namespace std;

int main()
{
    //includes one of the random number generators.
    std::minstd_rand0 generator (0);
    
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //Choose the RCP forcing scenario//////////////////////////////////////////////////////////////////////////////
    int Future_scenario;
    int Spatial_scenario;
    
    //This is target warming if on AMP scenario (for Future_scenario == 5)
    double Target1 = 3.0;
    double Target_threshold = 1.95; //The warming target above which generic≥2.0°C is used for spatial patterns on an AMP scenario
    
    //Set Future_scenario = 1 for RCP2.6, =2 for RCP4.5, =3 for RCP6.0 and =4 for RCP8.5
    Future_scenario = 5;
    
    //Set Spatial_scenario = 1 for arbitrary (average spatial pattern from RCPs 2.6, 4.5 and 8.5), =2 for close to meeting or meeting Paris Agreement (average of RCP2.6 'peak warming' and RCP4.5, =3 for exceeding Paris Agreement (average of RCP4.5 and RCP8.5), =4 for RCP2.6 'peak', =5 for RCP4.5, =6 for RCP8.5.
    Spatial_scenario = 1;
    
    //Choose the number of simulations in the initial ensemble/////////////////////////////////////////////////////
    int SCENARIOS; //This is the number of simulations in the initial ensemble - from which a smaller number will be extracted by the observational consistency tests
    
    SCENARIOS = 3000000; //Number of scenarios in observation-consistent ensemble.
    
    //defining and calculating the combined-scenario LGRTC and Stdev in LGRTC
    double LGRTC[64][128];
    double Stdev_LGRTC[64][128];
    
    double LGRTC_combined_arbitrary[64][128];
    double Stdev_LGRTC_combined_arbitrary[64][128];
    
    double LGRTC_combined_close_to_Paris[64][128];
    double Stdev_LGRTC_combined_close_to_Paris[64][128];
    
    double LGRTC_combined_exceeding_Paris[64][128];
    double Stdev_LGRTC_combined_exceeding_Paris[64][128];
    
    for (int k=0; k<64; k++)
    {
        for (int m=0; m<128; m++)
        {
            //For scenarios meeting or close to meeting the Paris Climate Agreement
            LGRTC_combined_close_to_Paris[k][m] = (LGRTC_RCP26_peak[k][m] + LGRTC_RCP45[k][m])/2.0;
            Stdev_LGRTC_combined_close_to_Paris[k][m] = sqrt( (Stdev_LGRTC_RCP26_peak[k][m]*Stdev_LGRTC_RCP26_peak[k][m]) + (Stdev_LGRTC_RCP45[k][m]*Stdev_LGRTC_RCP45[k][m]) );
            
            //This section claculates whether a point is insider the valid domain, and returns Not a Number if the point lies outside the domain
            if( abs(LGRTC_RCP26_peak[k][m] - LGRTC_RCP45[k][m]) > Stdev_LGRTC_combined_close_to_Paris[k][m])
            {
                LGRTC_combined_close_to_Paris[k][m] = sqrt(-1.0);
            }
            
            //For scenarios exceeding the Paris Climate Agreement
            LGRTC_combined_exceeding_Paris[k][m] = (LGRTC_RCP85[k][m] + LGRTC_RCP45[k][m])/2.0;
            Stdev_LGRTC_combined_exceeding_Paris[k][m] = sqrt( (Stdev_LGRTC_RCP85[k][m]*Stdev_LGRTC_RCP85[k][m]) + (Stdev_LGRTC_RCP45[k][m]*Stdev_LGRTC_RCP45[k][m]) );
            
            //This section claculates whether a point is insider the valid domain, and returns Not a Number if the point lies outside the domain
            if( abs(LGRTC_RCP85[k][m] - LGRTC_RCP45[k][m]) > Stdev_LGRTC_combined_exceeding_Paris[k][m])
            {
                LGRTC_combined_exceeding_Paris[k][m] = sqrt(-1.0);
            }
            
            //For arbitrary scenarios
            LGRTC_combined_exceeding_Paris[k][m] = (LGRTC_RCP26_peak[k][m] + LGRTC_RCP85[k][m] + LGRTC_RCP45[k][m])/3.0;
            Stdev_LGRTC_combined_arbitrary[k][m] = sqrt( (Stdev_LGRTC_RCP85[k][m]*Stdev_LGRTC_RCP85[k][m]) + (Stdev_LGRTC_RCP45[k][m]*Stdev_LGRTC_RCP45[k][m]) + (Stdev_LGRTC_RCP26_peak[k][m]*Stdev_LGRTC_RCP26_peak[k][m]));
            
            //This section claculates whether a point is insider the valid domain, and returns Not a Number if the point lies outside the domain
            if( abs(LGRTC_RCP85[k][m] - LGRTC_RCP26_peak[k][m]) > Stdev_LGRTC_combined_arbitrary[k][m])
            {
                LGRTC_combined_arbitrary[k][m] = sqrt(-1.0);
            }
            
            
        }
        
    }
    

    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
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
        cout << "Generating an initial ensemble of " << SCENARIOS <<  " simulations following an AMP scenario stabilising at " << Target1 << "degrees C warming" << endl;
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //The model must be run from a location containing a folder called 'RESULTS' for the output files below to be written into//////////////
    
    
    //This generates an array of spatial warming projections readable by MATLAB. Use ESM_4.m to then plot a figure in MATLAB
    ofstream results2 ("./RESULTS/Spatial_T_projection_2081_2100.m", std::ios::out);
    
    //Sets up the output files///////////////////////////////////////////////////////////////////
    
    //The model must be run from a location containing a folder called 'RESULTS' for the output files below to be written into//////////////
    
    //Outputs results of the model at the end of the year 2016, when the observational consistency tests are applied////////////////////
    ofstream results5 ("./RESULTS/Outputs.txt", std::ios::out);
    
    //Outputs the model input values at 2016 in the observation-consistent ensemble members/////////////////////////////////////////////
    
    //Format:
    ofstream inputs ("./RESULTS/Inputs.txt", std::ios::out);
    
    //Surface temperature anomnaly relative to the 1850-1900 average
    ofstream em_warming_T_2018 ("./RESULTS/Em_Warm_T_2018.txt", std::ios::out);
    
    //Cumulative emissions (this output of cumulative emissions is to check that the file above is working correctly, and clarify what the cumulative emissions relative to 2017 are in the file above)
    ofstream em_warming_I_2017 ("./RESULTS/Em_Warm_I_2018.txt", std::ios::out);
    
    //Surface temperature anomaly relative to preindustrial///////
    ofstream annual_warming ("./RESULTS/Warming.txt", std::ios::out);
    
    
    //Note in above files the un-weighted radiative forcing is used (i.e. no efficacy weighting is applied)
    
    ofstream Adjustment_TCRE ("./RESULTS/AdjustmentTCRE_INDC_SimTSL_AP2_coeff0.txt", std::ios::out);
    
    ofstream Adjustment_Emrate ("./RESULTS/AdjustmentEmrate_INDC_SimTSL_AP2_coeff0.txt", std::ios::out);
    
    ofstream Adjustment_DT ("./RESULTS/AdjustmentDT_INDC_SimTSL_AP2_coeff0.txt", std::ios::out);
    
    ofstream Adjustment_Emtime ("./RESULTS/AdjustmentEmtime_INDC_SimTSL_AP2_coeff0.txt", std::ios::out);
    
     /*
    //Outputs specific projections of temperature anonaly and carbon emissions for the 21st century/////////////////////////////////////
    
    //Format: ensemble-member number, mean temperature anonaly relative to preindustrial at 2081-2100, change in temperature anomaly 1986-2005 to 2081-2100, change in temperature anomaly 1986-2005 to 2046-2065, cumulative emissions 2017 to 2100, cumulative emissions 2012 to 2100
    ofstream projection ("./RESULTS/Projection.txt", std::ios::out);
    
    
    //Outputs monthly global surface surface temperature anonaly from years 1970 to 2016, values tab-seperated////////////////////////////////////////////////////////////////////
    //ofstream annual ("./RESULTS/Annual.txt", std::ios::out);
    
  
   
    //The following files output annual properties by year from 1765 to 2100 for all observation-consistent ensemble members////////////
    //All anomnalies are relative to the preindustrial//////////////////////////////////////////////////////////////////////////////////
    
    //Formats: Rows are values by year from 1765 to 2100, values are separated by a tab/////////////////////////////////////////////////
    
    //Surface temperature anomaly relative to preindustrial///////
    ofstream annual_warming ("./RESULTS/Warming.txt", std::ios::out);
    //Radiative forcing from CO2
    ofstream annual_RCO2 ("./RESULTS/RCO2.txt", std::ios::out);
    //efficacy-weighted radiative forcing from aerosols and other non-Kyoto agents
    ofstream annual_RnonKyoto ("./RESULTS/RnonKyoto.txt", std::ios::out);
    //Radiative forcing from CH4, N2O and Halogens (i.e. Kyoto agents other than CO2)
    ofstream annual_RKyoto_nonCO2 ("./RESULTS/RKyoto_nonCO2.txt", std::ios::out);
    //Sea Surface Temperature anomaly relative to preindustrial
    ofstream annual_SST ("./RESULTS/SST.txt", std::ios::out);
    //Whole-ocean heat content anomaly
    ofstream annual_heat ("./RESULTS/Heat.txt", std::ios::out);
    //Heat content anonaly of the upper 700m of the ocean
    ofstream annual_heat700m ("./RESULTS/Heat700m.txt", std::ios::out);
    //Ocean carbon storage anomaly
    ofstream annual_oceanC ("./RESULTS/OceanC.txt", std::ios::out);
    //Residual terrestrial carbon storage anonaly
    ofstream annual_landC ("./RESULTS/LandC.txt", std::ios::out);
    //Compatible cumulative emissions
    ofstream annual_Em ("./RESULTS/CumulativeEm.txt", std::ios::out);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    //The following files output properties every 5PgC emitted relative to the carbon emitted at the start of 2017/////////////////////
    
    //Format: values are given every 5Pg cumulative emissions relative to start of 2017, starting at -400PgC. Values are tab-separated//
    
    //Surface temperature anomnaly relative to the 1850-1900 average
    ofstream em_warming_T_2018 ("./RESULTS/Em_Warm_T_2018.txt", std::ios::out);
    
    //Cumulative emissions (this output of cumulative emissions is to check that the file above is working correctly, and clarify what the cumulative emissions relative to 2017 are in the file above)
    ofstream em_warming_I_2017 ("./RESULTS/Em_Warm_I_2017.txt", std::ios::out);
    
    
    //The sensitivities of temperature, radiative forcing and carbon emissions/////////////////////////////
    
    //Format: rows are values given by year from 1765 to 2100. Values are tab separated////////////////////
    
    //The sensitivity of temperature anonaly to radiative forcing by year from 1765 to 2100
    ofstream DT_DR ("./RESULTS/DT_DR.txt", std::ios::out);
    //The sensitivity of temperature anonaly to cumulative emissions by year from 1765 to 2100
    ofstream DT_DI ("./RESULTS/DT_DI.txt", std::ios::out);
    //The ratio of efficacy-weighted radiative forcing to radiative forcing from CO2 by year from 1765 to 2100
    ofstream DR_DRco2 ("./RESULTS/DR_DRco2.txt", std::ios::out);
    //The sensitivity of radiative forcing from CO2 to cumulative emissions by year from 1765 to 2100
    ofstream DRco2_DI ("./RESULTS/DRco2_DI.txt", std::ios::out);
    
    //Note in above files the un-weighted radiative forcing is used (i.e. no efficacy weighting is applied)
    
     */
   
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    

    //Definition of variables - while initial values are given, many parameters are re-assigned different values later
    
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
    double TCRE_safety = 0.0;
    
    double Restore_year = 2018.0;
    int Start_year = 1765;
    int Restore_year_int = 2018;
    
    int overshoot = 0;
    //////////////////////////////////////////////////////////////////////////////
    
    int count_annual = 0;
    
    int t_step_per_year = 12; //Number of timsteps per year
    
    int data_fit = 0; //Integer that counts every time a simulation passes an observation consistency test
    int data_fit2 = 0; //Integer that counts every time a simulation misses an observational consistency test
    double error_term = 0.0;//The total error outside the 90%/95% ranges of the observational consistency tests
    
    
    const int RESTORE = 1; //if RESTORE = 1, then emissions are set to restore CO2 to RCP scenarios. If RESTORE=0 then emissions are set regardless of CO_2 values.
    
    
    
    //Equilibrium climate variables/////////////
    double a_CO2 = 5.35; //W m-2 from Myhre et al [1998] 
    double lambda = 1.18;// 1.25; // W m-2 K-1 - this is equilibrium climate parameter, equal to 1/climate sensitivity, S
    double Ib = 3500.0;  //PgC - the buffered carbon inventory of the air-sea system
    
    double ratio = 0.8; //Ratio of global mean ocean warming to atmopsheric surface warming at equilibrium
    double ratio2 = 0.6; //Ratio of SST change to sub-surface ocean temperature change at equilibrium
    double S_ocean = (1.0/lambda)*ratio; //effective SST ocean climate sensitivity (K [Wm-2]-1)

    double epsilon = 1.28 ;//The relative efficacy of ocean heat uptake relative to radiative forcing from CO2
    double epsilon_aero = 1.0; //The relative efficacy of radiative forcing from aerosols compared to radiative forcing from CO2 (and other WMGHGs)
    
    
    //Required constants////////////////////////////////
    double CO2init = 590; //initial atmoispheric CO2 (PgC)
    double PgCtoppm = 278.0/CO2init; //constant converting CO2 in PgC to ppm
    double dt = 1.0/double (t_step_per_year); //time step in years
    double area = 5.1e8*1.0e6; //surface area of planet earth over which radiative forcing acts.
    double area_ocean = 3.5e8*1.0e6; //surface area of the ocean to allow conversion of N to effective ocean heat flux (in W m-2)
    
    double spatial_uncert = 0.0; //Scales the uncertainty in the patial warming from the CMIP5 ensemble
    double T_mean = 0.0;//
    
    //Parameters linked to radiative foricng
    double RnonCO2_Uncert=0.0;
    double RnonCO2_2011 = -0.6506;
    
    double RnonCO2_Kyoto_Uncert=0.0;
    double RnonCO2_Kyoto_2011 = 0.695;
    
    double RnonCO2_Kyoto[tmax]; //The non-CO2 from the Kyoto protocol gasses
    
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
    
    //Fraction of total ocean heat uptake going into the ocean
    double f_heat_ocean=0.93;
    
    double year[tmax]; //year
    
    //Carbon undersaturation inventories of the ocean boxes////////////////
    double I_Usat_mixed[tmax];
    double I_Usat_upper[tmax];
    double I_Usat_inter[tmax];
    double I_Usat_deep[tmax];
    double I_Usat_bottom[tmax];
    
    //Fluxes of carbon undersaturation between boxes///////////////////////
    double I_Usat_MixToDeep;
    double I_Usat_MixToUpper;
    double I_Usat_MixToInter;
    double I_Usat_MixToBottom;
    
    
    double N_mixed[tmax];      //Mixed layer heat uptake (W m-2 averaged across *whole planet's* surface area)
    
    double CO2_restore[tmax];  //CO2 concentration used to restore towards (e.g. to set to RCP pathway)
    
    double DeltaT[tmax];       //Mean surface temperature anonay relative to preindustrial
  
    double Heat_mixed[tmax]; //Mixed layer cumulative anthropogenic heat content
    double Heat_upper[tmax]; //Upper ocean cumulative anthropogenic heat content increase
    double Heat_inter[tmax]; //Intermediate box heat content anomaly
    double Heat_deep[tmax];  //deep ocean cumulative heat content increase
    double Heat_bottom[tmax];//bottom ocean cumulative anthrpogenic heat uptake
    
    double Em_heat = 0.0; //equivalent emission during time-step from temperature-CO2 solubility feedback (after Goodwin and Lenton, 2009 in GRL)
    
    double dCsatdT = -0.1/1.0e15; //PgC per K-1 m-3
    double I_eq_heat[tmax]; //Cumulative equivalent emissions from heating-CO2 solubility feedback
    
    double Steric_rise[tmax];   //To estimate steric contribution to sea level rise (m)
    double ice_melt_rise[tmax]; //To estimate ice_melt contribution to sea level rise (m)
    double Sea_surf_acid[tmax]; //to estimate sea surface acidification (Delta pH units)
    
    double I_veg[tmax];    //terrestrial cerbon stored in vegetation (PgC)
    double I_soil[tmax];   //terrestrial carbon stored in soil (PgC)
    double I_ter[tmax];    //total terrestrial carbon (PgC)
    
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
    double T_noise[tmax];
    T_noise[0] = 0.0;
    double gamma1_noise = 0.3;    //0.8//0.65
    double gamma2_noise = 0.4;
    double alpha_noise = 0.062;   //0.0075//0.01
    double z_noise = 0.0;
    ///////////////////////////////////////////////
    
    //////////////////////////////////////////////////////////////////////////////
    
    int T_restore_check = 0;
    int count_success = 0;
    
    for(int i = 0; i<tmax; i++) //Defines the year in terms of timesteps
    {
        year[i] = 1765.0 + double(i)*dt; 
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
    
    
    
    for (int j = 0; j<SCENARIOS; j++) //This loop runs 'SCENARIOS' number of ensemble members
    {
        output = 0;
        T_restore_check = 0;
        
        
        
        
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
                setRnonCO2_RCP3PD_K();
                convertRnonCO2(t_step_per_year); //Converts radiative forcing from other Kyoto agents to sub-annual timesteps////
                //Set RnonCO2_Kyoto gasses array
                for(int n=0; n<tmax; n++)
                {
                    RnonCO2_Kyoto[n] = RnonCO2[n]; //Radiative frocing from CH4, N2O and Halogens
                }
                setRnonCO2_RCP3PD();
                convertRnonCO2(t_step_per_year); //Converts to current sub-annual timesteps///////////////
                //Subtract the Kyoto gasses from RnonCO2, leaving only non-Kyoto gasses and aerosols + aerosol precursors
                for(int n=0; n<tmax; n++)
                {
                    RnonCO2[n] =  RnonCO2[n] - RnonCO2_Kyoto[n]; //Efficacy weighted radiative frocing from aerosols and other non-Kyoto agents
                }
                
            }
            
            if(Future_scenario == 2)
            {
                setCO2RCP45(PgCtoppm);
                convertCO2(t_step_per_year); //Converts the CO2 scenario to the sub-annual timesteps;
                setRnonCO2_RCP45_K();
                convertRnonCO2(t_step_per_year); //Converts radiative forcing from other Kyoto agents to sub-annual timesteps////
                //Set RnonCO2_Kyoto gasses array
                for(int n=0; n<tmax; n++)
                {
                    RnonCO2_Kyoto[n] = RnonCO2[n]; //Radiative frocing from CH4, N2O and Halogens
                }
                setRnonCO2_RCP45();
                convertRnonCO2(t_step_per_year); //Converts to current sub-annual timesteps///////////////
                //Subtract the Kyoto gasses from RnonCO2, leaving only non-Kyoto gasses and aerosols + aerosol precursors
                for(int n=0; n<tmax; n++)
                {
                    RnonCO2[n] =  RnonCO2[n] - RnonCO2_Kyoto[n]; //Efficacy weighted radiative frocing from aerosols and other non-Kyoto agents
                }
                
            }
            
            if(Future_scenario == 3)
            {
                setCO2RCP6(PgCtoppm);
                convertCO2(t_step_per_year); //Converts the CO2 scenario to the sub-annual timesteps;
                setRnonCO2_RCP60K();
                convertRnonCO2(t_step_per_year); //Converts radiative forcing from other Kyoto agents to sub-annual timesteps////
                //Set RnonCO2_Kyoto gasses array
                for(int n=0; n<tmax; n++)
                {
                    RnonCO2_Kyoto[n] = RnonCO2[n]; //Radiative frocing from CH4, N2O and Halogens
                }
                setRnonCO2_RCP6();
                convertRnonCO2(t_step_per_year); //Converts to current sub-annual timesteps///////////////
                //Subtract the Kyoto gasses from RnonCO2, leaving only non-Kyoto gasses and aerosols + aerosol precursors
                for(int n=0; n<tmax; n++)
                {
                    RnonCO2[n] =  RnonCO2[n] - RnonCO2_Kyoto[n]; //Efficacy weighted radiative frocing from aerosols and other non-Kyoto agents
                }
                
            }
            
            if(Future_scenario == 4)
            {
                setCO2RCP85(PgCtoppm);
                convertCO2(t_step_per_year); //Converts the CO2 scenario to the sub-annual timesteps;
                setRnonCO2_RCP85_K();
                convertRnonCO2(t_step_per_year); //Converts radiative forcing from other Kyoto agents to sub-annual timesteps////
                //Set RnonCO2_Kyoto gasses array
                for(int n=0; n<tmax; n++)
                {
                    RnonCO2_Kyoto[n] = RnonCO2[n]; //Radiative frocing from CH4, N2O and Halogens
                }
                setRnonCO2_RCP85();
                convertRnonCO2(t_step_per_year); //Converts to current sub-annual timesteps///////////////
                //Subtract the Kyoto gasses from RnonCO2, leaving only non-Kyoto gasses and aerosols + aerosol precursors
                for(int n=0; n<tmax; n++)
                {
                    RnonCO2[n] =  RnonCO2[n] - RnonCO2_Kyoto[n]; //Efficacy weighted radiative frocing from aerosols and other non-Kyoto agents
                }
                
            }
            
            if(Future_scenario == 5)
            {
                setCO2RCP3PD(PgCtoppm);
                convertCO2(t_step_per_year); //Converts the CO2 scenario to the sub-annual timesteps;
                
                if(Target1 < 2.01)
                {
                    setRnonCO2_RCP3PD_K();
                }
                else
                {
                    setRnonCO2_RCP45_K();
                }
                convertRnonCO2(t_step_per_year); //Converts radiative forcing from other Kyoto agents to sub-annual timesteps////
                //Set RnonCO2_Kyoto gasses array
                for(int n=0; n<tmax; n++)
                {
                    RnonCO2_Kyoto[n] = RnonCO2[n]; //Radiative frocing from CH4, N2O and Halogens
                }
                
                if(Target1 < Target_threshold)
                {
                    setRnonCO2_RCP3PD();
                }
                else
                {
                    setRnonCO2_RCP45();
                }
                convertRnonCO2(t_step_per_year); //Converts to current sub-annual timesteps///////////////
                //Subtract the Kyoto gasses from RnonCO2, leaving only non-Kyoto gasses and aerosols + aerosol precursors
                for(int n=0; n<tmax; n++)
                {
                    RnonCO2[n] =  RnonCO2[n] - RnonCO2_Kyoto[n]; //Efficacy weighted radiative frocing from aerosols and other non-Kyoto agents
                }
            }
            
            
            
            //Calculate 2011 values for the standard case with no uncertainty and scaling to match observations///
            RnonCO2_Kyoto_2011 = RnonCO2_Kyoto[(2011-1765)*t_step_per_year];
            RnonCO2_2011 = RnonCO2[(2011-1765)*t_step_per_year];
            //////////////////////////////////////////////////////////////////////////////////////////////////////
            
            ////////////////////////////////////////////////////////////
            
            //The section chooses pseudo-random values for internal model parameters, such that each of the 'SCENARIOS' simulations is different
            
            //Radnom normal value for radiative forcing coefficient
            a_CO2 = getRandomNormal(5.35,0.27);//ensemble_aCO2(j) ; //
            
            //lambda is the equilibrium climate parameter in W m-2 K-1
            //Randomly select values for lambda from the palaeo prior distribution of lambda from log-normal S distribution from Palaeosens
            double random_palaeosens1 = 0.0;
            random_palaeosens1 = getRandomLinear2(0.0, 1.0027623323, generator());
            lambda = get_lambda_palaeosens1(random_palaeosens1);//ensemble_lambda(j); //
            ///////////////////////////////////////////////////////////////////////////////////
            
            
            dNPPdT = getRandomLinear2(-5.0, 1.0, generator()); //ensemble_dNPPdT(j); //getRandomLinear2(-5.0, 1.0, generator());
            gamma_K = getRandomLinear2(0, 1.0, generator());// ensemble_gammaK(j); //
            dtaudT = getRandomLinear2(-2.0, 1.0, generator());//ensemble_dtaudT(j); //
            ratio = getRandomLinear2(0.30, 1.45, generator());//ensemble_ratio(j); ////(0.25,1.1);
            ratio2 = getRandomLinear2(0.01, 0.75, generator());//ensemble_ratio2(j); //
            epsilon = getRandomNormal(1.28,0.25); //ensemble_epsilon(j); //getRandomLinear2(0.33, 3.0, generator());
            epsilon_aero = getRandomLinear2(0.33, 3.0, generator());//ensemble_epsilon_aero(j); //
            
            f_heat_ocean = getRandomLinear2(0.9, 0.96, generator());//ensemble_f_heat_ocean(j); //
            
            tau_C_mixed = getRandomLinear2(0.1, 0.5, generator());//ensemble_tau_C_mixed(j); //
            tau_C_upper = getRandomLinear2(5.0, 40.0, generator());//ensemble_tau_C_upper(j); //
            tau_C_inter = getRandomLinear2(15.0, 60.0, generator());//ensemble_tau_C_inter(j); //
            tau_C_deep  = getRandomLinear2(75.0, 500.0, generator());//ensemble_tau_C_deep(j); //
            tau_C_bottom =getRandomLinear2(250.0, 1500.0, generator());//ensemble_tau_C_bottom(j); //
            Ib = getRandomLinear2(3100, 3900, generator());//ensemble_Ib(j);//
            
            //For random normal aerosol etc distribution to get total from aerosol
            RnonCO2_Uncert = getRandomNormal(-0.25,0.6079);//ensemble_R_nonCO2_Uncert(j);////
            RnonCO2_Kyoto_Uncert = getRandomNormal(0.32,0.06079);//ensemble_R_nonCO2_Kyoto_Uncert(j); //
            
            
            //A multiplier to give ensemble uncertainty in the LGRTC from the Standard Deviation in the LGRTC
            spatial_uncert = getRandomNormal(0.0,1.0);
            
            //If use sea level component then initalise random sea level coefficients here//////////////////
            //Sea level component described in Goodwin et al, Earth's Future 2017//////////////////////////
            
            //steric_coeff = getRandomNormal(1.5,0.3);//Random normal inputs for Williams et al (2012) in GRL approach to steric sea level rise (see Goodwin et al, Earth's Future 2017
            
            //ice_coeff = getRandomLinear(0.0,0.005); //For T, Rahmstorf (2007) coefficient in units mm K-1 yr-1
            
            //land_water = getRandomLinear(-0.01,0.09);
            
            
            //Now calculate the effective climate sensitivity of sea surface temperatures (in K [Wm-2]-1)
            S_ocean = (1.0/lambda)*ratio; //ocean climate sensitivity (K [Wm-2]-1)
            
            
            ////////////////////////////////////////////////////////////////
            
        }
        
        overshoot = 0;
        int i_last_check_year = 0;
        
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
        
        
        
        
        for (int i = 1; i<tmax; i++) 
             
        {
            
            
            //Scaling for non-CO2 forcing agents, scaled by uncertainties in 2011
            
            //First non-Kyoto and nonCO2 (e.g. aerosols)
            RnonCO2[i] = RnonCO2[i]*(1.0+(RnonCO2_Uncert/RnonCO2_2011));
           
            //Now Kyoto but nonCO2 (e.g. CH4, N2O and halogens)
            RnonCO2_Kyoto[i] = RnonCO2_Kyoto[i]*(1.0 + (RnonCO2_Kyoto_Uncert/RnonCO2_Kyoto_2011) );
                                                                                                                                
            
            
            //Noise in surface warming///////////////////////////////////////////////////
            z_noise = getRandomLinear2(-1.0, 1.0, generator());
            
            if(i<=2*t_step_per_year)
            {
                T_noise[i] = (gamma1_noise * T_noise[i-1]) + (alpha_noise * z_noise);
                
            }
            if(i>2*t_step_per_year)
            {
                
                T_noise[i] = (gamma1_noise * T_noise[i-1]) + (gamma2_noise*T_noise[i-2]) + (alpha_noise*z_noise);
            }
            /////////////////////////////////////////////////////////////////////////////
            
        
            //Carbon component////////////////////////////////////////////////////
            
            //Terrestrial Carbon Component////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            
            dIvegdt = getdCvdt(CO2[i-1]*PgCtoppm, CO2init*PgCtoppm, DeltaT[i-1], I_veg[i-1], I_veg_init, NPP_init, gamma_K, dNPPdT);
            
            dIsoildt = getdCsdt(I_veg[i-1], I_veg_init, I_soil[i-1], I_soil_init, NPP_init, DeltaT[i-1], dtaudT);
            
            I_veg[i] = I_veg[i-1] + dIvegdt*dt;
            
            I_soil[i] = I_soil[i-1] + dIsoildt*dt;
            
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
                
                
                
                EmRate = getRateINDCs(year[i]);
                
                
                
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
                        
                        if(data_fit>=16)
                        {
                            Adjustment_TCRE << j << '\t' << TCRE_estimate;
                            Adjustment_DT << j << '\t' << DT_ave2 - DT_ave_2003_2012 + 0.78;
                            if(I_remaining > 0.0)
                            {
                                Adjustment_Emtime << j << '\t' << year[i]+t_remaining ;
                                
                            }
                            Adjustment_Emrate << j << '\t' << EmRate_Assessment;
                        }
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
            
           
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            
            
            
            //Heat component//////////////////////////////////////////////////////
            
            //The following code advects heat content anomaly round the ocean, and heats the surface ocean due to radiative forcing
            
            //Calculates the equilibrium sea surface heat content change for the instantaneous efficacy-weighted radiative forcing
            H_mixed_equil = S_ocean * ( a_CO2*(log(CO2[i]/CO2init)) + epsilon_aero*RnonCO2[i] + RnonCO2_Kyoto[i] ) * (vol_total * vol_mixed) * c_p;
            
            
            double Heat_mixed_to_upper =  + (Heat_upper[i-1] - (ratio2*Heat_mixed[i-1]*vol_upper/vol_mixed) )*exp(-(dt/((tau_C_upper)))) + (ratio2*Heat_mixed[i-1]*vol_upper/vol_mixed) -  Heat_upper[i-1] ;
            
            double Heat_mixed_to_inter =  + (Heat_inter[i-1] - (ratio2*Heat_mixed[i-1]*vol_inter/vol_mixed) )*exp(-(dt/((tau_C_inter)))) + (ratio2*Heat_mixed[i-1]*vol_inter/vol_mixed) -  Heat_inter[i-1] ;
            
            double Heat_mixed_to_deep =  + (Heat_deep[i-1] - (ratio2*Heat_mixed[i-1]*vol_deep/vol_mixed) )*exp(-(dt/(tau_C_deep))) +  (ratio2*Heat_mixed[i-1]*vol_deep/vol_mixed) - Heat_deep[i-1] ;
            
            double Heat_mixed_to_bottom = + (Heat_bottom[i-1] - (ratio2*Heat_mixed[i-1]*vol_bottom/vol_mixed) )*exp(-(dt/(tau_C_bottom))) +  (ratio2*Heat_mixed[i-1]*vol_bottom/vol_mixed) - Heat_bottom[i-1]  ; 
            
            
            //Caclulate heat imbalance from distance from equilibrium SST due to efficacy weighted radiative forcing
            //Heat imbalance = total(non-efficacy weighted) radiative forcing times fraction of sea surface heat content anomaly towards equilibrium (See Goodwi, 2016 in Climate Dynamics)
            N_mixed[i] = ( a_CO2*(log(CO2[i]/CO2init)) + RnonCO2[i] + RnonCO2_Kyoto[i] ) * ((H_mixed_equil - (Heat_mixed[i-1]  ))/H_mixed_equil) ;
            
            //At low RF, the stochastic variability to surface temperatures can cause the ratio of ((H_mixed_equil - (Heat_mixed[i-1]  ))/H_mixed_equil) to vary leading to instability. The below code removes this instability:
            if(i>1)
            {
                if(N_mixed[i] - N_mixed[i-1] > 20.0 && fabs(N_mixed[i]) > fabs(a_CO2*(log(CO2[i]/CO2init)) + RnonCO2[i] + RnonCO2_Kyoto[i]))
                {
                    N_mixed[i] = a_CO2*(log(CO2[i]/CO2init)) + RnonCO2[i] + RnonCO2_Kyoto[i];
                }
                if(N_mixed[i] - N_mixed[i-1] < -20.0 && fabs(N_mixed[i]) > fabs(a_CO2*(log(CO2[i]/CO2init)) + RnonCO2[i] + RnonCO2_Kyoto[i]))
                {
                    N_mixed[i] = a_CO2*(log(CO2[i]/CO2init)) + RnonCO2[i] + RnonCO2_Kyoto[i];
                }
            }
            ////////////////////////////////////////////////////////////////////////////////////////////////
            
            //Insert f_heat_ocean fraction of energy imbalance into surface ocean mixed layer//////////////////////////
            Heat_mixed[i] = Heat_mixed[i-1] + (f_heat_ocean*N_mixed[i] * area*60.0*60.0*24.0*365.0*dt) + ((T_noise[i] - T_noise[i-1])*vol_total*vol_mixed*c_p)
            - Heat_mixed_to_upper - Heat_mixed_to_inter - Heat_mixed_to_deep - Heat_mixed_to_bottom ;
            
            
            
            Heat_upper[i] = Heat_upper[i-1] + Heat_mixed_to_upper ;   //Exchange of Heat with mixed layer
            
            Heat_inter[i] = Heat_inter[i-1] + Heat_mixed_to_inter ; 
            
            Heat_deep[i] = Heat_deep[i-1] + Heat_mixed_to_deep;  //Exchange of Heat with mixed layer
            
            
            Heat_bottom[i] = Heat_bottom[i-1] + Heat_mixed_to_bottom;          //Exchange of Heat with mixed layer
                        
        
            
            
            //Calculate warming/////////////////////////////////////////////////////////////////////////////////
           
            DeltaT[i] = (1.0/lambda)*(1.0 - (epsilon*N_mixed[i])/(a_CO2*log(CO2[i]/CO2init) + epsilon_aero*RnonCO2[i] + RnonCO2_Kyoto[i] ) ) * ( (a_CO2/Ib)* (I_em[i] + I_Usat_mixed[i] + I_Usat_deep[i] + I_Usat_upper[i] + I_Usat_inter[i] + I_Usat_bottom[i] + I_eq_heat[i]) + epsilon_aero*RnonCO2[i] + RnonCO2_Kyoto[i] ) + T_noise[i];
            
            
            
            
            //Calculate sea level rise//////////////////////////////////////////////////////////////////////////
            
            //Williams et al coefficient
            Steric_rise[i] = Steric_rise[i-1] + (steric_coeff/1000.0)*(area/area_ocean)*(f_heat_ocean*N_mixed[i])*dt; //steric_coeff [mm yr-1 (W m-2)-1] is the sea level rise per year per unit ocean heat flux from Williams et al (2012) in GRL
            
            
            
            //Rahmstorf (2007) T related ice_melt sea level rise
            ice_melt_rise[i] = ice_melt_rise[i-1] + ice_coeff*DeltaT[i]*dt; //Modelled after Rahmstorf (2007) and Vermeer and Rahmstorf (2009) with ice_coeff in m (K.yr)-1
            
            
            //Calculate surface ocean acidification/////////////////////////////////////////////////////////////
            
            Sea_surf_acid[i] = -0.387*log(CO2[i]/CO2init) - 2.41*(-(I_Usat_mixed[i]*1.0e15/12.0)/(vol_mixed*1.3e18));       //Using DpH = -0.387Dln (CO2) is at DIC saturation, then DpH = -2.41Delta DIC (DIC in moles C per m3) to account for distance from DIC saturation. Coefficients obtained from numerical carbonate chemistry solver for characteristic seasurface values.
            
            
            ///////////////////////////////////////////////////////////////////////
            
            
            
            
            
            
            
              
            //Data check at start of 2017////////////////////////////////////////////////////////////////////
            if (output == t_step_per_year-1 && year[i] >= 2017.0-dt/2.0 && year[i] <= 2017+dt/2.0  )
            {
                
                //Start with data_fit =0, add +2 to data_fit when within best estimate of observational constraint. Add +1 to data_fit when just outside best estimate range (up to 50% either way)
                //Accept simulation as observationally consistent based on value of data_fit once all observational consistency tests are conducted
                
                data_fit = 0;
                data_fit2 = 0;
                
                //Set the error term outside the 90/95% ranges equal to zero initially
                error_term = 0.0;
                double obs_mid = 0.0;
                double obs_range = 0.0;
                
                
                //Average warming 1850-1900 to 2003-2012 : data is 0.78 [0.72 to 0.85] K
                
                double DT_1850_1900 = 0.0;
                
                for (int n=0; n<51*t_step_per_year; n++)//Add up contributions
                {
                    DT_1850_1900 += DeltaT[i-167*t_step_per_year+n];
                    
                }
                DT_1850_1900  =  DT_1850_1900/(double(51*t_step_per_year));
                
                double DT_2003_2012 = 0.0;
                
                for (int n=0; n<10*t_step_per_year; n++)//Add up contributions
                {
                    DT_2003_2012 += DeltaT[i-14*t_step_per_year+n];
                    
                }
                DT_2003_2012  =  DT_2003_2012/(double(10*t_step_per_year));
                
                double DT_1986_2005 = 0.0;
                
                for (int n=0; n<20*t_step_per_year; n++)//Add up contributions
                {
                    DT_1986_2005 += DeltaT[i-31*t_step_per_year+n];
                    
                }
                DT_1986_2005  =  DT_1986_2005/(double(20*t_step_per_year));
                
                //Global surface temperature anomaly 1850-1900 through to 2003-2012
                double DT_1850_1900_to_2003_2012 = DT_2003_2012 - DT_1850_1900;
                /*
                 if(DT_1850_1900_to_2003_2012 >= 0.72 && DT_1850_1900_to_2003_2012 <= 0.85)
                 {
                 data_fit += 2;
                 }
                 obs_mid = 0.5*(0.72+0.85);
                 obs_range = 0.85-0.72;
                 if(DT_1850_1900_to_2003_2012 < 0.72)// && DT_1850_1900_to_2003_2012 >= 0.66)
                 {
                 data_fit2 += 2;
                 error_term += ( abs(obs_mid - DT_1850_1900_to_2003_2012)  / (0.5*obs_range) ) - 1.0;
                 
                 }
                 if(DT_1850_1900_to_2003_2012 > 0.85)// && DT_1850_1900_to_2003_2012 <= 0.92)
                 {
                 data_fit2 += 2;
                 error_term += ( abs(obs_mid - DT_1850_1900_to_2003_2012)  / (0.5*obs_range) ) - 1.0;
                 }*/
                
                //Global surface temperature anomaly 1850-1900 through to 1986-2005
                double DT_1850_1900_to_1986_2005 = DT_1986_2005 - DT_1850_1900;
                
                if(DT_1850_1900_to_1986_2005 >= 0.55 && DT_1850_1900_to_1986_2005 <= 0.67)
                {
                    data_fit += 2;
                }
                obs_mid = 0.5*(0.55+0.67);
                obs_range = 0.67-0.55;
                if(DT_1850_1900_to_1986_2005 < 0.55)//
                {
                    data_fit2 += 2;
                    error_term += ( abs(obs_mid - DT_1850_1900_to_1986_2005)  / (0.5*obs_range) ) - 1.0;
                    
                }
                if(DT_1850_1900_to_1986_2005 > 0.67)//
                {
                    data_fit2 += 2;
                    error_term += ( abs(obs_mid - DT_1850_1900_to_1986_2005)  / (0.5*obs_range) ) - 1.0;
                }
                
                
                
                
                double DT_1880_1900=0.0;//Average temperature 1880 to 1900
                for (int n =0; n<21*t_step_per_year; n++)
                {
                    DT_1880_1900 += DeltaT[i-137*t_step_per_year+n];
                }
                DT_1880_1900 = DT_1880_1900/double(21*t_step_per_year);
                
                double DT_2007_2016=0.0;//Average temperature 2007 to 2016
                for (int n = 0; n<10*t_step_per_year; n++)
                {
                    DT_2007_2016 += DeltaT[i-10*t_step_per_year+n];
                }
                DT_2007_2016 = DT_2007_2016/double(10*t_step_per_year);
                
                double DT_1971_1980=0.0;
                for (int n =0; n<10*t_step_per_year; n++)
                {
                    DT_1971_1980 += DeltaT[i-46*t_step_per_year+n];
                }
                DT_1971_1980 = DT_1971_1980/double(10*t_step_per_year);
                
                double DT_1951_1960=0.0;
                for (int n =0; n<10*t_step_per_year; n++)
                {
                    DT_1951_1960 += DeltaT[i-66*t_step_per_year+n];
                }
                DT_1951_1960 = DT_1951_1960/double(10*t_step_per_year);
                
                
                //Global surface temperature anomaly 1971-1980 through to 2007-2016
                if(DT_2007_2016 - DT_1971_1980 >= 0.56 && DT_2007_2016 - DT_1971_1980 <= 0.69)
                {
                    data_fit += 2;
                    
                }
                obs_mid = 0.5*(0.56+0.69);
                obs_range = 0.69-0.56;
                if( DT_2007_2016 - DT_1971_1980 < 0.56)
                {
                    data_fit2 +=2;
                    error_term += ( abs(obs_mid - (DT_2007_2016 - DT_1971_1980))  / (0.5*obs_range) ) - 1.0;
                }
                if( DT_2007_2016 - DT_1971_1980 > 0.69)
                {
                    data_fit2 +=2;
                    error_term += ( abs(obs_mid - (DT_2007_2016 - DT_1971_1980))  / (0.5*obs_range) ) - 1.0;
                }
                
                
                //Global surface temperature anomaly 1951-1960 through to 2007-2016
                if(DT_2007_2016 - DT_1951_1960 >= 0.54 && DT_2007_2016 - DT_1951_1960 <= 0.78)
                {
                    data_fit += 2;
                    
                }
                obs_mid = 0.5*(0.54+0.78);
                obs_range = 0.78-0.54;
                if(DT_2007_2016 - DT_1951_1960 > 0.78)// && DT_2007_2016 - DT_1951_1960 <= 0.96)
                {
                    data_fit2 +=2;
                    error_term += ( abs(obs_mid - (DT_2007_2016 - DT_1951_1960))  / (0.5*obs_range) ) - 1.0;
                }
                if(DT_2007_2016 - DT_1951_1960 < 0.54 )//&& DT_2007_2016 - DT_1951_1960 >= 0.38)
                {
                    data_fit2 +=2;
                    error_term += ( abs(obs_mid - (DT_2007_2016 - DT_1951_1960))  / (0.5*obs_range) ) - 1.0;
                }
                
                
                double SST_1850_1900=0.0;//Average SS temperature 1880 to 1900
                for (int n =0; n<51*t_step_per_year; n++)
                {
                    SST_1850_1900 += Heat_mixed[i-167*t_step_per_year+n]/(vol_mixed*vol_total*c_p);
                }
                SST_1850_1900 = SST_1850_1900/double(51*t_step_per_year);
                
                double SST_2003_2012 = 0.0;//average SS temperature
                
                for (int n=0; n<10*t_step_per_year; n++)//Add up contributions
                {
                    SST_2003_2012 += Heat_mixed[i-14*t_step_per_year+n]/(vol_mixed*vol_total*c_p);
                    
                }
                SST_2003_2012  =  SST_2003_2012/(double(10*t_step_per_year));
                
                
                //SST anomaly from 1850-1900 trhough to 2003-2012
                if(SST_2003_2012 - SST_1850_1900 >= 0.56 && SST_2003_2012 - SST_1850_1900 <= 0.68)
                {
                    data_fit += 2;
                    
                }
                obs_mid = 0.5*(0.56+0.68);
                obs_range = 0.68-0.56;
                if(SST_2003_2012 - SST_1850_1900 < 0.56)
                {
                    data_fit2 += 2;
                    error_term += ( abs(obs_mid - (SST_2003_2012 - SST_1850_1900))  / (0.5*obs_range) ) - 1.0;
                }
                if(SST_2003_2012 - SST_1850_1900 > 0.68 )//&& SST_2003_2012 - SST_1850_1900 <= 0.74)
                {
                    data_fit2 += 2;
                    error_term += ( abs(obs_mid - (SST_2003_2012 - SST_1850_1900))  / (0.5*obs_range) ) - 1.0;
                }
                
                
                
                
                
                //Whiole ocean heat content change from 1971 to 2010
                double Ocean_Heat_71_10 = 0.0;
                
                Ocean_Heat_71_10 = (Heat_mixed[i-6*t_step_per_year-1] + Heat_upper[i-6*t_step_per_year-1] + Heat_inter[i-6*t_step_per_year-1] + Heat_deep[i-6*t_step_per_year-1] + Heat_bottom[i-6*t_step_per_year-1])
                -  (Heat_mixed[i-46*t_step_per_year+t_step_per_year-1] + Heat_upper[i-46*t_step_per_year+t_step_per_year-1] + Heat_inter[i-46*t_step_per_year+t_step_per_year-1] + Heat_deep[i-46*t_step_per_year+t_step_per_year-1] + Heat_bottom[i-46*t_step_per_year+t_step_per_year-1]);
                
                
                if(Ocean_Heat_71_10 >=  1.17e23 && Ocean_Heat_71_10 <= 3.32e23)
                {
                    data_fit += 2;
                    
                    
                }
                obs_mid = 0.5*(1.17e23+3.32e23);
                obs_range = 3.32e23-1.17e23;
                if(Ocean_Heat_71_10 <  1.17e23 )//&& Ocean_Heat_71_10 >= 0.095e23)
                {
                    data_fit2 += 2;
                    error_term += ( abs(obs_mid - Ocean_Heat_71_10)  / (0.5*obs_range) ) - 1.0;
                    
                }
                if(Ocean_Heat_71_10 > 3.32e23)// && Ocean_Heat_71_10 <= 4.395e23)
                {
                    data_fit2 += 2;
                    error_term += ( abs(obs_mid - Ocean_Heat_71_10)  / (0.5*obs_range) ) - 1.0;
                    
                }
                
                
                
                //The change in heat content of the upper 700m of the ocean from 1971 to 2010
                double rate700m_heat_71_10 = 0.0;
                
                
                //Now DQ700m in ZJ
                rate700m_heat_71_10 =  Heat_mixed[i-6*t_step_per_year-1] + Heat_upper[i-6*t_step_per_year-1] - Heat_mixed[i-46*t_step_per_year+t_step_per_year-1] - Heat_upper[i-46*t_step_per_year+t_step_per_year-1];
                
                if(rate700m_heat_71_10  >= 98.0e21 && rate700m_heat_71_10  <= 170.0e21)
                {
                    data_fit += 2;
                    
                }
                obs_mid = 0.5*(98.0e21+170.0e21);
                obs_range = 170.0e21-98.0e21;
                if(rate700m_heat_71_10  < 98.0e21 )//&& rate700m_heat_71_10  >= 62.5e21)
                {
                    data_fit2 += 2;
                    error_term += ( abs(obs_mid - rate700m_heat_71_10)  / (0.5*obs_range) ) - 1.0;
                }
                if(rate700m_heat_71_10  > 170.0e21 )//&& rate700m_heat_71_10  <= 206.0e21)
                {
                    data_fit2 += 2;
                    error_term += ( abs(obs_mid - rate700m_heat_71_10)  / (0.5*obs_range) ) - 1.0;
                }
                
                
                double rate700m_heat_93_10 = 0.0;
                rate700m_heat_93_10 =  ( (Heat_mixed[i-6*t_step_per_year-1] + Heat_upper[i-6*t_step_per_year-1] - Heat_mixed[i-19*t_step_per_year] - Heat_upper[i-24*t_step_per_year]) / ((2010.0-1993.0)*(60.0*60.0*24.0*365.0))) ;//J per second
                rate700m_heat_93_10 = rate700m_heat_93_10/1.0e12; //convert to terra watts
                
                //Cumulative carbon emissions 1750 to 2011: data is 555 [470 to 640] PgC
                double Emissions_1750_2011 = 0.0;
                Emissions_1750_2011 =  I_em[i-7*t_step_per_year] + (I_ter[i-7*t_step_per_year] - I_ter[0]);
                
                
                
                //Terrestrial carbon uptake 1750 to 2011: data is 160 [70 to 250] PgC
                double Terrestrial_C_1750_2011 = 0.0;
                Terrestrial_C_1750_2011 = (I_ter[i-6*t_step_per_year] - I_ter[0]);
                
                if(Terrestrial_C_1750_2011 >= 70.0 && Terrestrial_C_1750_2011 <= 250.0)
                {
                    data_fit += 2;
                    
                }
                obs_mid = 0.5*(70.0+250.0);
                obs_range = 250.0-70.0;
                if(Terrestrial_C_1750_2011 < 70.0 )//&& Terrestrial_C_1750_2011 >= -20.0)
                {
                    data_fit2 += 2;
                    error_term += ( abs(obs_mid - Terrestrial_C_1750_2011)  / (0.5*obs_range) ) - 1.0;
                }
                if(Terrestrial_C_1750_2011 > 250.0 )//&& Terrestrial_C_1750_2011 <= 340.0)
                {
                    data_fit2 += 2;
                    error_term += ( abs(obs_mid - Terrestrial_C_1750_2011)  / (0.5*obs_range) ) - 1.0;
                }
                
                
                //Terrestrial carbon uptake
                double Terrestrial_C_2000_2009 = 0.0;
                Terrestrial_C_2000_2009 = (I_ter[i-8*t_step_per_year] - I_ter[i-18*t_step_per_year]);
                
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
                
                
                
                
                //Ocean carbon uptake 1750 to 2011: data is 155 [125 to 185] PgC
                double Ocean_C_1750_2011 = 0.0;
                Ocean_C_1750_2011 = I_em[i-5*t_step_per_year]  - (CO2[i-5*t_step_per_year] - CO2[0]);
                
                
                if(Ocean_C_1750_2011 >= 125.0 && Ocean_C_1750_2011 <= 185.0)
                {
                    data_fit += 2;
                    
                }
                obs_mid = 0.5*(125.0+185.0);
                obs_range = 185.0-125.0;
                if(Ocean_C_1750_2011 < 125.0 )//&& Ocean_C_1750_2011 >= 95.0)
                {
                    data_fit2 += 2;
                    error_term += ( abs(obs_mid - Ocean_C_1750_2011)  / (0.5*obs_range) ) - 1.0;
                    
                }
                if(Ocean_C_1750_2011 > 185.0 )//&& Ocean_C_1750_2011 <= 215.0)
                {
                    data_fit2 += 2;
                    error_term += ( abs(obs_mid - Ocean_C_1750_2011)  / (0.5*obs_range) ) - 1.0;
                    
                }
                
                
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
                
                //if data_fit2 + data_fit==18 is there to prevent runs with numerical instability passing the test
                if(error_term<0.1 && data_fit + data_fit2 == 18)
                {
                    //cout << data_fit << '\t' << error_term << endl;
                    data_fit = 18;
                }
                
                
                //Only accept the ensemble member into the data-consistent ensemble if data_fit is 17 or more
                //datafit <= m, where m should be tuned to '2*n - 1', where n is the number of observational constraints applied - see Goodwin (2016).
                if(data_fit >= 18)
                {
                    //prints to screen when observationally consistent ensemble member is found, and what number simulation is observationally consistent
                    cout << "Simulation" << '\t' <<  j << '\t' << "is observation consistent" << endl;
                    
                    //Writes the historic simulated climate properties out to file
                    results5 << DT_1850_1900_to_1986_2005 << '\t' << DT_2007_2016 - DT_1971_1980 << '\t' <<  DT_2007_2016 - DT_1951_1960 << '\t' << SST_2003_2012 - SST_1850_1900 << '\t' << rate700m_heat_71_10 << '\t'  << Ocean_Heat_71_10 << '\t' << Ocean_C_1750_2011 << '\t' << Terrestrial_C_1750_2011 << '\t' << Terrestrial_C_2000_2009/10.0 << endl;
                    
                    //Writes model input parameters to file/////////////////////////////////////////
                    inputs << lambda << '\t' <<  ratio << '\t' << ratio2 << '\t' <<  epsilon << '\t' << epsilon_aero << '\t' << tau_C_mixed << '\t' << tau_C_upper << '\t' << tau_C_inter <<'\t' << tau_C_deep << '\t' << tau_C_bottom << '\t' << Ib << '\t' << a_CO2 << '\t' << RnonCO2_Uncert << '\t' << RnonCO2_Kyoto_Uncert << '\t' << (a_CO2*log(CO2[i-6*t_step_per_year]/CO2init) + RnonCO2[i-6*t_step_per_year] + RnonCO2_Kyoto[i-6*t_step_per_year] ) << '\t' << a_CO2*log(CO2[i-6*t_step_per_year]/CO2init) << '\t' << RnonCO2_Kyoto[i-6*t_step_per_year] << '\t' << RnonCO2[i-6*t_step_per_year] << '\t' << f_heat_ocean << '\t' << dNPPdT << '\t' << dtaudT << '\t' << gamma_K << endl;
                    
                    
                    
                }
                
                
                
                
                
                //If not data-consistent then end simulation and move on to the next - no need to run past 2017 if not data-consistent
               if(data_fit <18)
               {
                   i = tmax - 1;
                   
               }
                     
               
                                         
            }
                                     
                                     
                                     
           
            T_mean = 0.0;
            
            //Output at 2100 only (for multi-parameter perturbations) with 20-year averages for quantities///////
            if (output == t_step_per_year-1 && year[i] >= 2100.0-0.05 && year[i] <= 2100+0.05 )
            {
                
                
                //Calculate 20-year averages for warming and sea level
                double DT_86_05=0.0;
                double DT_46_65=0.0;
                double DT_81_00=0.0;
                
                
                double DT_ave_2003_2012=0.0;
                
                for (int n=0; n<20*t_step_per_year; n++)//Add up contributions
                {
                    DT_86_05 += DeltaT[i-(114*t_step_per_year-n)];
                    
                    
                    DT_46_65 += DeltaT[i-(55*t_step_per_year-n)];
                    
                    DT_81_00 += DeltaT[i-(19*t_step_per_year-n)];
                    
                    
                    
                    if(n == 20*t_step_per_year-1)//At the end divide by 20 to get average
                    {
                        DT_86_05 = DT_86_05/(20.0*t_step_per_year);
                        
                        DT_46_65 = DT_46_65/(20.0*t_step_per_year);
                        
                        DT_81_00 = DT_81_00/(20.0*t_step_per_year);
                        
                    }
                    
                }
                
                for(int n = 0; n< 10*t_step_per_year; n++)
                {
                    DT_ave_2003_2012 += DeltaT[i-97*t_step_per_year+n];
                }
                DT_ave_2003_2012 = DT_ave_2003_2012/(10.0*double(t_step_per_year));
                
                int count_em = 0;
                
                double DT_50_00=0.0;
                
                for (int n=0; n<50*t_step_per_year; n++)
                {
                    DT_50_00 += DeltaT[i-(250*t_step_per_year-n)];
                    
                    if(n==50*t_step_per_year-1)//At the end divide by 50 to get average
                    {
                        DT_50_00 = DT_50_00/(50.0*t_step_per_year);
                    }
                }
                
                T_mean += DT_81_00 - DT_50_00;
                
                if( j == 0)
                {
                    results2 << "T_2081_2100= [";
                    
                }
                if(j>0)
                {
                    
                    results2 << "T_2081_2100(:,:," << count_success+1 << ")= [" ;
                }
                count_success ++;
                //results << j << '\t' <<  DeltaT[i]  ;
                for (int lat = 0; lat < 64; lat++)
                {
                    for(int longitude = 0; longitude < 128; longitude++)
                    {
                        if(Future_scenario == 5 && Target1 < Target_threshold)
                        {
                            //Use "Close to Paris" spatial pattern
                            LGRTC[lat][longitude] = LGRTC_combined_close_to_Paris[lat][longitude];
                            Stdev_LGRTC[lat][longitude] = Stdev_LGRTC_combined_close_to_Paris[lat][longitude];
                            
                        }
                        if(Future_scenario == 5 && Target1 > Target_threshold)
                        {
                            //Use "Exceeding Paris" spatial pattern
                            LGRTC[lat][longitude] = LGRTC_combined_exceeding_Paris[lat][longitude];
                            Stdev_LGRTC[lat][longitude] = Stdev_LGRTC_combined_exceeding_Paris[lat][longitude];
                    
                        }
                        if(Future_scenario == 4)
                        {
                            //Use RCP8.5 pattern
                            LGRTC[lat][longitude] = LGRTC_RCP85[lat][longitude];
                            Stdev_LGRTC[lat][longitude] = Stdev_LGRTC_RCP85[lat][longitude];
                    
                        }
                        if(Future_scenario == 2)
                        {
                            //Use RCP4.5 pattern
                            //Use RCP8.5 pattern
                            LGRTC[lat][longitude] = LGRTC_RCP45[lat][longitude];
                            Stdev_LGRTC[lat][longitude] = Stdev_LGRTC_RCP45[lat][longitude];
                        }
                        if(Future_scenario == 1)
                        {
                            //Use 'RCP2.6 peak' pattern
                            //Use RCP8.5 pattern
                            LGRTC[lat][longitude] = LGRTC_RCP26_peak[lat][longitude];
                            Stdev_LGRTC[lat][longitude] = Stdev_LGRTC_RCP26_peak[lat][longitude];
                        }
                        if(Future_scenario == 3)
                        {
                            //Use "Exceeding Paris" pattern
                            LGRTC[lat][longitude] = LGRTC_combined_exceeding_Paris[lat][longitude];
                            Stdev_LGRTC[lat][longitude] = Stdev_LGRTC_combined_exceeding_Paris[lat][longitude];
                            
                        }
                    }
                }
                
                for (int lat = 0; lat < 64; lat++)
                {
                    for(int longitude = 0; longitude < 128; longitude++)
                    {
                        
                        results2 << '\t' <<  (DT_81_00 - DT_50_00) * (LGRTC[lat][longitude] + spatial_uncert*Stdev_LGRTC[lat][longitude]) ;
                        
                        
                        if(longitude == 127 && lat < 63)
                        {
                            
                            results2 << ";";
                        }
                        if(longitude == 127 && lat == 63)
                        {
                            
                            results2 << '\t' << "];" ;
                        }
                    }
                    
                    
                    
                    
                }
                
                results2 << endl;
                
                                if(j== SCENARIOS-1)
                {
                    T_mean = T_mean/double(SCENARIOS);
                    cout << T_mean << endl;
                }
                
                //Output stored properties by year
                for (int i2 = -250; i2 <=0; i2 ++)
                {
                    double T_ave=0.0;
                   
                    for(int n=0; n<t_step_per_year; n++)
                    {
                        T_ave += DeltaT[i+i2*t_step_per_year - n]/double(t_step_per_year);
                        
                    }
                    
                    
                    annual_warming << T_ave << '\t' ;
                    
                    
                    if(i2 == 0)
                    {
                        annual_warming << endl ;
                    }
                }
                //}
                count_annual ++;
                
                count_em = 0;
                for(int i2 = -250*t_step_per_year; i2<=0; i2++)
                {
                    //Calculate emissions as of the START of 2018 (i.e. emissions as of the last time step in 2017).
                    double Em_2017 = I_em[i-83*t_step_per_year] + I_ter[i-83*t_step_per_year] - I_ter[0];
                    double T_ave = 0.0;
                    
                    if( (I_em[i+i2] + I_ter[i+i2] - I_ter[0]) - Em_2017 > 5.0*double(count_em) - 400.0 )
                    {
                        
                        double DT_1850_1900_analysis=0.0;
                        for (int n=0; n<51*t_step_per_year; n++)//Add up contributions
                        {
                            DT_1850_1900_analysis += DeltaT[i-250*t_step_per_year+n];
                            
                        }
                        DT_1850_1900_analysis  =  DT_1850_1900_analysis/(double(51*t_step_per_year));
                        
                        //Calculates year-average temperatures centred on January 1st each year
                        //This is so when plot year-average temperatures, plotted points are on Jan 1st.
                        for(int n=0; n<t_step_per_year; n++)
                        {
                            T_ave += DeltaT[i+i2 - n + 6 ]/double(t_step_per_year);
                        }
                        em_warming_T_2018 << T_ave - DT_1850_1900_analysis << '\t';
                        em_warming_I_2017 << I_em[i+i2] + I_ter[i+i2] - I_ter[0] - Em_2017 << '\t';
                        count_em ++;
                    }
                    if(i2 == 0)
                    {
                        em_warming_T_2018 << endl;
                        em_warming_I_2017 << endl;
                    }
                }
                
                
                
                
                
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
    
    
        
   
    

    
    //Close open files
    
    
    results2.close();
    
    results5.close();
    
    
    
    return 0;
}


