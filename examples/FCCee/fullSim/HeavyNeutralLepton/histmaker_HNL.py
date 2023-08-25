# list of processes (mandatory)
processList = {
    #'p8_ee_ZZ_ecm240':{'fraction':1},
    #'p8_ee_WW_ecm240':{'fraction':1}, 
    #'wzp6_ee_mumuH_ecm240':{'fraction':1},
    #'p8_ee_WW_mumu_ecm240':    {'fraction':1, 'crossSection': 0.25792}, 
    #'p8_ee_ZZ_mumubb_ecm240':  {'fraction':1, 'crossSection': 2 * 1.35899 * 0.034 * 0.152},
    #'p8_ee_ZH_Zmumu_ecm240':   {'fraction':1, 'crossSection': 0.201868 * 0.034},
    'HNL_50_v2_RECO_EDM4Hep_merged' : {'fraction':1, 'crossSection': 2.29*10**(-8)},

}

# Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics (mandatory)
#prodTag     = "FCCee/winter2023/IDEA/"

# Link to the dictonary that contains all the cross section informations etc... (mandatory)
procDict = "FCCee_procDict_winter2023_IDEA.json"

procDictAdd={"HNL_50_0_REC_EDM4Hep":{"numberOfEvents": 50000, "sumOfWeights": 50000, "crossSection": 2.29*10**(-8), "kfactor": 1.0, "matchingEfficiency": 1.0}}


includePaths = ["functions_HNL.h"]


# additional/custom C++ functions, defined in header files (optional)
#includePaths = ["functions.h"]

# Define the input dir (optional)
#inputDir    = "outputs/FCCee/higgs/mH-recoil/mumu/stage1"
#inputDir    = "/eos/home-j/jandrea/SampleFCCee_HNL/_HNL_50_RECO_EDM4Hep/merged"
inputDir    = "/eos/home-j/jandrea/SampleFCCee_HNL/_HNL_50_v2_RECO_EDM4Hep/merged"

#Optional: output directory, default is local running directory
outputDir   = "."

print(33)

# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = -1

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 150000000 # 150 /ab



# define some binning for various histograms
#bins_theta_el = (100, -90, 90) 
bins_theta_el = (100, -5, 5)
bins_p_el     = (100, 0, 100) # 100 MeV bins
bins_pt_el    = (100, 0, 100) # 100 MeV bins
bins_iso = (500, 0, 5)
bins_count = (10, -0.5, 9.5)
bins_invmass= (100, 0, 100)




# build_graph function that contains the analysis logic, cuts and histograms (mandatory)
def build_graph(df, dataset):
    
    results = []
    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")
    
    print(58)
    # define some aliases to be used later on
    
    df = df.Alias("ReconstructedParticles", "PandoraPFOs")
    #df = df.Alias("ReconstructedParticles", "LooseSelectedPandoraPFOs")
    
    #df = df.Define("RecoElectrons",    "ReconstructedParticle::sel_type(13, true) ( ReconstructedParticles )")
    df = df.Define("electrons_all", "ReconstructedParticle::sel_electrons(1) ( ReconstructedParticles )")
    
    df = df.Define("electrons_neg", "ReconstructedParticle::sel_charge_sign(-1) ( electrons_all )")
    df = df.Define("electrons_pos", "ReconstructedParticle::sel_charge_sign( 1) ( electrons_all )")
    
    
    #df = df.Alias("Electron0", "Electron#0.index")

    # get all the leptons from the collection
    #df = df.Define("electrons_all", "FCCAnalyses::ReconstructedParticle::get(Electron0, ReconstructedParticles)")

    
    # select leptons with momentum > 20 GeV
    #df = df.Define("electrons", "FCCAnalyses::ReconstructedParticle::sel_p(5)(electrons_all)")
    
     # select leptons with momentum > 5 GeV
    df = df.Define("electrons", "FCCAnalyses::ReconstructedParticle::sel_p(5)(electrons_all)")

    # compute the electron isolation and store electrons with an isolation cut of 0.25 in a separate column electrons_sel_iso
    df = df.Define("electrons_iso",     "FCCAnalyses::HNLfunctions::coneIsolation(0.01, 0.5)(electrons, ReconstructedParticles)")
    df = df.Define("electrons_sel_iso", "FCCAnalyses::HNLfunctions::sel_iso(0.25)(electrons, electrons_iso)")

    df = df.Define("electrons_p",     "FCCAnalyses::ReconstructedParticle::get_p(electrons)")
    df = df.Define("electrons_pt",    "FCCAnalyses::ReconstructedParticle::get_pt(electrons)")
    df = df.Define("electrons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(electrons)")
    df = df.Define("electrons_no",    "FCCAnalyses::ReconstructedParticle::get_n(electrons)")
    df = df.Define("electrons_q",     "FCCAnalyses::ReconstructedParticle::get_charge(electrons)")
        
    # baseline histograms, before any selection cuts (store with _cut0)
    results.append(df.Histo1D(("electrons_p_cut0",     "", *bins_p_el),     "electrons_p"    ))
    results.append(df.Histo1D(("electrons_pt_cut0",    "", *bins_pt_el),    "electrons_pt"   ))
    results.append(df.Histo1D(("electrons_theta_cut0", "", *bins_theta_el), "electrons_theta"))
    results.append(df.Histo1D(("electrons_iso_cut0",   "", *bins_iso),      "electrons_iso"))
    results.append(df.Histo1D(("electrons_no_cut0",    "", *bins_count),    "electrons_no"))


    #########
    ### CUT 0: all events
    #########
    df = df.Define("cut0", "0")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut0"))


    #########
    ### CUT 1: at least 1 electron with at least one isolated one
    #########
    df = df.Filter("electrons_no >= 1 && electrons_iso.size() > 0")
    df = df.Define("cut1", "1")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut1"))

    #########
    ### CUT 2 :at least 2 opposite-sign (OS) isolated leptons
    #########
    df = df.Filter("electrons_no == 2 && electrons_iso.size() == 2 && abs(Sum(electrons_q)) < electrons_q.size()")
    df = df.Define("cut2", "2")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut2"))
    
    df = df.Define("dilepton_system_build", "FCCAnalyses::HNLfunctions::dilepton_sys()(electrons)")
    
    df = df.Define("electrons_invmass", "FCCAnalyses::ReconstructedParticle::get_mass(dilepton_system_build)[0]") # recoil mass


    results.append(df.Histo1D(("electrons_invmass",    "", *bins_invmass),    "electrons_invmass"))

    return results, weightsum
