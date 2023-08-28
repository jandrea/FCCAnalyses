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

bin_mcpdg = (101, -50, 50)

bins_theta_el = (50, 0, 3.1415)
bins_p_el     = (100, 0, 60) # 100 MeV bins
bins_pt_el    = (100, 0, 60) # 100 MeV bins
bins_iso = (500, 0, 0.4)
bins_count = (10, -0.5, 9.5)
bins_invmass= (100, 0, 100)
bins_missingPT= (100, 0, 60)
bins_vertex_r= (100, 0, 100)
bins_vertex_z= (100, 0, 100)
bins_vertex_dist= (100, 0, 70)


bins_vertex_x_all= (100, 0, 70)
bins_vertex_Lxyz_all= (100, 0, 70)

# build_graph function that contains the analysis logic, cuts and histograms (mandatory)
def build_graph(df, dataset):
    
    results = []
    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")
    
    print(58)
    # define some aliases to be used later on
    
    #Pandora PF particles
    df = df.Alias("ReconstructedParticles", "PandoraPFOs")
    #MC particles particles
    df = df.Alias("GenParticles", "MCParticles")
    
    df = df.Define("MC_pdgs", "MCParticle::get_pdg(MCParticles)")

    results.append(df.Histo1D(("MC_pdgs_distrib",     "", *bin_mcpdg),     "MC_pdgs"    ))


    df = df.Define("MC_electrons", "FCCAnalyses::MCParticle::sel_pdgID(11, true) (MCParticles)")
    df = df.Define("MC_electrons_status1", "FCCAnalyses::MCParticle::sel_genStatus(1) (MC_electrons)")

    df = df.Define("MC_electrons_p",     "FCCAnalyses::MCParticle::get_p(MC_electrons_status1)")
    df = df.Define("MC_electrons_pt",    "FCCAnalyses::MCParticle::get_pt(MC_electrons_status1)")
    df = df.Define("MC_electrons_theta", "FCCAnalyses::MCParticle::get_theta(MC_electrons_status1)")

    results.append(df.Histo1D(("MC_electrons_p",     "", *bins_p_el),     "MC_electrons_p"    ))
    results.append(df.Histo1D(("MC_electrons_pt",    "", *bins_pt_el),    "MC_electrons_pt"   ))
    results.append(df.Histo1D(("MC_electrons_theta", "", *bins_theta_el), "MC_electrons_theta"))
    

    ## test of particule vertex position (initial for electrons)
    df = df.Define("gen_Vertex_x_distrib_all", "FCCAnalyses::MCParticle::get_vertex_x(MCParticles)")
    results.append(df.Histo1D(("gen_Vertex_x_distrib_all", "",    *bins_vertex_x_all),    "gen_Vertex_x_distrib_all")) 
    #df = df.Define("gen_Vertex_x_distrib_electrons_status3", "FCCAnalyses::MCParticle::get_vertex_x(MC_electrons_status1)")
    df = df.Define("MC_electrons_status23", "FCCAnalyses::MCParticle::sel_genStatus(23) (MC_electrons)")
    df = df.Define("gen_Vertex_x_distrib_electrons_status23", "FCCAnalyses::MCParticle::get_vertex_x(MC_electrons_status23)")
    results.append(df.Histo1D(("gen_Vertex_x_distrib_fromHNL", "",    *bins_vertex_x_all),    "gen_Vertex_x_distrib_electrons_status23")) 

    df = df.Define("gen_vertex_Lxyz", "FCCAnalyses::HNLfunctions::gen_vertex_Lxyz(MC_electrons_status23)")
    results.append(df.Histo1D(("gen_Vertex_Lxyz_distrib_fromHNL", "",    *bins_vertex_Lxyz_all),    "gen_vertex_Lxyz")) 

    
    #df = df.Define("MCParents_Ind", "_MCParticles_parents")
    #df = df.Define("MCElectrons_from_HNL", "FCCAnalyses::HNLfunctions::MCElectrons_from_HNL(MCParticles, MCParents_Ind)")
    #df = df.Define("gen_Vertex_x_distrib_fromHNL", "FCCAnalyses::MCParticle::get_vertex_x(MCElectrons_from_HNL)")
    
    #results.append(df.Histo1D(("gen_Vertex_x_distrib_fromHNL", "",    *bins_vertex_x_all),    "gen_Vertex_x_distrib_fromHNL")) 



    ## look at HNL particles
    df = df.Define("MC_HNL", "FCCAnalyses::MCParticle::sel_pdgID(9900012, true) (MCParticles)")

    df = df.Define("MC_HNL_p",     "FCCAnalyses::MCParticle::get_p(MC_HNL)")
    df = df.Define("MC_HNL_pt",    "FCCAnalyses::MCParticle::get_pt(MC_HNL)")
    df = df.Define("MC_HNL_theta", "FCCAnalyses::MCParticle::get_theta(MC_HNL)")

    results.append(df.Histo1D(("MC_HNL_p",     "", *bins_p_el),     "MC_HNL_p"    ))
    results.append(df.Histo1D(("MC_HNL_pt",    "", *bins_pt_el),    "MC_HNL_pt"   ))
    results.append(df.Histo1D(("MC_HNL_theta", "", *bins_theta_el), "MC_HNL_theta"))
    
    df = df.Define("MC_HNL_vertex_r", "FCCAnalyses::HNLfunctions::gen_vertex_r (MC_electrons_status1)")
    results.append(df.Histo1D(("gen_Vertex_r_distrib", "",    *bins_vertex_r),    "MC_HNL_vertex_r")) 


    df = df.Define("MC_HNL_vertex_x", "FCCAnalyses::MCParticle::get_vertex_x(MC_electrons_status1)")
    results.append(df.Histo1D(("gen_Vertex_x_distrib", "",    *bins_vertex_r),    "MC_HNL_vertex_x")) 



    ##look aat reco objects
    #df = df.Define("RecoElectrons",    "ReconstructedParticle::sel_type(13, true) ( ReconstructedParticles )")
    df = df.Define("electrons_all", "ReconstructedParticle::sel_electrons(0) ( ReconstructedParticles )")
    
    
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



    #########
    ### CUT 3: Z mass window
    #########  
    df = df.Filter("electrons_invmass < 86 || electrons_invmass > 96")
    df = df.Define("cut3", "3")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut3"))


    df = df.Define("missingEnergy",       "FCCAnalyses::HNLfunctions::missingEnergy(91., ReconstructedParticles)")
    df = df.Define("missingEnergy_pt",    "FCCAnalyses::ReconstructedParticle::get_pt(missingEnergy)")

    #results.append(df.Histo1D(("missingEnergy_pt", "", *bins_missingPT), "missingEnergy_pt")) # plot it before the cut
    results.append(df.Histo1D(("missingEnergy_pt", "", *bins_missingPT), "missingEnergy_pt")) # plot it before the cut



    #########
    ### CUT 4: missing Et Cuts
    #########  
    #df = df.Filter("missingEnergy_pt > 86")
    #df = df.Define("cut4", "4")
    #results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut4"))


    #df = df.Alias("ReconstructedTracks", "SiTracks_Refitted_trackStates")
    df = df.Alias("ReconstructedTracks", "_SiTracks_Refitted_trackStates")
    df = df.Define("trackStates_electrons", "FCCAnalyses::ReconstructedParticle2Track::getRP2TRK(electrons, ReconstructedTracks)")

    #vertex reconstruction, first attempts...
    df = df.Define("VertexObject_DiElectrons",  "FCCAnalyses::VertexFitterSimple::VertexFitter ( 1, electrons, trackStates_electrons) ")
    df = df.Define("Vertex",   "VertexingUtils::get_VertexData( VertexObject_DiElectrons )")   # primary vertex, in mm


    df = df.Define("Vertex_r",      "FCCAnalyses::HNLfunctions::vertex_r(Vertex)")   # primary vertex, in mm
    df = df.Define("Vertex_z",      "FCCAnalyses::HNLfunctions::vertex_z(Vertex)")   # primary vertex, in mm
    df = df.Define("Vertex_dist",   "FCCAnalyses::HNLfunctions::vertex_dist(Vertex)")   # primary vertex, in mm

    #results.append(df.Histo1D(("Vertex_r_distrib", "",    *bins_vertex_r),    "Vertex_r")) 
    #results.append(df.Histo1D(("Vertex_z_distrib", "",    *bins_vertex_z),    "Vertex_z")) 
    #results.append(df.Histo1D(("Vertex_dist_distrib", "", *bins_vertex_dist), "Vertex_dist")) 


    return results, weightsum
