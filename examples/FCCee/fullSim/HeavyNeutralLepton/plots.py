import ROOT

ROOT.gROOT.SetBatch(True)  # Run ROOT in batch mode to avoid displaying the plot

ROOT.gStyle.SetOptStat(0)

input_root_file  = ROOT.TFile("HNL_50_v2_RECO_EDM4Hep_merged.root")
c1 = ROOT.TCanvas ("c1", "", 800,800)


def plot_histo( distrib,  x_label, setlogy):
    hist = input_root_file.Get(distrib)
    hist.GetXaxis().SetTitle(x_label)
    c1.SetLogy(setlogy)
    hist.Draw()

    c1.SaveAs("plots/"+distrib+".png")

def plot_cutflow():

    hist = input_root_file.Get("cutFlow")
    hist.GetXaxis().SetTitle("selection step")
    hist.GetXaxis().SetBinLabel(1,"no sel")
    hist.GetXaxis().SetBinLabel(2,">0 iso el")
    hist.GetXaxis().SetBinLabel(3,"=2 OS el")
    hist.GetXaxis().SetBinLabel(4,"Z mass veto")
    hist.Draw()

    c1.SaveAs("plots/CutFlow.png")

plot_cutflow()
plot_histo("MC_HNL_pt", "Gen HNL p_{T} [GeV]", 0)
plot_histo("MC_HNL_theta", "Gen HNL #theta", 0)


plot_histo("MC_electrons_pt", "Gen electron p_{T} [GeV]", 0)
plot_histo("MC_electrons_theta", "Gen electron #theta", 0)


plot_histo("gen_Vertex_Lxyz_distrib_fromHNL", "HNL decay length [mm]", 1)




plot_histo("electrons_pt_cut0", "electron p_{T} [GeV]", 0)
plot_histo("electrons_theta_cut0", "electron #theta", 0)
plot_histo("electrons_iso_cut0", "electron isolation", 1)
plot_histo("electrons_no_cut0", "electrons multiplicity", 0)

plot_histo("electrons_invmass", "m_{ee} [GeV]",1)
plot_histo("missingEnergy_pt", "missing p_{T} [GeV]", 0)