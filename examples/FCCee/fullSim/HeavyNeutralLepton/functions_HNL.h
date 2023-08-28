#ifndef HNLfunctions_H
#define HNLfunctions_H

#include <cmath>
#include <vector>
#include <math.h>

#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/ParticleIDData.h"
#include "ReconstructedParticle2MC.h"


namespace FCCAnalyses { namespace HNLfunctions {

struct sel_iso {
    sel_iso(float arg_max_iso);
    float m_max_iso = .25;
    Vec_rp operator() (Vec_rp in, Vec_f iso);
};

sel_iso::sel_iso(float arg_max_iso) : m_max_iso(arg_max_iso) {};
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  sel_iso::operator() (Vec_rp in, Vec_f iso) {
    Vec_rp result;
    result.reserve(in.size());
    for (size_t i = 0; i < in.size(); ++i) {
        auto & p = in[i];
        if (iso[i] < m_max_iso) {
            result.emplace_back(p);
        }
    }
    return result;
}
// compute the cone isolation for reco particles
struct coneIsolation {

    coneIsolation(float arg_dr_min, float arg_dr_max);
    double deltaR(double eta1, double phi1, double eta2, double phi2) { return TMath::Sqrt(TMath::Power(eta1-eta2, 2) + (TMath::Power(phi1-phi2, 2))); };

    float dr_min = 0;
    float dr_max = 0.4;
    Vec_f operator() (Vec_rp in, Vec_rp rps) ;
};

coneIsolation::coneIsolation(float arg_dr_min, float arg_dr_max) : dr_min(arg_dr_min), dr_max( arg_dr_max ) { };
Vec_f coneIsolation::coneIsolation::operator() (Vec_rp in, Vec_rp rps) {
  
    Vec_f result;
    result.reserve(in.size());

    std::vector<ROOT::Math::PxPyPzEVector> lv_reco;
    std::vector<ROOT::Math::PxPyPzEVector> lv_charged;
    std::vector<ROOT::Math::PxPyPzEVector> lv_neutral;

    for(size_t i = 0; i < rps.size(); ++i) {

        ROOT::Math::PxPyPzEVector tlv;
        tlv.SetPxPyPzE(rps.at(i).momentum.x, rps.at(i).momentum.y, rps.at(i).momentum.z, rps.at(i).energy);
        
        if(rps.at(i).charge == 0) lv_neutral.push_back(tlv);
        else lv_charged.push_back(tlv);
    }
    
    for(size_t i = 0; i < in.size(); ++i) {

        ROOT::Math::PxPyPzEVector tlv;
        tlv.SetPxPyPzE(in.at(i).momentum.x, in.at(i).momentum.y, in.at(i).momentum.z, in.at(i).energy);
        lv_reco.push_back(tlv);
    }

    
    // compute the isolation (see https://github.com/delphes/delphes/blob/master/modules/Isolation.cc#L154) 
    for (auto & lv_reco_ : lv_reco) {
    
        double sumNeutral = 0.0;
        double sumCharged = 0.0;
    
        // charged
        for (auto & lv_charged_ : lv_charged) {
    
            double dr = coneIsolation::deltaR(lv_reco_.Eta(), lv_reco_.Phi(), lv_charged_.Eta(), lv_charged_.Phi());
            if(dr > dr_min && dr < dr_max) sumCharged += lv_charged_.P();
        }
        
        // neutral
        for (auto & lv_neutral_ : lv_neutral) {
    
            double dr = coneIsolation::deltaR(lv_reco_.Eta(), lv_reco_.Phi(), lv_neutral_.Eta(), lv_neutral_.Phi());
            if(dr > dr_min && dr < dr_max) sumNeutral += lv_neutral_.P();
        }
        
        double sum = sumCharged + sumNeutral;
        double ratio= sum / lv_reco_.P();
        result.emplace_back(ratio);
    }
    return result;
}




// build the Z resonance based on the available leptons. Returns the best lepton pair compatible with the Z mass and recoil at 125 GeV
// technically, it returns a ReconstructedParticleData object with index 0 the di-lepton system, index and 2 the leptons of the pair
struct dilepton_sys {
    dilepton_sys();
    Vec_rp operator()(Vec_rp legs) ;
};

dilepton_sys::dilepton_sys() {;}

Vec_rp dilepton_sys::dilepton_sys::operator()(Vec_rp legs) {

    Vec_rp result;
    
    rp reso;
    reso.charge = 0;
    TLorentzVector reso_lv; 

    TLorentzVector leg_lv_1, leg_lv_2;
    leg_lv_1.SetXYZM(legs[0].momentum.x, legs[0].momentum.y, legs[0].momentum.z, legs[0].mass);
    leg_lv_2.SetXYZM(legs[1].momentum.x, legs[1].momentum.y, legs[1].momentum.z, legs[1].mass);
    reso.charge = legs[0].charge + legs[1].charge;
    reso_lv = leg_lv_1+leg_lv_2;
    //std::cout << legs[0].momentum.x << "  " << legs[0].charge + legs[1].charge << std::endl;
    if(reso.charge != 0) return result; // neglect non-zero charge pairs
    reso.momentum.x = reso_lv.Px();
    reso.momentum.y = reso_lv.Py();
    reso.momentum.z = reso_lv.Pz();
    reso.mass = reso_lv.M();
    result.emplace_back(reso);
    //std::cout << reso.mass << endl;
    return result;
    
}    


// returns missing energy vector, based on reco particles
Vec_rp missingEnergy(float ecm, Vec_rp in, float p_cutoff = 0.0) {
    float px = 0, py = 0, pz = 0, e = 0;
    for(auto &p : in) {
        if (std::sqrt(p.momentum.x * p.momentum.x + p.momentum.y*p.momentum.y) < p_cutoff) continue;
        px += -p.momentum.x;
        py += -p.momentum.y;
        pz += -p.momentum.z;
        e += p.energy;
    }
    
    Vec_rp ret;
    rp res;
    res.momentum.x = px;
    res.momentum.y = py;
    res.momentum.z = pz;
    res.energy = ecm-e;
    ret.emplace_back(res);
    //std::cout << "missing pt " << std::sqrt(px*px + py*py) << std::endl;
    return ret;
}


float vertex_r(edm4hep::VertexData vertices){
    return (std::sqrt( vertices.position.x*vertices.position.x + vertices.position.y*vertices.position.y ));
}

float vertex_z(edm4hep::VertexData vertices){
    return  vertices.position.z;
}

float vertex_dist(edm4hep::VertexData vertices){
    return (std::sqrt( vertices.position.x*vertices.position.x + vertices.position.y*vertices.position.y + vertices.position.z*vertices.position.z));
}

float gen_vertex_r(ROOT::VecOps::RVec<edm4hep::MCParticleData> genpart){
    return (std::sqrt( genpart[0].vertex.x*genpart[0].vertex.x + genpart[0].vertex.y*genpart[0].vertex.y  ));
}


ROOT::VecOps::RVec<float> gen_vertex_Lxyz(ROOT::VecOps::RVec<edm4hep::MCParticleData> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    float Lxyz = std::sqrt( p.vertex.x*p.vertex.x +  p.vertex.y*p.vertex.y +p.vertex.z*p.vertex.z );
    result.push_back(Lxyz);
  }
  return result;
}


/*ROOT::VecOps::RVec<edm4hep::MCParticleData> MCElectrons_from_HNL(ROOT::VecOps::RVec<edm4hep::MCParticleData> genpart, const ROOT::VecOps::RVec<int> &ind){
    
    ROOT::VecOps::RVec<edm4hep::MCParticleData> output;
    for (size_t i = 0; i < genpart.size(); ++i) {   
        auto & p = genpart[i];
        if ( abs(p.PDG) != 11) continue; 
        for (unsigned j = p.parents_begin; j != p.parents_end; ++j) {
            int index = ind.at(j);
            int pdg_parent = genpart.at(index).PDG ;
            if( pdg_parent == abs(9900012) )  output.push_back(p);
        }
    }
    return output;
}*/

}}

#endif
