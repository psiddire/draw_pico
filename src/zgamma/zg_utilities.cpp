#include "zgamma/zg_utilities.hpp"
#include "zgamma/KinZfitter.h"

namespace ZgUtilities {
  using std::string;
  using std::to_string;
  using std::cout;
  using std::endl;
  // Returns negative lepton 4-momentum
  TLorentzVector AssignL1(const Baby &b, bool gen) {
    TLorentzVector l1;
    if(gen) {
      for(size_t i = 0; i < b.mc_id()->size(); i++)
        if(b.mc_id()->at(i) == 11 ||
           b.mc_id()->at(i) == 13 ||
           b.mc_id()->at(i) == 15) {
            l1.SetPtEtaPhiM(b.mc_pt()->at(i),
                            b.mc_eta()->at(i),
                            b.mc_phi()->at(i),
                            b.mc_mass()->at(i));
            if(b.mc_mom()->at(i) == 23) break;
          }
    }
    else{
      int il(-1);
      int indx = b.llphoton_ill()->at(0);
      if(b.ll_lepid()->at(indx) == 11){
        if(b.el_charge()->at(b.ll_i1()->at(indx)) < 0) il = b.ll_i1()->at(indx);
        else                                        il = b.ll_i2()->at(indx);
        l1.SetPtEtaPhiM(b.el_pt() ->at(il),
                        b.el_eta()->at(il),
                        b.el_phi()->at(il), 0.000511);
      }
      else if(b.ll_lepid()->at(indx) == 13){
        if(b.mu_charge()->at(b.ll_i1()->at(indx)) < 0) il = b.ll_i1()->at(indx);
        else                                        il = b.ll_i2()->at(indx);
        l1.SetPtEtaPhiM(b.mu_pt() ->at(il),
                        b.mu_eta()->at(il),
                        b.mu_phi()->at(il), 0.105);
      }
    }
    return l1;
  }

  // Returns positive lepton 4-momentum
  TLorentzVector AssignL2(const Baby &b, bool gen) {
    TLorentzVector l2;
    if(gen) {
      for(size_t i = 0; i < b.mc_id()->size(); i++)
        if(b.mc_id()->at(i) == -11 ||
           b.mc_id()->at(i) == -13 ||
           b.mc_id()->at(i) == -15) {
            l2.SetPtEtaPhiM(b.mc_pt()->at(i),
                            b.mc_eta()->at(i),
                            b.mc_phi()->at(i),
                            b.mc_mass()->at(i));
            if(b.mc_mom()->at(i) == 23) break;
          }
    }
    else{
      int il(-1);
      int indx = b.llphoton_ill()->at(0);
      if(b.ll_lepid()->at(indx) == 11){
        if(b.el_charge()->at(b.ll_i1()->at(indx)) < 0) il = b.ll_i2()->at(indx);
        else                                        il = b.ll_i1()->at(indx);
        l2.SetPtEtaPhiM(b.el_pt() ->at(il),
                        b.el_eta()->at(il),
                        b.el_phi()->at(il), 0.000511);
      }
      else if(b.ll_lepid()->at(indx) == 13){
        if(b.mu_charge()->at(b.ll_i1()->at(indx)) < 0) il = b.ll_i2()->at(indx);
        else                                        il = b.ll_i1()->at(indx);
        l2.SetPtEtaPhiM(b.mu_pt() ->at(il),
                        b.mu_eta()->at(il),
                        b.mu_phi()->at(il), 0.105);
      }
    }
    return l2;
  }

  // // Returns FSR photon
  // std::vector<TLorentzVector> AssignFSRPhoton(const Baby &b, bool gen) {
  //   TLorentzVector l1, l2, ph, null;
  //   null.SetPtEtaPhiM(0, 0, 0, 0);
  //   std::vector<TLorentzVector> fsrph;
  //   fsrph.append(null); fsrph.append(null);
  //   l1 = AssignL1(b);
  //   l2 = AssignL2(b);
  //   for(size_t i = 0; i < b.photon_pt()->size(); i++) {
  //     ph.SetPtEtaPhiM(b.photon_pt() ->at(i),
  // 			 b.photon_eta()->at(i),
  // 			 b.photon_phi()->at(i), 0.0);
  //     if (l1.DeltaR(ph) < 0.4)
  // 	fsrph[0] = ph;
  //     else if (l2.DeltaR(ph) < 0.4)
  // 	fsrph[1] = ph;
  //   }
  // }

  // // Returns FSR photon
  // std::vector<TLorentzVector> AssignFSRPhoton(const Baby &b, bool gen) {
  //   TLorentzVector l1, l2, ph, null;
  //   null.SetPtEtaPhiM(0, 0, 0, 0);
  //   std::vector<TLorentzVector> fsrph;
  //   fsrph.append(null); fsrph.append(null);
  //   l1 = AssignL1(b);
  //   l2 = AssignL2(b);
  //   for(size_t i = 0; i < b.nfsrphoton(); i++) {
  //     ph.SetPtEtaPhiM(b.fsrphoton_pt() ->at(i),
  // 		      b.fsrphoton_eta()->at(i),
  // 		      b.fsrphoton_phi()->at(i), 0.0);
  //     if (b.fsrphoton_muonidx->at(i)==b.ll_i1->at(b.llphoton_ill->at(0)) && l1.DeltaR(ph) < 0.5 && ph.Pt()/l1.Pt() < 0.4)
  // 	fsrph[0] = ph;
  //     else if (b.fsrphoton_muonidx->at(i)==b.ll_i2->at(b.llphoton_ill->at(0)) && l2.DeltaR(ph) < 0.5 && ph.Pt()/l2.Pt() < 0.4)
  // 	fsrph[1] = ph;
  //   }
  // }

  // Returns negative lepton 4-momentum error
  double AssignL1Error(const Baby &b) {
    double l1Err = 0.0;
    int il(-1);
    int indx = b.llphoton_ill()->at(0);
    TLorentzVector l1;
    if(b.ll_lepid()->at(indx) == 11){
      if(b.el_charge()->at(b.ll_i1()->at(indx)) < 0) il = b.ll_i1()->at(indx);
      else                                        il = b.ll_i2()->at(indx);
      l1.SetPtEtaPhiM(b.el_pt() ->at(il),
		      b.el_eta()->at(il),
		      b.el_phi()->at(il), 0.000511);
      // l1Err = b.el_energyErr()->at(il);
      // l1Err = b.el_energyErr()->at(il) / (b.el_etPt()->at(il) + 1.0);
      l1Err = b.el_energyErr()->at(il) * l1.Pt() / l1.P();
    }
    else if(b.ll_lepid()->at(indx) == 13){
      if(b.mu_charge()->at(b.ll_i1()->at(indx)) < 0) il = b.ll_i1()->at(indx);
      else                                        il = b.ll_i2()->at(indx);
      l1Err = b.mu_ptErr()->at(il);
    }
    return l1Err;
  }

  // Returns positive lepton 4-momentum error
  double AssignL2Error(const Baby &b) {
    double l2Err = 0.0;
    int il(-1);
    int indx = b.llphoton_ill()->at(0);
    TLorentzVector l2;
    if(b.ll_lepid()->at(indx) == 11){
      if(b.el_charge()->at(b.ll_i1()->at(indx)) < 0) il = b.ll_i2()->at(indx);
      else                                        il = b.ll_i1()->at(indx);
      l2.SetPtEtaPhiM(b.el_pt() ->at(il),
                      b.el_eta()->at(il),
                      b.el_phi()->at(il), 0.000511);
      // l2Err = b.el_energyErr()->at(il);
      // l2Err = b.el_energyErr()->at(il) / (b.el_etPt()->at(il) + 1.0);
      l2Err = b.el_energyErr()->at(il) * l2.Pt() / l2.P();
    }
    else if(b.ll_lepid()->at(indx) == 13){
      if(b.mu_charge()->at(b.ll_i1()->at(indx)) < 0) il = b.ll_i2()->at(indx);
      else                                        il = b.ll_i1()->at(indx);
      l2Err = b.mu_ptErr()->at(il);
    }
    return l2Err;
  }

  // Returns Z 4-momentum
  TLorentzVector AssignZ(const Baby &b, bool gen) {
    TLorentzVector ll;
    if(gen)
        ll = AssignL1(b,gen) + AssignL2(b,gen);
    else {
      int indx = b.llphoton_ill()->at(0);
      ll.SetPtEtaPhiM(b.ll_pt() ->at(indx),
                      b.ll_eta()->at(indx),
                      b.ll_phi()->at(indx),
                      b.ll_m()  ->at(indx));
    }
    return ll;
  }

  // Returns photon 4-momentum
  TLorentzVector AssignGamma(const Baby &b, bool gen) {
    TLorentzVector gamma;
    bool FoundGamma(false);
    if(gen)
      for(size_t i = 0; i < b.mc_id()->size(); i++)
        if(b.mc_id()->at(i) == 22 && b.mc_pt()->at(i) > 5) {
          gamma.SetPtEtaPhiM(b.mc_pt()->at(i),
                             b.mc_eta()->at(i),
                             b.mc_phi()->at(i),
                             b.mc_mass()->at(i));
          FoundGamma = true;
          if(b.mc_mom()->at(i) == 25) break;
	}
    if(!FoundGamma) {
      int indx = b.llphoton_iph()->at(0);
      gamma.SetPtEtaPhiM(b.photon_pt() ->at(indx),
                         b.photon_eta()->at(indx),
                         b.photon_phi()->at(indx), 0);
    }
    return gamma;
  }

  // Returns photon 4-momentum
  TLorentzVector AssignHGG(const Baby &b) {
    TLorentzVector g1;
    TLorentzVector g2;
    g1.SetPtEtaPhiM(b.photon_pt()->at(0), b.photon_eta()->at(0), b.photon_phi()->at(0), 0);
    g2.SetPtEtaPhiM(b.photon_pt()->at(1), b.photon_eta()->at(1), b.photon_phi()->at(1), 0);
    return g1+g2;
  }

  // Returns Higgs 4-momentum
  TLorentzVector AssignH(const Baby &b, bool gen) {
    TLorentzVector h;
    if(gen) {
      bool FoundH(false);
      for(size_t i = 0; i < b.mc_id()->size(); i++)
        if(b.mc_id()->at(i) == 25) {
          FoundH = true;
          h.SetPtEtaPhiM(b.mc_pt()->at(i),
                         b.mc_eta()->at(i),
                         b.mc_phi()->at(i),
                         b.mc_mass()->at(i));
        }
      if(!FoundH)
        h = AssignZ(b,gen) + AssignGamma(b,gen);
    }
    else {
      h.SetPtEtaPhiM(b.llphoton_pt()->at(0),
                     b.llphoton_eta()->at(0),
                     b.llphoton_phi()->at(0),
                     b.llphoton_m()->at(0));
    }
    return h;
  }

  //
  // Variables used for defining kinematic angles presented in https://arxiv.org/pdf/1108.2274.pdf
  //

  // Returns 4-momentum of q1 (quark from gluon-gluon fusion)
  //  Defined in Equation 4
  TLorentzVector AssignQ1(const Baby &b, bool gen) {
    TLorentzVector h = AssignH(b,gen);
    TVector3 htran = h.BoostVector();
    htran.SetZ(0);
    h.Boost(-1*htran);
    TLorentzVector k1;
    double pz, E;
    pz = h.Pz() + h.E();
    E  = h.E()  + h.Pz();
    k1.SetPxPyPzE(0,0,pz/2, E/2);
    k1.Boost(htran);
    return k1;
  }


  // Returns 4-momentum of q2 (quark from gluon-gluon fusion)
  //  Defined in Equation 5
  TLorentzVector AssignQ2(const Baby &b, bool gen) {
    TLorentzVector k2;
    TLorentzVector h = AssignH(b,gen);
    TVector3 htran = h.BoostVector();
    htran.SetZ(0);
    h.Boost(-1*htran);
    double pz, E;
    pz = h.Pz() - h.E();
    E  = h.E()  - h.Pz();
    k2.SetPxPyPzE(0,0,pz/2, E/2);
    k2.Boost(htran);
    return k2;
  }

  // Returns magnitude of Z candidate 3-momentum
  //  Defined in Equation 7
  double lambdaZ(const Baby &b, bool gen) {
    TLorentzVector P = AssignH(b, gen);
    TLorentzVector l1 = AssignL1(b, gen);
    TLorentzVector l2 = AssignL2(b, gen);
    TLorentzVector Z = AssignZ(b, gen);
    double M = P.M(), mll = Z.M();
    return sqrt(pow(P.Dot(Z)/M,2)-pow(mll,2));
  }

  // Cosine of angle between lepton 1 and parent Z in Higgs frame
  //  Defined in Equation 13
  double cos_theta(const Baby &b, bool gen) {
    TLorentzVector P =  AssignH(b,gen);
    TLorentzVector l1 = AssignL1(b,gen);
    TLorentzVector l2 = AssignL2(b,gen);
    double M = P.M();
    double lZ = lambdaZ(b,gen);
    double ctheta = P.Dot(l1-l2)/(M*lZ);
    if(ctheta > 1) ctheta = 0.999;
    if(ctheta <-1) ctheta = -0.999;
    return ctheta;
  }

  // Cosine of angle between incoming quarks and outgoing Zs in higgs frame
  //  Defined in Equation 8
  double cos_Theta(const Baby &b, bool gen) {
    TLorentzVector H  = AssignH(b,gen);
    TLorentzVector Z  = AssignZ(b,gen);
    TLorentzVector q1 = AssignQ1(b,gen);
    TLorentzVector q2 = AssignQ2(b,gen);
    double M = H.M();
    double lZ = lambdaZ(b,gen);
    double cosTheta = Z.Dot(q1-q2)/(M*lZ);
    if(abs(cosTheta) > 1.01) cout << "ERROR: cTheta = " << cosTheta <<  endl;
    return cosTheta;
  }

  // Angle of the Z decay plane from the z-axis (defined in Equation 1) in the higgs frame
  //  Defined in Equation 21+22
  double Getphi(const Baby &b, bool gen) {
    TVector3 l1 = AssignL1(b, gen).Vect();
    TVector3 l2 = AssignL2(b, gen).Vect();
    TVector3 q1 = AssignQ1(b, gen).Vect();
    TVector3 Z  = AssignZ( b, gen).Vect();
    double cosphi, sinphi;
    cosphi = -1*l1.Cross(l2).Dot(q1.Cross(Z))/l1.Cross(l2).Mag()/q1.Cross(Z).Mag();
    sinphi = -1*l1.Cross(l2).Dot(q1)/l1.Cross(l2).Mag()/q1.Mag();
    double phi(0);
    if(abs(cosphi) > 1.01) cout << "ERROR: cphi = " << cosphi <<  endl;
    if(cosphi > 1) cosphi = 1;
    if(cosphi < -1) cosphi = -1;
    if(sinphi < 0) phi = -1*acos(cosphi);
    else           phi = acos(cosphi);
    return phi;
  }

  double KinRefit(const Baby &b) {
    KinZfitter *kinZfitter;
    bool isData = false;
    kinZfitter = new KinZfitter(isData);

    std::map<unsigned int, TLorentzVector> selectedLeptons;
    selectedLeptons[0] = AssignL1(b);
    selectedLeptons[1] = AssignL2(b);
    std::map<unsigned int, double> errorLeptons;
    errorLeptons[0] = AssignL1Error(b);
    errorLeptons[1] = AssignL2Error(b);
    std::map<unsigned int, TLorentzVector> selectedFsrMap;
    // TLorentzVector nullFourVector(0, 0, 0, 0);
    // selectedFsrMap[0] = nullFourVector;
    // selectedFsrMap[1] = nullFourVector;

    kinZfitter->Setup(selectedLeptons, selectedFsrMap, errorLeptons);
    kinZfitter->KinRefitZ1();
    double massZ1REFIT = kinZfitter->GetRefitMZ1();
    return massZ1REFIT;
  }

  std::vector<TLorentzVector> RefitP4(const Baby &b) {
    KinZfitter *kinZfitter;
    bool isData = false;
    kinZfitter = new KinZfitter(isData);

    std::map<unsigned int, TLorentzVector> selectedLeptons;
    selectedLeptons[0] = AssignL1(b);
    selectedLeptons[1] = AssignL2(b);
    std::map<unsigned int, double> errorLeptons;
    errorLeptons[0] = AssignL1Error(b);
    errorLeptons[1] = AssignL2Error(b);
    std::map<unsigned int, TLorentzVector> selectedFsrMap;
    // TLorentzVector nullFourVector(0, 0, 0, 0);
    // selectedFsrMap[0] = nullFourVector;
    // selectedFsrMap[1] = nullFourVector;

    kinZfitter->Setup(selectedLeptons, selectedFsrMap, errorLeptons);
    kinZfitter->KinRefitZ1();
    std::vector<TLorentzVector> reFit = kinZfitter->GetRefitP4s();
    return reFit;
  }


}
