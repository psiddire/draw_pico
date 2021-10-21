#include <algorithm>
#include <stdlib.h>
#include <regex>
#include "zgamma/zg_utilities.hpp"
#include "core/utilities.hpp"
#include "core/palette.hpp"
#include "core/baby.hpp"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"

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
      if(b.ll_lepid()->at(0) == 11){
        if(b.el_charge()->at(b.ll_i1()->at(0)) < 0) il = b.ll_i1()->at(0);
        else                                        il = b.ll_i2()->at(0);
        l1.SetPtEtaPhiM(b.el_pt() ->at(il),
                        b.el_eta()->at(il),
                        b.el_phi()->at(il), 0.00511);
      }
      else if(b.ll_lepid()->at(0) == 13){
        if(b.mu_charge()->at(b.ll_i1()->at(0)) < 0) il = b.ll_i1()->at(0);
        else                                        il = b.ll_i2()->at(0);
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
      if(b.ll_lepid()->at(0) == 11){
        if(b.el_charge()->at(b.ll_i1()->at(0)) < 0) il = b.ll_i2()->at(0);
        else                                        il = b.ll_i1()->at(0);
        l2.SetPtEtaPhiM(b.el_pt() ->at(il),
                        b.el_eta()->at(il),
                        b.el_phi()->at(il), 0.00511);
      }
      else if(b.ll_lepid()->at(0) == 13){
        if(b.mu_charge()->at(b.ll_i1()->at(0)) < 0) il = b.ll_i2()->at(0);
        else                                        il = b.ll_i1()->at(0);
        l2.SetPtEtaPhiM(b.mu_pt() ->at(il),
                        b.mu_eta()->at(il),
                        b.mu_phi()->at(il), 0.105);
      }
    }
    return l2;
  }

  // Returns negative lepton 4-momentum error
  double AssignL1Error(const Baby &b) {
    double l1Err = 0.0;
    int il(-1);
    if(b.ll_lepid()->at(0) == 11){
      if(b.el_charge()->at(b.ll_i1()->at(0)) < 0) il = b.ll_i1()->at(0);
      else                                        il = b.ll_i2()->at(0);
      l1Err = b.el_energyErr()->at(il);
    }
    else if(b.ll_lepid()->at(0) == 13){
      if(b.mu_charge()->at(b.ll_i1()->at(0)) < 0) il = b.ll_i1()->at(0);
      else                                        il = b.ll_i2()->at(0);
      l1Err = b.mu_ptErr()->at(il);
    }
    return l1Err;
  }

  // Returns positive lepton 4-momentum error
  double AssignL2Error(const Baby &b) {
    double l2Err = 0.0;
    int il(-1);
    if(b.ll_lepid()->at(0) == 11){
      if(b.el_charge()->at(b.ll_i1()->at(0)) < 0) il = b.ll_i2()->at(0);
      else                                        il = b.ll_i1()->at(0);
      l2Err = b.el_energyErr()->at(il);
    }
    else if(b.ll_lepid()->at(0) == 13){
      if(b.mu_charge()->at(b.ll_i1()->at(0)) < 0) il = b.ll_i2()->at(0);
      else                                        il = b.ll_i1()->at(0);
      l2Err = b.mu_ptErr()->at(il);
    }
    return l2Err;
  }

  // Returns Z 4-momentum
  TLorentzVector AssignZ(const Baby &b, bool gen) {
    TLorentzVector ll;
    if(gen)
        ll = AssignL1(b,gen) + AssignL2(b,gen);
    else
      ll.SetPtEtaPhiM(b.ll_pt() ->at(0),
                      b.ll_eta()->at(0),
                      b.ll_phi()->at(0),
                      b.ll_m()  ->at(0));
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
    if(!FoundGamma)
      gamma.SetPtEtaPhiM(b.photon_pt() ->at(0),
                         b.photon_eta()->at(0),
                         b.photon_phi()->at(0), 0);
    return gamma;
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

  double getPhErr(double CorrectedEnergy, double eta) {
    double C, S, N;
    if (abs(eta) < 1.48) {
      C = 0.35 / 100;
      S = 5.51 / 100;
      N = 98. / 1000.;
    } else {
      C = 0;
      S = 12.8 / 100;
      N = 440. / 1000.;
    }
    double result = sqrt(C * C * CorrectedEnergy * CorrectedEnergy + S * S * CorrectedEnergy + N * N);
    return result;
  }

  double KinRefit(const Baby &b) {
    double l1, l2, lph1, lph2;
    l1 = 1.0; l2 = 1.0;
    lph1 = 1.0; lph2 = 1.0;
    TLorentzVector Z1_1 = AssignL1(b);
    TLorentzVector Z1_2 = AssignL2(b);
    double RECOpT1 = Z1_1.Pt();
    double RECOpT2 = Z1_2.Pt();
    double pTerrZ1_1 = AssignL1Error(b);
    double pTerrZ1_2 = AssignL2Error(b);

    TLorentzVector Z1_ph1, Z1_ph2;
    double pTerrZ1_ph1, pTerrZ1_ph2;
    double RECOpTph1, RECOpTph2;
    TLorentzVector nullFourVector(0, 0, 0, 0);
    Z1_ph1 = nullFourVector;
    Z1_ph2 = nullFourVector;
    RECOpTph1 = 0; RECOpTph2 = 0;
    pTerrZ1_ph1 = 0; pTerrZ1_ph2 = 0;
    if(b.nfsrphoton()>=1){
      Z1_ph1.SetPtEtaPhiM(b.fsrphoton_pt()->at(0), b.fsrphoton_eta()->at(0), b.fsrphoton_phi()->at(0), 0.0);
      pTerrZ1_ph1 = getPhErr(b.fsrphoton_pt()->at(0), b.fsrphoton_eta()->at(0));
      RECOpTph1 = Z1_ph1.Pt();
    }
    if(b.nfsrphoton()==2){
      Z1_ph2.SetPtEtaPhiM(b.fsrphoton_pt()->at(1), b.fsrphoton_eta()->at(1), b.fsrphoton_phi()->at(1), 0.0);
      pTerrZ1_ph2 = getPhErr(b.fsrphoton_pt()->at(1), b.fsrphoton_eta()->at(1));
      RECOpTph2 = Z1_ph2.Pt();
    }
    cout << "Reco Mass " << (Z1_1 + Z1_2 + Z1_ph1 + Z1_ph2).M() << endl;

    RooRealVar* pT1RECO = new RooRealVar("pT1RECO", "pT1RECO", RECOpT1, 5, 500);
    RooRealVar* pT2RECO = new RooRealVar("pT2RECO", "pT2RECO", RECOpT2, 5, 500);
    double RECOpT1min = std::max(5.0, RECOpT1-2*pTerrZ1_1);
    double RECOpT2min = std::max(5.0, RECOpT2-2*pTerrZ1_2);

    RooRealVar* pTph1RECO = new RooRealVar("pTph1RECO", "pTph1RECO", RECOpTph1, 5, 500);
    RooRealVar* pTph2RECO = new RooRealVar("pTph2RECO", "pTph2RECO", RECOpTph2, 5, 500);
    double RECOpTph1min = std::max(0.5, RECOpTph1-2*pTerrZ1_ph1);
    double RECOpTph2min = std::max(0.5, RECOpTph2-2*pTerrZ1_ph2);

    // observables pT1, pT2
    RooRealVar* pT1 = new RooRealVar("pT1", "pT1FIT", RECOpT1, RECOpT1min, RECOpT1+2*pTerrZ1_1 );
    RooRealVar* pT2 = new RooRealVar("pT2", "pT2FIT", RECOpT2, RECOpT2min, RECOpT2+2*pTerrZ1_2 );
    RooRealVar* m1 = new RooRealVar("m1", "m1", Z1_1.M());
    RooRealVar* m2 = new RooRealVar("m2", "m2", Z1_2.M());
    double Vtheta1, Vphi1, Vtheta2, Vphi2;
    Vtheta1 = (Z1_1).Theta();
    Vtheta2 = (Z1_2).Theta();
    Vphi1 = (Z1_1).Phi();
    Vphi2 = (Z1_2).Phi();
    RooRealVar* theta1 = new RooRealVar("theta1", "theta1", Vtheta1);
    RooRealVar* theta2 = new RooRealVar("theta2", "theta2", Vtheta2);
    RooRealVar* phi1   = new RooRealVar("phi1", "phi1", Vphi1);
    RooRealVar* phi2   = new RooRealVar("phi2", "phi2", Vphi2);
    // dot product to calculate (p1+p2).M()
    RooFormulaVar E1("E1", "TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1)))+@2*@2)", RooArgList(*pT1, *theta1, *m1));
    RooFormulaVar E2("E2", "TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1)))+@2*@2)", RooArgList(*pT2, *theta2, *m2));

    RooRealVar* pTph1 = new RooRealVar("pTph1", "pTph1FIT", RECOpTph1, RECOpTph1min, RECOpTph1+2*pTerrZ1_ph1 );
    RooRealVar* pTph2 = new RooRealVar("pTph2", "pTph2FIT", RECOpTph2, RECOpTph2min, RECOpTph2+2*pTerrZ1_ph2 );
    double Vthetaph1, Vphiph1, Vthetaph2, Vphiph2;
    Vthetaph1 = (Z1_ph1).Theta();
    Vthetaph2 = (Z1_ph2).Theta();
    Vphiph1 = (Z1_ph1).Phi();
    Vphiph2 = (Z1_ph2).Phi();
    RooRealVar* thetaph1 = new RooRealVar("thetaph1", "thetaph1", Vthetaph1);
    RooRealVar* phiph1   = new RooRealVar("phiph1", "phiph1", Vphiph1);
    RooRealVar* thetaph2 = new RooRealVar("thetaph2", "thetaph2", Vthetaph2);
    RooRealVar* phiph2   = new RooRealVar("phiph2", "phi2", Vphiph2);
    RooFormulaVar Eph1("Eph1", "TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1))))", RooArgList(*pTph1, *thetaph1));
    RooFormulaVar Eph2("Eph2", "TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1))))", RooArgList(*pTph2, *thetaph2));

    // 3-vector DOT
    RooFormulaVar* p1v3D2 = new RooFormulaVar("p1v3D2",
    					      "@0*@1*( ((TMath::Cos(@2))*(TMath::Cos(@3)))/((TMath::Sin(@2))*(TMath::Sin(@3))))",
    					      RooArgList(*pT1, *pT2, *theta1, *theta2));
    RooFormulaVar p1D2("p1D2", "@0*@1-@2", RooArgList(E1, E2, *p1v3D2));

    // lep DOT fsrPhoton1
    RooFormulaVar* p1v3Dph1 = new RooFormulaVar("p1v3Dph1",
						"@0*@1*( (TMath::Cos(@2)*TMath::Cos(@3))/(TMath::Sin(@2)*TMath::Sin(@3))+TMath::Cos(@4-@5))",
						RooArgList(*pT1, *pTph1, *theta1, *thetaph1, *phi1, *phiph1));
    RooFormulaVar p1Dph1("p1Dph1", "@0*@1-@2",RooArgList(E1,Eph1, *p1v3Dph1));
    // 3-vector DOT
    RooFormulaVar* p2v3Dph1 = new RooFormulaVar("p2v3Dph1",
						"@0*@1*( (TMath::Cos(@2)*TMath::Cos(@3))/(TMath::Sin(@2)*TMath::Sin(@3))+TMath::Cos(@4-@5))",
						RooArgList(*pT2, *pTph1, *theta2, *thetaph1, *phi2, *phiph1));
    RooFormulaVar p2Dph1("p2Dph1", "@0*@1-@2",RooArgList(E2,Eph1, *p2v3Dph1));
    // lep DOT fsrPhoton2
    RooFormulaVar* p1v3Dph2 = new RooFormulaVar("p1v3Dph2",
						"@0*@1*( (TMath::Cos(@2)*TMath::Cos(@3))/(TMath::Sin(@2)*TMath::Sin(@3))+TMath::Cos(@4-@5))",
						RooArgList(*pT1, *pTph2, *theta1, *thetaph2, *phi1, *phiph2));
    RooFormulaVar p1Dph2("p1Dph2", "@0*@1-@2",RooArgList(E1,Eph2, *p1v3Dph2));
    // 3-vector DOT
    RooFormulaVar* p2v3Dph2 = new RooFormulaVar("p2v3Dph2",
						"@0*@1*( (TMath::Cos(@2)*TMath::Cos(@3))/(TMath::Sin(@2)*TMath::Sin(@3))+TMath::Cos(@4-@5))",
						RooArgList(*pT2, *pTph2, *theta2, *thetaph2, *phi2, *phiph2));
    RooFormulaVar p2Dph2("p2Dph2", "@0*@1-@2",RooArgList(E2,Eph2, *p2v3Dph2));
    // fsrPhoton1 DOT fsrPhoton2
    RooFormulaVar* ph1v3Dph2 = new RooFormulaVar("ph1v3Dph2",
						 "@0*@1*( (TMath::Cos(@2)*TMath::Cos(@3))/(TMath::Sin(@2)*TMath::Sin(@3))+TMath::Cos(@4-@5))",
						 RooArgList(*pTph1, *pTph2, *thetaph1, *thetaph2, *phiph1, *phiph2));
    RooFormulaVar ph1Dph2("ph1Dph2", "@0*@1-@2",RooArgList(Eph1,Eph2, *ph1v3Dph2));

    // mZ1
    RooFormulaVar* mZ1;
    mZ1 = new RooFormulaVar("mZ1", "TMath::Sqrt(2*@0+@1*@1+@2*@2)", RooArgList(p1D2, *m1, *m2));
    if(b.nfsrphoton()==1)
      mZ1 = new RooFormulaVar("mZ1", "TMath::Sqrt(2*@0+2*@1+2*@2+@3*@3+@4*@4)",
			      RooArgList(p1D2, p1Dph1, p2Dph1, *m1, *m2));
    if(b.nfsrphoton()==2)
      mZ1 = new RooFormulaVar("mZ1", "TMath::Sqrt(2*@0+2*@1+2*@2+2*@3+2*@4+2*@5+@6*@6+@7*@7)",
			      RooArgList(p1D2, p1Dph1, p2Dph1, p1Dph2, p2Dph2, ph1Dph2, *m1, *m2));

    // pTerrs, 1, 2
    RooRealVar sigmaZ1_1("sigmaZ1_1", "sigmaZ1_1", pTerrZ1_1);
    RooRealVar sigmaZ1_2("sigmaZ1_2", "sigmaZ1_2", pTerrZ1_2);
    RooRealVar sigmaZ1_ph1("sigmaZ1_ph1", "sigmaZ1_ph1", pTerrZ1_ph1);
    RooRealVar sigmaZ1_ph2("sigmaZ1_ph2", "sigmaZ1_ph2", pTerrZ1_ph2);
    // resolution for decay products
    RooGaussian gauss1("gauss1", "gaussian PDF", *pT1RECO, *pT1, sigmaZ1_1);
    RooGaussian gauss2("gauss2", "gaussian PDF", *pT2RECO, *pT2, sigmaZ1_2);
    RooGaussian gaussph1("gaussph1", "gaussian PDF", *pTph1RECO, *pTph1, sigmaZ1_ph1);
    RooGaussian gaussph2("gaussph2", "gaussian PDF", *pTph2RECO, *pTph2, sigmaZ1_ph2);

    double alphaCB_, f1_, f2_, f3_, meanCB_, meanGauss1_, meanGauss2_, meanGauss3_, nCB_, sigmaCB_, sigmaGauss1_, sigmaGauss2_, sigmaGauss3_;
    alphaCB_ = 0.0; f1_ = 0.0; f2_ = 0.0; f3_ = 0.0; meanCB_ = 0.0; meanGauss1_ = 0.0; meanGauss2_ = 0.0; meanGauss3_ = 0.0;
    nCB_ = 0.0; sigmaCB_ = 0.0; sigmaGauss1_ = 0.0; sigmaGauss2_ = 0.0; sigmaGauss3_ = 0.0;;
    if(b.ll_lepid()->at(0) == 13){
      alphaCB_ =      1.11726e+00;
      f1_ =           8.82553e-01;
      f2_ =           5.15986e-01;
      f3_ =           6.51032e-01;
      meanCB_ =       9.08773e+01;
      meanGauss1_ =   1.01245e+02;
      meanGauss2_ =   9.11730e+01;
      meanGauss3_ =   9.11915e+01;
      nCB_ =          3.24956e+00;
      sigmaCB_ =      5.73125e+00;
      sigmaGauss1_ =   8.96820e+00;
      sigmaGauss2_ =   8.87225e-01;
      sigmaGauss3_ =   2.10692e+00;
    }
    else if(b.ll_lepid()->at(0) == 11){
      alphaCB_ =      5.30189e-01;
      f1_ =           7.91603e-01;
      f2_ =           5.74403e-01;
      f3_ =           6.51020e-01;
      meanCB_ =       8.71319e+01;
      meanGauss1_ =   9.37743e+01;
      meanGauss2_ =   9.10445e+01;
      meanGauss3_ =   9.06463e+01;
      nCB_ =          9.99998e+01;
      sigmaCB_ =      5.41560e+00;
      sigmaGauss1_ =   6.69235e+00;
      sigmaGauss2_ =   9.03736e-01;
      sigmaGauss3_ =   2.34845e+00;
    }
    RooRealVar meanCB("meanCB", "", meanCB_);
    RooRealVar sigmaCB("sigmaCB", "", sigmaCB_);
    RooRealVar alphaCB("alphaCB", "", alphaCB_);
    RooRealVar nCB("nCB", "", nCB_);
    RooRealVar meanGauss1("meanGauss1", "", meanGauss1_);
    RooRealVar sigmaGauss1("sigmaGauss1", "", sigmaGauss1_);
    RooRealVar f1("f1", "", f1_);
    RooRealVar meanGauss2("meanGauss2", "", meanGauss2_);
    RooRealVar sigmaGauss2("sigmaGauss2", "", sigmaGauss2_);
    RooRealVar f2("f2", "", f2_);
    RooRealVar meanGauss3("meanGauss3", "", meanGauss3_);
    RooRealVar sigmaGauss3("sigmaGauss3", "", sigmaGauss3_);
    RooRealVar f3("f3", "", f3_);
    RooCBShape* singleCB = new RooCBShape("singleCB", "", *mZ1, meanCB, sigmaCB, alphaCB, nCB);
    RooGaussian* gaussShape1 = new RooGaussian("gaussShape1", "", *mZ1, meanGauss1, sigmaGauss1);
    RooAddPdf* CBplusGauss = new RooAddPdf("CBplusGauss", "", *singleCB, *gaussShape1, f1);
    RooGaussian* gaussShape2 = new RooGaussian("gaussShape2", "", *mZ1, meanGauss2, sigmaGauss2);
    RooAddPdf* CBplus2Gauss = new RooAddPdf("CBplus2Gauss", "", *CBplusGauss, *gaussShape2, f2);
    RooGaussian* gaussShape3 = new RooGaussian("gaussShape3", "", *mZ1, meanGauss3, sigmaGauss3);
    RooAddPdf* CBplus3Gauss = new RooAddPdf("CBplus3Gauss", "", *CBplus2Gauss, *gaussShape3, f3);
    RooProdPdf *model = new RooProdPdf("model", "", RooArgList(gauss1, gauss2, *CBplus3Gauss) );
    if(b.nfsrphoton()==1)
      model = new RooProdPdf("model", "", RooArgList(gauss1, gauss2, gaussph1, *CBplus3Gauss) );
    if(b.nfsrphoton()==2)
      model = new RooProdPdf("model", "", RooArgList(gauss1, gauss2, gaussph1, gaussph2, *CBplus3Gauss) );

    // observable set
    RooArgSet *rastmp = new RooArgSet(*pT1RECO, *pT2RECO);
    if(b.nfsrphoton()==1)
      rastmp = new RooArgSet(*pT1RECO, *pT2RECO, *pTph1RECO);
    if(b.nfsrphoton()>=2)
      rastmp = new RooArgSet(*pT1RECO, *pT2RECO, *pTph1RECO, *pTph2RECO);

    RooDataSet* pTs = new RooDataSet("pTs", "pTs", *rastmp);
    pTs->add(*rastmp);
    RooFitResult* r = model->fitTo(*pTs, RooFit::Save(), RooFit::PrintLevel(-1));
    l1 = pT1->getVal()/RECOpT1;
    l2 = pT2->getVal()/RECOpT2;
    delete r;
    delete mZ1;
    delete pT1; delete pT2;
    delete pT1RECO; delete pT2RECO;
    delete model;
    delete pTs;
    delete rastmp;
    // pT scale after refitting w.r.t. reco pT
    TLorentzVector Z1_1_True(0, 0, 0, 0);
    Z1_1_True.SetPtEtaPhiM(l1*Z1_1.Pt(), Z1_1.Eta(), Z1_1.Phi(), Z1_1.M());
    TLorentzVector Z1_2_True(0, 0, 0, 0);
    Z1_2_True.SetPtEtaPhiM(l2*Z1_2.Pt(), Z1_2.Eta(), Z1_2.Phi(), Z1_2.M());

    TLorentzVector Z1_ph1_True(0, 0, 0, 0);
    if(b.nfsrphoton()>=1){
      lph1 = pTph1->getVal()/RECOpTph1;
      Z1_ph1_True.SetPtEtaPhiM(l1*Z1_ph1.Pt(), Z1_ph1.Eta(), Z1_ph1.Phi(), Z1_ph1.M());
    }
    TLorentzVector Z1_ph2_True(0, 0, 0, 0);
    if(b.nfsrphoton()==2){
      lph2 = pTph2->getVal()/RECOpTph2;
      Z1_ph2_True.SetPtEtaPhiM(l2*Z1_ph2.Pt(), Z1_ph2.Eta(), Z1_ph2.Phi(), Z1_ph2.M());
    }
    double trueMass = 0.0;
    trueMass = (Z1_1_True + Z1_2_True + Z1_ph1_True + Z1_ph2_True).M();
    cout << "True Mass " << trueMass << endl;
    return trueMass;
  }
}
