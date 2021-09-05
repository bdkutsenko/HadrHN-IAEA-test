
//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     HadrHN.cc
//
//      Author:        V.Ivanchenko
//
//      Creation date: 14 August 2018
//
//      Modifications: 20 Jule 2021 B.Kutsenko
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>
#include <algorithm>

#include "G4Material.hh"
#include "G4ElementVector.hh"
#include "G4ProcessManager.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"

#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"

#include "G4ForceCondition.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4GRSVolume.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4StateManager.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"

#include "Histo.hh"
#include "G4Timer.hh"

#include "G4HadronNucleonXsc.hh"
#include "G4ComponentSAIDTotalXS.hh"
#include "G4ParticleInelasticXS.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronElasticXS.hh"
#include "G4ChipsNeutronElasticXS.hh"
#include "G4ChipsNeutronInelasticXS.hh"
#include "G4ChipsProtonElasticXS.hh"
#include "G4ChipsProtonInelasticXS.hh"
#include "G4ChipsPionPlusElasticXS.hh"
#include "G4ChipsPionPlusInelasticXS.hh"
#include "G4ChipsPionMinusElasticXS.hh"
#include "G4ChipsPionMinusInelasticXS.hh"
#include "G4ChipsKaonPlusElasticXS.hh"
#include "G4ChipsKaonPlusInelasticXS.hh"
#include "G4ChipsKaonMinusElasticXS.hh"
#include "G4ChipsKaonMinusInelasticXS.hh"
#include "G4BGGNucleonElasticXS.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4BGGPionInelasticXS.hh"
#include "G4PhotoNuclearCrossSection.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4eCoulombScatteringModel.hh"
#include "G4EmParameters.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4GammaNuclearXS.hh"
#include "G4GammaNuclearIAEA.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"

#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

int main(int argc, char** argv)
{
  // Control on input
  if(argc < 4) {
    G4cout << "Input parameters are not specified! Exit" << G4endl;
    return 1;
  }

  G4cout << "=======================================================" << G4endl;
  G4cout << "======   Hadron Nucleon Cross Section Test     ========" << G4endl;
  G4cout << "=======================================================" << G4endl;

  // ------- Initialisation 

  G4int     verbose  = 1, ncomponents, natoms, Z , A=0;
  G4double  energy   = 1.*MeV, density;
  G4String name = argv[2],atomName = argv[3];
  G4MaterialCutsCouple* isocouple=new G4MaterialCutsCouple();
  G4ProductionCuts* pcut = new G4ProductionCuts();
  pcut->SetProductionCut(0.0, 1);
  pcut->SetProductionCut(0.0, 2);
  pcut->SetProductionCut(0.0, 3);
  pcut->SetProductionCut(DBL_MAX, 4);

  const G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_"+atomName);
  G4MaterialCutsCouple* couple = new G4MaterialCutsCouple(mat, pcut);
  Z = mat->GetZ();
  
  if(name == "ISOXS"||name == "ISOCHIPS" || name == "ISOIAEA" ){
  G4cout<<"Enter isotop A: "<<G4endl;
  G4cin>> A;
  if(Z<3||(Z>4&&Z<11)||(Z>13&&Z<19)||(Z>31&&Z<37)||(Z>50&&Z<55)||(Z>83&&Z<87)) natoms = 2;
  else natoms =1;
  G4Isotope* iso = new G4Isotope(atomName, Z, A);
  G4Element* elIso = new G4Element("Isotop_"+atomName ,atomName, ncomponents=1);
  elIso->AddIsotope(iso, 1.0);
  G4double    pressure = atmosphere, temperature = 293*kelvin, molar_constant = CLHEP::Avogadro*CLHEP::k_Boltzmann;
  density = (elIso->GetA()*pressure)/(temperature*molar_constant);
  G4Material* mIso = new G4Material(atomName+"_"+std::to_string(natoms), density, ncomponents=1);
  mIso->AddElement(elIso, natoms);
  isocouple->SetMaterial(mIso);
  isocouple->SetProductionCuts(pcut);
  isocouple->SetIndex(1); 
  }
  else 
  couple->SetIndex(0);
 
  const G4int nidx = 11;
  const G4String particle[nidx] = {"p", "n", "pi+", "pi-", "kaon+", "kaon-", "gamma",
                                   "deuteron","triton","He3","alpha"};
  const G4ParticleDefinition* pDef[nidx] = {
    G4Proton::Proton(),G4Neutron::Neutron(),G4PionPlus::PionPlus(),G4PionMinus::PionMinus(),
    G4KaonPlus::KaonPlus(),G4KaonMinus::KaonMinus(),G4Gamma::Gamma(),
    G4Deuteron::Deuteron(),G4Triton::Triton(),G4He3::He3(),G4Alpha::Alpha() };
  std::vector<G4String> p_mod_names = {"NS", "PDG16", "SAID", "CHIPS", "BGG", "XS"};
  std::vector<G4String> n_mod_names = {"NS", "PDG16", "SAID", "CHIPS", "BGG", "XS"};
  std::vector<G4String> pip_mod_names = {"NS", "PDG16", "SAID", "CHIPS", "BGG"};
  std::vector<G4String> pin_mod_names = {"NS", "PDG16", "SAID", "CHIPS", "BGG"};
  std::vector<G4String> kp_mod_names = {"kNS", "VG", "GG", "CHIPS"};
  std::vector<G4String> kn_mod_names = {"kNS", "VG", "GG", "CHIPS"};
  std::vector<G4String> gp_mod_names = {"CHIPS","XS","IAEA","ISOCHIPS","ISOXS","ISOIAEA"};
  std::vector<G4String> ip_mod_names = {"NN","XS","ION"};
  std::vector<G4String> part_mod_names;

  const G4ParticleDefinition* project = nullptr;  
  G4String partname = argv[1];
  G4int idx = 0;
  for (; idx < nidx; idx++) {
    if (partname == particle[idx]) {
      project = pDef[idx];
      if(0 == idx) { part_mod_names = p_mod_names; }
      else if(1 == idx) { part_mod_names = n_mod_names; }
      else if(2 == idx) { part_mod_names = pip_mod_names; }
      else if(3 == idx) { part_mod_names = pin_mod_names; }
      else if(4 == idx) { part_mod_names = kp_mod_names; }
      else if(5 == idx) { part_mod_names = kn_mod_names; }
      else if(6 == idx) { part_mod_names = gp_mod_names; }
      else { part_mod_names = ip_mod_names; }
      break;
    }
  }
  if(!project) {
    G4cout << "Projectile <" << partname << "> it's not implemented yet! Exit" << G4endl;
    return 1;
  }
  
  
  G4VCrossSectionDataSet* inel = nullptr;
  G4VCrossSectionDataSet* el = nullptr;
  const G4ParticleDefinition* target;
  if(part_mod_names != gp_mod_names &&  strcmp(atomName,"H")!=0)  {
    G4cout << "Collisions with different materials implemented only for gamma. For other particles use H!\nExit" << G4endl;
    return 0;
  }
  else target = pDef[0];


  G4HadronNucleonXsc hnxs;
  G4ComponentGGNuclNuclXsc* nnxs = new G4ComponentGGNuclNuclXsc();
  G4ComponentSAIDTotalXS* said = new G4ComponentSAIDTotalXS();
  G4eCoulombScatteringModel* ss = new G4eCoulombScatteringModel(false);


  
  if(std::find(part_mod_names.begin(),part_mod_names.end(),name)==part_mod_names.end()){
    G4cout << "For projectile <" << partname <<"" "> Model <"<<name<<"> is not implemented. Exit" << G4endl;
    return 1;
  }

  G4int  isHNX = -1;
  
  if (name == "CHIPS" || name == "ISOCHIPS") {
    if (0 == idx) {
      el = new G4ChipsProtonElasticXS();
      inel = new G4ChipsProtonInelasticXS();
    }
    else if (1 == idx) {
      el = new G4ChipsNeutronElasticXS();
      inel = new G4ChipsNeutronInelasticXS();
    }
    if (2 == idx) {
      el = new G4ChipsPionPlusElasticXS();
      inel = new G4ChipsPionPlusInelasticXS();
    }
    else if (3 == idx) {
      el = new G4ChipsPionMinusElasticXS();
      inel = new G4ChipsPionMinusInelasticXS();
    }
    if (4 == idx) {
      el = new G4ChipsKaonPlusElasticXS();
      inel = new G4ChipsKaonPlusInelasticXS();
    }
    else if (5 == idx) {
      el = new G4ChipsKaonMinusElasticXS();
      inel = new G4ChipsKaonMinusInelasticXS();
    }
    else if (6 == idx) {
      inel = new G4PhotoNuclearCrossSection();
    }
  } else if (name == "BGG") {
    if (0 == idx || 1 == idx) {
      el = new G4BGGNucleonElasticXS(project);
      inel = new G4BGGNucleonInelasticXS(project);
    } else if(2 == idx || 3 == idx) {
      el = new G4BGGPionElasticXS(project);
      inel = new G4BGGPionInelasticXS(project);
    }
  }else if (name == "XS"|| name == "ISOXS") {
    if (0 == idx) {
      inel = new G4ParticleInelasticXS(G4Proton::Proton());
      el = new G4BGGNucleonElasticXS(project);
    }
    else if (1 == idx) {
      el = new G4NeutronElasticXS();
      inel = new G4NeutronInelasticXS();
    }
    else if (6 == idx) {
      inel = new G4GammaNuclearXS();
    }
    else if (7 == idx) {
      inel = new G4ParticleInelasticXS(G4Deuteron::Deuteron());
    }
    else if (8 == idx) {
      inel = new G4ParticleInelasticXS(G4Triton::Triton());
    }
    else if (9 == idx) {
      inel = new G4ParticleInelasticXS(G4He3::He3());
    }
    else if (10 == idx) {
      inel = new G4ParticleInelasticXS(G4Alpha::Alpha());
    }
  }
  else if (name == "IAEA" ||name == "ISOIAEA" ) {
    if (6 == idx) {
      inel = new G4GammaNuclearIAEA();
    }
  }
  else if(name == "NS") { 
    isHNX = 0; 
  } else if(name == "PDG16") { 
    isHNX = 2;
  } else if(name == "SAID") { 
    isHNX = 3;
  } else if(name == "kNS") { 
    isHNX = 4;
  } else if(name == "VG") { 
    isHNX = 5;
  } else if(name == "GG") { 
    isHNX = 6;
  } else if(name == "NN") { 
    isHNX = 7;
  } 

  if(isHNX < 0 && !inel) {
    G4cout << "Cross-section model <" << name 
	   << "> does not exist! Exit" << G4endl;
    return 1;
  }
  if (inel) inel->BuildPhysicsTable(*project);
  if (el) el->BuildPhysicsTable(*project);
  G4DataVector dv;
  dv.resize(4, DBL_MAX);
  G4EmParameters::Instance()->SetMscThetaLimit(0.0);
  ss->SetPolarAngleLimit(0.0);
  ss->Initialise(project, dv);

  if(name == "ISOXS" || name == "ISOCHIPS" ||name == "ISOIAEA") ss->SetCurrentCouple(isocouple);
  else  ss->SetCurrentCouple(couple);
 
  
  G4ThreeVector aPosition = G4ThreeVector(0.,0.,0.);
  G4ThreeVector aDirection = G4ThreeVector(0.0,0.0,1.0);

  G4DynamicParticle dParticle(project,aDirection,energy);

  G4cout.setf( std::ios::scientific, std::ios::floatfield );
  G4int  prec = G4cout.precision(3);

  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  partTable->SetReadiness();

  // ------- Histograms name 
  Histo    histo;
  G4String hname;
  hname = "test/" + project->GetParticleName() + "_" 
	 + atomName + "_" + name + "_" +std::to_string(Z)+"_"+std::to_string(A) ;
  G4double mass = project->GetPDGMass()/MeV;
  const G4int lnbins = 900001;
  const G4int nbins = 50001;
  G4double lxmin = 0*MeV;
  G4double lxmax = 1e5*MeV;
  G4double dlx = (lxmax - lxmin)/G4double(lnbins-1);

  G4double lowxmin = 0*MeV;
  G4double lowxmax = 200*MeV;
  G4double dlowx = (lowxmax - lowxmin)/G4double(nbins-1);
  
  G4double pmin =  1*MeV;
  G4double pmax =  1e10*GeV;
  G4double xmin =  std::log10(pmin);
  G4double xmax =  std::log10(pmax);
  G4double dx = (xmax - xmin)/G4double(nbins-1);
  
  G4double xse[nbins],xsi[nbins],xst[nbins],em[nbins], lxsi[lnbins], lowxsi[nbins];
  
  histo.Add1D("0","Inelastic",nbins,xmin-dx*0.5,xmax+dx*0.5);
  histo.Add1D("1","Elastic",nbins,xmin-dx*0.5,xmax+dx*0.5);
  histo.Add1D("2","Total",nbins,xmin-dx*0.5,xmax+dx*0.5);
  histo.Add1D("3","EM Elastic",nbins,xmin-dx*0.5,xmax+dx*0.5);
  if(partname == "gamma") histo.Add1D("4","Photonuclear Low Energy Region",lnbins,lxmin-dlx*0.5,lxmax+dlx*0.5);
  if(partname == "gamma") histo.Add1D("5","Photonuclear Lowest Energy Region",nbins,lowxmin-dlowx*0.5,lowxmax+dlowx*0.5);
  histo.SetFileName(hname);
  histo.Book();
  G4cout << "Histograms are booked output file <" << hname << "> "
	 << G4endl;

  dParticle.SetDefinition(project);
  G4double lx0 = lxmin;
  
  for (G4int i=0; i<lnbins; ++i) {
    G4double lp = lx0;
    lx0 += dlx;
    G4double le = std::sqrt(lp*lp+mass*mass)-mass;
    dParticle.SetKineticEnergy(le);
      if(isHNX <0){
      if (idx == 6 && (name == "CHIPS" || name == "XS" || name == "IAEA")){
	lxsi[i] = inel->GetElementCrossSection(&dParticle,mat->GetZ(),0)/millibarn;
      } else if (idx == 6 && (name == "ISOCHIPS" || name == "ISOXS" ||name == "ISOIAEA")){
	lxsi[i] = inel->GetIsoCrossSection(&dParticle,Z,A)/millibarn;
      } 
    }
  if(idx == 6) histo.Fill(4,lx0,lxsi[i]);
  }
  
 

  G4int i;
  dParticle.SetDefinition(project);
  G4double x0 = xmin;
  G4double lowx0 = lxmin;
  
  if(verbose>1) { 
    G4cout << "----  "
	   << project->GetParticleName() << " + " << atomName
	   << " -------------" << G4endl;
    G4cout << "  PLAB        EKIN       XSIN        XSEL      XSTOT    EM    log(p/770)  log10(p)\n"
	   << "------------------------------------------------------------------------------------"
	   << G4endl; 
  }

  for (i=0; i<nbins; ++i) {
    G4double p = std::pow(10., x0)*MeV;   
    G4double lowp = lowx0;

    lowx0 += dlowx;    
    x0 += dx;
    
    G4double e = std::sqrt(p*p+mass*mass)-mass;
    G4double lowe = std::sqrt(lowp*lowp+mass*mass)-mass;
    //if(6 == idx && e > 20.*GeV && name == "XS") { e = 20.*GeV; }
    dParticle.SetKineticEnergy(e);
    xse[i] = 0.;

    if(isHNX < 0) {
      if (idx == 6 && (name == "CHIPS" || name == "XS" || name == "IAEA")){
	xsi[i] = inel->GetElementCrossSection(&dParticle,mat->GetZ(),0)/millibarn;
	dParticle.SetKineticEnergy(lowe);
	lowxsi[i] = inel->GetElementCrossSection(&dParticle,mat->GetZ(),0)/millibarn;
      } else if (idx == 6 && (name == "ISOCHIPS" || name == "ISOXS" ||name == "ISOIAEA")){
	xsi[i] = inel->GetIsoCrossSection(&dParticle,Z,A)/millibarn;
	dParticle.SetKineticEnergy(lowe);
	lowxsi[i] = inel->GetIsoCrossSection(&dParticle,Z,A)/millibarn;

      } else if (idx > 6) {
	xsi[i] = inel->GetElementCrossSection(&dParticle,1,0)/millibarn;
      } else if (name == "CHIPS") {
	xsi[i] = inel->GetIsoCrossSection(&dParticle,1,1)/millibarn;
	xse[i] = el->GetIsoCrossSection(&dParticle,1,1)/millibarn;
      } else {
	xsi[i] = inel->GetElementCrossSection(&dParticle,1,0)/millibarn;
	xse[i] = el->GetElementCrossSection(&dParticle,1,0)/millibarn;
      }
    }
    else if (isHNX == 0) {
      hnxs.HadronNucleonXscNS(project,target,e);
      xsi[i] = hnxs.GetInelasticHadronNucleonXsc()/millibarn;
      xse[i] = hnxs.GetElasticHadronNucleonXsc()/millibarn;
    }
    else if (isHNX == 2) {
      hnxs.HadronNucleonXscPDG(project,target,e);
      xsi[i] = hnxs.GetInelasticHadronNucleonXsc()/millibarn;
      xse[i] = hnxs.GetElasticHadronNucleonXsc()/millibarn;
    }
    else if (isHNX == 3) {
      hnxs.HadronNucleonXscPDG(project,target,e);
      xsi[i] = said->GetInelasticIsotopeCrossSection(project,e,1,1)/millibarn;
      xse[i] = said->GetElasticIsotopeCrossSection(project,e,1,1)/millibarn;
    }
    else if (isHNX == 4) {
      hnxs.KaonNucleonXscNS(project,target,e);
      xsi[i] = hnxs.GetInelasticHadronNucleonXsc()/millibarn;
      xse[i] = hnxs.GetElasticHadronNucleonXsc()/millibarn; 
    }
    else if (isHNX == 5) {
      hnxs.KaonNucleonXscVG(project,target,e);
      xsi[i] = hnxs.GetInelasticHadronNucleonXsc()/millibarn;
      xse[i] = hnxs.GetElasticHadronNucleonXsc()/millibarn;
    }
    else if (isHNX == 6) {
      hnxs.KaonNucleonXscGG(project,target,e);
      xsi[i] = hnxs.GetInelasticHadronNucleonXsc()/millibarn;
      xse[i] = hnxs.GetElasticHadronNucleonXsc()/millibarn;
    }
    else if (isHNX == 7) {
      xsi[i] = nnxs->GetInelasticElementCrossSection(project,e,1,1)/millibarn;
      xse[i] = nnxs->GetElasticElementCrossSection(project,e,1,1)/millibarn;
    }
    xst[i] = xsi[i] + xse[i];

    em[i] = 
      ss->ComputeCrossSectionPerAtom(project, e, 1.0, 1.0, 0.0, e)/millibarn;
    
    if(verbose>1) { 
      G4cout << i << ". " << p << "  ekin= " << e << "   " << xsi[i] 
             << "  " << xse[i] << "  " << xst[i] << "  " << em[i]
	     << "  " << G4Log(p/770.) 
	     << "  " << std::log10(p) << G4endl;
    }

    histo.Fill(0,x0,xsi[i]);
    histo.Fill(1,x0,xse[i]);
    histo.Fill(2,x0,xst[i]);
    //histo.Fill(3,x0,em[i]);
    if(idx == 6) histo.Fill(5,lowx0,lowxsi[i]);
  }
  G4cout << "------------------------------------------------------------------------------------"
	 << G4endl;
  G4cout.precision(prec);
  
  if(verbose > 0) { G4cout << "###### Save histograms" << G4endl; }
  histo.Save();
  
  if(verbose > 0) {
    G4cout << "###### End of run # " << G4endl;
  }
  G4double ma = G4PionPlus::PionPlus()->GetPDGMass();
  G4cout << std::sqrt(ma*ma  + 10000) - ma << G4endl;
}
