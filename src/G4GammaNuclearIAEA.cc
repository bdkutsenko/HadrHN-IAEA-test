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
//
// GEANT4 Class file
//
//
// File name:    G4GammaNuclearIAEA
//
// Author  
//
// Modifications:
//

#include "G4GammaNuclearIAEA.hh"
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"
#include "G4ElementTable.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4PhotoNuclearCrossSection.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4Gamma.hh"
#include "G4IsotopeList.hh"

#include <fstream>
#include <sstream>
#include <vector>

// factory
#include "G4CrossSectionFactory.hh"
//


G4_DECLARE_XS_FACTORY(G4GammaNuclearIAEA);

G4ElementData* G4GammaNuclearIAEA::data = nullptr;

G4double G4GammaNuclearIAEA::coeff[3][3];
G4String G4GammaNuclearIAEA::gDataDirectory = "";

#ifdef G4MULTITHREADED
  G4Mutex G4GammaNuclearIAEA::gNuclearIAEAMutex = G4MUTEX_INITIALIZER;
#endif

G4GammaNuclearIAEA::G4GammaNuclearIAEA() 
  : G4VCrossSectionDataSet(Default_Name()),
   gamma(G4Gamma::Gamma())
{
  //  verboseLevel = 0;
  if(verboseLevel > 0){
    G4cout  << "G4GammaNuclearIAEA::G4GammaNuclearIAEA Initialise for Z < " 
	    << MAXZGAMMAIAEA << G4endl;
  }
  ggXsection = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet("PhotoNuclearXS");
  if(ggXsection == nullptr) ggXsection = new G4PhotoNuclearCrossSection();
  theGamma = new G4DynamicParticle(gamma, G4ThreeVector(1,0,0), 0);
  SetForAllAtomsAndEnergies(true);
}

G4GammaNuclearIAEA::~G4GammaNuclearIAEA()
{
if(isMaster) { delete data; data = nullptr; }
}

void G4GammaNuclearIAEA::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4GammaNuclearIAEA calculates the gamma nuclear\n"
          << "cross-section on nuclei using data from the high precision\n"
          << "IAEA photonuclear database used on GDR energy region. Then liniear connection\n"
	  <<"implemented with previous CHIPS photonuclear model and connection \n with PDG parametrisation added at stil higher energies \n";
}

G4bool 
G4GammaNuclearIAEA::IsElementApplicable(const G4DynamicParticle*, 
                                      G4int, const G4Material*)
{
  return true;
}

G4bool G4GammaNuclearIAEA::IsIsoApplicable(const G4DynamicParticle*,
                                         G4int, G4int,
                                         const G4Element*, const G4Material*)
{
  return true;
}

G4double 
G4GammaNuclearIAEA::GetElementCrossSection(const G4DynamicParticle* aParticle,
                                         G4int ZZ, const G4Material*)
{
  return ElementCrossSection(aParticle->GetKineticEnergy(), ZZ);
}

G4double 
G4GammaNuclearIAEA::ElementCrossSection(G4double ekin, G4int ZZ){
    
  G4double xs = 0.0;
  G4double emax;
  G4int Z = (ZZ >= MAXZGAMMAIAEA) ? MAXZGAMMAIAEA - 1 : ZZ;
  theGamma->SetKineticEnergy(ekin);
  
  auto pv = GetPhysicsVector(Z);
  if(pv == nullptr) emax = 0.;
  else emax = pv->GetMaxEnergy();
  
  if(ekin <= emax) {
    xs = pv->Value(ekin*MeV);
  }
  else if(ekin <= rTransitionBound && emax !=0. ){
    theGamma->SetKineticEnergy(rTransitionBound);
    G4double rxs = ggXsection->GetElementCrossSection(theGamma, Z, 0);
    G4double lxs = pv->Value(emax);
    xs = lxs + (ekin - emax)*(rxs - lxs)/(rTransitionBound-emax);
  }
  else {
    xs = ggXsection->GetElementCrossSection(theGamma, Z, 0);
  }

#ifdef G4VERBOSE
  if(verboseLevel > 1) {
    G4cout  << "Z= " << Z << " Ekin(MeV)= " << ekin/CLHEP::MeV 
	    << ",  nElmIAEA(b)= " << xs/CLHEP::barn 
	    << G4endl;
  }
#endif
 
  return xs;

}

G4double G4GammaNuclearIAEA::GetIsoCrossSection(
         const G4DynamicParticle* aParticle, 
	 G4int ZZ, G4int A,
	 const G4Isotope*, const G4Element*,
	 const G4Material*)
{
  return IsoCrossSection(aParticle->GetKineticEnergy(), ZZ, A);
}

G4double 
G4GammaNuclearIAEA::IsoCrossSection(G4double ekin, G4int ZZ, G4int A)
{
  G4double xs = 0.0;
  G4double emax;
  G4int Z = (ZZ >= MAXZGAMMAIAEA) ? MAXZGAMMAIAEA - 1 : ZZ;
  
  /*
  G4cout << "IsoCrossSection  Z= " << Z << "  A= " << A 
         << "  Amin= " << amin[Z] << " Amax= " << amax[Z]
         << " E(MeV)= " << ekin << G4endl;
  */
  theGamma->SetKineticEnergy(ekin);
  auto pv = GetPhysicsVector(Z);
  if(pv == nullptr || Z == 1) emax = 0.; // Proton dont have Giant Dipole Resonance 
  else emax = pv->GetMaxEnergy();

  // compute isotope cross section if applicable
  if(amin[Z] > 0 && A >= amin[Z] && A <= amax[Z]){
    G4double emaxiso=0;
    auto pviso = data->GetComponentDataByIndex(Z, A - amin[Z]);
    if(pviso) emaxiso = pviso->GetMaxEnergy();
    if(pviso && ekin <= emaxiso ) {
      xs = pviso->Value(ekin*MeV);
  
#ifdef G4VERBOSE
    if(verboseLevel > 1) {
      G4cout  << "G4GammaNuclearIAEA::IsoIAEA: Z= " << Z << " A= " << A 
	    << " Ekin(MeV)= " << ekin/CLHEP::MeV 
	    << ", ElmXS(b)= " << xs/CLHEP::barn << G4endl;
  }
#endif
   return xs;
     }
   else if(pviso && ekin <= rTransitionBound) {
       theGamma->SetKineticEnergy(rTransitionBound);
       G4double rxs = ggXsection->GetIsoCrossSection(theGamma, Z, A);
       G4double lxs = pviso->Value(emaxiso);
       xs = lxs + (ekin - emaxiso)*(rxs - lxs)/(rTransitionBound-emaxiso);
   return xs;
     }   
  } 

 if(ekin <= emax) { 
    xs = pv->Value(ekin*MeV);
    xs *= A/aeff[Z];
  }
 else if(ekin <= rTransitionBound && emax !=0 ){
    theGamma->SetKineticEnergy(rTransitionBound);
    G4double rxs = ggXsection->GetIsoCrossSection(theGamma, Z, 0);
    G4double lxs = pv->Value(emax);
    xs = lxs + (ekin - emax)*(rxs - lxs)/(rTransitionBound-emax);
  }
 else {
   if(Z<=2 && ekin>10.*GeV) xs = coeff[Z][A - amin[Z]]*ggXsection->GetElementCrossSection(theGamma, Z, 0);
   else xs = ggXsection->GetIsoCrossSection(theGamma, Z, A);
 }

  return xs;

}


const G4Isotope* G4GammaNuclearIAEA::SelectIsotope(
       const G4Element* anElement, G4double kinEnergy, G4double)
{
  size_t nIso = anElement->GetNumberOfIsotopes();
  const G4Isotope* iso = anElement->GetIsotope(0);

  if(1 == nIso) { return iso; }

  const G4double* abundVector = anElement->GetRelativeAbundanceVector();
  G4double q = G4UniformRand();
  G4double sum = 0.0;
  size_t j;
  G4int Z = anElement->GetZasInt();

  if(0 == amin[Z] || Z >= MAXZGAMMAIAEA) {
  for (j=0; j<nIso; ++j) {
    sum += abundVector[j];
    if(q <= sum) {
      iso = anElement->GetIsotope(j);
      break;
    }
  }
  return iso;
  }
  // use isotope cross sections
  size_t nn = temp.size();
  if(nn < nIso) { temp.resize(nIso, 0.); }
  
  for (j=0; j<nIso; ++j) {
    //G4cout << j << "-th isotope " << (*isoVector)[j]->GetN() 
    //       <<  " abund= " << abundVector[j] << G4endl;
    sum += abundVector[j]*IsoCrossSection(kinEnergy, Z, anElement->GetIsotope(j)->GetN());
    temp[j] = sum;
  }
  sum *= q;
  for (j = 0; j<nIso; ++j) {
    if(temp[j] >= sum) {
      iso = anElement->GetIsotope(j);
      break;
    }
  }
  return iso;
}

void 
G4GammaNuclearIAEA::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if(verboseLevel > 0){
    G4cout << "G4GammaNuclearIAEA::BuildPhysicsTable for " 
	   << p.GetParticleName() << G4endl;
  }
  if(p.GetParticleName() != "gamma") { 
    G4ExceptionDescription ed;
    ed << p.GetParticleName() << " is a wrong particle type -"
       << " only gamma is allowed";
    G4Exception("G4GammaNuclearIAEA::BuildPhysicsTable(..)","had012",
		FatalException, ed, "");
    return; 
  }

    if(nullptr == data) { 
#ifdef G4MULTITHREADED
    G4MUTEXLOCK(&neutronInelasticXSMutex);
    if(nullptr == data) { 
#endif
      isMaster = true;
      data = new G4ElementData(); 
      data->SetName("PhotoNuclear");
      FindDirectoryPath();
#ifdef G4MULTITHREADED
    }
    G4MUTEXUNLOCK(&neutronInelasticXSMutex);
#endif
  }
    

  // it is possible re-initialisation for the second run
  // Upload data for elements used in geometry
  const G4ElementTable* table = G4Element::GetElementTable();
  if(isMaster) {
    for ( auto & elm : *table ) {
      G4int Z = std::max( 1, std::min( elm->GetZasInt(), MAXZGAMMAIAEA-1) );
      if ( nullptr == data->GetElementData(Z) ) { Initialise(Z); } 
    }
  }

    // prepare isotope selection
  size_t nIso = temp.size();
  for ( auto & elm : *table ) {
    size_t n = elm->GetNumberOfIsotopes();
    if(n > nIso) { nIso = n; }
  }
  temp.resize(nIso, 0.0);

   
}

const G4String& G4GammaNuclearIAEA::FindDirectoryPath()
{
  // check environment variable
  // build the complete string identifying the file with the data set
  if(gDataDirectory.empty()) {
    char* path = std::getenv("G4PARTICLEXSDATA");
    if (path) {
      std::ostringstream ost;
      ost << path << "/gamma/inel";
      gDataDirectory = ost.str();
    } else {
      G4Exception("G4GammaNuclearIAEA::Initialise(..)","had013",
		  FatalException,
		  "Environment variable G4PARTICLEIAEADATA is not defined");
    }
  }
  return gDataDirectory;
}

void G4GammaNuclearIAEA::InitialiseOnFly(G4int Z)
{
#ifdef G4MULTITHREADED
   G4MUTEXLOCK(&neutronInelasticXSMutex);
   if(nullptr == data->GetElementData(Z)) { 
#endif
     Initialise(Z);
#ifdef G4MULTITHREADED
   }
   G4MUTEXUNLOCK(&neutronInelasticXSMutex);
#endif
}

void G4GammaNuclearIAEA::Initialise(G4int Z)
{
  
  if(nullptr != data->GetElementData(Z)) { return; }

  // upload data from file
 
  std::ostringstream ost;
  //ost << FindDirectoryPath() << Z ;
  ost << "data/inel" << Z ;
  G4PhysicsVector* v = RetrieveVector(ost, true, Z);
  data->InitialiseForElement(Z, v);
 /*
  G4cout << "G4NeutronInelasticXS::Initialise for Z= " << Z 
	 << " A= " << Amean << "  Amin= " << amin[Z] 
	 << "  Amax= " << amax[Z] << G4endl;
  */
  // upload isotope data
  if(amin[Z] > 0) {
    size_t nmax = (size_t)(amax[Z]-amin[Z]+1);
    data->InitialiseForComponent(Z, nmax);
    for(G4int A=amin[Z]; A<=amax[Z]; ++A) {
      std::ostringstream ost1;
      //ost1 << gDataDirectory << Z << "_" << A;
      ost1 << "data/inel"<< Z << "_" << A;
      G4PhysicsVector* v1 = RetrieveVector(ost1, false, Z);
      data->AddComponent(Z, A, v1);
      if(Z<=2){
	theGamma->SetKineticEnergy(10.*GeV);
	G4double sig1 = ggXsection->GetIsoCrossSection(theGamma, Z, A);
	G4double sig2 = ggXsection->GetElementCrossSection(theGamma, Z, 0);
	//G4cout<< " sig1 = "<<sig1<<" "<<sig1/sig2<<G4endl;
	if(sig2 > 0.) coeff[Z][A-amin[Z]]=(sig1/sig2);
	else coeff[Z][A-amin[Z]]=1.;
      }
    }   
  }  
}


G4PhysicsVector* 
G4GammaNuclearIAEA::RetrieveVector(std::ostringstream& ost, G4bool warn, G4int Z)
{
  G4PhysicsVector* v = nullptr;
  
  std::ifstream filein(ost.str().c_str());
  if (!(filein)) {
    if(warn) {
      G4ExceptionDescription ed;
      ed << "Data file <" << ost.str().c_str()
	 << "> is not opened!";
      G4Exception("G4NeutronInelasticXS::RetrieveVector(..)","had014",
		  FatalException, ed, "Check G4PARTICLEXSDATA");
    }
  } else {
    if(verboseLevel > 1) {
      G4cout << "File " << ost.str() 
	     << " is opened by G4GammaNuclearIAEA" << G4endl;
    }
    // retrieve data from DB
    if(std::find(std::begin(freeVectorException), std::end(freeVectorException), Z ) == std::end(freeVectorException) && warn)  v = new G4PhysicsLinearVector();
    else v = new G4PhysicsVector();
    
    if(!v->Retrieve(filein, true)) {
      G4ExceptionDescription ed;
      ed << "Data file <" << ost.str().c_str()
	 << "> is not retrieved!";
      G4Exception("G4NeutronInelasticXS::RetrieveVector(..)","had015",
		  FatalException, ed, "Check G4PARTICLEXSDATA");
    }

  }
   
 return v;

}

