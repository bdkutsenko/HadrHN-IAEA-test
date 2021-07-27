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
#include "G4Gamma.hh"
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"
#include "G4ElementTable.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Gamma.hh"
#include "G4PhysicsLogVector.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4PhotoNuclearCrossSection.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
//#include "G4IsotopeList.hh"

#include <fstream>
#include <sstream>
#include <vector>

// factory
#include "G4CrossSectionFactory.hh"
//
G4_DECLARE_XS_FACTORY(G4GammaNuclearIAEA);

G4ElementData* G4GammaNuclearIAEA::data = nullptr;

G4double G4GammaNuclearIAEA::coeff[] = {0.0};
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
  SetForAllAtomsAndEnergies(true);
}

G4GammaNuclearIAEA::~G4GammaNuclearIAEA()
{
if(isMaster) { delete data; data = nullptr; }
}

void G4GammaNuclearIAEA::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4GammaNuclearIAEA calculates the gamma nuclear\n"
          << "cross section on nuclei using data from the high precision\n"
          << "IAEA gamma database.*bdk* INSERT DESCRIPTION LATER\n";
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
                                         G4int ZZ, const G4Material* mat)
{
  G4double xs = 0.0,rxs = 0.0, lxs = 0.0, lconnect = 130.*MeV, rconnect = 150.*MeV;
  G4double ekin = aParticle->GetKineticEnergy();
  G4DynamicParticle dParticle= *aParticle;
  G4int Z = (ZZ >= MAXZGAMMAIAEA) ? MAXZGAMMAIAEA - 1 : ZZ; 
  auto pv = GetPhysicsVector(Z);
  if(pv == nullptr) { return xs; }
  //if(ekin <= pv->GetMaxEnergy()) {
  if(ekin <= lconnect){
    xs = pv->Value(ekin*MeV);
    //xs = pv->LogVectorValue(ekin, aParticle->GetLogKineticEnergy());
  } else if(ekin > lconnect && ekin < rconnect){
    dParticle.SetKineticEnergy(rconnect);
    rxs = coeff[Z]*ggXsection->GetElementCrossSection(&dParticle, Z, mat);
    lxs = pv->Value(lconnect);
    xs = rxs + (rxs - lxs)/(rconnect-lconnect)*(ekin-rconnect); // Straight line equation xs(ekin) 
  }
  else {
    xs = coeff[Z]*ggXsection->GetElementCrossSection(aParticle, Z, mat);
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
//bdk// Add isotopes outside of IAEA
G4double G4GammaNuclearIAEA::GetIsoCrossSection(
         const G4DynamicParticle* aParticle, 
	 G4int ZZ, G4int A,
	 const G4Isotope*, const G4Element*,
	 const G4Material*)
{
  G4double xs = 0.0,rxs = 0.0, lxs = 0.0, lconnect = 130.*MeV, rconnect = 150.*MeV;
  G4int Z = (ZZ >= MAXZGAMMAIAEA) ? MAXZGAMMAIAEA - 1 : ZZ; 
  G4DynamicParticle dParticle= *aParticle;
  G4double ekin = aParticle->GetKineticEnergy();
  /*
  G4cout << "IsoCrossSection  Z= " << Z << "  A= " << A 
         << "  Amin= " << amin[Z] << " Amax= " << amax[Z]
         << " E(MeV)= " << ekin << G4endl;
  */
  //auto pv = GetPhysicsVector(Z);
  // if(pv == nullptr) { return xs; }
  
  // compute isotope cross section if applicable
  if(amin[Z] > 0 && A >= amin[Z] && A <= amax[Z]){
    auto pviso = data->GetComponentDataByIndex(Z, A - amin[Z]);
     if(pviso) {
        if(ekin <= lconnect){
    xs = pviso->Value(ekin*MeV);
    //xs = pv->LogVectorValue(ekin, aParticle->GetLogKineticEnergy());
  } else if(ekin > lconnect && ekin < rconnect){
    dParticle.SetKineticEnergy(rconnect);
    rxs = coeff[Z]*ggXsection->GetIsoCrossSection(&dParticle,Z,A);
    lxs = pviso->Value(lconnect);
    xs = rxs + (rxs - lxs)/(rconnect-lconnect)*(ekin-rconnect); // Straight line equation xs(ekin) 
  }
  else {
    xs = coeff[Z]*ggXsection->GetIsoCrossSection(&dParticle,Z,A);
  }
#ifdef G4VERBOSE
  if(verboseLevel > 1) {
    G4cout  << "G4GammaNuclearIAEA::IsoIAEA: Z= " << Z << " A= " << A 
	    << " Ekin(MeV)= " << ekin/CLHEP::MeV 
	    << ", ElmXS(b)= " << xs/CLHEP::barn << G4endl;
  }
#endif
   return xs;
     }
       /*//bdk// TODO Isotope cross section based on element cross-section
   auto pv = GetPhysicsVector(Z);
    if(pv) {
     if(ekin <= lconnect){
    xs = pviso->Value(ekin*MeV);
    //xs = pv->LogVectorValue(ekin, aParticle->GetLogKineticEnergy());
  } else if(ekin > lconnect && ekin < rconnect){
    dParticle.SetKineticEnergy(rconnect);
    rxs = coeff[Z]*ggXsection->GetIsoCrossSection(&dParticle,Z,A);
    lxs = pviso->Value(lconnect);
    xs = rxs + (rxs - lxs)/(rconnect-lconnect)*(ekin-rconnect); // Straight line equation xs(ekin) 
  }
  else {
    xs = coeff[Z]*ggXsection->GetIsoCrossSection(&dParticle,Z,A);
  }
  xs *= A/aeff[Z];   
#ifdef G4VERBOSE
  if(verboseLevel > 1) {
    G4cout  << "G4GammaNuclearIAEA::IsoIAEA: Z= " << Z << " A= " << A 
	    << " Ekin(MeV)= " << ekin/CLHEP::MeV 
	    << ", ElmXS(b)= " << xs/CLHEP::barn << G4endl;
  }
#endif
   return xs;
     }
       */
       
  } 


  xs = coeff[Z]*ggXsection->GetIsoCrossSection(&dParticle,Z,A);
  return xs;
  
  //return IsoCrossSection(aParticle->GetKineticEnergy(), 
  //    aParticle->GetLogKineticEnergy(), Z, A);
}


const G4Isotope* G4GammaNuclearIAEA::SelectIsotope(
       const G4Element* anElement, G4double kinEnergy, G4double)
{
  size_t nIso = anElement->GetNumberOfIsotopes();
  const G4Isotope* iso = anElement->GetIsotope(0);

  //G4cout << "SelectIsotope NIso= " << nIso << G4endl;
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
  
  G4ThreeVector aDirection = G4ThreeVector(0.0,0.0,1.0);
  G4DynamicParticle dParticle(G4Gamma::Gamma(),aDirection,kinEnergy);
  
  for (j=0; j<nIso; ++j) {
    //G4cout << j << "-th isotope " << (*isoVector)[j]->GetN() 
    //       <<  " abund= " << abundVector[j] << G4endl;
    sum += abundVector[j]*GetIsoCrossSection(&dParticle,  Z, anElement->GetIsotope(j)->GetN(),0,0,0);
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
    
  /*bdk//Coeff change//
  if(0. == coeff[0]) { 
#ifdef G4MULTITHREADED
    G4MUTEXLOCK(&gNuclearIAEAMutex);
    if(0. == coeff[0]) { 
#endif
      coeff[0] = 1.0;
      isMaster = true;
      FindDirectoryPath();
#ifdef G4MULTITHREADED
    }
    G4MUTEXUNLOCK(&gNuclearIAEAMutex);
#endif
  }
  */
  

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
  
  G4PhysicsVector* v = RetrieveVector(ost, true);
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
      ost1 << gDataDirectory << Z << "_" << A;
      G4PhysicsVector* v1 = RetrieveVector(ost1, false);
      data->AddComponent(Z, A, v1); 
    }
  }

  // smooth transition 
  G4double sig1 = (*v)[v->GetVectorLength()-1];
  //G4double ehigh= v->GetMaxEnergy();
  G4double sig2 = 0;//ggXsection->GetInelasticElementCrossSection(neutron, ehigh, Z, aeff[Z]);
  coeff[Z] = (sig2 > 0.) ? sig1/sig2 : 1.0; 
}


G4PhysicsVector* 
G4GammaNuclearIAEA::RetrieveVector(std::ostringstream& ost, G4bool warn)
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
	     << " is opened by G4NeutronInelasticXS" << G4endl;
    }
    // retrieve data from DB
    v = new G4PhysicsVector();
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





  
/*
void G4GammaNuclearIAEA::Initialise(G4int Z, G4int A)
{
  if(data[Z] != nullptr) { return; }
  for(int id=0; id<MAXZGAMMAIAEA; id++) {
    G4PhysicsLogVector* isodata[MAXZGAMMAIAEA]={nullptr};
    vecisodata.push_back(isodata[0]);
  }
  
  // upload data from file
  vecisodata[Z][A] = new G4PhysicsLogVector();
  
  std::ostringstream ost;
  ost << FindDirectoryPath() << Z ;
  std::ifstream filein(ost.str().c_str());
  if (!(filein)) {
    G4ExceptionDescription ed;
    ed << "Data file <" << ost.str().c_str()
       << "> is not opened!";
    G4Exception("G4GammaNuclearIAEA::Initialise(..)","had014",
                FatalException, ed, "Check G4PARTICLEIAEADATA");
    return;
  }
  if(verboseLevel > 1) {
    G4cout << "file " << ost.str() 
	   << " is opened by G4GammaNuclearIAEA" << G4endl;
  }
    
  // retrieve data from DB
  if(!data[Z]->Retrieve(filein, true)) {
    G4ExceptionDescription ed;
    ed << "Data file <" << ost.str().c_str()
       << "> is not retrieved!";
    G4Exception("G4GammaNuclearIAEA::Initialise(..)","had015",
		FatalException, ed, "Check G4PARTICLEIAEADATA");
    return;
  }
  // smooth transition 
}
*/
