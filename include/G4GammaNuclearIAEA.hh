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
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:    G4GammaNuclearIAEA
//
// Author  B.Kutsenko, V.Ivantchenko, Geant4, 1 August 2021
//

 
// Class Description:
// This is a base class for gamma-nuclear cross section based on
// data files from IAEA Evaluated Photonuclear Data Library (IAEA/PD-2019) https://www-nds.iaea.org/photonuclear/
// Database review - https://www.sciencedirect.com/science/article/pii/S0090375219300699?via%3Dihub
// Class Description - End

#ifndef G4GammaNuclearIAEA_h
#define G4GammaNuclearIAEA_h 1

#include "G4VCrossSectionDataSet.hh"
#include "globals.hh"
#include "G4ElementData.hh"
#include "G4PhysicsVector.hh"
#include "G4Threading.hh"
//#include "G4IsotopeList.hh"
#include <vector>

/*
static const G4int amin[] = {
  0,
  1,   3,   6,   9,  10,  12,  14,  16,  19,  20,  //1-10
 23,  24,  27,  28,  31,  32,  35,  36,  39,  40,  //11-20
 45,  46,  50,  50,  55,  54,  59,  58,  63,  64,  //21-30
 69,  70,  75,   74,   79,  78,   85,   84,  89,  90,  //31-40
  93,  92,   98,   96,   103, 102, 107, 106, 113, 112,  //41-50
  121,   120, 127,   124, 133, 130,   138,   136,   141,   142,  //51-60
  145,   144,   151, 152,   159,   156,   165,   162,   169, 168,  //61-70
  175,   174, 180, 180,   185,   184,   191, 190, 197,   196,  //71-80
203, 204, 209,   0,   0,   0,   0,   226,   0,   232,  //81-90
  231, 234, 237, 238 };

static const G4int amax[] = {
  0,
  3,   4,   7,   9,  11,  13,  15,  18,  19,  22,  //1-10
 23,  26,  27,  30,  31,  36,  37,  40,  41,  48,  //11-20
 45,  50,  51,  54,  55,  58,  59,  64,  65,  70,  //21-30
 71,  76,  75,   82,   81,  86,   87,   88,  89,  96,  //31-40
  93, 100,   98,   104,   103, 110, 109, 116, 115, 124,  //41-50
  123,   130, 127,   136, 133, 138,   139,   142,   141,   150,  //51-60
  145,   154,   153, 160,   159,   164,   165,   170,   169, 176,  //61-70
  176,   180, 181, 186,   187,   192,   193, 198, 197,   204,  //71-80
205, 208, 209,   0,   0,   0,   0,   226,   0,   232,  //81-90
  231, 238, 237, 241 };

static const G4double aeff[] = {
  0.,
 1.00794, 4.00264, 6.94003, 9.01218,  10.811, 12.0107, 14.0068, 15.9994, 18.9984,   20.18,  //1-10
 22.9898,  24.305, 26.9815, 28.0854, 30.9738, 32.0661, 35.4526, 39.9477, 39.0983,  40.078,  //11-20
 44.9559, 47.8667, 50.9415, 51.9961,  54.938, 55.8451, 58.9332, 58.6933, 63.5456, 65.3955,  //21-30
 69.7231, 72.6128, 74.9216, 78.9594, 79.9035, 83.7993, 85.4677, 87.6166, 88.9058, 91.2236,  //31-40
 92.9064, 95.9313, 97.9072, 101.065, 102.906, 106.415, 107.868, 112.411, 114.818,  118.71,  //41-50
  121.76, 127.603, 126.904, 131.292, 132.905, 137.327, 138.905, 140.115, 140.908, 144.236,  //51-60
 144.913, 150.366, 151.964, 157.252, 158.925, 162.497,  164.93, 167.256, 168.934, 173.038,  //61-70
 174.967, 178.485, 180.948, 183.842, 186.207, 190.225, 192.216, 195.078, 196.967, 200.599,  //71-80
 204.383, 207.217,  208.98, 208.982, 209.987, 222.018,  223.02, 226.025, 227.028, 232.038,  //81-90
 231.036, 238.029};
*/

class G4DynamicParticle;
class G4ParticleDefinition;
class G4Element;
class G4VComponentCrossSection;

class G4GammaNuclearIAEA final : public G4VCrossSectionDataSet
{
public: 

  explicit G4GammaNuclearIAEA();

  ~G4GammaNuclearIAEA() final;
    
  static const char* Default_Name() {return "G4GammaNuclearIAEA";}

  G4bool IsElementApplicable(const G4DynamicParticle*, 
			     G4int Z, const G4Material*) final;

  G4bool IsIsoApplicable(const G4DynamicParticle*, G4int Z, G4int A,
			 const G4Element*, const G4Material* mat) final;

  G4double GetElementCrossSection(const G4DynamicParticle*, 
			          G4int Z, const G4Material* mat) final; 

  G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A,
                              const G4Isotope* iso,
                              const G4Element* elm,
                              const G4Material* mat) final;

  const G4Isotope* SelectIsotope(const G4Element*, 
                                 G4double kinEnergy, G4double logE) final;

  void BuildPhysicsTable(const G4ParticleDefinition&) final;

  G4double IsoCrossSection(G4double ekin, G4int Z, G4int A);

  G4double ElementCrossSection(G4double ekin, G4int Z);

  void CrossSectionDescription(std::ostream&) const final;
      
  G4GammaNuclearIAEA & operator=(const G4GammaNuclearIAEA &right) = delete;
  G4GammaNuclearIAEA(const G4GammaNuclearIAEA&) = delete;
  
private: 

  void Initialise(G4int Z);

  void InitialiseOnFly(G4int Z);

  const G4String& FindDirectoryPath();

  inline G4PhysicsVector* GetPhysicsVector(G4int Z);

  G4PhysicsVector* RetrieveVector(std::ostringstream& in, G4bool warn, G4int Z);
  
  G4VCrossSectionDataSet* ggXsection = nullptr;

  G4DynamicParticle* theGamma = nullptr;
  std::vector<G4double> temp;
  const G4ParticleDefinition* gamma;
  
  G4bool isMaster = false;

  static const G4int MAXZGAMMAIAEA = 94;
  static const G4int rTransitionBound = 150.*CLHEP::MeV; // Upper limit of the linear transition between IAEA database and CHIPS model
  const G4int freeVectorException[11] = {4, 6, 7, 8, 27, 39, 45, 65, 67, 69, 73};
  static G4ElementData* data;
  static G4double coeff[3][3]; // Previous CHIPS photonuclear model had a problem with high energy parametrisation for isotopes of H and He, coefficient is needed to connect isotope cross section with element cross-section on high energies.
  static G4String gDataDirectory;
  
#ifdef G4MULTITHREADED
  static G4Mutex gNuclearIAEAMutex;
#endif
};

inline
G4PhysicsVector* G4GammaNuclearIAEA::GetPhysicsVector(G4int Z)
{
  G4PhysicsVector* pv = data->GetElementData(Z);
  if(pv == nullptr) { 
    InitialiseOnFly(Z);
    pv = data->GetElementData(Z);
  }
  return pv;
}

#endif
