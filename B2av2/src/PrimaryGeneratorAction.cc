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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the B2::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  G4int nofParticles = 1;
  particleGun = new G4ParticleGun(nofParticles);

  // default particle kinematic

  G4ParticleDefinition* particleDefinition
    = G4ParticleTable::GetParticleTable()->FindParticle("gamma");

  particleGun->SetParticleDefinition(particleDefinition);

  particleGun->SetParticleEnergy(662*keV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

///方式1 

 G4double conTheta = 2*G4UniformRand()-1., phi = 1.5*pi+(2*G4UniformRand()-1)*0.2*pi;//沿着长边方向30度立体角
G4double sinTheta = std::sqrt(1.-conTheta*conTheta);
G4double ux = sinTheta*std::cos(phi), uy = sinTheta*std::sin(phi), uz = conTheta;
particleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz));  
    
/* ///方式4pi立体角
auto theta = (G4UniformRand())*pi;
auto phi = (G4UniformRand())*2*pi;
G4double ux= std::sin(theta)*std::cos(phi);
G4double uy= std::sin(theta)*std::sin(phi);
G4double uz= std::cos(theta);
particleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz)); */

    G4double x0 = sourceX;
    G4double y0 = sourceY;
    G4double z0 = sourceZ;
    
particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
    
    


    

    particleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



