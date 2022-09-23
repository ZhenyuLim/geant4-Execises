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
/// \file PrimaryGeneratorAction.hh
/// \brief Definition of the B2::PrimaryGeneratorAction class

#ifndef B2PrimaryGeneratorAction_h
#define B2PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleGun.hh"

class G4ParticleGun;
class G4Event;



/// The primary generator action class with particle gum.
///
/// It defines a single particle which hits the Tracker
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class
/// (see the macros provided with this example).

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* ) override;

    G4ParticleGun* GetParticleGun() {return particleGun;}

    //reture source position
     G4double sourcepositionX (){ return sourceX;};

     G4double sourcepositionY (){ return sourceY;};

     G4double sourcepositionZ (){ return sourceZ;};

  private:
    G4ParticleGun* particleGun = nullptr; // G4 particle gun

    //source position

     DetectorConstruction* getsizeofroom = new DetectorConstruction;

     G4double roomX = getsizeofroom->getroominLength();

     G4double roomY = getsizeofroom->getroominWidth();

     G4double roomZ = getsizeofroom->getroominHighth();

     G4double sourceX = 0.*m;

     G4double sourceY = roomY/2-0.9*m;

     G4double sourceZ = -roomZ/2+1.6*m;

     G4double rBeam;

     G4double pi = 3.1415926;

};



#endif
