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
/// \file DetectorConstruction.hh
/// \brief Definition of the B2a::DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"



class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4GlobalMagFieldMessenger;





class DetectorMessenger;

/// Detector construction class to define materials, geometry
/// and global uniform magnetic field.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction() override;

  public:
    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    // Set methods
    void SetWorldMaterial (G4String );
    void SetRoomMaterial (G4String );    
    void SetCollimaterMaterial(G4String );


    
    G4double getroominLength(){ return  roominLength;};

     G4double getroominWidth(){return roominWidth;};

     G4double getroominHighth(){return roominHighth;};

    G4LogicalVolume* GetLogicDetector1(); 

    G4LogicalVolume* GetLogicDetector2(); 

  private:
    // methods
    void DefineMaterials();

    G4VPhysicalVolume* DefineVolumes();
    
    G4LogicalVolume* roomdetLV = nullptr;

    G4LogicalVolume* absorberdetLV = nullptr;
    
    G4Material* Air = nullptr; 
    G4Material* Pb = nullptr; 
    G4Material* Concrete = nullptr;
    G4Material* Fe = nullptr; 
    G4double roominLength = 4.5*m ; 

    G4double roominWidth  = 6*m ; 

    G4double roominHighth = 4.5*m; 

};


#endif
