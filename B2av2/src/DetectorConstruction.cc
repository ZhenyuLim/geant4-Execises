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
/// \file DetectorConstruction.cc
/// \brief Implementation of the B2a::DetectorConstruction class

#include "DetectorConstruction.hh"
/* #include "DetectorMessenger.hh" */
#include "TrackerSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"

#include "G4SubtractionSolid.hh"

#include "PrimaryGeneratorAction.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


DetectorConstruction::DetectorConstruction()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Material definition

  G4NistManager* nistManager = G4NistManager::Instance();

  // Air defined using NIST Manager
  Air = nistManager->FindOrBuildMaterial("G4_AIR");

  // Lead defined using NIST Manager
  Pb  = nistManager->FindOrBuildMaterial("G4_Pb");

  // Xenon gas defined using NIST Manager
  Concrete = nistManager->FindOrBuildMaterial("G4_CONCRETE");

  // Print materials
/*   G4cout << *(G4Material::GetMaterialTable()) << G4endl; */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
// world volume
  G4double worldLength = 7*m; // from chamber center to center!

  G4double worldWidth  = 10*m; // width of the chambers

  G4double worldHighth = 7*m; // width of the chambers

  G4Box* worldS = new G4Box("world",worldLength/2,worldWidth/2,worldHighth/2); 
  
  G4LogicalVolume* worldLV = new G4LogicalVolume(worldS,Air,"World");

  G4VPhysicalVolume* worldPV = new G4PVPlacement(0,G4ThreeVector(),worldLV,"World",0,true,0,true);

// room volume


  G4double roomoutLength = roominLength/2+20*cm; 

  G4double roomoutWidth  = roominWidth/2+20*cm; 

  G4double roomoutHighth = roominHighth/2+20*cm; 


  G4Box* innerroomS = new G4Box("innerroom",roominLength/2,roominWidth/2,roominHighth/2); 
  
  G4Box* outroomS = new G4Box("outroom",roomoutLength,roomoutWidth,roomoutHighth); 

  G4SubtractionSolid* roomsolid = new G4SubtractionSolid("Room",outroomS,innerroomS);

  G4LogicalVolume* roomLV = new G4LogicalVolume(roomsolid,Concrete,"roomLV");

  new G4PVPlacement(0,G4ThreeVector(),roomLV,"roomPV",worldLV,true,0,true);
  
  ////////////////////////////////get source position

  PrimaryGeneratorAction sourceposition;

  G4double sourceX = sourceposition.sourcepositionX();


  G4double sourceY = sourceposition.sourcepositionY();


  G4double sourceZ = sourceposition.sourcepositionZ();
 
  //construct source container
   
  G4double ScinLength = 2.5*cm; 

  G4double ScinWidth  = 2.5*cm; 

  G4double ScinHighth = 2.5*cm; 

  G4double ScoutLength = ScinLength+4*cm; 

  G4double ScoutWidth = ScinWidth+4*cm;

  G4double ScoutHighth = ScinHighth+4*cm;
  
    //Air Box
  G4Box* innerSc = new G4Box("innerSc",ScinLength,ScinWidth,ScinHighth); 
  
  G4Box* outSc = new G4Box("outSc",ScoutLength,ScoutWidth,ScoutHighth); 

  G4SubtractionSolid* Scsolid = new G4SubtractionSolid("SourceContain",outSc,innerSc);
  
    // window 2 --Tub
  G4double windowRadius = 2.5*cm; 

  G4double windwWidth = ScoutWidth-ScinWidth+0.1*cm;

  G4Tubs* windows = new G4Tubs("Outcoll_soild", 0, windowRadius, windwWidth, 0.*deg, 360.*deg);

  G4RotationMatrix* yrot = new G4RotationMatrix();
  
  yrot->rotateX(G4double(90*deg));

  G4ThreeVector yTrans(0,-ScinWidth-windwWidth,0);

  G4SubtractionSolid* subtraction = new G4SubtractionSolid("Box-window", Scsolid, windows, yrot, yTrans);

  G4LogicalVolume* ScLV = new G4LogicalVolume(subtraction,Pb,"SourceContainLV");

  new G4PVPlacement(0,G4ThreeVector(sourceX,sourceY,sourceZ),ScLV,"SourceContainPV",worldLV,true,0,true);



// collimater volume - type1
 //layer1
  G4double outcollimaterRadius1 = windowRadius+1*cm; 

  G4double outcollimaterHighth1 = 7.5*mm; 

  G4double inncollimaterRadius1 = windowRadius;

  G4double incollimaterHighth1 = 7.6*mm;  

  G4Tubs* outcoll1soild1  = new G4Tubs("Outcoll_soild", 0, outcollimaterRadius1, outcollimaterHighth1, 0.*deg, 360.*deg);

  G4Tubs* inncoll1soild1  = new G4Tubs("Innercoll_soild", 0, inncollimaterRadius1, incollimaterHighth1, 0.*deg, 360.*deg);

  G4SubtractionSolid* collimatersolid1 = new G4SubtractionSolid("Collimater",outcoll1soild1,inncoll1soild1);

  G4LogicalVolume* collimaterLV1 = new G4LogicalVolume(collimatersolid1,Pb,"Collimater1");

//rotate setting
  G4RotationMatrix* rot = new G4RotationMatrix();
  rot->rotateX(G4double(90*deg));

  new G4PVPlacement(rot,G4ThreeVector(0*m,sourceY-ScoutWidth-outcollimaterHighth1-0.5*mm,sourceZ),collimaterLV1,"Collimater1",worldLV,true,0,true);


 //layer2

  G4double inncollimaterRadius2 = outcollimaterRadius1;

  G4double incollimaterHighth2 = 7.6*mm;  

  G4double outcollimaterRadius2 = inncollimaterRadius2+1*cm; 

  G4double outcollimaterHighth2 = 7.5*mm; 

  G4Tubs* outcoll1soild2  = new G4Tubs("Outcoll_soild", 0, outcollimaterRadius2, outcollimaterHighth2, 0.*deg, 360.*deg);

  G4Tubs* inncoll1soild2  = new G4Tubs("Innercoll_soild", 0, inncollimaterRadius2, incollimaterHighth2, 0.*deg, 360.*deg);

  G4SubtractionSolid* collimatersolid2 = new G4SubtractionSolid("Collimater",outcoll1soild2,inncoll1soild2);

  G4LogicalVolume* collimaterLV2 = new G4LogicalVolume(collimatersolid2,Pb,"Collimater2");

  new G4PVPlacement(rot,G4ThreeVector(0*m,sourceY-ScoutWidth-3*outcollimaterHighth1-0.5*mm,sourceZ),collimaterLV2,"Collimater2",worldLV,true,0,true);




 //layer3

  G4double inncollimaterRadius3 = outcollimaterRadius2;

  G4double incollimaterHighth3 = 7.6*mm;  

  G4double outcollimaterRadius3 = inncollimaterRadius3+1*cm; 

  G4double outcollimaterHighth3 = 7.5*mm; 

  G4Tubs* outcoll1soild3  = new G4Tubs("Outcoll_soild", 0, outcollimaterRadius3, outcollimaterHighth3, 0.*deg, 360.*deg);

  G4Tubs* inncoll1soild3  = new G4Tubs("Innercoll_soild", 0, inncollimaterRadius3, incollimaterHighth3, 0.*deg, 360.*deg);

  G4SubtractionSolid* collimatersolid3 = new G4SubtractionSolid("Collimater",outcoll1soild3,inncoll1soild3);

  G4LogicalVolume* collimaterLV3 = new G4LogicalVolume(collimatersolid3,Pb,"Collimater3");

  new G4PVPlacement(rot,G4ThreeVector(0*m,sourceY-ScoutWidth-5*outcollimaterHighth1-0.5*mm,sourceZ),collimaterLV3,"Collimater3",worldLV,true,0,true);




 //layer4

  G4double inncollimaterRadius4 = outcollimaterRadius3;

  G4double incollimaterHighth4 = 7.6*mm;  

  G4double outcollimaterRadius4 = inncollimaterRadius4+1*cm; 

  G4double outcollimaterHighth4 = 7.5*mm; 

  G4Tubs* outcoll1soild4  = new G4Tubs("Outcoll_soild", 0, outcollimaterRadius4, outcollimaterHighth4, 0.*deg, 360.*deg);

  G4Tubs* inncoll1soild4  = new G4Tubs("Innercoll_soild", 0, inncollimaterRadius4, incollimaterHighth4, 0.*deg, 360.*deg);

  G4SubtractionSolid* collimatersolid4 = new G4SubtractionSolid("Collimater",outcoll1soild4,inncoll1soild4);

  G4LogicalVolume* collimaterLV4 = new G4LogicalVolume(collimatersolid4,Pb,"Collimater4");

  new G4PVPlacement(rot,G4ThreeVector(0*m,sourceY-ScoutWidth-7*outcollimaterHighth1-0.5*mm,sourceZ),collimaterLV4,"Collimater4",worldLV,true,0,true);



 //layer5

  G4double inncollimaterRadius5 = outcollimaterRadius4;

  G4double incollimaterHighth5 = 7.6*mm;  

  G4double outcollimaterRadius5 = inncollimaterRadius5+1*cm; 

  G4double outcollimaterHighth5 = 7.5*mm; 

  G4Tubs* outcoll1soild5  = new G4Tubs("Outcoll_soild", 0, outcollimaterRadius5, outcollimaterHighth5, 0.*deg, 360.*deg);

  G4Tubs* inncoll1soild5  = new G4Tubs("Innercoll_soild", 0, inncollimaterRadius5, incollimaterHighth5, 0.*deg, 360.*deg);

  G4SubtractionSolid* collimatersolid5 = new G4SubtractionSolid("Collimater",outcoll1soild5,inncoll1soild5);

  G4LogicalVolume* collimaterLV5 = new G4LogicalVolume(collimatersolid5,Pb,"Collimater5");

  new G4PVPlacement(rot,G4ThreeVector(0*m,sourceY-ScoutWidth-9*outcollimaterHighth1-0.5*mm,sourceZ),collimaterLV5,"Collimater5",worldLV,true,0,true);

 //detector volume 1

  G4double detecLength = roomoutLength+0.1*cm; 

  G4double detecWidth  = roomoutWidth+0.1*cm; 

  G4double detecHighth = roomoutHighth+0.1*cm; 

  G4Box* innerdetecS = new G4Box("innerdetec",detecLength,detecWidth,detecHighth);
  
  G4Box* outdetecS = new G4Box("outdetec",detecLength+2*cm,detecWidth+2*cm,detecHighth+2*cm); 

  G4SubtractionSolid* detecsolid = new G4SubtractionSolid("detecsd",outdetecS,innerdetecS);

  G4LogicalVolume* roomdetLV = new G4LogicalVolume(detecsolid,Concrete,"roomdetLV");

  new G4PVPlacement(0,G4ThreeVector(),roomdetLV,"detecpv",worldLV,true,0,true);  

  
/* detector volume 2   detection around source
  G4Box* outScdet = new G4Box("outScdet",ScinLength-0.1*cm,ScinWidth-0.1*cm,ScinHighth-0.1*cm); 

  G4Box* innerScdet = new G4Box("innerScdet",ScinLength-1.1*cm,ScinWidth-1.1*cm,ScinHighth-1.1*cm); 

  G4SubtractionSolid* Scdetsolid = new G4SubtractionSolid("SourceContain",outScdet,innerScdet);

  G4LogicalVolume* ScdetsolidLV = new G4LogicalVolume(Scdetsolid,Air,"ScdetsolidLV");

  new G4PVPlacement(0,G4ThreeVector(sourceX,sourceY,sourceZ),ScdetsolidLV,"ScdetsolidPV",worldLV,true,0,true);   */


// detector volume 3   A layer Pb absorber in the centre of room
/* G4Box* trackerS = new G4Box("tracker",50*cm,5*mm,50*cm);
G4LogicalVolume* trackerLV
    = new G4LogicalVolume(trackerS, Air, "Tracker",0,0,0);
  new G4PVPlacement(0,               // no rotation
                    G4ThreeVector(0,sourceY-1*m,sourceZ), // at (x,y,z)
                    trackerLV,       // its logical volume
                    "Tracker",       // its name
                    worldLV,         // its mother  volume
                    false,           // no boolean operations
                    0,               // copy number
                    true); // checking overlaps */
//////////////////



//////////////////

  G4Box* absorberdet = new G4Box("absorberdet",100*cm,4*cm,100*cm); 

  G4LogicalVolume* absorberdetLV = new G4LogicalVolume(absorberdet,Pb,"absorberdetLV");



  new G4PVPlacement(0,G4ThreeVector(0,sourceY-1*m,sourceZ),absorberdetLV,"absorberdetPV",worldLV,true,0,true); 







  return worldPV;
 


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  G4String Detect1SDname = "/Detect1SD";
  TrackerSD* aTrackerSD1 = new TrackerSD(Detect1SDname,
                                            "HitsCollection1");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD1);
  // Setting aTrackerSD to all logical volumes with the same name
  // of "Chamber_LV".
  SetSensitiveDetector("absorberdetLV", aTrackerSD1, true);

/////////////////////////////////////////////////

  G4String Detect2SDname = "/Detect2SD";
  TrackerSD* aTrackerSD2 = new TrackerSD(Detect2SDname,
                                            "HitsCollection2");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD2);
  // Setting aTrackerSD to all logical volumes with the same name
  // of "Chamber_LV".
  SetSensitiveDetector("roomdetLV", aTrackerSD2, true);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* DetectorConstruction::GetLogicDetector2()
{
  return roomdetLV;
}


G4LogicalVolume* DetectorConstruction::GetLogicDetector1()
{
  return absorberdetLV;
}


