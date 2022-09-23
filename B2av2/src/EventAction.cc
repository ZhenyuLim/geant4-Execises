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
/// \file EventAction.cc
/// \brief Implementation of the B2::EventAction class

#include "EventAction.hh"
#include "RunAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4AnalysisManager.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"

#include "G4HCofThisEvent.hh"
#include "TrackerHit.hh"
#include "TrackerSD.hh"

#include "Run.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction)
: fRunAction(runAction)
{
  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerHitsCollection*
EventAction::GetHitsCollection(G4int hcID,const G4Event* event) const
{

  auto hitsCollection
    = static_cast<TrackerHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));


return hitsCollection;
}

void EventAction::BeginOfEventAction(const G4Event*)
{

   fEdep1 = fEdep2 = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  // get number of stored trajectories

  G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  // periodic printing

  G4int eventID = event->GetEventID();
//  if ( eventID < 100 || eventID % 100 == 0) {
    G4cout << ">>> Event: " << eventID  << G4endl;
 //   if ( trajectoryContainer ) {
      G4cout << "    " << n_trajectories
             << " trajectories stored in this event." << G4endl;
//    }
    G4VHitsCollection* hc = event->GetHCofThisEvent()->GetHC(0);
    G4cout << "    "
           << hc->GetSize() << " hits stored in this event" << G4endl;


   fTkHCID1
      = G4SDManager::GetSDMpointer()->GetCollectionID("HitsCollection1");
    // Get hits collections
   auto TkHC1 = GetHitsCollection(fTkHCID1, event); 
   G4int n1 = TkHC1->entries();

   // Get hit with total values

       for (G4int i=0;i<n1;i++){
      auto TkHit = (*TkHC1)[i];
          fEdep1 += TkHit->GetEdep();
          
     }
         G4cout << "    "
           <<" sum in absorberdet is "<< fEdep1 << G4endl; 
     fRunAction->AddEdep1(fEdep1); 
///////////////////////////////////////////////////
   fTkHCID2
      = G4SDManager::GetSDMpointer()->GetCollectionID("HitsCollection1");
    // Get hits collections
   auto TkHC2 = GetHitsCollection(fTkHCID2, event); 
   G4int n2 = TkHC2->entries();

   // Get hit with total values

       for (G4int i=0;i<n2;i++){
      auto TkHit = (*TkHC2)[i];
          fEdep2 += TkHit->GetEdep();
          
     }
         G4cout << "    "
           <<" sum in roomdet is "<< fEdep2 << G4endl; 
     fRunAction->AddEdep2(fEdep2); 
//  }

    


///

   //填充Histogram
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillH1(0,fEdep1);
  analysisManager->FillH1(1,fEdep2);
  //填充 ntuple
  analysisManager->FillNtupleDColumn(0,fEdep1);
  analysisManager->FillNtupleDColumn(1,fEdep2);

  analysisManager->AddNtupleRow();



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



