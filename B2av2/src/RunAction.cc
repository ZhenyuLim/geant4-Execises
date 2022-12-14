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
/// \file RunAction.cc
/// \brief Implementation of the B2::RunAction class

#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4AccumulableManager.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{
  // set printing event number per each 100 events
  G4RunManager::GetRunManager()->SetPrintProgress(100);
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fEdep1);
  accumulableManager->RegisterAccumulable(fEdep2);

///
   auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
  analysisManager->SetFileName("B2a");


  //??????Histogram
 
  analysisManager->CreateH1("fEdep1","Edep in absorberdet",100,0.,700*keV);//ID 0

  analysisManager->CreateH1("fEdep2","Edep in roomdetLV",100,0.,700*keV);// ID 1


  //?????? ntuple
  analysisManager->CreateNtuple("Mydata","Energy deposit");
  analysisManager->CreateNtupleDColumn("fEdep1");
  analysisManager->CreateNtupleDColumn("fEdep2");
  analysisManager->FinishNtuple();
 


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{   /* delete G4AnalysisManager::Instance();  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();

  accumulableManager->Reset();

///
  //??????root??????
   auto analysisManager = G4AnalysisManager::Instance();
   analysisManager->Reset();
    analysisManager->OpenFile("B2a"); 


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* )
{
  
    G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge(); 
  
  
  G4double redep1 = fEdep1.GetValue();
  G4double redep2 = fEdep2.GetValue();
  G4cout << "Total energy in absorberdet is " << redep1 << G4endl;
  G4cout << "Total energy in roomdet is " << redep2 << G4endl;

  //?????????????????????
   auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile(); 

////

/*   G4cout
     << G4endl
     << " Energy deposit is " << fEdep << "MeV "
     << G4endl; */




}
void RunAction::AddEdep1(G4double edep)
{
  fEdep1  += edep;

}
void RunAction::AddEdep2(G4double edep)
{
  fEdep2  += edep;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......






