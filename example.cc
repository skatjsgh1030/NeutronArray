#include "globals.hh"

#include "G4RunManager.hh"
#include "G4VisExecutive.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"

#include "QGSP_BERT.hh"
#include "OTDetectorConstruction.hh"
#include "OTPrimaryGeneratorAction.hh"
#include "OTRunAction.hh"
#include "OTSteppingAction.hh"

int main(int argc, char** argv)
{
  G4RunManager* runManager = new G4RunManager;

  G4VModularPhysicsList* physicsList = new QGSP_BERT;
  runManager -> SetUserInitialization(physicsList);
  runManager -> SetUserInitialization(new OTDetectorConstruction());
  runManager -> SetUserAction(new OTPrimaryGeneratorAction());
  runManager -> SetUserAction(new OTRunAction());
  runManager -> SetUserAction(new OTSteppingAction());
  runManager -> Initialize();
  
  G4VisManager* visManager = new G4VisExecutive;
  visManager -> Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (argc != 1)
  {

    G4UIExecutive* ui = new G4UIExecutive(argc, argv);

    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager -> ApplyCommand(command+fileName);

    if (argc == 3)
    {

      G4String command1 = "/run/beamOn ";
      G4String fileName1 = argv[2];
      UImanager -> ApplyCommand(command1+fileName1);
    }
    //ui -> SessionStart();
    delete ui;
  }
  else 
  {
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    UImanager -> ApplyCommand("/control/execute vis.mac"); 
    ui -> SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;

  return 0;
}
