#include "OTDetectorConstruction.hh"
#include "G4SubtractionSolid.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UniformElectricField.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include <iostream>
#include "TMath.h"
//#include "G4RootAnalysisManager.hh"
//#include "B4Analysis.hh"//analysis manager header
//#include "g4root.hh"

/*B4RunAction::B4RunAction()
 : G4UserRunAction()
 {
 	//Create analysis manager
	auto analysisManager = G4AnalyshsanalysisManager::Instance();
	analysisManager -> SetVerboseLevel(1);
	analysisManager -> SetFirstHistoId(1);
}

B4RunAction::~B4RunAction()
{
	delete G4AnalysisManager::Instance();
}*/

using std::stringstream;
using namespace std;
using namespace CLHEP;
OTDetectorConstruction::OTDetectorConstruction()
: G4VUserDetectorConstruction()
{
}

OTDetectorConstruction::~OTDetectorConstruction()
{
}

G4VPhysicalVolume* OTDetectorConstruction::Construct()
{  
  G4NistManager* nist = G4NistManager::Instance();
  
 const G4double labTemp = CLHEP::STP_Temperature + 20.*kelvin;

  G4Element*	elC  = new G4Element("Carbon",   "C",  6,  12.011*g/mole);
  G4Element*	elH  = new G4Element("Hydrogen", "H",  1,  1.00794*g/mole);

  G4Material*	Methylene = new G4Material("Methylene", 6.262e-04*g/CLHEP::cm3, 2, kStateGas, labTemp);
  Methylene -> AddElement(elC, 1);
  Methylene -> AddElement(elH, 2);

  // -----------------------------------------------------
  // World

  G4Material* world_mat = nist -> FindOrBuildMaterial("G4_AIR");
  G4double world_size = 2000*mm;

  G4Box* solidWorld =    
    new G4Box("World",                       // its name
              5*world_size,                // half x
              5*world_size,                // half y
              5*world_size);               // half z
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      true);                 //overlaps checking


  // -----------------------------------------------------
  // Detector

  G4Material* scintillator_mat = nist -> FindOrBuildMaterial("G4_XYLENE");
  G4Material* veto_mat = nist -> FindOrBuildMaterial("G4_POLYETHYLENE");
  G4Material* tube_mat = nist -> FindOrBuildMaterial("G4_Pyrex_Glass");
  G4Material* target_mat = nist -> FindOrBuildMaterial("Methylene");
  G4double detector_size = 1*mm;
  G4double scintillator_offset_z = (63.5/2+3013)*mm;
  G4double veto_offset_z = (3005-420)*mm;
  G4double tube_offset_z = (63.5/2+3010)*mm;

  
  G4double sizeX = 2000 * detector_size;
  G4double sizeY = 2000 * detector_size;
  G4double sizeZ = 63.5 * detector_size;
//target
  G4Box* target =        
    new G4Box("target",
              150*mm/2,
              150*mm/2,
              50*mm/2);

 /* G4RotationMatrix* yRot = 
  	new G4RotationMatrix;
	yRot ->rotateY(M_PI/4.*rad);

  G4RotationMatrix invRot = yRot -> invert();
  G4Transform3D transform(invRot,0);*/
          
  G4LogicalVolume* logicalDetector3 =    
    new G4LogicalVolume(target,
                        target_mat,
                        "target");
    
    new G4PVPlacement(0,
                      G4ThreeVector(0,0,0),
                      logicalDetector3,
                      "target",
                      logicWorld,
                      false,
                      1,  
                      true);

//hollow pyrex glass tube
  G4Box* outerTube =    
    new G4Box("outerTube",
              sizeX/2,
              76.2*mm/2,
              sizeZ/2);

  G4Box* innerTube = 
  	new G4Box("innerTube",
			  (2000*mm-3*mm)/2,
			  (76.2*mm-3*mm)/2,
			  (63.5*mm-3*mm)/2);

  G4SubtractionSolid*hollowTube = 
  	new G4SubtractionSolid("hollowTube",
							outerTube,
							innerTube);
  
	G4double AngleNN = -27.07*CLHEP::degree;
	G4RotationMatrix* yRot_SS = new G4RotationMatrix;
								yRot_SS -> rotateY(AngleNN);
  G4LogicalVolume* logicalDetector = 
    new G4LogicalVolume(hollowTube,
                        tube_mat,
                        "hollowTube");
	for (G4int copyN =0;copyN < 24;copyN++)
	{
    new G4PVPlacement(yRot_SS,
                      G4ThreeVector(-1*tube_offset_z*TMath::Sin(AngleNN),(876.3-76.2*copyN)*mm,tube_offset_z*TMath::Cos(AngleNN)),
                      logicalDetector,
                      "hollowTube",
                      logicWorld,
                      false,
                      1,  
                      true);
	}
//liquide Scintillator
	G4double AngleN = -27.07*CLHEP::degree;
	G4RotationMatrix* yRot_S = new G4RotationMatrix;
								yRot_S -> rotateY(AngleN);
  G4Box* lScintillator =
    new G4Box("lScintillator",
              (2000*mm-3*mm)/2,
              (76.2*mm-3*mm)/2,
              (63.5*mm-3*mm)/2);     
 
  G4LogicalVolume* logicalDetector1 =                         
    new G4LogicalVolume(lScintillator,
                        scintillator_mat,
                        "lScintillator");
    for(G4int copyN_L =0; copyN_L < 24; copyN_L++)
	{
    new G4PVPlacement(yRot_S,
                      G4ThreeVector(-1*scintillator_offset_z*TMath::Sin(AngleN),(876.3-76.2*copyN_L)*mm,scintillator_offset_z*TMath::Cos(AngleN)),
                      logicalDetector1,
                      "lScintillator",
                      logicWorld,
                      false,
                      copyN_L+24,
                      true);
	}
//Veto Wall(front)
	G4double AngelF = -27.07*CLHEP::degree;
	const G4int vetoFN = 12;
	G4RotationMatrix* yRot = new G4RotationMatrix;
							yRot ->rotateY(AngelF);
  G4Box* veto =    
    new G4Box("veto",
              94*mm/2,
              sizeY/2,
              10*mm/2);
      
  G4LogicalVolume* logicalDetector2 =                         
    new G4LogicalVolume(veto,
                        veto_mat,
                        "veto");
	for(G4int copyN = 0; copyN < vetoFN; copyN++)
	{
                                   
   		 new G4PVPlacement(yRot,//pi,theta,psi
           		           G4ThreeVector(-1*(veto_offset_z*TMath::Sin(AngelF)+(1001-184*copyN)*mm*TMath::Cos(AngelF)),0,veto_offset_z*TMath::Cos(AngelF)-(1001-184*copyN)*mm*TMath::Sin(AngelF)),
                	       logicalDetector2,
                   		   "veto",
                   		   logicWorld,
                   		   false,
                    		  copyN*2,
                    		  true);
	}

//Veto Wall(back)
	G4double AngelB = -27.07*CLHEP::degree;
	const G4int vetoBN = 12;
	G4RotationMatrix* yBRot = new G4RotationMatrix;
							yBRot -> rotateY(AngelB);
  G4Box* vetoB =    
    new G4Box("vetoB",
              94*mm/2,
              sizeY/2,
              10*mm/2);
      
  G4LogicalVolume* logicalDetector4 =                         
    new G4LogicalVolume(vetoB,
                        veto_mat,
                        "vetoB");
   for(G4int copyN =0; copyN < vetoBN; copyN++)
   {
    new G4PVPlacement(yBRot,
                      G4ThreeVector(-1*((veto_offset_z+10*mm)*TMath::Sin(AngelB)+(817+94-184*copyN)*mm*TMath::Cos(AngelB)),0,(veto_offset_z+10*mm)*TMath::Cos(AngelB)-(817+94-184*copyN)*mm*TMath::Sin(AngelB)),
                      logicalDetector4,
                      "vetoB",
                      logicWorld,
                      false,
                      (copyN*2)+1,
                      true);
	}
/*G4VisAttributes* Red        = new G4VisAttributes( G4Colour(225/225. ,0/225.   , 0/225.   ));
G4VisAttributes* Yellow     = new G4VisAttributes( G4Colour(225/225. ,225/225. , 0/225.   ));
G4VisAttributes* LightBlue  = new G4VisAttributes( G4Colour(0/225.   ,204/225. , 204/225. ));
G4VisAttributes* LightGreen = new G4VisAttributes( G4Colour(153/225. ,255/225. , 153/225. ));*/

logicWorld -> SetVisAttributes (G4VisAttributes::GetInvisible());

G4VisAttributes* targetVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
targetVisAtt -> SetForceWireframe(true);
logicalDetector3 -> SetVisAttributes (targetVisAtt);
G4VisAttributes* vetoVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
vetoVisAtt -> SetForceWireframe(true);
logicalDetector2 -> SetVisAttributes (vetoVisAtt);
G4VisAttributes* vetoBVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
vetoBVisAtt -> SetForceWireframe(true);
logicalDetector4 -> SetVisAttributes (vetoBVisAtt);
G4VisAttributes* hollowTubeVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
hollowTubeVisAtt -> SetForceWireframe(true);
logicalDetector -> SetVisAttributes (hollowTubeVisAtt);

  return physWorld;
}
