#define vtkRenderingCore_AUTOINIT 2(vtkInteractionStyle, vtkRenderingOpenGL2)
//#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)
#include "vtkETHZGeoPhysicsHDF5SphericalReader.h"
#include "vtkProjMap.h"

#include "vtkGeometryFilter.h"
#include "vtkExtractGrid.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSmartPointer.h"
#include "vtkLookupTable.h"

#include <vtksys/SystemTools.hxx>
#include <vtksys/CommandLineArguments.hxx>
#define REST

#define VTK_CREATE(type, var) \
  vtkSmartPointer<type> var = vtkSmartPointer<type>::New();

int
main(int argc, char **argv)
{
  std::string filein, fileout;
  int compute=1;
  bool shuffle = false;
  bool full_vorticity = false;
  vtksys::CommandLineArguments args;
  args.Initialize(argc, argv);
  args.AddArgument(
    "-i", vtksys::CommandLineArguments::SPACE_ARGUMENT, &filein, "(the name of the HDF5 file used as input)");
/*
  args.AddArgument(
    "-o", vtksys::CommandLineArguments::SPACE_ARGUMENT, &fileout, "(the name of the VTI file used as output)");

  args.AddArgument(
    "--fv", vtksys::CommandLineArguments::NO_ARGUMENT, &full_vorticity, "(set if writing the full vorticity vector)");
*/
      //args.AddArgument(
    //"-c", vtksys::CommandLineArguments::SPACE_ARGUMENT, &compute, "(which field to evaluate, default is vorticity, 1 is vel_phi)");
    

  if ( !args.Parse() || argc == 1 || filein.empty() )
    {
    cerr << "\nTestJacksonReader: Written by Jean M. Favre at the Swiss National Supercomputing Center (CSCS)\n"
         << "options are:\n";
    cerr << args.GetHelp() << "\n";
    exit(1);
    }
  if(!vtksys::SystemTools::FileExists(filein.c_str()))
    {
    cerr << "\nFile " << filein.c_str() << " does not exist\n\n";
    exit(1);
    }
  VTK_CREATE(vtkETHZGeoPhysicsHDF5SphericalReader, reader);
    reader->DebugOff();
  cerr << "\nsetting File " << filein.c_str() << "\n";
  reader->SetFileName(filein.c_str());
 
  reader->UpdateInformation();
  // all should be available to computer Temperature
  //reader->DisableAllPointArrays();
  //reader->EnablePointArray("Entropy");
  reader->EnablePointArray("Velocity");
  //reader->Update();
  //exit(1);

  VTK_CREATE(vtkExtractGrid, voi);
  voi->SetInputConnection(reader->GetOutputPort(0));
  voi->SetVOI(
             0, 192,
             0, 193,
	    100, 100
             );

  VTK_CREATE(vtkProjMap, cart_map);
  cart_map->SetInputConnection(voi->GetOutputPort(0));
  cart_map->Update();
#ifdef REST
  VTK_CREATE(vtkGeometryFilter, geom2);
  geom2->SetInputConnection(cart_map->GetOutputPort(0));

  VTK_CREATE(vtkLookupTable, lut);
  lut->SetTableRange(1.32e-02, 4.9e-01);
  lut->SetHueRange(0.66,0.0);
  lut->SetNumberOfTableValues(256);
  lut->SetVectorModeToMagnitude();
  lut->Build();
  
  VTK_CREATE(vtkPolyDataMapper, ocMapper);
  ocMapper->SetInputConnection(geom2->GetOutputPort(0));
  ocMapper->ScalarVisibilityOn();
  ocMapper->SetColorModeToMapScalars();
  ocMapper->SetLookupTable(lut);
  ocMapper->UseLookupTableScalarRangeOn();
  ocMapper->SetScalarModeToUsePointFieldData();
  ocMapper->SelectColorArray("Velocity");

  ocMapper->ImmediateModeRenderingOn();
  ocMapper->DebugOff();

  VTK_CREATE(vtkActor, ocActor);
  ocActor->SetMapper(ocMapper);
  ocActor->DebugOff();

  VTK_CREATE(vtkRenderer, ren);
  VTK_CREATE(vtkRenderWindow, renWin);
  VTK_CREATE(vtkRenderWindowInteractor, iren);
  renWin->DebugOff();

  iren->SetRenderWindow(renWin);
  renWin->AddRenderer(ren);
  
  ren->AddActor(ocActor);

  renWin->SetSize(512, 512);
  renWin->Render();
  iren->Start();
#endif
}
