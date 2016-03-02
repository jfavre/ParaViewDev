/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkETHZBioMedicalFAIMReader.cxx,v $

  Copyright (c) Patrick Boenzli
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/


#include "vtkETHZBioMedicalFAIMReader.h"

#include "vtkETHZBioMedicalASCIIFile.cxx"
#include "vtkETHZBioMedicalASCIIFileFilter.cxx"
#include "vtkETHZBioMedicalASCIIFileSection.cxx"

//#include "vtkSystemIncludes.h"
#include <string>
#include <vtksys/SystemTools.hxx>

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDataArraySelection.h"
#include "vtkDataSetAttributes.h"
#include "vtkDoubleArray.h"
#include "vtkErrorCode.h"
#include "vtkFieldData.h"
#include "vtkFloatArray.h"
#include "vtkIdTypeArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnstructuredGrid.h"

#include <limits>
		   
//vtkCxxRevisionMacro(vtkETHZBioMedicalFAIMReader, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkETHZBioMedicalFAIMReader);

/**
 * request informations from files
 *
 *
 * this request asks the algorithm to provide as much information as it can about what the output data will look like once the algorithm has generated it.
 * get as much information as one can without actually reading in the entire data file and without taking up significant CPU time.
 *
 * @param request: the request chich specifies what the user is asking this algorithm to do. Eg. REQUEST_INFORMATION
 * @param inputVector: information vector for the input to this algorighm
 * @param outputVector: informaion vector for the output of this algorighm
 *
 * reader reads the header information from the file to get what inormation it can out of it.
 */
int vtkETHZBioMedicalFAIMReader::RequestInformation(vtkInformation* request,
    vtkInformationVector** vtkNotUsed(inputVector),
    vtkInformationVector* outputVector)
{
  // enable the displacement array
  this->PointDataArraySelection->AddArray("Nodal displacements");
  // enable the nodal forces array
  this->PointDataArraySelection->AddArray("Nodal forces");
  // enable border condition array
  this->PointDataArraySelection->AddArray("Boundary Conditions Fixed");
  this->PointDataArraySelection->DisableArray("Boundary Conditions Fixed");
  // enable border condition array
  this->PointDataArraySelection->AddArray("Boundary Conditions Translate");
  this->PointDataArraySelection->DisableArray("Boundary Conditions Translate");

  // add the strain and stress arrays
  for(int i=0; i < 7; i++)
  {
    this->CellDataArraySelection->AddArray(strain_var_names[i]);
  }
  for(int i=0; i < 7; i++)
  {
    this->CellDataArraySelection->AddArray(stress_var_names[i]);
  }
  
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(CAN_HANDLE_PIECE_REQUEST(), 1);
  //outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);
  return 1;
}


/**
 * request data from the file
 *
 * this function is called at the end and should return the informations the pipeline needs for further processing. the file is read in.
 */
int vtkETHZBioMedicalFAIMReader::RequestData(
  vtkInformation* request,
  vtkInformationVector** vtkNotUsed(inputVector),
  vtkInformationVector* outputVector)
{
  if (!this->FileName)
  {
    vtkErrorMacro(<< "no valid filename supported");
    return 0;
  }

  vtkDebugMacro( << "RequestData(BEGIN)");


  // ===============================================================================================================
  // initializing the file reader
  // ===============================================================================================================
  vtkETHZBioMedicalASCIIFile* f = new vtkETHZBioMedicalASCIIFile;
  //vtkETHZBioMedicalASCIIFile* f = vtkETHZBioMedicalASCIIFile::New();

  std::string   path = vtksys::SystemTools::GetFilenamePath(this->FileName);
  std::string   prefix = vtksys::SystemTools::GetFilenameWithoutLastExtension(this->FileName);
  std::string meshFileName = path + "/" + prefix + std::string(".mesh");
  cout << "mesh filename    = " << meshFileName<< endl;
  cout << "disp[s] filename = " << this->FileName<< endl;

  f->setFileName(meshFileName);
  if( f->processFile() == -1)
    {
    vtkErrorMacro(<< "Could not open mesh file " << meshFileName.c_str() << " for reading" << endl);
    return -1;
    }

  float progress = 0.;

  // check if there user canceled the operation
  this->UpdateProgress (++progress/13.0);
  if (this->GetAbortExecute())
    return -1;

  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()));

  int piece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  //int numPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  if(piece == 0)
    {  // for the moment, we only know how to ingest the data on one proc.
       // User should use D3 to distribute the data amongst processes
  int elementType;
  int nrOfInPnt;
  int nrOfComponents;

  // read in the information section
  f->selectSection( "Begin model");
  f->readLine(&elementType, &this->NbCells, &this->NbNodes, &nrOfInPnt, &nrOfComponents);
  cout << "nr of nodes: " << NbNodes << " - nr of cells: " << NbCells << endl;
  //cout << "nr of InPnt: " << nrOfInPnt << " - nr of Components: " << nrOfComponents << endl;

  double ElasticModulus, dummy, Poisson;
  int nprops, ntype;
  f->readLine(&dummy, &dummy, &dummy);
  f->readLine(&nprops, &ntype);
  f->readLine(&ElasticModulus, &Poisson);
  cout << "Elastic modulus = " << ElasticModulus << ", Poisson = " << Poisson << endl;
  // ===============================================================================================================
  // reading in NODES
  // ===============================================================================================================

  vtkDoubleArray *coords = vtkDoubleArray::New();
  coords->SetNumberOfComponents(nrOfComponents);
  coords->SetNumberOfTuples(this->NbNodes);

  double* nodes;                                               //!< the nodes from the file
  nodes = ((vtkDoubleArray *)coords)->GetPointer(0);

  // read in the nodes
  f->selectSection("Nodes");
  f->setWriteMode(vtkETHZBioMedicalASCIIFileSection::WM_WRITE);
  f->setPattern("%lf %lf %lf");
  f->read(&nodes);

  // write the result back to the information array
  vtkPoints *points = vtkPoints::New();
  points->SetData(coords);
  coords->Delete();
  output->SetPoints(points);
  points->Delete();

  // Call FastDelete() instead of Delete() to avoid garbage
  // collection checks. This improves the preformance significantly

  // check if there user canceled the operation
  this->UpdateProgress (++progress/13.0);
  if (this->GetAbortExecute())
    return -1;

  // ===============================================================================================================
  // reading in ELEMENTS
  // ===============================================================================================================

  vtkIdType* elements;
  //cerr << "sizeof(vtkIdType) = " << sizeof(vtkIdType) <<  endl;

  vtkIdTypeArray *vtkListCells = vtkIdTypeArray::New();
  vtkListCells->SetNumberOfValues(NbCells * (8 + 1));
  elements = vtkListCells->GetPointer(0);

  // now read out the elements
  f->selectSection("Elements");
  f->setWriteMode(vtkETHZBioMedicalASCIIFileSection::WM_WRITE);
  f->setPattern("%i %i %i %i %i %i %i %i");
  if(sizeof(vtkIdType) == 8)
     {
     long *lp = reinterpret_cast<long*>(elements);
     f->read(&lp);
     }
  else if(sizeof(vtkIdType) == 4)
     {
     int *lp = reinterpret_cast<int*>(elements);
     f->read(&lp);
     }

  // some structures for memory restructuring
  vtkIdType *inptr = &elements[8 * NbCells] - 1; // absolute last mem
  vtkIdType *destptr = &elements[9 * NbCells] - 1; // absolute last mem

  for(int i= NbCells - 1 ; i >= 0; i--)
  {
    for(int j = 7; j >= 0; j--)
    {
      *destptr-- = *inptr-- - 1;
      // substract 1 because VTK indices start at 0 not 1
    }
    *destptr-- = 8;
  }
/*
  destptr = elements; // absolute first mem
  for(int i=0; i < NbCells; i++)
  {
    for(int j = 0; j < 9; j++)
    {
      cout << *destptr++ << " ";
    }
    cout <<  "\n";
  }
*/

  vtkCellArray *cells = vtkCellArray::New();

  cells->SetCells(NbCells, vtkListCells);
  vtkListCells->Delete();
  // Call FastDelete() instead of Delete() to avoid garbage
  // collection checks. This improves the preformance significantly

  int *types = new int[NbCells];

  for(int i=0; i < NbCells; i++)
  {
    types[i] = VTK_HEXAHEDRON;
  }

  output->SetCells(types, cells);
  cells->Delete();

  // Call FastDelete() instead of Delete() to avoid garbage
  // collection checks. This improves the preformance significantly
  delete [] types;

  // check if there user canceled the operation
  this->UpdateProgress (++progress/13.0);
  if (this->GetAbortExecute())
    return -1;

  // ===============================================================================================================
  // reading in BOUNDARY CONDITIONS: Fixed
  // ===============================================================================================================

  if( this->PointDataArraySelection->ArrayIsEnabled("Boundary Conditions Fixed") ||
      this->PointDataArraySelection->ArrayIsEnabled("Boundary Conditions Translate"))
  {
    // now read out the elements
    f->selectSection("# Boundary conditions");

    int nrOfBoundaryConditions;     //!< the number of boundary conditions
    int nrOfBCPerLine = 5;          //!< this is the number of boundary conditions per line
    int* boundaries;                //!< temp buffer for the boundary condition
    // first read out the number of elements having boundary conditions
    f->setWriteMode(vtkETHZBioMedicalASCIIFileSection::WM_WRITE);
    f->readLine(&nrOfBoundaryConditions);

    // now read out the data, reader shall allocate the memory for itself
    f->setWriteMode(vtkETHZBioMedicalASCIIFileSection::WM_CREATE);
    f->setLinesToRead( (int)ceil(nrOfBoundaryConditions / nrOfBCPerLine + 0.5));
    f->read(&boundaries);

    // now creat the vtk container for the data
    vtkIntArray *bc = vtkIntArray::New();
    bc->SetNumberOfComponents(nrOfComponents);
    bc->SetNumberOfTuples(this->NbNodes);
    bc->SetName("Boundary Conditions Fixed");
    cerr << "allocating " << this->NbNodes << " * 3 * ints for Boundary Conditions" << endl;

    // init all elements with zero
    int empty[] = {0, 0, 0};
    for( int i = 0; i < this->NbNodes; i++)
      bc->SetTupleValue(i, empty);

    // now copy the values from the temp array to the member array
    for( int i = 0; i < nrOfBoundaryConditions * 4; i+=4)
    {
      // making 1 -> 0 and 0 -> 1
      boundaries[i+1] = -boundaries[i+1] + 1;
      boundaries[i+2] = -boundaries[i+2] + 1;
      boundaries[i+3] =  boundaries[i+3] - 1;

      bc->SetTupleValue(boundaries[i]-1, &boundaries[i+1]);
    }

    // check if the data is read correctly
    f->readLine(&nrOfBoundaryConditions);
    if( nrOfBoundaryConditions != 0)
      cerr << "error: the boundary conditions have an unknown format" << endl;

    // ===============================================================================================================
    // reading in BOUNDARY CONDITIONS: Translation
    // ===============================================================================================================

    // first read out the number of elements having translation boundary conditions
    f->setWriteMode(vtkETHZBioMedicalASCIIFileSection::WM_WRITE);
    f->readLine(&nrOfBoundaryConditions);

    double* translation;

    // now read out the data, reader shall allocate the memory for itself
    f->setWriteMode(vtkETHZBioMedicalASCIIFileSection::WM_CREATE);
    f->setPattern("%lf %lf %lf");
    f->read(&translation);

    //allocating necessary vtk container
    vtkDoubleArray *trans = vtkDoubleArray::New();
    trans->SetNumberOfComponents(3);
    trans->SetNumberOfTuples(this->NbNodes);
    trans->SetName("Boundary Conditions Translation");
     cerr << "allocating " << this->NbNodes << " * 3 * doubles for Boundary Conditions" << endl;

    double zerosDouble[3] = {0., 0., 0.};
    // init all elements with zero
    for( int i = 0; i < this->NbNodes; i++)
      trans->SetTupleValue(i, zerosDouble);

    // now copy the values from the temp array to the member array,
    // take care of all cases, where there are more translations than just one
    // for each dimension => read out value and write it back
    for( int i = 0; i < nrOfBoundaryConditions * 3; i+=3)
    {
      double coord[3];
      trans->GetTupleValue((int)translation[i]-1, coord);
      coord[ (int)translation[i+1] - 1] = translation[i+2];
      trans->SetTupleValue((int)translation[i]-1, coord);
    }

    // remove fixations with trailing translations
    for( int i = 0; i < this->NbNodes; i++)
    {
      // get the values
      int fixation[3];
      bc->GetTupleValue(i, fixation);
      double translation[3];
      trans->GetTupleValue(i, translation);

      // compare them
      if( translation[0] != 0. ||
          translation[1] != 0. ||
          translation[2] != 0.)
      {
        bc->SetTupleValue(i, empty);
      }
    }


    // pack the data to the array
    output->GetPointData()->AddArray(bc);
    // remove unused data
    bc->Delete();
    // now remove the tmp buffer again
    delete [] boundaries;

    // pack the data to the array
    output->GetPointData()->AddArray(trans);
    // remove unused data
    trans->Delete();
    // now remove the tmp buffer again
    delete [] translation;
  }

    // check if there user canceled the operation
    this->UpdateProgress (++progress/13.0);
    if (this->GetAbortExecute())
      return -1;

  // close the file again
  f->close();
  //f->Delete();
  delete f;

  // ===============================================================================================================
  // opening DISP FILE for further processing
  // ===============================================================================================================

  // create the file reader
  //f = vtkETHZBioMedicalASCIIFile::New();
  f = new vtkETHZBioMedicalASCIIFile;

  f->setFileName(this->FileName); // we now use the filename selected in browser
  // and read in the file
  if( f->processFile() == -1)
  {
    vtkErrorMacro(<< "Could not open file " << this->FileName << " for reading, only mesh informations loaded" << endl);
    return -1;
  }

  // check if there user canceled the operation
//   //  this->UpdateProgress (++progress/13.0);
  if (this->GetAbortExecute())
    return -1;

  // ===============================================================================================================
  // reading in NODAL DISPLACEMENTS
  // ===============================================================================================================
  if(this->GetPointArrayStatus("Nodal displacements"))
    {
    // allocate the memory for these vectors
    vtkDoubleArray *displacements = vtkDoubleArray::New();
    displacements->SetNumberOfComponents(3);
    displacements->SetNumberOfTuples(this->NbNodes);
    displacements->SetName("Nodal displacements");
    cerr << "allocating " << this->NbNodes << " * 3 * doubles for displacement\n";

    // get a simple pointer to the vtk data structure
    double* displacementsData = displacements->GetPointer(0);

    // read in the nodes
    f->selectSection("Nodal displacements");
    f->setWriteMode(vtkETHZBioMedicalASCIIFileSection::WM_WRITE);
    // set the pattern of the data (first entry is empty, the following tripple is used for coordinates)
    f->setPattern("%*d %lf %lf %lf");
    f->read(&displacementsData);

    output->GetPointData()->AddArray(displacements);
    displacements->Delete();
    output->GetPointData()->SetActiveAttribute("Nodal displacements", vtkDataSetAttributes::VECTORS);
    }

  // check if there user canceled the operation
   this->UpdateProgress (++progress/13.0);
  if (this->GetAbortExecute())
    return -1;

  // ===============================================================================================================
  // reading in NODAL FORCES
  // ===============================================================================================================
  if(this->GetPointArrayStatus("Nodal forces"))
    {
    vtkDoubleArray *forces = vtkDoubleArray::New();
    forces->SetNumberOfComponents(3);
    forces->SetNumberOfTuples(this->NbNodes);
    forces->SetName("Nodal forces");
    cerr << "allocating " << this->NbNodes << " * 3 * doubles for forces\n";

    // get a simple pointer to the data
    double* forcesData = forces->GetPointer(0);

    f->selectSection("# Nodal reaction forces");
    f->setWriteMode(vtkETHZBioMedicalASCIIFileSection::WM_WRITE);
    // set the pattern of the data (first entry is empty, the following tripple is used for vector coordinates)
    f->setPattern("%*d %lf %lf %lf");
    f->read(&forcesData);

    // clean up unused memory
    output->GetPointData()->AddArray(forces);
    forces->Delete();
    }

  // check if there user canceled the operation
  this->UpdateProgress (++progress/13.0);
  if (this->GetAbortExecute())
    return -1;

  // ===============================================================================================================
  // reading in strains
  // ===============================================================================================================

  // do we need to read in the stress variables (e11, e22, e33, e12, e23, e31, SED are 7 variables)
  bool NeedToReadStrain = false;
  for(int i = 0; i < 7; i++)
    {
    if(this->CellDataArraySelection->GetArraySetting(i))
      {
      NeedToReadStrain = true;
      }
    }
  if(NeedToReadStrain)
    {
    double* dataPointer = NULL;
    f->selectSection("# Element strain");
    // this time allocate the memory for yourself
    f->setWriteMode(vtkETHZBioMedicalASCIIFileSection::WM_CREATE);
    // set the pattern of the data (first two entries are empty, then we read in e11, e22, e33, e12, e23, e31, SED)
    f->setPattern("%*d %*d %lf %lf %lf %lf %lf %lf %lf");

    // now read in the data finally
    f->read(&dataPointer);

    // now add each array seperately if enabled in the GUI
    for( int i = 0; i < 6; i++)
      {
      if( this->CellDataArraySelection->GetArraySetting(i))
        {
        vtkDoubleArray *strain = vtkDoubleArray::New();
        strain->SetNumberOfComponents(1);
        strain->SetNumberOfTuples(this->NbCells);
        strain->SetName(this->CellDataArraySelection->GetArrayName(i));
        cerr << "allocating " << this->NbCells << " * 1 * doubles for " <<  this->CellDataArraySelection->GetArrayName(i)<< endl;

        // get a simple pointer to the vtk data structure
        double* copyArray;                                               //!< the nodes from the file
        copyArray = strain->GetPointer(0);

        // now copy the values in this array
        for( int j = 0; j < this->NbCells; j++)
        {
          copyArray[j] = dataPointer[i + j*7];
        }

        output->GetCellData()->AddArray(strain);
        strain->Delete();
        }
      }

    if(this->GetCellArrayStatus("SED"))
      {
      vtkDoubleArray *sed = vtkDoubleArray::New();
      sed->SetNumberOfComponents(1);
      sed->SetNumberOfTuples(this->NbCells);
      sed->SetName("SED");
      cerr << "allocating " << this->NbCells << " * 1 * doubles for SED" << endl;

      // get a simple pointer to the vtk data structure
      double* sedArray = sed->GetPointer(0);

      // now copy the values in this array
      for( int j = 0; j < this->NbCells; j++)
        {
        sedArray[j] = dataPointer[6 + j*7];
        }

      output->GetCellData()->AddArray(sed);

// sedArray still points to the SED values. We'll use it to evaluate EFF
      vtkDoubleArray *eff = vtkDoubleArray::New();
      eff->SetNumberOfComponents(1);
      eff->SetNumberOfTuples(this->NbCells);
      eff->SetName("EFF");
      cerr << "allocating " << this->NbCells << " * 1 * doubles for EFF" << endl;

        // get a simple pointer to the vtk data structure
      double* effArray = eff->GetPointer(0);

      for( int j = 0; j < this->NbCells; j++)
        {
        effArray[j] = sqrt(2.0*sedArray[j]/ElasticModulus);
        }

      output->GetCellData()->AddArray(eff);
      sed->Delete();
      eff->Delete();
      }
      delete [] dataPointer;
    }

  // check if there user canceled the operation
  this->UpdateProgress (++progress/13.0);
  if (this->GetAbortExecute())
    return -1;

  // ===============================================================================================================
  // reading in stress
  // ===============================================================================================================

  bool NeedToReadStress = false;
  for(int i = 0; i < 7; i++)
    {
    if(this->CellDataArraySelection->GetArraySetting(i+7))
      {
      NeedToReadStress = true;
      }
    }

  // do we need to read in the stress variables (s11,  s22 , s33  s12  s23  s31  vonMises)
  if(NeedToReadStress)
    {
    double* dataPointer = NULL;
    f->selectSection("# Element stress");
    // this time allocate the memory for yourself
    f->setWriteMode(vtkETHZBioMedicalASCIIFileSection::WM_CREATE);
    // set the pattern of the data (first two entries are empty, then we read in s11,  s22 , s33  s12  s23  s31  vonMises)
    f->setPattern("%*d %*d %lf %lf %lf %lf %lf %lf %lf");

    // now read in the data finally
    f->read(&dataPointer);

    // now add each array seperatley if enabled in the GUI
    for( int i = 0; i < 7; i++)
      {
      if( this->CellDataArraySelection->GetArraySetting(i+7))
        {
        vtkDoubleArray *stress = vtkDoubleArray::New();
        stress->SetNumberOfComponents(1);
        stress->SetNumberOfTuples(this->NbCells);
        stress->SetName(this->CellDataArraySelection->GetArrayName(i+7));
        cerr << "allocating " << this->NbCells << " * 1 * doubles for " <<  this->CellDataArraySelection->GetArrayName(i+7)<< endl;

        // get a simple pointer to the vtk data structure
        double* copyArray;                                               //!< the nodes from the file
        copyArray = stress->GetPointer(0);

        // now copy the values in this array
        for( int j = 0; j < this->NbCells; j++)
          {
          copyArray[j] = dataPointer[i + j*7];
          }

        output->GetCellData()->AddArray(stress);
        stress->Delete();
        }
      }
    delete [] dataPointer;
  }

  // check if there user canceled the operation
  this->UpdateProgress (++progress/13.0);

  // close file and delete it
  f->close();
  //f->Delete();
  delete f;
  } // only read on processor 0
  vtkDebugMacro( << "RequestData(END)");
  return 1;
}
