/*=========================================================================

  Program:   ParaView
  Module:    vtkGadgetReader.cxx


=========================================================================*/

#include "vtkGadgetReader.h"

#include "vtkDirectory.h"
#include "vtkCellType.h"
#include "vtkCellArray.h"
#include "vtkDataArray.h"
#include "vtkDataArraySelection.h"
#include "vtkFloatArray.h"
#include "vtkIdTypeArray.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#ifdef ALL_TYPES
#include "vtkMultiBlockDataSet.h"
#endif
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <vtksys/RegularExpression.hxx>
#include <vtksys/SystemTools.hxx>

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkGadgetReader, Controller, vtkMultiProcessController);
#endif

#include <algorithm>
#include <functional>
#include <map>

#include <hdf5.h>


static int ReadHDF5INT64Dataset(const char *name, hid_t mesh_id, long long *data)
{
  hid_t coords_id, filespace, attr1;
  herr_t  status;
  hsize_t dimsf[2];

  coords_id = H5Dopen(mesh_id, name, H5P_DEFAULT);
  filespace = H5Dget_space(coords_id);
  H5Sget_simple_extent_dims(filespace, dimsf, NULL);
  H5Sclose(filespace);

  if (H5Dread(coords_id, H5T_NATIVE_LONG, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT, data) < 0)
    return -1;

  H5Dclose(coords_id);
  return dimsf[0];
}

static int ReadHDF5Dataset(const char *name, hid_t mesh_id, float *data)
{
  hid_t coords_id, filespace, attr1;
  herr_t  status;
  hsize_t dimsf[2];
  //double CGSConversionFactor;
  //float aexpScaleExponent, hScaleExponent;

  coords_id = H5Dopen(mesh_id, name, H5P_DEFAULT);
  filespace = H5Dget_space(coords_id);
  H5Sget_simple_extent_dims(filespace, dimsf, NULL);
  H5Sclose(filespace);

  if (H5Dread(coords_id, H5T_NATIVE_FLOAT, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT, data) < 0)
    return -1;
/*
// read attributes of given dataset
        attr1 = H5Aopen_name(coords_id, "CGSConversionFactor");
        if (H5Aread(attr1, H5T_NATIVE_DOUBLE, &CGSConversionFactor) < 0)
    return -1;
        status = H5Aclose(attr1);

        attr1 = H5Aopen_name(coords_id, "aexp-scale-exponent");
        if (H5Aread(attr1, H5T_NATIVE_FLOAT, &aexpScaleExponent) < 0)
    return -1;
        status = H5Aclose(attr1);

        attr1 = H5Aopen_name(coords_id, "h-scale-exponent");
        if (H5Aread(attr1, H5T_NATIVE_FLOAT, &hScaleExponent) < 0)
    return -1;
        status = H5Aclose(attr1);
*/
// end of attributes read
// TODO Do coordinates have to be scaled with given attributes. Don't know yet.
    H5Dclose(coords_id);
    return dimsf[0];
}

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkGadgetReader);
//----------------------------------------------------------------------------
vtkGadgetReader::vtkGadgetReader()
{
  this->SetNumberOfInputPorts(0);
  this->DirectoryName = nullptr;
  this->CellType                 = 0;
  this->UpdatePiece              = 0;
  this->UpdateNumPieces          = 0;
  this->PointDataArraySelection  = vtkDataArraySelection::New();
  for (int i=0; i< 6; i++)
    {
    this->PartTypes[i] = false;
    this->NumPart_Total[i] = 0;
    }

  this->FieldArrays = {
    {"PartType0", {"Density", "Masses", "SmoothingLength", "Velocities", "ParticleIDs"} }, 
#ifdef ALL_TYPES
    {"PartType1", {"Masses", "Velocities", "ParticleIDs"} },
    {"PartType4", {"Masses", "Velocities", "ParticleIDs"} }
#endif
  };

#ifdef PARAVIEW_USE_MPI
  this->Controller = nullptr;
  this->SetController(vtkMultiProcessController::GetGlobalController());
#endif
}
//----------------------------------------------------------------------------
vtkGadgetReader::~vtkGadgetReader()
{
  this->CloseFile();

  delete [] this->DirectoryName;
  this->PointDataArraySelection->Delete();
  this->PointDataArraySelection = nullptr;

#ifdef PARAVIEW_USE_MPI
  this->SetController(NULL);
#endif
}
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
void vtkGadgetReader::CloseFile()
{

}
//----------------------------------------------------------------------------
int vtkGadgetReader::OpenFile()
{
  if (!this->DirectoryName)
    {
    vtkErrorMacro(<<"DirectoryName must be specified.");
    return 0;
    }

  if (FileModifiedTime>FileOpenedTime)
    {
    this->CloseFile();
    }

  return 1;
}

void vtkGadgetReader::SetDirectoryName(const char* dn)
{
  this->DirectoryName = strdup(vtksys::SystemTools::GetParentDirectory(dn).c_str());
}

//----------------------------------------------------------------------------
int vtkGadgetReader::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(CAN_HANDLE_PIECE_REQUEST(), 1);

#ifdef PARAVIEW_USE_MPI
  if (this->Controller)
    {
    this->UpdatePiece = this->Controller->GetLocalProcessId();
    this->UpdateNumPieces = this->Controller->GetNumberOfProcesses();
    }
#else
  this->UpdatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
#endif

    vtkDirectory* dir = vtkDirectory::New();
    int opened = dir->Open(this->DirectoryName);
    if (!opened)
    {
      vtkErrorMacro("Couldn't open " << this->DirectoryName);
      dir->Delete();
      return 0;
    }
    vtkIdType numFiles = dir->GetNumberOfFiles();

    this->GadgetFileNames.clear();

    for (vtkIdType i = 0; i < numFiles; i++)
    {
      if (strcmp(dir->GetFile(i), ".") == 0 ||
          strcmp(dir->GetFile(i), "..") == 0)
      {
        continue;
      }

      std::string fileString = this->DirectoryName;
      fileString += "/";
      fileString += dir->GetFile(i);

      if(fileString.find("snap") != std::string::npos)
        {
        this->GadgetFileNames.push_back(fileString);
        }
    }
 
// open only file0 and look what's inside
  hid_t    file_id = H5Fopen(this->GadgetFileNames[0].c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t    root_id, mesh_id, attr1;
  herr_t  status;
  root_id = H5Gopen(file_id, "/", H5P_DEFAULT);
  int j=0;

  H5E_auto_t func;
  void *client_data;
  H5Eget_auto2(H5E_DEFAULT, &func, &client_data);
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  for (auto const i: ParticleTypes)
    {
    if(H5Lexists(root_id, i.c_str(), H5P_DEFAULT))
      {
      mesh_id = H5Gopen(root_id, i.c_str(), H5P_DEFAULT);
      for( auto const s : this->FieldArrays[i])
        {
        if(H5Lexists(mesh_id, s.c_str(), H5P_DEFAULT))
          {
          this->PointDataArraySelection->AddArray((i + "/" + s).c_str());
          }
        }
      PartTypes[j] = true;
      H5Gclose(mesh_id);
      }
    j++;
    }
  H5Eset_auto2(H5E_DEFAULT, func, client_data);
  H5Gclose(root_id);

  root_id = H5Gopen(file_id, "/Header", H5P_DEFAULT);

  attr1 = H5Aopen_name(root_id, "NumPart_Total");
  if (H5Aread(attr1, H5T_NATIVE_INT, &this->NumPart_Total) < 0)
    vtkErrorMacro( << "cannot find the NumPart_Total");
  status = H5Aclose(attr1);

  H5Gclose(root_id);
  H5Fclose(file_id);
  dir->Delete();
  return 1;
}



int vtkGadgetReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkDataObject* doOutput = outInfo->Get(vtkDataObject::DATA_OBJECT());
#ifdef ALL_TYPES
  vtkMultiBlockDataSet* mb = vtkMultiBlockDataSet::SafeDownCast(doOutput);
  if (!mb)
    {
    return 0;
    }
#else
#ifdef OUTPUT_UG
  vtkUnstructuredGrid* output = vtkUnstructuredGrid::SafeDownCast(doOutput);
#else
  vtkPolyData* output = vtkPolyData::SafeDownCast(doOutput);
#endif
#endif

#ifdef PARAVIEW_USE_MPI
  if (this->Controller &&
      (this->UpdatePiece != this->Controller->GetLocalProcessId() ||
       this->UpdateNumPieces != this->Controller->GetNumberOfProcesses()))
  {
    vtkDebugMacro(<< "Parallel failure, Id's not right (ignore)");
    this->UpdatePiece = this->Controller->GetLocalProcessId();
    this->UpdateNumPieces = this->Controller->GetNumberOfProcesses();
  }
#else
  this->UpdatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
#endif

  int myType, load, MyNumber_of_Files;
  int nb_of_Files = this->GadgetFileNames.size();
  if(this->UpdateNumPieces == 1)
    {
    load = MyNumber_of_Files = nb_of_Files;
    }
  else
    {
    load = nb_of_Files / this->UpdateNumPieces;
    if (this->UpdatePiece < (this->UpdateNumPieces-1))
      {
      MyNumber_of_Files = load;
      }
    else
      {
      MyNumber_of_Files = nb_of_Files - (this->UpdateNumPieces-1) * load;
      }
    }

  hid_t   file_id, root_id, mesh_id, coords_id, dataset_id, filespace, attr1;
  herr_t  status;
  hsize_t dimsf[2];
  int LoadPart_Total[6] ={0,0,0,0,0,0};
  int lpT[6];
  for(int myFile = load*this->UpdatePiece; myFile< load*this->UpdatePiece + MyNumber_of_Files; myFile++)
    {
    file_id = H5Fopen(this->GadgetFileNames[myFile].c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    root_id = H5Gopen(file_id, "/Header", H5P_DEFAULT);

    attr1 = H5Aopen_name(root_id, "NumPart_ThisFile");
    if (H5Aread(attr1, H5T_NATIVE_INT, &lpT) < 0)
      vtkErrorMacro( << "cannot find the NumPart_ThisFile");

    for (int i=0; i< 6; i++)
      LoadPart_Total[i] += lpT[i];

    status = H5Aclose(attr1);

    H5Gclose(root_id);
    H5Fclose(file_id);
    }
#define PARALLEL_DEBUG 1
#ifdef PARALLEL_DEBUG
  std::ostringstream fname;
  fname << "/tmp/out." << this->UpdatePiece << ".txt" << ends;
  ofstream errs;
  errs.open(fname.str().c_str(),ios::app);
  errs << "piece " << this->UpdatePiece << " out of " << this->UpdateNumPieces << endl;

#endif

#ifdef ALL_TYPES
#ifdef OUTPUT_UG
  vtkUnstructuredGrid *output;
#else
  vtkPolyData *output;
#endif
#endif

  vtkFloatArray  *data;
  vtkIdTypeArray *uidata;

  int validPart;
  std::map<std::string, std::string> PVVarNames;
  PVVarNames["Velocities"] = "velocity";
  PVVarNames["Masses"] = "mass";
  PVVarNames["ParticleIDs"] = "id";
  PVVarNames["Density"] = "Density";
  PVVarNames["SmoothingLength"] = "SmoothingLength";

  std::map<std::string, std::string> HDF5VarNames;
  HDF5VarNames["velocity"] = "Velocities";
  HDF5VarNames["mass"] = "Masses";
  HDF5VarNames["id"] = "ParticleIDs";
  HDF5VarNames["Density"] = "Density";
  HDF5VarNames["SmoothingLength"] = "SmoothingLength";

  for(validPart=0, myType = Gas; myType<= Stars; myType++)
    {
    if(PartTypes[myType])
      {
#ifdef PARALLEL_DEBUG
      errs << " creating PolyData for PartType " << myType << " with " << LoadPart_Total[myType] << " points\n";
#endif

#ifdef ALL_TYPES
#ifdef OUTPUT_UG
      output = vtkUnstructuredGrid::New();
#else
      output = vtkPolyData::New();
#endif

      mb->SetBlock(validPart, output);
      mb->GetMetaData(validPart)->Set(vtkCompositeDataSet::NAME(), ParticleTypes[myType]);
      output->Delete();
#endif
      vtkFloatArray *coords = vtkFloatArray::New();
      coords->SetNumberOfComponents(3);
      coords->SetNumberOfTuples(LoadPart_Total[myType]);
      coords->SetName("coords");

      vtkPoints *points = vtkPoints::New();
      points->SetData(coords);
      output->SetPoints(points);
      coords->Delete();
      points->Delete();

      for(int i = 0 ; i < this->GetNumberOfPointArrays(); i++)
        {
        if(this->GetPointArrayStatus(this->GetPointArrayName(i)))
          {
          if(!strncmp(this->GetPointArrayName(i), ParticleTypes[myType].c_str(), 9)) // first 9 characters are equal
          {
          const char *name = &this->GetPointArrayName(i)[10];
#ifdef PARALLEL_DEBUG
      errs << " creating data array " << name << " for PartType " << myType << " with " << LoadPart_Total[myType] << " points\n";
#endif

          if(!strcmp(name, "Velocities"))
            {
            data = vtkFloatArray::New();
            data->SetNumberOfComponents(3);
            data->SetNumberOfTuples(LoadPart_Total[myType]);
            data->SetName(PVVarNames[name].c_str()); // use the mapped names; HDF5 reading will need the original
            output->GetPointData()->AddArray(data);
            data->Delete();
            }
          else if(!strcmp(name, "ParticleIDs"))
            {
            uidata = vtkIdTypeArray::New();
            uidata->SetNumberOfComponents(1);
            uidata->SetNumberOfTuples(LoadPart_Total[myType]);
            uidata->SetName(PVVarNames[name].c_str()); // use the mapped names; HDF5 reading will need the original
            output->GetPointData()->AddArray(uidata);
            uidata->Delete();
            }
          else
            {
            data = vtkFloatArray::New();
            data->SetNumberOfTuples(LoadPart_Total[myType]);
            data->SetName(PVVarNames[name].c_str()); // use the mapped names; HDF5 reading will need the original
            output->GetPointData()->AddArray(data);
            data->Delete();
            }
          }
          }
        }

      if (this->CellType == CellTypes::Vertex)
        {
        vtkCellArray *vertices =  vtkCellArray::New();
        vtkIdType* cells = vertices->WritePointer(LoadPart_Total[myType], 2 * LoadPart_Total[myType]);
        for (vtkIdType i = 0; i < LoadPart_Total[myType]; ++i)
          {
          cells[2 * i] = 1;
          cells[2 * i + 1] = i;
          }
       output->SetCells(VTK_VERTEX, vertices);
       vertices->Delete();
       }
     else if (this->CellType == CellTypes::PolyVertex)
        {
        vtkIdList *list = vtkIdList::New();
        list->SetNumberOfIds(LoadPart_Total[myType]);
        for(unsigned int i=0; i < LoadPart_Total[myType]; i++)
          list->SetId(i, i);
        output->Allocate(1);
        output->InsertNextCell(VTK_POLY_VERTEX, list);
        list->Delete();
        }
      validPart++;
      }
    }

  int offset, fileOffsetNodes[6] = {0,0,0,0,0,0};

  for(int myFile = load*this->UpdatePiece; myFile< load*this->UpdatePiece + MyNumber_of_Files; myFile++)
    {
    file_id = H5Fopen(this->GadgetFileNames[myFile].c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    root_id = H5Gopen(file_id, "/", H5P_DEFAULT);
    for(validPart=0, myType = Gas; myType<= Bndry; myType++)
      {
      if(PartTypes[myType])
        {
        char name[16];
        sprintf(name,"PartType%1d", myType);
#ifdef ALL_TYPES
#ifdef OUTPUT_UG
        output = static_cast<vtkUnstructuredGrid*>(mb->GetBlock(validPart));
#else
        output = static_cast<vtkPolyData*>(mb->GetBlock(validPart));
#endif
#endif
        mesh_id = H5Gopen(root_id, name, H5P_DEFAULT);

// insert coordinates read
        float *dptr = static_cast<vtkFloatArray *>(output->GetPoints()->GetData())->GetPointer(fileOffsetNodes[myType]*3);
        offset = ReadHDF5Dataset("Coordinates", mesh_id, dptr);
// end of coordinates read

// insert PointData here
      for(int i = 0 ; i < this->GetNumberOfPointArrays(); i++)
        {
        if(this->GetPointArrayStatus(this->GetPointArrayName(i)))
          {
          if(!strncmp(this->GetPointArrayName(i), ParticleTypes[myType].c_str(), 9)) // first 9 characters are equal
            {
            const char *name = PVVarNames[&this->GetPointArrayName(i)[10]].c_str();
#ifdef PARALLEL_DEBUG
      errs << this->GadgetFileNames[myFile] << ": reading data array " << name << " for PartType " << myType << " with " << this->NumPart_Total[myType] << " points\n";
#endif
          if(!strcmp(name, "id"))
            {
            uidata = static_cast<vtkIdTypeArray *>(output->GetPointData()->GetArray(name));
            vtkTypeInt64 *lptr = uidata->GetPointer( data->GetNumberOfComponents() * fileOffsetNodes[myType] );
            ReadHDF5INT64Dataset(HDF5VarNames[name].c_str(), mesh_id, lptr);
            }
          else
            {
            data = static_cast<vtkFloatArray *>(output->GetPointData()->GetArray(name));
            dptr = data->GetPointer( data->GetNumberOfComponents() * fileOffsetNodes[myType] );
            ReadHDF5Dataset(HDF5VarNames[name].c_str(), mesh_id, dptr);
            // split by component
            if(!strcmp(name, "velocity"))
              {
              double tuple[3];
              int NbTuples = data->GetNumberOfTuples();
              vtkFloatArray *vx = vtkFloatArray::New();
              vx->SetName("vx");
              vx->SetNumberOfTuples(NbTuples);
              output->GetPointData()->AddArray(vx);
              vx->Delete();

              vtkFloatArray *vy = vtkFloatArray::New();
              vy->SetName("vy");
              vy->SetNumberOfTuples(NbTuples);
              output->GetPointData()->AddArray(vy);
              vy->Delete();

              vtkFloatArray *vz = vtkFloatArray::New();
              vz->SetName("vz");
              vz->SetNumberOfTuples(NbTuples);
              output->GetPointData()->AddArray(vz);
              vz->Delete();
              for(vtkIdType i=0; i < NbTuples; i++)
                {
                data->GetTuple(i, tuple);
                vx->SetTuple1(i, tuple[0]);
                vy->SetTuple1(i, tuple[1]);
                vz->SetTuple1(i, tuple[2]);
                }
              }
            }
            }
          }
        }
// end of PointData read

        H5Gclose(mesh_id);
        validPart++;
        }
      fileOffsetNodes[myType] += offset;
      } // for all part types
    H5Gclose(root_id);
    H5Fclose(file_id);
    } // for all files in my load
  this->CloseFile();

#ifdef PARALLEL_DEBUG
  errs.close();
#endif
  return 1;
}


//----------------------------------------------------------------------------
const char* vtkGadgetReader::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}
//----------------------------------------------------------------------------
int vtkGadgetReader::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}
//----------------------------------------------------------------------------
void vtkGadgetReader::SetPointArrayStatus(const char* name, int status)
{
  if (status!=this->GetPointArrayStatus(name))
    {
    if (status)
      {
      this->PointDataArraySelection->EnableArray(name);
      }
    else
      {
      this->PointDataArraySelection->DisableArray(name);
      }
    this->Modified();
    }
}

void vtkGadgetReader::EnableAllPointArrays()
{
    this->PointDataArraySelection->EnableAllArrays();
}

//----------------------------------------------------------------------------
void vtkGadgetReader::DisableAllPointArrays()
{
    this->PointDataArraySelection->DisableAllArrays();
}

//----------------------------------------------------------------------------
int vtkGadgetReader::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}
//----------------------------------------------------------------------------
void vtkGadgetReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Directory Name: " <<
    (this->DirectoryName ? this->DirectoryName : "(none)") << "\n";
}
