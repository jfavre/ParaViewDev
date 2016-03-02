#include "vtkDataArraySelection.h"
#include "vtkDataSetAttributes.h"
#include "vtkDoubleArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkErrorCode.h"
#include "vtkExtentTranslator.h"
#include "vtkFieldData.h"
#include "vtkMath.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <math.h>
#include <set>
#include "vtkETHZGeoPhysicsHDF5SphericalReader.h"
#include "BaseGeoPhysicsHDF5SphericalReader.h"

#define DATA_LENGTH H5T_NATIVE_FLOAT

//vtkCxxRevisionMacro(vtkETHZGeoPhysicsHDF5SphericalReader, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkETHZGeoPhysicsHDF5SphericalReader);

vtkETHZGeoPhysicsHDF5SphericalReader::vtkETHZGeoPhysicsHDF5SphericalReader()
{
  this->FileName            = NULL; 
  this->DebugOff();
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
  this->PointDataArraySelection = vtkDataArraySelection::New();
  this->Radius = NULL;
  this->Phi = NULL;
  this->Theta = NULL;
  this->ghostLevel = 0;
}

vtkETHZGeoPhysicsHDF5SphericalReader::~vtkETHZGeoPhysicsHDF5SphericalReader()
{
  vtkDebugMacro(<< "cleaning up inside destructor");
  if (this->FileName)
    delete [] this->FileName;
  this->PointDataArraySelection->Delete();
 
  if (this->Radius)
    {
    //cerr << "this->Radius->Delete()\n";
    this->Radius->Delete();
    this->Radius = NULL;
    }
  if (this->Phi)
    {
    this->Phi->Delete();
    this->Phi = NULL;
    }
  if (this->Theta)
    {
    this->Theta->Delete();
    this->Theta = NULL;
    }
}

/*
int vtkETHZGeoPhysicsHDF5SphericalReader::FillOutputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkStructuredGrid");

  return 1;
}
*/

std::string BaseName;
std::set<std::string> Vector3GroupNames;
std::set<std::string> ScalarNames;

herr_t file_info2(hid_t loc_id, const char *name, void *opdata)
{
    H5G_stat_t statbuf;

    /*
     * Get type of the object and display its name and type.
     * The name of the object is passed to this function by
     * the Library. Some magic :-)
     */
    H5Gget_objinfo(loc_id, name, 0, &statbuf);
    switch (statbuf.type) {
    case H5G_GROUP:
         //cerr << "found group " << name << endl;
         Vector3GroupNames.insert(name);
         break;
    case H5G_DATASET:
         ScalarNames.insert(name);
         //cerr << "   found dataset " << name << endl;
         break;
    case H5G_TYPE:
         break;
    default:
         break;
    }
    return 0;
 }


int vtkETHZGeoPhysicsHDF5SphericalReader::RequestInformation(vtkInformation* request,
                         vtkInformationVector** inputVector,
                         vtkInformationVector* outputVector)
{
  hid_t err_id, filespace, dataset_id, file_id, root_id, geo_id, basegroup_id;
  int my_extents[6];
  hsize_t dimsf[3];
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  //outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);
//cerr << ".........H5Fopen(" << this->FileName << ").............."<< endl;
  file_id = H5Fopen(this->FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  root_id = H5Gopen(file_id, "/", H5P_DEFAULT);
  // before opening, we need to inquire if object really exists and turn off errors.
  Vector3GroupNames.clear();
  ScalarNames.clear();
  H5Giterate(file_id, "/", NULL, file_info2, NULL);
  //cerr << ".........opening the base group.............."<< endl;
  if(Vector3GroupNames.size() == 1)
    {
    BaseName = *Vector3GroupNames.begin();
    basegroup_id = H5Gopen(root_id, BaseName.c_str(), H5P_DEFAULT);
    H5Giterate(root_id, BaseName.c_str(), NULL, file_info2, NULL);
    }

  std::set<std::string>::iterator varit;
  for(varit = Vector3GroupNames.begin(); varit != Vector3GroupNames.end(); varit++)
    {
    if(*varit != BaseName)
      {
      //cerr << "........opening the 2nd level group... "<< *varit << endl;
      if(*varit != "grid")
        {
        this->PointDataArraySelection->AddArray(varit->c_str());
        //cerr << "add vector label " <<  *varit << endl;
        }
      }
    }
  for(varit = ScalarNames.begin(); varit != ScalarNames.end(); varit++)
    { // add any additional scalars which were found at the same level as the vector Groups.
    this->PointDataArraySelection->AddArray(varit->c_str());
    //cerr << "add scalar label " <<  *varit << endl;
    }
  herr_t error;
  H5E_auto_t func;
  void *client_data;
  H5Eget_auto2(H5E_DEFAULT, &func, &client_data);
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  std::string  var;

  this->Radius = vtkFloatArray::New();
  float temp, *p;
  //cerr << "this->Radius = vtkFloatArray::New()\n";
  this->Radius->SetName("R");
  this->Radius->SetNumberOfComponents(1);

  dataset_id = H5Dopen(basegroup_id, "grid/r axis", H5P_DEFAULT);
  if((dataset_id  < 0 ) && (dataset_id = H5Dopen(root_id, "raxis", H5P_DEFAULT)) < 0)  // does not exist
    {
    H5Fclose(file_id);
    return 0;
    }
  else
    {
    filespace = H5Dget_space(dataset_id);
    H5Sget_simple_extent_dims(filespace, dimsf, NULL);
    this->Radius->SetNumberOfTuples(dimsf[0]);
    H5Dread(dataset_id, DATA_LENGTH, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->Radius->GetPointer(0));
    H5Sclose(filespace);
    H5Dclose(dataset_id);
    my_extents[4] = 0;
    my_extents[5] = dimsf[0]-1;
/*
    p = this->Radius->GetPointer(0);
    for (int i=0; i < dimsf[0]/2 ; i++)
      {
      temp = p[i];
      p[i] = p[dimsf[0]-1-i];
      p[dimsf[0]-1-i] = temp;
      }
*/
    }

  this->Phi = vtkFloatArray::New();
  this->Phi->SetName("phi");
  this->Phi->SetNumberOfComponents(1);
  dataset_id = H5Dopen(basegroup_id, "grid/p axis", H5P_DEFAULT);
  if((dataset_id < 0) && (dataset_id = H5Dopen(root_id, "paxis", H5P_DEFAULT)) < 0)  // does not exist
    {
    H5Fclose(file_id);
    return 0;
    }
  else
    {
    filespace = H5Dget_space(dataset_id);
    H5Sget_simple_extent_dims(filespace, dimsf, NULL);
    this->Phi->SetNumberOfTuples(dimsf[0]);
    H5Dread(dataset_id, DATA_LENGTH, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->Phi->GetPointer(0));
    H5Sclose(filespace);
    H5Dclose(dataset_id);
    my_extents[0] = 0;
    my_extents[1] = dimsf[0]-1;
    // construct now the two arrays sin(phi) and cos(phi) which will be used a few of times
    this->SinPhi = new double[dimsf[0]];
    this->CosPhi = new double[dimsf[0]];
    for(int i=0; i < dimsf[0]; i++)
      {
      this->SinPhi[i] = sin(this->Phi->GetValue(i) - M_PI_2);
      this->CosPhi[i] = cos(this->Phi->GetValue(i) - M_PI_2);
      }
    }

  this->Theta = vtkFloatArray::New();
  this->Theta->SetName("theta");
  this->Theta->SetNumberOfComponents(1);

  dataset_id = H5Dopen(basegroup_id, "grid/t axis", H5P_DEFAULT);
  if((dataset_id < 0) && (dataset_id = H5Dopen(root_id, "taxis", H5P_DEFAULT)) < 0)  // does not exist
    {
    H5Fclose(file_id);
    return 0;
    }
  else
    {
    filespace = H5Dget_space(dataset_id);
    H5Sget_simple_extent_dims(filespace, dimsf, NULL);
    this->Theta->SetNumberOfTuples(dimsf[0]);
    H5Dread(dataset_id, DATA_LENGTH, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->Theta->GetPointer(0));
    //this->Theta->SetValue(0, 0.0); this->Theta->SetValue(1, 1.0); this->Theta->SetValue(2, 2.0);
    H5Sclose(filespace);
    H5Dclose(dataset_id);
    my_extents[2] = 0;
    my_extents[3] = dimsf[0]-1;

    this->SinTheta = new double[dimsf[0]];
    this->CosTheta = new double[dimsf[0]];
    for(int i=0; i < dimsf[0]; i++)
      {
      this->SinTheta[i] = sin(this->Theta->GetValue(i));
      this->CosTheta[i] = cos(this->Theta->GetValue(i));
      }
    }
  //cerr << "creating cos() and sin() arrays\n";
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), my_extents, 6);
#if VTK_MINOR_VERSION == 2
    outInfo->Set(CAN_PRODUCE_SUB_EXTENT(), 1);
#endif
  if(outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS()))
    {
    ghostLevel = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());
    }
/*
  cerr << "Set WHOLE_EXTENT  extents: " << my_extents[0] << ", " <<
                                                    my_extents[1] << ", " <<
                                                    my_extents[2] << ", " <<
                                                    my_extents[3] << ", " <<
                                                    my_extents[4] << ", " <<
                                                    my_extents[5] << endl;
*/
  H5Eset_auto2(H5E_DEFAULT, func, client_data);
  H5Gclose(basegroup_id);
  H5Gclose(root_id);
  H5Fclose(file_id);
  vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_FLOAT, 1);
// every pvserver process executes requestInformation but there is no notion of total # of processes
  return 1;
}

int vtkETHZGeoPhysicsHDF5SphericalReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  hid_t    mesh_id, dataset_id, dataset_id2, parameters_id, memspace;
  hid_t    aid1, attr1, file_id, root_id, dataspace_id, group_id, basegroup_id;
  herr_t   status;
  int ext[6];
  vtkDebugMacro( << "RequestData(BEGIN)");

  std::set<std::string>::iterator varit;
  std::vector<std::string>::iterator compit;
    // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  std::vector<std::string> comp;
  comp.clear();
  comp.push_back("p");
  comp.push_back("t");
  comp.push_back("r");

  this->UpdateProgress(0);
  int done=0;
  int tasks = 0;
  if(this->PointDataArraySelection->ArrayIsEnabled("vp"))
    tasks++;

  vtkStructuredGrid *output = vtkStructuredGrid::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), ext);
  output->SetExtent(ext);
  unsigned int updatePiece = 0;
  unsigned int updateNumPieces = 1;

    if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()) &&
      outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES()))
    {
    updatePiece = static_cast<unsigned int>(
        outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()));
    updateNumPieces =  static_cast<unsigned int>(
        outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES()));
    }

   if(outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS()))
    {
    this->ghostLevel = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());
    }

  if(!updatePiece)
  cerr << "GL = " << this->ghostLevel<< "\n" << "Set EXTENT  extents: " << ext[0] << ", " <<
                                                    ext[1] << ", " <<
                                                     ext[2] << ", " <<
                                                    ext[3] << ", " <<
                                                     ext[4] << ", " <<
                                                    ext[5] << endl;


  if( this->ghostLevel == 368 )
    {
    vtkUnsignedCharArray *cellGhostArray = vtkUnsignedCharArray::New();
    cellGhostArray->SetName("vtkGhostLevels");
    int size = (ext[1]-ext[0])*(ext[3]-ext[2])*(ext[5]-ext[4]);
    cerr << "allocate ghost cell array of size " << size << endl;
    cellGhostArray->SetNumberOfTuples(size);
    unsigned char *p = static_cast<unsigned char*>(cellGhostArray->GetVoidPointer(0));
    for(int i=0; i < size; i++){
      *p++ = 0;
    }
    //output->GetCellData()->AddArray(cellGhostArray);
    //cellGhostArray->Delete();
    }
      // Create the grid points from the grid image.
  vtkPoints *points = ConstructPoints(ext, this->SinTheta, this->SinPhi, this->CosTheta, this->CosPhi,
                                      this->Radius->GetPointer(0),
                                      this->Theta->GetPointer(0),
                                      this->Phi->GetPointer(0),
                                      true);
  output->SetPoints(points);
  points->Delete();
  points = NULL;

  if (!this->FileName)
    {
    vtkErrorMacro(<< "error reading header specified!");
    return 0;
    }

  //cerr << "H5Fopen " << this->FileName << "\n";
  file_id = H5Fopen(this->FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  root_id = H5Gopen(file_id, "/", H5P_DEFAULT);
  basegroup_id = H5Gopen(root_id, BaseName.c_str(), H5P_DEFAULT);

// read all 3-vector first
  for(varit = Vector3GroupNames.begin(); varit != Vector3GroupNames.end(); varit++)
    {
    const char *varname = varit->c_str();
    if(this->PointDataArraySelection->ArrayIsEnabled(varname))
      {
      bool NoErrorReadingAll3Comps = true;
      group_id = H5Gopen(basegroup_id, varname, H5P_DEFAULT);
      for(compit = comp.begin(); compit != comp.end(); compit++)
        {
	const char *comp_varname = strdup((varname + *compit).c_str());
        if((dataset_id = H5Dopen(group_id, comp_varname, H5P_DEFAULT)) < 0)
          {
	  cerr << "error opening dataset " << comp_varname << " for file " << this->FileName << endl;
	  NoErrorReadingAll3Comps = false;
	  }
        else
          {
          vtkFloatArray *tdata = ReadScalarVar(&group_id, ext, comp_varname);
          if(tdata)
            {
            output->GetPointData()->AddArray(tdata);
            tdata->Delete();
            }
          else
	    NoErrorReadingAll3Comps = false;
          H5Dclose(dataset_id);
          }
          free((void*)comp_varname);
        }// for all 3 components
      H5Gclose(group_id);
      if(NoErrorReadingAll3Comps)
        {
	//cout << "Composing the 3-vector " << varname << " after spherical-to-cartesian mapping\n";
	vtkFloatArray *fa = ReadVector(&basegroup_id, ext, varname,
                                       this->SinTheta, this->SinPhi, this->CosTheta, this->CosPhi,
                                       true);

        output->GetPointData()->AddArray(fa);
        fa->Delete();
	}
      }
    }

// finish reading all left-over scalars

  for(varit = ScalarNames.begin(); varit != ScalarNames.end(); varit++)
    {
    const char *varname = varit->c_str();
    if(this->PointDataArraySelection->ArrayIsEnabled(varname))
      {
      vtkFloatArray *tdata = ReadScalarVar(&basegroup_id, ext, varname);
      if(tdata)
        {
        output->GetPointData()->AddArray(tdata);
        tdata->Delete();
        }
      else cerr << "error reading scalar " << varname << endl;
      }
    }

  H5Gclose(basegroup_id);
  H5Gclose(root_id);
  H5Fclose(file_id);
  done++;
  output->GetFieldData()->AddArray(this->Radius); //this->Radius->Delete();
  output->GetFieldData()->AddArray(this->Phi);    //this->Phi->Delete();
  output->GetFieldData()->AddArray(this->Theta);  //this->Theta->Delete();
  //this->UpdateProgress((double)done/tasks);

  //cerr << "deleting cos() and sin() arrays\n";
  delete [] this->SinPhi;
  delete [] this->CosPhi;
  delete [] this->SinTheta;
  delete [] this->CosTheta;

  vtkDebugMacro( << "RequestData(END)");
  //cerr << "vtkETHZGeoPhysicsHDF5Reader::RequestData(END)\n";

  this->UpdateProgress(1);
  return 1;
}



void vtkETHZGeoPhysicsHDF5SphericalReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Reader:\n";
  os << indent << indent << "Filename      : " << this->FileName << endl;
}

void vtkETHZGeoPhysicsHDF5SphericalReader::EnablePointArray(const char* name)
{
  this->SetPointArrayStatus(name, 1);
}

void vtkETHZGeoPhysicsHDF5SphericalReader::DisablePointArray(const char* name)
{
  this->SetPointArrayStatus(name, 0);
}

void vtkETHZGeoPhysicsHDF5SphericalReader::EnableAllPointArrays()
{
  this->PointDataArraySelection->EnableAllArrays();
}

void vtkETHZGeoPhysicsHDF5SphericalReader::DisableAllPointArrays()
{
  this->PointDataArraySelection->DisableAllArrays();
}

const char* vtkETHZGeoPhysicsHDF5SphericalReader::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}

int vtkETHZGeoPhysicsHDF5SphericalReader::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}

void vtkETHZGeoPhysicsHDF5SphericalReader::SetPointArrayStatus(const char* name, int status)
{
  if(status)
    {
    this->PointDataArraySelection->EnableArray(name);
    }
  else
    {
    this->PointDataArraySelection->DisableArray(name);
    }
}

int vtkETHZGeoPhysicsHDF5SphericalReader::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}

#ifdef CELL
void vtkETHZGeoPhysicsHDF5Reader::EnableCellArray(const char* name)
{
  this->SetCellArrayStatus(name, 1);
}

void vtkETHZGeoPhysicsHDF5Reader::DisableCellArray(const char* name)
{
  this->SetCellArrayStatus(name, 0);
}

void vtkETHZGeoPhysicsHDF5Reader::EnableAllCellArrays()
{
  this->CellDataArraySelection->EnableAllArrays();
}

void vtkETHZGeoPhysicsHDF5Reader::DisableAllCellArrays()
{
  this->CellDataArraySelection->DisableAllArrays();
}

const char* vtkETHZGeoPhysicsHDF5Reader::GetCellArrayName(int index)
{
  return this->CellDataArraySelection->GetArrayName(index);
}

int vtkETHZGeoPhysicsHDF5Reader::GetCellArrayStatus(const char* name)
{
  return this->CellDataArraySelection->ArrayIsEnabled(name);
}

void vtkETHZGeoPhysicsHDF5Reader::SetCellArrayStatus(const char* name, int status)
{
  if(status)
    {
    this->CellDataArraySelection->EnableArray(name);
    }
  else
    {
    this->CellDataArraySelection->DisableArray(name);
    }
}

int vtkETHZGeoPhysicsHDF5Reader::GetNumberOfCellArrays()
{
  return this->CellDataArraySelection->GetNumberOfArrays();
}
#endif
