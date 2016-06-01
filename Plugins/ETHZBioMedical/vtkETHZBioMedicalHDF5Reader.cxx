// added reading the "Emoduli" contained in "Mesh"
// there was existing support for that, but it might break reading the older files.
// Thu Dec 19 13:43:44 CET 2013 
#include "vtkETHZBioMedicalHDF5Reader.h"
#include "BaseETHZBioMedicalHDF5Reader.h"

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
#include "vtkMultiBlockDataSet.h"
#include "vtkSelection.h"
#include "vtkSelectionNode.h"
#include "vtkExtractSelection.h"
#include "vtkXMLDataSetWriter.h"

#include <sys/time.h>
#include <algorithm>
#include <vector>
using namespace std;

#include "vtkPVConfig.h"
#ifdef PARAVIEW_USE_MPI
  #include "vtkMPI.h"
  #include "vtkMPIController.h"
  #include "vtkMPICommunicator.h"
#endif

#include <set>
#include <sstream>
hid_t    file_id;
hsize_t  array_dim[2];

herr_t file_info(hid_t loc_id, const char *name, void *opdata)
{
    H5G_stat_t statbuf;
    hsize_t dimsf[2];
    hid_t filespace;
    hid_t root_id, mesh_id, dataset_id;
    int N, *MeshSizes = (int *)opdata;
    /*
     * Get type of the object and display its name and type.
     * The name of the object is passed to this function by 
     * the Library. Some magic :-)
     */

    H5Gget_objinfo(loc_id, name, 0, &statbuf);
    switch (statbuf.type) {
    case H5G_GROUP: 
         printf(" Object with name %s is a group\n", name);
         break;
    case H5G_DATASET:

      root_id = H5Gopen(loc_id, "/", H5P_DEFAULT);

      mesh_id = H5Gopen(root_id, "Mesh", H5P_DEFAULT);

      dataset_id = H5Dopen(mesh_id, "Elements", H5P_DEFAULT);

      filespace = H5Dget_space(dataset_id);
      H5Sget_simple_extent_dims(filespace, dimsf, NULL);
      N = dimsf[0];
      H5Sclose(filespace);
 /*
         if (N == MeshSizes[0])
           printf("dataset for nodes");
         else if(N == MeshSizes[1])
           printf("dataset for cells");
         else
           printf("bad size on dataset");
*/
     H5Dclose(dataset_id);
      H5Gclose(mesh_id);
      H5Gclose(root_id);
         printf(" Object with name %s is a dataset of size %d\n", name, N);

         break;
    case H5G_TYPE: 
         printf(" Object with name %s is a named datatype\n", name);
         break;
    default:
         printf(" Unable to identify an object ");
    }
    return 0;
 }

vtkStandardNewMacro(vtkETHZBioMedicalHDF5Reader);

void vtkETHZBioMedicalHDF5Reader::Find_varnames(hid_t *loc_id)
{
  if(H5Lexists(*loc_id, "Nodal displacements", H5P_DEFAULT) || H5Lexists(*loc_id, "disp", H5P_DEFAULT))
    {
    this->PointDataArraySelection->AddArray("Nodal displacements");
    }
  if(H5Lexists(*loc_id, "Nodal forces", H5P_DEFAULT) || H5Lexists(*loc_id, "force", H5P_DEFAULT))
    {
    this->PointDataArraySelection->AddArray("Nodal forces");
    }
  if(H5Lexists(*loc_id, "Element strain", H5P_DEFAULT))
    {
    for(int i=0; i < 8; i++)
      {
      if(!this->CellDataArraySelection->AddArray(strain_var_names[i]))
        {
        //cerr << "nothing done. Variable already existed\n";
        }
      }
    }
  if(H5Lexists(*loc_id, "Element stress", H5P_DEFAULT))
    {
    for(int i=0; i < 7; i++)
      {
      if(!this->CellDataArraySelection->AddArray(stress_var_names[i]))
        {
        //cerr << "nothing done. Variable already existed\n";
        }
      }
    }

  const char* args[] ={"EFF","SED","DP", "DP_s", "DP_e","e_xx", "e_xy", "e_xz", "e_yy", "e_yz", "e_zz", "VonMises", "Emoduli"};
  int SizeOfargs = 13;
  for(int i=0; i < SizeOfargs; i++)
    {
    if(H5Lexists(*loc_id, args[i], H5P_DEFAULT))
      {
      this->CellDataArraySelection->AddArray(args[i]);
      }
    }
}

int vtkETHZBioMedicalHDF5Reader::RequestInformation(vtkInformation* request,
                         vtkInformationVector** inputVector,
                         vtkInformationVector* outputVector)
{
  file_id = H5Fopen(this->FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t g_id, root_id, mesh_id, dataset_id;
  hid_t filespace;
  hsize_t dimsf[2];
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(CAN_HANDLE_PIECE_REQUEST(), 1);
  //outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);

  root_id = H5Gopen(file_id, "/", H5P_DEFAULT);

  mesh_id = H5Gopen(root_id, "Mesh", H5P_DEFAULT);
  dataset_id = H5Dopen(mesh_id, "Elements", H5P_DEFAULT);
  filespace = H5Dget_space(dataset_id);
  H5Sget_simple_extent_dims(filespace, dimsf, NULL);
  this->NbCells = dimsf[0];
  H5Sclose(filespace);
  H5Dclose(dataset_id);

  dataset_id = H5Dopen(mesh_id, "Coordinates", H5P_DEFAULT);
  filespace = H5Dget_space(dataset_id);
  H5Sget_simple_extent_dims(filespace, dimsf, NULL);
  this->NbNodes = dimsf[0];
  H5Sclose(filespace);
  H5Dclose(dataset_id);
  H5Gclose(mesh_id);

  int MeshSizes[2] = {this->NbNodes, this->NbCells};
  if(H5Lexists(root_id, "Solution", H5P_DEFAULT))
    {
    hid_t solutions_id = H5Gopen(root_id, "Solution", H5P_DEFAULT);
    if(solutions_id >= 0)
      {
      this->Find_varnames(&solutions_id);
      H5Giterate(file_id, "/Solution", NULL, file_info, MeshSizes);

  //Get Number of time steps and the corresponding scale factors
  m_NumLoadSteps = 0;
  int status;

  if(H5Lexists(solutions_id, "#Load Steps", H5P_DEFAULT))
    {
    if(H5Lexists(solutions_id, "Load step 1", H5P_DEFAULT))
      {
      g_id = H5Gopen(solutions_id, "Load step 1", H5P_DEFAULT);
      this->Find_varnames(&g_id);
      H5Gclose(g_id);
      }
    hid_t dataset_id = H5Dopen(solutions_id, "#Load Steps", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_NumLoadSteps);
    H5Dclose(dataset_id);
    vtkDebugMacro( << "m_NumLoadSteps = " << m_NumLoadSteps);

    m_NumStateVariables = 0;
    char tmp_string[256];
    if (m_NumLoadSteps)
      {
      dataset_id = H5Dopen(solutions_id, "#State Variables", H5P_DEFAULT);
      status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_NumStateVariables);
      H5Dclose(dataset_id);
      for(int i=0; i < m_NumStateVariables; i++)
        {
	sprintf(tmp_string, "State Variable %d", i+1);
	if(!this->CellDataArraySelection->AddArray(tmp_string))
	  {
	    //cerr << "nothing done. Variable already existed\n";
	  }
        }
      }
  
    m_LoadStepScalingFactors = new double[m_NumLoadSteps];
    for(int i=0; i<m_NumLoadSteps; i++)
      {
      sprintf(tmp_string, "Load step %d/Boundary condition scaling factor", i+1);
      dataset_id = H5Dopen(solutions_id, tmp_string, H5P_DEFAULT);
      status = H5Dread(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
                       H5P_DEFAULT, &m_LoadStepScalingFactors[i]);
      H5Dclose(dataset_id);
      }

    double timeRange[2];
    timeRange[0]=m_LoadStepScalingFactors[0];
    timeRange[1]=m_LoadStepScalingFactors[m_NumLoadSteps-1];

    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
    for(int i=0;i<m_NumLoadSteps;i++)
      {    
      outInfo->Append(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), m_LoadStepScalingFactors[i]);
      }

      }
    H5Gclose(solutions_id);
    }
  }
  if(H5Lexists(root_id, "Mesh", H5P_DEFAULT))
    {
    hid_t mesh_id = H5Gopen(root_id, "Mesh", H5P_DEFAULT);
    if(mesh_id >= 0)
      {
      this->Find_varnames(&mesh_id);
      H5Giterate(file_id, "/Mesh", NULL, file_info, NULL);
      }
    H5Gclose(mesh_id);
    }

  H5Gclose(root_id);
  return 1;
}



int vtkETHZBioMedicalHDF5Reader::RequestData(
                vtkInformation* vtkNotUsed(request),
                vtkInformationVector** vtkNotUsed(inputVector),
                vtkInformationVector* outputVector)
{
  hid_t root_id, mesh_id, dataset_id, coords_id, solutions_id=0, dataset_id2, parameters_id, ids_id;
  hid_t filespace;
  hsize_t dimsf[2];
  herr_t   status;
  bool Materials_Found = false;
  int I;
  struct timeval tv0, tv1, res;

  vtkDebugMacro( << "RequestData(BEGIN)");
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkDataObject* doOutput = outInfo->Get(vtkDataObject::DATA_OBJECT());
  vtkMultiBlockDataSet* mb = vtkMultiBlockDataSet::SafeDownCast(doOutput);
  if (!mb)
    {
    return 0;
    }
//used temporarily, until we split with vtkSelection
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::New();
  //int rank, num_proc;
  //MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int piece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  int numPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
#define PARALLEL_DEBUG 1
#ifdef PARALLEL_DEBUG
  std::ostringstream fname;
  fname << "/tmp/out." << piece << ".txt" << ends;
  ofstream errs;
  errs.open(fname.str().c_str(),ios::app);
  errs << "piece " << piece << " out of " << numPieces << endl;
#endif

  int actualTimeStep = 0;
  double requestedTimeValue = 0.0;
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
    {
    requestedTimeValue = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

    double eps = 0.0001 * requestedTimeValue;
    for (int i=0; i<m_NumLoadSteps; i++)
      {
      if (fabs(requestedTimeValue - m_LoadStepScalingFactors[i]) < eps)
        {
	actualTimeStep = i+1;
	//cout << "actualTimeStep = " << actualTimeStep << endl;
	break;
        }
      }
    mb->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), requestedTimeValue);
    }
  //cout << "requestedTimeValue "  << requestedTimeValue << endl;
  //cout << "actualTimeStep "  << actualTimeStep << endl;

  vtkDebugMacro(<< "getting piece " << piece << " out of " << numPieces << " pieces");
  if (!this->FileName)
    {
    vtkErrorMacro(<< "error reading header specified!");
    return 0;
    }
  size_t size;
  //hid_t file_plist_id = H5Pcreate(H5P_FILE_ACCESS);
  //H5Pset_sieve_buf_size (file_plist_id, 64*1024);
  //file_id = H5Fopen(this->FileName, H5F_ACC_RDONLY, file_plist_id);
  file_id = H5Fopen(this->FileName, H5F_ACC_RDONLY, H5P_DEFAULT);

  //H5Pget_sieve_buf_size (file_plist_id, &size);
  //cerr << "default sieve size = " << size << endl;
  //H5Pclose(file_plist_id);

  root_id = H5Gopen(file_id, "/", H5P_DEFAULT);
  if(H5Lexists(root_id, "/Parameters", H5P_DEFAULT))
    {
    ParaSol = false;
/*
    dataset_id = H5Dopen(parameters_id, "#nodes", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->NbNodes);
    vtkDebugMacro( << "this->NbNodes = " << this->NbNodes);
    H5Dclose(dataset_id);

    dataset_id2 = H5Dopen(parameters_id, "#elements", H5P_DEFAULT);
    status = H5Dread(dataset_id2, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->NbCells);
    vtkDebugMacro( << "this->NbCells = " << this->NbCells);
    H5Dclose(dataset_id2);
*/
    }
  else
    {
    ParaSol = true;
    }

  mesh_id = H5Gopen(root_id, "/Mesh", H5P_DEFAULT);
  dataset_id = H5Dopen(mesh_id, "Elements", H5P_DEFAULT);
  coords_id = H5Dopen(mesh_id, "Coordinates", H5P_DEFAULT);
//printf("this->NbNodes = %d, this->NbCells = %d\n", this->NbNodes, this->NbCells);

// here we allocate the final list necessary for VTK. It includes an extra
// integer for every hexahedra to say that the next cell contains 8 nodes.
// to avoid allocating two lists and making transfers from one to the other
// we pre-allocate for 9 * NbCells.
// In the original code, we then we moved things in-memory, starting
// from the end of the list, so as to not overwrite things
// We now use HDF5's memory space to conveniently place the Nx8 array read from disk 
// in the Nx9 matrix, leaving the left-most column empty to add the value 8.
  long MyNumber_of_Cells, MyNumber_of_Nodes;
  long load;
  vtkIdType *vtkfinal_id_list, minId, maxId, *destptr;
  hsize_t count[2], offset[2];

  hid_t memspace, dataspace;

  if(numPieces == 1)
    {
    load = MyNumber_of_Cells = this->NbCells;
    }
  else
    {
    load = this->NbCells / numPieces;
    if (piece < (numPieces-1))
      {
      MyNumber_of_Cells = load;
      }
    else
      {
      MyNumber_of_Cells = this->NbCells - (numPieces-1) * load;
      }
    }
  vtkIdTypeArray *vtklistcells = vtkIdTypeArray::New();
  vtklistcells->SetNumberOfValues(MyNumber_of_Cells * (8 + 1));
  vtkfinal_id_list = vtklistcells->GetPointer(0);

  
#ifdef PARALLEL_DEBUG
  errs << "CPU " << piece << " aLLocating " << MyNumber_of_Cells  << " list of elements of size " << "*9*" << sizeof(vtkIdType) << " bytes = " << MyNumber_of_Cells * 9 *sizeof(vtkIdType) << " bytes\n";
#endif
  count[0] = MyNumber_of_Cells;
  count[1] = 8 + 1;
  memspace = H5Screate_simple(2, count, NULL);

  offset[0] = 0;
  offset[1] = 1;
  count[0] = MyNumber_of_Cells;
  count[1] = 8;
  H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset, NULL, count, NULL);
  hsize_t     dims_out[2];
  dataspace = H5Dget_space(dataset_id);
  int ndims  = H5Sget_simple_extent_ndims(dataspace);
  H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
  //printf("ndims %d, dimensions %lu x %lu \n", ndims, (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]));
  count[0] = MyNumber_of_Cells;
  count[1] = 8;
  offset[0] = piece * load;
  offset[1] = 0;
  //cerr <<  "CPU " << piece << ": count = " << count[0] << " x " <<  count[1] << endl;
  //cerr <<  "CPU " << piece << ": offset = " << offset[0] << " x " <<  offset[1] << endl;
  dataspace = H5Dget_space(dataset_id);
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

  hid_t plist_xfer = H5Pcreate(H5P_DATASET_XFER);
  status = H5Pset_buffer(plist_xfer, (hsize_t)1*1024*1024, NULL, NULL);
  // this is the default value
  if (sizeof(vtkIdType) == H5Tget_size(H5T_NATIVE_INT))
    {
    H5Dread(dataset_id, H5T_NATIVE_INT, memspace, dataspace,
            //H5P_DEFAULT,
            plist_xfer,
            vtkfinal_id_list);
    }
  else if (sizeof(vtkIdType) == H5Tget_size(H5T_NATIVE_LONG))
    {
    H5Dread(dataset_id, H5T_NATIVE_LONG, memspace, dataspace, 
            //H5P_DEFAULT,
            plist_xfer,
            vtkfinal_id_list);
    }
  else {
    cerr << "type&size error while reading element connectivity\n";
  }
  H5Pclose (plist_xfer);

  H5Dclose(dataset_id);
  H5Sclose(memspace);

  minId = vtkIdTypeArray::GetDataTypeValueMax();
  maxId = -1;
  destptr = vtklistcells->GetPointer(0);

  if(numPieces == 1) // serial job, do not renumber
    {
    for(int i= 0 ; i< MyNumber_of_Cells; i++)
      {
      *destptr++ = 8;
      for(int j = 0; j < 8; j++)
        {
        if(!ParaSol)
          (*destptr)--;
        destptr++;
        }
      }
      MyNumber_of_Nodes = this->NbNodes; minId = 0;
    }
  else
    {
    vtkIdTypeArray* originalPtIds = vtkIdTypeArray::New();
    originalPtIds->SetNumberOfComponents(1);
    originalPtIds->SetName("GlobalNodeIds");

    for(int i= 0 ; i< MyNumber_of_Cells; i++)
      {
      *destptr++ = 8;
      for(int j = 0; j < 8; j++)
        {
        if(!ParaSol)
          (*destptr)--;
        if((*destptr) > maxId) maxId = (*destptr);
        if((*destptr) < minId) minId = (*destptr);
        destptr++;
        }
      }
    destptr = vtklistcells->GetPointer(0);
    for(int i= 0 ; i< MyNumber_of_Cells; i++)
      {
      destptr++; // skip the first which contains the value 8
      for(int j = 0; j < 8; j++)
        {
        (*destptr) -= minId; // all nodes now start from 0
        destptr++;
        }
      }

    MyNumber_of_Nodes = maxId - minId + 1;
    originalPtIds->SetNumberOfTuples(MyNumber_of_Nodes);
#ifdef PARALLEL_DEBUG
  errs << "\nminId  = " << minId << " maxId = " <<  maxId << "\n";
  errs << "\naLLocating " << MyNumber_of_Nodes << " OriginalPointIds of size " <<  sizeof(vtkIdType) << " bytes = " << MyNumber_of_Nodes *sizeof(vtkIdType) << " bytes\n";
#endif
    for(int i= 0 ; i< MyNumber_of_Nodes; i++)
      {
      originalPtIds->SetTuple1(i, minId+i);
      }
    output->GetPointData()->SetGlobalIds(originalPtIds);
    originalPtIds->FastDelete();
    }

/*/ control print-out
  for(int i= 0 ; i < 10; i++)
    {
    cerr << i << ": ";
    for(int j = 0; j < 9; j++)
      {
      cerr << *(vtkfinal_id_list+9*i+j) << ", " ;
      }
    cerr << "\n";
    }
/*/
  vtkCellArray *cells = vtkCellArray::New();

  cells->SetCells(MyNumber_of_Cells, vtklistcells);
  vtklistcells->FastDelete();

  int *types = new int[MyNumber_of_Cells];
#ifdef PARALLEL_DEBUG
  errs << "\naLLocating " << MyNumber_of_Cells << " celltypes of size " <<  sizeof(int) << " bytes = " << MyNumber_of_Cells *sizeof(int) << " bytes\n";
#endif
  for(int i=0; i < MyNumber_of_Cells; i++)
    {
    types[i] = VTK_HEXAHEDRON;
    }

  output->SetCells(types, cells);
  cells->FastDelete();
  delete [] types;
  this->UpdateProgress(0.50);
  vtkFloatArray *coords = vtkFloatArray::New();
  coords->SetNumberOfComponents(3);
  coords->SetNumberOfTuples(MyNumber_of_Nodes);
#ifdef PARALLEL_DEBUG
  errs << "\naLLocating " << MyNumber_of_Nodes  << " coordinates of size " <<  3*sizeof(float) << " bytes = " << MyNumber_of_Nodes * 3 *sizeof(float) << " bytes\n";
#endif

  count[0] = MyNumber_of_Nodes;
  count[1] = 3;
  memspace = H5Screate_simple(2, count, NULL);

  offset[0] = 0;
  offset[1] = 0;
  count[0] = MyNumber_of_Nodes;
  count[1] = 3;
  H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset, NULL, count, NULL);

  count[0] = MyNumber_of_Nodes;
  count[1] = 3;
  offset[0] = minId; // read the array starting at minId and read only MyNumber_of_Nodes items
  offset[1] = 0;

  dataspace = H5Dget_space(coords_id);
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

  H5Dread(coords_id, H5T_NATIVE_FLOAT, memspace, dataspace, H5P_DEFAULT, static_cast<vtkFloatArray *>(coords)->GetPointer(0));

  H5Dclose(coords_id);
  H5Sclose(memspace);

  vtkPoints *points = vtkPoints::New();
  points->SetData(coords);
  coords->FastDelete();
  output->SetPoints(points);
  points->FastDelete();

typedef struct material
  {
  int ids;
  double properties[2];
  } material;

  std::vector<int> uniqueMaterial_v;
  std::vector<float> merged_material;
  if(H5Lexists(mesh_id, "Material IDs", H5P_DEFAULT))
    {
    ids_id = H5Dopen(mesh_id, "Material IDs", H5P_DEFAULT);

     // here we're not sorting by materials, but provide the global mesh with an int array called "materials"
    Materials_Found = false;
    hsize_t     count[2];          /* size of the hyperslab in the file */
    hsize_t     offset[2];         /* hyperslab offset in the file */ 

    vtkIntArray *materials = vtkIntArray::New();
    materials->SetNumberOfComponents(1);
    materials->SetNumberOfTuples(MyNumber_of_Cells);
    materials->SetName("Material IDs");
    cerr << setprecision(10) << "allocating " << MyNumber_of_Cells * sizeof(int) << " bytes for materials" << endl;
    dataspace = H5Dget_space (ids_id);

    offset[0] = piece * load;
    offset[1] = 0;
    count[0]  = MyNumber_of_Cells;
    count[1]  = 1;

    status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

    count[0] = MyNumber_of_Cells;
    count[1] = 1;
    memspace = H5Screate_simple (2, count, NULL);   
    offset[0] = 0;
    offset[1] = 0;
    count[0] = MyNumber_of_Cells;
    count[1] = 1;
    status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset, NULL,  count, NULL);

    status = H5Dread(ids_id, H5T_NATIVE_INT, memspace, dataspace, H5P_DEFAULT, materials->GetPointer(0));
    output->GetCellData()->AddArray(materials);
    materials->Delete();
    H5Sclose(dataspace);
    H5Sclose(memspace);
    if(this->SortByMaterials)
      {
      //std::vector<int> Material_v;
      //Material_v.resize(this->NbCells); // read the full mesh, in order to be able to sort to find all material ids

      parameters_id = H5Gopen(root_id, "/Parameters", H5P_DEFAULT);
      dataset_id2 = H5Dopen(parameters_id, "#material types", H5P_DEFAULT);
      status = H5Dread(dataset_id2, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->NbMaterials);
      H5Dclose(dataset_id2);

      //int *mats = &Material_v[0];

      //H5Dread(ids_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, mats);

      //vtkIntArray *cxa = vtkIntArray::New();
      //cxa->SetName("Material IDs");
      //cxa->SetNumberOfTuples(MyNumber_of_Cells);

      //for(int i=0; i < MyNumber_of_Cells; i++)
        //cxa->SetTuple1(i, mats[i+piece * load]);

      //output->GetCellData()->AddArray(cxa);
      //cxa->FastDelete();
// now we can sort and it does not matter if we mess up the order since we have just copied the materials to array cxa

      hid_t    mats_id, compound_id;
      hsize_t  adims[1] = {this->NbMaterials};
      hid_t loctype = H5Tarray_create (H5T_NATIVE_DOUBLE, 1, adims);
      compound_id = H5Tcreate(H5T_COMPOUND, sizeof(material));
      H5Tinsert(compound_id, "ids", HOFFSET(material, ids), H5T_NATIVE_INT);
      H5Tinsert(compound_id, "properties", HOFFSET(material, properties), loctype);

      material *materials = new material[this->NbMaterials];
      mats_id = H5Dopen(parameters_id, "materials", H5P_DEFAULT);
      if(mats_id >=0 )
        {
        status = H5Dread(mats_id, compound_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &materials[0]);
        }
      else
        {
        cerr << "error reading materials ids\n";
        }

      H5Tclose (loctype);
      H5Dclose(mats_id);
      H5Gclose(parameters_id);
      for(int i=0; i < this->NbMaterials; i++)
        uniqueMaterial_v.push_back(materials[i].ids);
      delete [] materials;

      Materials_Found = true;
      }
    else

      H5Dclose(ids_id);
    }
  else if(H5Lexists(mesh_id, "Emoduli", H5P_DEFAULT))
    {
    ids_id = H5Dopen(mesh_id, "Emoduli", H5P_DEFAULT);
    //cerr << setprecision(10) << "opening Emoduli" <<  endl;
    Materials_Found = false;
    hsize_t     count[2];          /* size of the hyperslab in the file */
    hsize_t     offset[2];         /* hyperslab offset in the file */ 

    vtkFloatArray *vtk_materials = vtkFloatArray::New();
    vtk_materials->SetNumberOfTuples(MyNumber_of_Cells);
    vtk_materials->SetName("Material IDs");

    std::vector<float> material_fv;
    std::vector<float>::iterator it;
    material_fv.resize(MyNumber_of_Cells);
    float *materials = material_fv.data();
    //float *materials = (float *)vtk_materials->GetVoidPointer(0);
#ifdef PARALLEL_DEBUG
  errs << "\naLLocating " << MyNumber_of_Cells  << " materials of size " <<  sizeof(float) << " bytes = " << MyNumber_of_Cells *sizeof(float) << " bytes\n";
#endif

    dataspace = H5Dget_space (ids_id);

    offset[0] = piece * load;
    offset[1] = 0;
    count[0]  = MyNumber_of_Cells;
    count[1]  = 1;

    status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

    count[0] = MyNumber_of_Cells;
    count[1] = 1;
    memspace = H5Screate_simple (2, count, NULL);   
    offset[0] = 0;
    offset[1] = 0;
    count[0] = MyNumber_of_Cells;
    count[1] = 1;
    status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset, NULL,  count, NULL);

    status = H5Dread(ids_id, H5T_NATIVE_FLOAT, memspace, dataspace, H5P_DEFAULT, materials);
    for(int i=0; i < MyNumber_of_Cells; i++)
      vtk_materials->SetTuple1(i, materials[i]);
    output->GetCellData()->AddArray(vtk_materials);
    vtk_materials->Delete();
    H5Sclose(dataspace);
    H5Sclose(memspace);
    //material_fv.assign(materials, materials+MyNumber_of_Cells);

    if(this->SortByMaterials)
      if(numPieces==1)
      { // should find the unique values by sorting the full array

      std::vector<float>::iterator it;
      merged_material.resize(this->NbCells);
      float *mats = &merged_material[0];
      H5Dread(ids_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, mats);
      std::sort(merged_material.begin(), merged_material.end());
      it = std::unique(merged_material.begin(), merged_material.end());
      merged_material.resize(std::distance(merged_material.begin(), it));

      Materials_Found = true;
      }
#ifdef PARAVIEW_USE_MPI
    else
      { // sort the local array, unique() it, send to master process for gathering, sorting and re-unique() fnal result, then scatter final result
      std::sort(material_fv.begin(), material_fv.end());
      it = std::unique(material_fv.begin(), material_fv.end());
      material_fv.resize(std::distance(material_fv.begin(), it));
      int num = material_fv.size();
      float *ptr = material_fv.data();
#ifdef PARALLEL_DEBUG
      errs << "local list of sorted unique elements\n";
      for(int i=0; i < num; i++)
        {
        errs <<  ptr[i] << " ";
        }
      errs <<  "\n===========================\n";
#endif
/*

see http://www.mpi-forum.org/docs/mpi-1.1/mpi-11-html/node70.html
process "i" must send "num" floats to process root. the various values of num are not known to root,
so a separate gather must first be run to find these out. The data is placed contiguously at the receiving end.
*/
      int *rcounts = NULL;
      rcounts = (int *)malloc(numPieces*sizeof(int));
      MPI_Gather(&num, 1, MPI_INT, rcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
/*
  process root now has correct rcounts, using these we set displs[] so 
  that data is placed contiguously (or concatenated) at receive end
*/

      int i, *displs = NULL;
      int total_s = rcounts[numPieces-1];
      if(piece == 0)
        {
        displs = (int *)malloc(numPieces*sizeof(int)); 
        displs[0] = 0; 
        for (i=1; i<numPieces; ++i)
          { 
          displs[i] = displs[i-1]+rcounts[i-1]; 
          total_s += rcounts[i-1]; 
          } 

#ifdef PARALLEL_DEBUG
        errs << "root process will gather " << total_s << "  elements\n";
        for(i=0; i < numPieces; i++)
          {
          errs <<  displs[i] << " " << rcounts[i] <<"\n";
          }
#endif
        }
/*
http://www.mpi-forum.org/docs/mpi-1.1/mpi-11-html/node69.html#Node69
And, create receive buffers
*/ 
      float *rbuf = NULL;

      if(piece == 0)
        {
        merged_material.resize(total_s);
        rbuf = merged_material.data();  
        }
      MPI_Gatherv((void *)material_fv.data(), num, MPI_FLOAT, (void*)rbuf, (int *)rcounts, (int *)displs, MPI_FLOAT, 0, MPI_COMM_WORLD);
      if(piece == 0)
        {
        std::vector<float>::iterator it;
        std::sort(merged_material.begin(), merged_material.end());
        it = std::unique(merged_material.begin(), merged_material.end());
        merged_material.resize(std::distance(merged_material.begin(), it));
        rbuf = merged_material.data();
        total_s = merged_material.size();
        }
      MPI_Bcast((void*)&total_s, 1, MPI_INT, 0, MPI_COMM_WORLD);

      merged_material.resize(total_s);

      MPI_Bcast((void*)merged_material.data(), total_s, MPI_FLOAT, 0, MPI_COMM_WORLD);

#ifdef PARALLEL_DEBUG
      rbuf = merged_material.data();
      for(i=0; i < total_s; i++)
        {
        errs <<  rbuf[i] <<" ";
        }
#endif
      Materials_Found = true;
      }
#endif
      H5Dclose(ids_id);
    } // Emoduli case

  H5E_auto_t func;
  void *client_data;

  if(H5Gget_objinfo(root_id, "Solution", 0, NULL) >= 0)
    {
    solutions_id = H5Gopen(root_id, "/Solution", H5P_DEFAULT);
    }
  if(m_NumLoadSteps && actualTimeStep >= 0) {
    char tmp_string[256];
    sprintf(tmp_string, "/Solution/Load step %d", actualTimeStep);
    solutions_id = H5Gopen(root_id, tmp_string, H5P_DEFAULT);
  }
  if(this->PointDataArraySelection->ArrayIsEnabled("Nodal displacements"))
    {
    vtkFloatArray *displacements;

    if(H5Lexists(solutions_id, "Nodal displacements", H5P_DEFAULT) )
      displacements = ReadVectorAtNodes(solutions_id, "Nodal displacements", "Ux Uy Uz", minId, MyNumber_of_Nodes,
#ifdef COLLECTIVE
                                              plist_xfer);
#else
                                              H5P_DEFAULT);
#endif
    else if(H5Lexists(solutions_id, "disp", H5P_DEFAULT))
      displacements = ReadVectorAtNodes(solutions_id, "disp", "Ux Uy Uz", minId, MyNumber_of_Nodes,
#ifdef COLLECTIVE
                                              plist_xfer);
#else
                                              H5P_DEFAULT);
#endif
    else
      cerr << "No displacement vector found\n";
    output->GetPointData()->AddArray(displacements);
    displacements->Delete();
    output->GetPointData()->SetActiveAttribute("Nodal displacement", vtkDataSetAttributes::VECTORS);
    }

  if(this->PointDataArraySelection->ArrayIsEnabled("Nodal forces"))
    {
    vtkFloatArray *forces;

    if(H5Lexists(solutions_id, "Nodal forces", 0))
      forces = ReadVectorAtNodes(solutions_id, "Nodal forces", "Fx Fy Fz", minId, MyNumber_of_Nodes,
#ifdef COLLECTIVE
                                              plist_xfer);
#else
                                              H5P_DEFAULT);
#endif
    else if(H5Lexists(solutions_id, "force", 0))
      forces = ReadVectorAtNodes(solutions_id, "force", "Fx Fy Fz", minId, MyNumber_of_Nodes,
#ifdef COLLECTIVE
                                              plist_xfer);
#else
                                              H5P_DEFAULT);
#endif
    else
      cerr << "No force vector found\n";
    output->GetPointData()->AddArray(forces);
    forces->Delete();
    }

  this->UpdateProgress(0.70);
  for(int cell_var_index =0; cell_var_index < this->GetNumberOfCellArrays(); cell_var_index++)
    {
    if(this->CellDataArraySelection->GetArraySetting(cell_var_index))
      {
      if(strcmp(this->CellDataArraySelection->GetArrayName(cell_var_index), "Emoduli"))
        {
        vtkFloatArray *data = ReadScalarAtCells(solutions_id,
                                              cell_var_index,
                                              piece * load,
                                              m_NumStateVariables,
                                              this->CellDataArraySelection->GetArrayName(cell_var_index),
                                              MyNumber_of_Cells,
#ifdef COLLECTIVE
                                              plist_xfer);
#else
                                              H5P_DEFAULT);
#endif
        output->GetCellData()->AddArray(data);
        data->Delete();
        }
      else
        {
        vtkFloatArray *data = ReadScalarAtCells(mesh_id,
                                              cell_var_index,
                                              piece * load,
                                              m_NumStateVariables,
                                              this->CellDataArraySelection->GetArrayName(cell_var_index),
                                              MyNumber_of_Cells,
#ifdef COLLECTIVE
                                              plist_xfer);
#else
                                              H5P_DEFAULT);
#endif
        output->GetCellData()->AddArray(data);
        data->Delete();
        }
      }  // if cell was activated
    } // for all cell arrays

  H5Gclose(mesh_id);
  if(Materials_Found)
    {
    vtkSelectionNode *selectionNode = vtkSelectionNode::New();
    selectionNode->SetFieldType(vtkSelectionNode::CELL);

// see http://www.vtk.org/pipermail/vtkusers/2013-February/127538.html
// to know why selecting by VALUES is not fast
    selectionNode->SetContentType(vtkSelectionNode::THRESHOLDS);
    vtkDoubleArray *thismat = vtkDoubleArray::New();
    thismat->SetNumberOfTuples(2);

    thismat->SetName("Material IDs");

    std::vector<float>::iterator ck = merged_material.begin();
    char name[64];

    vtkSelection *selection = vtkSelection::New();
    selection->AddNode(selectionNode);

    gettimeofday(&tv0, NULL);
    for (int i=0; ck != merged_material.end(); ck++, i++)
      {
      //cerr << "Extracting material id " << *ck << "...";
      thismat->SetTuple1(0, *ck); //selection based on threshold requirs a min and a max. Here they're equal
      thismat->SetTuple1(1, *ck);

      sprintf(name, "%d", (int)*ck);
      //selectionNode->Initialize();

      selectionNode->SetSelectionList(thismat);

      vtkExtractSelection *extThres = vtkExtractSelection::New();
      extThres->SetInputData(1, selection);
      extThres->SetInputData(0, output);
      extThres->PreserveTopologyOff();
      extThres->DebugOn();
      extThres->Update();

      mb->SetBlock(i, vtkUnstructuredGrid::SafeDownCast(extThres->GetOutput()));
      mb->GetMetaData(i)->Set(vtkCompositeDataSet::NAME(), name);
      extThres->FastDelete();
      //cerr << "done\n";
      I = i+1;
      }
  gettimeofday(&tv1, NULL);
  timersub(&tv1, &tv0, &res);
  //cerr << "ExtractionByMaterials time elapsed = " << res.tv_sec  << "." << res.tv_usec << endl;
    thismat->Delete();
    selectionNode->Delete();
    selection->Delete();
    }
  else
    {
    I=0;
    mb->SetBlock(I, output);
    mb->GetMetaData(I)->Set(vtkCompositeDataSet::NAME(), "mesh");
    }
  output->Delete();

  if(solutions_id)
    H5Gclose(solutions_id);
  H5Gclose(root_id);
  H5Fclose(file_id);
  this->UpdateProgress(1.0);
  vtkDebugMacro( << "RequestData(END)");
#ifdef PARALLEL_DEBUG
  errs.close();
#endif

  return 1;
}
/*
h5ls -r /local/data/BioMedical/jus_f10_1032r_intact_bone_120um_r1vox.h5 
/                        Group
/Image_Data              Group
/Image_Data/Fixed_Displacement_Coordinates Dataset {38778, 4}
/Image_Data/Fixed_Displacement_Values Dataset {38778}
/Image_Data/Image        Dataset {851, 357, 337}
/Image_Data/Poisons_ratio Dataset {SCALAR}
/Image_Data/Voxelsize    Dataset {SCALAR}
/Mesh                    Group
/Mesh/Coordinates        Dataset {14201004, 3}
/Mesh/Elements           Dataset {10200221, 8}
/Mesh/Emoduli            Dataset {10200221, 1}
/Solution                Group
/Solution/EFF            Dataset {10200221, 1}
/Solution/disp           Dataset {14201004, 3}
/Solution/e_vol          Dataset {10200221, 1}
/Solution/force          Dataset {14201004, 3}

piece 0 out of 2
CPU 0 aLLocating 5100110 list of elements of size *9*8 bytes = 367207920 bytes

aLLocating 10559297 OriginalPointIds of size 8 bytes = 84474376 bytes

aLLocating 5100110 celltypes of size 4 bytes = 20400440 bytes

aLLocating 10559297 coordinates of size 12 bytes = 126711564 bytes
aLLocating 5100110 materials of size 4 bytes = 20400440 bytes
local list of sorted unique elements
18000 
===========================
root process will gather 2  elements
0 1
1 1
18000
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
piece 1 out of 2
CPU 1 aLLocating 5100111 list of elements of size *9*8 bytes = 367207992 bytes

aLLocating 5722929 OriginalPointIds of size 8 bytes = 45783432 bytes

aLLocating 5100111 celltypes of size 4 bytes = 20400444 bytes

aLLocating 5722929 coordinates of size 12 bytes = 68675148 bytes
aLLocating 5100111 materials of size 4 bytes = 20400444 bytes
local list of sorted unique elements
18000 
===========================
18000


*/
