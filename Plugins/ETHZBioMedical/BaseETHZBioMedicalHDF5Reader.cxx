#include "BaseETHZBioMedicalHDF5Reader.h"
#ifdef ETHZBioMedical_EXPORTS
#include "vtkFloatArray.h"
#endif

#ifdef ETHZBioMedical_EXPORTS
vtkFloatArray *
#else
float *
#endif
ReadScalarAtCells(hid_t solutions_id, int cell_var_index, int offset0, int m_NumStateVariables, const char *label1, int size, hid_t plist_xfer)
{
  hid_t dataset_id;
  hid_t memspace, dataspace;
  hsize_t count[2], offset[2];
  herr_t   status;
  bool     ParaSol = false;
  vtkFloatArray *data = vtkFloatArray::New();
  data->SetNumberOfComponents(1);
  data->SetNumberOfTuples(size);
  data->SetName(label1);
  //cerr << setprecision(10) << "allocating " << size * sizeof(float) << " bytes for " <<  label1 << endl;

//check if datasets exists; this first test could only be the paraSol case if it is true
  if(H5Lexists(solutions_id, label1, H5P_DEFAULT))
    {
    dataset_id = H5Dopen(solutions_id, label1, H5P_DEFAULT);
    ParaSol = true;
    }
// from here down, the older paraFE format
  else if(cell_var_index < 8)
    {
    if(H5Gget_objinfo(solutions_id, "Element strain", 0, NULL) >= 0)
      dataset_id = H5Dopen(solutions_id, "Element strain", H5P_DEFAULT);
    else return NULL;
    }
  else if (cell_var_index < 15)
    {
    if(H5Gget_objinfo(solutions_id, "Element stress", 0, NULL) >= 0)
      dataset_id = H5Dopen(solutions_id, "Element stress", H5P_DEFAULT);
    else return NULL;
    }
  else if (cell_var_index < 15 + m_NumStateVariables)
    {
    if(H5Gget_objinfo(solutions_id, "State Variables", 0, NULL) >= 0)
      dataset_id = H5Dopen(solutions_id, "State Variables", H5P_DEFAULT);
    else return NULL;
    }

  dataspace = H5Dget_space (dataset_id);

  offset[0] = offset0;
  if(!ParaSol)
    {
    if(cell_var_index < 8)
      offset[1] = cell_var_index;
    else if (cell_var_index < 15)
      offset[1] = cell_var_index - 8;
    else if (cell_var_index <= 15+m_NumStateVariables)
      offset[1] = cell_var_index - 15;
    }
  else
    offset[1] = 0;

  count[0]  = size;
  count[1]  = 1;

  status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

  count[0] = size;
  count[1] = 1;
  memspace = H5Screate_simple (2, count, NULL);   
  offset[0] = 0;
  offset[1] = 0;
  count[0] = size;
  count[1] = 1;
  status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset, NULL,  count, NULL);

  status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, memspace, dataspace, plist_xfer, data->GetPointer(0));

  H5Dclose(dataset_id);
  return data;
}

#ifdef ETHZBioMedical_EXPORTS
vtkFloatArray *
#else
float *
#endif
ReadVectorAtNodes(hid_t solutions_id, const char *label1, const char *label2, int minId, int size, hid_t plist_xfer)
{
  hid_t s1_tid, s2_tid, s3_tid, dataset_id;
  hid_t memspace, dataspace;
  hsize_t count[2], offset[2];

  vtkFloatArray *data = vtkFloatArray::New();
  data->SetNumberOfComponents(3);
  data->SetNumberOfTuples(size);
  data->SetName(label1);
  //cerr << setprecision(10) << "allocating " << size * 3 * sizeof(float) << " bytes for " << label1 << endl;

  dataset_id = H5Dopen(solutions_id, label1, H5P_DEFAULT);
  s1_tid     = H5Dget_type(dataset_id);
  
  hsize_t     dims_out[2];
  dataspace = H5Dget_space(dataset_id);
  int rank  = H5Sget_simple_extent_ndims(dataspace);
  H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
  //printf("rank %d, dimensions %lu x %lu \n", rank, (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]));

  if(H5Tget_class(s1_tid) == H5T_COMPOUND)
    {
    //* superfluous?
    count[0] = size;
    count[1] = 1;
    memspace = H5Screate_simple(2, count, NULL);

    offset[0] = 0;
    offset[1] = 0;
    count[0] = size;
    count[1] = 1;
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    //*/
    count[0] = size;
    count[1] = 1;
    offset[0] = minId;
    offset[1] = 0;

    dataspace = H5Dget_space(dataset_id);
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
  
    hsize_t    array_dim[] = {3};
    int array_rank = 1;
    s2_tid = H5Tarray_create(H5T_NATIVE_FLOAT, array_rank, array_dim);
    /* Create the memory data type                                    */
    s3_tid = H5Tcreate (H5T_COMPOUND, 3 * sizeof(float));
    H5Tinsert(s3_tid, label2, 0, s2_tid);

    H5Dread(dataset_id, s3_tid, memspace, dataspace, plist_xfer, data->GetPointer(0));
    H5Tclose(s2_tid);
    H5Tclose(s3_tid);
    }
  else
    {
    //* not superfluous since the parallel reads breaks right here
    count[0] = size;
    count[1] = 3;
    memspace = H5Screate_simple(2, count, NULL);

    offset[0] = 0;
    offset[1] = 0;
    count[0] = size;
    count[1] = 3;
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    //*/ not superfluous
    count[0] = size;
    count[1] = 3;
    offset[0] = minId;
    offset[1] = 0;

    dataspace = H5Dget_space(dataset_id);
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

    H5Dread(dataset_id, H5T_NATIVE_FLOAT, memspace, dataspace, plist_xfer, data->GetPointer(0));
    }
  H5Dclose(dataset_id);
  return data;
}
