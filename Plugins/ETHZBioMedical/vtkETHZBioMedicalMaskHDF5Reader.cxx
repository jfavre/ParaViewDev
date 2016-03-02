#include "vtkETHZBioMedicalMaskHDF5Reader.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDataArraySelection.h"
#include "vtkDataSetAttributes.h"
#include "vtkDoubleArray.h"
#include "vtkErrorCode.h"
#include "vtkFieldData.h"
#include "vtkFloatArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkCharArray.h"
#include "vtkShortArray.h"
#include "vtkTimerLog.h"
#include "vtkSelection.h"
#include "vtkSelectionNode.h"
#include "vtkExtractSelection.h"
#include "vtkIdTypeArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnstructuredGrid.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkUniformGrid.h"

#include <algorithm>
#include <vector>
using namespace std;
#include <vtksys/SystemTools.hxx>


//vtkCxxRevisionMacro(vtkETHZBioMedicalMaskHDF5Reader, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkETHZBioMedicalMaskHDF5Reader);

int vtkETHZBioMedicalMaskHDF5Reader::FillOutputPortInformation(int port, vtkInformation* info)
{
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  if(port == 1)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUniformGrid");
  return 1;
}


vtkETHZBioMedicalMaskHDF5Reader::vtkETHZBioMedicalMaskHDF5Reader()
{
  this->AIMReader = vtkAIMReader::New();
  this->HDF5Reader= vtkETHZBioMedicalHDF5Reader::New();
  this->Masks = vtkDataArraySelection::New();
  this->SetNumberOfInputPorts(0);
#ifdef MULTI_PORT
  this->SetNumberOfOutputPorts(2);
#else
  this->SetNumberOfOutputPorts(1);
#endif
  this->HDF5FileName = NULL;
  this->FileName =  NULL;
  
  //vtkImageData *ug = vtkImageData::New();
  //ug->ReleaseData();
  //this->GetExecutive()->SetOutputData(1, ug);
  //ug->Delete();

}

vtkETHZBioMedicalMaskHDF5Reader::~vtkETHZBioMedicalMaskHDF5Reader()
{
  this->AIMReader->Delete();
  this->HDF5Reader->Delete();
  this->Masks->Delete();
}

int vtkETHZBioMedicalMaskHDF5Reader::RequestInformation(vtkInformation* request,
                         vtkInformationVector** inputVector,
                         vtkInformationVector* outputVector)
{
  vtkDebugMacro( << "RequestInformation(BEGIN)");
  
  this->AIMReader->SetFileName(this->FileName);
  this->AIMReader->DebugOff();
  this->AIMReader->Update();
  vtkCellData *cd = this->AIMReader->GetOutput()->GetCellData();
  vtkDataArray *da = cd->GetArray("aim_data");
  if(da->IsA("vtkShortArray"))
    {
    cerr << "short array: not implemented yet\n";
    }
  else if(da->IsA("vtkCharArray"))
    {
    char *data = static_cast<vtkCharArray*>(da)->GetPointer(0);
    cerr << "copying array of size " << da->GetNumberOfTuples() <<endl;
    vtkTimerLog *tl = vtkTimerLog::New();

    std::vector<char> tags2;
    std::vector<char>::iterator tagit;
    int tagit2, CurrentSize;
    tl->StartTimer();
#define    USING_STL
#ifdef USING_STL
    std::vector<char> tags(data, data + da->GetNumberOfTuples());
    sort(tags.begin(), tags.end());
    tags2.assign( tags.begin(), std::unique(tags.begin(), tags.end()));
#else
    cerr << "0: found new tag: " << (int)(*data) << endl;
    tags2.push_back(*data);
    CurrentSize = 1;
    for (int t = 0; t < da->GetNumberOfTuples(); t++, data++)
      {
      char test = *data;
      bool NotFound = true;
      tagit2 = 0; 
      while(NotFound && (tagit2 < CurrentSize))
        {
	cerr << "t = " << t << ", " << (int)tags2[tagit2] << ", " << (int)test << " ! ";
        if (tags2[tagit2] == test);
	  {
	  //cerr << "equal tag: " << tagname << endl;
          NotFound = false;
	  }
	tagit2++;
        } 
	cerr  << "out of w loop\n";
      if((tagit2 == CurrentSize) && (150==t))
        {
	cerr << "NotFound = " <<NotFound << endl;
	cerr << "found new tag: " << (int)(*data) << endl;
        tags2.push_back(*data);
	CurrentSize++;
	}
      }
#endif
    tl->StopTimer();
    cerr << "time used to find tags: " << tl->GetElapsedTime() << endl;
    tl->Delete();
    
    char tagname[16];
    for (tagit = tags2.begin(); tagit != tags2.end(); tagit++)
      {
      sprintf(tagname,"%d", (int)(*tagit));
      cerr << "adding tag to GUI: " << tagname << endl;
      this->Masks->AddArray(tagname);
      }
    }
  else if(da->IsA("vtkFloatArray"))
    {
    cerr << "float array: not implemented yet\n";
    }
  this->HDF5Reader->PointDataArraySelection->DebugOff();
  this->HDF5Reader->PointDataArraySelection->AddArray("Nodal displacements");
  this->HDF5Reader->PointDataArraySelection->AddArray("Nodal forces");
  
  this->HDF5Reader->CellDataArraySelection->DebugOff();
  for(int i=0; i < 7; i++)
    {
    this->HDF5Reader->CellDataArraySelection->AddArray(vtkETHZBioMedicalReader::strain_var_names[i]);
    }
  this->HDF5Reader->CellDataArraySelection->AddArray("EFF");
  for(int i=0; i < 7; i++)
    {
    this->HDF5Reader->CellDataArraySelection->AddArray(vtkETHZBioMedicalReader::stress_var_names[i]);
    }
  int my_extents[6], dims[3];
#ifdef MULTI_PORT
  vtkInformation* outInfo = outputVector->GetInformationObject(1);
  this->AIMReader->GetOutput()->GetDimensions(dims);
  my_extents[0] = 0;
  my_extents[1] = dims[0]-1;
  my_extents[2] = 0;
  my_extents[3] = dims[1]-1;
  my_extents[4] = 0;
  my_extents[5] = dims[2]-1;
  cerr << "Set WHOLE_EXTENT  extents: " << my_extents[0] << ", " <<
                                                    my_extents[1] << ", " <<
                                                    my_extents[2] << ", " <<
                                                    my_extents[3] << ", " <<
                                                    my_extents[4] << ", " <<
                                                    my_extents[5] << endl;
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), my_extents, 6);
#endif
  vtkDebugMacro( << "\nRequestInformation(END)");
  return 1;
}

int vtkETHZBioMedicalMaskHDF5Reader::RequestData(vtkInformation* request,
                         vtkInformationVector** inputVector,
                         vtkInformationVector* outputVector)
{
  vtkDebugMacro( << "RequestData(BEGIN)");
  
  vtkInformation* info = outputVector->GetInformationObject(0);
  vtkDataObject* doOutput = info->Get(vtkDataObject::DATA_OBJECT());
  vtkMultiBlockDataSet* mb = vtkMultiBlockDataSet::SafeDownCast(doOutput);
  if (!mb)
    {
    return 0;
    }
#ifdef MULTI_PORT
  vtkInformation* info1 = outputVector->GetInformationObject(1);
  vtkDataObject* doOutput1 = info1->Get(vtkDataObject::DATA_OBJECT());
  vtkUniformGrid* ug = vtkUniformGrid::SafeDownCast(doOutput1);
  if (!ug)
    {
    return 0;
    }
  else
    {
    int dimensions[3];
    this->AIMReader->GetOutput()->GetDimensions(dimensions);
    ug->SetSpacing(.02,.02,.02);
    //ug->SetWholeExtent(0, dimensions[0]-1, 0, dimensions[1]-1, 0, dimensions[2]-1);
    ug->SetOrigin(0.0, 0.0, 0.0);
    ug->SetDimensions(dimensions[0], dimensions[1], dimensions[2]);
    }
#endif
  vtkDataArray *da = this->AIMReader->GetOutput()->GetCellData()->GetArray("aim_data");
  
  vtkUnsignedCharArray *vis = vtkUnsignedCharArray::New();
  vis->SetNumberOfTuples(da->GetNumberOfTuples());
#ifdef MULTI_PORT
  ug->SetCellVisibilityArray(vis);
  ug->GetCellData()->SetScalars(vis); 
  vis->Delete();
#endif
  char *data = static_cast<vtkCharArray*>(da)->GetPointer(0);
  for(int k=0; k < this->GetNumberOfMaskArrays(); k++)
    {
    unsigned char t = static_cast<vtkCharArray*>(da)->GetValue(k);
    vis->SetValue(k, t);
    }
    vis->Delete();


  data = static_cast<vtkCharArray*>(da)->GetPointer(0);
  if (! this->HDF5FileName || ! vtksys::SystemTools::FileExists(this->HDF5FileName))
     {
     //cerr << "no HDF5 Filename was set\n";
     std::string   path = vtksys::SystemTools::GetFilenamePath(this->FileName);
     std::string   prefix = vtksys::SystemTools::GetFilenameWithoutLastExtension(this->FileName);
     std::string meshFileName = path + "/" + prefix + std::string(".h5");
     if(! vtksys::SystemTools::FileExists(meshFileName.c_str()))
       {
       cerr << "No HDF5 filename was set, and no matching file was found. Can't read data\n";
       return 0;
       }
     else
       this->HDF5Reader->SetFileName(meshFileName.c_str());
     }
  else
    this->HDF5Reader->SetFileName(this->HDF5FileName);

  this->HDF5Reader->DebugOff();
  //this->HDF5Reader->DisableAllCellArrays();
  //this->HDF5Reader->DisableAllPointArrays();
  this->HDF5Reader->Update();
  unsigned int valid_masks = 0;
  for(int k=0; k < this->GetNumberOfMaskArrays(); k++)
    {
    if(this->GetMaskArrayStatus(this->GetMaskArrayName(k)))
      {
      //cout << "found tag turned ON: " << this->GetMaskArrayName(k) << endl;
      data = static_cast<vtkCharArray*>(da)->GetPointer(0);
      vtkIdTypeArray *masks = vtkIdTypeArray::New();
      masks->SetNumberOfComponents(1);
      masks->SetName("masks");
      int H5Mesh_Cell_Number=0; // cell_number in the HDF5 mesh file
      cerr << "search for tag " << atoi(this->GetMaskArrayName(k)) << endl;
      for(int i=0; i < da->GetNumberOfTuples(); i++, data++)
        {
        if(*data == static_cast<char>(atoi(this->GetMaskArrayName(k))))
          {
        //cerr << "inserted " << *data << endl;
          masks->InsertNextTuple1(H5Mesh_Cell_Number);
          }
        if(*data != static_cast<char>(0))
          H5Mesh_Cell_Number++; // cell number has increased since we have a non-null cell
        }
/*
    ofstream *fs = new ofstream("/tmp/mask.txt", ios::out);
    if (fs->fail())
      {
      cerr << "Unable to open mask file: "<< endl;
      delete fs;
      fs = NULL;
      return 0;
      }
    else
      {
      int i=0;
      //for(; i < masks->GetNumberOfTuples(); i++){(*fs) << masks->GetValue(i) << endl;}
      i=1;
      int start, previous;
      start = previous = masks->GetValue(0);
      while(i < masks->GetNumberOfTuples())
        {
	int current = masks->GetValue(i);
	while(current == (previous + 1))
	  {
	  previous = current;
	  i++;
	  current = masks->GetValue(i);
	  } // once we break out, we have just finished a sequence
	(*fs) << start << " " << previous << endl;
	start = current; 
	previous = current - 1; 
	}
      }
    delete fs;
      */
      vtkSelection *sel = vtkSelection::New();
      vtkSelectionNode *selnode = vtkSelectionNode::New();
      selnode->SetFieldType(vtkSelection::CELL);
      selnode->SetContentType(vtkSelectionNode::INDICES);
      selnode->SetSelectionList(masks);
      sel->AddNode(selnode);
      masks->Delete();
  
      vtkExtractSelection *extThres = vtkExtractSelection::New();
      extThres->SetInputData(1, sel);
      extThres->SetInputConnection(0, this->HDF5Reader->GetOutputPort(0));
      extThres->Update();
      //output->ShallowCopy(extThres->GetOutput());
      mb->SetBlock(valid_masks, extThres->GetOutput());
      mb->GetMetaData(valid_masks)->Set(vtkCompositeDataSet::NAME(), this->GetMaskArrayName(k));

      extThres->Delete();
      selnode->Delete();
      sel->Delete();
      valid_masks++;
      }
    }
  vtkDebugMacro( << "\nRequestData(END)");
  return 1;
}
