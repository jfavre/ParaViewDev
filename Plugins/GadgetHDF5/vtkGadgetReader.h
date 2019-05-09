// .NAME vtkGadgetReader
// .SECTION Description
// vtkGadgetReader reads collection of "Gadget" HDF5 files 
// 
// .SECTION Thanks
// Jean M. Favre
// CSCS - Swiss National Supercomputing Centre for creating and contributing
// this class.

#ifndef vtkGadgetReader_h
#define vtkGadgetReader_h

#include <string>
#include <vector>

#define ALL_TYPES 1
#define OUTPUT_UG 1

#ifdef ALL_TYPES
#include "vtkMultiBlockDataSetAlgorithm.h"
static std::vector<std::string> ParticleTypes = {"PartType0", "PartType1", "PartType2", "PartType3", "PartType4", "PartType5"};
#else
#ifdef OUTPUT_UG
#include "vtkUnstructuredGridAlgorithm.h"
#else
#include "vtkPolyDataAlgorithm.h"
#endif
static std::vector<std::string> ParticleTypes = {"PartType0"};
#endif

#include <map>
#include <sstream>

class vtkDataArraySelection;
class vtkStdString;
class vtkMultiProcessController;

enum  ParticleType : int {Gas=0, Halo=1, Disk=2, Bulge=3, Stars=4, Bndry=5};
enum  CellTypes : int {None=0, Vertex=1, PolyVertex=2};

#ifdef ALL_TYPES
class vtkGadgetReader : public vtkMultiBlockDataSetAlgorithm
#else
#ifdef OUTPUT_UG
class vtkGadgetReader : public vtkUnstructuredGridAlgorithm
#else
class vtkGadgetReader : public vtkPolyDataAlgorithm
#endif
#endif
{
public:
  static vtkGadgetReader *New();
#ifdef ALL_TYPES
  vtkTypeMacro(vtkGadgetReader,vtkMultiBlockDataSetAlgorithm);
#else
#ifdef OUTPUT_UG
  vtkTypeMacro(vtkGadgetReader,vtkUnstructuredGridAlgorithm);
#else
  vtkTypeMacro(vtkGadgetReader,vtkPolyDataAlgorithm);
#endif
#endif
  void PrintSelf(ostream& os, vtkIndent indent);

  void SetDirectoryName(const char* dn);
  vtkGetStringMacro(DirectoryName);

  // Description:
  // When using ParaView, cell generation is recommended, without them
  // many filter operations are unavailable
  // Note that Point Gaussian Rendering does not require cells. 
  // When CellType == Vertex, the reader will generate one vertex cell
  // for each point/particle read.
  // When CellType == PolyVertex, the reader will generate a single cell
  // for each ParticleType
  vtkSetMacro(CellType, int);
  vtkGetMacro(CellType, int);

  int         GetNumberOfPointArrays();
  const char* GetPointArrayName(int index);
  int         GetPointArrayStatus(const char* name);
  void        SetPointArrayStatus(const char* name, int status);
  void        DisableAllPointArrays();
  void        EnableAllPointArrays();
  //
  int         GetNumberOfPointArrayStatusArrays() { return GetNumberOfPointArrays(); }
  const char* GetPointArrayStatusArrayName(int index) { return GetPointArrayName(index); }
  int         GetPointArrayStatusArrayStatus(const char* name) { return GetPointArrayStatus(name); }
  void        SetPointArrayStatusArrayStatus(const char* name, int status) { SetPointArrayStatus(name, status); }

#ifdef PARAVIEW_USE_MPI

    // Description:
    // Set/Get the controller use in compositing (set to
    // the global controller by default)
    // If not using the default, this must be called before any
    // other methods.
    virtual void SetController(vtkMultiProcessController* controller);
    vtkGetObjectMacro(Controller, vtkMultiProcessController);

#endif

protected:
   vtkGadgetReader();
  ~vtkGadgetReader();
  //
  int   RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int   RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int   OpenFile();
  void  CloseFile();

  //
  // Internal Variables
  //
  char*         DirectoryName;
  int           CellType;
  vtkTimeStamp  FileModifiedTime;
  vtkTimeStamp  FileOpenedTime;
  int           UpdatePiece;
  int           UpdateNumPieces;
  bool          PartTypes[6];
  int           NumPart_Total[6];
  typedef std::vector<std::string>  stringlist;
  std::map<std::string, stringlist> FieldArrays;
  stringlist                        GadgetFileNames;
  // To allow paraview gui to enable/disable scalar reading
  vtkDataArraySelection* PointDataArraySelection;

  #ifdef PARAVIEW_USE_MPI
  vtkMultiProcessController* Controller;
  #endif

private:
  vtkGadgetReader(const vtkGadgetReader&) VTK_DELETE_FUNCTION;
  void operator=(const vtkGadgetReader&) VTK_DELETE_FUNCTION;
};

#endif
