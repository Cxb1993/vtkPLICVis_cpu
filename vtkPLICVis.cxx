#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkSmartPointer.h"
#include "vtkRectilinearGrid.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"
#include "vtkMPIController.h"

#include "vtkPLICVis.h"
#include "plicvis_impl.h"

#include <vector>
#include <array>

vtkStandardNewMacro(vtkPLICVis);

//----------------------------------------------------------------------------
vtkPLICVis::vtkPLICVis():
  NumGhostLevels(1)
{
  this->SetNumberOfInputPorts(1);
}
//----------------------------------------------------------------------------
int vtkPLICVis::FillInputPortInformation(int port, vtkInformation* info)
{
  if (!this->Superclass::FillInputPortInformation(port, info)) {
    return 0;
  }
  if (port == 0) {
    info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkRectilinearGrid");
    return 1;
  }
  return 0;
}

//----------------------------------------------------------------------------
void vtkPLICVis::RemoveGhostCellsFromExtent(const int extent[6],
					    int fieldExtent[6])
{
  vtkMPIController *controller = vtkMPIController::New();
  if (controller->GetCommunicator() == 0) {
    controller->Delete();
    return;
  }  

  const int numSides = 6;
  int processId = controller->GetLocalProcessId();
  int numProcesses = controller->GetNumberOfProcesses();

  // prepare buffers for communication ---------------------------------------
  std::vector<vtkIdType> recvLengths(numProcesses);
  std::vector<vtkIdType> recvOffsets(numProcesses);
  for (int i = 0; i < numProcesses; ++i) {
    recvLengths[i] = numSides;
    recvOffsets[i] = i*numSides;
  }
  
  // find global extent ------------------------------------------------------
  std::vector<int> allExtents(numSides*numProcesses);
  controller->AllGatherV(&extent[0], &allExtents[0], numSides, 
			 &recvLengths[0], &recvOffsets[0]);

  int globalExtent[6];

  globalExtent[0] = globalExtent[2] = globalExtent[4] = std::numeric_limits<int>::max();
  globalExtent[1] = globalExtent[3] = globalExtent[5] = - globalExtent[0];

  for (int i = 0; i < allExtents.size()/6; ++i) {
    if (globalExtent[0] > allExtents[i*6+0]) globalExtent[0] = allExtents[i*6+0];
    if (globalExtent[1] < allExtents[i*6+1]) globalExtent[1] = allExtents[i*6+1];
    if (globalExtent[2] > allExtents[i*6+2]) globalExtent[2] = allExtents[i*6+2];
    if (globalExtent[3] < allExtents[i*6+3]) globalExtent[3] = allExtents[i*6+3];
    if (globalExtent[4] > allExtents[i*6+4]) globalExtent[4] = allExtents[i*6+4];
    if (globalExtent[5] < allExtents[i*6+5]) globalExtent[5] = allExtents[i*6+5];
  }

  // reduce field extent to one without ghost cells --------------------------------  
  if (extent[0] > globalExtent[0]) fieldExtent[0] += NumGhostLevels;
  if (extent[1] < globalExtent[1]) fieldExtent[1] -= NumGhostLevels;
  if (extent[2] > globalExtent[2]) fieldExtent[2] += NumGhostLevels;
  if (extent[3] < globalExtent[3]) fieldExtent[3] -= NumGhostLevels;
  if (extent[4] > globalExtent[4]) fieldExtent[4] += NumGhostLevels;
  if (extent[5] < globalExtent[5]) fieldExtent[5] -= NumGhostLevels;

  controller->Delete();
}

//----------------------------------------------------------------------------
int vtkPLICVis::RequestUpdateExtent(vtkInformation *vtkNotUsed(request),
				    vtkInformationVector **inputVector,
				    vtkInformationVector *outputVector)
{
  // set one ghost level -----------------------------------------------------
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

  int ngl = 0;
  ngl = inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());
  if (ngl == 0) {
    NumGhostLevels = 1;
    inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), NumGhostLevels);
  }
  else {
    NumGhostLevels = ngl;
  }
  
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), NumGhostLevels);

  return 1;
}
//----------------------------------------------------------------------------
int vtkPLICVis::RequestData(vtkInformation *vtkNotUsed(request),
			    vtkInformationVector **inputVector,
			    vtkInformationVector *outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkRectilinearGrid *input = vtkRectilinearGrid::
    SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkPolyData *output = vtkPolyData::
      SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkDataArray *coords[3] = {input->GetXCoordinates(),
			     input->GetYCoordinates(),
			     input->GetZCoordinates()};
  std::vector<float> dx;
  std::vector<float> dy;
  std::vector<float> dz;

  nodeCoordsToEdgeLengths(coords[0], dx);
  nodeCoordsToEdgeLengths(coords[1], dy);
  nodeCoordsToEdgeLengths(coords[2], dz);

  std::vector<float3> vertices(0);
  std::vector<int> indices(0);

  int nodeRes[3];
  input->GetDimensions(nodeRes);
  int cellRes[3] = {nodeRes[0]-1,nodeRes[1]-1,nodeRes[2]-1};

  
  int extent[6];
  input->GetExtent(extent);    
  int fieldExtent[6] = {0, cellRes[0], 0, cellRes[1], 0, cellRes[2]};
  RemoveGhostCellsFromExtent(extent, fieldExtent);
    
  int imin = fieldExtent[0];
  int imax = fieldExtent[1];
  int jmin = fieldExtent[2];
  int jmax = fieldExtent[3];
  int kmin = fieldExtent[4];
  int kmax = fieldExtent[5];

  if (cellRes[0] < 2) {
    imin = 0;
    imax = 1;
  }
  if (cellRes[1] < 2) {
    jmin = 0;
    jmax = 1;
  }
  if (cellRes[2] < 2) {
    kmin = 0;
    kmax = 1;
  }


  int memSize = 0;  
  vtkDataArray *data = input->GetCellData()->GetArray("Data");
  memSize = data->GetActualMemorySize();

  std::cout << "memSize = " << memSize << std::endl;

  const int memAval = 128*256*256*4*13;

  std::vector<std::array<int,6>> dataChunks;
  dataChunks.resize(1);
  dataChunks[0][0] = imin;
  dataChunks[0][1] = imax;
  dataChunks[0][2] = jmin;
  dataChunks[0][3] = jmax;
  dataChunks[0][4] = kmin;
  dataChunks[0][5] = kmax;
  
  int vertexID = 0;
  for (int c = 0; c < dataChunks.size(); ++c) {
  
    for (int k = dataChunks[c][4]; k < dataChunks[c][5]; ++k) {
      float oz = coords[2]->GetComponent(k,0);
      for (int j = dataChunks[c][2]; j < dataChunks[c][3]; ++j) {
	float oy = coords[1]->GetComponent(j,0);
	for (int i = dataChunks[c][0]; i < dataChunks[c][1]; ++i) {
	  float ox = coords[0]->GetComponent(i,0);

	  int idx = i + j*cellRes[0] + k*cellRes[0]*cellRes[1];
	  float f = data->GetComponent(idx, 0);
	  if (f <= EMF0) {
	    continue;
	  }
	  if (f >= EMF1 && !interfaceCell(data, i, j, k, cellRes)) {
	    continue;
	  }

	  float3 grad = computeGradient(data, i, j, k, cellRes, dx, dy, dz);
	  generatePLIC(f, grad, dx[i], dy[j], dz[k], 
		       ox, oy, oz, vertices, indices, vertexID);
	}
      }
    }
  }
  
  std::vector<std::vector<int>> borders(0);
  
  //  extractPLICBorders(vertices, indices, borders);

  const int numPoints = vertices.size();
  const int numTriangles = indices.size()/3;

  vtkPoints *points = vtkPoints::New();
  points->SetNumberOfPoints(numPoints);
  for (int i = 0; i < numPoints; ++i) {
    float p[3] = {vertices[i].x, vertices[i].y, vertices[i].z};
    points->SetPoint(i, p);
  }

  vtkIdTypeArray *cellIndices = vtkIdTypeArray::New();
  cellIndices->SetNumberOfComponents(1);
  cellIndices->SetNumberOfTuples(numTriangles*4);
  for (int i = 0; i < numTriangles; ++i) {
    cellIndices->SetValue(i*4, 3);
    cellIndices->SetValue(i*4+1, indices[i*3+0]);
    cellIndices->SetValue(i*4+2, indices[i*3+1]);
    cellIndices->SetValue(i*4+3, indices[i*3+2]);
  }
  vtkCellArray *cells = vtkCellArray::New();
  cells->SetCells(numTriangles, cellIndices);

  output->SetPoints(points);
  output->SetPolys(cells);

  return 1;
}

////////// External Operators /////////////
void vtkPLICVis::PrintSelf(ostream &os, vtkIndent indent)
{
}
