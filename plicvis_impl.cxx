#include "vtkDataArray.h"
#include <cstdlib>
#include <map>
#include <cmath>
#include <limits>
#include <set>
#include <list>
#include "mc_tables.h"
#include "plicvis_impl.h"

float3 computeGradient(vtkDataArray *data, int cell_i, int cell_j, int cell_k, int res[3],
		       std::vector<float> &dx,std::vector<float> &dy,std::vector<float> &dz)
{
  float3 cellGradient = make_float3(0.0f);
  float dfm1, dfm2;

  // compute gradient on nodes
  for (int k = 0; k < 2; ++k) {
    int km = std::max(0,           k+cell_k-1);
    int kp = std::min(res[2]-1,k+cell_k);
    float dzc = (dz[km] + dz[kp])*0.5f;

    for (int j = 0; j < 2; ++j) {
      int jm = std::max(0,           j+cell_j-1);
      int jp = std::min(res[1]-1,j+cell_j);
      float dyc = (dy[jm] + dy[jp])*0.5f;

      for (int i = 0; i < 2; ++i) {
	int im = std::max(0,           i+cell_i-1);
	int ip = std::min(res[0]-1,i+cell_i);
	float dxc = (dx[im] + dx[ip])*0.5f;

	float f[8] = {data->GetComponent(im + jm*res[0] + km*res[0]*res[1],0),
		      data->GetComponent(ip + jm*res[0] + km*res[0]*res[1],0),
		      data->GetComponent(im + jp*res[0] + km*res[0]*res[1],0),
		      data->GetComponent(ip + jp*res[0] + km*res[0]*res[1],0),
		      data->GetComponent(im + jm*res[0] + kp*res[0]*res[1],0),
		      data->GetComponent(ip + jm*res[0] + kp*res[0]*res[1],0),
		      data->GetComponent(im + jp*res[0] + kp*res[0]*res[1],0),
		      data->GetComponent(ip + jp*res[0] + kp*res[0]*res[1],0)};

	dfm1 = (f[7] - f[6])*dz[km] + (f[3] - f[2])*dz[kp];
	dfm2 = (f[5] - f[4])*dz[km] + (f[1] - f[0])*dz[kp];	    
	float nx = 0.25f*(dfm1*dy[jm]+dfm2*dy[jp]) / (dxc*dyc*dzc);

	dfm1 = (f[7] - f[5])*dz[km] + (f[3] - f[1])*dz[kp];
	dfm2 = (f[6] - f[4])*dz[km] + (f[2] - f[0])*dz[kp];	    
	float ny = 0.25f*(dfm1*dx[im]+dfm2*dx[ip]) / (dxc*dyc*dzc);

	dfm1 = (f[7] - f[3])*dy[jm] + (f[5] -  f[1])*dy[jp];
	dfm2 = (f[6] - f[2])*dy[jm] + (f[4] -  f[0])*dy[jp];	    
	float nz = 0.25f*(dfm1*dx[im]+dfm2*dx[ip]) / (dxc*dyc*dzc);

	cellGradient += make_float3(-nx,-ny,-nz);// normal points from f outwards
      }
    }
  }

  return cellGradient/8.0f;
}
void vertexInterp(float isolevel, float3 p0, float3 p1, float f0, float f1, float3 &p)
{
  float t = (isolevel - f0) / (f1 - f0);
  p = lerp(p0, p1, t);
}

class compare_float3 {
public:
  bool operator()(const float3 a, const float3 b) const {
    return (a.x < b.x || (a.x == b.x && (a.y < b.y || (a.y == b.y && (a.z < b.z)))));
  }
};

static void mergeTriangles(std::vector<float3> &pos, int totalVerts,
			   std::vector<int>& indices, std::vector<float3>& vertices,
			   int &vertexID)
{
  std::map<float3, int, compare_float3> vertexInfo;

  for (int t = 0; t < totalVerts/3; t++) {

    float3 vrt[3] = {pos[t*3+0], pos[t*3+1], pos[t*3+2]};

    for (int v = 0; v < 3; v++) {

      float3 &key = vrt[v];

      if (vertexInfo.find(key) == vertexInfo.end()) {

    	vertexInfo[key] = vertexID;

    	vertices.push_back(pos[t*3+v]);
    	indices.push_back(vertexID);

    	vertexID++;
      }
      else {

      	indices.push_back(vertexInfo[key]);
      }
    }
  }
}

void cross(double a[3], double b[3], double c[3])
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

double dot(double a[3], double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double computeVolume(std::vector<float3> &vertices)
{
  int numCells = vertices.size()/3;
  float volume = 0.0;

  for (int c = 0; c < numCells; c++) {
    int id0 = c*3+0;
    int id1 = c*3+1;
    int id2 = c*3+2;

    double pos[3][3] = {vertices[id0].x, vertices[id0].y, vertices[id0].z,
			vertices[id1].x, vertices[id1].y, vertices[id1].z,
			vertices[id2].x, vertices[id2].y, vertices[id2].z};
    double cr[3];
    cross(pos[1],pos[2],cr);
    volume += dot(pos[0], cr)/6.0f;
  }
  return volume;
}

static int edge[12][2] = { 	// hex
  {0,1},{2,3},{4,5},{6,7},{0,2},{1,3},{4,6},{5,7},{0,4},{1,5},{2,6},{3,7}
};	


void generateSurface(float f, float3 grad, 
		     float dx, float dy, float dz,
		     float ox, float oy, float oz,
		     const float isoValue, 
		     int closed,
		     std::vector<float3> &vertices)
{
  const int fSize = 4;
  int b0, b1;

  std::vector<float> nodeVals(64, -1.0f*std::numeric_limits<float>::max());
  float3 d2 = make_float3(dx/2.0f,dy/2.0f,dz/2.0f);
  nodeVals[21] = f + dot(grad,make_float3(-d2.x,-d2.y,-d2.z));
  nodeVals[22] = f + dot(grad,make_float3( d2.x,-d2.y,-d2.z));
  nodeVals[25] = f + dot(grad,make_float3(-d2.x, d2.y,-d2.z));
  nodeVals[26] = f + dot(grad,make_float3( d2.x, d2.y,-d2.z));
  nodeVals[37] = f + dot(grad,make_float3(-d2.x,-d2.y, d2.z));
  nodeVals[38] = f + dot(grad,make_float3( d2.x,-d2.y, d2.z));
  nodeVals[41] = f + dot(grad,make_float3(-d2.x, d2.y, d2.z));
  nodeVals[42] = f + dot(grad,make_float3( d2.x, d2.y, d2.z));

  float3 p0 = make_float3(ox-dx,oy-dy,oz-dz);

  if (closed) {
    b0 = 0;
    b1 = 3;    
  }
  else {
    b0 = 1;
    b1 = 2;     
  }
  float3 cellCenter = make_float3(ox+dx/2.0f,oy+dy/2.0f,oz+dz/2.0f);

  for (int k = b0; k < b1; k++) {
    for (int j = b0; j < b1; j++) {
      for (int i = b0; i < b1; i++) {
  	int ip = i+1;
  	int jp = j+1;
  	int kp = k+1;
	float3 v[8];

	float3 org = {p0.x+i*dx,p0.y+j*dy,p0.z+k*dz};

	v[0] = make_float3(org.x,   org.y,   org.z);
	v[1] = make_float3(org.x+dx,org.y,   org.z);
	v[2] = make_float3(org.x+dx,org.y+dy,org.z);
	v[3] = make_float3(org.x,   org.y+dy,org.z);
	v[4] = make_float3(org.x,   org.y,   org.z+dz);
	v[5] = make_float3(org.x+dx,org.y,   org.z+dz);
	v[6] = make_float3(org.x+dx,org.y+dy,org.z+dz);
	v[7] = make_float3(org.x,   org.y+dy,org.z+dz);

  	int ids[8] = {i  + j*fSize  + k*fSize*fSize,
		      ip + j*fSize  + k*fSize*fSize,
		      ip + jp*fSize + k*fSize*fSize,
		      i  + jp*fSize + k*fSize*fSize,
		      i  + j*fSize  + kp*fSize*fSize,
		      ip + j*fSize  + kp*fSize*fSize,
		      ip + jp*fSize + kp*fSize*fSize,
		      i  + jp*fSize + kp*fSize*fSize};

	float field[8] = {nodeVals[ids[0]],
			  nodeVals[ids[1]],
			  nodeVals[ids[2]],
			  nodeVals[ids[3]],
			  nodeVals[ids[4]],
			  nodeVals[ids[5]],
			  nodeVals[ids[6]],
			  nodeVals[ids[7]]};

	// calculate flag indicating if each vertex is inside or outside isosurface
	unsigned int cubeindex = 0;
	cubeindex += uint(field[0] < isoValue);
	cubeindex += uint(field[1] < isoValue)*2;
	cubeindex += uint(field[2] < isoValue)*4;
	cubeindex += uint(field[3] < isoValue)*8;
	cubeindex += uint(field[4] < isoValue)*16;
	cubeindex += uint(field[5] < isoValue)*32;
	cubeindex += uint(field[6] < isoValue)*64;
	cubeindex += uint(field[7] < isoValue)*128;

	float3 vertlist[12];
	
	vertexInterp(isoValue, v[0], v[1], field[0], field[1], vertlist[0]);
	vertexInterp(isoValue, v[1], v[2], field[1], field[2], vertlist[1]);
	vertexInterp(isoValue, v[2], v[3], field[2], field[3], vertlist[2]);
	vertexInterp(isoValue, v[3], v[0], field[3], field[0], vertlist[3]);
	
	vertexInterp(isoValue, v[4], v[5], field[4], field[5], vertlist[4]);
	vertexInterp(isoValue, v[5], v[6], field[5], field[6], vertlist[5]);
	vertexInterp(isoValue, v[6], v[7], field[6], field[7], vertlist[6]);
	vertexInterp(isoValue, v[7], v[4], field[7], field[4], vertlist[7]);
	
	vertexInterp(isoValue, v[0], v[4], field[0], field[4], vertlist[8]);
	vertexInterp(isoValue, v[1], v[5], field[1], field[5], vertlist[9]);
	vertexInterp(isoValue, v[2], v[6], field[2], field[6], vertlist[10]);
	vertexInterp(isoValue, v[3], v[7], field[3], field[7], vertlist[11]);

	int numVerts = numVertsTable[cubeindex];

	for (int i=0; i<numVerts; i++) {
	  uint edge = triTable[cubeindex][i];
	  float3 v = vertlist[edge];
	  // todo: set parameter in gui
	  float3 dist = v - cellCenter;
	  v = cellCenter + dist*0.999f;
	  vertices.push_back(v);
	}
      }
    }
  }
}

float computeIsoValue(float f, float3 grad, 
		      float dx, float dy, float dz,
		      float ox, float oy, float oz)
{
  float3 d2 = make_float3(dx/2.0f,dy/2.0f,dz/2.0f);
  float err = dx*dy*dz*0.0001f;
  float minVal = f+dot(-1.0f*fabs(grad),d2)-err;
  float maxVal = f+dot(fabs(grad),d2)+err;
  float isoValue = 0.5;
  float volume;
  int i = 0;
        
  isoValue = (maxVal+minVal)/2.0;
  std::vector<float3> vertices;
  while (i < 10) {
      
    vertices.clear();

    generateSurface(f, grad, dx, dy, dz, ox, oy, oz, 
		    isoValue, 1, vertices);

    volume = fabs(computeVolume(vertices));
    
    if (fabs(volume-f) < err)
	break;

    if (volume < f)
      maxVal = isoValue;
    else
      minVal = isoValue;
    isoValue = (maxVal+minVal)/2.0;

    i++;
  }
  return isoValue;
}

void generatePLIC(float f, float3 grad, 
		  float dx, float dy, float dz,
		  float ox, float oy, float oz, 
		  std::vector<float3> &vertices,
		  std::vector<int> &indices,
		  int &vertexID)
{
  float fScaled = f*dx*dy*dz;
  float isoValue = computeIsoValue(fScaled, -grad, dx, dy, dz, ox, oy, oz);

  std::vector<float3> verticesTmp(0);

  generateSurface(fScaled, -grad, dx, dy, dz, ox,  oy, oz,
		  isoValue, 0, verticesTmp);
  
  mergeTriangles(verticesTmp, verticesTmp.size(), indices, vertices, vertexID);
}

void nodeCoordsToEdgeLengths(vtkDataArray *coords, std::vector<float> &dx)
{
  int cellRes = coords->GetNumberOfTuples()-1;
  dx.resize(cellRes);
  float p0 = coords->GetComponent(0,0);
  for (int i = 0; i < cellRes; ++i) {
    float p1 = coords->GetComponent(i+1,0);
    dx[i] = p1 - p0;
    p0 = p1;
  }
}

int interfaceCell(vtkDataArray *data, int cell_i, int cell_j, int cell_k, int cellRes[3])
{
  int im = std::max(0,cell_i-1);
  int ip = std::min(cellRes[0]-1,cell_i+1);
  int jm = std::max(0,cell_j-1);
  int jp = std::min(cellRes[1]-1,cell_j+1);
  int km = std::max(0,cell_k-1);
  int kp = std::min(cellRes[2]-1,cell_k+1);

  for (int k = km; k <= kp; ++k) {
    for (int j = jm; j <= jp; ++j) {
      for (int i = im; i <= ip; ++i) {
	
	if (i == cell_i && j == cell_j && k == cell_k) 
	  continue;

	int idx = i + j*cellRes[0] + k*cellRes[0]*cellRes[1];
	if (data->GetComponent(idx,0) < EMF1) {
	  return 1;
	}
      }
    }
  }
  return 0;
}

//
//http://stackoverflow.com/questions/14108553/get-border-edges-of-mesh-in-winding-order
//
void extractPLICBorders(std::vector<float3> &vertices,
			std::vector<int> &indices,
			std::vector<std::vector<int>> &borders)
{
  std::map<std::pair<int,int>,int> edges;
  edges.clear();

  const int numTriangles = indices.size()/3;
  for (int i = 0; i < numTriangles; ++i) {
    int ids[3] = {indices[i*3+0], indices[i*3+1], indices[i*3+2]};
    
    for (int j = 0; j < 3; ++j) {

      int e0 = ids[j];//std::min(ids[j],ids[(j+1)%3]);
      int e1 = ids[(j+1)%3];//std::max(ids[j],ids[(j+1)%3]);
      std::pair<int,int> e = std::pair<int,int>(e0,e1);

      if (edges.find(e) != edges.end()) {
	edges[e] += 1;
      }
      else {
	edges[e] = 1;
      }
    }
  }

  std::list<std::pair<int,int>> edge_list;
  edge_list.clear();

  for (const auto &edge : edges) {
    if (edge.second == 1) { // edge with only one traingle attached
      edge_list.push_front(edge.first);
    }
  }

  while (!edge_list.empty()) {
    
    std::pair<int,int> e = edge_list.front();
    edge_list.pop_front();
    const int i_first = e.first;
    int i_next = e.second;

    std::vector<int> border(0);
    border.push_back(i_first);
    border.push_back(i_next);

    std::list<std::pair<int,int>>::iterator it = edge_list.begin();
    while (it != edge_list.end() && i_next != i_first) {
      if (it->first == i_next) {
	border.push_back(it->second);
	i_next = it->second;
	edge_list.erase(it);	
      }
      else if (it->second == i_next) {
	border.push_back(it->first);
	i_next = it->first;
	edge_list.erase(it);
      }
      else {
	++it;
      }
    }
    borders.push_back(border);
  }
}
