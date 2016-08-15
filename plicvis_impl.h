#ifndef PLICVIS_IMPL_H
#define PLICVIS_IMPL_H

#include "vtkDataArray.h"
#include <vector>

#define EMF0 0.000001f
#define EMF1 0.999999f

//----------------------------------------------------------------------------
// float3 structure and related functions taken from NVIDIA CUDA SDK
//----------------------------------------------------------------------------
struct float3
{
    float x, y, z;
};

inline float3 make_float3(float x, float y, float z)
{
  return {x,y,z};
}

inline float3 make_float3(float s)
{
    return make_float3(s, s, s);
}

inline void operator+=(float3 &a, float3 b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
}

inline float3 operator/(float3 a, float b)
{
    return make_float3(a.x / b, a.y / b, a.z / b);
}

inline void operator/=(float3 &a, float b)
{
    a.x /= b;
    a.y /= b;
    a.z /= b;
}

inline float3 operator+(float3 a, float3 b)
{
    return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline float3 operator-(float3 a, float3 b)
{
    return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline float3 operator-(float3 &a)
{
    return make_float3(-a.x, -a.y, -a.z);
}

inline float3 operator*(float b, float3 a)
{
    return make_float3(b * a.x, b * a.y, b * a.z);
}

inline float3 operator*(float3 a, float b)
{
    return make_float3(a.x * b, a.y * b, a.z * b);
}

inline float dot(float3 a, float3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline float3 lerp(float3 a, float3 b, float t)
{
    return a + t*(b-a);
}

inline float3 fabs(float3 v)
{
    return make_float3(fabs(v.x), fabs(v.y), fabs(v.z));
}

inline float length(float3 v)
{
    return sqrtf(dot(v, v));
}

inline float3 normalize(float3 v)
{
  float len = length(v);
  if (len == 0.0f) {
    return make_float3(0.0f);
  }
  return v/len;
}

inline float3 cross(float3 a, float3 b)
{
  return make_float3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

inline float clamp(float f, float a, float b)
{
  return std::max(a, std::min(f, b));
}

//----------------------------------------------------------------------------

float3 computeGradient(vtkDataArray *data, int cell_i, int cell_j, int cell_k, int cellRes[3],
		       std::vector<float> &dx,std::vector<float> &dy,std::vector<float> &dz);

void generatePLIC(float f, float3 grad, 
		  float dx, float dy, float dz,
		  float ox, float oy, float oz, 
		  std::vector<float3> &vertices,
		  std::vector<int> &indices,
		  int &vertexID);

void nodeCoordsToEdgeLengths(vtkDataArray *coords, std::vector<float> &dx);

int interfaceCell(vtkDataArray *data, int cell_i, int cell_j, int cell_k, int cellRes[3]);

void extractPLICBorders(std::vector<float3> &vertices,
			std::vector<int> &indices,
			std::vector<std::vector<int>> &borders);

#endif//PLICVIS_IMPL_H
