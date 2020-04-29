#include "DPmesh.h"

#include <cmath>

#include "Box.h"
#include "Triangle.h"

DPmesh::DPmesh(int n)
    : mBowlFaceCount(4*n+4) {
  mVertices.resize(48+8*n);
  mNormals.resize(48+8*n);
  mFaces.resize(20+8*n);

  // d/p common stem  (v/n: 0 .. 23, f: 0 .. 11)
  const glm::vec3 S0(0.1f,1,1),
                  S1(0.1f,-1,1),
                  S2(-0.1f,-1,1),
                  S3(-0.1f,1,1),
                  S4(0.1f,1,-1),
                  S5(0.1f,-1,-1),
                  S6(-0.1f,-1,-1),
                  S7(-0.1f,1,-1),
                  EX(1,0,0),
                  EY(0,1,0),
                  EZ(0,0,1);
  mVertices[0] = S0;     mNormals[0] = EZ;     mFaces[0] = Face(0,2,1);
  mVertices[1] = S1;     mNormals[1] = EZ;     mFaces[1] = Face(0,3,2);
  mVertices[2] = S2;     mNormals[2] = EZ;
  mVertices[3] = S3;     mNormals[3] = EZ;
  mVertices[4] = S4;     mNormals[4] = -EZ;    mFaces[2] = Face(4,5,6);
  mVertices[5] = S5;     mNormals[5] = -EZ;    mFaces[3] = Face(4,6,7);
  mVertices[6] = S6;     mNormals[6] = -EZ;
  mVertices[7] = S7;     mNormals[7] = -EZ;

  mVertices[8] = S0;     mNormals[8] = EX;     mFaces[4] = Face(8,9,10);
  mVertices[9] = S1;     mNormals[9] = EX;     mFaces[5] = Face(8,10,11);
  mVertices[10] = S5;    mNormals[10] = EX;
  mVertices[11] = S4;    mNormals[11] = EX;
  mVertices[12] = S3;    mNormals[12] = -EX;   mFaces[6] = Face(12,15,14);
  mVertices[13] = S2;    mNormals[13] = -EX;   mFaces[7] = Face(12,14,13);
  mVertices[14] = S6;    mNormals[14] = -EX;
  mVertices[15] = S7;    mNormals[15] = -EX;

  mVertices[16] = S0;    mNormals[16] = EY;    mFaces[8] = Face(16,17,18);
  mVertices[17] = S4;    mNormals[17] = EY;    mFaces[9] = Face(16,18,19);
  mVertices[18] = S7;    mNormals[18] = EY;
  mVertices[19] = S3;    mNormals[19] = EY;
  mVertices[20] = S1;    mNormals[20] = -EY;   mFaces[10] = Face(20,23,22);
  mVertices[21] = S5;    mNormals[21] = -EY;   mFaces[11] = Face(20,22,21);
  mVertices[22] = S6;    mNormals[22] = -EY;
  mVertices[23] = S2;    mNormals[23] = -EY;

  // p bowl (v/n: 24 .. 35+4n, f: 12 .. 16+4n)
  const glm::vec3 B0(0.2f,1,1),
                  B1(0.2f,0,1),
                  B2(0.2f,0,-1),
                  B3(0.2f,1,-1),
                  B4(0.44f,1,1),
                  B5(0.44f,0,1),
                  B6(0.44f,0,-1),
                  B7(0.44f,1,-1);
  mVertices[24] = B0;    mNormals[24] = -EX;   mFaces[12] = Face(24,26,25);
  mVertices[25] = B1;    mNormals[25] = -EX;   mFaces[13] = Face(24,27,26);
  mVertices[26] = B2;    mNormals[26] = -EX;
  mVertices[27] = B3;    mNormals[27] = -EX;

  mVertices[28] = B1;  mNormals[28] = EZ;  // v/n: 28 .. 29+n
  for (int i=0; i < n; ++i) {
    float t = PI*(i/float(n-1)-0.5f);
    mVertices[29+i] = 0.5f*(B4 + B5) + 0.5f*cos(t)*EX + 0.5f*sin(t)*EY;
    mNormals[29+i] = EZ;
  }
  mVertices[29+n] = B0;  mNormals[29+n] = EZ;
  for (int i=0; i < n; ++i)              // f: 14 .. 13+n
    mFaces[14+i] = Face(28,29+i,30+i);

  for (int i=0; i < n+2; ++i) {   // v/n: 30+n .. 31+2n
    mVertices[30+n+i] = mVertices[28+i] - 2.0f*EZ;
    mNormals[30+n+i] = -EZ;
  }
  for (int i=0; i < n; ++i)     // f: 14+n .. 13+2n
    mFaces[14+n+i] = Face(30+n,31+n+i,32+n+i);

  for (int i=0; i < n+2; ++i) {  // v/n: 32+2n .. 33+3n & 34+3n .. 35+4n
    mVertices[32+2*n+i] = mVertices[28+i];
    mVertices[34+3*n+i] = mVertices[30+n+i];
  }
  mNormals[32+2*n] = mNormals[34+3*n] = -EY;
  for (int i=0; i < n; ++i) {
    float t = PI*(i/float(n-1)-0.5f);
    mNormals[33+2*n+i] = mNormals[35+3*n+i]
                      = cos(t)*EX + sin(t)*EY;
  }
  mNormals[33+3*n] = mNormals[35+4*n] = EY;
  for (int i=0; i < n+1; ++i) {  // f: 14+2n .. 15+4n
    mFaces[14+2*n+2*i+0] = Face(32+2*n+i,33+2*n+i,35+3*n+i);
    mFaces[14+2*n+2*i+1] = Face(32+2*n+i,35+3*n+i,34+3*n+i);
  }

  // d bowl
  for (int i=0; i < 4*n+12; ++i) {  // v/n: 36+4n .. 47+8n
    mVertices[36+4*n+i] = glm::vec3(-1,-1,1) * mVertices[24+i];
    mNormals[36+4*n+i] = glm::vec3(-1,-1,1) * mNormals[24+i];
  }
  for (int i=0; i < 4*n+4; ++i) { // f: 16+4n .. 19+8n
    const Face &f = mFaces[12+i];
    mFaces[16+4*n+i] = Face(f.index1+4*n+12,f.index2+4*n+12,f.index3+4*n+12);
  }

}

bool DPmesh::IsIntersection(const Ray& r, IntersectionData& data)
{
	// Check for intersection with standard cube
	data.depth = std::numeric_limits<float>::max();

	// Take inverse of model matrix
	glm::mat3 modelMatrix = { mUAxis, mVAxis, mWAxis };
	glm::mat3 modelMatrixInverse = glm::inverse(modelMatrix);

	// Transform ray into object space to do unit sphere intersection
	// R'(t) = M-1(P0 - C) + tM-1d0
	glm::vec3 rayDirection = modelMatrixInverse * r.direction;
	glm::vec3 cubeRay = modelMatrixInverse * (r.position - mCenter);
	Ray objectSpaceRay{ cubeRay, rayDirection };

	// Check for intersection with unit cube
	Box unitCube;
	IntersectionData cubeIntersection;
	if (!unitCube.IsIntersection(objectSpaceRay, cubeIntersection))
	{
		return false;
	}

	// Check for intersection with all the triangles in the mesh
	bool hasIntersection = false;
	IntersectionData faceIntersection = IntersectionData();
	glm::vec3 barycentricCoordinates(0.0f);
	for (Face face : mFaces)
	{
		// Create a triangle for each face to test for intersection
		Triangle tri(mVertices[face.index1], mVertices[face.index2], mVertices[face.index3]);
		// Use transformed ray
		if (tri.IsIntersection(objectSpaceRay, faceIntersection, &barycentricCoordinates))
		{
			hasIntersection = true;
			if (faceIntersection.depth < data.depth)
			{
				data.depth = faceIntersection.depth;
				// Interpolate normal
				glm::vec3 interpolatedNormal = glm::transpose(glm::inverse(modelMatrix)) * (mNormals[face.index1] * barycentricCoordinates.x + mNormals[face.index2] * barycentricCoordinates.z + mNormals[face.index3] * barycentricCoordinates.y);
				data.normal = glm::normalize(interpolatedNormal);
				data.object = this;
				data.view = glm::normalize(-r.direction);
				data.point = r.position + data.depth * r.direction;
			}
		}
	}

	return hasIntersection;
}
