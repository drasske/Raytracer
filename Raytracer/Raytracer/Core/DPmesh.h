// DPmesh.h
#pragma once

#include "Object.h"

#include "IntersectionData.h"

class DPmesh : public Object
{
public:
	DPmesh(int n = 20);

	bool IsIntersection(const Ray& r, IntersectionData& data);

	int StemStartIndex() const {return 0;}
	int StemIndexCount() const { return 12; }
	int Bowl1StartIndex() const { return 12; }
	int Bowl1IndexCount() const { return mBowlFaceCount; }
	int Bowl2StartIndex() const { return 12 + mBowlFaceCount; }
	int Bowl2IndexCount() const { return mBowlFaceCount; }

public:
  // Each face is indices into vertices that gives a triangle
	struct Face {
		Face(int i = -1, int j = -1, int k = -1) : index1(i), index2(j), index3(k) {}
		int index1, index2, index3;
	};

public:
	std::vector<glm::vec3> mVertices;
	std::vector<glm::vec3> mNormals;
	std::vector<Face> mFaces;

	int mBowlFaceCount;

	glm::vec3 mCenter = glm::vec3(0.0f, 0.0f, 0.0f);
	glm::vec3 mUAxis = glm::vec3(1.0f, 0.0f, 0.0f);
	glm::vec3 mVAxis = glm::vec3(0.0f, 1.0f, 0.0f);
	glm::vec3 mWAxis = glm::vec3(0.0f, 0.0f, 1.0f);
private:
	const float PI = 4.0f * atan(1.0f);
};