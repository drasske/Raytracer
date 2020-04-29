#pragma once

#include <vector>

#include "Object.h"
#include "Triangle.h"

#include "IntersectionData.h"

class Polygon : public Object
{
public:
	bool IsIntersection(const Ray& r, IntersectionData& data) override;
public:
	int mNumberOfVertices;
	std::vector<glm::vec3> mVertices;
	std::vector<Triangle*> mTriangles;
};