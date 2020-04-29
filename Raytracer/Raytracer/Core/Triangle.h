#pragma once

#include "Object.h"

#include "IntersectionData.h"

class Triangle : public Object
{
public:
	Triangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
	{
		mVertices[0] = v0;
		mVertices[1] = v1;
		mVertices[2] = v2;
	};
	bool IsIntersection(const Ray& r, IntersectionData& data) override;
	bool IsIntersection(const Ray& r, IntersectionData& data, glm::vec3* coordinates);
public:
	glm::vec3 mVertices[3];
};