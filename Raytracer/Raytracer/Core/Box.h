// Box.h
#pragma once

#include "Object.h"

#include "IntersectionData.h"

class Box : public Object
{
public:
	bool IsIntersection(const Ray& r, IntersectionData& data) override;
public:
	glm::vec3 mCenter = glm::vec3(0.0f, 0.0f, 0.0f);
	glm::vec3 mUVector = glm::vec3(1.0f, 0.0f, 0.0f);
	glm::vec3 mVVector = glm::vec3(0.0f, 1.0f, 0.0f);
	glm::vec3 mWVector = glm::vec3(0.0f, 0.0f, 1.0f);
};