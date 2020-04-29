// Ellipsoid.h
#pragma once

#include "Object.h"

#include "IntersectionData.h"

class Ellipsoid : public Object
{
public:
	bool IsIntersection(const Ray& r, IntersectionData& data) override;

public:
	glm::vec3 mCenter = glm::vec3(0.0f, 0.0f, 0.0f);
	glm::vec3 mUAxis = glm::vec3(1.0f, 0.0f, 0.0f);
	glm::vec3 mVAxis = glm::vec3(0.0f, 1.0f, 0.0f);
	glm::vec3 mWAxis = glm::vec3(0.0f, 0.0f, 1.0f);
};