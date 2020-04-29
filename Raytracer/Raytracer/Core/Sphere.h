// Sphere.h
#pragma once

#include "Object.h"

#include "IntersectionData.h"

class Sphere : public Object
{
public:
	bool IsIntersection(const Ray& r, IntersectionData& data) override;

public:
	glm::vec3 mCenter = glm::vec3(0.0f, 0.0f, 0.0f);
	float mRadius = 1.0f;
};