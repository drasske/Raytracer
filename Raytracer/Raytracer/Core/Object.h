// Object.h
#pragma once

#include <iostream>
#include <vector>
#include <glm/glm.hpp>

#include "Ray.h"

#define sq(x) x*x

#define Log(x) std::cout << x << std::endl;

struct IntersectionData;

class Object
{
public:
	virtual ~Object() {}
	virtual bool IsIntersection(const Ray& r, IntersectionData& data) = 0;

public:
	glm::vec3 mDiffuseCoefficient = glm::vec3(0.0f);
	float mSpecularCoefficient = 0.0f;
	float mSpecularExponent = 1.0f;
	glm::vec3 mAttenuiationCoefficient = glm::vec3(0.0f);
	float mPermittivity = 0.0f;
	float mPermeability = 0.0f;
	bool mIsLight = false;

};