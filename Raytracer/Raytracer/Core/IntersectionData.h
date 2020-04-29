// IntersectionData.h
#pragma once

#include <glm/glm.hpp>

#include "Object.h"

struct IntersectionData 
{
	float depth;
	glm::vec3 point;
	glm::vec3 normal;
	glm::vec3 view;
	Object* object;
};