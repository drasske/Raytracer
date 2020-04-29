#pragma once

#include <glm/glm.hpp>

#include "Sphere.h"

class Light : public Sphere
{
public:
	glm::vec3 mColor;
};