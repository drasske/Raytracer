// Camera.h
#pragma once

#include <glm/glm.hpp>

#include "Ray.h"

class Camera 
{
public:
	Ray GetRay(float x, float y) const;

public:
	glm::vec3 mCenter	= glm::vec3(0.0f, 0.0f, 0.0f);
	glm::vec3 mRight	= glm::vec3(1.0f, 0.0f, 0.0f);
	glm::vec3 mUp		= glm::vec3(0.0f, 1.0f, 0.0f);
	glm::vec3 mEye		= glm::vec3(0.0f, 0.0f, 1.0f);
};