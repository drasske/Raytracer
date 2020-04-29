#include "Distribution.h"
#include <iostream>

std::random_device Distribution::mRandomDevice;
std::mt19937 Distribution::mGenerator(mRandomDevice());
std::uniform_real_distribution<float> Distribution::mDistribution(0.0, 1.0);

float Beckmann::probability(const glm::vec3& v)
{
	float cosTheta = glm::dot(direction, v);
	if (cosTheta <= 0.0f)
	{
		return 0.0f;
	}

	float theta = acosf(cosTheta);
	float p =	expf(-(tanf(theta) * tanf(theta)) / (roughness * roughness)) /
				(PI * roughness * roughness * cosTheta * cosTheta * cosTheta);
	
	return p;
}

glm::vec3 Beckmann::generate()
{
	// Generate two samples on uniform distribution [0, 1]
	float u = mDistribution(mGenerator);
	float v = mDistribution(mGenerator);

	float tau = -(roughness * roughness) * logf(1.0f - u);
	float sqrtTau = sqrtf(tau / (1.0f + tau));

	float x = sqrtTau*cosf(2.0f*PI*v);
	float y = sqrtTau*sinf(2.0f*PI*v);
	float z = 1.0f / sqrtf(1.0f + tau);

	// Apply transformation to get sample about L
	// Choose r vector not parallel to L
	glm::vec3 rVec = glm::vec3(direction.x + 1.0f, direction.y - 1.0f, direction.z + 1.0f);
	// Should put a check here just in case direction points in (1, -1, 1)
	glm::vec3 uVec = glm::normalize(glm::cross(direction, rVec));
	glm::vec3 vVec = glm::cross(direction, uVec);

	glm::vec3 sample = x * uVec + y * vVec + z * direction;

	return sample;
}