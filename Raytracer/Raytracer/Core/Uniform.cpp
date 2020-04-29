#include "Distribution.h"

float Uniform::probability(const glm::vec3& v)
{
	return glm::dot(v, direction) >= 0 ? 1.0f / (2.0f * PI) : 0.0f;
}

glm::vec3 Uniform::generate()
{
	// Generate two samples on uniform distribution [0, 1]
	float u = mDistribution(mGenerator);
	float v = mDistribution(mGenerator);
	//glm::vec3 Client::UniformHemisphereDistribution(float u, float v, const glm::vec3& N)
	// Compute cartesian coordinates
	float x = sqrtf(u*(2.0f - u))*cosf(2.0f*PI*v);
	float y = sqrtf(u*(2.0f - u))*sinf(2.0f*PI*v);
	float z = 1.0f - u;

	// Apply transformation to get sample about N
	// Choose r vector not parallel to N
	glm::vec3 rVec = glm::vec3(direction.x + 1.0f, direction.y - 1.0f, direction.z + 1.0f);
	// Should put a check here just in case normal points in (1, -1, 1)
	glm::vec3 uVec = glm::normalize(glm::cross(direction, rVec));
	glm::vec3 vVec = glm::cross(direction, uVec);

	glm::vec3 sample = x * uVec + y * vVec + z * direction;

	return sample;
}