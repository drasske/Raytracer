#include "Ellipsoid.h"

bool Ellipsoid::IsIntersection(const Ray& r, IntersectionData& data)
{
	// Set default distance to the maximum value
	data.depth = std::numeric_limits<float>::max();

	// Take inverse of model matrix
	glm::mat3 modelMatrix = { mUAxis, mVAxis,mWAxis };
	glm::mat3 modelMatrixInverse = glm::inverse(modelMatrix);

	// Transform ray into object space to do unit sphere intersection
	// R'(t) = M-1(P0 - C) + tM-1d0
	glm::vec3 rayDirection = modelMatrixInverse * r.direction;
	glm::vec3 ellipsoidRay = modelMatrixInverse * (r.position - mCenter);
	float centerOffset = glm::length(ellipsoidRay);

	float a = sq(glm::length(rayDirection));
	float b = 2.0f * glm::dot(rayDirection, ellipsoidRay);
	float c = sq(centerOffset) - 1.0f;

	float discriminant = sq(b) - (4.0f * a * c);
	if (discriminant < 0.0f)
	{
		return false;
	}

	float tp = (-b + sqrtf(discriminant)) / (2.0f * a);
	float tn = (-b - sqrtf(discriminant)) / (2.0f * a);
	float intersectionDepth = 0.0f;

	if (tp >= 0.0f)
	{
		intersectionDepth = (tn < 0.0f) ? tp : tn;
	}
	else
	{
		return false;
	}

	// Use the depth to get the point of intersection on the unit sphere
	glm::vec3 intersectionPointUnitSphere = ellipsoidRay + intersectionDepth * rayDirection;
	glm::vec3 sphereNormal = intersectionPointUnitSphere;

	// Use the point of intersection to get the ellipsoid surface normal
	glm::vec3 n = glm::transpose(modelMatrixInverse) * sphereNormal;

	// Set intersection data
	data.depth = intersectionDepth;
	data.object = this;
	data.normal = glm::normalize(n);
	data.view = glm::normalize(-r.direction);
	data.point = r.position + intersectionDepth * r.direction;

	return true;
}