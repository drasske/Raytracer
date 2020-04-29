#include "Sphere.h"

bool Sphere::IsIntersection(const Ray& r, IntersectionData& data)
{
	// Set default depth to max in case of no intersection
	data.depth = std::numeric_limits<float>::max();
	
	// Check for sphere ray intersection
	glm::vec3 sphereRay = (r.position - mCenter);
	float positionCenter = glm::length(sphereRay);

	float a = sq(glm::length(r.direction));
	float b = 2.0f * glm::dot(r.direction, sphereRay);
	float c = sq(positionCenter) - sq(mRadius);

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

	// Use the depth to get the point of intersection
	glm::vec3 intersectionPoint = r.position + intersectionDepth * r.direction;

	// Use the point of intersection to get the sphere surface normal
	glm::vec3 n = intersectionPoint - mCenter;

	// Set intersection data
	data.depth = intersectionDepth;
	data.object = this;
	data.normal = glm::normalize(n);
	data.view = glm::normalize(-r.direction);
	data.point = intersectionPoint;

	return true;
}