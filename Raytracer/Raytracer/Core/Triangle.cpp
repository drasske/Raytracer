#include "Triangle.h"

bool Triangle::IsIntersection(const Ray& r, IntersectionData& data)
{
	return IsIntersection(r, data, nullptr);
}

bool Triangle::IsIntersection(const Ray& r, IntersectionData& data, glm::vec3* coordinates)
{
	data.depth = std::numeric_limits<float>::max();

	// Define barycentric coordiantes of triangle
	glm::vec3 C = mVertices[0];
	glm::vec3 a = mVertices[2] - mVertices[0];
	glm::vec3 b = mVertices[1] - mVertices[0];

	// Get plane triangle lies on
	glm::vec3 n = glm::cross(b, a);

	glm::vec3 triangleRay = r.position - C;
	float ddotn = glm::dot(r.direction, n);

	float tIntersection = -(glm::dot(triangleRay, n) / ddotn);

	if (tIntersection < 0.0f)
	{
		return false;
	}

	// Get point of intersection on plane
	glm::vec3 intersectionPoint = r.position + tIntersection * r.direction;
	const glm::vec3& P = intersectionPoint;

	// Check if plane itersection point lies on plane
	float adotb = glm::dot(a, b);
	float alength = glm::length(a) * glm::length(a);
	float blength = glm::length(b) * glm::length(b);

	float scaleDet = 1.0f / (alength * blength - sq(adotb));

	glm::mat2 determinant{ blength, -adotb, -adotb, alength };
	determinant *= scaleDet;
	glm::vec2 position{ glm::dot(intersectionPoint - C, a), glm::dot(intersectionPoint - C, b) };
	glm::vec2 barycentricCoordinates = determinant * position;

	float mu = barycentricCoordinates.x;
	float nu = barycentricCoordinates.y;
	float lambda = 1.0f - mu - nu;

	if (lambda < 0.0f || mu < 0.0f || nu < 0.0f)
	{
		return false;
	}

	// Actual intersection, store intersection data
	data.depth = tIntersection;
	data.normal = n;
	data.object = this;
	data.point = intersectionPoint;
	data.view = -r.direction;

	// Set barycentric coordinates
	if (coordinates)
	{
		coordinates->x = lambda;
		coordinates->y = mu;
		coordinates->z = nu;
	}

	return true;
}
