#include "Box.h"

bool Box::IsIntersection(const Ray& r, IntersectionData& data)
{
	data.depth = std::numeric_limits<float>::max();
	// TODO: Move this into constructor
	// Generate half space points and normals
	glm::vec3 points[6];
	glm::vec3 normals[6];

	// Front
	points[0] = mCenter + mWVector;
	normals[0] = glm::cross(mUVector, mVVector);

	// Back 
	points[1] = mCenter - mWVector;
	normals[1] = -normals[0];

	// Left
	points[2] = mCenter + mUVector;
	normals[2] = glm::cross(mVVector, mWVector);

	// Right
	points[3] = mCenter - mUVector;
	normals[3] = -normals[2];

	// Top
	points[4] = mCenter + mVVector;
	normals[4] = glm::cross(mWVector, mUVector);

	// Bottom
	points[5] = mCenter - mVVector;
	normals[5] = -normals[4];

	for (int i = 0; i < 6; ++i)
	{
		normals[i] = glm::normalize(normals[i]);
	}

	// Determine intersection interval
	float intervalMin = 0.0f;
	float intervalMax = std::numeric_limits<float>::max();
	unsigned int numberOfFaces = 6;
	unsigned int faceIndex = 0;
	float intersectionT = 0.0f;
	glm::vec3 intersectionNormalMin = glm::vec3(0.0f, 0.0f, 0.0f);
	glm::vec3 intersectionNormalMax = glm::vec3(0.0f, 0.0f, 0.0f);

	while ((faceIndex < numberOfFaces) && (intervalMin <= intervalMax))
	{
		glm::vec3 point = points[faceIndex];
		glm::vec3 normal = normals[faceIndex];
		float ddotn = glm::dot(r.direction, normal);
		if (ddotn < 0.0f)
		{
			intersectionT = -(glm::dot((r.position - point), normal) / ddotn);
			if (intervalMin < intersectionT)
			{
				intersectionNormalMin = normal;
			}
			intervalMin = std::fmax(intersectionT, intervalMin);
		}
		else if (ddotn > 0.0f)
		{
			intersectionT = -(glm::dot((r.position - point), normal) / ddotn);
			if (intervalMax > intersectionT)
			{
				intersectionNormalMax = normal;
			}
			intervalMax = std::fmin(intersectionT, intervalMax);
		}
		else if(glm::dot((r.position - point), normal) > 0.0f)
		{
			// No intersection, exit loop
			intervalMin = intervalMax + 1.0f;
		}
		++faceIndex;
	}

	glm::vec3 intersectionNormal = glm::vec3(0.0f, 0.0f, 0.0f);

	if (intervalMin > intervalMax)
	{
		return false;
	}
	else if (intervalMin == 0.0f)
	{
		intersectionT = intervalMax;
		intersectionNormal = intersectionNormalMax;
	}
	else
	{
		intersectionT = intervalMin;
		intersectionNormal = intersectionNormalMin;
	}

	// Populate intersection data
	data.depth = intersectionT;
	data.normal = intersectionNormal;
	data.object = this;
	data.point = r.position + intersectionT * r.direction;
	data.view = glm::normalize(-r.direction);

	return true;
}
