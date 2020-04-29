#include "Polygon.h"

bool Polygon::IsIntersection(const Ray& r, IntersectionData& data)
{
	// Check for intersection with every triangle
	data.depth = std::numeric_limits<float>::max();
	IntersectionData polyData;
	bool hasIntersection = false;
	for (int vertexIndex = 2; vertexIndex < mNumberOfVertices; ++vertexIndex)
	{
		// Create a triangle
		Triangle tri(mVertices[0], mVertices[vertexIndex - 1], mVertices[vertexIndex]);
		if (tri.IsIntersection(r, polyData))
		{
			hasIntersection = true;
		}
		if (polyData.depth < data.depth)
		{
			data = polyData;
			data.object = this;
		}
	}
	
	return hasIntersection;
}
