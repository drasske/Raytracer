#include "Scene.h"

#include <fstream>

#include "Sphere.h"
#include "Ellipsoid.h"
#include "Box.h"
#include "Polygon.h"
#include "DPmesh.h"

// Deletes all objects from the scene
void Scene::Clear()
{
	for (Object* pObject : mObjects)
	{
		if (!pObject->mIsLight)
		{
			delete pObject;
		}
	}
	mObjects.clear();
	for (Light* pLight : mLights)
	{
		delete pLight;
	}
	mLights.clear();
	mAir.mAttenuationFactor = glm::vec3(1.0f, 1.0f, 1.0f);
	mAir.mPermeability = 1.0f;
	mAir.mPermittivity = 1.0f;
	mAmbientLight = glm::vec3(0.0f);
}

// Loads a .txt file into the scene
void Scene::LoadFile(const char* fileName)
{
	std::ifstream in(fileName);
	mAmbientLight = glm::vec3(0.0f, 0.0f, 0.0f);
	while (in)
	{
		std::string start;
		in >> start;
		// Ignore comments in the file
		if (start[0] == '#')
		{
			in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}
		else
		{
			if (start == "CAMERA")
			{
				ParseCamera(in);
			}
			else if (start == "SPHERE")
			{
				ParseSphere(in);
			}
			else if (start == "ELLIPSOID")
			{
				ParseEllipsoid(in);
			}
			else if (start == "BOX")
			{
				ParseBox(in);
			}
			else if (start == "POLYGON")
			{
				ParsePolygon(in);
			}
			else if (start == "DPLOGO")
			{
				ParseDPLogo(in);
			}
			// TODO: Add more shapes
			else if (start == "LIGHT")
			{
				ParseLight(in);
			}
			else if (start == "AMBIENT")
			{
				ParseAmbient(in);
			}
			else if (start == "AIR")
			{
				ParseAir(in);
			}
			else if (in)
			{
				#ifdef VERBOSE
				std::cout << "Unrecognized symbol: " << start << std::endl;
				#endif
				in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			}
		}
	}
}

bool Scene::ParseTriple(std::istream& in, float& x, float& y, float& z)
{
	//   ( x , y , z )
	char lp, c1, c2, rp;
	in >> lp >> x >> c1 >> y >> c2 >> z >> rp;
	return bool(in) && lp == '(' && c1 == ',' && c2 == ',' && rp == ')';
}

bool Scene::ParseTriple(std::istream& in, glm::vec3& vec)
{
	return ParseTriple(in, vec.x, vec.y, vec.z);
}

bool Scene::ParseFloat(std::istream& in, float& x)
{
	in >> x;
	return bool(in);
}

bool Scene::ParseInt(std::istream& in, int& i)
{
	in >> i;
	return bool(in);
}

bool Scene::ParseMaterial(std::istream& in, Object& o)
{
	bool okay = ParseTriple(in, o.mDiffuseCoefficient)
			&&	ParseFloat(in, o.mSpecularCoefficient)
			&&	ParseFloat(in, o.mSpecularExponent)
			&&	ParseTriple(in, o.mAttenuiationCoefficient)
			&&	ParseFloat(in, o.mPermittivity)
			&&	ParseFloat(in, o.mPermeability);	

	return okay;
}

void Scene::ParseCamera(std::istream& in)
{
	Camera cam;
	bool okay =	ParseTriple(in, cam.mCenter)
			&&	ParseTriple(in, cam.mRight)
			&&	ParseTriple(in, cam.mUp)
			&&	ParseTriple(in, cam.mEye);
	if (!okay)
	{
		throw std::runtime_error("Failed to parse camera");
	}
	mCamera = cam;
	#ifdef VERBOSE
	std::cout << "Camera parsed" << std::endl;
	#endif
}

void Scene::ParseSphere(std::istream& in)
{
	Sphere* pSphere = new Sphere();
	bool okay = ParseTriple(in, pSphere->mCenter)
			&& ParseFloat(in, pSphere->mRadius)
			&& ParseMaterial(in, *pSphere);
	if (!okay)
	{
		throw std::runtime_error("Failed to parse sphere");
	}
	#ifdef VERBOSE
	std::cout << "Sphere parsed" << std::endl;
	#endif
	mObjects.push_back(pSphere);
}

void Scene::ParseEllipsoid(std::istream& in)
{
	Ellipsoid* pEllipsoid = new Ellipsoid();
	bool okay = ParseTriple(in, pEllipsoid->mCenter)
		&& ParseTriple(in, pEllipsoid->mUAxis)
		&& ParseTriple(in, pEllipsoid->mVAxis)
		&& ParseTriple(in, pEllipsoid->mWAxis)
		&& ParseMaterial(in, *pEllipsoid);
	if (!okay)
	{
		throw std::runtime_error("Failed to parse ellipsoid");
	}
	#ifdef VERBOSE
	std::cout << "Ellipsoid parsed" << std::endl;
	#endif
	mObjects.push_back(pEllipsoid);
}

void Scene::ParseBox(std::istream& in)
{
	Box* pBox = new Box();
	bool okay = ParseTriple(in, pBox->mCenter)
			&&	ParseTriple(in, pBox->mUVector)
			&&	ParseTriple(in, pBox->mVVector)
			&&	ParseTriple(in, pBox->mWVector)
			&&	ParseMaterial(in, *pBox);
	if (!okay)
	{
		throw std::runtime_error("Failed to parse box");
	}
#ifdef VERBOSE
	std::cout << "Box parsed" << std::endl;
#endif
	mObjects.push_back(pBox);
}

void Scene::ParsePolygon(std::istream& in)
{
	Polygon* pPolygon = new Polygon();
	// Parse vertices, sent as triples
	bool okay = ParseInt(in, pPolygon->mNumberOfVertices);
	for (int vertexIndex = 0; vertexIndex < pPolygon->mNumberOfVertices; ++vertexIndex)
	{
		glm::vec3 vertex;
		okay = okay && ParseTriple(in, vertex);
		if (okay)
		{
			pPolygon->mVertices.push_back(vertex);
		}
	}
	okay = okay && ParseMaterial(in, *pPolygon);
	if (!okay)
	{
		throw std::runtime_error("Failed to parse polygon");
	}
#ifdef VERBOSE
	std::cout << "Polygon parsed" << std::endl;
#endif
	mObjects.push_back(pPolygon);
}

void Scene::ParseDPLogo(std::istream& in)
{
	DPmesh* pDPMesh = new DPmesh();
	bool okay = ParseTriple(in, pDPMesh->mCenter)
		&& ParseTriple(in, pDPMesh->mUAxis)
		&& ParseTriple(in, pDPMesh->mVAxis)
		&& ParseTriple(in, pDPMesh->mWAxis)
		&& ParseMaterial(in, *pDPMesh);
	if (!okay)
	{
		throw std::runtime_error("Failed to parse DPLogo");
	}
#ifdef VERBOSE
	std::cout << "DPLogo parsed" << std::endl;
#endif
	mObjects.push_back(pDPMesh);
}

void Scene::ParseLight(std::istream& in)
{
	Light* pLight = new Light;
	bool okay = ParseTriple(in, pLight->mCenter)
		&& ParseTriple(in, pLight->mColor)
		&& ParseFloat(in, pLight->mRadius);
	pLight->mIsLight = true;
	pLight->mDiffuseCoefficient = pLight->mColor;
	if (!okay)
	{
		throw std::runtime_error("Failed to parse Light");
	}
#ifdef VERBOSE
	std::cout << "Light parsed" << std::endl;
#endif
	mLights.push_back(pLight);
	mObjects.push_back(pLight);
}

void Scene::ParseAmbient(std::istream& in)
{
	bool okay = ParseTriple(in, mAmbientLight);
	if (!okay)
	{
		throw std::runtime_error("Failed to parse Ambient");
	}
#ifdef VERBOSE
	std::cout << "Ambient parsed" << std::endl;
#endif
}

void Scene::ParseAir(std::istream& in)
{
	bool okay = ParseFloat(in, mAir.mPermittivity)
		&& ParseFloat(in, mAir.mPermeability)
		&& ParseTriple(in, mAir.mAttenuationFactor);
	if (!okay)
	{
		throw std::runtime_error("Failed to parse Air");
	}
#ifdef VERBOSE
	std::cout << "Air parsed" << std::endl;
#endif
}
