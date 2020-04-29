// Scene.h
#pragma once

#include <vector>
#include <iostream>

#include "Object.h"
#include "Camera.h"
#include "Light.h"

struct Air
{
	float mPermittivity = 1.0f;
	float mPermeability = 1.0f;
	glm::vec3 mAttenuationFactor = glm::vec3(1.0f, 1.0f, 1.0f);
};

class Scene 
{
public:
	void Clear();
	void LoadFile(const char* fileName);

public:
	std::vector<Object*> mObjects;
	std::vector<Light*> mLights;
	glm::vec3 mAmbientLight = glm::vec3(0.0f, 0.0f, 0.0f);
	Camera mCamera;
	Air mAir;

private:
	// File parsing 
	bool ParseTriple(std::istream& in, float& x, float& y, float& z);
	bool ParseTriple(std::istream& in, glm::vec3& vec);
	bool ParseFloat(std::istream& in, float& x);
	bool ParseInt(std::istream& in, int& i);
	bool ParseMaterial(std::istream& in, Object& o);
	void ParseCamera(std::istream& in);
	void ParseSphere(std::istream& in);
	void ParseEllipsoid(std::istream& in);
	void ParseBox(std::istream& in);
	void ParsePolygon(std::istream& in);
	void ParseDPLogo(std::istream& in);
	void ParseLight(std::istream& in);
	void ParseAmbient(std::istream& in);
	void ParseAir(std::istream& in);
};