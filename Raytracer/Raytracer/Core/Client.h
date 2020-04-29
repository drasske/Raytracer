// Client.h
#pragma once

#include <thread>
#include <random>

#include <SDL2/SDL.h>

#include "Scene.h"
#include "Bitmap.h"
#include "Distribution.h"

#define EPSILON 0.001f

class Client 
{
public:
	Client(SDL_Window* pWindow);
	~Client();

	void Load(const char* fname);
	
	void RenderStart();
	void RenderStop();

	void Resize();

	void Blit();

	void Save();

private:
	void Draw();
	void RenderScene();
	bool CastRay(const Ray& r, IntersectionData& closestObjectData);
	bool CastRay(const Ray& r);
	glm::vec3 RayTrace(const IntersectionData& d, const Ray& r, float ni, unsigned int depth);
	glm::vec3 PathTrace(const IntersectionData& d, const Ray& r, unsigned int depth);
	float CalculateReflectionCoefficient(const glm::vec3& i, const glm::vec3& N, float ni, float nt, float ui, float ut);
	glm::vec3 LocalLighting(const IntersectionData& d, float R, const glm::vec3& N);
	glm::vec3 GenerateRandomSample(const glm::vec3& N, const glm::vec3& point, float& probability);
	glm::vec3 UniformHemisphereDistribution(float u, float v, const glm::vec3& N);

private:
	SDL_Window* mpWindow;
	Scene mScene;
	Bitmap* mpBitmap;
	std::thread* mpRenderThread;
	bool mKillFlag;
	const float PI = 4.0f * atan(1.0f);
	Combined mCombinedDistribution;
};