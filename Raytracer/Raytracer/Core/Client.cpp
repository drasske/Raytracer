#include "Client.h"

#include <iostream>
#include <fstream>

#include <GL/glew.h>

#include "IntersectionData.h"
#include "Distribution.h"

// Uncomment to apply supersampling anti-aliasing
//#define ANTIALIAS
#define SAMPLESPERPIXEL 10

// Client constructor, needs an SDL window
Client::Client(SDL_Window* pWindow) : mpRenderThread(nullptr), mKillFlag(false), mpWindow(pWindow)
{
	// Create original bitmap
	mpBitmap = new Bitmap(500, 500);
	for (int i = 0; i < mpBitmap->Height(); ++i)
	{
		for (int j = 0; j < mpBitmap->Width(); ++j)
		{
			mpBitmap->SetPixel(i, j, 100, 0, 100);
		}
	}
}

// Client destructor, clears scene and deletes bitmap
Client::~Client()
{
	RenderStop();
	mScene.Clear();
	delete mpBitmap;
}

// Client load, loads a file give filename
void Client::Load(const char* fname) 
{
	mScene.Clear();
	try 
	{
		mScene.LoadFile(fname);
	}
	catch (std::exception& e) 
	{
		std::cout << e.what() << std::endl;
	}
}

// Draws the scene to the screen as a bitmap
void Client::Draw() 
{
	SDL_SetWindowTitle(mpWindow, "Rendering");
	RenderScene();
	SDL_SetWindowTitle(mpWindow, "Finished Rendering");
}

void Client::RenderScene()
{
	bool defkill = false;
	bool* killptr = (&mKillFlag == nullptr) ? &defkill : &mKillFlag;
	bool& kill = *killptr;
	for (int j = 0; !kill && j < mpBitmap->Height(); ++j)
	{
		for (int i = 0; !kill && i < mpBitmap->Width(); ++i) 
		{
			// Convert to NDC
			float x = 2.0f * (i + 0.5f) / mpBitmap->Width() - 1.0f;
			float y = 2.0f * (j + 0.5f) / mpBitmap->Height() - 1.0f;
			
			glm::vec3 color(0.0f, 0.0f, 0.0f);

			Ray r = mScene.mCamera.GetRay(x, y);
			unsigned int depth = 4;
			for (int i = 0; i < SAMPLESPERPIXEL; ++i)
			{
				IntersectionData data = IntersectionData();
				CastRay(r, data);
				color += PathTrace(data, r, depth);
			}
			color = color * (1.0f / SAMPLESPERPIXEL);

#ifdef ANTIALIAS
			//// Create additional rays if anti-aliasing is desired
			//float N = 9.0f; // Set as a square
			//float halfSquare = 0.5f * sqrtf(N);
			//float hOffset = 1.0f / mpBitmap->Width();
			//float vOffset = 1.0f / mpBitmap->Height();

			//for (float i = -halfSquare; i < halfSquare; ++i)
			//{
			//	for (float j = -halfSquare; j < halfSquare; ++j)
			//	{
			//		Ray r = mScene.mCamera.GetRay(x+hOffset*i, y+vOffset*j);
			//		IntersectionData data = IntersectionData();
			//		CastRay(r, data);
			//		color += RayTrace(data, r, 1.0, 10);
			//	}
			//}

			//color = color / N;
#else
			//// Calculate the ray for the pixel
			//Ray r = mScene.mCamera.GetRay(x, y);
			//IntersectionData data = IntersectionData();
			//if (CastRay(r, data))
			//{
			//	color = RayTrace(data, r, 1.0, 10);
			//}
#endif

			// Gamma correction
			//float igamma = 1.0f/2.2f;
			//color = glm::vec3(pow(color.r,igamma),
			//                  pow(color.g,igamma),
			//                  pow(color.b,igamma));
			color = glm::clamp(255.0f * color, 0.0f, 255.0f);
			mpBitmap->SetPixel(i, j, color.r, color.g, color.b);
		}
	}
}

// Determines if cast ray hits an object and stores intersection data
bool Client::CastRay(const Ray& r, IntersectionData& closestObjectData)
{
	// Set the smallest depth value
	float smallestDepth = std::numeric_limits<float>::max();
	closestObjectData.object = nullptr;

	// Create intersection data
	IntersectionData data = IntersectionData();
	Object* pCurrentObject = nullptr;

	// Loop through all the objects in the scene
	for (Object* pObject : mScene.mObjects)
	{
		// Check for ray object intersection
		if (pObject->IsIntersection(r, data))
		{
			// Update closest object if a closer depth is found
			if (data.depth < smallestDepth)
			{
				smallestDepth = data.depth;
				closestObjectData = data;
				closestObjectData.object = pObject;
			}
		}
	}

	// Get color information from closest object
	if (closestObjectData.object == nullptr)
	{
		return false;
	}

	return true;
}

// Used for casting rays when object intersection data isnt wanted
bool Client::CastRay(const Ray& r)
{
	IntersectionData temp;
	return CastRay(r, temp);
}

// Apply lighting given intersection data
glm::vec3 Client::RayTrace(const IntersectionData& d, const Ray& r, float ni, unsigned int depth)
{
	// Check for no interseciton and end of recursion depth
	if (depth == 0 || d.object == nullptr)
	{
		return glm::vec3(0.0f, 0.0f, 0.0f);
	}

	// Get object information
	glm::vec3 N = (ni == 1.0f) ?  glm::normalize(d.normal) : -glm::normalize(d.normal);
	glm::vec3 eps = EPSILON * N;
	float ks = d.object->mSpecularCoefficient;

	// Get index of refraction for transmitted material
	float nt = (ni == 1.0f) ? sqrtf(d.object->mPermeability * d.object->mPermittivity) : 1.0f;

	// Get light attenuation coefficients
	float t = d.depth;
	glm::vec3 A = (ni == 1.0) ? mScene.mAir.mAttenuationFactor : d.object->mAttenuiationCoefficient;

	// Get magnetic permeability
	float ui = (ni == 1.0f) ? 1.0f : d.object->mPermeability;
	float ut = (ui == 1.0f) ? d.object->mPermeability : 1.0f;

	// Compute reflection coefficient R
	float R = CalculateReflectionCoefficient(glm::normalize(d.view), N, ni, nt, ui, ut);
	float T = 1.0f - R;
	R = R * ks;

	// Go through all lights and apply phong local lighting for now
	glm::vec3 color = LocalLighting(d, R, N);

	// Generate reflected ray if reflection coefficient not zero
	if (R > 0.0f)
	{
		// Generate reflected ray
		glm::vec3 V = glm::normalize(d.view);
		float NdotV = glm::dot(N, V);

		Ray reflected{ d.point + eps, glm::normalize((2.0f * (NdotV)* N) - V) };
		IntersectionData reflectData = IntersectionData();
		if (CastRay(reflected, reflectData))
		{
			color += R * RayTrace(reflectData, reflected, ni, depth - 1);
		}
	}

	// Get transmission coefficient from reflection coefficient
	//float T = 1.0f - R;
	T = ks * T;
	
	// Generate transmitted ray
	if (T > 0.0f)
	{
		// Compute transmission direction
		const glm::vec3& i = glm::normalize(-r.direction);
		float idotn = glm::dot(i, N);

		float ni_nt2 = (ni / nt) * (ni / nt);
		float radicand = 1.0f - (ni_nt2 * (1.0f - idotn * idotn));
		// Check for total interal reflection
		if (radicand >= 0.0f)
		{
			float cosT = sqrtf(radicand);

			// Caveat if (idotn) >= 0, tidotn = -costhetat
			//		  if (idotn) < 0, tidotn = costhetat
			float tdotn = (idotn >= 0.0f) ? -cosT : cosT;

			Ray transmitted{ d.point - eps, glm::normalize((tdotn + (ni / nt) * idotn) * N - (ni / nt) * i) };
			IntersectionData transmitData = IntersectionData();
			if (CastRay(transmitted, transmitData))
			{
				color += T * RayTrace(transmitData, transmitted, nt, depth - 1);
			}
		}
	}

	// Apply light attenuation
	glm::vec3 attenuationFactor = glm::vec3(powf(A.r, t), powf(A.g, t), powf(A.b, t));
	color = color * attenuationFactor;

	return color;
}

glm::vec3 Client::PathTrace(const IntersectionData& d, const Ray& r, unsigned int depth)
{
	// Check for no interseciton and end of recursion depth
	if (depth == 0 || d.object == nullptr)
	{
		return glm::vec3(0.0f, 0.0f, 0.0f);
	}

	// Check if object is a light, return radiance if so
	if (d.object->mIsLight)
	{
		return d.object->mDiffuseCoefficient * 100.0f;
	}

	// Clear distributions
	for (Distribution* dist : mCombinedDistribution.distributions)
	{
		delete dist;
	}
	mCombinedDistribution.distributions.clear();

	// Set up uniform distribution
	glm::vec3 N = glm::normalize(d.normal);
	Uniform* u = new Uniform(N);
	mCombinedDistribution.distributions.push_back(u);

	// Set up distributions based on number of lights
	for (Light* pLight : mScene.mLights)
	{
		glm::vec3 L = glm::normalize(pLight->mCenter - d.point);
		float distanceToLight = glm::distance(pLight->mCenter, d.point);
		float m = sqrtf(2.0f) * (pLight->mRadius / distanceToLight);
		Beckmann* b = new Beckmann(L, m);
		mCombinedDistribution.distributions.push_back(b);
	}
	
	glm::vec3 omega = mCombinedDistribution.generate();
	float probability = mCombinedDistribution.probability(omega);
	//Log(probability);
	// Create new ray
	glm::vec3 eps = EPSILON * N;
	Ray path{ d.point + eps, omega };
	IntersectionData data = IntersectionData();
	CastRay(path, data);

	// Return color and cast ray bounces
	return (d.object->mDiffuseCoefficient)*(1.0f/(2.0f*PI))*glm::dot(N, omega)*(1.0f/probability)*PathTrace(data, path, depth-1);
}

float Client::CalculateReflectionCoefficient(const glm::vec3& i, const glm::vec3& N, float ni, float nt, float ui, float ut)
{
	float R = 1.0f;
	float idotn = glm::dot(i, N);

	float sinI = sqrtf(1.0f - idotn * idotn);

	// If total interal reflection the transmission isn't computed, reflection coefficient is unity
	if (sinI >= (nt / ni))
	{
		return R;
	}

	// Compute reflection coefficient
	float cosI = idotn;
	float ni_nt2 = (ni / nt) * (ni / nt);
	float cosT = sqrtf(1.0f - (ni_nt2 * (1.0f - cosI * cosI)));

	// Calculate light polarized in perpendicular direction
	float Eperp = ((ni / nt) * cosI - (ui / ut) * cosT)
				/ ((ni / nt) * cosI + (ui / ut) * cosT);

	// Calculate light polarized in parallel direction
	float Eparl = ((ui / ut) * cosI - (ni / nt) * cosT)
				/ ((ui / ut) * cosI + (ni / nt) * cosT);

	// Reflection coefficient from Fresnel Equations
	R = 0.5f * (Eperp * Eperp + Eparl * Eparl);

	if (R > 1.0f)
	{
		R = 1.0f;
	}

	return R;
}

glm::vec3 Client::LocalLighting(const IntersectionData& d, float R, const glm::vec3& N)
{
	// Get object lighting information
	glm::vec3 V = glm::normalize(d.view);
	glm::vec3 kd = d.object->mDiffuseCoefficient;
	float ns = d.object->mSpecularExponent;

	glm::vec3 color(0.0f, 0.0f, 0.0f);

	for (Light* pLight : mScene.mLights)
	{
		// Get L vector
		glm::vec3 L = glm::normalize(pLight->mCenter - d.point);

		// Check if point is in shadow
		Ray shadowFeeler{ d.point + EPSILON * N, L };
		IntersectionData shadowData = IntersectionData();
		CastRay(shadowFeeler, shadowData);

		// Skip lighting if point is in shadow
		if (shadowData.depth > 0.0f && shadowData.depth <= glm::distance(pLight->mCenter, d.point))
		{
			continue;
		}

		float NdotL = glm::dot(N, L);
		if (NdotL < 0.0f)
		{
			continue;
		}

		// Calculate reflection vector
		glm::vec3 Ref = glm::normalize((2.0f * NdotL * N) - L);
		float VdotR = glm::clamp(glm::dot(V, Ref), 0.0f, 1.0f);

		// Get light color
		glm::vec3 lightColor = pLight->mColor;

		glm::vec3 diffuseColor = kd * NdotL * lightColor;
		glm::vec3 specularColor = R * powf(VdotR, ns) * lightColor;

		color += diffuseColor + specularColor;
	}

	// Apply ambient light
	color = color + kd * mScene.mAmbientLight;
	return color;
}

// Creates a thread to start rendering
void Client::RenderStart(void) 
{
	mKillFlag = false;
	mpRenderThread = new std::thread(&Client::Draw, std::ref(*this));
	//mpRenderThread->join();
}

// Sets kill flag to true, deletes thread
void Client::RenderStop(void) 
{
	if (mpRenderThread != nullptr)
	{
		mKillFlag = true;
		mpRenderThread->join();
		delete mpRenderThread;
		mpRenderThread = nullptr;
	}
}

// Draws bitmap to window
void Client::Blit(void) 
{
	glDrawPixels(mpBitmap->Width(), mpBitmap->Height(), GL_RGB, GL_UNSIGNED_BYTE, mpBitmap->Bytes());
	SDL_GL_SwapWindow(mpWindow);
}

// Handles adjusting window size
void Client::Resize(void) 
{
	int W, H;
	SDL_GetWindowSize(mpWindow, &W, &H);
	float aspect = glm::length(mScene.mCamera.mRight) / glm::length(mScene.mCamera.mUp);
	W = (int)round(H * aspect);
	SDL_SetWindowSize(mpWindow, W, H);
	delete mpBitmap;
	mpBitmap = new Bitmap(W, H);
}

// Saves bitmap data to file
void Client::Save(void) {
	// Write bitmap file
	int dataSize = mpBitmap->Height() * mpBitmap->Stride();
	std::fstream out("RayCastTest.bmp", std::ios_base::binary | std::ios_base::out);
	
	// Write bitmap header
	char header[54] = { 'B', 'M', 0 };
	*reinterpret_cast<unsigned*>(header + 2) = 54 + dataSize;
	*reinterpret_cast<unsigned*>(header + 10) = 54;
	*reinterpret_cast<unsigned*>(header + 14) = 40;
	*reinterpret_cast<int*>(header + 18) = mpBitmap->Width();
	*reinterpret_cast<int*>(header + 22) = mpBitmap->Height();
	*reinterpret_cast<unsigned short*>(header + 26) = 1;
	*reinterpret_cast<unsigned short*>(header + 28) = 24;
	*reinterpret_cast<unsigned*>(header + 34) = dataSize;
	out.write(header, 54);

	// Write bitmap data
	char* data = new char[dataSize];
	for (int j = 0; j < mpBitmap->Height(); ++j)
	{
		for (int i = 0; i < mpBitmap->Width(); ++i) 
		{
			int index = j * mpBitmap->Stride() + i * 3;
			// Write in BGR order
			data[index + 0] = mpBitmap->Bytes()[index + 2];
			data[index + 1] = mpBitmap->Bytes()[index + 1];
			data[index + 2] = mpBitmap->Bytes()[index + 0];
		}
	}
	out.write(data, dataSize);

	delete[] data;
}