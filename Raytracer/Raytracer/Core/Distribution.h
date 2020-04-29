#pragma once
#include <random>

#include <glm/glm.hpp>


class Distribution
{
public:
	virtual ~Distribution() {}
	virtual float probability(const glm::vec3& v) = 0;
	virtual glm::vec3 generate() = 0;

protected:
	const float PI = 4.0f * atan(1.0f);
	static std::random_device mRandomDevice;
	static std::mt19937 mGenerator;
	static std::uniform_real_distribution<float> mDistribution;
};

class Uniform : public Distribution
{
public:
	Uniform(const glm::vec3 N) : direction(N) {}

	float probability(const glm::vec3& v);
	glm::vec3 generate();

	glm::vec3 direction;
};

class Beckmann : public Distribution
{
public:
	Beckmann(const glm::vec3 L, float m) : direction(L), roughness(m) {}
	float probability(const glm::vec3& v) override;
	glm::vec3 generate() override;

public:
	glm::vec3 direction;
	float roughness;
};

class Combined : public Distribution
{
public:
	std::vector<Distribution*> distributions;
	float probability(const glm::vec3& v) override;
	glm::vec3 generate() override;
};