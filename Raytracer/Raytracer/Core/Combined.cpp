#include "Distribution.h"

float Combined::probability(const glm::vec3& v)
{
	float prob = 0.0f;
	// Take the average of the probabilities of all distributions
	for (Distribution* pDistribution : distributions)
	{
		prob += pDistribution->probability(v);
	}
	return (1.0f / distributions.size()) * prob;
}

glm::vec3 Combined::generate()
{
	float u = mDistribution(mGenerator);
	// Choose a random sample to generate from the possible distributions
	for (unsigned int index = 0; index < distributions.size(); ++index)
	{
		if (u < ((index + 1.0f) / distributions.size()))
		{
			return distributions[index]->generate();
		}
	}
	// Default return from first distribution
	return distributions[0]->generate();
}