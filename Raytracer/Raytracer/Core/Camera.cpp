#include "Camera.h"

Ray Camera::GetRay(float x, float y) const
{
	// x, y are NDC coordinates of the pixel
	// Need to find the A = aU + bV
	// Direction isnt just x, y should be Q - Pe
	// Where Q is Center + aRIGHT + bUP
	// And Pe = CENTER + EYE
	glm::vec3 position = mCenter + mEye;
	glm::vec3 Q = mCenter + x * mRight + y * mUp;
	glm::vec3 direction = glm::normalize(Q - position);

	return Ray{position, direction};
}
