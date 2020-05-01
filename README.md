# Raytracer
Ray tracing implementation done for CS 500 at Digipen.
Has the options to run either a naive ray tracing algorithm or a stochastic path tracing algorithm.

### Naive Ray Tracing
For the naive ray tracing approach, the output of each pixel is calculated starting with one ray. This ray is shot into the scene and creates a transmissive and reflective ray if an object hit. This process continues recursively until a light is hit or a maximum depth is reached.

#### Antialising
Antialising can be enabled to smooth out jagged lines. It is implemented as super sampling each pixel. Instead of only using one ray per pixel, each pixel is subdivided into nine smaller pixels that are then averaged to get the final result.

#### Transmission And Reflection
Transmission and reflection coefficients are calculated using the Fresnel equations. They account for transparent or reflective objects. 

### Stochastic Path Tracing
For the stochastic path tracing approach, each pixel is sampled multiple times. A ray is shot in a random direction and if an object is hit, continues to bounce until it hits a light or nothing. Each time an object is hit, a new ray is generated at that point. This causes the images to have a large amount of noise in them unless the number of samples per pixel is large ~4000.

#### Importance Sampling
Importance sampling is used to influence the randomness used by the stochastic path tracing. Instead of creating a ray in a random direction using a uniform distribution, a combined distribution, in this case combined Beckmann distributions for each light, is used. This greatly reduces the amount of noise with the same number of samples per pixel. 
![Image](https://github.com/drasske/Raytracer/blob/master/Images/ImportanceSampling/Ptest_Importance_1.bmp)\
![Image](https://github.com/drasske/Raytracer/blob/master/Images/ImportanceSampling/Ptest_Importance_10.bmp)
![Image](https://github.com/drasske/Raytracer/blob/master/Images/ImportanceSampling/Ptest_Importance_100.bmp)
![Image](https://github.com/drasske/Raytracer/blob/master/Images/ImportanceSampling/Ptest_Importance_1000.bmp)

#### Potential Updates
The path tracer does not account for transmission and reflection like the naive ray tracer does.
