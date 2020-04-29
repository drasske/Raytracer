// Bitmap.h
#pragma once

#include <glm/glm.hpp>

class Bitmap {
public:
	Bitmap(int W, int H);
	~Bitmap();

	char* Bytes();
	int Width() const;
	int Height() const;
	int Stride() const;
	void SetPixel(int i, int j, int r, int g, int b);
		
private:
	int mWidth;
	int mHeight;
	int mStride;
	unsigned char* mpData;
};
