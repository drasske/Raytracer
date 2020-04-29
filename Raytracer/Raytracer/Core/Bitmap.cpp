#include "Bitmap.h"

// Bitmap constructor
Bitmap::Bitmap(int W, int H) : mWidth(W), mHeight(H), mStride(4 * int(ceil(0.75f * W))), mpData(new unsigned char[H * mStride])
{
}

// Bitmap destructor, deletes data
Bitmap::~Bitmap()
{
	delete mpData;
}

char* Bitmap::Bytes()
{
	return reinterpret_cast<char*>(mpData);
}

int Bitmap::Width() const
{
	return mWidth;
}

int Bitmap::Height() const
{
	return mHeight;
}

int Bitmap::Stride() const
{
	return mStride;
}

void Bitmap::SetPixel(int i, int j, int r, int g, int b)
{
	mpData[mStride * j + 3 * i + 0] = r;
	mpData[mStride * j + 3 * i + 1] = g;
	mpData[mStride * j + 3 * i + 2] = b;
}
