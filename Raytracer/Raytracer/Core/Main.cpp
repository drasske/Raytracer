// Main.cpp

#include <SDL2/SDL.h>
#include <GL/glew.h>
#include <GL/gl.h>

#include "Client.h"

int main(int argc, char* argv[]) 
{
	// Initialize and create a window
	SDL_Init(SDL_INIT_VIDEO);
	const char* title = "Press 'p' to write bitmap to file";
	int width = 500;
	int height = 500;
	SDL_Window* pWindow = SDL_CreateWindow(title, SDL_WINDOWPOS_UNDEFINED,
		SDL_WINDOWPOS_UNDEFINED, width, height,
		SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);
	SDL_GLContext context = SDL_GL_CreateContext(pWindow);

	// Create client
	Client* pClient = new Client(pWindow);

	bool done = false;

	while (!done) 
	{
		// Poll events
		SDL_Event event;
		while (SDL_PollEvent(&event)) 
		{
			switch (event.type) 
			{
				case SDL_QUIT:
					done = true;
					break;
				case SDL_KEYDOWN:
					if (event.key.keysym.sym == SDLK_ESCAPE)
						done = true;
					else if (event.key.keysym.sym == SDLK_p)
						pClient->Save();
					break;
				case SDL_DROPFILE:
					pClient->RenderStop();
					pClient->Load(event.drop.file);
					pClient->Resize();
					pClient->RenderStart();
					SDL_free(event.drop.file);
					break;
			}
		}
		pClient->Blit();
	}

	// Clean up
	delete pClient;
	SDL_GL_DeleteContext(context);
	SDL_DestroyWindow(pWindow);
	SDL_Quit();

	return 0;
}