#include <SDL.h>
#include<iostream>
#include<vector>
#include <list>

using namespace std;

const int WIDTH = 640;
const int HEIGHT = 640;
SDL_Window* window = NULL;
SDL_Renderer* renderer = NULL;



float fNear = 300;
float fFar = 0;

struct vec2d
{
	float x, y;
};

struct vec3d
{
	float x, y, z;
};



struct triangle
{
	vec2d p[3];
};

struct triangle3d
{
	vec3d p[3];
};



struct mesh
{
	triangle tris[8];
};

struct mesh3d
{

	triangle3d tris[8];

};

//***********

struct vertex {
	vec3d p;
	int id;
};


struct face {
	vector<int> points;
};

struct meshv
{


private:
	int index = 0;

public:
	vector<vertex> vertices;

	vector<face> faces;

	vec3d origin;

	vec3d position;

	void addVertex(vector<vertex> verticesA) {
		for (size_t i = 0; i < verticesA.size(); i++)
		{
			vertex v;
			v.p = verticesA[i].p;
			v.id = index;
			vertices.push_back(v);
			index++;
		}
	}

	void addFace(vector<vector<int>> facesA) {
		for (size_t i = 0; i < facesA.size(); i++)
		{
			face f;
			f.points = facesA[i];
			faces.push_back(f);
		}
	}

};


void createMesh() {

	float cube[3][2] = { {0.0f, 0.0f}, {0.0f, 1.0f}, {1.0f, 1.0f} };

	triangle tri1 = { vec2d {0.0f, 0.0f}, vec2d {0.0f, 1.0f}, vec2d {1.0f, 1.0f} };

	triangle mesh[] = { tri1 };


}

vector<int> a = { 5 };


void drawTris(float tri[][2]) {


	SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
	int j;
	for (size_t i = 0; i < 3; i++)
	{
		if (i != 2) {
			j = i + 1;
		}
		else {
			j = 0;
		}


		SDL_Point points[2] = {
		{tri[i][0], tri[i][1]},
		{tri[j][0], tri[j][1]},
		};

		SDL_RenderDrawLines(renderer, points, 2);
	}

}


void drawTris(mesh m) {


	SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);



	for (auto tri : m.tris) {

		int j;
		for (size_t i = 0; i < 3; i++)
		{
			if (i != 2) {
				j = i + 1;
			}
			else {
				j = 0;
			}


			SDL_Point points[2] = {
			{tri.p[i].x + WIDTH / 2, tri.p[i].y + HEIGHT / 2},
			{tri.p[j].x + WIDTH / 2, tri.p[j].y + HEIGHT / 2},
			};

			SDL_RenderDrawLines(renderer, points, 2);
		}
	}

}

void drawTris(mesh3d m) {


	SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);



	for (auto tri : m.tris) {

		int j;
		for (size_t i = 0; i < 3; i++)
		{
			if (i != 2) {
				j = i + 1;
			}
			else {
				j = 0;
			}


			SDL_Point points[2] = {
			{tri.p[i].x + WIDTH / 2, tri.p[i].y + HEIGHT / 2},
			{tri.p[j].x + WIDTH / 2, tri.p[j].y + HEIGHT / 2},
			};

			SDL_RenderDrawLines(renderer, points, 2);
		}
	}

}

mesh3d& project(mesh3d& m) {

	for (auto& tri : m.tris) {


		for (size_t i = 0; i < 3; i++)
		{

			float z = tri.p[i].z;

			float zRate = z / (fFar - fNear);

			tri.p[i].x = tri.p[i].x * zRate;
			tri.p[i].y = tri.p[i].y * zRate;

		}
	}

	return m;
}

meshv& project(meshv& m) {

	for (auto& vert : m.vertices) {

		float z = m.position.z + vert.p.z;

		float zRate = z / (fFar - fNear);

		//zRate = 1;

		vert.p.x = (m.position.x + vert.p.x) * zRate;
		vert.p.y = (m.position.y + vert.p.y) * zRate;
	}

	return m;
}

void drawTris(meshv m) {


	project(m);

	SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);


	vector<vec3d> pos;

	for (auto face : m.faces) {

		int j;
		pos.clear();

		for (size_t i = 0; i < 3; i++)
		{

			int index = face.points[i];



			for (auto v : m.vertices) {

				if (v.id == index) {
					pos.push_back(v.p);
				}

			}

		}

		for (size_t i = 0; i < 3; i++) {
			if (i != 2) {
				j = i + 1;
			}
			else {
				j = 0;
			}


			SDL_Point points[2] = {
			{pos[i].x + WIDTH / 2,pos[i].y + HEIGHT / 2},
			{pos[j].x + WIDTH / 2, pos[j].y + HEIGHT / 2},
			};

			SDL_RenderDrawLines(renderer, points, 2);
		}


	}


}





bool fillTriangle(vec3d p, vec3d a, vec3d b, vec3d c) {
	// Compute vectors
	vec2d v0 = { c.x - a.x, c.y - a.y };
	vec2d v1 = { b.x - a.x, b.y - a.y };
	vec2d v2 = { p.x - a.x, p.y - a.y };

	// Compute dot products
	float dot00 = v0.x * v0.x + v0.y * v0.y;
	float dot01 = v0.x * v1.x + v0.y * v1.y;
	float dot02 = v0.x * v2.x + v0.y * v2.y;
	float dot11 = v1.x * v1.x + v1.y * v1.y;
	float dot12 = v1.x * v2.x + v1.y * v2.y;

	// Compute barycentric coordinates
	float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
	float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	// Check if point is in triangle
	return (u >= 0) && (v >= 0) && (u + v <= 1);
}

// Cross product
vec3d cross(const vec3d& a, const vec3d& b) {
	return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
}

// Dot product
float dot(const vec3d& a, const vec3d& b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

vec3d multiplayVector(vec3d vec, float f) {
	vec3d multiplied;
	multiplied.x = vec.x * f;
	multiplied.y = vec.y * f;
	multiplied.z = vec.z * f;

	return multiplied;
}

vec3d addVector(vec3d vec1, vec3d vec2) {
	vec3d added;
	added.x = vec1.x + vec2.x;
	added.y = vec1.y + vec2.y;
	added.z = vec1.z + vec2.z;

	return added;
}

// Rotate point around an axis
vec3d rotateAroundPoint(const vec3d& point, const vec3d mCenter, const vec3d& axis, float angle) {
	// Translate point to origin
	vec3d center;
	center.x = mCenter.x;
	center.y = mCenter.y;
	center.z = mCenter.z - 450;

	vec3d p = { point.x - center.x, point.y - center.y, point.z - center.z };

	// Normalize the axis
	float length = std::sqrt(axis.x * axis.x + axis.y * axis.y + axis.z * axis.z);
	vec3d u = { axis.x / length, axis.y / length, axis.z / length };

	// Rodrigues' rotation formula
	vec3d p_rot =
		addVector(addVector(multiplayVector(p, std::cos(angle)),
			multiplayVector(cross(u, p), std::sin(angle))),
			multiplayVector(multiplayVector(u, dot(u, p)), (1 - std::cos(angle))));

	// Translate point back
	return { p_rot.x + center.x, p_rot.y + center.y, p_rot.z + center.z };
}

void move(meshv& m, vec3d v) {
	m.position.x += v.x;
	m.position.y += v.y;
	m.position.z += v.z;
}

void rotate(meshv& m, const vec3d& axis, float angle) {


	for (auto& vert : m.vertices) {

		vec3d rotated = rotateAroundPoint(vert.p, m.position, axis, angle);

		vert.p.x = rotated.x;
		vert.p.y = rotated.y;
		vert.p.z = rotated.z;
	}

}



int main(int argc, char* argv[])
{




	SDL_Init(SDL_INIT_VIDEO);
	window = SDL_CreateWindow("Graphics Engine", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WIDTH, HEIGHT, SDL_WINDOW_SHOWN);
	renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);



	createMesh();

	int close = 0;




	SDL_SetRenderDrawColor(renderer, 255 * 0.227, 255 * 0.227, 255 * 0.227, 255);
	SDL_RenderClear(renderer);


	//SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
	//SDL_RenderClear(renderer);
	SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
	//for (int i = 0; i < 640; ++i)
	//	SDL_RenderDrawPoint(renderer, i, i);



	SDL_Point points[2] = {
	{0, 0},
	{640, 640},
	};

	//SDL_RenderDrawLines(renderer, points, 2);

	SDL_Point points2[2] = {
	{0, 640},
	{640, 0},
	};

	//SDL_RenderDrawLines(renderer, points2, 2);


	//float tri[3][2] = { {50.0f, 50.0f}, {100.0f, 50.0f}, {50.0f, 100.0f} };
	triangle tri = { vec2d {50.0f, 50.0f}, {100.0f, 50.0f}, {50.0f, 100.0f} };

	mesh cube;
	cube.tris[0] = { vec2d {50.0f, 50.0f}, {100.0f, 50.0f}, {50.0f, 100.0f} };
	cube.tris[1] = { vec2d {100.0f, 50.0f}, {100.0f, 100.0f}, {50.0f, 100.0f} };

	cube.tris[2] = { vec2d {50.0f, 50.0f}, {100.0f, 50.0f}, {50.0f, 100.0f} };
	cube.tris[3] = { vec2d {100.0f, 50.0f}, {100.0f, 100.0f}, {50.0f, 100.0f} };





	SDL_RenderPresent(renderer);



	int a = 0;
	float time = 0;



#pragma region cube

	SDL_SetRenderDrawColor(renderer, 255 * 0.227, 255 * 0.227, 255 * 0.227, 255);
	SDL_RenderClear(renderer);

	//
			//*********
	float z1 = -50;
	float z2 = 50;

	meshv cubeV;
	cubeV.position = { 0,0,500 };
	cubeV.addVertex({
		{-50.0f, 50.0f, z1}, {50.0f, 50.0f, z1}, {-50.0f, -50.0f, z1},{50.0f, -50.0f, z1},
		{-50.0f, 50.0f, z2}, {50.0f, 50.0f, z2}, {-50.0f, -50.0f, z2},{50.0f, -50.0f, z2},
		}
	);
	cubeV.addFace({ {0,1,2},{1,3,2},{4,5,6},
		{5,7,6},{4,0,6},{6,0,2},
		{4,5,0},{5,1,0},{1,5,3},
		{3,5,7},{2,3,6},{3,7,6} });

	//move(cubeV, { 100,100,25 });


	//rotate(cubeV, { 1, 0, 0 }, 3.14159f / 4);
	//

#pragma endregion



	// animation loop
	while (!close) {

		SDL_SetRenderDrawColor(renderer, 255 * 0.227, 255 * 0.227, 255 * 0.227, 255);
		SDL_RenderClear(renderer);

		SDL_Event event;


		rotate(cubeV, { 1, 1, 1 }, 3.14159f / 250);

		drawTris(cubeV);


		while (SDL_PollEvent(&event)) {
			switch (event.type) {

			case SDL_QUIT:
				// handling of close button
				close = 1;
				break;

			case SDL_KEYDOWN:
				// keyboard API for key pressed
				switch (event.key.keysym.scancode) {
				case SDL_SCANCODE_W:
				case SDL_SCANCODE_UP:
					move(cubeV, { 0,-10,0 });
					//dest.y -= speed / 30;
					break;
				case SDL_SCANCODE_A:
				case SDL_SCANCODE_LEFT:
					move(cubeV, { -10,0,0 });
					break;
				case SDL_SCANCODE_S:
				case SDL_SCANCODE_DOWN:
					move(cubeV, { 0,10,0 });
					break;
				case SDL_SCANCODE_D:
				case SDL_SCANCODE_RIGHT:
					move(cubeV, { 10,0,0 });
					break;
				case SDL_SCANCODE_G:

					break;
				default:
					break;
				}
			}

		}






		mesh3d cube3d;
		cube3d.tris[0] = { vec3d {50.0f, 50.0f, z1}, {100.0f, 50.0f, z1}, {50.0f, 100.0f, z1} };
		cube3d.tris[1] = { vec3d {100.0f, 50.0f, z1}, {100.0f, 100.0f, z1}, {50.0f, 100.0f, z1} };

		cube3d.tris[2] = { vec3d {50.0f, 50.0f, z2}, {100.0f, 50.0f, z2}, {50.0f, 100.0f, z2} };
		cube3d.tris[3] = { vec3d {100.0f, 50.0f, z2}, {100.0f, 100.0f, z2}, {50.0f, 100.0f, z2} };

		//project(cube3d);

		//drawTris(cube);
		//drawTris(cube3d);






		SDL_RenderPresent(renderer);

		SDL_Delay(1000 / 60);
		a++;
	}





	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();
	return 0;
}