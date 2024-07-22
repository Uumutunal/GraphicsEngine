#include <SDL.h>
#include<iostream>
#include<vector>

#include"methods.h"

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

struct vertex {
	vec3d p;
	int id;
	vec3d normal;
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


vec3d subVector(vec3d vec1, vec3d vec2) {
	vec3d subbed;
	subbed.x = vec1.x - vec2.x;
	subbed.y = vec1.y - vec2.y;
	subbed.z = vec1.z - vec2.z;

	return subbed;
}

bool fillTriangle2(vec3d p, vec3d a, vec3d b, vec3d c) {
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


void fillTriangle(vector<vec3d> vertices) {

	if (vertices.size() < 3) {
		return;
	}
	vertices.push_back(vertices[0]);
	vertices.push_back(vertices[1]);


	vec3d longestEdge[2];
	vec3d thirdVertex;
	float edgeLengthSquared = -1;

	for (size_t i = 0; i < vertices.size() - 2; i++)
	{
		vec3d edge = subVector(vertices[i + 1], vertices[i]);

		float temp = edge.x * edge.x + edge.y * edge.y;
		if (edgeLengthSquared < temp) {
			edgeLengthSquared = temp;

			longestEdge[0] = vertices[i];
			longestEdge[1] = vertices[i + 1];
			thirdVertex = vertices[i + 2];

		}

	}


	//Equation for the longest edge, Ax + By + C = 0
	//A = y2 - y1
	//B = x1 - x2
	//C = (x2 * y1) - (x1 * y2)
	float a = longestEdge[1].y - longestEdge[0].y;
	float b = longestEdge[0].x - longestEdge[1].x;
	float c = (longestEdge[1].x * longestEdge[0].y) - (longestEdge[0].x * longestEdge[1].y);

	//Equation for the parallel line that passes through the third point, Ax + By + D = 0
	float d = -a * thirdVertex.x - b * thirdVertex.y;


	//Equations for the other edges
	float a2 = longestEdge[1].y - thirdVertex.y;
	float b2 = thirdVertex.x - longestEdge[1].x;
	float c2 = (longestEdge[1].x * thirdVertex.y) - (thirdVertex.x * longestEdge[1].y);

	float a3 = thirdVertex.y - longestEdge[0].y;
	float b3 = longestEdge[0].x - thirdVertex.x;
	float c3 = (thirdVertex.x * longestEdge[0].y) - (longestEdge[0].x * thirdVertex.y);

	//Loop through the points on the longest edge


	for (float x = longestEdge[0].x; x < longestEdge[1].x; x++)
	{
		//Calculate the y for x
		float y = -(a * x + c) / b;

		//Calculate perpendicular line
		float aP = -b;
		float bP = a;

		float cP = -(aP * x + bP * y);


		//Find intersection point
		float determinant = a * bP - aP * b;

		float x2 = (b * cP - bP * d) / determinant;
		float y2 = (aP * d - a * cP) / determinant;

		int sign2 = (a2 * longestEdge[0].x + b2 * longestEdge[0].y + c2);

		int sign3 = (a3 * longestEdge[1].x + b3 * longestEdge[1].y + c3);

		int s2;
		int s3;

		if (sign2 < 0) {
			s2 = 1;
		}
		else {
			s2 = -1;
		}

		if (sign3 < 0) {
			s3 = 1;
		}
		else {
			s3 = -1;
		}

		for (float i = x; i < x2; i++)	
		{
			vec3d p;
			p.x = i;
			p.y = (aP * x + cP) / -bP;

			//Fix <> signs
			if ((a2 * p.x + b2 * p.y + c2)*s2 < 0 && (a3 * p.x + b3 * p.y + c3)*s3 < 0) {

				SDL_RenderDrawPoint(renderer, p.x + 320, p.y + 320);

			}
			else {
				break;
			}
		}

	}


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

		fillTriangle(pos);

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
	center.z = 50;

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



	int close = 0;



	SDL_SetRenderDrawColor(renderer, 255 * 0.227, 255 * 0.227, 255 * 0.227, 255);
	SDL_RenderClear(renderer);


	SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);


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

	//cubeV.addFace({ {0,1,2},{1,3,2},{4,5,6},
	//	{5,7,6},{4,0,6},{6,0,2},
	//	{4,5,0},{5,1,0},{1,5,3},
	//	{3,5,7},{2,3,6},{3,7,6} });


	cubeV.addFace({ {0,1,2} });

	//move(cubeV, { 100,100,25 });


	//rotate(cubeV, { 1, 0, 0 }, 3.14159f / 4);
	//

#pragma endregion



	// animation loop
	while (!close) {

		SDL_SetRenderDrawColor(renderer, 255 * 0.227, 255 * 0.227, 255 * 0.227, 255);
		SDL_RenderClear(renderer);

		SDL_Event event;


		//rotate(cubeV, { 1, 1, 1 }, 3.14159f / 250);

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
					move(cubeV, { 0,10,0 });
					//dest.y -= speed / 30;
					break;
				case SDL_SCANCODE_A:
				case SDL_SCANCODE_LEFT:
					move(cubeV, { 10,0,0 });
					break;
				case SDL_SCANCODE_S:
				case SDL_SCANCODE_DOWN:
					move(cubeV, { 0,-10,0 });
					break;
				case SDL_SCANCODE_D:
				case SDL_SCANCODE_RIGHT:
					move(cubeV, { -10,0,0 });
					break;
				case SDL_SCANCODE_G:

					break;
				default:
					break;
				}
			}

		}








		SDL_RenderPresent(renderer);

		SDL_Delay(1000 / 60);
		a++;
	}





	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();
	return 0;
}