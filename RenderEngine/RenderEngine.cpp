#include <SDL.h>
#include<iostream>
#include<vector>
#include<math.h>
#include <chrono>
#include <map>
#include <algorithm>
#include <fstream>
#include <strstream>

using namespace std;


#include"methods.h"
#include"include/structs.h"
#include"include/vectorCalc.h"
#include"mesh.h"
#include"importer.h"



const int WIDTH = 1000;
const int HEIGHT = 1000;
SDL_Window* window = NULL;
SDL_Renderer* renderer = NULL;

vector<vector<pixel>> zBuffer(WIDTH, vector<pixel>(HEIGHT, { -std::numeric_limits<float>::infinity() }));

void reinitializeZBuffer(std::vector<std::vector<pixel>>& zBuffer, int width, int height) {
	for (int x = 0; x < width; ++x) {
		for (int y = 0; y < height; ++y) {
			zBuffer[x][y].depth = -std::numeric_limits<float>::infinity();
		}
	}
}


float fNear = 300;
float fFar = 0;



void drawLine(SDL_Renderer* renderer, const SDL_Point* points, float depth) {

	int startX = points[0].x;
	int endX = points[0].x + 1;
	int startY = points[0].y;
	int endY = points[1].y;

	if (startY > endY) {
		std::swap(startY, endY);
	}

	for (int j = startX - 1; j <= endX; j++) {
		for (int i = startY; i <= endY; i++) {
			// Boundary checks
			if (i >= 0 && i < HEIGHT && j >= 0 && j < WIDTH) {
				if (zBuffer[i][j].depth < depth) {
					SDL_RenderDrawPoint(renderer, j, i);
					zBuffer[i][j].depth = depth;
				}
			}
		}
	}

}













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




void project(meshv& m) {

	//for (auto& vert : m.vertices) {

	//	float z = m.position.z + vert.p.z;

	//	float zRate = z / (fFar - fNear);

	//	//zRate = 1;

	//	vert.p.x = (m.position.x + vert.p.x) * zRate;
	//	vert.p.y = (m.position.y + vert.p.y) * zRate;
	//}
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		float z = m.position.z + m.vertices[i].p.z;

		float zRate = z / (fFar - fNear);


		m.verticesProjected[i].p.x = (m.position.x + m.vertices[i].p.x) * zRate;
		m.verticesProjected[i].p.y = (m.position.y + m.vertices[i].p.y) * zRate;
	}

	//TODO: check this
	//for (auto& face : m.faces) {
	//for (int j = 0; j < m.faces.size(); j++) {
	//	for (size_t i = 0; i < 3; i++)
	//	{
	//		float z = m.position.z + m.faces[j].pos[0].z;

	//		float zRate = z / (fFar - fNear);

	//		//zRate = 1;

	//		m.faces[j].pos[i].x = (m.position.x + m.faces[j].pos[i].x) * zRate;
	//		m.faces[j].pos[i].y = (m.position.y + m.faces[j].pos[i].y) * zRate;
	//	}
	//}

	for (auto& f : m.faces) {

		vec3d v1 = subVector(m.vertices[f.points[1]].p, m.vertices[f.points[0]].p);
		vec3d v2 = subVector(m.vertices[f.points[2]].p, m.vertices[f.points[0]].p);
		vec3d normal = cross(v1, v2);
		normalize(normal);
		f.normal = normal;
	}
	int t = 0;

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

	vec3d top;
	vec3d bot;


	float max = -99999999;
	float min = 99999999;
	for (int i = 0; i < vertices.size(); i++)
	{
		if (vertices[i].x > max) {
			max = vertices[i].x;
			top = vertices[i];
		}
		if (vertices[i].x < min) {
			min = vertices[i].x;
			bot = vertices[i];
		}

	}

	float depth = (vertices[0].z + vertices[1].z + vertices[2].z) / 3.0;

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
	float b = (longestEdge[0].x - longestEdge[1].x);
	float c = (longestEdge[1].x * longestEdge[0].y) - (longestEdge[0].x * longestEdge[1].y);

	//Equation for the parallel line that passes through the third point, Ax + By + D = 0
	float d = -a * thirdVertex.x - b * thirdVertex.y;


	//Equations for the other edges
	float a2 = longestEdge[1].y - thirdVertex.y;
	float b2 = (thirdVertex.x - longestEdge[1].x);
	float c2 = (longestEdge[1].x * thirdVertex.y) - (thirdVertex.x * longestEdge[1].y);

	float a3 = thirdVertex.y - longestEdge[0].y;
	float b3 = (longestEdge[0].x - thirdVertex.x);
	float c3 = (thirdVertex.x * longestEdge[0].y) - (longestEdge[0].x * thirdVertex.y);

	//Loop through the points on the longest edge



	for (float x = min; x < max; x++)
	{
		float p1 = (-a * x - c) / b;
		float p2 = (-a2 * x - c2) / b2;
		float p3 = (-a3 * x - c3) / b3;


		if (b == 0) {
			SDL_Point points[2] = {
			{x + WIDTH / 2,p2 + HEIGHT / 2},
			{x + WIDTH / 2, p3 + HEIGHT / 2},
			};
			//SDL_RenderDrawLines(renderer, points, 2);
			drawLine(renderer, points, depth);
			continue;
		}
		else if (b2 == 0) {
			SDL_Point points[2] = {
			{x + WIDTH / 2,p1 + HEIGHT / 2},
			{x + WIDTH / 2, p3 + HEIGHT / 2},
			};
			//SDL_RenderDrawLines(renderer, points, 2);
			drawLine(renderer, points, depth);
			continue;
		}
		else if (b3 == 0) {
			SDL_Point points[2] = {
			{x + WIDTH / 2,p2 + HEIGHT / 2},
			{x + WIDTH / 2, p1 + HEIGHT / 2},
			};
			//SDL_RenderDrawLines(renderer, points, 2);
			drawLine(renderer, points, depth);
			continue;
		}


		float l1 = (top.x - x) * (top.x - x) + (top.y - p1) * (top.y - p1);
		float l2 = (top.x - x) * (top.x - x) + (top.y - p2) * (top.y - p2);
		float l3 = (top.x - x) * (top.x - x) + (top.y - p3) * (top.y - p3);

		float A;
		float B;
		float C;

		A = bot.y - top.y;
		B = (top.x - bot.x);
		C = (bot.x * top.y) - (top.x * bot.y);

		if (abs(A * x + B * p1 + C) <= 0.01) {
			float lH1 = abs(p1 - p2);
			float lH2 = abs(p1 - p3);


			if (lH1 < lH2) {
				SDL_Point points[2] = {
				{x + WIDTH / 2,p1 + HEIGHT / 2},
				{x + WIDTH / 2, p2 + HEIGHT / 2},
				};


				//SDL_RenderDrawLines(renderer, points, 2);
				drawLine(renderer, points, depth);
			}
			else {
				SDL_Point points[2] = {
				{x + WIDTH / 2,p1 + HEIGHT / 2},
				{x + WIDTH / 2, p3 + HEIGHT / 2},
				};
				//SDL_RenderDrawLines(renderer, points, 2);
				drawLine(renderer, points, depth);

			}

		}
		else if (abs(A * x + B * p2 + C) <= 0.01) {
			float lH1 = abs(p2 - p1);
			float lH2 = abs(p2 - p3);


			if (lH1 < lH2) {
				SDL_Point points[2] = {
				{x + WIDTH / 2,p2 + HEIGHT / 2},
				{x + WIDTH / 2, p1 + HEIGHT / 2},
				};
				//SDL_RenderDrawLines(renderer, points, 2);
				drawLine(renderer, points, depth);

			}
			else {
				SDL_Point points[2] = {
				{x + WIDTH / 2,p2 + HEIGHT / 2},
				{x + WIDTH / 2, p3 + HEIGHT / 2},
				};
				//SDL_RenderDrawLines(renderer, points, 2);
				drawLine(renderer, points, depth);

			}
		}
		else if (abs(A * x + B * p3 + C) <= 0.01) {
			float lH1 = abs(p3 - p2);
			float lH2 = abs(p3 - p1);


			if (lH1 < lH2) {
				SDL_Point points[2] = {
				{x + WIDTH / 2,p3 + HEIGHT / 2},
				{x + WIDTH / 2, p2 + HEIGHT / 2},
				};
				//SDL_RenderDrawLines(renderer, points, 2);
				drawLine(renderer, points, depth);

			}
			else {
				SDL_Point points[2] = {
				{x + WIDTH / 2,p3 + HEIGHT / 2},
				{x + WIDTH / 2, p1 + HEIGHT / 2},
				};
				//SDL_RenderDrawLines(renderer, points, 2);
				drawLine(renderer, points, depth);

			}
		}
	}

	return;

}


void drawTris(meshv& m) {


	project(m);


	for (auto& face : m.faces) {

		SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);

		int j;
		//pos.clear();
		vec3d abc = { 0,0,1 };
		float d = dot(face.normal, abc);
		//float d = face.normal.z;

		if (d > 0) {
			continue;
			//SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);
			int tt = 0;
		}
		if (face.normal.y > 0) {
			//SDL_SetRenderDrawColor(renderer, 0, 125, 125, 125);
			float col = -face.normal.y * 255;
			SDL_SetRenderDrawColor(renderer, col, col, col, 255);

		}



		//SDL_SetRenderDrawColor(renderer, 0, 255, 255, 255);
		//fillTriangle(pos);

		vector<vec3d> tri = { m.verticesProjected[face.points[0]].p ,m.verticesProjected[face.points[1]].p ,m.verticesProjected[face.points[2]].p };

		fillTriangle(tri);

		for (size_t i = 0; i < 3; i++) {
			if (i != 2) {
				j = i + 1;
			}
			else {
				j = 0;
			}

			SDL_Point points[2] = {
			{m.verticesProjected[face.points[i]].p.x + WIDTH / 2,m.verticesProjected[face.points[i]].p.y + HEIGHT / 2},
			{m.verticesProjected[face.points[j]].p.x + WIDTH / 2, m.verticesProjected[face.points[j]].p.y + HEIGHT / 2},
			};
			SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);

			//SDL_RenderDrawLines(renderer, points, 2);
		}
		//break;

	}

}








// Rotate point around an axis
vec3d rotateAroundPoint(const vec3d& point, const vec3d mCenter, const vec3d& axis, float angle, bool camera) {
	// Translate point to origin
	vec3d center;

	if (camera) {
		center = mCenter;
	}
	else {
		center.x = 0;
		center.y = 0;
		center.z = 0;
	}



	vec3d p = { point.x - center.x, point.y - center.y, point.z - center.z };

	// Normalize the axis
	// sqrt
	float length = (axis.x * axis.x + axis.y * axis.y + axis.z * axis.z);
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
	//meshv mP = m;
	//project(mP);
	//for (auto& vert : m.vertices) {

	//	vec3d rotated = rotateAroundPoint(vert.p, m.position, axis, angle);


	//	vert.p.x = rotated.x;
	//	vert.p.y = rotated.y;
	//	vert.p.z = rotated.z;

	//	//m.verticesProjected[vert.id].p.x = rotated.x;
	//	//m.verticesProjected[vert.id].p.y = rotated.y;
	//	//m.verticesProjected[vert.id].p.z = rotated.z;

	//	int t = 0;
	//}

	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		vec3d rotated = rotateAroundPoint(m.vertices[i].p, m.position, axis, angle, false);


		if (m.vertices[i].id == 1840) {
			int a = 0;
		}

		m.vertices[i].p.x = rotated.x;
		m.vertices[i].p.y = rotated.y;
		m.vertices[i].p.z = rotated.z;

		//m.verticesProjected[i].p.x = rotated.x;
		//m.verticesProjected[i].p.y = rotated.y;
		//m.verticesProjected[i].p.z = rotated.z;

		m.verticesProjected[i].p.x = m.vertices[i].p.x;
		m.verticesProjected[i].p.y = m.vertices[i].p.y;
		m.verticesProjected[i].p.z = m.vertices[i].p.z;

	}

	for (auto& f : m.faces) {

		//for (size_t i = 0; i < 3; i++)
		//{
		//	vec3d rotated = rotateAroundPoint(m.vertices[f.points[i]].p, m.position, axis, angle);
		//	m.vertices[f.points[i]].p.x = rotated.x;
		//	m.vertices[f.points[i]].p.y = rotated.y;
		//	m.vertices[f.points[i]].p.z = rotated.z;
		//}

		vec3d v1 = subVector(m.verticesProjected[f.points[1]].p, m.verticesProjected[f.points[0]].p);
		vec3d v2 = subVector(m.verticesProjected[f.points[2]].p, m.verticesProjected[f.points[0]].p);

		vec3d normal = cross(v1, v2);


		normalize(normal);
		f.normal = normal;


	}
	int ta = 0;

}

void rotateCamera(meshv& m, camera cam, const vec3d& axis, float angle) {
	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		vec3d focus = { 0,0,500 };

		vec3d rotated = rotateAroundPoint(m.vertices[i].p, cam.position, axis, angle, true);


		if (m.vertices[i].id == 1840) {
			int a = 0;
		}

		m.vertices[i].p.x = rotated.x;
		m.vertices[i].p.y = rotated.y;
		m.vertices[i].p.z = rotated.z;


		m.verticesProjected[i].p.x = m.vertices[i].p.x;
		m.verticesProjected[i].p.y = m.vertices[i].p.y;
		m.verticesProjected[i].p.z = m.vertices[i].p.z;

	}

	for (auto& f : m.faces) {

		vec3d v1 = subVector(m.verticesProjected[f.points[1]].p, m.verticesProjected[f.points[0]].p);
		vec3d v2 = subVector(m.verticesProjected[f.points[2]].p, m.verticesProjected[f.points[0]].p);

		vec3d normal = cross(v1, v2);


		normalize(normal);
		f.normal = normal;


	}
}




int main(int argc, char* argv[])
{

	SDL_Init(SDL_INIT_VIDEO);
	window = SDL_CreateWindow("Graphics Engine", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WIDTH, HEIGHT, SDL_WINDOW_SHOWN);
	renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);

	camera mainCamera;

	int close = 0;

	meshv importM;
	importMesh("../models/teapot.obj", importM);
	importM.position = { 0,-50,800 };

	SDL_SetRenderDrawColor(renderer, 255 * 0.227, 255 * 0.227, 255 * 0.227, 255);
	SDL_RenderClear(renderer);


	SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);


	SDL_RenderPresent(renderer);



	int a = 0;
	float time = 0;

	print();
	print2();

#pragma region cube

	SDL_SetRenderDrawColor(renderer, 255 * 0.227, 255 * 0.227, 255 * 0.227, 255);
	SDL_RenderClear(renderer);

	//
			//*********
	float z1 = 50;
	float z2 = -50;

	meshv cubeV;
	cubeV.position = { 0,0,500 };
	cubeV.addVertex({
		{-50.0f, 50.0f, z1}, {50.0f, 50.0f, z1}, {-50.0f, -50.0f, z1},{50.0f, -50.0f, z1},
		{-50.0f, 50.0f, z2}, {50.0f, 50.0f, z2}, {-50.0f, -50.0f, z2},{50.0f, -50.0f, z2},
		}
	);

	//cubeV.addVertex({
	//{-50.0f, 50.0f, z1}, {50.0f, 50.0f, z1}, {-50.0f, -50.0f, z1},
	//	}
	//);

	//cubeV.addFace({
	//	{0,1,2},{5,4,6} });


	cubeV.addFace({ {0,1,2},{1,3,2},{5,4,6},
		{7,5,6},{4,0,6},{6,0,2},
		{4,5,0},{5,1,0},{1,5,3},
		{3,5,7},{2,3,6},{3,7,6} });

	//move(cubeV, { 100,100,25 });


	//rotate(cubeV, { 1, 1, 1 }, 3.14159f * 193 / 250);
	//
	//rotate(importM, { 0, 1, 0 }, 3.14159f*50 / 50);

#pragma endregion


	int mouseX = 0;
	int mouseY = 0;
	bool gPressed = false;
	bool rPressed = false;
	bool shiftPressed = false;



	//meshv& renderedMesh = cubeV;
	meshv& renderedMesh = importM;

	// animation loop
	while (!close) {

		auto start = std::chrono::high_resolution_clock::now();

		reinitializeZBuffer(zBuffer, WIDTH, HEIGHT);

		SDL_SetRenderDrawColor(renderer, 255 * 0.227, 255 * 0.227, 255 * 0.227, 255);
		SDL_RenderClear(renderer);

		SDL_Event event;

		drawTris(renderedMesh);
		//drawTris(importM);

		//rotate(renderedMesh, { 1, 1, 1 }, 3.14159f / 250);
		//rotate(importM, { 0, 1, 0 }, 3.14159f / 50);



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
					move(renderedMesh, { 0,10,0 });
					//dest.y -= speed / 30;
					break;
				case SDL_SCANCODE_A:
				case SDL_SCANCODE_LEFT:
					move(renderedMesh, { 10,0,0 });
					break;
				case SDL_SCANCODE_S:
				case SDL_SCANCODE_DOWN:
					move(renderedMesh, { 0,-10,0 });
					break;
				case SDL_SCANCODE_D:
				case SDL_SCANCODE_RIGHT:
					move(renderedMesh, { -10,0,0 });
					break;

				case SDL_SCANCODE_G:
					if (!gPressed) {
						gPressed = true;
						SDL_GetMouseState(&mouseX, &mouseY);
					}
					cout << "G" << endl;
					//rotate(cubeV, { 0, 1, 0 }, 3.14159f / 250);

					break;

				case SDL_SCANCODE_R:
					if (!rPressed) {
						rPressed = true;
						SDL_GetMouseState(&mouseX, &mouseY);
					}
					cout << "R" << endl;

					break;
				case SDL_SCANCODE_ESCAPE:
					gPressed = false;
					rPressed = false;
					break;
				case SDL_SCANCODE_LSHIFT:
				case SDL_SCANCODE_RSHIFT:
					std::cout << "Shift key down" << std::endl;
					shiftPressed = true;
					break;
				default:
					break;
				}
				break;

			case SDL_KEYUP:

				switch (event.key.keysym.scancode) {

				case SDL_SCANCODE_G:
					cout << "G up" << endl;
					break;
				case SDL_SCANCODE_R:
					cout << "R up" << endl;
					break;
				case SDL_SCANCODE_LSHIFT:
				case SDL_SCANCODE_RSHIFT:
					std::cout << "Shift key up" << std::endl;
					shiftPressed = false;
					break;
				}

				break;

			case SDL_MOUSEWHEEL:
				cout << "mouse wheel " << event.wheel.y << endl;

				vec3d moveDir = { 0,0,30 * event.wheel.y };

				move(renderedMesh, moveDir);
				break;

			}

			
			switch (event.type) {
			case SDL_MOUSEBUTTONDOWN:
				cout << "mouse wheel down" << endl;
				SDL_GetMouseState(&mouseX, &mouseY);
				break;

			case SDL_MOUSEBUTTONUP:
				cout << "mouse wheel up" << endl;
				break;
			}

			if (event.type == SDL_MOUSEBUTTONDOWN) {
			}
		}

		//Move
		if (gPressed) {
			int x = 0;
			int y = 0;
			SDL_GetMouseState(&x, &y);

			float z = renderedMesh.position.z;
			float zRate = z / (fFar - fNear);

			vec3d moveDir = { (mouseX - x) / -zRate , (mouseY - y) / -zRate, 0 };
			mouseX = x;
			mouseY = y;
			move(renderedMesh, moveDir);

		}

		//Rotate
		if (rPressed) {
			int x = 0;
			int y = 0;
			SDL_GetMouseState(&x, &y);

			float z = renderedMesh.position.z;
			float zRate = z / (fFar - fNear);

			vec3d moveDir = { (mouseX - x) / -zRate , (mouseY - y) / -zRate, 0 };
			mouseX = x;
			mouseY = y;

			rotate(renderedMesh, { 0,1,0 }, 3.14 * (moveDir.x) / 50);
			//rotate(renderedMesh, { 0,0,1 }, 3.14 * (moveDir.y) / 50);

			//rotate(renderedMesh, { 0,1,0 }, 3.14 / 50);
		}

		if (SDL_GetMouseState(NULL, NULL) == 1) {
			gPressed = false;
			rPressed = false;
			cout << "mouse down" << endl;
		}

		//Pan camera
		if (SDL_GetMouseState(NULL, NULL) == 2 && shiftPressed) {
			int x = 0;
			int y = 0;
			SDL_GetMouseState(&x, &y);

			float z = renderedMesh.position.z;
			float zRate = z / (fFar - fNear);

			vec3d moveDir = { (mouseX - x) / -zRate , (mouseY - y) / -zRate, 0 };
			mouseX = x;
			mouseY = y;

			mainCamera.position = moveDir;

			move(renderedMesh, moveDir);
			
			//cout << "wheel pressed" << endl;
		}

		//Rotate camera
		if (SDL_GetMouseState(NULL, NULL) == 2) {
			int x = 0;
			int y = 0;
			SDL_GetMouseState(&x, &y);

			float z = renderedMesh.position.z;
			float zRate = z / (fFar - fNear);

			vec3d moveDir = { (mouseX - x) / -zRate , (mouseY - y) / -zRate, 0 };
			mouseX = x;
			mouseY = y;

			mainCamera.position = moveDir;

			rotateCamera(renderedMesh, mainCamera, { 0,1,0 }, 3.14 * (moveDir.x) / 50);
			rotateCamera(renderedMesh, mainCamera, { 1,0,0 }, -3.14 * (moveDir.y) / 50);
			
			//cout << "wheel pressed" << endl;
		}

		SDL_SetRenderDrawColor(renderer, 255 * 1, 255 * 1, 255 * 1, 255);

		//SDL_RenderDrawPoint(renderer, 0+a, 0+a);

		auto end = std::chrono::high_resolution_clock::now();

		// Calculate the duration
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		//rotate(importM, { 0, 1, 0 }, 3.14159f * duration.count() / 4000);
		// Output the duration in milliseconds
		//std::cout << "Time taken: " << duration.count() << " milliseconds" << std::endl;

		int x = 0;
		int y = 0;
		SDL_GetMouseState(&x, &y);

		//cout << y << endl;

		SDL_RenderPresent(renderer);

		SDL_Delay(1000 / 60);
		a++;
	}





	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();
	return 0;
}