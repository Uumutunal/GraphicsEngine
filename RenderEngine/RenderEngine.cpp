#include <SDL.h>
#include<iostream>
#include<vector>
#include<math.h>
#include <chrono>
#include <fstream>
#include <strstream>

using namespace std;


#include"methods.h"
#include"include/structs.h"
#include"include/vectorCalc.h"
#include"include/mesh.h"
#include"include/importer.h"



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


float fNear = 500;
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



void project(meshv& m) {

	for (size_t i = 0; i < m.vertices.size(); i++)
	{
		float z = m.position.z + m.vertices[i].p.z;

		float zRate = z / (fFar - fNear);


		m.verticesProjected[i].p.x = (m.position.x + m.vertices[i].p.x) * zRate;
		m.verticesProjected[i].p.y = (m.position.y + m.vertices[i].p.y) * zRate;
	}


	for (auto& f : m.faces) {

		vec3d v1 = subVector(m.vertices[f.points[1]].p, m.vertices[f.points[0]].p);
		vec3d v2 = subVector(m.vertices[f.points[2]].p, m.vertices[f.points[0]].p);
		vec3d normal = cross(v1, v2);
		normalize(normal);
		f.normal = normal;
	}
	int t = 0;

}



//TODO: rewrite
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


	}

}





void move(meshv& m, vec3d v) {


	m.position.x += v.x;
	m.position.y += v.y;
	m.position.z += v.z;


}



// Rotate point around an axis
vec3d rotateAroundPoint(vec3d& point, const vec3d mCenter, const vec3d& axis, float angle, bool camera) {
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

	point.x = p_rot.x + center.x;
	point.y = p_rot.y + center.y;
	point.z = p_rot.z + center.z;


	return { p_rot.x + center.x, p_rot.y + center.y, p_rot.z + center.z };
}


void rotate(meshv& m, const vec3d& axis, float angle) {


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

//Fix camera clipping

void drawGrid(camera cam) {

	SDL_SetRenderDrawColor(renderer, 100, 100, 100, 255);


	int range = cam.position.x / 100;

	for (int i = -WIDTH - range * 100 - 1000; i < WIDTH - range * 100 + 1000; i += 100)
	{

		vec3d lineStart = { i, 0, fFar + 200 };
		vec3d lineEnd = { i, 0, fNear + 4000 };

		lineStart.x += cam.position.x;
		lineStart.y += cam.position.y;
		lineStart.z += cam.position.z;

		lineEnd.x += cam.position.x;
		lineEnd.y += cam.position.y;
		lineEnd.z += cam.position.z;


		rotateAroundPoint(lineStart, cam.position, { 0,1,0 }, -cam.rotation.y, true);
		rotateAroundPoint(lineStart, { 0,0,500 }, { 1,0,0 }, -cam.rotation.x, true);

		rotateAroundPoint(lineEnd, cam.position, { 0,1,0 }, -cam.rotation.y, true);
		rotateAroundPoint(lineEnd, { 0,0,500 }, { 1,0,0 }, -cam.rotation.x, true);


		float z1 = lineEnd.z;
		float z2 = lineStart.z;

		float zRate1 = (fFar - fNear) / (z1);
		float zRate2 = (fFar - fNear) / (z2);

		if (zRate1 > 0) {
			continue;
		}
		if (zRate2 > 0) {
			continue;
		}

		lineEnd.x = lineEnd.x * zRate1;
		lineEnd.y = lineEnd.y * zRate1;

		lineStart.x = lineStart.x * zRate2;
		lineStart.y = lineStart.y * zRate2;



		SDL_RenderDrawLine(renderer, lineStart.x + WIDTH / 2, lineStart.y + HEIGHT / 2,
			lineEnd.x + WIDTH / 2, lineEnd.y + HEIGHT / 2);


	}

	for (int j = 0; j < 5000; j += 100)
	{
		vec3d lineStart = { -2000 - range * 10, 0, j };
		vec3d lineEnd = { 2000 - range * 10, 0, j };

		lineStart.x += cam.position.x;
		lineStart.y += cam.position.y;

		lineEnd.x += cam.position.x;
		lineEnd.y += cam.position.y;


		rotateAroundPoint(lineStart, cam.position, { 0,1,0 }, -cam.rotation.y, true);
		rotateAroundPoint(lineStart, { 0,0,500 }, { 1,0,0 }, -cam.rotation.x, true);

		rotateAroundPoint(lineEnd, cam.position, { 0,1,0 }, -cam.rotation.y, true);
		rotateAroundPoint(lineEnd, { 0,0,500 }, { 1,0,0 }, -cam.rotation.x, true);


		if (lineStart.z < 0) {
			float A;
			float B;
			float C;

			A = lineEnd.y - lineStart.y;
			B = (lineStart.x - lineEnd.x);
			C = (lineEnd.x * lineStart.y) - (lineStart.x * lineEnd.y);

			vec3d direction = { lineEnd.x - lineStart.x, lineEnd.y - lineStart.y, lineEnd.z - lineStart.z };

			float t = -lineStart.z / direction.z;

			float x = lineStart.x + t * direction.x;
			float y = lineStart.y + t * direction.y;
			float z = 0;

			//lineStart.x = x;
			//lineStart.y = y;

			//lineStart.x *= -1;
			//lineStart.y *= -1;
			lineStart.z *= -1;

		}

		float z1 = lineEnd.z;
		float z2 = lineStart.z;

		float zRate1 = (fFar - fNear) / (z1);
		float zRate2 = (fFar - fNear) / (z2);

		if (zRate1 > 0 || zRate2 > 0) {
			//continue;
		}



		lineEnd.x = lineEnd.x * zRate1;
		lineEnd.y = lineEnd.y * zRate1;

		lineStart.x = lineStart.x * zRate2;
		lineStart.y = lineStart.y * zRate2;



		//SDL_RenderDrawLine(renderer, lineStart.x, lineStart.y, lineEnd.x, lineEnd.y);
		SDL_RenderDrawLine(renderer, lineStart.x + WIDTH / 2, lineStart.y + HEIGHT / 2,
			lineEnd.x + WIDTH / 2, lineEnd.y + HEIGHT / 2);


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


#pragma region cube

	SDL_SetRenderDrawColor(renderer, 255 * 0.227, 255 * 0.227, 255 * 0.227, 255);
	SDL_RenderClear(renderer);

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
					//std::cout << "Shift key down" << std::endl;
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

				//Zoom
			case SDL_MOUSEWHEEL:
				cout << "mouse wheel " << event.wheel.y << endl;

				vec3d moveDir = { 0,0,30 * event.wheel.y };

				move(renderedMesh, moveDir);
				mainCamera.position.z = moveDir.z;
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

			mainCamera.position.x += moveDir.x;
			mainCamera.position.y += moveDir.y;

			cout << "pan camera" << endl;


			move(renderedMesh, moveDir);
		}

		//Rotate camera
		if (SDL_GetMouseState(NULL, NULL) == 2 && !shiftPressed) {
			int x = 0;
			int y = 0;
			SDL_GetMouseState(&x, &y);

			float z = renderedMesh.position.z;
			float zRate = z / (fFar - fNear);

			vec3d moveDir = { (mouseX - x) / -zRate , (mouseY - y) / -zRate, 0 };
			mouseX = x;
			mouseY = y;

			//mainCamera.position.x += moveDir.x;

			//cout << mainCamera.position.x << endl;

			rotateCamera(renderedMesh, mainCamera, { 0,1,0 }, 3.14 * (moveDir.x) / 50);
			rotateCamera(renderedMesh, mainCamera, { 1,0,0 }, -3.14 * (moveDir.y) / 50);

			mainCamera.rotation.x += -3.14 * (moveDir.y) / 50;
			mainCamera.rotation.y += -3.14 * (moveDir.x) / 50;

			cout << "rotate camera" << endl;
		}

		SDL_SetRenderDrawColor(renderer, 255 * 1, 255 * 1, 255 * 1, 255);


		auto end = std::chrono::high_resolution_clock::now();

		// Calculate the duration
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

		// Output the duration in milliseconds
		//std::cout << "Time taken: " << duration.count() << " milliseconds" << std::endl;

		int x = 0;
		int y = 0;
		SDL_GetMouseState(&x, &y);

		//drawGrid(mainCamera);


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