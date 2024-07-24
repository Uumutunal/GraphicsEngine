#include <SDL.h>
#include<iostream>
#include<vector>
#include<math.h>
#include <chrono>
#include"methods.h"
#include <map>
#include <algorithm>

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
	vector<vec3d> pos;
	vec3d normal;
};



// Cross product
vec3d cross(const vec3d& a, const vec3d& b) {
	return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
}
vec3d subVector(vec3d vec1, vec3d vec2) {
	vec3d subbed;
	subbed.x = vec1.x - vec2.x;
	subbed.y = vec1.y - vec2.y;
	subbed.z = vec1.z - vec2.z;

	return subbed;
}
// Dot product
float dot(const vec3d& a, const vec3d& b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

struct meshv
{


private:
	int index = 0;

public:
	vector<vertex> vertices;

	map<int, vertex> verts;

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

			verts[index] = verticesA[i];

			index++;

		}
	}

	void addFace(vector<vector<int>> facesA) {
		for (size_t i = 0; i < facesA.size(); i++)
		{
			face f;
			f.points = facesA[i];
			vector<vec3d> posV;
			posV.push_back(verts[facesA[i][0]].p);
			posV.push_back(verts[facesA[i][1]].p);
			posV.push_back(verts[facesA[i][2]].p);

			vec3d v1 = subVector(verts[facesA[i][1]].p, verts[facesA[i][0]].p);
			vec3d v2 = subVector(verts[facesA[i][2]].p, verts[facesA[i][0]].p);

			vec3d normal = cross(v1, v2);

			verts[facesA[i][0]].normal = normal;
			verts[facesA[i][1]].normal = normal;
			verts[facesA[i][2]].normal = normal;

			f.normal = normal;

			f.pos = posV;
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

	//for (auto& face : m.faces) {
	for (int j = 0; j < m.faces.size(); j++) {
		for (size_t i = 0; i < 3; i++)
		{
			float z = m.position.z + m.faces[j].pos[0].z;

			float zRate = z / (fFar - fNear);

			//zRate = 1;

			m.faces[j].pos[i].x = (m.position.x + m.faces[j].pos[i].x) * zRate;
			m.faces[j].pos[i].y = (m.position.y + m.faces[j].pos[i].y) * zRate;

		}

	}
	int t = 0;

	return m;
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


		float l1 = (top.x - x) * (top.x - x) + (top.y - p1) * (top.y - p1);
		float l2 = (top.x - x) * (top.x - x) + (top.y - p2) * (top.y - p2);
		float l3 = (top.x - x) * (top.x - x) + (top.y - p3) * (top.y - p3);

		float A;
		float B;
		float C;

		if ((a * top.x + b * top.y + c) == 0 && (a * bot.x + b * bot.y + c) == 0) {
			A = a;
			B = b;
			C = c;
		}
		else if ((a2 * top.x + b2 * top.y + c2) == 0 && (a2 * bot.x + b2 * bot.y + c2) == 0) {
			A = a2;
			B = b2;
			C = c2;
		}
		else if ((a3 * top.x + b3 * top.y + c3) == 0 && (a3 * bot.x + b3 * bot.y + c3) == 0) {
			A = a3;
			B = b3;
			C = c3;
		}
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
				SDL_RenderDrawLines(renderer, points, 2);
			}
			else {
				SDL_Point points[2] = {
				{x + WIDTH / 2,p1 + HEIGHT / 2},
				{x + WIDTH / 2, p3 + HEIGHT / 2},
				};
				SDL_RenderDrawLines(renderer, points, 2);
			}

		}
		if (abs(A * x + B * p2 + C) <= 0.01) {
			float lH1 = abs(p2 - p1);
			float lH2 = abs(p2 - p3);


			if (lH1 < lH2) {
				SDL_Point points[2] = {
				{x + WIDTH / 2,p2 + HEIGHT / 2},
				{x + WIDTH / 2, p1 + HEIGHT / 2},
				};
				SDL_RenderDrawLines(renderer, points, 2);
			}
			else {
				SDL_Point points[2] = {
				{x + WIDTH / 2,p2 + HEIGHT / 2},
				{x + WIDTH / 2, p3 + HEIGHT / 2},
				};
				SDL_RenderDrawLines(renderer, points, 2);
			}
		}
		if (abs(A * x + B * p3 + C) <= 0.01) {
			float lH1 = abs(p3 - p2);
			float lH2 = abs(p3 - p1);


			if (lH1 < lH2) {
				SDL_Point points[2] = {
				{x + WIDTH / 2,p3 + HEIGHT / 2},
				{x + WIDTH / 2, p2 + HEIGHT / 2},
				};
				SDL_RenderDrawLines(renderer, points, 2);
			}
			else {
				SDL_Point points[2] = {
				{x + WIDTH / 2,p3 + HEIGHT / 2},
				{x + WIDTH / 2, p1 + HEIGHT / 2},
				};
				SDL_RenderDrawLines(renderer, points, 2);
			}
		}

		//***************

		//if (l1 >= l2 && l1 >= l3) {

		//	if (p2 < -99999999 || p2 > 99999999) {
		//		p2 = p1;
		//	}
		//	if (p3 < -99999999 || p3 > 99999999) {
		//		p3 = p1;
		//	}

		//	SDL_Point points[2] = {
		//	{x + WIDTH / 2,p2 + HEIGHT / 2},
		//	{x + WIDTH / 2, p3 + HEIGHT / 2},
		//	};
		//	SDL_RenderDrawLines(renderer, points, 2);
		//}
		//else if (l2 >= l3 && l2 >= l1) {
		//	if (p3 < -999999 || p3 > 999999) {
		//		p3 = p2;
		//	}
		//	if (p1 < -999999 || p1 > 999999) {
		//		p1 = p2;
		//	}
		//	SDL_Point points[2] = {
		//	{x + WIDTH / 2,p3 + HEIGHT / 2},
		//	{x + WIDTH / 2, p1 + HEIGHT / 2},
		//	};
		//	SDL_RenderDrawLines(renderer, points, 2);
		//}
		//else if (l3 >= l2 && l3 >= l1) {
		//	if (p2 < -999999 || p2 > 999999) {
		//		p2 = p3;
		//	}
		//	if (p1 < -999999 || p1 > 999999) {
		//		p1 = p3;
		//	}
		//	SDL_Point points[2] = {
		//	{x + WIDTH / 2,p2 + HEIGHT / 2},
		//	{x + WIDTH / 2, p1 + HEIGHT / 2},
		//	};
		//	SDL_RenderDrawLines(renderer, points, 2);
		//}
	}

	return;



#pragma region old
	float xS = longestEdge[0].x;
	float xF = longestEdge[1].x;

	if (longestEdge[0].x > longestEdge[1].x) {
		xS = longestEdge[1].x;
		xF = longestEdge[0].x;
	}

	for (float x = xS; x < xF; x++)
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


		float xS2 = x;
		float xF2 = thirdVertex.x;

		if (x > thirdVertex.x) {
			xS2 = thirdVertex.x;
			xF2 = x;
		}


		for (float i = xS2; i < xF2; i++)
		{
			vec3d p;
			p.x = i;
			p.y = (aP * x + cP) / -bP;
			SDL_RenderDrawPoint(renderer, p.x + WIDTH / 2, p.y + HEIGHT / 2);
			//Fix <> signs
			if ((a2 * p.x + b2 * p.y + c2) * s2 < 0 && (a3 * p.x + b3 * p.y + c3) * s3 < 0) {

				SDL_RenderDrawPoint(renderer, p.x + WIDTH / 2, p.y + HEIGHT / 2);

			}
			else {
				break;
			}
		}

	}
#pragma endregion




}

//TODO: fix face-vertex
void drawTris(meshv m) {


	int t = 0;
	project(m);
	int t2 = 0;

	map<int, vec3d> f;

	for (auto face : m.faces) {


		for (size_t i = 0; i < 3; i++)
		{

			int index = face.points[i];

			for (auto v : m.vertices) {

				if (v.id == index) {
					f[index] = v.p;
					int t = 0;
					break;
				}

			}
		}
	}

	for (auto& face : m.faces) {

		face.pos[0] = f[face.points[0]];
		face.pos[1] = f[face.points[1]];
		face.pos[2] = f[face.points[2]];

	}

	sort(m.faces.begin(), m.faces.end(), [](face &f1, face &f2) {

		float z1 = (f1.pos[0].z + f1.pos[1].z + f1.pos[2].z) / 3.0;
		float z2 = (f2.pos[0].z + f2.pos[1].z + f2.pos[2].z) / 3.0;
		return z1 < z2;
	});


	

	//SDL_SetRenderDrawColor(renderer, 0, 255, 255, 255);

	vector<vec3d> pos;

	for (auto face : m.faces) {

		SDL_SetRenderDrawColor(renderer, 0, 255, 255, 255);

		int j;
		pos.clear();
		vec3d abc = { 0,0,-1 };
		//float d = dot(face.normal, abc);
		float d = face.normal.z;

		if (d < 0) {
			continue;
			SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);
			int tt = 0;
		}
		if (face.normal.y > 0) {
			SDL_SetRenderDrawColor(renderer, 0, 125, 125, 125);

		}

		for (size_t i = 0; i < 3; i++)
		{

			int index = face.points[i];

			//m.verts[index].p;
			//pos.push_back(m.verts[index].p);

			for (auto v : m.vertices) {

				if (v.id == index) {
					//v.p.y = -v.p.y;
					pos.push_back(v.p);
					int t = 0;
					break;
				}

			}

		}



		//SDL_SetRenderDrawColor(renderer, 0, 255, 255, 255);
		fillTriangle(face.pos);
		//fillTriangle(pos);

		for (size_t i = 0; i < 3; i++) {
			if (i != 2) {
				j = i + 1;
			}
			else {
				j = 0;
			}


			//SDL_Point points[2] = {
			//{pos[i].x + WIDTH / 2,pos[i].y + HEIGHT / 2},
			//{pos[j].x + WIDTH / 2, pos[j].y + HEIGHT / 2},
			//};

			SDL_Point points[2] = {
			{face.pos[i].x + WIDTH / 2,face.pos[i].y + HEIGHT / 2},
			{face.pos[j].x + WIDTH / 2, face.pos[j].y + HEIGHT / 2},
			};
			SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);

			SDL_RenderDrawLines(renderer, points, 2);
		}

	}

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
	meshv mP = m;
	project(mP);
	for (auto& vert : m.vertices) {

		vec3d rotated = rotateAroundPoint(vert.p, m.position, axis, angle);

		vert.p.x = rotated.x;
		vert.p.y = rotated.y;
		vert.p.z = rotated.z;
	}
	int t = 0;
	for (auto& f : m.faces) {

		vec3d v1 = subVector(mP.vertices[f.points[1]].p, mP.vertices[f.points[0]].p);
		vec3d v2 = subVector(mP.vertices[f.points[2]].p, mP.vertices[f.points[0]].p);
		vec3d normal = cross(v1, v2);
		f.normal = normal;
	}
	int ta = 0;

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

	//cubeV.addVertex({
	//{-50.0f, 50.0f, z1}, {50.0f, 50.0f, z1}, {-50.0f, -50.0f, z1},
	//	}
	//);

	//cubeV.addFace({ 
	//	{3,5,7} });


	cubeV.addFace({ {0,1,2},{1,3,2},{5,4,6},
		{7,5,6},{4,0,6},{6,0,2},
		{4,5,0},{5,1,0},{1,5,3},
		{3,5,7},{2,3,6},{3,7,6} });

	//move(cubeV, { 100,100,25 });


	//rotate(cubeV, { 1, 1, 1 }, 3.14159f * 193 / 250);
	//

#pragma endregion



	// animation loop
	while (!close) {

		auto start = std::chrono::high_resolution_clock::now();

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
					rotate(cubeV, { 1, 1, 1 }, 3.14159f / 250);
					break;
				default:
					break;
				}
			}

		}



		auto end = std::chrono::high_resolution_clock::now();

		// Calculate the duration
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

		// Output the duration in milliseconds
		//std::cout << "Time taken: " << duration.count() << " milliseconds" << std::endl;


		SDL_RenderPresent(renderer);

		SDL_Delay(1000 / 60);
		a++;
	}





	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();
	return 0;
}