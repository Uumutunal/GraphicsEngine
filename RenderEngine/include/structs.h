#pragma once


using namespace std;

#include <map>
#include<vector>



#ifndef STRUCTS_H
#define STRUCTS_H



void print2() {
	std::cout << "header structs" << std::endl;
}


struct pixel
{
	float depth;
};

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

struct camera {
	vec3d position = { 0,0,0 };
	vec3d rotation = { 0,0,0 };
};











#endif