#pragma once


#ifndef VECTORCALC_H
#define VECTORCALC_H

// Cross product
vec3d cross(const vec3d& a, const vec3d& b) {
	return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
}
vec3d subVector(const vec3d& vec1, const vec3d& vec2) {
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

void normalize(vec3d& v) {

	float length = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);

	if (length == 0) {
		return;
	}

	v.x /= length;
	v.y /= length;
	v.z /= length;

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

#endif