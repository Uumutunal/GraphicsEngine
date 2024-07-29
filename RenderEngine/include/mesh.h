#pragma once



#ifndef MESH_H
#define MESH_H


struct meshv
{

private:
	int index = 0;

public:
	vector<vertex> vertices;
	vector<vertex> verticesProjected;

	map<int, vertex> verts;

	vector<face> faces;

	vec3d cameraPos;

	vec3d position;


	void addVertex(vector<vertex> verticesA) {
		for (size_t i = 0; i < verticesA.size(); i++)
		{
			vertex v;
			v.p = verticesA[i].p;
			v.id = index;
			vertices.push_back(v);
			verticesProjected.push_back(v);
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

			normalize(normal);

			verts[facesA[i][0]].normal = normal;
			verts[facesA[i][1]].normal = normal;
			verts[facesA[i][2]].normal = normal;

			f.normal = normal;

			f.pos = posV;
			faces.push_back(f);
		}
	}

};


#endif