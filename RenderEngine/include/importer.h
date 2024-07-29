#pragma once




#ifndef IMPORTER_H
#define IMPORTER_H



bool importMesh(string sFilename, meshv& importedMesh) {

	//meshv importedMesh;


	ifstream f(sFilename);
	if (!f.is_open())
		return false;

	// Local cache of verts
	vector<vertex> verts;
	vector<vec3d> vertsP;
	vector<face> faces;
	map<int, vertex> vv;

	vector<vector<int>> facesA;

	int index = 0;
	while (!f.eof())
	{
		char line[128];
		f.getline(line, 128);

		strstream s;
		s << line;

		char junk;

		if (line[0] == 'v')
		{
			vertex vert;
			vec3d v;
			s >> junk >> v.x >> v.y >> v.z;
			vert.p = v;
			vert.id = index;
			verts.push_back(vert);
			vertsP.push_back(v);
			vv[index] = vert;
			index++;
		}

		if (line[0] == 'f')
		{
			face ff;

			int f[3];
			s >> junk >> f[0] >> f[1] >> f[2];

			ff.points.push_back(f[2] - 1);
			ff.points.push_back(f[1] - 1);
			ff.points.push_back(f[0] - 1);

			facesA.push_back(ff.points);

			ff.pos.push_back(vertsP[f[2] - 1]);
			ff.pos.push_back(vertsP[f[1] - 1]);
			ff.pos.push_back(vertsP[f[0] - 1]);

			vec3d v1 = subVector(vertsP[f[1] - 1], vertsP[f[0] - 1]);
			vec3d v2 = subVector(vertsP[f[2] - 1], vertsP[f[0] - 1]);
			vec3d normal = cross(v1, v2);
			normalize(normal);
			ff.normal = normal;

			faces.push_back(ff);
		}
	}

	importedMesh.addVertex(verts);
	importedMesh.addFace(facesA);

	importedMesh.position = { 0,0,500 };

	return true;
}


#endif