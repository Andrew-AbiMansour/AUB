#include <fstream>
#include<iomanip>
#include <limits>
#include "Main.h"

#ifndef DATASTRUCT_H
#define DATASTRUCT_H

class DataStructure {
private:
	int ** tris, ** neigh, Ntris, Nverts;
	double ** verts, * area;

public:
	DataStructure(char* fname);
	int GetNtris() const { return Ntris; }
	int GetNverts() const { return Nverts; }
	double GetVerts(const int i, const int j) const { return verts[i][j]; }
	int GetTris(const int i, const int j) const { return tris[i][j]; }
	int GetNeigh(const int i, const int j) const { return neigh[i][j]; }
	double GetArea(const int i) const { return area[i]; }
};

#endif
