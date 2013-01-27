#ifndef DATASTRUCT_H
#define DATASTRUCT_H

#include <fstream>
#include <iomanip>
#include <limits>
#include "Common.h"

class DataStructure {
	int ** tris, ** neigh, Ntris, Nverts;
	int GRIDPTS_X, GRIDPTS_Y, GRIDPTS, *is, *js, **id;
	double ** verts, * area;

public:
	// Structureless data structure
	DataStructure(char* fname);
	int GetNtris() const { return Ntris; }
	int GetNverts() const { return Nverts; }
	double GetVerts(const int i, const int j) const { return verts[i][j]; }
	int GetTris(const int i, const int j) const { return tris[i][j]; }
	int GetNeigh(const int i, const int j) const { return neigh[i][j]; }
	double GetArea(const int i) const { return area[i]; }
	
	// Structured (grid) data structure
	DataStructure(int Npts1, int Npts2);
	int GetNpts() const { return GRIDPTS; }
	int GetXNpts() const { return GRIDPTS_X; }
	int GetYNpts() const { return GRIDPTS_Y; }
	int operator()(int i, int j) const { return id[i][j]; }
	int operator[](int i) const { return is[i]; }
	int operator()(int j) const { return js[j]; }
	
	// Missing destructor
};

#endif
