#include "DataStruct.h"

DataStructure::DataStructure(char* fname) {
	ifstream ReadFile;
	ReadFile.open(fname);
	
	ReadFile >> Nverts >> Ntris;
	ReadFile.ignore(numeric_limits<streamsize>::max(), '\n');
	
	tris = new int*[Ntris];
	neigh = new int*[Nverts];
	verts = new double*[Nverts];
	area = new double[Ntris];

	for(int i = 0; i < Nverts; i++)
		verts[i] = new double[3];

	for(int i = 0; i < Ntris; i++)
		tris[i] = new int[3];

	for(int i = 0; i < Nverts; i++) {
		for(int j = 0; j < 3; j++)
			ReadFile >> verts[i][j];
		ReadFile.ignore(numeric_limits<streamsize>::max(), '\n');
	}
	
	for(int i = 0; i < Ntris; i++) {
		for(int j = 0; j < 3; j++)
			ReadFile >> tris[i][j];
		ReadFile.ignore(numeric_limits<streamsize>::max(), '\n');
	}
	
	for(int i = 0; i < Nverts; i++) {
		int length; 
		ReadFile >> length;
		neigh[i] = new int[length+1];
		neigh[i][0] = length;

		for(int j = 1; j <= neigh[i][0]; j++)
			ReadFile >> neigh[i][j];
		ReadFile.ignore(numeric_limits<streamsize>::max(), '\n');
	}
	
	for(int i = 0; i < Ntris; i++) {
		double x1 = verts[tris[i][0]][0],
			   x2 = verts[tris[i][1]][0],
			   x3 = verts[tris[i][2]][0],

			   y1 = verts[tris[i][0]][1],
			   y2 = verts[tris[i][1]][1],
			   y3 = verts[tris[i][2]][1];

		area[i] = 0.5 * fabs( x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2 );
	}

}

DataStructure::DataStructure(int Npts1, int Npts2) : GRIDPTS_X(Npts1), GRIDPTS_Y(Npts2) {
		GRIDPTS = (GRIDPTS_X - 2) * (GRIDPTS_Y - 2);
		is = new int[GRIDPTS];
		js = new int[GRIDPTS];
		id = new int*[GRIDPTS_Y];

		for(int i = 0; i < GRIDPTS_Y; i++)
			id[i] = new int[GRIDPTS_X];

		int count = 0;

		for(int i = 0; i < GRIDPTS_Y; i++)
			for(int j = 0; j < GRIDPTS_X; j++)
				if(i == 0)
					id[i][j] = -1;
				else 
					if(j == 0)
						id[i][j] = -2;
					else
						if(i == (GRIDPTS_Y - 1))
							id[i][j] = -3;
						else
							if(j == (GRIDPTS_X - 1))
								id[i][j] = -4;
							else {
								is[count] = i;
								js[count] = j;
								id[i][j] = count++;
							}
}
