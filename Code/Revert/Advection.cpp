#include "Sparse.h"
#include "DataStruct.h"
#include "Vector.h"

extern double t, dt;

typedef Sparse<double> sparse;
typedef Vector<int> vector;
typedef Vector<double> vec;

template <typename T>
int find(Vector<T>& v, int c) {
	for(int i = 0; i < v.GetLength(); i++)
		if(v(i) == c)
			return i;

	return -1;
}

bool CheckBC(int index, const Vector<int>& BoundPoints) {
	bool IsBC = false;
	
	for(int i = 0; i < BoundPoints.GetLength(); i++)
		if(index == BoundPoints(i)) {
			IsBC = true;
			break;
		}
		
	return IsBC;
}

sparse AssembleAdvection(DataStructure& Data, double velocity, char* BCType) {
	/* Dirichlet BC is implemented inefficiently here */
	
	int Nverts = Data.GetNverts(), nnz = 0, NNZMAX = 0, nBC = 0, uBC;
	bool DirichletIsTrue = false;
	
	if(!strcmp(BCType,"Dirichlet"))
		DirichletIsTrue = true;
		
	for(int i = 0; i < Data.GetNverts(); i++)
		NNZMAX += Data.GetNeigh(i,0);

	sparse Advection(Nverts,Nverts,2*NNZMAX);
	ifstream file;
	
	if(DirichletIsTrue) {	
		file.open("Dirichlet.dat");
		file >> nBC >> uBC;
	}
		 
	vector BoundPts(nBC);
	
	if(DirichletIsTrue) 
		for(int i = 0; i < nBC; i++)
			file >> BoundPts(i);
			 
	file.close();
	
	for(int i = 0; i < Nverts; i++) {
		int n = Data.GetNeigh(i,0);
		double Area_C = .0;
		
		for(int k = 1; k <= n; k++)
			Area_C += Data.GetArea(Data.GetNeigh(i,k)); 

		Area_C *= 1.0/3.0;
		
		for(int j = 1; j <= n; j++) {
			double x[3], y[3]; vector nodes(3); 

			for(int k = 0; k < 3; k++)
				nodes(k) = Data.GetTris(Data.GetNeigh(i,j),k);

			int index1 = find(nodes,i), index2 = 1, index3 = 2;
	
			if(index1 == 1) {
				index2 = 2;
				index3 = 0;
			}

			if(index1 == 2) {
				index2 = 0;
				index3 = 1;
			}

			x[0] = Data.GetVerts(nodes(index1),0);
			y[0] = Data.GetVerts(nodes(index1),1);

			x[1] = Data.GetVerts(nodes(index2),0);
			y[1] = Data.GetVerts(nodes(index2),1);

			x[2] = Data.GetVerts(nodes(index3),0);
			y[2] = Data.GetVerts(nodes(index3),1);

			bool IsBoundary1 = CheckBC(nodes(index1),BoundPts),
				 IsBoundary2 = CheckBC(nodes(index2),BoundPts),
				 IsBoundary3 = CheckBC(nodes(index3),BoundPts);
				 
			double delta_x1 = x[2]/3.0 - x[1]/6.0 - x[0]/6.0,
				   delta_y1 = y[2]/3.0 - y[1]/6.0 - y[0]/6.0,

				   delta_x2 = x[2]/6.0 - x[1]/3.0 + x[0]/6.0,
				   delta_y2 = y[2]/6.0 - y[1]/3.0 + y[0]/6.0;
					
			double qf1 = velocity * (delta_y1 - delta_x1),
			       qf2 = velocity * (delta_y2 - delta_x2);
			
			if(qf1 > 0) {
				   Advection.Row(nnz) = i;
				   Advection.Col(nnz) = Data.GetTris(Data.GetNeigh(i,j),index1);
				   if(IsBoundary1)
						Advection.Val(nnz) = .0;
				   else
						Advection.Val(nnz) = - qf1 / Area_C;
				   nnz++;
			}
			else {
				   Advection.Row(nnz) = i;
				   Advection.Col(nnz) = Data.GetTris(Data.GetNeigh(i,j),index2);
				   if(IsBoundary2)
						Advection.Val(nnz) = .0;
				   else
						Advection.Val(nnz) = - qf1 / Area_C;
				   nnz++;
			}
			
			if(qf2 > 0) {
				   Advection.Row(nnz) = i;
				   Advection.Col(nnz) = Data.GetTris(Data.GetNeigh(i,j),index1);
				   if(IsBoundary1)
						Advection.Val(nnz) = .0;
				   else
						Advection.Val(nnz) = - qf2 / Area_C;
				   nnz++;
        	}
        	else {
				   Advection.Row(nnz) = i;
				   Advection.Col(nnz) = Data.GetTris(Data.GetNeigh(i,j),index3);
				   if(IsBoundary3)
						Advection.Val(nnz) = .0;
				   else
						Advection.Val(nnz) = - qf2 / Area_C;
				   nnz++;
        	}
        }
	}
	Advection.CreateSparse();
	return Advection;
}

sparse BuildAdvection(const DataStructure& Domain, double velocity) {
	
	int Npts = Domain.GetNpts();
	sparse Advection(Npts,Npts,3*Npts);

	int nnz = 0;

	if(velocity > .0)
		for(int i = 0; i < Npts; i++) {
		
			Advection.Val(nnz) = - 1.0 / hx - 1.0 / hy;
			Advection.Row(nnz) = i;
			Advection.Col(nnz) = i;
			nnz++;

			if(Domain(Domain[i],Domain(i)+1) != -4) {
				Advection.Val(nnz) = 1.0 / hx;
				Advection.Row(nnz) = i;
				Advection.Col(nnz) = Domain(Domain[i],Domain(i)+1);
				nnz++;
			}
			else {
				Advection.Val(nnz) = 1.0 / hx;
				Advection.Row(nnz) = i;
				Advection.Col(nnz) = i;
				nnz++;
			}

			if(Domain(Domain[i] + 1,Domain(i)) != -3) {
				Advection.Val(nnz) = 1.0/hy;
				Advection.Row(nnz) = i;
				Advection.Col(nnz) = Domain(Domain[i]+1,Domain(i));
				nnz++;
			}
			else {
				Advection.Val(nnz) = 1.0/hy;
				Advection.Row(nnz) = i;
				Advection.Col(nnz) = i;
				nnz++;
			}
	}
	else
		for(int i = 0; i < Npts; i++) {
		
			Advection.Val(nnz) = 1.0 / hx + 1.0 / hy;
			Advection.Row(nnz) = i;
			Advection.Col(nnz) = i;
			nnz++;
			
		if(Domain(Domain[i],Domain(i) - 1) != -2) {
			Advection.Val(nnz) = - 1.0 / hx;
			Advection.Row(nnz) = i;
			Advection.Col(nnz) = Domain(Domain[i],Domain(i)-1);
			nnz++;
		}
		else {
			Advection.Val(nnz) = - 1.0 / hx;
			Advection.Row(nnz) = i;
			Advection.Col(nnz) = i;
			nnz++;
		}

		if(Domain(Domain[i]-1,Domain(i)) != -1) {
			Advection.Val(nnz) = - 1.0/hy;
			Advection.Row(nnz) = i;
			Advection.Col(nnz) = Domain(Domain[i]-1,Domain(i));
			nnz++;
		}
		else {
			Advection.Val(nnz) = - 1.0/hy;
			Advection.Row(nnz) = i;
			Advection.Col(nnz) = i;
			nnz++;
		}
	}

	Advection.CreateSparse();	
	return Advection;
}
