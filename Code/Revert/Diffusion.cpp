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

sparse AssembleIsotropicDiffusion(const DataStructure& Data) {
	int Nverts = Data.GetNverts(), nnz = 0, NNZMAX = 0;

	for(int i = 0; i < Data.GetNverts(); i++)
		NNZMAX += Data.GetNeigh(i,0);

	sparse Diffusion(Nverts,Nverts,3*NNZMAX);

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

			// Tris(:,:) contains nodes in an anti-clockwise geomtric pattern
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

			double delta_x1 = x[2]/3.0 - x[1]/6.0 - x[0]/6.0,
				   delta_y1 = y[2]/3.0 - y[1]/6.0 - y[0]/6.0,

				   delta_x2 = x[2]/6.0 - x[1]/3.0 + x[0]/6.0,
				   delta_y2 = y[2]/6.0 - y[1]/3.0 + y[0]/6.0,

				   Area = Data.GetArea(Data.GetNeigh(i,j)),

				   div_x[3] = {(y[1] - y[2]) / (2.0*Area), 
					           (y[2] - y[0]) / (2.0*Area), 
					           (y[0] - y[1]) / (2.0*Area)},

        		   div_y[3] = {(x[2] - x[1]) / (2.0*Area),
					           (x[0] - x[2]) / (2.0*Area),
					           (x[1] - x[0]) / (2.0*Area)},

				   a1 = div_x[0] * (delta_y1 + delta_y2) - div_y[0] * (delta_x1 + delta_x2),
				   a2 = div_x[1] * (delta_y1 + delta_y2) - div_y[1] * (delta_x1 + delta_x2),
				   a3 = div_x[2] * (delta_y1 + delta_y2) - div_y[2] * (delta_x1 + delta_x2);

				   Diffusion.Row(nnz) = i;
				   Diffusion.Col(nnz) = Data.GetTris(Data.GetNeigh(i,j),index1);
				   Diffusion.Val(nnz) = - a1 / Area_C; 
				   nnz++;

				   Diffusion.Row(nnz) = i;
				   Diffusion.Col(nnz) = Data.GetTris(Data.GetNeigh(i,j),index2);
				   Diffusion.Val(nnz) = - a2 / Area_C;
				   nnz++;

				   Diffusion.Row(nnz) = i;
				   Diffusion.Col(nnz) = Data.GetTris(Data.GetNeigh(i,j),index3);
				   Diffusion.Val(nnz) = - a3 / Area_C;
				   nnz++;
        	}
	}

	Diffusion.CreateSparse();
	return Diffusion;
}

sparse BuildLaplacian(const DataStructure& Domain) {
	
	int Npts = Domain.GetNpts();
	sparse Laplace(Npts,Npts,5*Npts);

	int nnz = 0;

	for(int i = 0; i < Npts; i++) {
	
		Laplace.Val(nnz) = - 2.0 / (hx*hx) - 2.0 / (hy*hy);
		Laplace.Row(nnz) = i;
		Laplace.Col(nnz) = i;
		nnz++;

		if(Domain(Domain[i],Domain(i)+1) != -4) {
			Laplace.Val(nnz) = 1.0 / (hx*hx);
			Laplace.Row(nnz) = i;
			Laplace.Col(nnz) = Domain(Domain[i],Domain(i)+1);
			nnz++;
		}
		else {
			Laplace.Val(nnz) = 1.0 / (hx*hx);
			Laplace.Row(nnz) = i;
			Laplace.Col(nnz) = i;
			nnz++;
		}

		if(Domain(Domain[i],Domain(i) - 1) != -2) {
			Laplace.Val(nnz) = 1.0 / (hx*hx);
			Laplace.Row(nnz) = i;
			Laplace.Col(nnz) = Domain(Domain[i],Domain(i)-1);
			nnz++;
		}
		else {
			Laplace.Val(nnz) = 1.0 / (hx*hx);
			Laplace.Row(nnz) = i;
			Laplace.Col(nnz) = i;
			nnz++;
		}

		if(Domain(Domain[i] + 1,Domain(i)) != -3) {
			Laplace.Val(nnz) = 1.0/(hy*hy);
			Laplace.Row(nnz) = i;
			Laplace.Col(nnz) = Domain(Domain[i]+1,Domain(i));
			nnz++;
		}
		else {
			Laplace.Val(nnz) = 1.0/(hy*hy);
			Laplace.Row(nnz) = i;
			Laplace.Col(nnz) = i;
			nnz++;
		}

		if(Domain(Domain[i]-1,Domain(i)) != -1) {
			Laplace.Val(nnz) = 1.0/(hy*hy);
			Laplace.Row(nnz) = i;
			Laplace.Col(nnz) = Domain(Domain[i]-1,Domain(i));
			nnz++;
		}
		else {
			Laplace.Val(nnz) = 1.0/(hy*hy);
			Laplace.Row(nnz) = i;
			Laplace.Col(nnz) = i;
			nnz++;
		}
	}

	Laplace.CreateSparse();	
	return Laplace;
}
