#include "Main.h"
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

void AssembleTransient(sparse &Op, double Coef) {
	int index1 = 0, index2 = 0;

	for(int k = 0; k < Op.COLS() ; k++) {
		index2 = Op.ReturnOp()->p[k+1];
		for(int j = index1; j < index2; j++)
			if(Op.ReturnOp()->i[j] == k)
				Op.ReturnOp()->x[j] += Coef;
		index1 = index2;
	}
}

sparse AssembleSource(sparse &Op, Vector<double>& v) {
	int index1 = 0, index2 = 0;
	sparse Oper = Op;

	for(int k = 0; k < Op.COLS() ; k++) {
		index2 = Op.ReturnOp()->p[k+1];
		for(int j = index1; j < index2; j++)
			if(Op.ReturnOp()->i[j] == k)
				Oper.ReturnOp()->x[j] += kr * v(k);
		index1 = index2;
	}

	return Oper;
}

sparse AssembleNucleation(sparse &Op, vec& n, vec& p, vec& c, double Coef) {
	int index1 = 0, index2 = 0;
	sparse Oper = Op;

	for(int k = 0; k < Op.COLS(); k++) {
		bool theta1 = c(k) >= c1, theta2 = c(k) >= c2, theta3 = c(k) >= c3;

		double D1 = Coef + beta * c2 - beta * c(k) + (c(k) - c2) * theta2 * beta
					+ gama * theta2 * (c(k) - c2),

			   D2 = beta * n(k) * (theta2 - 1.0) + gama * theta2 * n(k) - alfa * theta3,
			   D3 = gama * theta2 * (c2 - c(k)),
			   D4 = Coef + delta * theta1 * (c1 - c(k)),
			   D5 = - delta * theta1 * p(k) - gama * theta2 * n(k),
			   D6 = beta * c(k) * (1.0 - theta2) + beta * c2 * (theta2 - 1.0),
			   D7 = delta * theta1 * (c(k) - c1),
			   L  = alfa * theta3 + delta * theta1 * p(k) + beta * n(k) * (1.0 - theta2),

			   V = - D6 * D2/D1 - D5 * D7/D4 + D2*D3*D7 / (D1*D4) + L;

		index2 = Op.ReturnOp()->p[k+1];
		for(int j = index1; j < index2; j++)
			if(Op.ReturnOp()->i[j] == k)
				Oper.ReturnOp()->x[j] += V;
		index1 = index2;
	}

	return Oper;
}

sparse AssembleNucleation_Picard(sparse &Op, vec& n, vec& p, vec& c) {
	int index1 = 0, index2 = 0;
	sparse Oper = Op;

	for(int k = 0; k < Op.COLS(); k++) {
		bool theta1 = c(k) >= c1, theta2 = c(k) >= c2, theta3 = c(k) >= c3;

		double V = alfa * theta3 + delta * theta1 * p(k) + beta * n(k) * (1.0 - theta2);

		index2 = Op.ReturnOp()->p[k+1];
		for(int j = index1; j < index2; j++)
			if(Op.ReturnOp()->i[j] == k)
				Oper.ReturnOp()->x[j] += V;
		index1 = index2;
	}

	return Oper;
}

sparse AssembleIsotropicDiffusion(DataStructure& Data) {

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
