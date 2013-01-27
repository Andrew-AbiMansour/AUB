#include "Main.h"
#include "Matrix.h"
#include "DataStruct.h"

typedef Matrix<double> matrix;
typedef Vector<Matrix<double> > tensor;

void LoadIC(tensor& x) {

	ifstream fp("input.dat");
	
	if(fp.is_open()) {
		for(int i = 0; i < Nc; i++) {
			for(int j = 0; j < x(i).GetCols(); j++)
				fp >> x(i)[0](j);

			x(i)[1] = x(i)[0];
		}
		fp.close();
	}
}

void LoadIC(matrix& c, double c0, double c1, DataStructure& Data, char* TYPE) {
	
	if(!strcmp(TYPE,"Circle"))
		for(int i = 0; i < c.GetCols(); i++) {
			double x = Data.GetVerts(i,0), y = Data.GetVerts(i,1),
			       x0 = 0.0, y0 = 0.0, r0 = sqrt(200.0);

			if(pow(x-x0,2.0) + pow(y-y0,2.0) <= pow(r0,2.0))
				c[0](i) = c0;
			else
				c[0](i) = c1;
		}

	if(!strcmp(TYPE,"Square"))
                for(int i = 0; i < c.GetCols(); i++) {
                        double x = Data.GetVerts(i,0), y = Data.GetVerts(i,1),
                               l0 = sqrt(pi * 200.0);

                        if(fabs(x) <= l0/2.f && fabs(y) <= l0/2.f)
                                c[0](i) = c0;
                        else
                                c[0](i) = c1;
                }

	if(!strcmp(TYPE,"Two-Circles"))
		for(int i = 0; i < c.GetCols(); i++) {
			double x = Data.GetVerts(i,0), y = Data.GetVerts(i,1),
				   x0 = 3.0, y0 = 2.0, r0 = sqrt(1.0);
		
			if(pow(x-x0,2.0) + pow(y+y0,2.0) <= pow(r0,2.0) || 
			   pow(x+x0,2.0) + pow(y-y0,2.0) <= pow(r0,2.0))
				c[0](i) = c0;
			else
				c[0](i) = c1;
		}
		 
	if(!strcmp(TYPE,"Multiple-Circles")) {

		int NumCirclesX = 2, NumCirclesY = 2;
		double fact = 1.0;
		double r0 = sqrt(150.0), x0 = - 100.0 * fact, y0 = x0, 

 			   SpacingX = (2.0 * fabs(x0) - 2.0 * NumCirclesX * r0) / (NumCirclesX + 1.0),
			   SpacingY = (2.0 * fabs(y0) - 2.0 * NumCirclesY * r0) / (NumCirclesY + 1.0);
		
		for(int i = 0; i < c.GetCols(); i++) {		
			bool IsInterior = false;
			double x = Data.GetVerts(i,0), y = Data.GetVerts(i,1);
	
			for(int j = 0; j < NumCirclesY; j++)
				for(int k = 0; k < NumCirclesX; k++)
if( (pow(x -(x0  +  SpacingX+r0+k*(2*r0+SpacingX)),2.0)
    +pow(y -(y0  + SpacingY+r0+j*(2*r0+SpacingY)),2.0)) <= pow(r0,2.0) )
						IsInterior = true;

			if(IsInterior)
				c[0](i) = c0;
			else
				c[0](i) = c1;
		}
	}

	if(!strcmp(TYPE,"Triangle"))
		for(int i = 0; i < c.GetCols(); i++) {
			double x = Data.GetVerts(i,0), y = Data.GetVerts(i,1),
				   x0 = sqrt(pi*100.0), y0 = x0 * cos(pi/6.0), trans = 0.0;
			
			if( (x >= -x0 && x <= x0) && ( y >= (-y0 - trans) && y <= (y0 - trans)) )
				if( (y + y0 + trans) <= (x0 - fabs(x)) * tan(pi/3.0))
					c[0](i) = c0;
				else
					c[0](i) = c1;
			else
				c[0](i) = c1;
		}
	

	if(!strcmp(TYPE,"Multiple-Triangles")) {

		int NumTrisX = 1, NumTrisY = 1;

		double L = 50.0, l = sqrt(20.0),
			   SpacingX = (L - NumTrisX * l) / (NumTrisX + 1.f),
			   SpacingY = (L - NumTrisY * sqrt(3.f) * l / 2.f) / (NumTrisY + 1.f);
		
		for(int i = 0; i < c.GetCols(); i++) {		
			bool IsInterior = false;
			double x = Data.GetVerts(i,0), y = Data.GetVerts(i,1), dx;
	
			for(int j = 0; j < NumTrisY; j++)
				for(int k = 0; k < NumTrisX; k++)
					if(x >= (SpacingX*(k+1) + k*l) && x <= (SpacingX*(k+1) + (k+1)*l)) {
						if(x >= (SpacingX*(k+1) + k*l))
							dx = ((SpacingX*(k+1) + (k+1)*l) - x);
						else
							dx = x - (SpacingX*(k+1) + k*l);

						if(y >= (SpacingY*(j+1) + j*l*sqrt(3.f)/2.f) && y <= (tan(pi/4.f) * dx 
							+ SpacingY*(j+1) + j*l*sqrt(3.f)/2.f))
							IsInterior = true;
						}
			if(IsInterior)
				c[0](i) = c0;
			else
				c[0](i) = c1;
		}
	}

	c[1] = c[0];
}
