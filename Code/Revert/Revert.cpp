#include "Common.h"
#include "Matrix.h"
#include "Vector.h"
#include "DataStruct.h"
#include "Sparse.h"

double t(0.0), dt0 = 0.0001, PERIOD;
const double tmax = 1.0, tmin = tmax, dc = 0.05;
typedef Matrix<double> matrix;
typedef Vector<double> vector;
typedef Sparse<double> sparse;
typedef Vector<Matrix<double> > tensor;

vector ConstructBDF(int, const vector&);
sparse AssembleIsotropicDiffusion(const DataStructure&);
sparse BuildLaplacian(const DataStructure&);
sparse BuildAdvection(const DataStructure&, double velocity);
void LoadIC(matrix& c, double c0, double c1, DataStructure& Data,char* TYPE);
void ReLoadIC(tensor& x);
sparse AssembleAdvection(DataStructure& Data, double velocity, char*);
void FixRHS(vector& RHS);
double GlobRate(const DataStructure&, const vector&);

vector ConstructTheta(double cr, vector& c) {
	vector theta(c.GetLength(),.0);

	for(int i = 0; i < c.GetLength(); i++)
		if(c(i) >= cr)
			theta(i) = 1.0;

	return theta;
}

double FixedPointStep(const DataStructure& Mesh, const sparse& Diff, const sparse& Adv1, const sparse& Adv2, const matrix& RHS, 
					  vector& x1, vector& x2, vector& x3, vector& x4, vector& x5, double coef,char* BCType, int Npts)
{
	vector  theta1 = ConstructTheta(c1,x3),
			theta2 = ConstructTheta(c2,x3),
			theta3 = ConstructTheta(c3,x3),

		RHS1 = RHS[0],
		RHS2 = RHS[1],
		RHS3 = RHS[2] + kr * x1 * x2 + alfa * c3 * theta3 + delta * c1 * theta1 * x5 + beta * (1.0 - theta2) * c2 * x4,
		RHS4 = RHS[3] + (x3 - c3) * alfa * theta3,
		RHS5 = RHS[4] + gama * (x3 - c2) * theta2 * x4;
	
	if(!strcmp(BCType,"Dirichlet")) {
		FixRHS(RHS1);
		FixRHS(RHS2);
	}
		
	sparse Op1 = - (Diff * Da) + (Adv1 * velocity1), Op2 = - (Diff * Db) + (Adv2 * velocity2), Op3 = - Diff * Dc;

	vector v1 = coef + kr * x2,
	       v2 = coef + kr * x1,
	       v3 = coef + alfa * theta3 + delta * theta1 * x5 + beta * (1.0 - theta2) * x4;
 
	Op1.Diag(v1);
	Op2.Diag(v2);
	Op3.Diag(v3);

	vector sol1 = Op1 / RHS1, //Op1.PCCG(RHS1,x1),
	       sol2 = Op2 / RHS2, //Op2.PCCG(RHS2,x2),
	       sol3 = Op3 / RHS3, //Op3.PCCG(RHS3,x3),
	       sol4(Npts,.0),
	       sol5(Npts,.0),

	       Op4 = coef + gama * (x3 - c2) * theta2 + beta * (1.0 - theta2) * (c2 - x3),
	       Op5 = coef - delta * (x3 - c1) * theta1;

	sol4 = RHS4 / Op4;
	sol5 = RHS5 / Op5;

	double error = max((sol1 - x1).Amax(),max((sol2 - x2).Amax(),max((sol3 - x3).Amax(),max((sol4 - x4).Amax(),(sol5 - x5).Amax()))));

	x1 = sol1; x2 = sol2; x3 = sol3; x4 = sol4; x5 = sol5;

	return error;
}

double FPI_RD(const DataStructure& Mesh, sparse& Diff, matrix RHS, vector& x1, vector& x2, double coef, double dt)
{
	int Npts = Mesh.GetNverts();

	vector	RHS1 = RHS[0],
		RHS2 = RHS[1];

	sparse Op1 = Diff * (Da * dt), Op2 = Diff * (Db * dt);

	vector v1 = coef + (dt * kr) * x2,
	       v2 = coef + (dt * kr) * x1;
 
	Op1.Diag(v1);
	Op2.Diag(v2);

	vector sol1 = Op1.PCCG(RHS1,x1),
	       sol2 = Op2.PCCG(RHS2,x2);

	double error = max((sol1 - x1).Amax(),(sol2 - x2).Amax());

	x1 = sol1; x2 = sol2;

	return error;
}

double NTR_NU(const DataStructure& Mesh, sparse& Diff, matrix RHS, vector& x1, vector& x2, vector& x3, vector& x4, 
	vector& x5, double coef, double dt)
{
	int Npts = Mesh.GetNverts();
	sparse Opc = Diff * (Dc * dt);

	vector  theta1 = ConstructTheta(c1,x3),
		theta2 = ConstructTheta(c2,x3),
		theta3 = ConstructTheta(c3,x3),
		z(Npts,.0),

		RHS3 = RHS[2] - coef * x3 - Opc * x3 + dt * (kr * x1 * x2 - alfa * theta3 * (x3 - c3) 
		     - delta * theta1 * (x3 - c1) * x5 + beta * (1.0 - theta2) * (c2 - x3) * x4),
		RHS4 = RHS[3] - coef * x4 + dt * ((x3 - c3) * alfa * theta3 - beta * (1.0 - theta2) * (x3 - c2) * x4
		     - gama * theta2 * (x3 - c2) * x4),
		RHS5 = RHS[4] - coef * x5 + dt * (gama * (x3 - c2) * theta2 * x4 + delta * (x3 - c1) * theta1 * x5);

 	vector  J11 = coef + dt * (alfa * theta3 + delta * theta1 * x5 + beta * (1.0 - theta2) * x4),
	        J12 = dt * (- beta * (c2 - x3) * (1.0 - theta2)),
	        J13 = dt * (delta * (x3 - c1) * theta1),
		J21 = dt * (- alfa * theta3 + beta * (1.0 - theta2) * x4 + gama * theta2 * x4),
	        J22 = coef + dt * (beta * (c2 - x3) * (1.0 - theta2) + gama * theta2 * (x3 - c2)),
		J31 = dt * (- gama * theta2 * x4 - delta * theta1 * x5),
		J32 = dt *(- gama * (x3 - c2) * theta2),
		J33 = coef - dt * (delta * (x3 - c1) * theta1),
 
		 V = J11 - J12 * J21 / J22 - J13 * J31 / J33 + J13 * J32 * J21 / (J22 * J33);
			   
	//cout << (Diff.Diag()).max() << " " << (Diff.Diag()).min() << endl;
	Opc.Diag(V);
	//Opc = PC * Opc;

	vector   RHSc = (RHS3 + (J12 / J22 - J13 * J32 / (J22*J33)) * RHS4 + J13 * RHS5 / J33),
		 dc = Opc.PCCG(RHSc,z),
		 dn = -  (- RHS4 + J21 * dc) / J22,
		 dp = -  (- RHS5 + J31 * dc + J32 * dn) / J33;

	double error = max(dc.Amax(),max(dn.Amax(),dp.Amax()));
		//cout << x3.Norm() << " " << x4.Norm() << " " << x5.Norm() << endl;

	x3 += dc; x4 += dn; x5 += dp;
	//cout << "Cond num = " << (Opc.Diag()).max() / (Opc.Diag()).min() << endl;
	return error;
}

bool NewtonSolver(const DataStructure& Mesh, const sparse& Diff, const sparse& Adv1, const sparse& Adv2, tensor& x, vector& dt, int Order,
				  char* BCType, int Npts) 
{
	int theta, GOrder, max_iters = 20;

	if(t >= (Order-1) * dt0)
		GOrder = Order;
	else {
		theta = (t + dt0)/dt0;
		GOrder = theta;
	}

	vector BDF1 = ConstructBDF(GOrder,dt), BDF2(GOrder,.0);

	matrix RHS1(Nc,Npts), RHS2(Nc,Npts);

	for(int j = 0; j < Nc; j++)
		for(int i = 1; i < GOrder + 1; i++)
			RHS1[j] -= BDF1(i) * x(j)[i];

	if(GOrder > 1) {
		 BDF2 = ConstructBDF(GOrder-1,dt);

		 for(int j = 0; j < Nc; j++)
			for(int i = 1; i < GOrder; i++)
				 RHS2[j] -= BDF2(i) * x(j)[i];
	}
	
	vector x1 = x(0)[0], x2 = x(1)[0], x3 = x(2)[0], x4 = x(3)[0], x5 = x(4)[0];

	double error = 1.0, tol = pow(10.0,-5.0), est = .0, abs_error = 0.001, rel_error = 0.001;
	int iters = 0;
	
	while(error > tol) {
		iters++;
		error = FixedPointStep(Mesh,Diff,Adv1,Adv2,RHS1,x1,x2,x3,x4,x5,BDF1(0),BCType,Npts);
		if(iters > max_iters)
			break;
	}
	
	x(0)[Order+2] = x1; x(1)[Order+2] = x2; x(2)[Order+2] = x3; x(3)[Order+2] = x4; 
	x(4)[Order+2] = x5;
	
	if(GOrder > 1 && iters <= max_iters) {
		vector err(Nc*Npts,.0);
		matrix deriv1 = RHS1, deriv2 = RHS2;
		
		for(int i = 0; i < Nc; i++) {
			deriv1[i] -= BDF1(0) * x(i)[Order+2];
			deriv2[i] -= BDF2(0) * x(i)[Order+2];
		}
			
		for(int i = 0; i < Nc; i++)
			for(int j = 0; j < Npts; j++)
				if(fabs(deriv1[i](j)) > pow(10.0,-8.0))
					err(i*Npts+j) = (fabs(deriv1[i](j)) - fabs(deriv2[i](j))); // / fabs(deriv1[i](j));

		est = err.Norm(); double deno = .0;
		
		for(int i = 0; i < Nc; i++)
			deno += deriv1[i].Norm(); est = est / deno; 
	}

	if(iters > max_iters) {
		cout << "Failure to converge with attempted timestep = " << dt(0) << endl;
		est = rel_error * 2.0; }

	if( est < rel_error) {

		for(int i = 0; i < Nc; i++)
			x(i)[0] = x(i)[Order+2];

		for(int i = 0; i < Nc; i++)
			for(int j = Order; j > 0; j--)
				x(i)[j] = x(i)[j-1];

		t += dt(0);
		vector rate = x(0)[0]*x(1)[0];
		vector time(1,t+dt(0)), rg(1,GlobRate(Mesh,rate)), cm(1,(rate.Imax()-(Ny-1)/2)*hy);
		
		time.Write("adaptive.dat");
		cm.Write("cmax.dat");
		rg.Write("GlobRate.dat");
		
		for(int i = Order - 1; i > 0; i--)
			dt(i) = dt(i-1);

		cout << "Newton Adaptive BDF Solver Converged with c_max = " << cm << " at time " << t << endl; 
	}

	if(GOrder > 1)
		dt(0) = min(tmax,0.8 * dt(0) *  pow(rel_error/est,1.0/(GOrder)));

	return 0;
}

int main(int argc, char* argv[]) {
	
	if(argc < 4) {

		cout << "====================================================================" << endl;
		cout << "====================================================================" << endl;
		cout << "CVFEM PDE solver for Lisegang systems, written by Andrew Abi Mansour" << endl;
		cout << "------------------------------------------------------------------" << endl;
		cout << "Program in Computational Science, American University of Beirut" << endl;
		cout << "====================================================================" << endl;
		cout << "====================================================================" << endl;
		cout << "Usage: ./PDES InitialCond Period InitialDist [OutputFname] [Mesh] [Order]" << endl << endl;
		cout << "Required Args:" << endl << endl;
		cout << "InitialCond = a string which initializes the simulation if" << endl 
		     << "              'Initialize' is selected, otherwise the simulation" << endl
		     << "              is resumed based on an input file input.dat" << endl;
		cout << "Period      = a floating point number that specifies the time at" << endl
		     << "              which the simulation ends." << endl;
		cout << "InitialDist = a string that specifies how each components" << endl
		     << "              is initially distributed." << endl << endl;	
		cout << "Optional Args:" << endl << endl;
		cout << "Order       = an integer specifying the maximum BDF order." << endl; 
		cout << "OutputFname = a string that specifies the output filename," << endl
		     << "              by defaul it is output.dat" << endl;
		cout << "Mesh        = a string that specifies the mesh data structure" << endl
			 << "              file (for FEM/CVFEM only)" << endl; 
		exit(1);
	}
	else
		{
			char* InitCond = argv[1];
			double PERIOD = atof(argv[2]);
			char* InitialDist = argv[3];
			
			int BDFOrder = 4;
			if(argc >= 5)
				BDFOrder = atoi(argv[4]);
				
			char* OutputFname = "output.dat";
			if(argc >= 6)
				OutputFname = argv[5];

			DataStructure Mesh(Nx,Ny);
			
			int Npts = Mesh.GetNpts(); //Mesh.GetNverts();
			char* BCType = "Neumann";
	
			/*sparse Diff = AssembleIsotropicDiffusion(Mesh), Adv1 = AssembleAdvection(Mesh,velocity1,BCType),
				   Adv2 = AssembleAdvection(Mesh,velocity2,BCType);*/
				   
			sparse Diff = BuildLaplacian(Mesh), Adv1 = BuildAdvection(Mesh,velocity1),
				   Adv2 = BuildAdvection(Mesh,velocity2);
			
			cout << "BDF order : " << BDFOrder << endl << "Simulation time : " << PERIOD << endl;

			if(!strcmp(InitCond,"Initialize"))
				cout << "Initializing simulation ..." << endl;
			else
				cout << "Resuming simulation ..." << endl;
		
			vector dt(BDFOrder,dt0);

			matrix m(BDFOrder+3,Npts,.0);
			tensor x(Nc,m);
	
			if(strcmp(InitCond,"Initialize"))
				ReLoadIC(x);
			else {
				LoadIC(x(0),a0,.0,Mesh,InitialDist);
				LoadIC(x(1),.0,b0,Mesh,InitialDist);
			}

			while(t <= PERIOD)
				NewtonSolver(Mesh,Diff,Adv1,Adv2,x,dt,BDFOrder,BCType,Npts);

			for(int i = 0; i < Nc; i++)
				x(i)[0].Write(OutputFname);
		}

	return 0;
}
