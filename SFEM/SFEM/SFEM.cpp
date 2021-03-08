// SFEM.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include "shapequad.h"
#include "gridmesh.h"
#include "elastmat2D.h"
#include "KLGalerkinRF.h"
#include "vonmises.h"
#include "druckerprager.h"
#include "elastoplastic2D.h"
#include "ludcmp.h"
#include <math.h>
#include <cmath>
#include "mins_ndim.h"

using namespace std;

Doub beamproblem(Int nx, Int ny, Int order);
void stochasticproblem(MatDoub &HHAT, gridmesh & grid);
void stochasticproblemtalude(MatDoub &HHAT);
void MonteCarlo(MatDoub &HHAT, gridmesh & grid);
void OutPutFile(MatDoub & postdata, std::ofstream &file);
void OutPutFile1var(MatDoub & postdata, std::ofstream &file);
void OutPutFile4var(MatDoub & postdata, std::ofstream &file);
void convergencemeshstudy();
void OutPutPost(vector<vector<double>> & postdata, std::ofstream &file);
void OutPutPost(MatDoub & postdata, std::ofstream &file);
void OutPutPost(MatInt & postdata, std::ofstream &file);
void ReadMesh(vector<vector< vector<Doub > > > &allcoords, MatDoub & meshcoords, MatInt & meshtopology, string filenameel, string filenamecoord);
vector<Int> vecstr_to_vecint(vector<string> vs);
vector<Doub> vecstr_to_vecdoub(vector<string> vs);
void IterativeProcess();
vector<vector<int>>  LineTopology(vector<int> ids, Int order);
void ToMatInt(vector<vector<int>> in, MatInt & out);
void IterativeProcessPressure();
void IterativeSlopeStability();
template <class T>
vector<T> vecstr_to_vec(vector<string> vs);


bool flag = false;

vector<vector<int>>  LineTopology(vector<int> ids, Int order)
{
	Int k = 0;

	vector<vector<int>> vg;
	for (int j = 0;j < ids.size() / order;j++)
	{
		vector<int> v;
		for (int i = 0;i < order + 1;i++)
		{
			v.push_back(ids[i + k]);
		}
		vg.push_back(v);
		k += order;
	}
	return vg;
}

void ToMatInt(vector<vector<int>> in, MatInt & out)
{
	Int  rows = in.size();
	Int cols = in[0].size();
	out.assign(rows, cols, 0.);
	for (Int i = 0;i < rows;i++)for (Int j = 0;j < cols;j++)out[i][j] = in[i][j];
}
//LineTopology[ids_, order_] : = Block[{k, vg, v, j, i},
//k = 0;
//vg = {};
//For[j = 1, j < Length[ids] / (order), j++,
//	v = {};
//For[i = 1, i <= order + 1, i++,
//AppendTo[v, ids[[i + k]]];
//];
//AppendTo[vg, v];
//k += order;
//];
//vg
//];
//

//VecDoub p = ...;
//Func func;

Doub func(VecDoub_I &x)
{
	return pow(x[0], 4) - 3 * pow(x[0], 2) - x[0];
}

Doub func2(VecDoub_I &x)
{
	return sin(tan(x[0]) / 2.);
}

Doub func2d(VecDoub_I &x)
{
	return  (x[0]*x[0] + x[1]*x[1] - 16)*(x[0] * x[0] + x[1] * x[1] - 16);
}


//struct Funcd {
//	Doub operator() (VecDoub_I &x)
//	{
//		return (x[0] * x[0] + x[1] * x[1] - 16)*(x[0] * x[0] + x[1] * x[1] - 16);
//	}
//	void df(VecDoub_I &x, VecDoub_O &deriv) 
//	{
//		deriv[0] = 4 * x[0] * (-16 + pow(x[0], 2) + pow(x[1], 2));
//		deriv[1] = 4 * x[1] * (-16 + pow(x[0], 2) + pow(x[1], 2));
//	}
//};

struct Funcd2 {

	Doub young = 20000.;
	Doub nu = 0.49;
	Doub c = 50;
	Doub phi = M_PI / 9.;
	Doub sig1 = 105.627, sig2 = 78.0876, sig3 = 53.7359;
	Doub tanphi = (3.*tan(phi)) / sqrt(9. + 12.* tan(phi) *tan(phi));
	Doub apex = c * 1. / tan(phi);
	Doub a = (c / (sqrt(3)*tan(phi)) - apex);
	Doub b = a*tanphi;
	Doub K = (young) / (3.* (1. - 2.* nu));
	Doub G = young / (2.* (1. + nu));
	Doub beta = atan((sqrt(3)*(-sig2 + sig3)) / (-2 * sig1 + sig2 + sig3));

	Doub operator() (VecDoub_I &x)
	{
		Doub xiint = x[0];
		return((4 * pow(sqrt(3)*sig1 + sqrt(3)*sig2 + sqrt(3)*sig3 - 3 * xiint, 2)) / K + (9 * pow(-2 * sig1 + sig2 + sig3 +
			2 * sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xiint + pow(xiint, 2))) / pow(a, 2))*cos(beta), 2)) / G +
			(3 * pow(-3 * sig2 + 3 * sig3 + 2 * sqrt(3)*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xiint + pow(xiint, 2))) / pow(a, 2))*sin(beta), 2)) / G) / 108.;
	}
	void df(VecDoub_I &x, VecDoub_O &deriv)
	{
		Doub xi = x[0];
		deriv[0] = ((-24 * (sqrt(3)*sig1 + sqrt(3)*sig2 + sqrt(3)*sig3 - 3 * xi)) / K + (18 * pow(b, 2)*(-2 * sqrt(3)*apex + 2 * xi)*cos(beta)*
			(-2 * sig1 + sig2 + sig3 + 2 * sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*cos(beta))) /
			(pow(a, 2)*G*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))) +
			(6 * sqrt(3)*pow(b, 2)*(-2 * sqrt(3)*apex + 2 * xi)*sin(beta)*(-3 * sig2 + 3 * sig3 +
				2 * sqrt(3)*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*sin(beta))) /
			(pow(a, 2)*G*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2)))) / 108.;
	}
};

Doub distfunddp(VecDoub_I &x){

	Doub young = 20000.;
	Doub nu = 0.49;
	Doub c = 50;
	Doub phi = M_PI / 9.;
	Doub sig1 = 105.627, sig2 = 78.0876, sig3 = 53.7359;
	Doub tanphi = (3.*tan(phi)) / sqrt(9. + 12.* tan(phi) *tan(phi));
	Doub apex = c * 1. / tan(phi);

	Doub a = (c / (sqrt(3)*tan(phi)) - apex);
	Doub b = a*tanphi;
	Doub K = (young) / (3.* (1. - 2.* nu));
	Doub G = young / (2.* (1. + nu));
	Doub xiint = x[0];

	Doub beta = atan((sqrt(3)*(-sig2 + sig3)) / (-2 * sig1 + sig2 + sig3));

	Doub func = ((4 * pow(sqrt(3)*sig1 + sqrt(3)*sig2 + sqrt(3)*sig3 - 3 * xiint,2)) / K + (9 * pow(-2 * sig1 + sig2 + sig3 +
	2 * sqrt((pow(b,2)*(-3 * pow(a,2) + 3 * pow(apex,2) - 2 * sqrt(3)*apex*xiint + pow(xiint,2))) / pow(a,2))*cos(beta),2)) / G +
		(3 * pow(-3 * sig2 + 3 * sig3 + 2 * sqrt(3)*sqrt((pow(b,2)*(-3 * pow(a,2) + 3 * pow(apex,2) - 2 * sqrt(3)*apex*xiint + pow(xiint,2))) / pow(a,2))*sin(beta),2)) / G) / 108.;
	return func;
}


struct Func {
	Doub operator()(VecDoub_I &x);
};

Doub Func::operator()(VecDoub_I &x)
{
	Doub young = 20000.;
	Doub nu = 0.49;
	Doub c = 50;
	Doub phi = M_PI / 9.;
	Doub sig1 = 105.627, sig2 = 78.0876, sig3 = 53.7359;
	Doub tanphi = (3.*tan(phi)) / sqrt(9. + 12.* tan(phi) *tan(phi));
	Doub apex = c * 1. / tan(phi);

	Doub a = (c / (sqrt(3)*tan(phi)) - apex);
	Doub b = a*tanphi;
	Doub K = (young) / (3.* (1. - 2.* nu));
	Doub G = young / (2.* (1. + nu));
	Doub xiint = x[0];

	Doub beta = atan((sqrt(3)*(-sig2 + sig3)) / (-2 * sig1 + sig2 + sig3));

	Doub func = ((4 * pow(sqrt(3)*sig1 + sqrt(3)*sig2 + sqrt(3)*sig3 - 3 * xiint, 2)) / K + (9 * pow(-2 * sig1 + sig2 + sig3 +
		2 * sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xiint + pow(xiint, 2))) / pow(a, 2))*cos(beta), 2)) / G +
		(3 * pow(-3 * sig2 + 3 * sig3 + 2 * sqrt(3)*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xiint + pow(xiint, 2))) / pow(a, 2))*sin(beta), 2)) / G) / 108.;
	return func;
}




int main()
{


	IterativeSlopeStability();

	system("PAUSE");
	return 0;

	Powell<Doub(VecDoub_I &)> powell(distfunddp);
	VecDoub initialguess(1,  (105.627+ 78.0876+ 53.7359)/sqrt(3));
	VecDoub p = powell.minimize(initialguess);
	p.Print();


	Func func;
	Powell<Func> powell2(func);
	p = powell2.minimize(initialguess);
	p.Print();


	Funcd2 funcd;
	Frprmn<Funcd2> frprmn(funcd);
	VecDoub pinit(1,1.);
	VecDoub p0 = frprmn.minimize(initialguess);
	p0.Print();



	//T &func, const Doub ftoll = 3.0e-8



	//Powell<Func> powell(func);

	//Powell<Doub(VecDoub_I &)> powell(func2, 0.0001);

	///VecDoub initialguess(1, -2.);

	//VecDoub p = powell.minimize(initialguess);

	//p.Print();

	//Powell<Doub(VecDoub_I &)> powell2(func2d2, 0.0001);

	//VecDoub initialguess2d(2, 0.);
	//initialguess2d[0] = 3.6;
	//initialguess2d[1] = 1.6;

	//VecDoub p2 = powell2.minimize(initialguess2d);

	//p2.Print();



	IterativeProcessPressure();

	system("PAUSE");
	return 0;
	MatDoub  meshcoords, elcoords;
	MatInt meshtopology;
	vector<vector<vector<Doub>>> allcoords;
	string  elsstr = "elements-pressure.txt";
	string nodestr = "nodes-pressure.txt";
	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);

	std::ofstream filemesh1("meshcoords.txt");
	OutPutPost(meshcoords, filemesh1);

	std::ofstream filemesh2("meshtopology.txt");

	OutPutPost(meshtopology, filemesh2);

	cout << " \n number of elements " << allcoords.size() << endl;
	MatDoub bodyforce(2, 1, 0.);
	elastoplastic2D< vonmises > *  material = new elastoplastic2D< vonmises >( 1., bodyforce, 0, 2);
	MatDoub KG, FG;
	MatInt  linetopology;

	Int ndivs = 10000;
	MatDoub path1;
	vector<int>  idpath1, idpath2, idpath3, iddisplace;
	ndivs = 1000;
	Doub delta;
	path1.assign(ndivs+1, 2, 0.);
	Int i = 0;
	delta = (M_PI / 2.) /( Doub(ndivs));
	for (Doub theta = 0;theta < M_PI / 2.; theta+= delta) {
		path1[i][0] = 100. * cos(theta);
		path1[i][1] = 100. * sin(theta);
		i++;
	}

	gridmesh::FindIdsInPath(path1, allcoords, meshtopology, idpath1);
	for (int i = 0;i < idpath1.size();i++)std::cout << " ID  = " <<idpath1[i] << endl;

	vector<vector<int>>  linetopol = LineTopology(idpath1, 2);
	ToMatInt(linetopol, linetopology);
	linetopology.Print();
	Doub force = 0.19209;
	material->ContributeCurvedLine(KG, FG, meshcoords, linetopology, force);

	IterativeProcess();
	system("PAUSE");
	return 0;
	NRmatrix<NRvector<Doub>> teste(3,4);
	for (int i = 0;i < teste.nrows();i++)
	{
		for (int j = 0;j < teste.ncols();j++) {
			teste[i][j].assign(2, 0.);
		}
	}
	//TESTE DE MUDANCA NO CODIGO
	//Matrix3DDoub teste2(3,3,3);
	//mat.assign(3, 3, 4, 0.);
	//teste.Print();
	//NRmatrix3D<Doub> teste3(3, 3, 4);
	//Doub young = 210., nu = 0.3,sigy=0.24;
	//vonmises * vm = new vonmises(young, nu, sigy);
	//TensorDoub  epst, epsp, projstress, projstrain;
	//MatDoub Dep;
	//Doub  projgamma;
	//epst.XX() = 0.004; epst.YY() = 0.001; epst.XY() = 0.002;
	////epsp.Print();
	////epst.Print();

	//vm->closestpointproj(epst, epsp, projstress, projstrain, Dep, projgamma);

	//std::cout << projgamma << std::endl;
	//Dep.Print();

	system("PAUSE");
	return 0;
	//beamproblem(3,2, 1);
	//convergencemeshstudy();
	//return 0;
	std::clock_t start;
	double duration;

	start = std::clock();

	cout << "\n starting simulation " << endl;
	Doub L = 400.;
	Doub h = 40.;
	Int nx = 10;
	Int ny = 5;
	Int order = 2;
	gridmesh grid = gridmesh(L, h, nx, ny, order);
	MatDoub  HHAT;
	//stochasticproblem(HHAT, grid);
	stochasticproblemtalude(HHAT);



	//MonteCarlo(HHAT, grid);

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;

	std::cout << "\n total simualtion time  " << duration << '\n';

	system("PAUSE");
	return 0;
}

void IterativeSlopeStability()
{
	MatDoub  meshcoords, elcoords;
	MatInt meshtopology;
	vector<vector<vector<Doub>>> allcoords;
	//string  elsstr = "els-grossa-very.txt";
	//string nodestr = "nodes-grossa-very.txt";
	string  elsstr = "mesh-talude-els.txt";
	string nodestr = "mesh-talude-nodes.txt";
	//string nodestr = "mesh-talude-nodes-goodmesh.txt";
	//string elsstr = "mesh-talude-els-goodmesh.txt";
	
	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);

	std::ofstream filemesh1("meshcoords.txt");
	OutPutPost(meshcoords, filemesh1);

	std::ofstream filemesh2("meshtopology.txt");

	OutPutPost(meshtopology, filemesh2);

	cout << " \n number of elements " << allcoords.size() << endl;

	Doub thickness = 1.;
	Doub young = 20000.;
	Doub nu = 0.49;
	Doub c = 50;
	Doub phi = 20. * M_PI / 180.;
	Int planestress = 0;
	MatDoub bodyforce(2, 1, 0.),newbodyforce;
	bodyforce[1][0] = -20.;

	Int ndivs = 1000;
	MatDoub pathbottom, pathleft, pathright,pathdisplace;
	vector<int>  idsbottom, idsleft, idsright, iddisplace;
	VecDoub a(2), b(2);
	a[0] = 0.;a[1] = 0.;
	b[0] = 75.;b[1] = 0;
	gridmesh::Line(a, b, ndivs, pathbottom);
	gridmesh::FindIdsInPath(pathbottom, allcoords, meshtopology, idsbottom);
	//cout << "IDS BOTTOM " << endl;
	//for (Int i = 0;i < idsbottom.size();i++)cout << idsbottom[i] << endl;

	a[0] = 0.;a[1] = 0.;
	b[0] = 0.;b[1] = 40.;
	gridmesh::Line(a, b, ndivs, pathleft);
	gridmesh::FindIdsInPath(pathleft, allcoords, meshtopology, idsleft);
	//cout << "IDS idsleft " << endl;
	//for (Int i = 0;i < idsleft.size();i++)cout << idsleft[i] << endl;
	a[0] = 75.;a[1] = 0.;
	b[0] = 75.;b[1] = 30;
	gridmesh::Line(a, b, ndivs, pathright);
	gridmesh::FindIdsInPath(pathright, allcoords, meshtopology, idsright);
	//cout << "IDS idsright " << endl;
	//for (Int i = 0;i < idsright.size();i++)cout << idsright[i] << endl;
	a[0] = 34.99;a[1] = 39.99;
	b[0] = 35.;b[1] = 40.;
	gridmesh::Line(a, b, ndivs, pathdisplace);
	gridmesh::FindIdsInPath(pathdisplace, allcoords, meshtopology, iddisplace);
	cout << "IDS iddisplace " << endl;
	for (Int i = 0;i < iddisplace.size();i++)cout << iddisplace[i] << endl;
	Int sz = 2 * meshcoords.nrows();
	MatDoub KG(sz, sz, 0.), FBODY(sz, 1, 0.),FINT(sz, 1, 0.), ptsweigths;

	int order = 2;
	shapequad shape = shapequad(order, 1);
	shape.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();
	Int nglobalpts = meshtopology.nrows()* npts;

	MatDoub displace;
	displace.assign(sz, 1, 0.);

	elastoplastic2D< druckerprager > *  material = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order);
	material->fYC.setup(young, nu, c,phi);
	material->SetMemory(nglobalpts, sz);
	material->UpdateBodyForce(bodyforce);
	vector<double>  fac = {0.1,1., 2., 2.135, 2.27, 2.405, 2.54, 2.675, 2.81, 2.945, 3.08, 3.215, 3.35,
		3.485, 3.62, 3.755, 3.9, 4.1, 4.2, 4.3,4.4,4.5,4.525,4.55,4.56,4.57,4.58,4.59,4.6,4.62};
	Doub l = 5, lamb = 1.,lambn,diff=100;
	Int steps = fac.size();
	Int counterout = 0;
	MatDoub solpost(1000, 2, 0.);
	while (counterout < 20 && fabs(diff)> 1.e-2)
	{
		std::cout << "load step = " << counterout << std::endl;
		Int counter = 0, maxcount = 300;
		Doub err1 = 10., err2 = 10., tol = 10e-6;
		MatDoub dws(sz, 1, 0.),dwb(sz,1,0.), dww(sz, 1, 0.),   R;
		lambn = lamb;
		//lamb = 1.;
		Doub dlamb = 0.;
		MatDoub dw(sz, 1, 0.);
		while (counter <  maxcount && err1 > tol)
		{
			
			//newbodyforce = bodyforce;
			//newbodyforce *= lamb;
			//material->UpdateBodyForce(newbodyforce);

			material->Assemble(KG, FINT, FBODY,allcoords, meshcoords, meshtopology);


	
			R = FBODY;
			R *= lamb;
			R -= FINT;
			
			//FBODY *= 1. / lamb;
			Int dir, val;
			dir = 1;
			val = 0;
			material->DirichletBC(KG, R, idsbottom, dir, val);
			dir = 0;
			val = 0;
			material->DirichletBC(KG, R, idsright, dir, val);
			dir = 0;
			val = 0;
			material->DirichletBC(KG, R, idsleft, dir, val);

			dir = 1;
			val = 0;
			material->DirichletBC(KG, FBODY, idsbottom, dir, val);
			dir = 0;
			val = 0;
			material->DirichletBC(KG, FBODY, idsright, dir, val);
			dir = 0;
			val = 0;
			material->DirichletBC(KG, FBODY, idsleft, dir, val);

			MatDoub invKG, sol;
			Cholesky * chol = new Cholesky(KG);
			VecDoub x(sz,0.), b(sz,0.);
			for (int i = 0;i < sz;i++)b[i] = R[i][0];
			chol->solve(b, x);
			if (chol->fail)
			{
				LUdcmp *lu = new LUdcmp(KG);
				lu->solve(b, x);
			}

			for (int i = 0;i < sz;i++)dws[i][0] = x[i];

			x.assign(sz, 0.);
			b.assign(sz, 0.);
			for (int i = 0;i < sz;i++)b[i] = FBODY[i][0];
			chol->solve(b, x);
			if (chol->fail)
			{
				LUdcmp *lu = new LUdcmp(KG);
				lu->solve(b, x);
			}
			for (int i = 0;i < sz;i++)dwb[i][0] = x[i];


			Doub aa=0.;
			for (int i = 0;i < sz;i++)aa += dwb[i][0] * dwb[i][0];

		//	dws.Print();

			Doub bb=0.;
			MatDoub dwcopy = dw;
			dwcopy += dws;
			for (int i = 0;i < sz;i++)bb += dwb[i][0] * dwcopy[i][0];
			bb *= 2;
			//bb = 2 (dw + dws).dwb;

			Doub cc = 0.;
			for (int i = 0;i < sz;i++)cc += dwcopy[i][0] * dwcopy[i][0];
			cc -= l*l;
			//ccc = (dw + dws).(dw + dws) - l ^ 2;
			Doub delta = bb*bb - 4.*aa*cc;
			dlamb = (-bb + sqrt(delta)) / (2. * aa);
		//	cout << "a = " << aa << endl;
		//	cout << "b = " << bb << endl;
		//	cout << "c = " << cc << endl;
		//	cout <<"dlamb = " <<dlamb << endl;
			//dws.Print();
			dww = dwb;
			dww *= dlamb;
			dww += dws;
			dw += dww;
		//	cout << "  dlam = " << dlamb << endl;
			//cout << "  b lamb = " << lamb << endl;
			lamb += dlamb;
		//	cout << " a lamb = " << lamb << endl;

			


			//dww = dws + dlamb dwb;
			//dw += dww;
			//lamb += dlamb;
			//res = Norm[dww];
			//displace += dww;

			//dw.Print();
			displace += dww;
			//displace.Print();
			material->UpdateDisplacement(displace);
			Doub rnorm = 0., normdw = 0., normfg = 0., unorm = 0.;
			rnorm = R.NRmatrixNorm();
			normdw = dww.NRmatrixNorm();
			unorm = displace.NRmatrixNorm();
			FBODY *= lamb;

			err1 = rnorm / FBODY.NRmatrixNorm();
			err2 = normdw / unorm;
			//if (fabs(dlamb) > 0.1)err1 = 10;

			std::cout << " Iteration number = " << counter << " |  |R|/FE = " << err1 <<" |  |R| = " << rnorm <<" | Unrom  = " << unorm << " | lamb  = " << lamb  << " |  dlamb "<< dlamb <<std::endl;
			counter++;
		}
		diff = lamb - lambn;
		cout << "diff = " << diff <<endl ;
		material->UpdatePlasticStrain();
		solpost[counterout][0] = fabs(displace[2 * iddisplace[0]+1][0]);
		solpost[counterout][1] = fabs(lamb*bodyforce[1][0]/20.);
		counterout++;

	}

	std::ofstream file8("loadvsdisplacementlu.txt");
	OutPutFile(solpost, file8);
	vector<vector<double>> solx, soly;
	material->PostProcess(allcoords, meshtopology, displace, solx, soly);
	std::ofstream file("soly.txt");
	OutPutPost(soly, file);

	vector<vector<double>> epsppost;
	material->PostProcessIntegrationPointVar(allcoords, meshtopology, displace, epsppost);
	std::ofstream file2("epsppost.txt");
	OutPutPost(epsppost, file2);

	//return sol[2 * idpath1[0] + 1][0];

}





void IterativeProcessPressure()
{
	MatDoub  meshcoords, elcoords;
	MatInt meshtopology;
	vector<vector<vector<Doub>>> allcoords;
	string  elsstr = "elements-pressure-fino.txt";
	string nodestr = "nodes-pressure-fino.txt";
	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);

	std::ofstream filemesh1("meshcoords.txt");
	OutPutPost(meshcoords, filemesh1);

	std::ofstream filemesh2("meshtopology.txt");

	OutPutPost(meshtopology, filemesh2);

	cout << " \n number of elements " << allcoords.size() << endl;

	Doub young = 210., nu = 0.3, thickness = 1., sigy = 0.24;
	Int planestress = 0;

	MatDoub  bodyforce(2,1,0.);
	Int ndivs = 10000;
	MatDoub path1, path2, path3;
	vector<int>  idpath1, idpath2, iddisplace;
	VecDoub a(2), b(2);
	a[0] = 100.;a[1] = 0.;
	b[0] = 200.;b[1] = 0;
	gridmesh::Line(a, b, ndivs, path1);
	gridmesh::FindIdsInPath(path1, allcoords, meshtopology, idpath1);

	a[0] = 0.;a[1] = 100.;
	b[0] = 0.;b[1] = 200.;
	gridmesh::Line(a, b, ndivs, path2);
	gridmesh::FindIdsInPath(path2, allcoords, meshtopology, idpath2);

	a[0] = 200.;a[1] = 0.;
	b[0] = 200.;b[1] = 0.01;
	gridmesh::Line(a, b, ndivs, path3);
	gridmesh::FindIdsInPath(path3, allcoords, meshtopology, iddisplace);

	Int sz = 2 * meshcoords.nrows();
	MatDoub KG(sz, sz, 0.), FG(sz, 1, 0.), ptsweigths;

	int order = 2;
	shapequad shape = shapequad(order, 1);
	shape.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();
	Int nglobalpts = meshtopology.nrows()* npts;

	MatDoub displace;
	displace.assign(sz, 1, 0.);

	elastoplastic2D< vonmises > *  material = new elastoplastic2D< vonmises >(thickness, bodyforce, planestress, order);
	material->fYC.setup(young, nu, sigy);
	material->SetMemory(nglobalpts, sz);

	Doub finalload = 0.19209;
	Doub fac[] = { 0.1 / finalload, 0.14 / finalload, 0.18 / finalload, 0.19 / finalload, 1. };



	MatInt  linetopology;
	vector<int> idpathcirc;
	MatDoub pathcirc;
	ndivs = 1000;
	Doub delta;
	pathcirc.assign(ndivs + 1, 2, 0.);
	Int i = 0;
	delta = (M_PI / 2.) / (Doub(ndivs));
	for (Doub theta = 0;theta < M_PI / 2.; theta += delta) {
		pathcirc[i][0] = 100. * cos(theta);
		pathcirc[i][1] = 100. * sin(theta);
		i++;
	}

	gridmesh::FindIdsInPath(pathcirc, allcoords, meshtopology, idpathcirc);
	//for (int i = 0;i < idpathcirc.size();i++)std::cout << " ID  = " << idpathcirc[i] << endl;

	vector<vector<int>>  linetopol = LineTopology(idpathcirc, 2);
	ToMatInt(linetopol, linetopology);
	//linetopology.Print();
	material->ContributeCurvedLine(KG, FG, meshcoords, linetopology, finalload);

	//FG.Print();

	Int steps = 5;
	Int counterout = 1;
	MatDoub solpost(1000, 2, 0.);
	for (Int iload = 0; iload < steps; iload++)
	{
		std::cout << "load step = " << iload << std::endl;
		Int counter = 0, maxcount = 30;
		Doub err1 = 10., err2 = 10., tol = 10.e-5;
		MatDoub dw(sz, 1, 0.), res(sz, 1, 0.), FINT,FBODY, R;
		while (counter <  maxcount && err1 > tol)
		{
			MatDoub FGint = FG;
			material->Assemble(KG, FINT, FBODY, allcoords, meshcoords, meshtopology);

			//KG.Print();
			//FINT.Print();

			FGint *= fac[iload];
			FGint -= FINT;
			R = FGint;

			//R.Print();

			Int dir, val;
			dir = 1;
			val = 0;
			material->DirichletBC(KG, R, idpath1, dir, val);
			dir = 0;
			val = 0;
			material->DirichletBC(KG, R, idpath2, dir, val);

			MatDoub invKG, sol;
			Cholesky * chol = new Cholesky(KG);
			chol->inverse(invKG);

			//LUdcmp * LU = new LUdcmp(KG);
			//LU->inverse(invKG);

			invKG.Mult(R, dw);
			//dw.Print();
			displace += dw;
			material->UpdateDisplacement(displace);
			Doub rnorm = 0., normdw = 0., normfg = 0., unorm = 0.;
			for (Int i = 0;i < R.nrows();i++)rnorm += R[i][0] * R[i][0];
			for (Int i = 0;i < dw.nrows();i++)normdw += dw[i][0] * dw[i][0];
			for (Int i = 0;i < FG.nrows();i++)normfg += FG[i][0] * fac[iload] * FG[i][0] * fac[iload];
			for (Int i = 0;i < displace.nrows();i++)unorm += displace[i][0] * displace[i][0];
			rnorm = sqrt(fabs(rnorm));
			normdw = sqrt(fabs(normdw));
			normfg = sqrt(fabs(normfg));
			unorm = sqrt(fabs(unorm));
			err1 = rnorm / normfg;
			err2 = normdw / unorm;
			std::cout << " Iteration number = " << counter << " |  |R|/|FE| = " << err1 << " | deltau/u " << err2 << std::endl;
			counter++;
		}
		material->UpdatePlasticStrain();
		counterout++;
		solpost[iload][0] = fabs(displace[2 * iddisplace[0]][0]);
		solpost[iload][1] = fabs(fac[iload] * finalload);

	}

	std::ofstream file8("loadvsdisplacementlu.txt");
	OutPutFile(solpost, file8);
	vector<vector<double>> solx, soly;
	material->PostProcess(allcoords, meshtopology, displace, solx, soly);
	std::ofstream file("soly.txt");
	OutPutPost(soly, file);
	//return sol[2 * idpath1[0] + 1][0];

}







void IterativeProcess()
{
	Int steps = 23;
	Doub factor = 1. / Doub(steps);
	Doub L = 400;
	Doub h = 40;
	Int nx = 10;
	Int ny = 5;
	Int order = 2;

	gridmesh grid = gridmesh(L, h, nx, ny, order);

	MatDoub  meshcoords;
	MatInt meshtopology;
	vector<vector<vector<Doub>>> allcoords;
	grid.CreateMesh(allcoords, meshcoords, meshtopology);

	std::ofstream filemesh1("meshcoords.txt");
	OutPutPost(meshcoords, filemesh1);

	std::ofstream filemesh2("meshtopology.txt");

	OutPutPost(meshtopology, filemesh2);

	cout << " \n number of elements " << allcoords.size() << endl;

	Doub young = 3000., nu = 0.2, thickness = 1., sigy = 300.;
	Int planestress = 0;


	Int ndivs = 10000;
	MatDoub path1, path2, path3,path4;
	vector<int>  idpath1, idpath2, idpath3,iddisplace;
	VecDoub a(2), b(2);
	a[0] = L;a[1] = 0.;
	b[0] = L;b[1] = h;
	grid.Line(a, b, ndivs, path1);
	grid.FindIdsInPath(path1, allcoords, meshtopology, idpath1);

	a[0] = 0.;a[1] = 0.;
	b[0] = 0.01;b[1] = 0.01;
	grid.Line(a, b, ndivs, path2);
	grid.FindIdsInPath(path2, allcoords, meshtopology, idpath2);


	a[0] = 199;a[1] = 40;
	b[0] = 201;b[1] = 40;
	grid.Line(a, b, ndivs, path3);
	grid.FindIdsInPath(path3, allcoords, meshtopology, idpath3);

	a[0] = L ;a [1] = 0.;
	b[0] = L ;b[1] = 0.01;
	grid.Line(a, b, ndivs, path4);
	grid.FindIdsInPath(path4, allcoords, meshtopology, iddisplace);

	Int sz = 2 * meshcoords.nrows();
	MatDoub elcoords, KG(sz,sz,0.), FG(sz,1,0.), ptsweigths;
	Doub initialload = -100;
	FG[2 * idpath3[0] + 1][0] = initialload;


	shapequad shape = shapequad(order, 1);
	shape.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();
	Int nglobalpts = meshtopology.nrows()* npts;

	MatDoub displace;
	displace.assign(sz, 1, 0.);
	MatDoub  bodyforce(2, 1, 0.);
	elastoplastic2D< vonmises > *  material = new elastoplastic2D< vonmises >(thickness, bodyforce, planestress, order);
	material->fYC.setup(young, nu, sigy);
	material->SetMemory(nglobalpts, sz);

	Doub fac[] = { 3., 3.06555, 3.1311, 3.19665, 3.2622, 3.32775, 3.3933, 3.45885,
		3.5244, 3.58995, 3.6555, 3.72105, 3.7866, 3.85215, 3.9177, 3.98325,
		4.0488, 4.11435, 4.1799, 4.24545, 4.311, 4.37655, 4.4421, 4.50765,
		4.5732, 4.63875, 4.7043, 4.76985, 4.8354, 4.90095, 4.9665, 5.03205,
		5.0976, 5.16315, 5.2287, 5.29425, 5.3598, 5.42535, 5.4909, 5.55645,
		5.622, 5.68755, 5.7531, 5.81865, 5.8842, 5.94975, 6.0153, 6.08085,
		6.1464, 6.21195, 6.2775, 6.34305, 6.4086, 6.47415, 6.5397, 6.60525,
		6.6708, 6.73635, 6.8019, 6.86745, 6.933 };
	steps = 61;
	Int counterout = 1;
	MatDoub solpost(steps, 2, 0.);
	for (Int iload = 0; iload < steps; iload++)
	{
		std::cout << "load step = " << iload << std::endl;
		Int counter = 0, maxcount = 30;
		Doub err1=10., err2=10., tol=10.e-5;
		MatDoub dw(sz, 1, 0.),res(sz,1,0.), FINT,FBODY,R;
		while (counter <  maxcount && err1 > tol)
		{
			MatDoub FGint = FG;
			material->Assemble(KG, FINT, FBODY, allcoords, meshcoords, meshtopology);

			//KG.Print();
			//FINT.Print();

			FGint *= fac[iload];
			FGint -= FINT;
			R = FGint;

			//R.Print();

			Int dir, val;
			dir = 0;
			val = 0;
			material->DirichletBC(KG, R, idpath1, dir, val);
			dir = 1;
			val = 0;
			material->DirichletBC(KG, R, idpath2, dir, val);

			MatDoub invKG, sol;
			//Cholesky * chol = new Cholesky(KG);
			//chol->inverse(invKG);

			LUdcmp * LU = new LUdcmp(KG);
			LU->inverse(invKG);

			invKG.Mult(R, dw);
			//dw.Print();
			displace += dw;
			material->UpdateDisplacement(displace);
			Doub rnorm = 0.,normdw=0., normfg=0.,unorm=0.;
			for (Int i = 0;i < R.nrows();i++)rnorm += R[i][0]* R[i][0];
			for (Int i = 0;i < dw.nrows();i++)normdw += dw[i][0]* dw[i][0];
			for (Int i = 0;i < FG.nrows();i++)normfg += FG[i][0] * fac[iload]* FG[i][0]*fac[iload];
			for (Int i = 0;i < displace.nrows();i++)unorm += displace[i][0]* displace[i][0];
			rnorm = sqrt(fabs(rnorm));
			normdw = sqrt(fabs(normdw));
			normfg = sqrt(fabs(normfg));
			unorm = sqrt(fabs(unorm));
			err1 = rnorm / normfg;
			err2 = normdw / unorm;
			std::cout << " Iteration number = " << counter << " |  |R|/|FE| = " << err1 << " | deltau/u " << err2 << std::endl;
			counter++;
		}
		material->UpdatePlasticStrain();
		counterout++;
		solpost[iload][0] = fabs(displace[2 * iddisplace[0] + 1][0]);
		solpost[iload][1] = fabs(fac[iload]*initialload);

	}
	
	std::ofstream file8("loadvsdisplacementlu.txt");
	OutPutFile(solpost, file8);
	vector<vector<double>> solx, soly;
	material->PostProcess(allcoords, meshtopology, displace, solx, soly);
	std::ofstream file("soly.txt");
	OutPutPost(soly, file);
	//return sol[2 * idpath1[0] + 1][0];

}

void stochasticproblemtalude(MatDoub &HHAT)
{


	MatDoub  meshcoords, elcoords;
	MatInt meshtopology;
	vector<vector<vector<Doub>>> allcoords; 
    string  elsstr= "mesh-talude-els.txt";
    string nodestr = "mesh-talude-nodes.txt";
	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);

	std::ofstream filemesh1("meshcoords.txt");
	OutPutPost(meshcoords, filemesh1);

	std::ofstream filemesh2("meshtopology.txt");

	OutPutPost(meshtopology, filemesh2);

	cout << " \n number of elements " << allcoords.size() << endl;

	Doub young = 20000., nu = 0.49, thickness = 1., bodyforce = 0.;
	Int planestress = 0;

	//subst2 = { young -> 20000, nu -> 0.49, G->young / (2 (1 + nu)),
	//	K -> (young) / (3 (1 - 2 nu)), a -> (c / (Sqrt[3] Tan[phi]) - apex),
	//	b->a tanphi, apex->c Cot[phi] , phi->Pi / 9., c -> 50.,
	//	tanphi -> - (3 Tan[phi]) / Sqrt[9 + 12 Tan[phi] ^ 2] };

	//Doub Lx = 4 * grid.fL / grid.fnx;//(*Correlation length in x direction*)
	//Doub Ly = 4 * grid.fh / grid.fny;//(*Correlation length in y direction*)

	//cout << " \n Correlation Lx " << Lx << endl;
	//cout << " \n Correlation Ly " << Ly << endl;

	Doub Lx = 70;//(*Correlation length in x direction*)ww
	Doub Ly = 70;//(*Correlation length in y direction*)
				  //Lx = 80.;
	Int samples = 1000, expansionorder = 30;
	Doub sig = 0.3;
	Int type = 1, order = 2;
	KLGalerkinRF *objKLGalerkinRF = new KLGalerkinRF(young, nu, thickness, bodyforce, planestress, order, Lx, Ly, sig, type, samples, expansionorder);

	VecComplex  val; MatDoub  vec;
	objKLGalerkinRF->SolveGenEigValProblem(allcoords, meshcoords, meshtopology, val, vec, HHAT);

	//elastmat2D *mat = new elastmat2D();
	//autovetores nas colunas de vec
	for (Int iveccol = 0;iveccol < vec.ncols(); iveccol++)
	{
		vector<vector<double>> eigenfuncx;
		string name = "eigenfunction";
		string ext = ".txt";
		auto s = std::to_string(iveccol);
		name += s;
		name += ext;


		MatDoub eigenfunc(vec.nrows(), 1, 0.);
		for (Int ivecrow = 0;ivecrow < vec.nrows();ivecrow++) {
			eigenfunc[ivecrow][0] = vec[ivecrow][iveccol];
		}
		objKLGalerkinRF->PostProcess(allcoords, meshtopology, eigenfunc, eigenfuncx);
		std::ofstream file(name);
		OutPutPost(eigenfuncx, file);
	}

}



void stochasticproblem(MatDoub &HHAT, gridmesh & grid)
{


	MatDoub  meshcoords, elcoords;
	MatInt meshtopology;
	vector<vector<vector<Doub>>> allcoords;
	grid.CreateMesh(allcoords, meshcoords, meshtopology);

	std::ofstream filemesh1("meshcoords.txt");
	OutPutPost(meshcoords, filemesh1);

	std::ofstream filemesh2("meshtopology.txt");

	OutPutPost(meshtopology, filemesh2);

	cout << " \n number of elements " << allcoords.size() << endl;

	Doub young = 3000., nu = 0.2, thickness = 1., bodyforce = 0.;
	Int planestress = 0;

	//Doub Lx = 4 * grid.fL / grid.fnx;//(*Correlation length in x direction*)
	//Doub Ly = 4 * grid.fh / grid.fny;//(*Correlation length in y direction*)

	//cout << " \n Correlation Lx " << Lx << endl;
	//cout << " \n Correlation Ly " << Ly << endl;

	Doub Lx = 50.;//(*Correlation length in x direction*)ww
	Doub Ly = 20.;//(*Correlation length in y direction*)
	//Lx = 80.;
	Int samples = 1000, expansionorder = 30;
	Doub sig = 0.3;
	Int type =1, order = grid.forder;
	KLGalerkinRF *objKLGalerkinRF = new KLGalerkinRF(young, nu, thickness, bodyforce, planestress, order, Lx, Ly, sig, type, samples, expansionorder);

	VecComplex  val; MatDoub  vec;
	objKLGalerkinRF->SolveGenEigValProblem(allcoords, meshcoords, meshtopology, val, vec, HHAT);

	//elastmat2D *mat = new elastmat2D();
	//autovetores nas colunas de vec
	for (Int iveccol = 0;iveccol < vec.ncols(); iveccol++)
	{
		vector<vector<double>> eigenfuncx;
		string name = "eigenfunction";
		string ext = ".txt";
		auto s = std::to_string(iveccol);
		name += s;
		name += ext;


		MatDoub eigenfunc(vec.nrows(), 1, 0.);
		for (Int ivecrow = 0;ivecrow < vec.nrows();ivecrow++) {
			eigenfunc[ivecrow][0] = vec[ivecrow][iveccol];
		}
		objKLGalerkinRF->PostProcess(allcoords, meshtopology, eigenfunc, eigenfuncx);
		std::ofstream file(name);
		OutPutPost(eigenfuncx, file);
	}

}

void MonteCarlo(MatDoub &HHAT, gridmesh & grid)
{
	std::clock_t start;
	double duration;

	vector<vector< vector<Doub > > > allcoords;
	MatDoub  meshcoords;
	MatInt meshtopology;
	grid.GetData(allcoords, meshcoords, meshtopology);

	Doub young = 3000., nu = 0.2, thickness = 1., bodyforce = 0., L = grid.fL, h = grid.fh;
	Int planestress = 0, postprintfreq = 50;

	Int order = grid.forder;

	Int nsamples = HHAT.ncols();

	Int ndivs = 10000;
	MatDoub path1, path2, path3;
	vector<int>  idpath1, idpath2, idpath3;
	VecDoub a(2), b(2);
	a[0] = L;a[1] = 0.;
	b[0] = L;b[1] = h;
	grid.Line(a, b, ndivs, path1);
	grid.FindIdsInPath(path1, allcoords, meshtopology, idpath1);

	a[0] = 0.;a[1] = 0.;
	b[0] = 0.01;b[1] = 0.01;
	grid.Line(a, b, ndivs, path2);
	grid.FindIdsInPath(path2, allcoords, meshtopology, idpath2);

	a[0] = 199;a[1] = 40;
	b[0] = 201;b[1] = 40;
	grid.Line(a, b, ndivs, path3);
	grid.FindIdsInPath(path3, allcoords, meshtopology, idpath3);

	for (int i = 0; i < idpath3.size();i++)
	{
		std::cout << "\n" << idpath3[i] << std::endl;
	}

	MatDoub hhatinho(HHAT.nrows(), 1, 0.);
	elastmat2D*  material = new elastmat2D(young, nu, thickness, bodyforce, planestress, order, hhatinho);

	for (Int i = 0;i < HHAT.nrows();i++)hhatinho[i][0] = HHAT[i][3];
	vector<vector<double>> hhat, hhat2;
	material->PostProcess(allcoords, meshtopology, hhatinho, hhat);
	std::ofstream file("hhat.txt");
	OutPutPost(hhat, file);

	MatDoub solpost(nsamples, 2, 0.), solpost2(nsamples, 1, 0.);
	Doub sum = 0.;
	for (Int imc = 0;imc < nsamples;imc++)
	{
		MatDoub hhatinho(HHAT.nrows(), 1, 0.);
		Doub mean = 0.;
		Doub var = 0.;
		for (Int i = 0;i < HHAT.nrows();i++) {
			hhatinho[i][0] = HHAT[i][imc];
			mean += hhatinho[i][0];
			mean /= (i+1);
			var += (mean - hhatinho[i][0])*(mean - hhatinho[i][0]);
		}

	

		material->fHHAT = hhatinho;

		MatDoub KG, FG;
		material->Assemble(KG, FG, allcoords, meshcoords, meshtopology);

		if (imc == 0)
		{
			cout << "\n inverting matrix size " << KG.nrows() << endl;
		}

		Int dir1, val1;
		dir1 = 0;
		val1 = 0;
		material->DirichletBC(KG, FG, idpath1, dir1, val1);
		dir1 = 1;
		val1 = 0;
		material->DirichletBC(KG, FG, idpath2, dir1, val1);
		FG[2 * idpath3[0] + 1][0] = -100;
		MatDoub invKG, sol;
		Cholesky * chol = new Cholesky(KG);


		if (imc % postprintfreq == 0)
		{
			start = std::clock();

		}
		chol->inverse(invKG);
		bool fail = chol->fail;
		if (fail == true) {
			start = std::clock();
			LUdcmp * LU = new LUdcmp(KG);
			MatDoub invKG2;
			LU->inverse(invKG2);
			invKG2.Mult(FG, sol);
			duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
			//std::cout << " time to solve with (LU) " << duration << '\n';
		}
		else {
			invKG.Mult(FG, sol);
		}

		if (imc % postprintfreq == 0)
		{
			std::cout << " mean = " << mean << std::endl;
			std::cout << " sdev = " << sqrt(var) << std::endl;

			duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
			std::cout << " time to solve with (Cholesky) " << duration << '\n';

			vector<vector<double>> hhatx;
			string name = "hhat";
			string ext = ".txt";
			auto s = std::to_string(imc);
			name += s;
			name += ext;
			material->PostProcess(allcoords, meshtopology, hhatinho, hhatx);
			std::ofstream file(name);
			OutPutPost(hhatx, file);

		}

		sum += sol[2 * idpath1[0] + 1][0];
		if (imc % postprintfreq == 0)
		{
			cout << " MC sample:  " << imc << ".  Mean vertical displacement in the midle of the beam spam: " << sum / (imc + 1) << endl;

		}

		solpost[imc][0] = imc;
		solpost[imc][1] += sum / (imc + 1);

		solpost2[imc][0] = sol[2 * idpath1[0] + 1][0];
	}

	//solpost.Print();
	std::ofstream file8("saidafina.txt");
	OutPutFile(solpost, file8);

	std::ofstream file2("saidafina2.txt");
	OutPutFile1var(solpost2, file2);

}




Doub beamproblem(Int nx, Int ny, Int order)
{
	Doub L = 400;
	Doub h = 40;

	gridmesh grid = gridmesh(L, h, nx, ny, order);

	MatDoub  meshcoords;
	MatInt meshtopology;
	vector<vector<vector<Doub>>> allcoords;
	grid.CreateMesh(allcoords, meshcoords, meshtopology);

	Doub young = 3000., nu = 0.2, thickness = 1., bodyforce = 0.;
	Int planestress = 0;
	elastmat2D*  material = new elastmat2D(young, nu, thickness, bodyforce, planestress, order);
	MatDoub elcoords, KG, FG;
	material->Assemble(KG, FG, allcoords, meshcoords, meshtopology);

	KG.Print();

	Int ndivs = 10000;
	MatDoub path1, path2, path3;
	vector<int>  idpath1, idpath2, idpath3;
	VecDoub a(2), b(2);
	a[0] = L;a[1] = 0.;
	b[0] = L;b[1] = h;
	grid.Line(a, b, ndivs, path1);
	grid.FindIdsInPath(path1, allcoords, meshtopology, idpath1);

	a[0] = 0.;a[1] = 0.;
	b[0] = 0.01;b[1] = 0.01;
	grid.Line(a, b, ndivs, path2);
	grid.FindIdsInPath(path2, allcoords, meshtopology, idpath2);


	a[0] = 199;a[1] = 40;
	b[0] = 201;b[1] = 40;
	grid.Line(a, b, ndivs, path3);
	grid.FindIdsInPath(path3, allcoords, meshtopology, idpath3);

	std::cout << "\n PATH 1" << std::endl;

	for (int i = 0; i < idpath1.size();i++)
	{
	std::cout << "\n" << idpath1[i] << std::endl;
	}
	cout << "\n PATH 2" << endl;
	for (int i = 0; i < idpath2.size();i++)
	{
	std::cout << "\n" << idpath2[i] << std::endl;
	}
	cout << "\n PATH 3" << endl;
	for (int i = 0; i < idpath3.size();i++)
	{
	std::cout << "\n" << idpath3[i] << std::endl;
	}

	Int dir, val;
	dir = 0;
	val = 0;
	material->DirichletBC(KG, FG, idpath1, dir, val);
	dir = 1;
	val = 0;
	material->DirichletBC(KG, FG, idpath2, dir, val);

	KG.Print();

	if (idpath3.size() == 0)return 0.;
	FG[2 * idpath3[0] + 1][0] = -100;

	FG.Print();
	MatDoub invKG, sol;
	Cholesky * chol = new Cholesky(KG);


	chol->inverse(invKG);
	invKG.Mult(FG, sol);
	sol.Print();
	vector<vector<double>> solx, soly;
	material->PostProcess(allcoords, meshtopology, sol, solx, soly);
	std::ofstream file("soly.txt");
	OutPutPost(soly, file);
	return sol[2 * idpath1[0] + 1][0];


}

void convergencemeshstudy() {
	Int nx = 5, ny = 2, order = 1, nxmax = 30, nymax = 10, ordermax = 2;
	order = 1;
	Doub sol = 0.;
	MatDoub post((nxmax - 5)*(nymax - 2)*ordermax, 4, 0.);
	Int counter = 0;
	sol = beamproblem(6, 2, 1);
	for (nx = 5;nx < nxmax;nx++)
	{
		for (ny = 2;ny < nymax;ny++)
		{
			for (order = 1;order <= ordermax;order++)
			{
				cout << " *sim* " << endl;
				sol = beamproblem(nx, ny, order);
				post[counter][0] = nx;
				post[counter][1] = ny;
				post[counter][2] = order;
				post[counter][3] = sol;
				counter++;
			}

		}
	}
	std::ofstream file2("convergenciadamalha.txt");
	OutPutPost(post, file2);
}


void OutPutFile(MatDoub & postdata, std::ofstream &file)
{

	file.clear();

	for (Int i = 0;i < postdata.nrows(); i++)
	{
		file << postdata[i][0] << " " << postdata[i][1] << endl;
	}

	file.close();
}

void OutPutFile1var(MatDoub & postdata, std::ofstream &file)
{

	file.clear();

	for (Int i = 0;i < postdata.nrows(); i++)
	{
		file << postdata[i][0] << endl;
	}

	file.close();
}

void OutPutFile4var(MatDoub & postdata, std::ofstream &file)
{

	file.clear();

	for (Int i = 0;i < postdata.nrows(); i++)
	{
		file << postdata[i][0] << " " << postdata[i][1] << " " << postdata[i][2] << " " << postdata[i][3] << endl;
	}

	file.close();
}


void OutPutPost(vector<vector<double>> & postdata, std::ofstream &file)
{
	file.clear();
	for (Int i = 0;i < postdata.size(); i++)
	{
		for (Int j = 0;j < postdata[0].size();j++)
		{
			file << postdata[i][j] << " ";
		}
		file << endl;
	}
	file.close();
}

void OutPutPost(MatDoub & postdata, std::ofstream &file)
{
	file.clear();
	for (Int i = 0;i < postdata.nrows(); i++)
	{
		for (Int j = 0;j < postdata.ncols();j++)
		{
			file << postdata[i][j] << " ";
		}
		file << endl;
	}
	file.close();
}
void OutPutPost(MatInt & postdata, std::ofstream &file)
{
	file.clear();
	for (Int i = 0;i < postdata.nrows(); i++)
	{
		for (Int j = 0;j < postdata.ncols();j++)
		{
			file << postdata[i][j] << " ";
		}
		file << endl;
	}
	file.close();
}

void ReadMesh(vector<vector< vector<Doub > > > &allcoords, MatDoub & meshcoords, MatInt & meshtopology, string filenameel, string filenamecoord)
{
	std::vector<vector<Int>> topol;
	string line,temp;

	ifstream myfile(filenameel);
	//
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			vector<string> tokens;
			istringstream iss(line);
			while (iss >> temp)
				tokens.push_back(temp);
			vector<Int> input_int = vecstr_to_vecint(tokens);
			for (int k = 0;k < input_int.size();k++)
			{
				input_int[k] = input_int[k]-1;
			}
			topol.push_back(input_int);
		}
		myfile.close();
	}
	else cout << "Unable to open file";

	meshtopology.CopyFromVector(topol);
	//meshtopology.Print();

	std::vector<vector<Doub>> coords;
	string line2, temp2;
	ifstream myfile2(filenamecoord);
	if (myfile2.is_open())
	{
		while (getline(myfile2, line2))
		{
			vector<string> tokens;
			istringstream iss(line2);
			while (iss >> temp2)
				tokens.push_back(temp2);
			vector<Doub> input_doub = vecstr_to_vecdoub(tokens);

			//vector<Doub> input_doub2(input_doub.size() - 1);
			//for (int k = 1;k < input_doub.size();k++)
			//{
			//	input_doub2[k] = input_doub[k];
			//}

			coords.push_back(input_doub);
		}
		myfile2.close();
	}
	else cout << "Unable to open file";

	meshcoords.CopyFromVector(coords);
	//meshcoords.Print();


	vector<Doub> temp33(2);
	for (Int i = 0; i < meshtopology.nrows();i++)
	{
		vector< vector<Doub> > temp22;
		for (Int j = 0; j < meshtopology.ncols(); j++)
		{
			Int top = meshtopology[i][j];
			temp33[0] = meshcoords[top][0];
			temp33[1] = meshcoords[top][1];
			temp22.push_back(temp33);
		}
		allcoords.push_back(temp22);
	}



}
vector<Int> vecstr_to_vecint(vector<string> vs)
{
	vector<Int> ret;
	for (vector<string>::iterator it = vs.begin()+1;it != vs.end();++it)
	{
		istringstream iss(*it);
		Int temp;
		iss >> temp;
		ret.push_back(temp);
	}
	return ret;
}

vector<Doub> vecstr_to_vecdoub(vector<string> vs)
{
	vector<Doub> ret;
	for (vector<string>::iterator it = vs.begin()+1;it != vs.end()-1;++it)
	{
		istringstream iss(*it);
		Doub temp;
		iss >> temp;
		ret.push_back(temp);
	}
	return ret;
}

template <class T>
vector<T> vecstr_to_vec(vector<string> vs)
{
	vector<T> ret;
	for (vector<string>::iterator it = vs.begin();it != vs.end();++it)
	{
		istringstream iss(*it);
		T temp;
		iss >> temp;
		ret.push_back(temp);
	}
	return ret;
}


