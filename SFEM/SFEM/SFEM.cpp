// SFEM.cpp : Defines the entry point for the console application.
//

#include <iostream>
//#include <boost/filesystem.hpp>
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/vector_proxy.hpp>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/triangular.hpp>
//#include <boost/numeric/ublas/lu.hpp>
//#include <boost/numeric/ublas/io.hpp>
//#include <iostream>
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/matrix_sparse.hpp>
//#include <boost/numeric/ublas/vector_sparse.hpp>
//#include <boost/numeric/ublas/lu.hpp>
//#include <boost/numeric/ublas/io.hpp>

#include <direct.h>
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
#include <algorithm>
//#include <boost/filesystem.hpp>

//bool globalfail = false;
using namespace std;

MatDoub boostsolve(MatDoub KG, MatDoub FG);
Doub beamproblem(Int nx, Int ny, Int order);
void stochasticproblem(MatDoub &HHAT, gridmesh & grid);
void stochasticproblemtalude(MatDoub &HHAT);
void MonteCarlo(MatDoub &HHAT, gridmesh & grid);
void OutPutFile(MatDoub & postdata, std::ofstream &file);
void OutPutFile1var(MatDoub & postdata, std::ofstream &file);
void OutPutFile4var(MatDoub & postdata, std::ofstream &file);
void convergencemeshstudy();
void OutPutPost(std::vector<std::vector<double>> & postdata, std::ofstream &file);
void OutPutPost(MatDoub & postdata, std::ofstream &file);
void OutPutPost(MatInt & postdata, std::ofstream &file);
void ReadMesh(std::vector<std::vector< std::vector<Doub > > > &allcoords, MatDoub & meshcoords, MatInt & meshtopology, string filenameel, string filenamecoord);
std::vector<Int> vecstr_to_vecint(std::vector<string> vs);
std::vector<Doub> vecstr_to_vecdoub(std::vector<string> vs);
void IterativeProcess();
std::vector<std::vector<int>>  LineTopology(std::vector<int> ids, Int order);
void ToMatInt(std::vector<std::vector<int>> in, MatInt & out);
void IterativeProcessPressure();
std::vector<std::vector<double>>  IterativeSlopeStability(std::vector<std::vector< std::vector<Doub > > > &allcoords, MatDoub & meshcoords, MatInt & meshtopology, MatDoub  hhatinho,elastoplastic2D< druckerprager > *  material);
template <class T>
std::vector<T> vecstr_to_vec(std::vector<string> vs);
std::vector<std::vector<double>>  IterativeSlopeStabilityNew(std::vector<std::vector< std::vector<Doub > > > &allcoords, MatDoub & meshcoords, MatInt & meshtopology, MatDoub  hhatinho, elastoplastic2D< druckerprager > *  material);

void MonteCarloTalude2();
void ReadMatDoub(MatDoub & matdoub, std::string  file);


bool flag = false;

std::vector<std::vector<int>>  LineTopology(std::vector<int> ids, Int order)
{
	Int k = 0;

	std::vector<std::vector<int>> vg;
	for (int j = 0;j < ids.size() / order;j++)
	{
		std::vector<int> v;
		for (int i = 0;i < order + 1;i++)
		{
			v.push_back(ids[i + k]);
		}
		vg.push_back(v);
		k += order;
	}
	return vg;
}

void ToMatInt(std::vector<std::vector<int>> in, MatInt & out)
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


void MonteCarloTalude();


//#include <boost/locale.hpp>
//#include <iostream>
//
//#include <ctime>
//
//int main()
//{
//	using namespace boost::locale;
//	using namespace std;
//	generator gen;
//	locale loc = gen("");
//	// Create system default locale
//
//	locale::global(loc);
//	// Make it system global
//
//	cout.imbue(loc);
//	// Set as default locale for output
//
//	cout << format("Today {1,date} at {1,time} we had run our first localization example") % time(0)
//		<< endl;
//
//	cout << "This is how we show numbers in this locale " << as::number << 103.34 << endl;
//	cout << "This is how we show currency in this locale " << as::currency << 103.34 << endl;
//	cout << "This is typical date in the locale " << as::date << std::time(0) << endl;
//	cout << "This is typical time in the locale " << as::time << std::time(0) << endl;
//	cout << "This is upper case " << to_upper("Hello World!") << endl;
//	cout << "This is lower case " << to_lower("Hello World!") << endl;
//	cout << "This is title case " << to_title("Hello World!") << endl;
//	cout << "This is fold case " << fold_case("Hello World!") << endl;
//
//}


//#include "boost/numeric/ublas/matrix.hpp"
//#include "boost/numeric/ublas/io.hpp"

//using std::vector;
//using boost::numeric::ublas::matrix;

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#ifdef _DEBUG
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
// Replace _NORMAL_BLOCK with _CLIENT_BLOCK if you want the
// allocations to be of _CLIENT_BLOCK type
#else
#define DBG_NEW new
#endif

#include <sys/types.h>
#include <sys/stat.h>

void repruducecriticalcase();

//using namespace boost::ublas;

//bool InvertMatrix(const matrix<double>& input, matrix<double>& inverse);

int main()
{

	repruducecriticalcase();

	system("PAUSE");
	return 0;

	//_CrtMemState sOld;
	//_CrtMemState sNew;
	//_CrtMemState sDiff;
	//_CrtMemCheckpoint(&sOld); //take a snapchot


	//boost::filesystem::create_directories("/tmp/a/b/c");
	////string nodestr = "nodes-mais-fina.txt";
	///string elsstr = "elements-mais-fina.txt";
	////string nodestr = "nodes-fina.txt";
	////string elsstr = "elements-fina.txt";
	//string nodestr = "nodes-extrema-fina.txt";
	//string elsstr = "elements-extrema-fina.txt";
	//string nodestr = "concentrated-nos.txt";
	//string elsstr = "concentrated-els.txt";
	//string nodestr = "nodes-good.txt";
	//string elsstr = "els-good.txt";
	//string nodestr = "nos-110.txt";//ruim
	//string elsstr = "els-110.txt";//ruim
	//string nodestr = "nos-74.txt";
	//string elsstr = "els-74.txt";

	string nodestr = "nos-75.txt";
	string elsstr = "els-75.txt";

	//string nodestr = "nos-94.txt";
	//string elsstr = "els-94.txt";

	MatDoub  meshcoords, elcoords;
	MatInt meshtopology;
	std::vector<std::vector<std::vector<Doub>>> allcoords;
	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);

	std::ofstream filemesh1("meshcoords.txt");
	OutPutPost(meshcoords, filemesh1);

	std::ofstream filemesh2("meshtopology.txt");

	OutPutPost(meshtopology, filemesh2);

	std::clock_t start;
	double duration;
	start = std::clock();

	cout << "\n starting stochastic simulation " << endl;


	Doub thickness = 1.;
	Doub young = 20000.;
	Doub nu = 0.49;
	Doub c = 16.25;
	Doub phi = 20 * M_PI / 180.;
	Int planestress = 0;
	MatDoub bodyforce(2, 1, 0.), newbodyforce;
	bodyforce[1][0] = -20.;
	MatDoub ptsweigths;
	int order = 2;
	shapequad shape = shapequad(order, 1);
	shape.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();
	Int nglobalpts = meshtopology.nrows()* npts;
	Int sz = 2 * meshcoords.nrows();
	MatDoub hhatinho;

	elastoplastic2D< druckerprager > *  material = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
	material->fYC.setup(young, nu, c, phi);
	material->SetMemory(nglobalpts, sz);
	material->UpdateBodyForce(bodyforce);

	std::vector<std::vector<double>>  sol = IterativeSlopeStabilityNew(allcoords, meshcoords, meshtopology, hhatinho, material);//x = desloc y = loadfactor

	delete material;
	////MonteCarloTalude();
	MonteCarloTalude2();
	//_CrtMemCheckpoint(&sNew); //take a snapchot 
	//if (_CrtMemDifference(&sDiff, &sOld, &sNew)) // if there is a difference
	//{
	//	cout << "-----------_CrtMemDumpStatistics ---------"<<endl;
	//	_CrtMemDumpStatistics(&sDiff);
	//	cout << "-----------_CrtMemDumpAllObjectsSince ---------"<<endl;
	//	_CrtMemDumpAllObjectsSince(&sOld);
	//	cout << "-----------_CrtDumpMemoryLeaks ---------" <<endl;
	//	_CrtDumpMemoryLeaks();
	//}

	
	_CrtDumpMemoryLeaks();
	system("PAUSE");
	return 0;

}

void repruducecriticalcase()
{
	MatDoub hhatinho;
	ReadMatDoub(hhatinho, "D:\\DClib\\results-75els-011996\\hhatinho1.145088.txt");
	hhatinho.Print();

	string nodestr = "nos-75.txt";
	string elsstr = "els-75.txt";
	MatDoub  meshcoords, elcoords;
	MatInt meshtopology;
	std::vector<std::vector<std::vector<Doub>>> allcoords;
	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);

	std::ofstream filemesh1("meshcoords.txt");
	OutPutPost(meshcoords, filemesh1);

	std::ofstream filemesh2("meshtopology.txt");

	OutPutPost(meshtopology, filemesh2);

	std::clock_t start;
	double duration;
	start = std::clock();

	cout << "\n starting stochastic simulation " << endl;


	Doub thickness = 1.;
	Doub young = 20000.;
	Doub nu = 0.49;
	Doub c = 16.25;
	Doub phi = 20 * M_PI / 180.;
	Int planestress = 0;
	MatDoub bodyforce(2, 1, 0.), newbodyforce;
	bodyforce[1][0] = -20.;
	MatDoub ptsweigths;
	int order = 2;
	shapequad shape = shapequad(order, 1);
	shape.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();
	Int nglobalpts = meshtopology.nrows()* npts;
	Int sz = 2 * meshcoords.nrows();

	elastoplastic2D< druckerprager > *  material = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
	material->fYC.setup(young, nu, c, phi);
	material->SetMemory(nglobalpts, sz);
	material->UpdateBodyForce(bodyforce);

	std::vector<std::vector<double>>  sol = IterativeSlopeStabilityNew(allcoords, meshcoords, meshtopology, hhatinho, material);//x = desloc y = loadfactor

	int check;
	auto s = std::to_string(rand() % 30 + 1985);
	string namefolder = "D:/DClib/criticalcases";
	namefolder += s;

	char *cstr = new char[namefolder.length() + 1];
	strcpy(cstr, namefolder.c_str());

	check = mkdir(cstr);

	string datafile = namefolder;

	string filename = namefolder;
	std::vector<std::vector<double>> hhatx;
	string name = "/Coesao";
	string ext = ".txt";
	filename += name;
	filename += ext;
	material->PostProcess(0, allcoords, meshtopology, hhatinho, hhatx);
	std::ofstream file(filename);
	OutPutPost(hhatx, file);


	filename = namefolder;
	std::vector<std::vector<double>> hhatx2;
	string namesss = "/Phi";
	string extsss = ".txt";
	filename += namesss;
	filename += ext;
	material->PostProcess(1, allcoords, meshtopology, hhatinho, hhatx2);
	std::ofstream filesss(filename);
	OutPutPost(hhatx2, filesss);



	filename = namefolder;
	std::vector<std::vector<double>> solx, soly;
	material->PostProcess(allcoords, meshtopology, material->fdisplace, solx, soly);
	string name2 = "/soly";
	string ext2 = ".txt";
	filename += name2;
	filename += ext2;
	std::ofstream file2(filename);
	OutPutPost(soly, file2);


	filename = namefolder;
	name2 = "/solx";
	ext2 = ".txt";
	filename += name2;
	filename += ext2;
	std::ofstream file22(filename);
	OutPutPost(solx, file22);


	filename = namefolder;
	name2 = "/hhatinho";
	ext2 = ".txt";
	filename += name2;
	filename += ext2;
	std::ofstream file222(filename);
	OutPutPost(hhatinho, file222);


	filename = namefolder;
	std::vector<std::vector<double>> epsppost;
	material->PostProcessIntegrationPointVar(allcoords, meshtopology, material->fdisplace, epsppost);
	string name3 = "/plasticsqrtj2";
	string ext3 = ".txt";
	filename += name3;
	filename += ext3;
	std::ofstream file3(filename);
	OutPutPost(epsppost, file3);

}

//MatDoub boostsolve(MatDoub KG, MatDoub FG)
//{
//	int sz = KG.nrows();
//	matrix<double> A;
//	boost::numeric::ublas::vector<double> b(sz);
//	KG.Toboost(A);
//	for (int i = 0;i < sz;i++)b(i) = FG[i][0];
//	//compressed_matrix<double, column_major, 0> A(sz,sz,sz*sz);
//	//A = m;
//	permutation_matrix<size_t> pm(A.size1());
//	lu_factorize(A, pm);
//	lu_substitute(A, pm, b);
//	MatDoub x(sz,1, 0.);
//	for (int i = 0;i < sz;i++)x[i][0] = b(i);
//	return x;
//}

std::vector<std::vector<double>>  IterativeSlopeStabilityNew(std::vector<std::vector< std::vector<Doub > > > &allcoords, MatDoub & meshcoords, MatInt & meshtopology, MatDoub  hhatinho, elastoplastic2D< druckerprager > *  material)
{




	Int ndivs = 1000;
	MatDoub pathbottom, pathleft, pathright, pathdisplace;
	std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
	VecDoub a(2), b(2);
	a[0] = 0.;a[1] = 0.;
	b[0] = 30.;b[1] = 0;
	gridmesh::Line(a, b, ndivs, pathbottom);
	gridmesh::FindIdsInPath(pathbottom, allcoords, meshtopology, idsbottom);
	//cout << "IDS BOTTOM " << endl;
	//for (Int i = 0;i < idsbottom.size();i++)cout << idsbottom[i] << endl;

	a[0] = 0.;a[1] = 0.;
	b[0] = 0.;b[1] = 15.;
	gridmesh::Line(a, b, ndivs, pathleft);
	gridmesh::FindIdsInPath(pathleft, allcoords, meshtopology, idsleft);
	//cout << "IDS idsleft " << endl;
	//for (Int i = 0;i < idsleft.size();i++)cout << idsleft[i] << endl;
	a[0] = 30.;a[1] = 0.;
	b[0] = 30.;b[1] = 5;
	gridmesh::Line(a, b, ndivs, pathright);
	gridmesh::FindIdsInPath(pathright, allcoords, meshtopology, idsright);
	//cout << "IDS idsright " << endl;
	//for (Int i = 0;i < idsright.size();i++)cout << idsright[i] << endl;
	a[0] = 9.99;a[1] = 14.99;
	b[0] = 10.;b[1] = 15.;
	gridmesh::Line(a, b, ndivs, pathdisplace);
	gridmesh::FindIdsInPath(pathdisplace, allcoords, meshtopology, iddisplace);
	//cout << "IDS iddisplace " << endl;
	//for (Int i = 0;i < iddisplace.size();i++)cout << iddisplace[i] << endl;
	Int sz = 2 * meshcoords.nrows();
	MatDoub KG(sz, sz, 0.), FBODY(sz, 1, 0.), FINT(sz, 1, 0.), ptsweigths;


	material->ResetPlasticStrain();
	material->ResetDisplacement();

	MatDoub displace;
	displace.assign(sz, 1, 0.);
	//Doub l = 10., lamb = 1., lambn=0, lamb3, diff = 100;
	Doub l =0.5, lamb =0.5, lambn = 0, lamb3, diff = 100,diff2=100;
	Int counterout = 0;
	std::vector<double> solcount(2, 0.);
	std::vector<std::vector<double>> solpost;
	solpost.push_back(solcount);
	//	while (counterout < 15 && fabs(diff)>0.1)
	while (counterout < 10 && fabs(diff2)>0.005)
	{
		std::cout << "load step = " << counterout << std::endl;
		Int counter = 0, maxcount = 30;
		Doub err1 = 10., err2 = 10., tol = 1.e-3;
		MatDoub dws(sz, 1, 0.), dwb(sz, 1, 0.), dww(sz, 1, 0.), R;
		
		//lamb = 1.;
		Doub dlamb = 1., lambn0 = lamb;
		MatDoub dw(sz, 1, 0.);
		//while (counter <  maxcount && err1 > tol && fabs(dlamb)>1.e-3)
		std::cout << "diff = " << diff << std::endl;
		diff = 10;
		//while (counter <  maxcount && err1 > tol && fabs(diff)>0.001)
		while (counter <  maxcount && err1 > tol && fabs(dlamb) > tol)
		{
			lambn = lamb;
			//newbodyforce = bodyforce;
			//newbodyforce *= lamb;
			//material->UpdateBodyForce(newbodyforce);

			material->Assemble(KG, FINT, FBODY, allcoords, meshcoords, meshtopology);



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
			material->DirichletBC(KG, FBODY, idsbottom, dir, val);

			dir = 0;
			val = 0;
			material->DirichletBC(KG, FBODY, idsright, dir, val);
			dir = 0;
			val = 0;
			material->DirichletBC(KG, FBODY, idsleft, dir, val);

			//MatDoub invKG, sol;
			Cholesky * chol = new Cholesky(KG);
			VecDoub x(sz, 0.), b(sz, 0.);
			//MatDoub x(sz,1, 0.), b(sz,1, 0.);
			for (int i = 0;i < sz;i++)b[i] = R[i][0];
			//x = boostsolve(KG, R);
			chol->solve(b, x);
			
			//if (chol->fail)
			//{
			//	LUdcmp *lu = new LUdcmp(KG);
			//	lu->solve(b, x);
		//	}

			for (int i = 0;i < sz;i++)dws[i][0] = x[i];

			x.assign(sz,0.);
			b.assign(sz, 0.);
			for (int i = 0;i < sz;i++)b[i]= FBODY[i][0];
			//x = boostsolve(KG, FBODY);
			chol->solve(b, x);
			delete chol;
			for (int i = 0;i < sz;i++)dwb[i][0] = x[i];
			Doub aa = 0.;
			for (int i = 0;i < sz;i++)aa += dwb[i][0] * dwb[i][0];
			Doub bb = 0.;
			MatDoub dwcopy = dw;
			dwcopy += dws;
			for (int i = 0;i < sz;i++)bb += dwb[i][0] * dwcopy[i][0];
			bb *= 2;
			Doub cc = 0.;
			for (int i = 0;i < sz;i++)cc += dwcopy[i][0] * dwcopy[i][0];
			cc -= l*l;
			Doub delta = bb*bb - 4.*aa*cc;
			dlamb = (-bb + sqrt(delta)) / (2. * aa);
			dww = dwb;
			dww *= dlamb;
			dww += dws;
			dw += dww;
			lamb += dlamb;
			displace += dww;
			material->UpdateDisplacement(displace);
			Doub rnorm = 0., normdw = 0., normfg = 0., unorm = 0.;
			rnorm = R.NRmatrixNorm();
			normdw = dww.NRmatrixNorm();
			unorm = displace.NRmatrixNorm();
			FBODY *= lamb;
			err1 = rnorm / FBODY.NRmatrixNorm();
			err2 = normdw / unorm;
			std::cout << " Iteration number = " << counter << " |  |du|/|u| = " << err2 << " |  |R| = " << " |  |R|/FE = " << err1 << " |  |R| = " << rnorm << " | Unrom  = " << unorm << " | lambn  = " << lambn << " | lamb  = " << lamb << " |  dlamb " << dlamb << std::endl;
			counter++;
			diff = fabs(lamb) - fabs(lambn);
		}
		diff2 = fabs(lambn0 -lamb)  ;
		cout << " diff2 = " << diff2 << endl;
		material->UpdatePlasticStrain();
		solcount[0] = fabs(displace[2 * iddisplace[0] + 1][0]);
		solcount[1] = lamb;
		solpost.push_back(solcount);
		counterout++;
		//std::cout << " Iteration number = " << counter << " |  |R|/FE = " << err1 << " |  |R| = " << rnorm << " | Unrom  = " << unorm << " | lambn  = " << lambn << " | lamb  = " << lamb << " |  dlamb " << dlamb << std::endl;
	}
	//MatDoub out(2, 1, 0.);
	//out[0][0] = solpost[counterout - 1][0];
	//out[1][0] = solpost[counterout - 1][1];

	MatDoub solpost23;
	solpost23.CopyFromVector(solpost);
	string names = "loadvsdisplacementnew";
	string exts = ".txt";
	names += exts;
	std::ofstream file8(names);
	OutPutFile(solpost23, file8);


	std::vector<std::vector<double>> epsppost;
	material->PostProcessIntegrationPointVar(allcoords, meshtopology, material->fdisplace, epsppost);
	string name3 = "epsppostnew";
	string ext3 = ".txt";
	name3 += ext3;
	std::ofstream file3(name3);
	OutPutPost(epsppost, file3);


	return solpost;

}

std::vector<std::vector<double>>  IterativeProcessSlope(std::vector<std::vector< std::vector<Doub > > > &allcoords, MatDoub & meshcoords, MatInt & meshtopology, MatDoub  hhatinho, elastoplastic2D< druckerprager > *  material);



std::vector<std::vector<double>>  IterativeProcessSlope(std::vector<std::vector< std::vector<Doub > > > &allcoords, MatDoub & meshcoords, MatInt & meshtopology, MatDoub  hhatinho, elastoplastic2D< druckerprager > *  material)
{
	Int ndivs = 1000;
	MatDoub pathbottom, pathleft, pathright, pathdisplace;
	std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
	VecDoub a(2), b(2);
	a[0] = 0.;a[1] = 0.;
	b[0] = 30.;b[1] = 0;
	gridmesh::Line(a, b, ndivs, pathbottom);
	gridmesh::FindIdsInPath(pathbottom, allcoords, meshtopology, idsbottom);

	a[0] = 0.;a[1] = 0.;
	b[0] = 0.;b[1] = 15.;
	gridmesh::Line(a, b, ndivs, pathleft);
	gridmesh::FindIdsInPath(pathleft, allcoords, meshtopology, idsleft);
	a[0] = 30.;a[1] = 0.;
	b[0] = 30.;b[1] = 5;
	gridmesh::Line(a, b, ndivs, pathright);
	gridmesh::FindIdsInPath(pathright, allcoords, meshtopology, idsright);
	a[0] = 9.99;a[1] = 14.99;
	b[0] = 10.;b[1] = 15.;
	gridmesh::Line(a, b, ndivs, pathdisplace);
	gridmesh::FindIdsInPath(pathdisplace, allcoords, meshtopology, iddisplace);
	Int sz = 2 * meshcoords.nrows();
	MatDoub KG(sz, sz, 0.), FBODY(sz, 1, 0.), FINT(sz, 1, 0.), ptsweigths;


	material->ResetPlasticStrain();
	material->ResetDisplacement();
	Doub postres = 0;
	MatDoub displace;
	displace.assign(sz, 1, 0.);
	//Doub l = 10., lamb = 1., lambn=0, lamb3, diff = 100;
	Doub l = 0.5, lamb = 0.5, lambn = 0, lamb3, diff = 100, diff2 = 100;
	Int counterout = 0;
	std::vector<double> solcount(7, 0.);
	std::vector<std::vector<double>> solpost;
	solpost.push_back(solcount);
	//	while (counterout < 15 && fabs(diff)>0.1)
    while (counterout <10  && fabs(diff2)>0.005)
	//while (counterout < 2)
	{
		//std::cout << "load step = " << counterout << std::endl;
		Int counter = 0, maxcount = 30;
		Doub err1 = 10., err2 = 10., tol = 1.e-3;
		MatDoub dws(sz, 1, 0.), dwb(sz, 1, 0.), dww(sz, 1, 0.), R;

		//lamb = 1.;
		Doub dlamb = 1., lambn0 = lamb;
		MatDoub dw(sz, 1, 0.);
		//while (counter <  maxcount && err1 > tol && fabs(dlamb)>1.e-3)
		//std::cout << "diff = " << diff << std::endl;
		diff = 10;
		//while (counter <  maxcount && err1 > tol && fabs(diff)>0.001)
		while (counter <  maxcount && err1 > tol && fabs(dlamb) > tol)
		{
			lambn = lamb;
			//newbodyforce = bodyforce;
			//newbodyforce *= lamb;
			//material->UpdateBodyForce(newbodyforce);

			material->Assemble(KG, FINT, FBODY, allcoords, meshcoords, meshtopology);



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
			material->DirichletBC(KG, FBODY, idsbottom, dir, val);

			dir = 0;
			val = 0;
			material->DirichletBC(KG, FBODY, idsright, dir, val);
			dir = 0;
			val = 0;
			material->DirichletBC(KG, FBODY, idsleft, dir, val);

			MatDoub invKG, sol;
			Cholesky * chol = new Cholesky(KG);
			VecDoub x(sz, 0.), b(sz, 0.);
			for (int i = 0;i < sz;i++)b[i] = R[i][0];
			chol->solve(b, x);

			if (chol->fail)
			{
				//globalfail = true;
				throw ("cholesky fail ");
				return solpost;
			}

			for (int i = 0;i < sz;i++)dws[i][0] = x[i];

			x.assign(sz, 0.);
			b.assign(sz, 0.);
			for (int i = 0;i < sz;i++)b[i] = FBODY[i][0];
			chol->solve(b, x);
			if (chol->fail)
			{
				//globalfail = true;
				throw ("cholesky fail ");
				return solpost;
			}
			delete chol;

			for (int i = 0;i < sz;i++)dwb[i][0] = x[i];
			Doub aa = 0.;
			for (int i = 0;i < sz;i++)aa += dwb[i][0] * dwb[i][0];
			Doub bb = 0.;
			MatDoub dwcopy = dw;
			dwcopy += dws;
			for (int i = 0;i < sz;i++)bb += dwb[i][0] * dwcopy[i][0];
			bb *= 2;
			Doub cc = 0.;
			for (int i = 0;i < sz;i++)cc += dwcopy[i][0] * dwcopy[i][0];
			cc -= l*l;
			Doub delta = bb*bb - 4.*aa*cc;
			dlamb = (-bb + sqrt(delta)) / (2. * aa);
			dww = dwb;
			dww *= dlamb;
			dww += dws;
			dw += dww;
			lamb += dlamb;

			displace += dww;
			material->UpdateDisplacement(displace);
			Doub rnorm = 0., normdw = 0., normfg = 0., unorm = 0.;
			rnorm = R.NRmatrixNorm();
			normdw = dww.NRmatrixNorm();
			unorm = displace.NRmatrixNorm();
			FBODY *= lamb;
			err1 = rnorm / FBODY.NRmatrixNorm();
			err2 = normdw / unorm;
			//std::cout << " Iteration number = " << counter << " |  |du|/|u| = " << err2 << " |  |R| = " << " |  |R|/FE = " << err1 << " |  |R| = " << rnorm << " | Unrom  = " << unorm << " | lambn  = " << lambn << " | lamb  = " << lamb << " |  dlamb " << dlamb << std::endl;
			counter++;
			diff = fabs(lamb) - fabs(lambn);
			postres = err1;
		}
		diff2 = fabs(lambn0 - lamb);
		//std::cout << " Iteration number = " << counter;
		//cout << " diff2 = " << diff2 << endl;
		//cout << " diff2 = " << diff2 << endl;
		
		//cout << "diff = " << diff << endl;
		material->UpdatePlasticStrain();
		counterout++;
		solcount[0] = diff;
		solcount[1] = lamb;
		solcount[2] = err1;
		solcount[3] = err2;
		solcount[4] = diff;
		solcount[5] = diff2;
		solcount[6] = counterout;
		solpost.push_back(solcount);
		//std::cout << " Iteration number = " << counter << "   |R|/FE = " << err1 << " | lambn  = " << lambn << " | lamb  = " << lamb << std::endl;
	}
	std::cout << " Iteration number = " << counterout << "   |R|/FE = " << postres << " | lambn  = " << lambn << " | lamb  = " << lamb << std::endl;
	
	return solpost;

}






void MonteCarloTalude2()
{
	


	//string nodestr = "nodes-mais-fina.txt";
	//string elsstr = "elements-mais-fina.txt";
	//string nodestr = "nodes-fina.txt";
	//string elsstr = "elements-fina.txt";
	//string nodestr = "nodes-extrema-fina.txt";
	//string elsstr = "elements-extrema-fina.txt";
	//string nodestr = "concentrated-nos.txt";
	//string elsstr = "concentrated-els.txt";
	//string nodestr = "nodes-good.txt";
	//string elsstr = "els-good.txt";

	//string nodestr = "nos-110.txt";
	//string elsstr = "els-110.txt";

	//string nodestr = "nos-74.txt";
	//string elsstr = "els-74.txt";

	string nodestr = "nos-75.txt";
	string elsstr = "els-75.txt";

	//string nodestr = "nos-94.txt";
	//string elsstr = "els-94.txt";

	MatDoub  meshcoords, elcoords;
	MatInt meshtopology;
	std::vector<std::vector<std::vector<Doub>>> allcoords;
	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);


	
	std::clock_t start;
	double duration;
	start = std::clock();

	cout << "\n starting stochastic simulation " << endl;


	Doub thickness = 1.;
	Doub young = 20000.;
	Doub nu = 0.49;
	Doub c = 16.25;
	Doub phi = 20 * M_PI / 180.;
	Int planestress = 0;
	MatDoub bodyforce(2, 1, 0.), newbodyforce;
	bodyforce[1][0] = -20.;
	MatDoub ptsweigths;
	int order = 2;
	shapequad shape = shapequad(order, 1);
	shape.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();
	Int nglobalpts = meshtopology.nrows()* npts;
	Int sz = 2 * meshcoords.nrows();

	Doub Lx = 25;//(*Correlation length in x direction*)
	Doub Ly = 2.5;//(*Correlation length in y direction*)

	Int samples = 2000, expansionorder = 40;
	Doub sig = 0.1;
	Int type = 1;
	KLGalerkinRF *objKLGalerkinRF = new KLGalerkinRF(young, nu, thickness, bodyforce[1][0], planestress, order, Lx, Ly, sig, type, samples, expansionorder);

	//system("mkdir -p D:\DClib\results");

	//boost::filesystem::create_directories("D:\DClib\results");
	int check;
	//char* dirname = "D:\DClib\results";
	auto s = std::to_string( rand() % 30 + 1985);
	//auto s = std::to_string(8);
	string namefolder = "D:/DClib/results-75els-01";

	namefolder += s;

	char *cstr = new char[namefolder.length() + 1];
	strcpy(cstr, namefolder.c_str());

	check = mkdir(cstr);

	string datafile = namefolder;
	datafile += "/DATA.txt";
	std::ofstream file(datafile);
	file << " Young = " << young << " | nu = " << nu << endl;
	file << " c = " << c << " | phi = " << phi << endl;
	file << " bodyforce = " << bodyforce[1][0]  << endl;
	file << " Mesh  = " << nodestr  << endl;
	file << " samples = " << samples << " | expansion order = " << expansionorder << " | func type = "<< type << endl;
	file <<"Lx = "<<Lx<< " | Ly = "<< Ly << " variance = " << sig  << endl;
	VecComplex  val; MatDoub  vec, HHAT;
	NRmatrix<MatDoub> randomfield;
	objKLGalerkinRF->SolveGenEigValProblem(allcoords, meshcoords, meshtopology, val, vec, randomfield);

	delete objKLGalerkinRF;

	datafile = namefolder;
	datafile += "/vec.txt";
	std::ofstream filevec(datafile);
	OutPutPost(vec, filevec);

	datafile = namefolder;
	datafile += "/val.txt";
	std::ofstream fileval(datafile);
	for (Int j = 0;j < val.size();j++)
	{
		fileval << val[j].real() << endl;
	}

	datafile = namefolder;
	datafile += "/meshcoords.txt";
	std::ofstream filemesh1(datafile);
	OutPutPost(meshcoords, filemesh1);

	datafile = namefolder;
	datafile += "/meshtopology.txt";
	std::ofstream filemesh2(datafile);

	OutPutPost(meshtopology, filemesh2);


	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "\n  simualtion time  " << duration << '\n';

	Int  fail = 0;
	cout << "\n starting Monte Carlo " << endl;
	start = std::clock();
	Int postprintfreq = 50;
	Doub sum = 0.;
	MatDoub solpost(samples, 2, 0.), solpost2(samples, 1, 0.);
	Doub soldatamin = 10.;
	Doub soldatamax =-10;
	std::vector<double> solvec;
	for (Int imc = 0;imc < samples;imc++)
	{
		
		//globalfail = false;
		Int nrandomvars = randomfield.nrows();
		Int nrowss = randomfield[0][0].nrows();
		VecDoub mean(nrandomvars, 0.), var(nrandomvars, 0.);
		MatDoub hhatinho(randomfield[0][0].nrows(), randomfield.nrows(), 0.);
		MatDoub coesmat(randomfield[0][0].nrows(),1,0.),phimat(randomfield[0][0].nrows(),1,0.);
		for (Int ivar = 0;ivar < nrandomvars;ivar++)
		{
			for (Int isample = 0;isample < nrowss;isample++)
			{
				hhatinho[isample][ivar] = randomfield[ivar][0][isample][imc];
				//hhatinho[i][0] = 0.;
				mean[ivar] += hhatinho[isample][ivar];
				mean[ivar] /= (isample + 1);
				var[ivar] += (mean[ivar] - hhatinho[isample][ivar])*(mean[ivar] - hhatinho[isample][ivar]);
			}
		}
		for (int i = 0;i < randomfield[0][0].nrows();i++)coesmat[i][0] = hhatinho[i][0];
		for (int i = 0;i < randomfield[0][0].nrows();i++)phimat[i][0] = hhatinho[i][1];



		elastoplastic2D< druckerprager > *  material = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
		material->fYC.setup(young, nu, c, phi);
		material->SetMemory(nglobalpts, sz);
		material->UpdateBodyForce(bodyforce);

		std::vector<std::vector<double>>  sol = IterativeProcessSlope(allcoords, meshcoords, meshtopology, hhatinho, material);//x = desloc y = loadfactor
		//if (globalfail)
		//{
		//	throw ("fail ");
		//	continue;
	//	}
	//	else {


			MatDoub solpost23;
			solpost23.CopyFromVector(sol);

			Int last = solpost23.nrows() - 1;
			Doub data = solpost23[last][1];
			solvec.push_back(data);

			double min = *min_element(solvec.begin(), solvec.end());
			double max = *max_element(solvec.begin(), solvec.end());
			

			//solpost.Print();

			if (soldatamin>min || soldatamax<max)
			{
				string  filename = namefolder;
				datafile = namefolder;
				datafile += "/information.txt";
				std::ofstream fileinfo(datafile);
				fileinfo << "Monte Carlo Sample = " << imc << std::endl;
				fileinfo << "Safety Factor = " << data << std::endl;
				fileinfo << "r/normr = " << solpost23[last][2] << std::endl;
				fileinfo << "u/normu = " << solpost23[last][3] << std::endl;
				fileinfo << "diff= " << solpost23[last][4] << std::endl;
				fileinfo << "diff2 = " << solpost23[last][5] << std::endl;
				fileinfo << "counterout = " << solpost23[last][6] << std::endl;

			//	solcount[0] = diff;
			//	solcount[1] = lamb;
			//	solcount[2] = err1;
			//	solcount[3] = err2;
			//	solcount[4] = diff;
			///	solcount[5] = diff2;
			//	solcount[6] = counterout;

				filename = namefolder;
				std::vector<std::vector<double>> hhatx;
				string name = "/Coesao";
				string ext = ".txt";
				filename += name;
				auto s = std::to_string(data);
				filename += s;
				filename += ext;
				material->PostProcess(0, allcoords, meshtopology, hhatinho, hhatx);
				std::ofstream file(filename);
				OutPutPost(hhatx, file);


				filename = namefolder;
				std::vector<std::vector<double>> hhatx2;
				string namesss = "/Phi";
				string extsss = ".txt";
				filename += namesss;
				auto sss = std::to_string(data);
				filename += sss;
				filename += ext;
				material->PostProcess(1, allcoords, meshtopology, hhatinho, hhatx2);
				std::ofstream filesss(filename);
				OutPutPost(hhatx2, filesss);

				filename = namefolder;
				string names = "/FxU";
				string exts = ".txt";
				filename += names;
				auto ss = std::to_string(data);
				filename += ss;
				filename += exts;
				std::ofstream file8(filename);
				OutPutFile(solpost23, file8);

				filename = namefolder;
				std::vector<std::vector<double>> solx, soly;
				material->PostProcess(allcoords, meshtopology, material->fdisplace, solx, soly);
				string name2 = "/soly";
				string ext2 = ".txt";
				filename += name2;
				auto s2 = std::to_string(data);
				filename += s2;
				filename += ext2;
				std::ofstream file2(filename);
				OutPutPost(soly, file2);


				filename = namefolder;
				name2 = "/solx";
				ext2 = ".txt";
				filename += name2;
				s2 = std::to_string(data);
				filename += s2;
				filename += ext2;
				std::ofstream file22(filename);
				OutPutPost(solx, file22);


				//filename = namefolder;
				//std::vector<std::vector<double>> hc, hphi;
				//material->PostProcessIntegrationPointVar(allcoords, meshtopology, coesmat, hc);
				//material->PostProcessIntegrationPointVar(allcoords, meshtopology, phimat, hphi);
				//name2 = "/CoesionIntPoint";
				//ext2 = ".txt";
				//filename += name2;
				//s2 = std::to_string(data);
				//filename += s2;
				//filename += ext2;
				//std::ofstream filehc(filename);
				//OutPutPost(hc, filehc);

				//filename = namefolder;
				//name2 = "/PhiIntPoint";
				//ext2 = ".txt";
				//filename += name2;
				//s2 = std::to_string(data);
				//filename += s2;
				//filename += ext2;
				//std::ofstream filephi(filename);
				//OutPutPost(hphi, filephi);


				filename = namefolder;
				name2 = "/hhatinho";
				ext2 = ".txt";
				filename += name2;
				s2 = std::to_string(data);
				filename += s2;
				filename += ext2;
				std::ofstream file222(filename);
				OutPutPost(hhatinho, file222);


				filename = namefolder;
				std::vector<std::vector<double>> epsppost;
				material->PostProcessIntegrationPointVar(allcoords, meshtopology, material->fdisplace, epsppost);
				string name3 = "/plasticsqrtj2";
				string ext3 = ".txt";
				filename += name3;
				auto s3 = std::to_string(data);
				filename += s3;
				filename += ext3;
				std::ofstream file3(filename);
				OutPutPost(epsppost, file3);



				std::cout << "min = " << min << std::endl;
				std::cout << "max = " << max << std::endl;
				std::cout << "soldatamax = " << soldatamax << std::endl;
				std::cout << "soldatamin = " << soldatamin << std::endl;

				if (soldatamin > min)
				{
					soldatamin = min;
				}
				else {
					soldatamax = max;
				}
				
			}

			string filename = namefolder;
			filename += "/montecarlosafetyfactor.txt";
			std::ofstream file23(filename);
			OutPutFile1var(solpost2, file23);

			if (data <= 1.)
			{
				fail++;
			}
			cout << " mc it = " << imc << " | Current safety fator = " << data << endl;
			solpost2[imc][0] = data;
			delete material;
		}
	//}
	string filename = namefolder;
	filename += "/montecarlosafetyfactor.txt";
	std::ofstream file23(filename);
	OutPutFile1var(solpost2, file23);
	file << "failue probability = " << Doub(fail) / Doub(samples) << endl;
	cout << "failue probability = " << Doub(fail) / Doub(samples);
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "\n Monte Carlo simualtion time  " << duration << '\n';
	file << "\n Monte Carlo simualtion time  " << duration << '\n';



}



void MonteCarloTalude()
{
	MatDoub  meshcoords, elcoords;
	MatInt meshtopology;
	std::vector<std::vector<std::vector<Doub>>> allcoords;
	//string  elsstr = "els-grossa-very.txt";
	//string nodestr = "nodes-grossa-very.txt";
	//string  elsstr = "mesh-talude-els.txt";
	//string nodestr = "mesh-talude-nodes.txt";
	//string nodestr = "mesh-talude-nodes-goodmesh.txt";
	//string elsstr = "mesh-talude-els-goodmesh.txt";
//	string nodestr = "nodes-intermediario.txt";
	//string elsstr = "els-intermediario.txt";
	//string nodestr = "nodes-intermediario181.txt";
	//string elsstr = "els-intermediario181.txt";
	//string nodestr = "nodes-talude-tora.txt";
	//string elsstr = "els-talude-tora.txt";
	//string elsstr = "direcional-els.txt";
	//string nodestr = "direcional-nodes.txt";

	string nodestr = "nos-sz-v2.txt";
	string elsstr = "els-sz-v2.txt";

	//string nodestr = "nos-sz-v4.txt";
	//string elsstr = "els-sz-v4.txt";

	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);

	cout << "elements = " << meshtopology.nrows() << endl;

	std::ofstream filemesh1("meshcoords.txt");
	OutPutPost(meshcoords, filemesh1);

	std::ofstream filemesh2("meshtopology.txt");

	OutPutPost(meshtopology, filemesh2);

	std::clock_t start;
	double duration;
	start = std::clock();

	cout << "\n starting stochastic simulation " << endl;


	Doub thickness = 1.;
	Doub young = 20000.;
	Doub nu = 0.49;
	Doub c = 50;
	Doub phi = 20. * M_PI / 180.;
	Int planestress = 0;
	MatDoub bodyforce(2, 1, 0.), newbodyforce;
	bodyforce[1][0] = -20.;
	MatDoub ptsweigths;
	int order = 2;
	shapequad shape = shapequad(order, 1);
	shape.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();
	Int nglobalpts = meshtopology.nrows()* npts;
	Int sz = 2 * meshcoords.nrows();

	Doub Lx = 10;//(*Correlation length in x direction*)
	Doub Ly = 10;//(*Correlation length in y direction*)

	Int samples = 190, expansionorder = 30;
	Doub sig = 0.2;
	Int type = 1;
	KLGalerkinRF *objKLGalerkinRF = new KLGalerkinRF(young, nu, thickness, bodyforce[1][0], planestress, order, Lx, Ly, sig, type, samples, expansionorder);

	VecComplex  val; MatDoub  vec,HHAT;
	NRmatrix<MatDoub> randomfield;
	objKLGalerkinRF->SolveGenEigValProblem(allcoords, meshcoords, meshtopology, val, vec, randomfield);
	

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "\n  simualtion time  " << duration << '\n';


	cout << "\n starting Monte Carlo " << endl;
	start = std::clock();
	Int nsamples = HHAT.ncols();
	Int postprintfreq = 100;
	Doub sum = 0.;
	MatDoub solpost(nsamples, 2, 0.), solpost2(nsamples, 1, 0.);
	for (Int imc = 0;imc < nsamples;imc++)
	{
		MatDoub hhatinho(HHAT.nrows(), randomfield.nrows(), 0.);
		
		Int nrandomvars = randomfield.nrows();
		Int nsamples= randomfield[0][0].nrows();
		VecDoub mean(nrandomvars,0.), var(nrandomvars, 0.);
		for (Int ivar = 0;ivar < nrandomvars;ivar++)
		{
			for (Int isample = 0;isample < nsamples;isample++)
			{
				hhatinho[isample][ivar] = randomfield[ivar][0][isample][imc];
				//hhatinho[i][0] = 0.;
				mean[ivar] += hhatinho[isample][ivar];
				mean[ivar] /= (isample + 1);
				var[ivar] += (mean[ivar] - hhatinho[isample][ivar])*(mean[ivar] - hhatinho[isample][ivar]);
			}
		}

		elastoplastic2D< druckerprager > *  material = new elastoplastic2D< druckerprager >(thickness, bodyforce, planestress, order, hhatinho);
		material->fYC.setup(young, nu, c, phi);
		material->SetMemory(nglobalpts, sz);
		material->UpdateBodyForce(bodyforce);

		std::vector<std::vector<double>>  sol = IterativeSlopeStability(allcoords, meshcoords, meshtopology, hhatinho,material);//x = desloc y = loadfactor

		MatDoub solpost23;
		solpost23.CopyFromVector(sol);

		//solpost.Print();

		if (imc % postprintfreq == 0)
		{
			//std::cout << " mean = " << mean << std::endl;
			//std::cout << " sdev = " << sqrt(var) << std::endl;

			//duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
			//std::cout << " time to solve with (Cholesky) " << duration << '\n';

			std::vector<std::vector<double>> hhatx;
			string name = "hhat3";
			string ext = ".txt";
			auto s = std::to_string(imc);
			name += s;
			name += ext;
			material->PostProcess(0,allcoords, meshtopology, hhatinho, hhatx);
			std::ofstream file(name);
			OutPutPost(hhatx, file);

			string names = "loadvsdisplacementlu3";
			string exts = ".txt";
			auto ss = std::to_string(imc);
			names += ss;
			names += exts;
			std::ofstream file8(names);
			OutPutFile(solpost23, file8);

			//std::vector<std::vector<double>> solx, soly;
			//material->PostProcess(allcoords, meshtopology, material->fdisplace, solx, soly);
			//string name2 = "soly";
			//string ext2 = ".txt";
			//auto s2 = std::to_string(imc);
			//name2 += s2;
			//name2 += ext2;
			//std::ofstream file2(name2);
			//OutPutPost(soly, file2);

			std::vector<std::vector<double>> epsppost;
			material->PostProcessIntegrationPointVar(allcoords, meshtopology, material->fdisplace, epsppost);
			string name3 = "epsppost3";
			string ext3 = ".txt";
			auto s3 = std::to_string(imc);
			name3 += s3;
			name3 += ext3;
			std::ofstream file3(name3);
			OutPutPost(epsppost, file3);

		}
		Int last = solpost23.nrows()-1;
		//cout << "last = " << last;
		Doub data = solpost23[last][1];
		sum += data;
		//if (imc % 100 == 0)
		//{
			cout << " MC sample:  " << imc << " Slope load factor " << sum / (imc + 1) <<"  Nload steps = "<< last << endl;
		//}

		solpost[imc][0] = imc;
		solpost[imc][1] += sum / (imc + 1);
		solpost2[imc][0] = sol[1][0];
		std::ofstream file8("saidafina.txt");
		OutPutFile(solpost, file8);
	}
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "\n Monte Carlo simualtion time  " << duration << '\n';
	//solpost.Print();
	//std::ofstream file8("saidafina.txt");
	//OutPutFile(solpost, file8);

	std::ofstream file2("saidafina2.txt");
	OutPutFile1var(solpost2, file2);

}

std::vector<std::vector<double>>  IterativeSlopeStability(std::vector<std::vector< std::vector<Doub > > > &allcoords, MatDoub & meshcoords, MatInt & meshtopology, MatDoub  hhatinho, elastoplastic2D< druckerprager > *  material)
{




	Int ndivs = 1000;
	MatDoub pathbottom, pathleft, pathright, pathdisplace;
	std::vector<int>  idsbottom, idsleft, idsright, iddisplace;
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
	//cout << "IDS iddisplace " << endl;
	//for (Int i = 0;i < iddisplace.size();i++)cout << iddisplace[i] << endl;
	Int sz = 2 * meshcoords.nrows();
	MatDoub KG(sz, sz, 0.), FBODY(sz, 1, 0.), FINT(sz, 1, 0.), ptsweigths;


	material->ResetPlasticStrain();
	material->ResetDisplacement();

	MatDoub displace;
	displace.assign(sz, 1, 0.);
	//Doub l = 10., lamb = 1., lambn=0, lamb3, diff = 100;
	Doub l =2., lamb = 1., lambn = 0, lamb3, diff = 100;
	Int counterout = 0;
	std::vector<double> solcount(2,0.);
	std::vector<std::vector<double>> solpost;
	
//	while (counterout < 15 && fabs(diff)>0.1)
	while (counterout < 30 && fabs(diff)>1.e-3)
	{
		std::cout << "load step = " << counterout << std::endl;
		Int counter = 0, maxcount = 30;
		Doub err1 = 10., err2 = 10., tol = 10e-8;
		MatDoub dws(sz, 1, 0.), dwb(sz, 1, 0.), dww(sz, 1, 0.), R;
		lambn = lamb;
		//lamb = 1.;
		Doub dlamb = 1.;
		MatDoub dw(sz, 1, 0.);
		//while (counter <  maxcount && err1 > tol && fabs(dlamb)>1.e-3)
		while (counter <  maxcount && err2 > tol && fabs(dlamb)>1.e-20)
		{
			lamb3 = lamb;
			//newbodyforce = bodyforce;
			//newbodyforce *= lamb;
			//material->UpdateBodyForce(newbodyforce);

			material->Assemble(KG, FINT, FBODY, allcoords, meshcoords, meshtopology);



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
			VecDoub x(sz, 0.), b(sz, 0.);
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
			Doub aa = 0.;
			for (int i = 0;i < sz;i++)aa += dwb[i][0] * dwb[i][0];
			Doub bb = 0.;
			MatDoub dwcopy = dw;
			dwcopy += dws;
			for (int i = 0;i < sz;i++)bb += dwb[i][0] * dwcopy[i][0];
			bb *= 2;
			Doub cc = 0.;
			for (int i = 0;i < sz;i++)cc += dwcopy[i][0] * dwcopy[i][0];
			cc -= l*l;
			Doub delta = bb*bb - 4.*aa*cc;
			dlamb = (-bb + sqrt(delta)) / (2. * aa);
			dww = dwb;
			dww *= dlamb;
			dww += dws;
			dw += dww;
			lamb += dlamb;
			displace += dww;
			material->UpdateDisplacement(displace);
			Doub rnorm = 0., normdw = 0., normfg = 0., unorm = 0.;
			rnorm = R.NRmatrixNorm();
			normdw = dww.NRmatrixNorm();
			unorm = displace.NRmatrixNorm();
			FBODY *= lamb;
			err1 = rnorm / FBODY.NRmatrixNorm();
			err2 = normdw / unorm;
			std::cout << " Iteration number = " << counter << " |  |du|/|u| = " << err2 << " |  |R| = " << " |  |R|/FE = " << err1 << " |  |R| = " << rnorm << " | Unrom  = " << unorm << " | lambn  = " << lambn << " | lamb  = " << lamb << " |  dlamb " << dlamb << std::endl;
			counter++;
		}

		diff = fabs(lamb) - fabs(lambn);
		//cout << "diff = " << diff << endl;
		material->UpdatePlasticStrain();
		solcount[0] = fabs(displace[2 * iddisplace[0] + 1][0]);
		solcount[1] = lamb;
		solpost.push_back(solcount);
		counterout++;
		//std::cout << " Iteration number = " << counter << " |  |R|/FE = " << err1 << " |  |R| = " << rnorm << " | Unrom  = " << unorm << " | lambn  = " << lambn << " | lamb  = " << lamb << " |  dlamb " << dlamb << std::endl;
	}
	//MatDoub out(2, 1, 0.);
	//out[0][0] = solpost[counterout - 1][0];
	//out[1][0] = solpost[counterout - 1][1];

	MatDoub solpost23;
	solpost23.CopyFromVector(solpost);
	string names = "loadvsdisplacementluphiiguala1";
	string exts = ".txt";
	names += exts;
	std::ofstream file8(names);
	OutPutFile(solpost23, file8);


	std::vector<std::vector<double>> epsppost;
	material->PostProcessIntegrationPointVar(allcoords, meshtopology, material->fdisplace, epsppost);
	string name3 = "epsppostphiiguala1";
	string ext3 = ".txt";
	name3 += ext3;
	std::ofstream file3(name3);
	OutPutPost(epsppost, file3);


	return solpost;

}




void stochasticproblemtalude(MatDoub &HHAT)
{


	MatDoub  meshcoords, elcoords;
	MatInt meshtopology;
	std::vector<std::vector<std::vector<Doub>>> allcoords;
	//string  elsstr= "mesh-talude-els.txt";
	//string nodestr = "mesh-talude-nodes.txt";
	//string nodestr = "nodes-intermediario.txt";
	//string elsstr = "els-intermediario.txt";
	//string  elsstr = "els-grossa-very.txt";
	//string nodestr = "nodes-grossa-very.txt";

	//string nodestr = "mesh-talude-nodes-goodmesh.txt";
	//string elsstr = "mesh-talude-els-goodmesh.txt";
	//string nodestr = "nodes-intermediario.txt";
	//string elsstr = "els-intermediario.txt";

	//string elsstr = "direcional-els.txt";
	//string nodestr = "direcional-nodes.txt";

	//string nodestr = "nodes-talude-tora.txt";
	//string elsstr = "els-talude-tora.txt";

	string nodestr = "nos-sz.txt";
	string elsstr = "els-sz.txt";


	ReadMesh(allcoords, meshcoords, meshtopology, elsstr, nodestr);

	std::ofstream filemesh1("meshcoords.txt");
	OutPutPost(meshcoords, filemesh1);

	std::ofstream filemesh2("meshtopology.txt");

	OutPutPost(meshtopology, filemesh2);

	cout << " \n number of elements " << allcoords.size() << endl;

	Doub young = 20000., nu = 0.49, thickness = 1., bodyforce = 0.;
	Int planestress = 0;

	Doub Lx = 70;//(*Correlation length in x direction*)ww
	Doub Ly = 70;//(*Correlation length in y direction*)

	Int samples = 50, expansionorder = 30;
	Doub sig = 0.2;
	Int type = 1, order = 2;
	KLGalerkinRF *objKLGalerkinRF = new KLGalerkinRF(young, nu, thickness, bodyforce, planestress, order, Lx, Ly, sig, type, samples, expansionorder);

	VecComplex  val; MatDoub  vec;
	NRmatrix<MatDoub> randomfield;
	objKLGalerkinRF->SolveGenEigValProblem(allcoords, meshcoords, meshtopology, val, vec, randomfield);
	HHAT = randomfield[0][0];

	for (Int iveccol = 0;iveccol < vec.ncols(); iveccol++)
	{
		std::vector<std::vector<double>> eigenfuncx;
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





void IterativeProcessPressure()
{
	MatDoub  meshcoords, elcoords;
	MatInt meshtopology;
	std::vector<std::vector<std::vector<Doub>>> allcoords;
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
	std::vector<int>  idpath1, idpath2, iddisplace;
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
	std::vector<int> idpathcirc;
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

	std::vector<std::vector<int>>  linetopol = LineTopology(idpathcirc, 2);
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
	std::vector<std::vector<double>> solx, soly;
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
	std::vector<std::vector<std::vector<Doub>>> allcoords;
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
	std::vector<int>  idpath1, idpath2, idpath3,iddisplace;
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
	std::vector<std::vector<double>> solx, soly;
	material->PostProcess(allcoords, meshtopology, displace, solx, soly);
	std::ofstream file("soly.txt");
	OutPutPost(soly, file);
	//return sol[2 * idpath1[0] + 1][0];

}





void stochasticproblem(MatDoub &HHAT, gridmesh & grid)
{


	MatDoub  meshcoords, elcoords;
	MatInt meshtopology;
	std::vector<std::vector<std::vector<Doub>>> allcoords;
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
	NRmatrix<MatDoub> randomfield;
	objKLGalerkinRF->SolveGenEigValProblem(allcoords, meshcoords, meshtopology, val, vec, randomfield);
	HHAT = randomfield[0][0];

	//elastmat2D *mat = new elastmat2D();
	//autovetores nas colunas de vec
	for (Int iveccol = 0;iveccol < vec.ncols(); iveccol++)
	{
		std::vector<std::vector<double>> eigenfuncx;
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

	std::vector<std::vector< std::vector<Doub > > > allcoords;
	MatDoub  meshcoords;
	MatInt meshtopology;
	grid.GetData(allcoords, meshcoords, meshtopology);

	Doub young = 3000., nu = 0.2, thickness = 1., bodyforce = 0., L = grid.fL, h = grid.fh;
	Int planestress = 0, postprintfreq = 50;

	Int order = grid.forder;

	Int nsamples = HHAT.ncols();

	Int ndivs = 10000;
	MatDoub path1, path2, path3;
	std::vector<int>  idpath1, idpath2, idpath3;
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
	std::vector<std::vector<double>> hhat, hhat2;
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

			std::vector<std::vector<double>> hhatx;
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
	std::vector<std::vector<std::vector<Doub>>> allcoords;
	grid.CreateMesh(allcoords, meshcoords, meshtopology);

	Doub young = 3000., nu = 0.2, thickness = 1., bodyforce = 0.;
	Int planestress = 0;
	elastmat2D*  material = new elastmat2D(young, nu, thickness, bodyforce, planestress, order);
	MatDoub elcoords, KG, FG;
	material->Assemble(KG, FG, allcoords, meshcoords, meshtopology);

	KG.Print();

	Int ndivs = 10000;
	MatDoub path1, path2, path3;
	std::vector<int>  idpath1, idpath2, idpath3;
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
	std::vector<std::vector<double>> solx, soly;
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


void OutPutPost(std::vector<std::vector<double>> & postdata, std::ofstream &file)
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

void ReadMesh(std::vector<std::vector< std::vector<Doub > > > &allcoords, MatDoub & meshcoords, MatInt & meshtopology, string filenameel, string filenamecoord)
{
	std::vector<std::vector<Int>> topol;
	string line,temp;

	ifstream myfile(filenameel);
	//
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			std::vector<string> tokens;
			istringstream iss(line);
			while (iss >> temp)
				tokens.push_back(temp);
			std::vector<Int> input_int = vecstr_to_vecint(tokens);
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

	std::vector<std::vector<Doub>> coords;
	string line2, temp2;
	ifstream myfile2(filenamecoord);
	if (myfile2.is_open())
	{
		while (getline(myfile2, line2))
		{
			std::vector<string> tokens;
			istringstream iss(line2);
			while (iss >> temp2)
				tokens.push_back(temp2);
			std::vector<Doub> input_doub = vecstr_to_vecdoub(tokens);

			//std::vector<Doub> input_doub2(input_doub.size() - 1);
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


	std::vector<Doub> temp33(2);
	for (Int i = 0; i < meshtopology.nrows();i++)
	{
		std::vector< std::vector<Doub> > temp22;
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
std::vector<Int> vecstr_to_vecint(std::vector<string> vs)
{
	std::vector<Int> ret;
	for (std::vector<string>::iterator it = vs.begin()+1;it != vs.end();++it)
	{
		istringstream iss(*it);
		Int temp;
		iss >> temp;
		ret.push_back(temp);
	}
	return ret;
}

std::vector<Doub> vecstr_to_vecdoub(std::vector<string> vs)
{
	std::vector<Doub> ret;
	for (std::vector<string>::iterator it = vs.begin()+1;it != vs.end()-1;++it)
	{
		istringstream iss(*it);
		Doub temp;
		iss >> temp;
		ret.push_back(temp);
	}
	return ret;
}

std::vector<Doub> vecstr_to_vecdoub2(std::vector<string> vs)
{
	std::vector<Doub> ret;
	for (std::vector<string>::iterator it = vs.begin() ;it != vs.end() ;++it)
	{
		istringstream iss(*it);
		Doub temp;
		iss >> temp;
		ret.push_back(temp);
	}
	return ret;
}

template <class T>
std::vector<T> vecstr_to_vec(std::vector<string> vs)
{
	std::vector<T> ret;
	for (std::vector<string>::iterator it = vs.begin();it != vs.end();++it)
	{
		istringstream iss(*it);
		T temp;
		iss >> temp;
		ret.push_back(temp);
	}
	return ret;
}

void ReadMatDoub(MatDoub & matdoub, std::string  file)
{
	std::vector<std::vector<Doub>> coords;
	string line2, temp2;
	ifstream myfile2(file);
	if (myfile2.is_open())
	{
		while (getline(myfile2, line2))
		{
			std::vector<string> tokens;
			istringstream iss(line2);
			while (iss >> temp2)
				tokens.push_back(temp2);
			std::vector<Doub> input_doub = vecstr_to_vecdoub2(tokens);
			coords.push_back(input_doub);
		}
		myfile2.close();
	}
	else cout << "Unable to open file";
	//for (int i = 0;i < coords.size();i++)
	//{
	//	for (int j = 0;j < coords[0].size();j++)
	//	{
	//		cout << coords[i][j] << endl;
	//	}
	//	cout << endl;
	//}

	matdoub.CopyFromVector(coords);
}


//// REMEMBER to update "lu.hpp" header includes from boost-CVS
//
//
///* Matrix inversion routine.
//Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
//bool InvertMatrix(const matrix<double>& input, matrix<double>& inverse) {
//	using namespace boost::numeric::ublas;
//	typedef permutation_matrix<std::size_t> pmatrix;
//	// create a working copy of the input
//	matrix<double> A(input);
//	// create a permutation matrix for the LU-factorization
//	pmatrix pm(A.size1());
//
//	// perform LU-factorization
//	int res = lu_factorize(A, pm);
//	if (res != 0) return false;
//
//	// create identity matrix of "inverse"
//	inverse.assign(identity_matrix<double>(A.size1()));
//
//	// backsubstitute to get the inverse
//	lu_substitute(A, pm, inverse);
//
//	return true;
//}