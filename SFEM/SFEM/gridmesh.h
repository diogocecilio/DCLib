#pragma once
#include "nr3.h"

class gridmesh
{
public:
	gridmesh(Doub L, Doub h, Int nx, Int ny, Int order);
	~gridmesh();

	void CreateMesh(vector<vector< vector<Doub > > > &allcoords, MatDoub & meshcoords, MatInt & meshtopology);
	void PrintAllCoords(vector<vector< vector<Doub > > > allcoords, MatInt  meshtopology);
	static void GetElCoords(vector<vector< vector<Doub > > > allcoords, Int el, MatDoub & elcoords);
	static void FindIdsInPath(const MatDoub & path, vector<vector< vector<Doub > > > &allcoords, MatInt & meshtopology, vector<int> & idpath);
	static void Line(VecDoub a, VecDoub b, Int ndivs, MatDoub & path);
	void PrintGMeshVTK(vector<vector< vector<Doub > > >  allcoords, MatInt meshtopology, std::ofstream &file);

	inline void GetData(vector<vector< vector<Doub > > > &allcoords, MatDoub  &meshcoords, MatInt &meshtopology)
	{
		allcoords = fallcoords;
		meshcoords = fmeshcoords;
		meshtopology = fmeshtopology;
	}


private:




public:

	vector<vector< vector<Doub > > > fallcoords;
	MatDoub  fmeshcoords;
	MatInt  fmeshtopology;
	Int forder;
	Doub fL;
	Doub fh;
	Int fnx;
	Int fny;

};

