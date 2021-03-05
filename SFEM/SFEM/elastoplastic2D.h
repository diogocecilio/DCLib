#pragma once

#include "elastmat2D.h"


template <class YC>
class elastoplastic2D : shapequad{
public:
	elastoplastic2D( Doub thickness, MatDoub bodyforce, Int planestress, Int order);
	elastoplastic2D(Doub young, Doub nu, Doub sigy, Doub thickness, MatDoub bodyforce, Int planestress, Int order, MatDoub  HHAT);
	elastoplastic2D(elastoplastic2D & copy);
	elastoplastic2D();
	~elastoplastic2D();

//	elastmat2D(Doub young, Doub nu, Doub thickness, Doub bodyforce, Int planestress, Int order);
//	elastmat2D(Doub young, Doub nu, Doub thickness, Doub bodyforce, Int planestress, Int order, MatDoub  HHAT);
//	elastmat2D();
//	~elastmat2D();
//
	void Contribute(MatDoub &ek, MatDoub &ef, Doub xi, Doub eta, Doub w, MatDoub elcoords, MatDoub eldisplace);
	void CacStiff(MatDoub &ek, MatDoub &ef, const MatDoub  &elcoords, MatDoub eldisplace);
	void Assemble(MatDoub &KG, MatDoub &FG, const vector<vector< vector<Doub > > > &allcoords, const MatDoub &meshnodes, const MatInt meshtopology);
	void assembleBandN(MatDoub &B, MatDoub &N, const MatDoub &psis, const MatDoub &GradPhi);
	void assembleConstitutiveMatrix(MatDoub &C, Doub mult);
	void GetElCoords(vector<vector< vector<Doub > > > allcoords, Int el, MatDoub & elcoords);
	void DirichletBC(MatDoub &KG, MatDoub & FG, vector<int> ids, Int  dir, Int val);
	void ContributeLineNewan(MatDoub &KG, MatDoub & FG, vector<int> ids, Int  dir, Int val);
	void ContributeCurvedLine(MatDoub &KG, MatDoub &FG, MatDoub meshnodes, MatInt linetopology, Doub force);
	void SolPt(const vector<vector< vector<Doub > > > &allcoords,const MatInt &meshtopology, const Int &el, const  MatDoub &solG, const Doub &xi, const Doub &eta, MatDoub &xycoords, MatDoub &sol);
//
	void PostProcess(const vector<vector< vector<Doub > > > &allcoords,const MatInt &meshtopology, const MatDoub & nodalsol, vector<vector<double>> &solx, vector<vector<double>> &soly);
    void PostProcess(const vector<vector< vector<Doub > > > &allcoords,const MatInt &meshtopology, const MatDoub & nodalsol, vector<vector<double>> &sol);

	//void SetMemory(MatDoub displace, NRvector<TensorDoub> epspvec, NRvector<TensorDoub>  epspsolitern, Int globalcounter);
	void SetMemory(Int ngloblapoints, Int systemsize);
	void UpdateDisplacement(MatDoub displace);
	void UpdatePlasticStrain();
	void UpdateBodyForce(MatDoub newbodyforce);

//
//
private:

	MatDoub fbodyforce;
	Int fplanestress;
	Doub fthickness;
	Int fOrder;



public:
	MatDoub fHHAT;
	MatDoub fhhatvel;
	YC fYC;

	MatDoub fdisplace;
	NRvector<TensorDoub> fepspvec, fepspsolitern;
	Int fglobalcounter;


};

