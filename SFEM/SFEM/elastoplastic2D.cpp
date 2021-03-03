#include "elastoplastic2D.h"
#include "vonmises.h"


template <class YC>
elastoplastic2D<YC>::elastoplastic2D(Doub young, Doub nu, Doub sigy, Doub thickness, Doub bodyforce, Int planestress, Int order, MatDoub  HHAT) : shapequad(order, 1)
{
	fyoung = young;
	fnu = nu;
	fbodyforce = bodyforce;
	fplanestress = planestress;
	fthickness = thickness;
	fOrder = order;
	fHHAT = HHAT;
}

template <class YC>
elastoplastic2D<YC>::elastoplastic2D(Doub young, Doub nu, Doub sigy, Doub thickness, Doub bodyforce, Int planestress, Int order) : shapequad(order, 1)
{
	fyoung = young;
	fnu = nu;
	fbodyforce = bodyforce;
	fplanestress = planestress;
	fthickness = thickness;
	fOrder = order;
}

template <class YC>
elastoplastic2D<YC>::elastoplastic2D(elastoplastic2D & copy)
{
}


template <class YC>
elastoplastic2D<YC>::~elastoplastic2D()
{
}

template <class YC>
void elastoplastic2D<YC>::SetMemory(Int nglobalpts, Int systemsize)
{
	fdisplace.assign(systemsize,1, 0.);
	fepspvec.assign(nglobalpts, 0.);
	fepspsolitern.assign(nglobalpts, 0.);
	fglobalcounter = 0;
}

template <class YC>
void elastoplastic2D<YC>::UpdateDisplacement(MatDoub displace)
{
	fdisplace = displace;
}
template <class YC>
void elastoplastic2D<YC>::UpdatePlasticStrain()
{
	fepspsolitern = fepspvec;
	fglobalcounter = 0;
}

template <class YC>
void elastoplastic2D<YC>::Contribute(MatDoub &ek, MatDoub &ef, Doub xi, Doub eta, Doub w, MatDoub elcoords,MatDoub eldisplace)
{
	MatDoub psis, GradPsi, elcoordst, xycoords, Jac, InvJac(2, 2), GradPhi, B, BT, N, NT, psist, C, BC, BCS, stress(3, 1, 0.), bodyforce(2, 1), temp, CS, KSt;
	bodyforce[0][0] = 0;
	bodyforce[1][0] = -fbodyforce;
	shapes(psis, GradPsi, xi, eta);
	psis.Transpose(psist);
	psist.Mult(elcoords, xycoords);
	MatDoub hhat;
	if (fhhatvel.nrows() != 0)
	{
		psist.Mult(fhhatvel, hhat);
	}
	GradPsi.Mult(elcoords, Jac);
	Int nnodes = psis.nrows();
	Doub DetJ = -Jac[0][1] * Jac[1][0] + Jac[0][0] * Jac[1][1];
	if (DetJ <= 0)
	{

		std::cout << "\n DetJ < 0 " << endl;
		std::cout << "\n xi " << xi << endl;
		std::cout << "\n eta " << eta << endl;
		GradPsi.Print();
		elcoords.Print();
		GradPsi.Mult(elcoords, Jac);
		xycoords.Print();

		psis.Print();

		shapes(psis, GradPsi, xi, eta);
		return;
	}
	InvJac[0][0] = Jac[1][1] / DetJ;   InvJac[0][1] = -Jac[0][1] / DetJ;
	InvJac[1][0] = -Jac[1][0] / DetJ;	InvJac[1][1] = Jac[0][0] / DetJ;
	InvJac.Mult(GradPsi, GradPhi);
	MatDoub gradu;
	GradPhi.Mult(eldisplace, gradu);

	//{ {dudx, dudy}, { dvdx, dvdy }} = gradprevsol;
	Doub ex = gradu[0][0];// dudx;
	Doub ey = gradu[1][1];
	Doub exy = (gradu[0][1] + gradu[1][0]);


	TensorDoub  epst(0.), epsp(0.), projstress(0.), projstrain(0.), epspeint(0.);
	MatDoub Dep;
	Doub  projgamma=0.;
	epst.XX() = ex;epst.YY() = ey;epst.XY() = exy;
	//epsp = epspsoliternGLOBAL[globalcounter];
	epsp = fepspsolitern[fglobalcounter];

	//std::cout << "epst = " << std::endl;
	//epst.Print();
	//std::cout << "epsp = " << std::endl;
	//epsp.Print();
	fYC.closestpointproj(epst,epsp,projstress,projstrain,Dep,projgamma);

	//std::cout << "projstress = " << std::endl;
	//projstress.Print();

	//epspeint = epst - (epse);
	epspeint = epst;
	epspeint = epspeint - projstrain;

	fepspvec[fglobalcounter] = epspeint;
	//epspvecGLOBAL[globalcounter] = epspeint;

	fglobalcounter++;
	assembleBandN(B, N, psis, GradPhi);
	N.Transpose(NT);
	B.Transpose(BT);
	
	//Dep.Print();
	//BT.Print();
	BT.Mult(Dep, BC);
	BC.Mult(B, ek);

	if (fhhatvel.nrows() != 0)
	{
		Doub mult = hhat[0][0];
		BT.Mult(Dep, BCS);
		BCS.Mult(B, KSt);
		ek += KSt;
	}
	else {

	}
	stress[0][0] = projstress.XX();stress[1][0] = projstress.YY();stress[2][0] = projstress.XY();
	ek *= w*DetJ*fthickness;
	//std::cout << " stress "<< std::endl;
	//stress.Print();
	//ek.Print();
	BT.Mult(stress, ef);
	NT.Mult(bodyforce, temp);
	ef -= temp;
	ef *= w*DetJ;
}
template <class YC>
void elastoplastic2D<YC>::CacStiff(MatDoub &ek, MatDoub &ef, const MatDoub  &elcoords,MatDoub eldisplace)
{
	MatDoub ptsweigths, ekt, eft;
	Doub xi, eta, w;
	Int nnodes = elcoords.nrows();
	ek.assign(nnodes * 2, nnodes * 2, 0.);
	ef.assign(nnodes * 2, 1, 0.);
	shapequad shape = shapequad(fOrder, 1);
	shape.pointsandweigths(ptsweigths);
	Int npts = ptsweigths.nrows();

	for (Int ipt = 0;ipt < npts;ipt++)
	{
		xi = ptsweigths[ipt][0];
		eta = ptsweigths[ipt][1];
		w = ptsweigths[ipt][2];
		//std::cout << "integration point :" << " xi = " << xi << "  eta = "<<  eta << std::endl;
		Contribute(ekt, eft, xi, eta, w, elcoords, eldisplace);
		ek += ekt;
		ef += eft;
	}
}
template <class YC>
void elastoplastic2D<YC>::Assemble(MatDoub &KG, MatDoub &FG, const vector<vector< vector<Doub > > > &allcoords, const MatDoub &meshnodes, const MatInt meshtopology)
{
	MatDoub ek, ef, elcoords, eltopology;
	GetElCoords(allcoords, 0, elcoords);
	Int nnodes = meshnodes.nrows();
	Int rows = elcoords.nrows();
	Int sz = 2 * nnodes;
	Int cols = rows;
	KG.assign(sz, sz, 0.);
	FG.assign(sz, 1, 0.);
	Int nels = allcoords.size();
	
	//uglob = Table[
		//Table[{displacement[[2 topol[[k, j]] - 1]],
		//	displacement[[2 topol[[k, j]]]]}, { j, 1,
			//Length[topol[[k]]] }], { k, 1, nels }];
	NRmatrix<NRvector<Doub>> uglob;
	uglob.resize(nels, nnodes);
	for (int i = 0;i < nels;i++)
	{
		for (int j = 0;j < nnodes;j++) {
			uglob[i][j].assign(2, 0.);
		}
	}
	
	for (Int iel = 0; iel < nels;iel++)
	{
		for (Int node = 0;node < rows;node++)
		{
			uglob[iel][node][0] = fdisplace[2* meshtopology[iel][node]  ][0];
			uglob[iel][node][1]= fdisplace[2* meshtopology[iel][node]+1][0];
		}
	}

	//uglob.Print2();

	Int fu = 0;

	for (Int iel = 0;iel < nels;iel++)
	{
		if (fHHAT.nrows() != 0)
		{
			fhhatvel.assign(rows, 1, 0.);
			for (Int inode = 0;inode < rows;inode++)
			{
				fhhatvel[inode][0] = fHHAT[meshtopology[iel][inode]][0];
			}
		}
		MatDoub elementdisplace(elcoords.nrows(), 2,0.);
		for (Int i = 0;i < elcoords.nrows(); i++)for (Int j = 0;j < 2;j++)elementdisplace[i][j] = uglob[iel][i][j];
		GetElCoords(allcoords, iel, elcoords );
		CacStiff(ek, ef, elcoords, elementdisplace);
		for (Int irow = 0;irow < rows;irow++)
		{
			Int rowglob = meshtopology[iel][irow];
			for (Int icol = 0;icol < cols;icol++)
			{
				Int colglob = meshtopology[iel][icol];
				KG[2 * rowglob + fu][2 * colglob + fu] += ek[2 * irow + fu][2 * icol + fu];
				KG[2 * rowglob + fu][2 * colglob + 1 + fu] += ek[2 * irow + fu][2 * icol + 1 + fu];
				KG[2 * rowglob + 1 + fu][2 * colglob + fu] += ek[2 * irow + 1 + fu][2 * icol + fu];
				KG[2 * rowglob + 1 + fu][2 * colglob + 1 + fu] += ek[2 * irow + 1 + fu][2 * icol + 1 + fu];

			}
			FG[2 * rowglob + fu][0] += ef[2 * irow + fu][0];
			FG[2 * rowglob + 1 + fu][0] += ef[2 * irow + 1 + fu][0];
		}

	}
	fglobalcounter = 0;
}
template <class YC>
void elastoplastic2D<YC>::assembleBandN(MatDoub &B, MatDoub &N, const MatDoub &psis, const MatDoub &GradPhi) 
{
	B.assign(3, psis.nrows() * 2, 0.);
	N.assign(2, psis.nrows() * 2, 0.);

	Int j = 0, k = 0;
	for (Int i = 0;i < psis.nrows() * 2;i++)
	{
		if (i % 2 == 0)
		{
			B[0][i] = GradPhi[0][j];
			B[1][i] = 0;
			B[2][i] = GradPhi[1][j];

			N[0][i] = psis[j][0];
			N[1][i] = 0;
			j++;
		}
		else {
			B[0][i] = 0;
			B[1][i] = GradPhi[1][k];
			B[2][i] = GradPhi[0][k];

			N[0][i] = 0;
			N[1][i] = psis[k][0];

			k++;
		}

	}
}
template <class YC>
void elastoplastic2D<YC>::assembleConstitutiveMatrix(MatDoub &C, Doub mult)
{
	Doub nusqr = fnu*fnu;
	Doub young = mult*fyoung, nu = fnu;
	C.assign(3, 3, 0.);
	C[0][0] = young / (1 - nusqr);   C[0][1] = nu*young / (1 - nusqr);C[0][2] = 0.;
	C[1][0] = nu*young / (1 - nusqr);C[1][1] = young / (1 - nusqr);   C[1][2] = 0.;
	C[2][0] = 0.;                    C[2][1] = 0.;                    C[2][2] = young / (2 * (1 + nu));
}
template <class YC>
void elastoplastic2D<YC>::GetElCoords(vector<vector< vector<Doub > > > allcoords, Int el, MatDoub & elcoords) 
{
	elcoords.assign(allcoords[el].size(), 2, 0.);
	for (Int j = 0; j < allcoords[el].size(); j++)
	{
		Doub x = allcoords[el][j][0];
		Doub y = allcoords[el][j][1];
		elcoords[j][0] = x;
		elcoords[j][1] = y;
	}
}
template <class YC>
void elastoplastic2D<YC>::DirichletBC(MatDoub &KG, MatDoub & FG, vector<int> ids, Int  dir, Int val)
{
	Int nodes = ids.size();
	Int sz = KG.nrows();
	Int fu = 0, cu = 0;
	for (Int i = 0;i < nodes;i++)
	{
		Int pso = ids[i];
		if (dir == 0)
		{
			for (Int j = 0; j < sz;j++)
			{
				KG[2 * pso][j] = 0;
				KG[j][2 * pso] = 0;
			}
			KG[2 * pso][2 * pso] = 1;
			FG[2 * pso][0] = val;
		}
		else
		{
			for (Int j = 0; j < sz;j++)
			{
				KG[2 * pso + 1][j] = 0;
				KG[j][2 * pso + 1] = 0;
			}
			KG[2 * pso + 1][2 * pso + 1] = 1;
			FG[2 * pso + 1][0] = val;
		}


	}
}
template <class YC>
void elastoplastic2D<YC>::ContributeLineNewan(MatDoub &KG, MatDoub & FG, vector<int> ids, Int  dir, Int val)
{
	MatDoub psis, gradpsis, ptsws;

	pointsandweigths1D(ptsws);
	Int npts = ptsws.nrows();

	Doub xi, w;
	MatDoub DetJ;

	for (Int ipt = 0;ipt < npts;ipt++)
	{
		xi = ptsws[ipt][0];
		w = ptsws[ipt][2];
		shapes1D(psis, gradpsis, xi);
		psis.Mult(gradpsis, DetJ);

	}

}

template <class YC>
void elastoplastic2D<YC>::ContributeCurvedLine(MatDoub &KG, MatDoub &FG, MatDoub meshnodes, MatInt linetopology,Doub force)
{
	Int sz = 2*meshnodes.nrows();
	FG.assign(sz, 1, 0.);
	//std::cout << "sz = " << sz << std::endl;
	MatDoub psis, gradpsis, ptsws, DetJ, psist;

	pointsandweigths1D(ptsws);

	Int npts = 1000,els=linetopology.nrows(),nodes= forder+1;
	Doub xi, w;
	for (Int iel = 0;iel < els;iel++)
	{
		MatDoub xy(1, 2, 0.),elcoords(nodes,2,0.), diff(1,2,0.),temp(1, 2, 0.), temp2(1, 2, 0.), xycoords(1, 2, 0.);
		Doub x=0., y=0.;
		for (Int inode = 0;inode < nodes;inode++)
		{
			Doub xmesh = meshnodes[linetopology[iel][inode]][0];
			Doub ymesh = meshnodes[linetopology[iel][inode]][1];
			elcoords[inode][0] = xmesh;
			elcoords[inode][1] = ymesh;
		}
		Doub delta = 2. / npts;
		xi = -1.;
		for (Int ipt = 0;ipt <= npts+1;ipt++)
		{
			shapes1D(psis, gradpsis, xi);
			psis.Transpose(psist);
			temp = xycoords;
			psist.Mult(elcoords, xycoords);
			if (ipt > 0)
			{
				temp2 = xycoords;
				temp2 -= temp;
				diff += temp2;
			}
			xi += delta;
		}
		Doub jac=sqrt(diff[0][0]* diff[0][0] + diff[0][1] * diff[0][1])/2.;
		Int npts = ptsws.nrows();
		MatDoub integral(forder+1, 1, 0.) ;
		for (Int inode = 0;inode < forder + 1;inode++)
		{
			for (Int ipt = 0;ipt < npts;ipt++)
			{
				xi = ptsws[ipt][0];
				w = ptsws[ipt][1];
				shapes1D(psis, gradpsis, xi);
				integral[inode][0] += psis[inode][0] * force*jac*w;
			
			}

			VecDoub normal(2,0.);
			Doub norm = sqrt(elcoords[inode][0] * elcoords[inode][0] + elcoords[inode][1] * elcoords[inode][1]);
			normal[0]= elcoords[inode][0] / norm;
			normal[1] = elcoords[inode][1] / norm;
			//normal.Print();
			//std::cout << "iel " << iel << std::endl;
			//std::cout << "inode " << inode << std::endl;
			//std::cout << " linetopology[iel][inode] * 2 " <<linetopology[iel][inode] * 2 << std::endl;
			FG[ linetopology[iel][inode] * 2 ][0] += integral[inode][0]*normal[0];
			FG[ linetopology[iel][inode] * 2 + 1 ][0] += integral[inode][0]*normal[1];
		}
		
	}

	//FG.Print();
}

template <class YC>
void elastoplastic2D<YC>::SolPt(const vector<vector< vector<Doub > > > &allcoords, const MatInt &meshtopology, const Int &el, const  MatDoub &solG, const Doub &xi, const Doub &eta, MatDoub &xycoords, MatDoub &sol)
{
	MatDoub psis, GradPsi, elcoords, psist, solel;
	GetElCoords(allcoords, el, elcoords);
	shapes(psis, GradPsi, xi, eta);
	Int nodes = psis.nrows();
	Int nstatevars = 1;
	solel.assign(nodes, 1, 0);
	psis.Transpose(psist);
	psist.Mult(elcoords, xycoords);

	for (Int inode = 0; inode < nodes;inode++)
	{
		solel[inode][0] = solG[meshtopology[el][inode]][0];
	}
	psist.Mult(solel, sol);
}

template <class YC>
void elastoplastic2D<YC>::PostProcess(const vector<vector< vector<Doub > > > &allcoords, const MatInt &meshtopology, const MatDoub & nodalsol, vector<vector<double>> &solx, vector<vector<double>> &soly)
{
	MatDoub elcoords, eltopology, psis, gradpsis, xycoords, psist;
	GetElCoords(allcoords, 0, elcoords);
	Int rows = elcoords.nrows();
	Int cols = rows;
	Int nels = allcoords.size();
	Doub refine = 0.1;

	for (Int iel = 0;iel < nels;iel++)
	{
		GetElCoords(allcoords, iel, elcoords);
		for (Doub xi = -1.;xi < 1 - refine;xi += refine)
		{
			vector<double> sol(3);
			for (Doub eta = -1.; eta < 1 - refine;eta += refine)
			{
				Doub approx = 0., approy = 0.;
				shapes(psis, gradpsis, xi, eta);
				psis.Transpose(psist);
				psist.Mult(elcoords, xycoords);
				sol[0] = xycoords[0][0];
				sol[1] = xycoords[0][1];
				for (Int inode = 0;inode < elcoords.nrows();inode++)
				{
					approx += psis[inode][0] * nodalsol[meshtopology[iel][inode] * 2][0];
					approy += psis[inode][0] * nodalsol[meshtopology[iel][inode] * 2 + 1][0];
				}
				sol[2] = approx;
				solx.push_back(sol);
				sol[2] = approy;
				soly.push_back(sol);

			}
		}
	}
}
template <class YC>
void elastoplastic2D<YC>::PostProcess(const vector<vector< vector<Doub > > > &allcoords, const MatInt &meshtopology, const MatDoub & nodalsol, vector<vector<double>> &sol)
{
	MatDoub elcoords, eltopology, psis, gradpsis, xycoords, psist;
	GetElCoords(allcoords, 0, elcoords);
	Int rows = elcoords.nrows();
	Int cols = rows;
	Int nels = allcoords.size();
	Doub refine = 0.05;

	for (Int iel = 0;iel < nels;iel++)
	{
		GetElCoords(allcoords, iel, elcoords);
		for (Doub xi = -1.;xi < 1 - refine;xi += refine)
		{
			vector<double> soli(3);
			for (Doub eta = -1.; eta < 1 - refine;eta += refine)
			{
				Doub approx = 0., approy = 0.;
				shapes(psis, gradpsis, xi, eta);
				psis.Transpose(psist);
				psist.Mult(elcoords, xycoords);
				soli[0] = xycoords[0][0];
				soli[1] = xycoords[0][1];
				for (Int inode = 0;inode < elcoords.nrows();inode++)
				{
					approx += psis[inode][0] * nodalsol[meshtopology[iel][inode]][0];
				}
				soli[2] = approx;
				sol.push_back(soli);

			}
		}
	}
}

template class elastoplastic2D<vonmises>;