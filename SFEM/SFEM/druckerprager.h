#pragma once
#include "nr3.h"
#include <math.h>
#include <cmath>
#include "elastoplasticbase.h"
#include "roots.h"
class druckerprager
{
public:
	druckerprager(Doub young, Doub nu, Doub coesion, Doub frictionangle);
	druckerprager();
	~druckerprager();

	void closestpointproj(TensorDoub epst, TensorDoub epsp, TensorDoub & projstress, TensorDoub & projstrain, MatDoub & Dep, Doub & projgamma);
	Doub yield(Doub xi, Doub rho);
	MatDoub GetElasticMatrix();
	MatDoub GetInverseElasticMatrix();
	MatDoub F1HWCylDruckerPragerSmoothPSMATCH(Doub xisol, Doub rho, Doub betasol);
	MatDoub HW(MatDoub sig);
	MatDoub stressrecosntruction(MatDoub val, MatDoub vec);
	MatDoub dadsig(TensorDoub sigprojvoigt);
	MatDoub P();
	MatDoub avec(TensorDoub sigprojvoigt);

	Doub FindMinimum(MatDoub pt,Doub xitrial,bool flag);
	MatDoub ComputeQ(MatDoub fulltensorproj, MatDoub tempepsemat, TensorDoub & projstress, TensorDoub & projstrain, Doub & projgamma,MatDoub & nvec);
	
	// Doub func2(Doub xi)
	//{
	//	Doub sig1 = pt[0][0], sig2 = pt[0][1], sig3 = pt[0][2];
	//	Doub beta = atan((sqrt(3)*(-sig2 + sig3)) / (-2 * sig1 + sig2 + sig3));

	//	Doub fun = ((4. * pow(sqrt(3.)*sig1 + sqrt(3)*sig2 + sqrt(3.)*sig3 - 3. * xi, 2)) / fK + (9. * pow(-2. * sig1 + sig2 + sig3 +
	//		2 * sqrt((pow(fb, 2)*(-3. * pow(fa, 2) + 3. * pow(fapex, 2) - 2. * sqrt(3.)*fapex*xi + pow(xi, 2))) / pow(fa, 2.))*cos(beta), 2.)) / fG +
	//		(3 * pow(-3. * sig2 + 3. * sig3 + 2 * sqrt(3.)*sqrt((pow(fb, 2)*(-3. * pow(fa, 2) + 3. * pow(fapex, 2) - 2. * sqrt(3.)*fapex*xi + pow(xi, 2))) / pow(fa, 2))*sin(beta), 2)) / fG) / 108.;
	//	return fun;
	//}
	//Doub   distfunddp(MatDoub &pt, VecDoub_I &x);
	inline void setup(Doub young, Doub nu, Doub coesion, Doub frictionangle) {
		fyoung = young;
		fnu = nu;
		fcoesion = coesion;
		fphi = frictionangle;
		fG = young / (2.* (1. + nu));
		fK = young / (3.* (1. - 2.* nu));
		ftanphi = (3.*tan(fphi)) / sqrt(9. + 12.* tan(fphi) *tan(fphi));
		fapex = fcoesion * 1. / tan(fphi);
		fa = (fcoesion / (sqrt(3)*tan(fphi)) - fapex);
		fb = fa*ftanphi;

	}
	 Doub fyoung;
	 Doub fnu;

private:
	Doub fcoesion;
	Doub fphi;
	Doub fK;
	Doub fG;
	Doub ftanphi;
	Doub fapex;
	Doub fa, fb;
	bool fflag = true;
};


struct Funcstruct:druckerprager {

	MatDoub pt;
	Doub fyoung;
	Doub fnu;
	Doub fcoesion;
	Doub fphi;
	Doub fK;
	Doub fG;
	Doub ftanphi;
	Doub fapex;
	Doub fa, fb;



	Doub operator()(VecDoub_I &x)
	{
		Doub sig1 = pt[0][0], sig2 = pt[0][1], sig3 = pt[0][2];
		Doub xi = x[0];
		if (xi > fapex) {
			//xi = fapex;
		}
		Doub beta = atan((sqrt(3)*(-sig2 + sig3)) / (-2 * sig1 + sig2 + sig3));

		Doub func2 = ((4. * pow(sqrt(3.)*sig1 + sqrt(3)*sig2 + sqrt(3.)*sig3 - 3. * xi, 2)) / fK + (9. * pow(-2. * sig1 + sig2 + sig3 +
			2 * sqrt((pow(fb, 2)*(-3. * pow(fa, 2) + 3. * pow(fapex, 2) - 2. * sqrt(3.)*fapex*xi + pow(xi, 2))) / pow(fa, 2.))*cos(beta), 2.)) / fG +
			(3 * pow(-3. * sig2 + 3. * sig3 + 2 * sqrt(3.)*sqrt((pow(fb, 2)*(-3. * pow(fa, 2) + 3. * pow(fapex, 2) - 2. * sqrt(3.)*fapex*xi + pow(xi, 2))) / pow(fa, 2))*sin(beta), 2)) / fG) / 108.;
		return func2;
	}

	Doub operator()(const Doub &x)
	{
		Doub sig1 = pt[0][0], sig2 = pt[0][1], sig3 = pt[0][2];
		Doub xi = x;
		if (xi > fapex) {
			//xi = fapex;
		}
		Doub beta = atan((sqrt(3)*(-sig2 + sig3)) / (-2 * sig1 + sig2 + sig3));
		Doub G = fG, K = fK, a = fa, b = fb, apex = fapex;
		Doub func2 = ((-24 * (sqrt(3)*sig1 + sqrt(3)*sig2 + sqrt(3)*sig3 - 3 * xi)) / K + (18 * pow(b, 2)*(-2 * sqrt(3)*apex + 2 * xi)*cos(beta)*
			(-2 * sig1 + sig2 + sig3 + 2 * sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*cos(beta))) /
			(pow(a, 2)*G*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))) +
			(6 * sqrt(3)*pow(b, 2)*(-2 * sqrt(3)*apex + 2 * xi)*sin(beta)*(-3 * sig2 + 3 * sig3 +
				2 * sqrt(3)*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*sin(beta))) /
			(pow(a, 2)*G*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2)))) / 108.;
		return func2;
	}

	void df(VecDoub_I &x, VecDoub_O &deriv)
	{
		Doub sig1 = pt[0][0], sig2 = pt[0][1], sig3 = pt[0][2];
		Doub xi = x[0];
		if (xi > fapex) {
			//xi = fapex;
		}
		Doub beta = atan((sqrt(3)*(-sig2 + sig3)) / (-2 * sig1 + sig2 + sig3));
		Doub G = fG,K = fK, a = fa, b = fb, apex = fapex;
		deriv[0] = ((-24 * (sqrt(3)*sig1 + sqrt(3)*sig2 + sqrt(3)*sig3 - 3 * xi)) / K + (18 * pow(b, 2)*(-2 * sqrt(3)*apex + 2 * xi)*cos(beta)*
			(-2 * sig1 + sig2 + sig3 + 2 * sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*cos(beta))) /
			(pow(a, 2)*G*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))) +
			(6 * sqrt(3)*pow(b, 2)*(-2 * sqrt(3)*apex + 2 * xi)*sin(beta)*(-3 * sig2 + 3 * sig3 +
				2 * sqrt(3)*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*sin(beta))) /
			(pow(a, 2)*G*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2)))) / 108.;
	}


	void setpt(MatDoub ptext, Doub young, Doub nu, Doub coesion, Doub frictionangle)
	{
		pt = ptext;
		setup(young,nu,coesion,frictionangle);
	}
	inline void setup(Doub young, Doub nu, Doub coesion, Doub frictionangle) {
		fyoung = young;
		fnu = nu;
		fcoesion = coesion;
		fphi = frictionangle;
		fG = young / (2.* (1. + nu));
		fK = young / (3.* (1. - 2.* nu));
		ftanphi = (3.*tan(fphi)) / sqrt(9. + 12.* tan(fphi) *tan(fphi));
		fapex = fcoesion * 1. / tan(fphi);
		fa = (fcoesion / (sqrt(3)*tan(fphi)) - fapex);
		fb = fa*ftanphi;

	}
};


struct Funcd :druckerprager {
	MatDoub pt;
	Doub fyoung;
	Doub fnu;
	Doub fcoesion;
	Doub fphi;
	Doub fK;
	Doub fG;
	Doub ftanphi;
	Doub fapex;
	Doub fa, fb;
	Doub operator() (const Doub x) {
		Doub sig1 = pt[0][0], sig2 = pt[0][1], sig3 = pt[0][2];
		Doub xi = x;
		if (xi > fapex) {
			//xi = fapex;
		}
		Doub beta = atan((sqrt(3)*(-sig2 + sig3)) / (-2 * sig1 + sig2 + sig3));
		Doub G = fG, K = fK, a = fa, b = fb, apex = fapex;
		return ((-24 * (sqrt(3)*sig1 + sqrt(3)*sig2 + sqrt(3)*sig3 - 3 * xi)) / K + (18 * pow(b, 2)*(-2 * sqrt(3)*apex + 2 * xi)*cos(beta)*
			(-2 * sig1 + sig2 + sig3 + 2 * sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*cos(beta))) /
			(pow(a, 2)*G*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))) +
			(6 * sqrt(3)*pow(b, 2)*(-2 * sqrt(3)*apex + 2 * xi)*sin(beta)*(-3 * sig2 + 3 * sig3 +
				2 * sqrt(3)*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*sin(beta))) /
			(pow(a, 2)*G*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2)))) / 108.;
	}
	Doub df(const Doub x) {
		Doub sig1 = pt[0][0], sig2 = pt[0][1], sig3 = pt[0][2];
		Doub xi = x;
		if (xi > fapex) {
			//xi = fapex;
		}
		Doub beta = atan((sqrt(3)*(-sig2 + sig3)) / (-2 * sig1 + sig2 + sig3));
		Doub G = fG, K = fK, a = fa, b = fb, apex = fapex;
		return (72 / K + (72 * pow(b, 2)*pow(-(sqrt(3)*apex) + xi, 2)*pow(cos(beta), 2)) / (pow(a, 2)*G*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) -
			(9 * pow(b, 4)*pow(2 * sqrt(3)*apex - 2 * xi, 2)*cos(beta)*(-2 * sig1 + sig2 + sig3 + 2 * sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*cos(beta))) /
			(pow(a, 4)*G*pow((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2), 1.5)) +
			(36 * pow(b, 2)*cos(beta)*(-2 * sig1 + sig2 + sig3 + 2 * sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*cos(beta))) /
			(pow(a, 2)*G*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))) +
			(72 * pow(b, 2)*pow(-(sqrt(3)*apex) + xi, 2)*pow(sin(beta), 2)) / (pow(a, 2)*G*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) -
			(3 * sqrt(3)*pow(b, 4)*pow(2 * sqrt(3)*apex - 2 * xi, 2)*sin(beta)*(-3 * sig2 + 3 * sig3 +
				2 * sqrt(3)*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*sin(beta))) /
			(pow(a, 4)*G*pow((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2), 1.5)) +
			(12 * sqrt(3)*pow(b, 2)*sin(beta)*(-3 * sig2 + 3 * sig3 + 2 * sqrt(3)*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*sin(beta))) /
			(pow(a, 2)*G*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2)))) / 108.;
	}

	void setpt(MatDoub ptext, Doub young, Doub nu, Doub coesion, Doub frictionangle)
	{
		pt = ptext;
		setup(young, nu, coesion, frictionangle);
	}
	inline void setup(Doub young, Doub nu, Doub coesion, Doub frictionangle) {
		fyoung = young;
		fnu = nu;
		fcoesion = coesion;
		fphi = frictionangle;
		fG = young / (2.* (1. + nu));
		fK = young / (3.* (1. - 2.* nu));
		ftanphi = (3.*tan(fphi)) / sqrt(9. + 12.* tan(fphi) *tan(fphi));
		fapex = fcoesion * 1. / tan(fphi);
		fa = (fcoesion / (sqrt(3)*tan(fphi)) - fapex);
		fb = fa*ftanphi;

	}
};


