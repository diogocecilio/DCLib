#include "druckerprager.h"
#include <math.h>
#include "ludcmp.h"
#include "eigen_sym.h"
#include <math.h>
#include <cmath>
#include "cholesky.h"
#include "mins_ndim.h"
druckerprager::druckerprager()
{
}

druckerprager::druckerprager(Doub young, Doub nu, Doub coesion, Doub frictionangle)
{
	setup(young, nu,coesion,frictionangle);
}


druckerprager::~druckerprager()
{
}

void druckerprager::closestpointproj(TensorDoub epst, TensorDoub epsp, TensorDoub & projstress, TensorDoub & projstrain, MatDoub & Dep, Doub & projgamma)
{
	TensorDoub epse = epst - epsp;
	MatDoub C = GetElasticMatrix();
	//C.Print();
	//epse.Print();

	MatDoub tempepsemat, stresstrial, Dept;
	epse.FromTensorToNRmatrix(tempepsemat);
	//std::cout << "tempepsemat (elastic strain MatDoub) = " << std::endl;
	//tempepsemat.Print();
	C.Mult(tempepsemat, stresstrial);
	//stresstrial.Print();
	TensorDoub stresstrialtensor;
	epse.FromNRmatrixToTensor(stresstrial, stresstrialtensor);
	Doub I1, J2;
	J2 = stresstrialtensor.J2();
	I1 = stresstrialtensor.I1();
	Doub xi = I1 / sqrt(3);
	Doub rho = sqrt(2 * J2);
	Doub yieldcr = yield(xi, rho);
	if (yieldcr < 0) {

		projstress = stresstrialtensor;
		projstrain = epse;
		Dept = C;
		Dep.assign(3, 3, 0.);
		Dep[0][0] = Dept[0][0];Dep[0][1] = Dept[0][1];Dep[0][2] = Dept[0][5];
		Dep[1][0] = Dept[1][0];Dep[1][1] = Dept[1][1];Dep[1][2] = Dept[1][5];
		Dep[2][0] = Dept[5][0];Dep[2][1] = Dept[5][1];Dep[2][2] = Dept[5][5];
		projgamma = 0.;
	}
	else {
		Doub xitrial = I1 / sqrt(3.);
		Doub rhotrial = sqrt(2. * J2);
		//std::cout << "rhotrial = " << rhotrial << std::endl;
		MatDoub pt, vec, vect;
		stresstrialtensor.EigenSystem(pt, vec);
		//std::cout << "autovalores  \n";
		//pt.Print();
		//std::cout << "autovetores  \n";
		//vec.Print();
		Doub sig1 = pt[0][0], sig2 = pt[0][1], sig3 = pt[0][2];
		Doub betasol = atan((sqrt(3.)*(-sig2 + sig3)) / (-2.* sig1 + sig2 + sig3));
		//std::cout << "betasol  = "<< betasol <<std::endl;
		
		VecDoub initguess(1,0.);
		initguess[0] = xitrial;
		Funcstruct func;
		func.setpt(pt);
		Powell<Funcstruct> powell(func);
		VecDoub xisolvec = powell.minimize(initguess);
		Doub xisol = xisolvec[0];


		MatDoub sig = HW(F1HWCylDruckerPragerSmoothPSMATCH(xisol, rhotrial, betasol));

		//std::cout << "sig  = "  << std::endl;
		//sig.Print();

		MatDoub fulltensorproj = stressrecosntruction(sig, vec);
		//std::cout << "fullprojtensor  = "  << std::endl;
		//fulltensorproj.Print();
		MatDoub projVoigtMat, invCe = GetInverseElasticMatrix(), epsemat;
		fulltensorproj.FromFullToVoigt(projVoigtMat);
		//std::cout << "projVoigtMat  = " << std::endl;
		//projVoigtMat.Print();

		TensorDoub voigtProjTensor, nvec, epsptensor;
		nvec.FromNRmatrixToTensor(projVoigtMat, voigtProjTensor);
		projstress = voigtProjTensor;

		voigtProjTensor.ComputeS(nvec);

		//std::cout << "ComputeS  = " << std::endl;
		//nvec.Print();

		nvec *= sqrt(3.) / (2.* sqrt(voigtProjTensor.J2()));

		//std::cout << "nvec  = " << std::endl;
		//nvec.Print();

		invCe.Mult(projVoigtMat, epsemat);
		MatDoub epspmat = tempepsemat;
		epspmat -= epsemat;
		nvec.FromNRmatrixToTensor(epspmat, epsptensor);

		nvec.FromNRmatrixToTensor(epsemat, projstrain);
		Doub gamma = epsptensor.Norm() / nvec.Norm();
		projgamma = gamma;
		//std::cout << "gamma  = " << gamma  <<std::endl;

		MatDoub dnvecdsig = dadsig(voigtProjTensor);

		//std::cout << "dnvecdsig  = " << std::endl;

		//dnvecdsig.Print();

		MatDoub Q(6, 6, 0.), invQ, R;for (Int i = 0;i < 6;i++)Q[i][i] = 1.;
		MatDoub sec;
		C.Mult(dnvecdsig, sec);
		sec *= gamma;
		Q += sec;
		//std::cout << "Q  = " << std::endl;
		//Q.Print();

		LUdcmp *lu = new LUdcmp(Q);
		lu->inverse(invQ);
		//std::cout << "invQ  = " << std::endl;
		//invQ.Print();

		invQ.Mult(C, R);

		//std::cout << "R  = " << std::endl;
		//R.Print();

		Dept = R;
		MatDoub asol, tempprod, tempprodT, temp2;
		nvec.FromTensorToNRmatrix(asol);
		R.Mult(asol, tempprod);
		Doub sum = 0.;
		for (Int i = 0;i < asol.nrows();i++)sum += asol[i][0] * tempprod[i][0];

		//std::cout << "sum = " << sum  <<std::endl;
		tempprod.Transpose(tempprodT);
		tempprod.Mult(tempprodT, temp2);
		//std::cout << "Outer[a.R,a.R] = " << std::endl;
		//temp2.Print();

		temp2 *= 1. / sum;
		Dept -= temp2;

		Dep.assign(3, 3, 0.);
		Dep[0][0] = Dept[0][0];Dep[0][1] = Dept[0][1];Dep[0][2] = Dept[0][5];
		Dep[1][0] = Dept[1][0];Dep[1][1] = Dept[1][1];Dep[1][2] = Dept[1][5];
		Dep[2][0] = Dept[5][0];Dep[2][1] = Dept[5][1];Dep[2][2] = Dept[5][5];

		/*Dep2D = { { Dep[[1, 1]], Dep[[1, 2]], Dep[[1, 6]] },{ Dep[[2, 1]],
		Dep[[2, 2]], Dep[[2, 6]] },{ Dep[[6, 1]], Dep[[6, 2]],
		Dep[[6, 6]] } };
		sigproj2D = { sigprojvoigth[[1]], sigprojvoigth[[2]],
		sigprojvoigth[[6]] };
		*/
		//std::cout << "DEP = " << std::endl;
		//Dep.Print();
	}
}



Doub druckerprager::yield(Doub xi,Doub rho)
{
	Doub xiint;
	if (xi > fapex)
	{
		 xiint = fapex;
	}
	else {
		xiint = xi;
	}
	return 1 + pow(rho, 2) / (2.*pow(fb, 2)) - pow(fapex - xiint/ sqrt(3), 2) / pow(fa, 2);

}

MatDoub druckerprager::GetElasticMatrix()
{
	MatDoub C(6, 6, 0.);
	Doub G = fG, K = fK;
	C[0][0] = (4 * G) / 3 + K;    C[0][1] = -((2 * G) / 3) + K;C[0][2] = -((2 * G) / 3) + K;C[0][3] = 0.;C[0][4] = 0.;C[0][5] = 0.;
	C[1][0] = -((2 * G) / 3) + K; C[1][1] = (4 * G) / 3 + K;C[1][2] = -((2 * G) / 3) + K;C[1][3] = 0.;C[1][4] = 0.;C[1][5] = 0.;
	C[2][0] = -((2 * G) / 3) + K; C[2][1] = -((2 * G) / 3) + K;C[2][2] = (4 * G) / 3 + K;C[2][3] = 0.;C[2][4] = 0.;C[2][5] = 0.;
	C[3][0] = 0;                  C[3][1] = 0;                 C[3][2] = 0;                 C[3][3] = G; C[3][4] = 0.;C[3][5] = 0.;
	C[4][0] = 0;                  C[4][1] = 0;                 C[4][2] = 0;                 C[4][3] = 0.;C[4][4] = G; C[4][5] = 0.;
	C[5][0] = 0;                  C[5][1] = 0;                 C[5][2] = 0;                 C[5][3] = 0.;C[5][4] = 0.;C[5][5] = G;
	return C;
}
MatDoub druckerprager::GetInverseElasticMatrix()
{
	MatDoub C(6, 6, 0.);
	Doub G = fG, K = fK;
	C[0][0] = (G + 3 * K) / (9.*G*K);    C[0][1] = -1 / (6.*G) + 1 / (9.*K);C[0][2] = -1 / (6.*G) + 1 / (9.*K);C[0][3] = 0.;C[0][4] = 0.;C[0][5] = 0.;
	C[1][0] = -1 / (6.*G) + 1 / (9.*K);  C[1][1] = (G + 3 * K) / (9.*G*K);  C[1][2] = -1 / (6.*G) + 1 / (9.*K);C[1][3] = 0.;C[1][4] = 0.;C[1][5] = 0.;
	C[2][0] = -1 / (6.*G) + 1 / (9.*K);  C[2][1] = -1 / (6.*G) + 1 / (9.*K);C[2][2] = (G + 3 * K) / (9.*G*K);  C[2][3] = 0.;C[2][4] = 0.;C[2][5] = 0.;
	C[3][0] = 0;                         C[3][1] = 0;                       C[3][2] = 0;                       C[3][3] = 1. / G;C[3][4] = 0.;C[3][5] = 0.;
	C[4][0] = 0;                         C[4][1] = 0;                       C[4][2] = 0;                       C[4][3] = 0.;C[4][4] = 1. / G;C[4][5] = 0.;
	C[5][0] = 0;                         C[5][1] = 0;                       C[5][2] = 0;                       C[5][3] = 0.;C[5][4] = 0.;C[5][5] = 1. / G;
	return C;
}
MatDoub druckerprager::F1HWCylDruckerPragerSmoothPSMATCH(Doub xi, Doub rho, Doub beta)
{
	MatDoub solproj(3, 1);


	Doub xiint;
	if (xi > fapex)
	{
		xiint = fapex;
	}
	else {
		xiint = xi;
	}
	Doub rhoint = sqrt(0.6666666666666666)*sqrt(-((pow(fb, 2)*(3 * pow(fa, 2) - 3 * pow(fapex, 2) + 2 * sqrt(3)*fapex*xiint - pow(xiint, 2))) / pow(fa, 2)));
	solproj[0][0] = xiint;
	solproj[1][0] = rhoint;
	solproj[2][0] = beta;
	return solproj;
}
MatDoub druckerprager::HW(MatDoub sig)
{
	MatDoub sol(3, 1);
	Doub xi, rho, beta;
	xi = sig[0][0];
	rho = sig[1][0];
	beta = sig[2][0];
	sol[0][0] = xi / sqrt(3.) + sqrt(2. / 3.)*rho*cos(beta);
	sol[1][0] = xi / sqrt(3.) + sqrt(2. / 3.)*rho*cos(beta - 2. *  M_PI / 3.);
	sol[2][0] = xi / sqrt(3.) + sqrt(2. / 3.)*rho*cos(beta + 2. *  M_PI / 3.);
	return sol;
}
MatDoub druckerprager::stressrecosntruction(MatDoub val, MatDoub vec)
{
	MatDoub sol(3, 3, 0.);
	for (Int i = 0;i < 3;i++)
	{
		MatDoub colvec(3, 1, 0.), colvect, temp;
		for (Int j = 0;j < 3;j++)
		{
			colvec[j][0] = vec[j][i];
		}

		colvec.Transpose(colvect);
		colvec.Mult(colvect, temp);
		temp *= val[i][0];
		sol += temp;
	}
	//sol.Print();
	return sol;
}
MatDoub druckerprager::dadsig(TensorDoub sigprojvoigt)
{
	//(P / b ^ 2) - (2 Outer[Times, Ii, Ii] / (9 a ^ 2))
	MatDoub  first = P();
	first *= 1/ (fb*fb);
	MatDoub second, I(6,1,0.),It(1,6,0.);
	I[0][0] = 1.;I[1][0] = 1.;I[2][0] = 1.;
	I.Transpose(It);
	I.Mult(It, second);
	second *= 2. / (9 * fa*fa);
	first -= second;
	return first;
}

void druckerprager::avec(MatDoub & avec)
{
	//nvec=(6 apex b ^ 2 Ii - 2 b ^ 2 ComputeI1[sigma] Ii + 9 a ^ 2 ComputeS[sigma]) / (9 a ^ 2 b ^ 2) //. subst2
}

MatDoub druckerprager::P()
{
	MatDoub P(6, 6);
	P[0][0] = 2. / 3.;    P[0][1] = -1. / 3.; P[0][2] = -1. / 3.; P[0][3] = 0.;P[0][4] = 0.;P[0][5] = 0.;
	P[1][0] = -1. / 3.;   P[1][1] = 2. / 3.;  P[1][2] = -1. / 3.; P[1][3] = 0.;P[1][4] = 0.;P[1][5] = 0.;
	P[2][0] = -1. / 3.;   P[2][1] = -1. / 3.;  P[2][2] = 2. / 3.; P[2][3] = 0.;P[2][4] = 0.;P[2][5] = 0.;
	P[3][0] = 0.;      P[3][1] = 0.;     P[3][2] = 0.;     P[3][3] = 2.;P[3][4] = 0.;P[3][5] = 0.;
	P[4][0] = 0.;       P[4][1] = 0.;      P[4][2] = 0.;      P[4][3] = 0.;P[4][4] = 2.;P[4][5] = 0.;
	P[5][0] = 0.;       P[5][1] = 0.;      P[5][2] = 0.;      P[5][3] = 0.;P[5][4] = 0.;P[5][5] = 2.;

	return P;
}
