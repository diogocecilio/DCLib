#include "druckerprager.h"



druckerprager::druckerprager()
{
}

druckerprager::druckerprager(Doub young, Doub nu, Doub coesion, Doub frictionangle)
{
	fyoung = young;
	fnu = nu;
	fcoesion = coesion;
	ffrictionangle = frictionangle;
	fG = young / (2.* (1. + nu));
	fK = young / (3.* (1. - 2.* nu));
}


druckerprager::~druckerprager()
{
}

void druckerprager::closestpointproj(TensorDoub epst, TensorDoub epsp, TensorDoub & projstress, TensorDoub & projstrain, MatDoub & Dep, Doub & projgamma)
{

}
Doub druckerprager::yield(Doub J2)
{

}
MatDoub druckerprager::GetElasticMatrix()
{

}
MatDoub druckerprager::GetInverseElasticMatrix()
{

}
MatDoub druckerprager::F1HWCylVonMises(Doub xisol, Doub rho, Doub betasol)
{

}
MatDoub druckerprager::HW(MatDoub sig)
{

}
MatDoub druckerprager::stressrecosntruction(MatDoub val, MatDoub vec)
{

}
MatDoub druckerprager::dadsig(TensorDoub sigprojvoigt)
{

}
MatDoub druckerprager::P()
{

}
