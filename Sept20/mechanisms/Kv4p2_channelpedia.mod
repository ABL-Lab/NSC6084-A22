:[$URL: https://bbpteam.epfl.ch/svn/analysis/trunk/IonChannel/xmlTomod/CreateMOD.c $]
:[$Revision: 1499 $]
:[$Date: 2012-01-28 10:45:44 +0100 (Sat, 28 Jan 2012) $]
:[$Author: rajnish $]
:Comment :
:Reference :Properties of voltage-gated potassium currents in nucleated patches from large layer 5 cortical pyramidal neurons of the rat. J. Physiol. (Lond.), 2000, 525 Pt 3, 593-609
: EM - added q10 to correct for temperature

NEURON	{
	SUFFIX Kv4_2_0016
	USEION k READ ek WRITE ik
	RANGE gKv4_2bar, gKv4_2, ik, BBiD, q10 
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	celsius
	gKv4_2bar = 0.00001 (S/cm2) 
	BBiD = 40 
	q10=3
}

ASSIGNED	{
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	gKv4_2	(S/cm2)
	mInf
	mTau
	hInf
	hTau
}

STATE	{ 
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gKv4_2 = gKv4_2bar*m*m*m*h
	ik = gKv4_2*(v-ek)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
	h' = (hInf-h)/hTau
}

INITIAL{
	rates()
	m = mInf
	h = hInf
}

PROCEDURE rates(){
	LOCAL qt
	qt=q10^((celsius-23)/10)
	UNITSOFF 
		mInf = (1/(1 + exp((v- -18.8)/-16.6)))
		if(v < -50){
			mTau = 1.0/((0.026* exp(-0.026*v)) + (35* exp(0.136*v)))
		}
		if(v >= -50){
			mTau = 1.7/(1+ exp((-42 - v)/-26)) + 0.34
		} 
        mTau = mTau/qt
		hInf = 1/(1 + exp((v - -81.6)/6.7)) 
		hTau = (0.01*v + 6.7)/qt
	UNITSON
}
