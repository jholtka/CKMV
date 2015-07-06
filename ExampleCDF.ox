#include <oxstd.h>
#include <oxfloat.h> //needed for M_PI
#include <oxprob.h> // for probmvn

/* Example for dimensions in CDF of Proposition 2 */

main()
{
//	decl m_Ti = 3;
	decl m_Ti = 2;  //max. dim. for 'probmvn' is 3 = 2+1
	decl mA, mV, mSigma;
/* Matrix A */
	mA = -(ones(m_Ti,1)~unit(m_Ti));
//	println("Matrix A", mA);
/* Matrix V (including Psi=susq*unit(m_Ti)) */
	decl susq = 2; 	//conv. sigma_u^2 (time-varying)
	decl sufix = 10; //persistent
	mV = (sufix|zeros(m_Ti,1)) ~ (zeros(1,m_Ti)|(susq*unit(m_Ti)));	
//	println("Matrix V", mV);
/* Matrix Sigma */
	decl sesq = 5; //sigma_e^2 (noise)
	decl sbsq = 15; //sigma_b^2 (random effect)
	mSigma = sesq*unit(m_Ti) + sbsq*ones(m_Ti,m_Ti);
//	println("Matrix Sigma", mSigma);
/* Matrix Lambda = covariance in CDF */
	decl mLambda;
	mLambda = invertgen(invertgen(mV) + mA' * invertgen(mSigma) * mA);
//	println("Matrix mLambda", mLambda);
/* Matrix R -- skewness */
	decl mR;
	mR = mLambda * mA' * invertgen(mSigma);
//	println("Matrix mR", mR);

/* random error, vector(Tx1) */
	decl veps;
	veps = rann(m_Ti, 1);
	println("veps ", veps);
/* input to CDF */
	decl vx = mR * veps;  // R* (y_i - X_i \beta - 1_T \beta_0)

	println("vx' = (mR * veps)' ", vx'); //here, transpose. row vector as input
	println("result from probmvn \n", probmvn(vx', mLambda));
}