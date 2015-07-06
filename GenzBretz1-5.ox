/*Example,Genz&Bretz:p.4*/
#include <oxstd.h>
#include <oxfloat.h>
#include <oxprob.h>

main()
{
	//(X1,X2,X3)
	decl mx = <1,4,2>; //rann(const r, const c);
	decl msigma = <
	  1,   3/5,   1/3;
	3/5,     1, 11/15;
	1/3, 11/15,     1;>;
	
	print(mx);  //row vector (input!)

	print(msigma);

	println("probmvn ",probmvn(mx, msigma));

}