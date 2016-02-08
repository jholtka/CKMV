#include <oxstd.h>
#include <packages/gnudraw/gnudraw.h>
 //#import  <packages/sfamb/sfamb>
#include  <sfambRE.ox>

main(){
/*new object of class 'Sfa', Load data*/

	decl fob = new Sfa(); fob.Load("DataFrontierRiceProdPhilSort.xlsx");
//	fob.Info(); exit(1);

/**********************************************************************************
The output is annual rice production (measured in tons).
We include the
-area planted (measured in hectares),
-labor employed (measured in days per worker), and
-the amount of fertilizer (measured in kilograms) used as inputs.
**********************************************************************************/

	fob.SetConstant();
	
/*Identification of panel structure*/

	fob.Ident(fob.GetVar("FMERCODE"), fob.GetVar("YEARDUM"));


/*Data is already in logs, renew the names*/
	fob.Renew(fob.PrepData(fob.GetVar("PROD"), 0), "lny");
	
	fob.Renew(fob.PrepData(fob.GetVar("LABOR"), 0), "lnlab");
	fob.Renew(fob.PrepData(fob.GetVar("AREA"), 0), "lnarea");
	fob.Renew(fob.PrepData(fob.GetVar("NPK"), 0), "lnfert");
//	fob.Info(); exit(2);

/*Set up model*/
    fob.Select(Y_VAR, {"lny",0,0});			// Select dependent variable

    fob.Select(X_VAR, {						// Select regressors
						"Constant",0,0, 
						"lnlab",0,0, 
						"lnarea",0,0,
						"lnfert",0,0
						});

/**********************************************************************
Model specification
log-likelihood    -86.4304278
see Colombi et al., Table 6: TFT -86.430
**********************************************************************/
						
	// Select estimation sample
	fob.SetSelSample(-1, 1, -1, 1);  // full sample
	fob.SetPrintSfa(TRUE);
	MaxControl(1000,10,TRUE);

	fob.SetPrintDetails(1);
	
	fob.Estimate();

	delete fob;
}