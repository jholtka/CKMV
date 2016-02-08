#include <oxstd.h>
#include <packages/gnudraw/gnudraw.h>
 //#import  <packages/sfamb/sfamb>
#include  <sfambGTREml.ox>

main(){
/*new object of class 'Sfa', Load data*/

	decl fob = new Sfa(); fob.Load("demo_panel.xls");
//	fob.Info(); exit(1);

	fob.SetConstant();
	
/*Identification of panel structure*/

	fob.Ident(fob.GetVar("id"), fob.GetVar("time"));


/*Data is already in logs, renew the names*/
	
	fob.Renew(fob.GetVar("y"), "lny");

	fob.Renew(fob.GetVar("x1"), "lnx1");
	fob.Renew(fob.GetVar("x2"), "lnx2");
//	fob.Info(); exit(2);

/*Set up model*/
    fob.Select(Y_VAR, {"lny",0,0});			// Select dependent variable

    fob.Select(X_VAR, {						// Select regressors
						"Constant",0,0, 
						"lnx1",0,0, 
						"lnx2",0,0
//						"time",0,0
						});

	// Select estimation sample
	fob.SetSelSample(-1, 1, -1, 1);  // full sample
	fob.SetPrintSfa(TRUE);
	MaxControl(1000,10,TRUE);

	fob.SetPrintDetails(1);
	
	fob.Estimate();

	delete fob;
}