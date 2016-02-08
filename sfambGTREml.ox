#include <oxstd.h>
#include <oxfloat.h>
#import <modelbase>
#import <maximize>
#import <maxsqp>
#include "sfambGTREml.h"

//////-----------Template for GTRE model (Colombi et al 2014) ------------------


/*   01-07-2015

	 Feb-2016
*/


/*-------------------------- Sfa : Database ---------------------------*/
Sfa::Sfa()
 {
    this->Modelbase();           // intialize base class
    m_iModel = 0;               // 1 when estimated
    m_bUse_maxSQP = m_bUse_maxSQPF = m_fPrintDetails = FALSE;
    m_iResult = -1;
    m_dLoglikOLS = 0;
    m_par = 0;
    m_dAlpha = 0.05;
    m_vStart = 0;
    m_iTl = 0;
    m_bUseRobustStdErr = TRUE;
    m_vUpper = m_vLower = <>;
    m_fcon_ge0 = m_fcon_eq0 = 0;
    if (m_fPrintDetails)
        print("Sfa class, object created on ",date(), ".\n\n");
}

Sfa::Ident(const vID, const vPer)
{
//----------------- Firms ID -----------------------
	decl vtid, i, vid;
	vtid = zeros(sizer(vID),1);//empty matrix
	for (i = 1; i < sizer(vID); ++i)
	{
		vtid[i][] = vID[i][] - vID[i-1][];//to find the first obs. of each firm	(note: value of the first firm is zero!)
	}									  
	vid = selectifr(vtid, vtid .!=0);//get rid of the zeros
	decl vpid = zeros(sizer(vid),1);
	for (i = 0; i < sizer(vid); ++i)
	{
		vpid[i][] = i+2;//"real" ID, i.e. the identifiers 1,2,3...  
	}
	m_ID = 1 | vpid;//[#firms][1]
//--------- ID vector for all observations (expansion)
	decl vFirm = vID[0][] | selectifr(vID, vtid .!= 0);//vector of original firm variable (FirmA, FirmB...)
	decl mData = vFirm ~ m_ID;
	decl vEmp, vRiden = zeros(sizer(vID),1);
	for (i=0; i < sizer(vID); i++)//LoopL[#obs]
	{
	vEmp = selectifr(mData, mData[][0] .== vID[i][]);
	if(sizer(vEmp)>1) {println("\nPanel structure incorrectly specified!",
							   "\nOrder of groups");exit(1);}
	vRiden[i][] = vEmp[][1];
	}
	m_vID = vRiden;//[#obs][1]
	m_vPer = vPer;//orig. data
//---------get individual T-i and expand-------------
	decl mFirm, mTi = zeros(maxc(m_vID), 1);
	for (i = 0; i < maxc(m_vID); i++)//LoopL[#firms] 
	{
	mFirm = selectifr(m_vPer, m_vID .== i+1);
	mTi[i][] = sizer(mFirm);//individual T-i
	if(mTi[i][]<2){
		println("\nPanel structure incorrectly specified!",
				"\nOnly one observation for group ", i+1);exit(1);}
	}
	m_Ti = mTi;
	decl vEm, vAllTi = zeros(sizer(m_vID),1);//expand it to m_vTi
	decl mData2 = m_ID ~ m_Ti;//#firms
	for (i=0; i < sizer(m_vID); i++)//LoopL[#obs], expansion
	{
	vEm = selectifr(mData2, mData2[][0] .== m_vID[i][]);
	vAllTi[i][] = vEm[][1];
	}
	m_vTi = vAllTi;//T-i as a [#obs][1] vector
}

Sfa::DropGroupIf(const mifr)
{
	decl mSel = selectifr(IDandPer(), mifr);

	decl i, vDiff, vDrop = zeros(sizer(mSel),1);

	for (i = 1; i < sizer(mSel); i++)
			{vDiff = mSel[i][0] - mSel[i-1][0];

			if(vDiff .!= 0)
			{vDrop[i][] = mSel[i][0];}}//find IDs

	vDrop = mSel[0][0] | vDrop[1:][];//append first ID

	vDrop = deleteifr(vDrop, vDrop .== 0);//delete zeros
	
	//expand:
	decl vIndex = zeros(sizer(m_vID), 1);
	
	for (i = 0; i < sizer(m_vID); i++){
		vIndex[i][] = selectifr(vDrop, vDrop .== m_vID[i][]);}

	//shorten the vectors of identification
	decl mvidnew = deleteifr(m_vID, vIndex .!= 0);
	decl mvpernew = deleteifr(m_vPer, vIndex .!= 0);

	//shorten the database
	m_mData = deleteifr(m_mData, vIndex .!= 0);

	println("\nRestricted Sample:");
	
	Ident(mvidnew, mvpernew);//new call to 'Ident()' to renew the corresp. data members
}

Sfa::SetPrintSfa(const iPri)
{
	print("\n\t  #groups:   #periods(max):  avg.T-i:");
	print(maxc(m_vID) ~ maxc(m_vPer) - minc(m_vPer) + 1 ~ sumc(m_Ti)/maxc(m_ID));

	SetPrint(iPri);
}

Sfa::IDandPer()	 
{
	return m_vID ~ m_vTi;//Firm number ~ individual T-i; [#obs][2]
}

Sfa::PrepData(const mSel, iNorm) 
{
	decl mVar;
	
	if (iNorm == TRUE){
		mVar = log(mSel/meanc(mSel));}
	
	else {mVar = log(mSel);}

	return mVar;
}

Sfa::GetPackageName(){return "Sfa";}

Sfa::GetPackageVersion(){return "1.0";}

Sfa::SetConstant() 
{
	Deterministic(-1);// (-1) == no seasonals
	
}

Sfa::Select(const iGroup, const aSel)  // RE model uses "Constant" !
{
	decl aSelnew;

//	if (iGroup != Y_VAR)
//	{
//		decl iIndex = strfind(aSel, "Constant");//relevant index of "Constant"
//
//		if (iIndex == -1) //in case "Constant" is not specified
//		{
//		aSelnew = aSel;
//		}
//		else
//		{
//		aSelnew = dropr(aSel, iIndex~iIndex+1~iIndex+2);//drop "Constant" and related zeros
//		}
//	}
//	else
//	{
		aSelnew = aSel;
//	}

	decl asVars = Database::Select(iGroup, aSelnew);

	return asVars;
}

Sfa::SetStart(const vStart)
{
        m_vStart = vStart;
        return TRUE;
}

Sfa::SetTranslog(const iTl){
        m_iTl = iTl;
        return TRUE;
}

Sfa::SetConfidenceLevel(const alpha){
        m_dAlpha = alpha;
        return TRUE;
}

Sfa::SetLowerBounds(const vBounds){
  m_bUse_maxSQP=TRUE;
  if (m_vUpper != <>){
    if (m_vUpper < vBounds){
      eprint("At least one upper bound is below the lower bound\n");
      eprint("no lower bounds have been applied");
      return FALSE;
    }
  }
  m_vLower = vBounds;
  return FALSE;
}

Sfa::SetUpperBounds(const vBounds){
  m_bUse_maxSQP=TRUE;
  if (m_vLower != <>){
    if (m_vLower > vBounds){
      eprint("At least one lower bound is below the lower bound\n");
      eprint("no upper bounds have been applied");
      return FALSE;
    }
  }
  m_vUpper = vBounds;
  return FALSE;
}

Sfa::fEqRest(const fIn){
  m_bUse_maxSQP=TRUE;
  m_fcon_eq0=fIn;
  return TRUE;
}

Sfa::fIneqRest(const fIn){
  m_bUse_maxSQP=TRUE;
  m_fcon_ge0=fIn;
  return TRUE;
}

/*
LLF in Proposition 2 of CKMV (2014)
*/

Sfa::fSfa(const vP, const prob, const avScore, const amHessian)
{
	decl epsi  = m_var[][0]-m_var[][1:m_cX]*vP[0:(m_cX-1)]; //NTx1 (including beta_0 )

	decl svsq = exp(vP[m_cX]);//sig-v^2 == sigma_e^2 (noise)  (note!)
	decl susq = exp(vP[m_cX+1]);//sig-u^2

	decl su2fix = exp(vP[m_cX+2]);//sig-u^2(fix, i.e. persistent component)
	
	decl sbsq = exp(vP[m_cX+3]);//sigma_b^2 (random effect) 

	
	decl i, veps, sumsqeps, epsbar, mustar, sigstarsq, lnL;

	decl lnfx = zeros(sizer(m_ID),1);

	lnL = zeros(sizer(m_ID),1);

/* matrices A, V, Sigma, Gamma */	
	decl mA, mV, mSigma, mGamma, mLambda, mR, vOut, n, j, k, iDim;

	decl lnFX = zeros(sizer(m_ID),1);

	
	for (i = 0; i < sizer(m_ID); i++)//loopL[#firms]
	{
	veps = selectifr(epsi, m_vID .== i+1);//Tx1, e_it (full error: eta + epsilon)

	mA = -(ones(m_Ti[i][],1) ~ unit(m_Ti[i][]));

	mV = (su2fix|zeros(m_Ti[i][],1)) ~ (zeros(1,m_Ti[i][])|(susq*unit(m_Ti[i][])));	

	mSigma = svsq*unit(m_Ti[i][]) + sbsq*ones(m_Ti[i][],m_Ti[i][]);

	mGamma = mSigma + mA*mV*mA';
	
	mLambda = invertgen(invertgen(mV) + mA' * invertgen(mSigma) * mA); /* covariance in CDF */

	mR = mLambda * mA' * invertgen(mSigma); /* skewness */

	iDim = rows(mLambda);             // 'iDim' could be replaced by data member 'm_Ti'

 	vOut = zeros(1, iDim*(iDim-1)/2);
	k=0;

  		for (n=0; n<(iDim-1); ++n){									//cf. 'veclower()'
    			for (j=n+1; j<iDim; ++j){
      			vOut[k] = mLambda[j][n]; ++k;}}
      			
/********** option, using m_Ti ******************
  vOut = zeros(1, m_Ti[i][]*(m_Ti[i][]-1)/2);
  k=0;

  		for (n=0; n<(m_Ti[i][]-1); ++n){									//cf. 'veclower()'
    			for (j=n+1; j<m_Ti[i][]; ++j){
      			vOut[k] = mLambda[j][n]; ++k;}}
*************************************/

//******************LLF*****************************
	lnfx[i][] =
				-(m_Ti[i][]/2) * log(2*M_PI)

				-0.5*log(determinant(mGamma))

				-0.5 * veps' * invertgen(mGamma) * veps;//pdf

				
	lnFX[i][] = log(probmvtdst((mR * veps)', vOut));//cdf

	}

	decl iC = sizer(m_ID)*(meanc(m_Ti)+1) * log(2);
	
	prob[0] = double(iC + sumc(lnfx) + sumc(lnFX));	
	
	return (!(isnan(prob[0])||isdotinf(prob[0]) ));      // 1 indicates success
}







Sfa::TE()
{
	decl my=0;
	decl u  = GetResiduals();
	decl svsq = exp(m_par[m_cX]);
	decl susq = exp(m_par[m_cX+1]);
	decl s  = sqrt(svsq + susq);
	decl lam =  sqrt(susq) / sqrt(svsq);
	decl gam = 1 ./ (1 + lam .^ 2);

    decl mystar = gam .*  - my + (1 - gam) .* ( - u);
    decl s2star = gam .* susq;

	return( (probn(mystar ./ sqrt(s2star) - sqrt(s2star))
         ./ probn(mystar ./ sqrt(s2star)) .* exp( - mystar + 0.5 .* s2star)));
}

Sfa::Ineff()
{
	decl my=0;
	decl u  = GetResiduals();//eps
	decl svsq = exp(m_par[m_cX]);
	decl susq = exp(m_par[m_cX+1]);
	decl s  = sqrt(svsq + susq);//sigma
	decl lam =  sqrt(susq) / sqrt(svsq);
	decl gam = 1 ./ (1 + lam .^ 2); //attention: gam is s_v^2/s, not s_u^2/s!!
	decl mystar = gam .*  - my + (1 - gam) .* ( - u);//cf. mustar of B&C1988,eq(9)
	decl s2star = gam .* susq;//sigma-star^2 = sigma_u^2*sigma_v^2/sigma^2; cf. B&C1988,eq(10)

	decl vjlms = sqrt(s2star) .* (				  //cf. K&L2000,eq(3.2.50)
								  mystar./sqrt(s2star)+
								  densn(mystar./(sqrt(s2star))) ./ (1-probn(-mystar./(sqrt(s2star))))
								  );
	return vjlms;
}

Sfa::GetResiduals()	
{
	return m_var[][0]-m_var[][1:m_cX]*m_par[0:(m_cX-1)];// returns w= v-u = y - \hat y
}

Sfa::InitData()
{
    decl  cp, vp, mh, i;

    m_iT1est = m_iT1sel;  m_iT2est = m_iT2sel;

	decl my = GetGroup(Y_VAR);
	decl mx = GetGroup(X_VAR);	
    //m_cX = columns(mx);
	decl cX = columns(mx);
//-------------------SetTranslog-----------------------------
	if (m_iTl)
    {   print("Constructing Squares and Cross-Products...");
        if (m_iTl>1) cX=m_iTl;
	mx=mx~0.5*(mx[][0:(cX-1)] .* mx[][0:(cX-1)]);//squares

	for (i = 0 ; i < cX-1 ; i++)
    {										
		mx=mx~(mx[][i+1:cX-1] .* mx[][i]);//cross-products
	}
        print("done.\n");

	m_mTldata = my~mx;//constructed data
	}
	
	m_cT = m_iT2est - m_iT1est + 1;
    if (m_cT <= 2)
    {   eprint("Only ", m_cT, " observations. This is not enough.\n");
        return FALSE;
    }

//////---------------------dummies-----------------
////	decl j, mDum = zeros(sizer(m_vID),sizer(m_ID));//NT x N
////	
////	for (i = 0; i < sizer(m_vID); i++){//[NT]rows 
////		for (j = 0; j < sizer(m_ID); j++)//[N]columns
////		{
////			if   (m_vID[i][] == m_ID[j][])
////			mDum[i][j]=mDum[i][j]+1;						
////		}}

////	m_var = my~mx~mDum;
////	m_cX = columns(mx)+columns(mDum);

	m_var = my~mx;
	m_cX = columns(mx);
	
	m_mY = my;
	m_cY = 1;

	decl mZ = GetGroup(Z_VAR);
    decl cZ = columns(mZ);

    decl mU = GetGroup(U_VAR);
    decl cU = columns(mU);

	if (cZ || cU){
		print("This model: neither Z-VAR nor U-VAR \n");
		return FALSE;
	}
		
	if (!m_cX)
    {   eprint("Need some regressors\n");
        return FALSE;
    }
	
    m_iModelStatus = MS_DATA;

    return TRUE;
}

Sfa::GridSearch(const beta, const s)
{
    decl i, fnow, fbest, fr, par_buffer;
    fnow = fbest = -10000;

    for (i = 0.05 ; i < 0.99 ; i += 0.01)
    {   fr = s/(1 - i * 0.6366198); //2/pi

			par_buffer =				 
			
//			(beta[0]+0.7978846*sqrt(fr*i))|			 
			beta[0:m_cX-1]
			|log((1-i)*fr)
//			|log(fr*i);
			|log(fr*i)
            |zeros(2, 1); //account for \sig_u,fix^2  and \sig_w^2

		if (!fSfa(par_buffer,&fnow,0,0))
        {   print("function evaluation failed at gamma ",i," !\n");
            fnow = -1000;
        }
        else if (fnow > fbest)
        {
            fbest = fnow;
            m_par = par_buffer;
        }
    }
	return TRUE;
}

// fcon(const avF, const vP){
//       avF[0]=matrix(1-vP[1]-vP[2]-vP[3]-vP[4]);
//       return 1;
//     }

Sfa::DoEstimation(vPar)
{
  decl  temp,cp,sd,gam,time;

  if (m_vStart==0)
    ols2c(m_var[][0], m_var[][1:m_cX], &temp);

  else 
    temp=m_vStart;

  sd=
    ((m_var[][0]-m_var[][1:m_cX]*temp[:m_cX-1])'(m_var[][0]-m_var[][1:m_cX]*temp[:m_cX-1]))/m_cT;//'
  m_dLoglikOLS = -m_cT/2*(log(sd)+2.837877066409); /* log-likelihood OLS*/

  if (sizeof(temp) == m_cX) {  
// check for length, if TRUE, only function pars are given, hence  perform gridsearch
	if (!GridSearch(temp,sd))
	  print("Error in GridSearch");
      }
  else m_par = temp;
  
////  if (m_fPrintDetails) 
////    print("starting values: ",
////			m_par[0:(m_cX-sizer(m_ID)-1)] | m_par[m_cX:(m_cX+1)]); //TFE-spec.

  if (m_fPrintDetails) 
    print("starting values: ", m_par);

			
  cp = sizer(m_par);
  SetParCount(cp);
  m_mCovP = unit(cp);
    
  time=timer();
  SetResult(5);

  SetResult(MaxBFGS(fSfa, &m_par, &m_dLogLik, &m_mCovP, TRUE));


//	  m_dLogLik *= m_cT; //not here!

	  if (m_fPrintDetails)
	    print("Elapsed time: ",timespan(time),"\n");
	  vPar=m_par;

		//SetResult(MAX_CONV); //uncomment to get SE's even after non-convergence
	  return vPar;
}

Sfa::Covar()				 
{
  m_mCovar = diag(ones(m_cPar,1));
    if (GetResult() == MAX_CONV || GetResult() == MAX_WEAK_CONV)
    {
        decl mXprod, dfunc;
        //fSfa(m_par, &dfunc, &mXprod, 0);
        //println(rows(mXprod)," ",columns(mXprod));
////////        mXprod = m_mScore * m_mScore' /m_cT; //'
        //print(rows(mXprod)," ",columns(mXprod));exit(0);
        if (Num2Derivative(fSfa, m_par, &m_mCovP))
        {
            m_mCovP = invertgen(-m_mCovP,30) ; //println(m_fPrintDetails);
			println("");
                if (m_fPrintDetails)
                    println("standard errors from Hessian.");
        }
        else
        {   m_mCovP = m_mCovP ;
                if (m_fPrintDetails)
                    println("standard errors from last update:");
        }

//		m_mCovar = m_mCovP / m_cT;//original sfamb
		m_mCovar = m_mCovP;//adjustment
		
//	if (m_iMethod==0)//orig.sfamb:
//		{
//		if (m_bUseRobustStdErr)
//            m_mCovarRobust = m_cT * m_mCovar * mXprod * m_mCovar ;
//		}			
    }
}

Sfa::GetParNames()
{
    decl i, j, iSize, asy = {}, asx = {};
	decl asvu = {"ln{\\sigma_v^2}","ln{\\sigma_u^2}"};
	decl asinvariant = {"ln{\\sig_u,fix^2}","ln{\\sigma_w^2}"};

	GetGroupLagNames(Y_VAR, 1, 100, &asy);
    if (m_cX) // regressors
        GetGroupNames(X_VAR, &asx);
    if (m_iTl){  // print squares and crossproducts, if necessary; adjusted for WT model
        iSize = (m_iTl > 1 ? m_iTl : sizeof(asx));
        for (i = 0 ; i < iSize; i++)    //print squares
                asx ~= sprint(".5*",asx[i],"^2");
        for (i = 0 ; i < iSize-1 ; i++)    //print cross terms
        {
				for (j = i; j < iSize-1; j++)
						asx ~= sprint(asx[i],"*",asx[j+1]);
        }
	}

	asx ~= asvu ~= asinvariant;
	
	for (i = sizeof(asx); i < m_cPar; ++i)
        asx ~= sprint("Par ", "-%2d", i + 1);
    return asy ~ asx;
}

Sfa::Output()
{
	println("-new model-");

    Modelbase::Output();

	decl sv = sqrt(exp(m_par[m_cX]));
	decl su = sqrt(exp(m_par[m_cX+1]));
	decl lambda = double(su/sv);
	println("lambda  ","%21.4g", lambda);	

}

Sfa::TestGraphicAnalysis()
{
    decl mEff = TE();
    // Histogram and boxplot
    DrawDensity(0, (mEff[][0])', "Technical efficiency", 0, 1, 0);
	DrawBoxPlot(1, (mEff[][0])', "Technical efficiency");
	ShowDrawWindow();
}

Sfa::SetPrintDetails(const bool)
{
  m_fPrintDetails=bool;
  return 1;
}

Sfa::SetUse_maxSQPF(const bool)
{
  if (bool)
    m_bUse_maxSQP=FALSE;
  m_bUse_maxSQPF=bool;
  return 1;
}

Sfa::SetRobustStdErr(const bool)
{
  m_bUseRobustStdErr = bool;
  return 1;
}

Sfa::AiHat()
{
	return m_par[(m_cX-sizer(m_ID)):m_cX-1];
}

Sfa::Elast(const sXname)
{
  decl i, maxParIdx, sTmp, vIdx, vData; //asX = {"Constant",sXname}, 
  decl asX = {sXname};
  decl asNames;
  asX = {sXname};
  asNames=GetParNames();
// get index of sigma_v^2 to exclude re-occuring X in Z or U
// when integrating this into the package, use a better var
  maxParIdx = strfind(asNames,"ln{\sigma_v^2}");
  if (maxParIdx == -1){
    eprint("Error with finding sigma_v in the list of vars, exiting...");
	exit(1);
  }
// we need the index of the first order term
  vIdx = strfind(asNames,sXname);
// find indices of second order terms
  for (i=vIdx+1;i< maxParIdx; ++i){
  	vIdx ~= strfind(asNames[i],sXname) == -1 ? <> : i;
  }
  for (i=2;i<sizerc(vIdx);++i){
    sTmp=asNames[vIdx[i]];//println(sTmp);
  	asX ~= strfind(sTmp,sXname) == 0
		? sTmp[sizeof(sXname)+1:]
		: sTmp[:strfind(sTmp,"*")-1];
  }
  vData = ones(GetcT(),1)~GetVar(asX);//vector of ones added
  decl vEstimates = GetPar()[vIdx];
//  decl sigma2, vCovar = GetCovarRobust()[vIdx][vIdx];//not for WT

	//04-06-2015

  decl sigma2;
  decl vCovar = GetCovar()[vIdx][vIdx];
  decl vTvalue=zeros(GetcT(),1);
	for (i=0;i<GetcT();i++){
	  sigma2 = (vData[i][] * vCovar * (vData[i][])');
	  vTvalue[i]=((vData * vEstimates)[i]/sqrt(sigma2));
	}
  
  return (vData * vEstimates)~(vTvalue);
}

Sfa::GetTldata()
{
  return m_mTldata;
}

Sfa::GetLLFi()
{
	decl i, vLLFi = zeros(maxc(m_vID),1);

	for (i = 0; i < maxc(m_vID); i++){ 
		vLLFi[i][] = m_dLogLik/m_cT * m_Ti[i][];}
	return vLLFi;
}

Sfa::GetMeans()
{
  return m_mMeans;
}

Sfa::GetWithins()
{
  return m_mWithins;
}