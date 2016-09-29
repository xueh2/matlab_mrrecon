#include <string.h>
#include <omp.h>
#include "mex.h" 

#define DIMS 6
#define IDX(step) el[0]*step[0]+el[1]*step[1]+el[2]*step[2]+el[3]*step[3]+el[4]*step[4]+el[5]*step[5]
#define IDX_OMP(step) el0 * step[order[0]] + el1 * step[order[1]] + el2 * step[order[2]] + el3 * step[order[3]] + el4 * step[order[4]] + el5 * step[order[5]]
#define SUMDIM(x) x[0]+x[1]+x[2]+x[3]+x[4]+x[5]
#define PRODDIM(x) x[0]*x[1]*x[2]*x[3]*x[4]*x[5]
#define changeOrder(o1,o2,o3,o4,o5,o6) order[0]=o1;order[1]=o2;order[2]=o3;order[3]=o4;order[4]=o5;order[5]=o6;

double wavSoftTreshold(double* out_r,double* out_i,double* in_r,double* in_i,int* wavDim,double* weights,double step,int* order);

//to be mexed by compile_and_test_redundantHaar6D.m
//the use of the function is explained there
//
//soft tresholds the coefficients of wavelet coefficients
//inputs: - x a matrix up to 6D representing the wavelet coefficents of an image
//        - lambdas: a [2][2][2][2][2][2] dimensioned matrix giving the tresholds to apply for each quadrant for each combination of low and high frequencies
//output: - y, the result of the softtresholding of x in the wavelet dimension (before reprojection on the image domain
//
//CAREFUL: the weights for the last elements in high frequency is put to 0 as it loops with the first element
//it may (or may not) be necessary to do the same for the low frequencies
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    //input pointers
    const double* in_r;
    const double* in_i;
    const double* weights;
    double step;

    //output pointers
    double* out_r;
    double* out_i;
    double L1norm;
    double* normPtr;
    


    //input/outputdimensions
    mwSize numDim;
    mwSize numWeights;
    const mwSize* wavDimMatlab;
    const mwSize* weightsDimMatlab;
    int wavDim[DIMS];           //dimension of the wavelet data
    int order[DIMS]={5, 2, 3, 4, 0, 1};
    
    //loop
    int k;
    
    //data
    int complex;

    
    //use all the processors if possible
    omp_set_num_threads(omp_get_num_procs());

    //initialization of variables
    for (k=0;k<DIMS;++k)
        wavDim[k]=1;
    
    
    //Error management
    if ((nrhs<2)||(nrhs>3))
    	mexErrMsgTxt ("Either 2 or 3 inputs are accepted!\n"); 
    if (nlhs>3)
    	mexErrMsgTxt ("Only 1 or 2 output argument is allowed!\n"); 
    
    if (nrhs==3)
        step = mxGetScalar(prhs[2]);
    else
        step = 1.0;

    //check and prepare dimensions
    numDim=mxGetNumberOfDimensions (prhs[0]);
    wavDimMatlab=mxGetDimensions (prhs[0]);
    
    numWeights=mxGetNumberOfDimensions (prhs[1]);
    weightsDimMatlab=mxGetDimensions (prhs[1]);
    if (numWeights!=6)
        mexErrMsgTxt ("Weights need to be [2][2][2][2][2][2]!\n");
    for (k=0;k<6;++k)
        if (weightsDimMatlab[k]!=2)
            mexErrMsgTxt ("Weights need to be [2][2][2][2][2][2]!\n");
    
    complex=mxIsComplex(prhs[0]);
    for (k=0;k<numDim;k++)
	{
        wavDim[k]=wavDimMatlab[k];
		if (wavDimMatlab[k]>1)
			if (wavDimMatlab[k]%2)
				mexErrMsgTxt ("Wavelets dimensions are supposed to be even!\n");
	}
    
    //read input data
    in_r = mxGetPr(prhs[0]);
    in_i = (complex? mxGetPi(prhs[0]):0);
    weights = mxGetPr(prhs[1]);
        
    //create output data
    plhs[0] = (mxArray *) mxCreateNumericArray (numDim, wavDimMatlab, mxGetClassID (prhs[0]), complex);
    out_r=mxGetPr(plhs[0]);
    out_i=(complex? mxGetPi(plhs[0]):0);
    


    //soft tresholding
    L1norm=wavSoftTreshold(out_r,out_i,in_r,in_i,wavDim,weights,step,order);
    if (nlhs==2)
    {
        plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
        normPtr=mxGetPr(plhs[1]);
        normPtr[0]=L1norm;;
    }
    
}




///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
#define WDIM(hl,num) hl*weightStep[num]
#define SUBIND(hl,num) hl * wavDimSteps[order[num]]*dataDim[num]
#define COPYCOND(el,num,hl) ((el>(dataDim[num]-HLdim[num]))&&(hl==1))//from this condition comes the fact that wavelet weigths are 0 on the high frequency border
double wavSoftTreshold(double* out_r,double* out_i,double* in_r,double* in_i,int* wavDim,double* weights,double step,int* order)
{
    
    double L1norm=0;
    double coeff;
    double lambda;
    double tmp;
    int el0, el1, el2, el3, el4, el5,ind,k,weiInd,subindex,hl0,hl1,hl2,hl3,hl4,hl5;
    int weightStep[DIMS];
    int wavDimSteps[DIMS];
    int dataDim[DIMS];
    int HLdim[DIMS];

    wavDimSteps[0]=1;
    for (k=1;k<DIMS;++k)
        wavDimSteps[k]   = wavDimSteps[k-1]*wavDim[k-1];
    for (k=0;k<DIMS;++k)
    {
        dataDim[k]   = (wavDim[order[k]]>1?wavDim[order[k]]/2:1);
        weightStep[k]=(int)(pow(2,order[k]));
        HLdim[k]=(wavDim[order[k]]>1?2:1);
    }

    if (in_i)
    {
        
        
        
        for (hl0=0;hl0<HLdim[0];++hl0)
            for (hl1=0;hl1<HLdim[1];++hl1)
                for (hl2=0;hl2<HLdim[2];++hl2)
                    for (hl3=0;hl3<HLdim[3];++hl3)
                        for (hl4=0;hl4<HLdim[4];++hl4)
                            for (hl5=0;hl5<HLdim[5];++hl5)
                            {
                                subindex=SUBIND(hl0,0)+SUBIND(hl1,1)+SUBIND(hl2,2)+SUBIND(hl3,3)+SUBIND(hl4,4)+SUBIND(hl5,5);
                                lambda=weights[WDIM(hl0,0)+WDIM(hl1,1)+WDIM(hl2,2)+WDIM(hl3,3)+WDIM(hl4,4)+WDIM(hl5,5)];
                                tmp = lambda*step;
                                
                                //mexPrintf("lambda = %f\n",lambda);
								for (el0=0;el0<dataDim[0];++el0) //0 weight is given to the last element has the wavelet is circular (the last element is a combination of the first and the last, it doesnt make sense to weight it)
									for (el1=0;el1<dataDim[1];++el1)
										for (el2=0;el2<dataDim[2];++el2)
											for (el3=0;el3<dataDim[3];++el3)
											{
												#pragma omp parallel for private(el4, el5,coeff,ind)
												for (el4 = 0; el4 < dataDim[4]; ++el4)
												{
													for (el5 = 0; el5 < dataDim[5]; ++el5)
													{
														ind=subindex+IDX_OMP(wavDimSteps);

														if ((lambda==0)||COPYCOND(el0,0,hl0)||COPYCOND(el1,1,hl1)||COPYCOND(el2,2,hl2)||COPYCOND(el3,3,hl3)||COPYCOND(el4,4,hl4)||COPYCOND(el5,5,hl5))
														{
															out_r[ind]=in_r[ind];
															out_i[ind]=in_i[ind];
														}
														else
														{
															coeff=sqrt(in_r[ind]*in_r[ind]+in_i[ind]*in_i[ind]);
															if (coeff-tmp>0)
															{
																L1norm+=(coeff-tmp)*lambda;
																coeff=1-tmp/coeff;
																out_r[ind]=in_r[ind]*coeff;
																out_i[ind]=in_i[ind]*coeff;
															}
															else
															{
																out_r[ind]=0;
																out_i[ind]=0;
															}
														}


													}
												}
											}
                            }
    }
    else
    {
    	for (hl0=0;hl0<HLdim[0];++hl0)
            for (hl1=0;hl1<HLdim[1];++hl1)
                for (hl2=0;hl2<HLdim[2];++hl2)
                    for (hl3=0;hl3<HLdim[3];++hl3)
                        for (hl4=0;hl4<HLdim[4];++hl4)
                            for (hl5=0;hl5<HLdim[5];++hl5)
                            {
                                subindex=SUBIND(hl0,0)+SUBIND(hl1,1)+SUBIND(hl2,2)+SUBIND(hl3,3)+SUBIND(hl4,4)+SUBIND(hl5,5);
                                lambda=weights[WDIM(hl0,0)+WDIM(hl1,1)+WDIM(hl2,2)+WDIM(hl3,3)+WDIM(hl4,4)+WDIM(hl5,5)];
                                tmp=lambda*step;

								for (el0=0;el0<dataDim[0];++el0)
									for (el1=0;el1<dataDim[1];++el1)
										for (el2=0;el2<dataDim[2];++el2)
											for (el3=0;el3<dataDim[3];++el3)
											{
												#pragma omp parallel for private(el4, el5,coeff,ind)
												for (el4 = 0; el4 < dataDim[4]; ++el4)
												{
													for (el5 = 0; el5 < dataDim[5]; ++el5)
													{

														ind=subindex+IDX_OMP(wavDimSteps);
														//mexPrintf("index is %i .\n",ind);
														//mexPrintf("weight indices %i %i %i %i %i %i.\n",WDIM(el0,0),WDIM(el1,1),WDIM(el2,2),WDIM(el3,3),WDIM(el4,4),WDIM(el5,5));
														//mexPrintf("tresh is %f .\n",lambda);														
                                                        if ((lambda==0)||COPYCOND(el0,0,hl0)||COPYCOND(el1,1,hl1)||COPYCOND(el2,2,hl2)||COPYCOND(el3,3,hl3)||COPYCOND(el4,4,hl4)||COPYCOND(el5,5,hl5))
														{
															out_r[ind]=in_r[ind];
														}
														else
														{
															coeff=fabs(in_r[ind]);//sqrt(wav_r[ind]*wav_r[ind]);

															if (coeff-tmp>0)
															{
																out_r[ind]=in_r[ind]*(1-tmp/coeff);
																L1norm+=lambda*(coeff-tmp);
															}
															else
																out_r[ind]=0;
														}
													}
												}
											}
                            }
    }
    
    return L1norm;
    
}