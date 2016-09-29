#include <string.h>
#include <omp.h>
#include "mex.h" 

#define DIMS 6
#define IDX(step) el[0]*step[0]+el[1]*step[1]+el[2]*step[2]+el[3]*step[3]+el[4]*step[4]+el[5]*step[5]
#define IDX_OMP(step) el0 * step[order[0]] + el1 * step[order[1]] + el2 * step[order[2]] + el3 * step[order[3]] + el4 * step[order[4]] + el5 * step[order[5]]
#define SUMDIM(x) x[0]+x[1]+x[2]+x[3]+x[4]+x[5]
#define PRODDIM(x) x[0]*x[1]*x[2]*x[3]*x[4]*x[5]
#define changeOrder(o1,o2,o3,o4,o5,o6) order[0]=o1;order[1]=o2;order[2]=o3;order[3]=o4;order[4]=o5;order[5]=o6;

double wavL1norm(double* in_r,double* in_i,int* wavDim,double* weights,int* order);


//to be mexed by compile_and_test_redundantHaar6D.m
//the use of the function is explained there
//
//computes the weighted L1 norm, given the redundant Haar wavelet coefficients
//inputs: - x a matrix up to 6D representing the wavelet coefficents of an image
//        - lambdas: a [2][2][2][2][2][2] dimensioned matrix giving the tresholds to apply for each quadrant for each combination of low and high frequencies
//output: - L1 norm: the results of ||Lambda \otimes x||_1
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
    


    //soft tresholding
    L1norm=wavL1norm(in_r,in_i,wavDim,weights,order);
    plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
    normPtr=mxGetPr(plhs[0]);
    normPtr[0]=L1norm;;
    
}




///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
#define WDIM(hl,num) hl*weightStep[num]
#define SUBIND(hl,num) hl * wavDimSteps[order[num]]*dataDim[num]
#define COPYCOND(el,num,hl) ((el>(dataDim[num]-HLdim[num]))&&(hl==1)) //from this condition comes the fact that wavelet weigths are 0 on the high frequency border
double wavL1norm(double* in_r,double* in_i,int* wavDim,double* weights,int* order)
{
    
    double L1norm=0;
    double* L1norm_temp;
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
    
    L1norm_temp = (double*)malloc(dataDim[4] * sizeof(double));

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
                                
                                //mexPrintf("lambda = %f\n",lambda);
                                if (lambda)
                                    for (el0=0;el0<dataDim[0]-hl0;++el0) //0 weight is given to the last element has the wavelet is circular (the last element is a combination of the first and the last, it doesnt make sense to weight it)
                                        for (el1=0;el1<dataDim[1]-hl1;++el1)
                                            for (el2=0;el2<dataDim[2]-hl2;++el2)
                                                for (el3=0;el3<dataDim[3]-hl3;++el3)
                                                {
                                                    #pragma omp parallel for private(el4, el5,ind)
                                                    for (el4 = 0; el4 < dataDim[4]-hl4; ++el4)
                                                    {
                                                        L1norm_temp[el4]=0;
                                                        for (el5 = 0; el5 < dataDim[5]-hl5; ++el5)
                                                        {
                                                            ind=subindex+IDX_OMP(wavDimSteps);
                                                            L1norm_temp[el4]+=lambda*sqrt(in_r[ind]*in_r[ind]+in_i[ind]*in_i[ind]);



                                                        }
                                                    }
                                                    
                                                    for (el4 = 0; el4 < dataDim[4]-hl4; ++el4)
                                                        L1norm+=L1norm_temp[el4];
                                                    
                                                    
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

                                if (lambda)
                                    for (el0=0;el0<dataDim[0]-hl0;++el0)
                                        for (el1=0;el1<dataDim[1]-hl1;++el1)
                                            for (el2=0;el2<dataDim[2]-hl2;++el2)
                                                for (el3=0;el3<dataDim[3]-hl3;++el3)
                                                {
                                                    #pragma omp parallel for private(el4, el5,coeff,ind)
                                                    for (el4 = 0; el4 < dataDim[4]-hl4; ++el4)
                                                    {
                                                        L1norm_temp[el4]=0;
                                                        for (el5 = 0; el5 < dataDim[5]-hl5; ++el5)
                                                        {

                                                            ind=subindex+IDX_OMP(wavDimSteps);
                                                            L1norm_temp[el4]+=lambda*fabs(in_r[ind]);//sqrt(wav_r[ind]*wav_r[ind]);
                                                        }
                                                    }
                                                    
                                                    for (el4 = 0; el4 < dataDim[4]-hl4; ++el4)
                                                        L1norm+=L1norm_temp[el4];
                                                }
                            }
    }
    
    free(L1norm_temp);
    
    return L1norm;
    
}