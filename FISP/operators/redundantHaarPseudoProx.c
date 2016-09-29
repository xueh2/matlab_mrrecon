#include <string.h>
#include <omp.h>
#include "mex.h" 

#define DIMS 6
#define IDX(step) el[0]*step[0]+el[1]*step[1]+el[2]*step[2]+el[3]*step[3]+el[4]*step[4]+el[5]*step[5]
#define IDX_OMP(step) el0 * step[order[0]] + el1 * step[order[1]] + el2 * step[order[2]] + el3 * step[order[3]] + el4 * step[order[4]] + el5 * step[order[5]]
#define SUMDIM(x) x[0]+x[1]+x[2]+x[3]+x[4]+x[5]
#define PRODDIM(x) x[0]*x[1]*x[2]*x[3]*x[4]*x[5]
#define changeOrder(o1,o2,o3,o4,o5,o6) order[0]=o1;order[1]=o2;order[2]=o3;order[3]=o4;order[4]=o5;order[5]=o6;

double wavSoftTreshold(double* wav_r,double* wav_i,int* wavDim,double* weights,double step,int* order);


//to be mexed by compile_and_test_redundantHaar6D.m
//the use of the function is explained there
//inputs: - x a matrix up to 6D
//        - lambdas: a [2][2][2][2][2][2] dimensioned matrix giving the tresholds to apply for each quadrant for each combination of low and high frequencies
//output: - y, the result of the softtresholding of x in the wavelet dimension, and reprojected to image space. W'*softtreshold(W*x,lambdas)
//        - L1norm: the weighted wavelet L1 norm of softtreshold(W*x,lambdas), be careful, it is slightly different from the the weighted L1 norm of W*y as W*W' is not identity
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
    
    //wavelet pointer
    double* wav_r;
    double* wav_i;

    //input/outputdimensions
    mwSize numDim;
    mwSize numWeights;
    const mwSize* inDimMatlab;
    const mwSize* weightsDimMatlab;
    int wavDim[DIMS];           //dimension of the wavelet data
    int inDim[DIMS];           //current dimensions (is changed at each wavelet computation)
    int order[DIMS]={5, 2, 3, 4, 0, 1};
    
    //loop
    int k;
    
    //data
    int complex;

    
    //use all the processors if possible
    omp_set_num_threads(omp_get_num_procs());

    //initialization of variables
    for (k=0;k<DIMS;++k)
    {
        wavDim[k]=1;
        inDim[k]=1;
    }
    
    
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
    inDimMatlab=mxGetDimensions (prhs[0]);
    
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
        wavDim[k]=(inDimMatlab[k]>1? 2*inDimMatlab[k]:1);
        inDim[k]=inDimMatlab[k];
    }
    
    //read input data
    in_r = mxGetPr(prhs[0]);
    in_i = (complex? mxGetPi(prhs[0]):0);
    weights = mxGetPr(prhs[1]);
        
    //create output data
    plhs[0] = (mxArray *) mxCreateNumericArray (numDim, inDimMatlab, mxGetClassID (prhs[0]), complex);
    out_r=mxGetPr(plhs[0]);
    out_i=(complex? mxGetPi(plhs[0]):0);
    

    //allocation of the data
    wav_r = (double*)malloc(PRODDIM(wavDim) * sizeof(double));
    if (complex)
        wav_i = (double*)malloc(PRODDIM(wavDim) * sizeof(double));
    else
        wav_i=0;

    //forward
    computeAllForwardWavelets(in_r,wav_r,inDim,wavDim);
    if (complex)
        computeAllForwardWavelets(in_i,wav_i,inDim,wavDim);

    //soft tresholding
    L1norm=wavSoftTreshold(wav_r,wav_i,wavDim,weights,step,order);
    if (nlhs==2)
    {
        plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
        normPtr=mxGetPr(plhs[1]);
        normPtr[0]=L1norm;;
    }

    //backward
    computeAllBackwardWavelets(out_r,wav_r,inDim,wavDim);
    if (complex)
        computeAllBackwardWavelets(out_i,wav_i,inDim,wavDim);
    
    //free memory
    free(wav_r);
    if (complex)
        free(wav_i);
    
}

int computeAllForwardWavelets(double* data,double* wav,const int dataDim[DIMS],int wavDim[DIMS])
{
    //order of computation to copy the data (and later for wavelets)
    //from outer to inner loop
    //usual data are [x][y][coil][z][t][repetition]
    //largest dimension are preferably in the innermost loops
    int order[DIMS]={5, 2, 3, 4, 1, 0};
    
    
    int isFirst=1;          //if it is the first step then we need to copy the data 
    int loopDim[DIMS];      //dimension over which to lool (evolves with the successive computations)
    int k;                  //loop
    for (k=0;k<DIMS;k++)
        loopDim[k]=dataDim[k];
    
    //copy the data to the wavelet output
    //copyData(data,wav,dataDim,wavDim,order);

    //run the wavelet computation for each existing dimension
    if (dataDim[0]>1)
    {
        //order needs to be optimized
        //the loop before the last is parrallelized
        changeOrder(5, 2, 3, 4, 1, 0) 
        copyAndComputeForwardWavelets(data,wav,dataDim,wavDim,loopDim,order);
        isFirst=0;
        loopDim[0]=loopDim[0]*2; //twice as many data after we compute wavelets along one dimension
    }
    if (dataDim[1]>1)
    {
        changeOrder(5, 2, 3, 4, 0, 1)
        if (isFirst) //either copy the data and compute wavelets
            copyAndComputeForwardWavelets(data,wav,dataDim,wavDim,loopDim,order);
        else        //or just compute wavelets (faster)
            computeForwardWavelets(wav,wavDim,loopDim,order);
        isFirst=0;
        loopDim[1]=loopDim[1]*2;
    }
    if (dataDim[2]>1)
    {
        changeOrder(5, 3, 4, 1, 0, 2)
        if (isFirst)
            copyAndComputeForwardWavelets(data,wav,dataDim,wavDim,loopDim,order);
        else
            computeForwardWavelets(wav,wavDim,loopDim,order);
        isFirst=0;
        loopDim[2]=loopDim[2]*2;
    }
    if (dataDim[3]>1)
    {
        changeOrder(5, 2, 4, 1, 0, 3)
        if (isFirst)
            copyAndComputeForwardWavelets(data,wav,dataDim,wavDim,loopDim,order);
        else
            computeForwardWavelets(wav,wavDim,loopDim,order);
        isFirst=0;
        loopDim[3]=loopDim[3]*2;
    }
    if (dataDim[4]>1)
    {
        changeOrder(5, 2, 3, 1, 0, 4)
        if (isFirst)
            copyAndComputeForwardWavelets(data,wav,dataDim,wavDim,loopDim,order);
        else
            computeForwardWavelets(wav,wavDim,loopDim,order);
        isFirst=0;
        loopDim[4]=loopDim[4]*2;
    }
    if (dataDim[5]>1)
    {
        changeOrder(2, 3, 4, 1, 0, 5)
        if (isFirst)
            copyAndComputeForwardWavelets(data,wav,dataDim,wavDim,loopDim,order);
        else
            computeForwardWavelets(wav,wavDim,loopDim,order);
        isFirst=0;
        loopDim[5]=loopDim[5]*2;
    }
    
    return 1;
}












int copyAndComputeForwardWavelets(double* inWav,double* outWav, const int inWavDim[DIMS],const int outWavDim[DIMS],const int loopDim[DIMS],const int order[DIMS])
{
   
    //prepare dimensions
    int inDimStep[DIMS];  //step to go to next element for each dimension
    int outDimStep[DIMS];  //step to go to next element for each dimension
    int dataStep;       //size to swith from mean wavelets to difference wavelets

    int inCurDimStep,outCurDimStep;     //step to increment for the current computed dimension
    int loopDim0, loopDim1, loopDim2, loopDim3, loopDim4, loopDim5; //better
    int el0, el1, el2, el3, el4,el5;        //positions
    
    int* inInd;
    int* outInd;
    

    double* temp;        //to keep track of the 1st element (the wavelet being looped, it is needed for the last element computation)
    
    int k;              //loop variable
    
    //mexPrintf("Called for dimension %i.\n",order[DIMS-1]);

    //populate the variables (steps and length of dimensions)
    inDimStep[0]=1;
    outDimStep[0]=1;
    for (k=1;k<DIMS;++k)
    {
        inDimStep[k]   = inDimStep[k-1]*inWavDim[k-1];
        outDimStep[k]   = outDimStep[k-1]*outWavDim[k-1];
    }
    dataStep=outDimStep[order[DIMS-1]]*loopDim[order[DIMS-1]];
    inCurDimStep=inDimStep[order[DIMS-1]];
    outCurDimStep=outDimStep[order[DIMS-1]];
    
    //save loop size in separate variables (better optimized)
    loopDim0=loopDim[order[0]]; loopDim1=loopDim[order[1]]; loopDim2=loopDim[order[2]]; loopDim3=loopDim[order[3]]; loopDim4=loopDim[order[4]]; loopDim5=loopDim[order[5]];
    inInd = (int*)malloc(loopDim4 * sizeof(int));
    outInd = (int*)malloc(loopDim4 * sizeof(int));
    temp = (double*)malloc(loopDim4 * sizeof(double));
    
    //Real values
    for (el0=0;el0<loopDim0;++el0)
        for (el1=0;el1<loopDim1;++el1)
            for (el2=0;el2<loopDim2;++el2)
                for (el3=0;el3<loopDim3;++el3)
                {
                    #pragma omp parallel for private(el4, el5)
                    for (el4=0;el4<loopDim4;++el4)
                    {
                        //record the first element
                        el5=0;
                        inInd[el4]=IDX_OMP(inDimStep);
                        if (inWav==outWav)
                            outInd[el4]=inInd[el4];
                        else
                            outInd[el4]=IDX_OMP(outDimStep);

                        
                        temp[el4]=inWav[inInd[el4]];
                        
                        for (el5=0;el5<loopDim5-1;++el5)
                        {
                            inInd[el4]=IDX_OMP(inDimStep);
                            if (inWav==outWav)
                                outInd[el4]=inInd[el4];
                            else
                                outInd[el4]=IDX_OMP(outDimStep);
                            
                            //mexPrintf("el index is %i %i %i %i %i %i .\n",el[0]*dimStep[0],el[1]*dimStep[1],el[2]*dimStep[2],el[3]*dimStep[3],el[4]*dimStep[4],el[5]*dimStep[5]);
                            //high frequency
                            outWav[outInd[el4]+dataStep] = (inWav[inInd[el4]]-inWav[inInd[el4]+inCurDimStep])/2;
                            //low frequency
                            outWav[outInd[el4]]          = inWav[inInd[el4]]-outWav[outInd[el4]+dataStep];
                            

                        }
                        
                        inInd[el4]=IDX_OMP(inDimStep);
                        if (inWav==outWav)
                            outInd[el4]=inInd[el4];
                        else
                            outInd[el4]=IDX_OMP(outDimStep);
                        
                        //loop with first element to compute last coefficients
                        //ind=IDX_OMP(inDdimStep);//WARNING: we use here the value of el[order[5]] AFTER the loop
                        outWav[outInd[el4]+dataStep]  = (inWav[inInd[el4]]-temp[el4])/2;
                        outWav[outInd[el4]]           = inWav[inInd[el4]]-outWav[outInd[el4]+dataStep];
                    }
                }
    
    free(inInd);
    free(outInd);
    free(temp);
    return 1;
}



int computeForwardWavelets(double* wav,const int wavDim[DIMS],const int loopDim[DIMS],const int order[DIMS])
{
   
    //prepare dimensions
    int dimStep[DIMS];  //step to go to next element for each dimension
    int dataStep;       //size to swith from mean wavelets to difference wavelets
    int ind;            //current index
    int curDimStep;     //step to increment for the current computed dimension
    int loopDim0, loopDim1, loopDim2, loopDim3, loopDim4, loopDim5;
    int el0, el1, el2, el3, el4, el5;        //positions

    double temp;        //to keep track of the 1st element (the wavelet being looped, it is needed for the last element computation)
    
    int k;              //loop variable
    
    //mexPrintf("Called for dimension %i.\n",order[DIMS-1]);

    //populate the variables (steps and length of dimensions)
    dimStep[0]=1;
    for (k=1;k<DIMS;++k)
        dimStep[k]   = dimStep[k-1]*wavDim[k-1];
    dataStep=dimStep[order[DIMS-1]]*loopDim[order[DIMS-1]];
    curDimStep=dimStep[order[DIMS-1]];
    loopDim0=loopDim[order[0]]; loopDim1=loopDim[order[1]]; loopDim2=loopDim[order[2]]; loopDim3=loopDim[order[3]]; loopDim4=loopDim[order[4]]; loopDim5=loopDim[order[5]];

    //Real values
    for (el0=0;el0<loopDim0;++el0)
        for (el1=0;el1<loopDim1;++el1)
            for (el2=0;el2<loopDim2;++el2)
                for (el3=0;el3<loopDim3;++el3)
                {
                    #pragma omp parallel for private(el4, el5, ind, temp)
                    for (el4=0;el4<loopDim4;++el4)
                    {
                        //record the first element
                        el5=0;
                        ind=IDX_OMP(dimStep);

                        
                        temp=wav[ind];
                        
                        for (el5=0;el5<loopDim5-1;++el5)
                        {
                            //current possition
                            ind=IDX_OMP(dimStep);
                            
                            //mexPrintf("el index is %i %i %i %i %i %i .\n",el[0]*dimStep[0],el[1]*dimStep[1],el[2]*dimStep[2],el[3]*dimStep[3],el[4]*dimStep[4],el[5]*dimStep[5]);
                            //high frequency
                            wav[ind+dataStep] = (wav[ind]-wav[ind+curDimStep])/2;
                            //low frequency
                            wav[ind]          = wav[ind]-wav[ind+dataStep];

                        }
                        
                        //loop with first element to compute last coefficients
                        ind=IDX_OMP(dimStep);//WARNING: we use here the value of el[order[5]] AFTER the loop
                        wav[ind+dataStep]  = (wav[ind]-temp)/2;
                        wav[ind]           = wav[ind]-wav[ind+dataStep];
                    }
                }
    return 1;
}





//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


int computeAllBackwardWavelets(double* data,const double* wav,int dataDim[DIMS],int wavDim[DIMS])
{
    //order of computation to copy the data (and later for wavelets)
    //from outer to inner loop
    //usual data are [x][y][coil][z][t][repetition]
    //largest dimension are preferably in the innermost loops
    int order[DIMS]={5, 2, 3, 4, 1, 0};
    
    int current=0;          //current dimension being processed
    int loopDim[DIMS];      //dimension over which to lool (evolves with the successive computations)
    int k;                  //loop
    int wavSteps=0;         //how many backward steps to perform (last needs to copy the data)
    int computedSteps=1;    //steps that have already been computed
    //pointers to the correct structures for computation, and to temp data
    double* in;
    double* out;
    double* temp;  
    
    int* inDim;
    int* outDim;
    int* tempDim;
    
    for (k=0;k<DIMS;k++)
    {
        if (wavDim[k]>1)
            ++wavSteps;
        loopDim[k]=wavDim[k];
    }
    
    //mexPrintf("number of steps: %i\n",wavSteps);
    //allocate temp memory
    temp = wav;
    tempDim=wavDim;
    
    //copy the data to the wavelet output
    //copyData(data,wav,dataDim,wavDim,order);

    //run the wavelet computation for each existing dimension
    in=wav;
    out=temp;
    inDim=wavDim;
    outDim=tempDim;
    
    current=0;
    if (dataDim[current]>1)
    {
        //mexPrintf("step number: %i\n",computedSteps);
        //order needs to be optimized
        //the loop before the last is parrallelized
        changeOrder(5, 2, 3, 4, 1, 0) 
        
        //reduce the dimension
        loopDim[current]=loopDim[current]/2;

        
        //if last loop then save into data
        if (computedSteps==wavSteps)
        {
            out=data;
            outDim=dataDim;
        }
        copyAndComputeBackwardWavelets(in,out,inDim,outDim,loopDim,order);
        
        //increment to next number of computed steps
        ++computedSteps;
        in=temp;
        inDim=tempDim;

    }
    
    current=1;
    if (dataDim[current]>1)
    {
        //mexPrintf("step number: %i out of %i.\n",computedSteps,wavSteps);
        changeOrder(5, 2, 3, 4, 0, 1)
        
        //reduce the dimension
        loopDim[current]=loopDim[current]/2;

        
        //if last loop then save into data
        if (computedSteps==wavSteps)
        {
            out=data;
            outDim=dataDim;
        }
        
        copyAndComputeBackwardWavelets(in,out,inDim,outDim,loopDim,order);
        
        //increment to next number of computed steps
        ++computedSteps;
        in=temp;
        inDim=tempDim;

    }
    
    current=2;
    if (dataDim[current]>1)
    {
        //mexPrintf("step number: %i\n",computedSteps);
        changeOrder(5, 3, 4, 1, 0, 2)
        
        //reduce the dimension
        loopDim[current]=loopDim[current]/2;

        
        //if last loop then save into data
        if (computedSteps==wavSteps)
        {
            out=data;
            outDim=dataDim;
        }
        
        
        copyAndComputeBackwardWavelets(in,out,inDim,outDim,loopDim,order);
        
        //increment to next number of computed steps
        ++computedSteps;
        in=temp;
        inDim=tempDim;

    }
    
    current=3;
    if (dataDim[current]>1)
    {
        //mexPrintf("step number: %i\n",computedSteps);
        changeOrder(5, 2, 4, 1, 0, 3)
        
        //reduce the dimension
        loopDim[current]=loopDim[current]/2;

        
        //if last loop then save into data
        if (computedSteps==wavSteps)
        {
            out=data;
            outDim=dataDim;
        }
        copyAndComputeBackwardWavelets(in,out,inDim,outDim,loopDim,order);
        
        //increment to next number of computed steps
        ++computedSteps;
        in=temp;
        inDim=tempDim;

    }
    
    current=4;
    if (dataDim[current]>1)
    {
        //mexPrintf("step number: %i\n",computedSteps);
        changeOrder(5, 2, 3, 1, 0, 4)
        
        //reduce the dimension
        loopDim[current]=loopDim[current]/2;

        
        //if last loop then save into data
        if (computedSteps==wavSteps)
        {
            out=data;
            outDim=dataDim;
        }
        copyAndComputeBackwardWavelets(in,out,inDim,outDim,loopDim,order);
        
        //increment to next number of computed steps
        ++computedSteps;
        in=temp;
        inDim=tempDim;

    }
    
    current=5;
    if (dataDim[current]>1)
    {
        //mexPrintf("step number: %i\n",computedSteps);
        changeOrder(2, 3, 4, 1, 0, 5)
        
        //reduce the dimension
        loopDim[current]=loopDim[current]/2;

        
        //if last loop then save into data
        if (computedSteps==wavSteps)
        {
            out=data;
            outDim=dataDim;
        }
        copyAndComputeBackwardWavelets(in,out,inDim,outDim,loopDim,order);
        
        //increment to next number of computed steps
        ++computedSteps;
        in=temp;
        inDim=tempDim;

    }
    

    return 1;
}



int copyAndComputeBackwardWavelets(double* inWav,double* outWav, const int inWavDim[DIMS],const int outWavDim[DIMS],const int loopDim[DIMS],const int order[DIMS])
{
   
    //prepare dimensions
    int inDimStep[DIMS];  //step to go to next element for each dimension
    int outDimStep[DIMS];  //step to go to next element for each dimension
    int dataStep;       //size to swith from mean wavelets to difference wavelets

    int inCurDimStep,outCurDimStep;     //step to increment for the current computed dimension
    int loopDim0, loopDim1, loopDim2, loopDim3, loopDim4, loopDim5; //better
    int el0, el1, el2, el3, el4,el5;        //positions
    
    int* inInd;
    int* outInd;
    

    double* temp;        //to keep track of the 1st element (the wavelet being looped, it is needed for the last element computation)
    
    int k;              //loop variable
    
    //mexPrintf("Called for dimension %i.\n",order[DIMS-1]);

    //populate the variables (steps and length of dimensions)
    inDimStep[0]=1;
    outDimStep[0]=1;
    for (k=1;k<DIMS;++k)
    {
        inDimStep[k]   = inDimStep[k-1]*inWavDim[k-1];
        outDimStep[k]   = outDimStep[k-1]*outWavDim[k-1];
    }
    dataStep=inDimStep[order[DIMS-1]]*loopDim[order[DIMS-1]];
    inCurDimStep=inDimStep[order[DIMS-1]];
    outCurDimStep=outDimStep[order[DIMS-1]];
    
    //save loop size in separate variables (better optimized)
    loopDim0=loopDim[order[0]]; loopDim1=loopDim[order[1]]; loopDim2=loopDim[order[2]]; loopDim3=loopDim[order[3]]; loopDim4=loopDim[order[4]]; loopDim5=loopDim[order[5]];
    inInd = (int*)malloc(loopDim4 * sizeof(int));
    outInd = (int*)malloc(loopDim4 * sizeof(int));
    temp = (double*)malloc(loopDim4 * sizeof(double));
    
    //Real values
    for (el0=0;el0<loopDim0;++el0)
        for (el1=0;el1<loopDim1;++el1)
            for (el2=0;el2<loopDim2;++el2)
                for (el3=0;el3<loopDim3;++el3)
                {
                    #pragma omp parallel for private(el4, el5)
                    for (el4=0;el4<loopDim4;++el4)
                    {
                        //record the first element
                        el5=loopDim5-1;
                        inInd[el4]=IDX_OMP(inDimStep);
                        if (inWav==outWav)
                            outInd[el4]=inInd[el4];
                        else
                            outInd[el4]=IDX_OMP(outDimStep);

                        
                        temp[el4]=inWav[inInd[el4]];
                        
                        for (el5=loopDim5-1;el5>0;--el5)
                        {
                            inInd[el4]=IDX_OMP(inDimStep);
                            if (inWav==outWav)
                                outInd[el4]=inInd[el4];
                            else
                                outInd[el4]=IDX_OMP(outDimStep);
                            
                            //mexPrintf("el index is %i %i %i %i %i %i .\n",el[0]*dimStep[0],el[1]*dimStep[1],el[2]*dimStep[2],el[3]*dimStep[3],el[4]*dimStep[4],el[5]*dimStep[5]);

                            outWav[outInd[el4]] = (inWav[inInd[el4]]            +inWav[inInd[el4]-inCurDimStep]
                                                  +inWav[inInd[el4]+dataStep]   -inWav[inInd[el4]-inCurDimStep+dataStep])/2;


                        }
                        
                        inInd[el4]=IDX_OMP(inDimStep);
                        if (inWav==outWav)
                            outInd[el4]=inInd[el4];
                        else
                            outInd[el4]=IDX_OMP(outDimStep);
                        
                        //loop with first element to compute last coefficients
                        //ind=IDX_OMP(inDdimStep);//WARNING: we use here the value of el[order[5]] AFTER the loop
                        outWav[outInd[el4]] = (inWav[inInd[el4]]            +temp[el4]
                                              +inWav[inInd[el4]+dataStep]   -inWav[inInd[el4]+(loopDim5-1)*inCurDimStep+dataStep])/2;
                    }
                }
    
    free(inInd);
    free(outInd);
    free(temp);
    return 1;
}









//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//not used anymore
int copyData(double* in,double* out,const int inDim[DIMS],const int outDim[DIMS],const int order[DIMS])
{
    //declarations
    //int el[DIMS];           //element
    int inDimSteps[DIMS];   //step to increment element for each dimension (input data)
    int outDimSteps[DIMS];  //step to increment element for each dimension (output data)
    int k, el0, el1, el2, el3, el4, el5;
    
    //populate
    inDimSteps[0] =1;
    outDimSteps[0]=1;
    for (k=1;k<DIMS;++k)
    {
        inDimSteps[k]   = inDimSteps[k-1]   *inDim[k-1];
        outDimSteps[k]  = outDimSteps[k-1]  *outDim[k-1];
    }
    

    //copy all data for real values
    for (el0=0;el0<inDim[order[0]];++el0)
        for (el1=0;el1<inDim[order[1]];++el1)
            for (el2=0;el2<inDim[order[2]];++el2)
                for (el3=0;el3<inDim[order[3]];++el3)
                {
                    #pragma omp parallel for private(el4, el5)
                    for (el4 = 0; el4 < inDim[order[4]]; ++el4)
                    {
                        for (el5 = 0; el5 < inDim[order[5]]; ++el5)
                        	out[ IDX_OMP(outDimSteps) ] = in[ IDX_OMP(inDimSteps) ];
                    }
                }
    
    
//                        for (el[order[4]]=0;el[order[4]]<inDim[order[4]];++el[order[4]])
    return 1;
}








///////////////////////////////////////////////////////////////////////////
#define WDIM(hl,num) hl*weightStep[num]
#define SUBIND(hl,num) hl * wavDimSteps[order[num]]*dataDim[num]
double wavSoftTreshold(double* wav_r,double* wav_i,int* wavDim,double* weights,double step,int* order)
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

    if (wav_i)
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
                                if (lambda)
                                    for (el0=0;el0<dataDim[0]-hl0;++el0) //0 weight is given to the last element has the wavelet is circular (the last element is a combination of the first and the last, it doesnt make sense to weight it)
                                        for (el1=0;el1<dataDim[1]-hl1;++el1)//from this kind of "aborted" comes the fact that wavelet weigths are 0 on the high frequency border
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
                                                            coeff=sqrt(wav_r[ind]*wav_r[ind]+wav_i[ind]*wav_i[ind]);

                                                            if (coeff-tmp>0)
                                                            {
                                                                L1norm_temp[el4]+=(coeff-tmp)*lambda;
                                                                coeff=1-tmp/coeff;
                                                                wav_r[ind]=wav_r[ind]*coeff;
                                                                wav_i[ind]=wav_i[ind]*coeff;
                                                            }
                                                            else
                                                            {
                                                                wav_r[ind]=0;
                                                                wav_i[ind]=0;
                                                            }


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
                                tmp=lambda*step;
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
                                                            //mexPrintf("index is %i .\n",ind);
                                                            //mexPrintf("weight indices %i %i %i %i %i %i.\n",WDIM(el0,0),WDIM(el1,1),WDIM(el2,2),WDIM(el3,3),WDIM(el4,4),WDIM(el5,5));
                                                            //mexPrintf("tresh is %f .\n",lambda);
                                                            coeff=fabs(wav_r[ind]);//sqrt(wav_r[ind]*wav_r[ind]);

                                                            if (coeff-tmp>0)
                                                            {
                                                                wav_r[ind]=wav_r[ind]*(1-tmp/coeff);
                                                                L1norm_temp[el4]+=lambda*(coeff-tmp);
                                                            }
                                                            else
                                                                wav_r[ind]=0;
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
