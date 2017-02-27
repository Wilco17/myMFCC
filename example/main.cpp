#include <stdio.h>
#include <windows.h>
#include <include/myFCC.h>

int main()
{
	size_t 				nSignal 					    = 1536;
    vector<double> 		signal(nSignal,1.0);
    double 				freqRange[2] 					= {200,8000};
	size_t			    numFilterBanks 					= 32;
	size_t 				numMFCC 						= 12;
	double 				preEmphasisCoef 			    = 0.98;
	double 				samplingFreq 					= 48E3;
	double 				liftCoef 						= 20;
	vector<double> 		coefs(numMFCC,0.0);
	
	/** Object creation **/
    myFCC superMFCC(numFilterBanks,numMFCC,freqRange,preEmphasisCoef,signal,samplingFreq,liftCoef);
    
	/** MFCC computation **/
    DWORD start_time = GetTickCount();
    for (size_t ii=0; ii<1000;ii++)
        superMFCC.getMFCC(signal,coefs);
    DWORD end_time = GetTickCount();
	/** Output **/
	cout << "Copyrigth (c) 2017 -- Ingacio Jauregui Novo. University of Vigo, SPAIN"<<endl;
	cout << "Software under GPL"<<endl;
	cout<<"--------------------------------------------------------------------------"<<endl;
	cout <<"FFTW is copyright (c) 1997--1999 Massachusetts Institute of Technology" << endl;
    cout <<"Copyright (c) 2003, 2007-14 Matteo Frigo"<<endl;
    cout <<"Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology"<<endl;
	cout<<"--------------------------------------------------------------------------"<<endl;
    cout <<"MFCC time: " <<(double)(end_time - start_time)/1000 << " miliseconds"<<endl;
	for (size_t i=0;i<numMFCC;i++)
		cout <<"Coef "<<i<<": "<<coefs[i]<<endl;
	
	
	return 0;
}