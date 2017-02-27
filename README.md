# myMFCC
Mel Frequency Cepstral Coefficient computation in C++

[COPYRIGHT]

myFCC:
Copyright (c) 2017 - Ignacio Jáuregui Novo
Multimedia Technologies Group. University of Vigo, Spain
email: ignaciojauregui@gts.uvigo.es
 
This software is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

FFTW:
FFTW is copyright © 1997--1999 Massachusetts Institute of Technology.
Copyright (c) 2003, 2007-14 Matteo Frigo
Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology

The following statement of license applies *only* to this header file,
and *not* to the other files distributed with FFTW or derived therefrom:

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

[DESCRIPTION]

This C++ software computes MFCC from a vector signal.

This software uses OMP. You can deactivate OMP by commenting line #59 in /include/myFCC.h

DCT computation is not optimized, although it might not be necessary for only a few MFCC coefficients. Contributions in this aspect would be nicely received

Workflow:

1- First order pre-emphasis and Hamming windowing

2- PSD computation

3- MEL triangular filter bank application to PSD

4- LOG of MEL filtered PSD

5- DCT decorrelation of Log-MEL filtered PSD

6- LIFTERING (sinusoidal) of MFCC coefficients


[USAGE]
1 - Create myFCC object: myFCC(size_t numMel, size_t numMFCC, double* freqRange, double alpha, vector<double> signal,double fs, double liftparam);

	Example:
	
	myFCC superMFCC(numFilterBanks,numMFCC,freqRange,preEmphasisCoef,signal,samplingFreq,liftCoef);	
	
	Where:
	
		numFilterBanks: 		number of MEL filter banks applied to FFT
		
		numMFCC:			number of MFFC coefficients
		
		freqRange:		        Frequency range of computation in Hertz.
		
		preEmphasisCoef:		1st. order FIR filter tap for pre-emphasis  
		
		signal:				vector of source signal
		
		samplingFreq:			Sampling frequency
		
		liftCoef:			Sinusoidal lifting coefficient 
		
2- Get MFCC from myFCC-object method getMFCC:     void  getMFCC(vector<double> signal, vector<double>& mffcCoefs)

	Example:
	
	superMFCC.getMFCC(signal,coefs);
	
	Where:
	
		signal:					vector of source signal
		
		coefs:					pointer to vector of coefficients
		
EXAMPLE INCLUDED IN /example FOLDER 

[ISNTALLATION]

This software is cross-platform.

Add /include/myFCC.h and myFCC.cpp to your project and link properly.

IMPORTANT:

THIS SOFTWARE USES FFTW FOR FFT COMPUTATION SO YOU HAVE TO INSTALL IT. http://www.fftw.org/

[ISSUES]

You can report any issue, commentary or contribution to jaureguinovo@gmail.com or ignaciojauregui@gts.uvigo.es
