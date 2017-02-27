/*  Copyright (c) 2017 - Ignacio Jáuregui Novo
 *  Multimedia Technologies Group. University of Vigo, Spain
 *  email: ignaciojauregui@gts.uvigo.es
 *
 *  This file (myFCC.cpp) is  part of free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 *  FFTW is copyright © 1997--1999 Massachusetts Institute of Technology.
 *  Copyright (c) 2003, 2007-14 Matteo Frigo
 *  Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
 *
 * The following statement of license applies *only* to this header file,
 * and *not* to the other files distributed with FFTW or derived therefrom:
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * */

#include "include/myFCC.h"

using namespace std;
myFCC::myFCC(size_t numMel, size_t numMFCC, double* freqRange, double Alpha, vector<double> signal,double fs, double liftparam)
{
      numMelComponents  =   numMel;
      numMFCCcoefs      =   numMFCC;
      signalL           =   static_cast<size_t>(signal.size());
      nFFT              =   pow(2,ceil(log2(signalL)));
      halfFFt           =   static_cast<size_t>(nFFT / 2) + 1;
      alpha             =   Alpha;
      liftParam         =   liftparam;
      Fs                =   fs;
      procSignal.resize(signalL,0.0);
      /** Resizing matrixes **/
      melFilterBank.resize(numMelComponents);
      for(size_t i=0;i<numMelComponents;i++)
          melFilterBank[i].resize(halfFFt,0.0);
      liftVector.resize(numMFCCcoefs);
      dctMatrix.resize(numMFCCcoefs);
      for(size_t i=0;i<numMFCCcoefs;i++)
      {
          dctMatrix[i].resize(numMelComponents,0.0);
          liftVector[i] = 1+0.5*liftParam*sin(M_PI*i/liftParam);
      }
      /** Mel-FilterBank matrix generation **/
      for(size_t jj=0; jj < numMelComponents;jj++)
      {
            double lowMelf  =   1127*log(1+freqRange[0]/700);
            double highMelf =   1127*log(1+freqRange[1]/700);
            double aux_c_0  =   700*exp((lowMelf+jj*(((highMelf -lowMelf)/(numMelComponents+1))))/1127)-700;
            double aux_c_1  =   700*exp((lowMelf+(jj+1)*(((highMelf -lowMelf)/(numMelComponents+1))))/1127)-700;
            double aux_c_2  =   700*exp((lowMelf+(jj+2)*(((highMelf -lowMelf)/(numMelComponents+1))))/1127)-700;
            vector<double> aux_mel(halfFFt,0.0);
            for (size_t nn = 0; nn < halfFFt; nn++)
            {
                double f      = 0.5*Fs/(halfFFt-1)*nn;
                if(f>=aux_c_0 && f<=aux_c_1)
                    aux_mel[nn] = (f-aux_c_0)/(aux_c_1-aux_c_0);
                else if (f>=aux_c_1 && f<=aux_c_2)
                    aux_mel[nn] = (aux_c_2-f)/(aux_c_2-aux_c_1);
            }
            melFilterBank[jj] = aux_mel;
      }
      /** DCT Matrix generation **/
      for (size_t ii = 0; ii < numMFCCcoefs;ii++)
      {
          vector<double> aux(numMelComponents,0.0);
          for (size_t jj = 0; jj < numMelComponents; jj++)
          {
                aux[jj]=(sqrt(2.0/static_cast<double>(numMelComponents))*cos((ii)*(jj+0.5)/static_cast<double>(numMelComponents)*M_PI));
          }
          dctMatrix[ii]=aux;
      }
};
void myFCC::getMFCC(vector<double> signal, vector<double>& mfccCoefs)
{
    fftw_plan           planForward;
    fftw_complex*       fftOut;
	vector<double>      melFilteredPSD(numMelComponents,0.0);
    fftOut =static_cast<fftw_complex*>( fftw_malloc ( sizeof ( fftw_complex ) * halfFFt ));

    /** Pre-emphasis & windowing  **/
    procSignal[0] = signal[0]*0.08;
    for(size_t i=1;i<signalL;i++)
        procSignal[i] = (signal[i]-alpha*signal[i-1])*(0.54 - 0.46*cos(2*M_PI*i/(signalL-1)));

    /** PSD **/
    planForward = fftw_plan_dft_r2c_1d ( static_cast<int>(nFFT),static_cast<double*>(&procSignal[0]), fftOut, FFTW_ESTIMATE);
    fftw_execute ( planForward );

    /** Apply MEL-filterbank to PSD **/
#ifdef MYFCC_USE_OMP
    #pragma omp parallel for
#endif
    for (int i=0;i<static_cast<int>(numMelComponents);i++)
	{
       vector<double> aux = melFilterBank[static_cast<size_t>(i)];
       for (int j=0, jlen = static_cast<int>(halfFFt); j < jlen;j++)
           melFilteredPSD[static_cast<size_t>(i)]+=(fftOut[j][0]*fftOut[j][0]+fftOut[j][1]*fftOut[j][1])*aux[static_cast<size_t>(j)];
       melFilteredPSD[static_cast<size_t>(i)] = static_cast<double>(log(melFilteredPSD[static_cast<size_t>(i)]));
    }
    fftw_free(fftOut);

    /** DCT and sinusoidal LIFTERING**/
#ifdef MYFCC_USE_OMP
    #pragma omp parallel for
#endif
    for (int i=0;i<static_cast<int>(numMFCCcoefs);i++)
	{
        double auxm =0.0;
        for (int j=0, jlen =static_cast<int>(numMelComponents);j<jlen;j++)
            auxm +=dctMatrix[static_cast<size_t>(i)][static_cast<size_t>(j)]*melFilteredPSD[static_cast<size_t>(j)]*liftVector[static_cast<size_t>(i)];
        mfccCoefs[static_cast<size_t>(i)] = auxm;
    }
};
