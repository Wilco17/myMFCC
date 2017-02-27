/*  Copyright (c) 2017 - Ignacio Jáuregui Novo
 *  Multimedia Technologies Group. University of Vigo, Spain
 *  email: ignaciojauregui@gts.uvigo.es
 *
 *  This file (myFCC.h) is  part of free software: you can redistribute it and/or modify
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
 * FFTW is copyright © 1997--1999 Massachusetts Institute of Technology.
 * Copyright (c) 2003, 2007-14 Matteo Frigo
 * Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
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
 */

#ifndef MYFCC_H
#define MYFCC_H
#include <stdlib.h>
#include <fftw3.h>
#include <vector>
#include <omp.h>
#define M_PI 3.141592653589793
/** Comment the following line if OMP is not desired **/
#define MYFCC_USE_OMP
using namespace std;
class myFCC{
private:
    size_t                           numMelComponents;							/** Number of MEL components**/
    size_t                           numMFCCcoefs;								/** Number of MFCC coefs **/
    double                           nFFT;										/** FFT points **/
    size_t                           halfFFt;									/** Length of unique part of FFT **/
    size_t                           signalL;									/** Source signal length **/
    vector<double>                   procSignal;								/** Processed signal (filtered and windowed) **/
    double                           alpha;										/** Pre-emphasis 1st. order filter tap **/
    double                           liftParam;									/** Sinusoidal lifting parameter **/
    double                           Fs;										/** Sampling frequency **/

    vector<double>                   liftVector;								/** Sinusoidal lifting vector **/
    vector<vector<double>>           melFilterBank;								/** Mel triangular filter-bank matrix	 numMelComponents-by-halfFFt **/
    vector<vector<double>>           dctMatrix;									/** DCT Matrix numMFCCcoefs-by-numMelComponents **/
public:
    void        getMFCC(vector<double> signal, vector<double>& mffcCoefs);
    myFCC(size_t numMel, size_t numMFCC, double* freqRange, double alpha, vector<double> signal,double fs, double liftparam);
};
#endif // MYFCC_H
