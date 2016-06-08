#include <string>    // string
#include <vector>    // vector
#include <cmath>     // pow
#include <fstream>   // fstream, ifstream, ofstream
#include <iostream>  // cout, endl
#include <fftw3.h>   // FFTW3
#include <algorithm> // swap
#include <cstdlib>   // EXIT_FAILURE, EXIT_SUCCESS
#include <chrono>    // time
#include <random>    // random numbers

#include <stdlib.h>  /* srand, rand */
#include <time.h>    /* time */

using namespace std;

//--------------------------------------------------------------------
// WRITE TO BINARY FILE
//--------------------------------------------------------------------

void binWrite(vector<int8_t> &data, string filePath) {
    clock_t begin = clock();
    cout << "[BINWRITE] Writing to binary file... " << flush;

    ofstream binData;
    binData.open(filePath, ios::out | ios::binary);
    size_t dataLength = data.size();
    vector<int8_t> dataINT8(dataLength);
    for (size_t i{0}; i < dataLength; ++i) {
        binData.write(reinterpret_cast<const char *>(&data[i]), 1);
    }
    binData.close();
    clock_t end = clock();
    cout << "Done (" << double(end - begin) / CLOCKS_PER_SEC << "s)" << endl;
}

//--------------------------------------------------------------------
// SIMULATION DATA
//--------------------------------------------------------------------

void simulateRealData(vector<double> &simData, vector<double> pulseProfile, double carrierFreq, double step, double timeSamples) {
    clock_t begin = clock();
    cout << "[SIMULATE] Simulating data... " << flush; 

    /* initialize random seed: */
    srand (time(NULL));

    // construct a trivial random generator engine from a time-based seed:
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    normal_distribution<double> distribution (0.0, 1.0);

    size_t x{0};
    int timeSeg = 32000000;
    double carrierPeriod = 1 / carrierFreq;
    for (size_t j{}; j < timeSamples/timeSeg; ++j) {
      	for (size_t i{0}; i < timeSeg; ++i) {
            simData[i] = 127 * distribution(generator);
            // simple cosine wave
      	    //simData[i] = 127 * cos(2 * M_PI * carrierFreq * step * i);
            simData[i] = floor(simData[i] + 0.5);
      	}
    	  clock_t end = clock();
    	  cout << "Done (" << double(end - begin) / CLOCKS_PER_SEC << "s)" << endl;
  	}
}

//--------------------------------------------------------------------
// MODULATE SIGNAL
//--------------------------------------------------------------------

void modulateData(vector<double> &simData, vector<double> pulseProfile, double timeSamples) {
    clock_t begin = clock();
    cout << "[MODULATE] Modulating data... " << flush; 

    size_t x{0};
    for (size_t i{0}; i < timeSamples; ++i) {
        if (x == pulseProfile.size()) {
            x = 0; // pulsar is off
        }  
        // simulate signal - pulsar modulated normal distribution
        simData[i] *= pulseProfile[x];  
        simData[i] = floor(simData[i] + 0.5);
        ++x;
    }

    clock_t end = clock();
    cout << "Done (" << double(end - begin) / CLOCKS_PER_SEC << "s)" << endl;
}

//--------------------------------------------------------------------
// SQUARE PULSE PROFILE
//--------------------------------------------------------------------

vector<double> getSquareProfile(vector<double> &profile, double dutyCycle) {
  clock_t begin = clock();
  cout << "[SQUARE PROFILE] Doing square profile... " << flush;

  size_t someLength = profile.size() * (1 - (dutyCycle / 100));
  for (size_t i{0}; i < someLength; ++i) {
      profile[i] = 0.0;
  }
  for (size_t i{someLength}; i < profile.size(); ++i) {
      profile[i] = 1.0;
  }

  clock_t end = clock();
  cout << "Done (" << double(end - begin) / CLOCKS_PER_SEC << "s)" << endl;
  return profile;
}

//--------------------------------------------------------------------
// GAUSSIAN PULSE PROFILE
//--------------------------------------------------------------------

void getGaussProfile(vector<double> &profile, double FWHM, double step) {
    clock_t begin = clock();
    cout << "[GAUSS PROFILE] Doing Gauss profile... " << flush;

    // number of samples in FWHM
    double fwhmSamples = FWHM * 1000000 / step; 
    //cout << endl << "FWHM samples = " << fwhmSamples << endl;
    double sigma = fwhmSamples / (2. * sqrt(2. * log(2.)));
    //cout << "sigma = " << sigma << endl;

    int middlePoint = profile.size() / 2;
    //cout << "middle point = " << middlePoint << endl;
    for (int i{0 - middlePoint}; i <= middlePoint; ++i) {
        profile[i + middlePoint] = exp(-1 * (pow (i, 2.) / (2. * pow (sigma, 2.))));
    }

    clock_t end = clock();
    cout << "Done (" << double(end - begin) / CLOCKS_PER_SEC << "s)" << endl;
}

//--------------------------------------------------------------------
// HILBERT TRANSFORM
//--------------------------------------------------------------------

void doHilbertTransform(vector<double> &dataImagParts, vector<double> data) { 
    
    clock_t begin = clock();
    //cout << "[HILBERT] Doing Hilbert transform... " << endl;
    size_t dataLength = data.size();

    // do FFT
    cout << "[HILBERT] Doing FFT... " << flush;
    fftw_complex *in, *out, *in2;
    fftw_plan plan_forward;

    in  = (fftw_complex*) fftw_malloc(dataLength * sizeof(fftw_complex)); // for input
    in2 = (fftw_complex*) fftw_malloc(dataLength * sizeof(fftw_complex)); // for output
    out = (fftw_complex*) fftw_malloc(dataLength * sizeof(fftw_complex)); // for comparison

    for (size_t i{0}; i < dataLength; ++i) {
        in[i][0] = data[i];
        in[i][1] = 0.;
    }   

    plan_forward = fftw_plan_dft_1d(dataLength, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_forward);  
    fftw_destroy_plan(plan_forward);
    clock_t end = clock();
    cout << "Done (" << double(end - begin) / CLOCKS_PER_SEC << "s)" << endl; 

    begin = clock();
    cout << "[HILBERT] Doing Hilbert transform... " << flush;  

    // do phase shift
    double l2 = (dataLength + 1) / 2;
    for (size_t i = 1; i < l2; i++) {
        out[i][0] *= 2.;
        out[i][1] *= 2.;
    }
    l2 = dataLength / 2 + 1;
    for (size_t i = l2; i < dataLength; i++) {
        out[i][0] = 0.;
        out[i][1] = 0.;
    }
    end = clock();
    cout << "Done (" << double(end - begin) / CLOCKS_PER_SEC << "s)" << endl; 
      
    // do inverse FFT
    begin = clock();
    cout << "[HILBERT] Doing inverse FFT... " << flush;  
    fftw_plan pb = fftw_plan_dft_1d(data.size(), out, in2, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(pb);
    fftw_destroy_plan(pb);
        
    for (size_t i{0}; i < dataLength; ++i) {    
        dataImagParts[i] = floor((in2[i][1]) / dataLength + 0.5);
    }
        
    free (in);
    free (in2);
    fftw_free (out);
    end = clock();
    cout << "Done (" << double(end - begin) / CLOCKS_PER_SEC << "s)" << endl; 
}  

//--------------------------------------------------------------------
// MAIN PROGRAM
//--------------------------------------------------------------------

int main() { 
    // define parameters
    // all time variables multiplied by 2, because Hilbert transform
    // returns half data
    double pulsePeriod = 0.3587384107696 * 2; // seconds
    double dutyCycle = 10; // only needed for square-shape pulsar pulse
    double step = 0.0625; // useconds
    double FWHM = 0.006 * 2; // seconds
    double duration = 1 * 2; // seconds // or 10 * 2 for full dataset
    double carrierFreq = 4; // MHz // only needed for monochromatic wave case
    

    // TEST Hilbert with cosine and sine
    
    /*
    double pulsePeriod = 0.02; // seconds
    double dutyCycle = 10; 
    double step = 10; // useconds
    double FWHM = 0.005; // seconds
    double duration = 1; // seconds
    double carrierFreq = 0.002; // MHz
    */

    /*
    double pulsePeriod = 0.02; // seconds
    double dutyCycle = 10; 
    double step = 10; // useconds
    double FWHM = 0.005; // seconds
    double duration = 1; // seconds
    double carrierFreq = 0.000002; // MHz
    */
    
    // number of samples in pulse period
    double pulsePeriodSamples = pulsePeriod * 1000000 / step; 
    vector<double> pulseProfile(pulsePeriodSamples);

    // produce pulsar profile
    // getSquareProfile(pulseProfile, dutyCycle);
    getGaussProfile(pulseProfile, FWHM, step);

    // convert duration to number of samples
    size_t timeSamples = duration * 1000000 / step; 

    // simulate real data components
    vector<double> simDataRealR(timeSamples); // polarisation 1
    vector<double> simDataRealL(timeSamples); // polarisation 2 
    simulateRealData(simDataRealR, pulseProfile, carrierFreq, step, timeSamples);
    simulateRealData(simDataRealL, pulseProfile, carrierFreq, step, timeSamples);
    cout << "[MAIN] Number of data points simulated = " << simDataRealR.size() << endl;
    
    // do Hilbert transform - produce imaginary data components
    cout << "[MAIN] Getting imaginary parts for first polarisation..." << endl;
    vector<double> simDataImagR(simDataRealR.size()); // polarisation 1
    doHilbertTransform(simDataImagR, simDataRealR);
    cout << "[MAIN] Getting imaginary parts for second polarisation..." << endl;
    vector<double> simDataImagL(simDataRealL.size());  // polarisation 2
    doHilbertTransform(simDataImagL, simDataRealL);

    // print first 20 real data points
    size_t N{20};
    cout << "[MAIN] First " << N << " real parts:" << endl;
    for (size_t i{0}; i < N; ++i) {
        cout << simDataRealR[i] << endl;
    }
    
    // print firs 20 imaginary data points
    cout << "[MAIN] First " << 20 << " imaginary parts:" << endl;
    for (size_t i{0}; i < 20; ++i) {
        cout << simDataImagR[i] << endl;
    }

    // Modulate data with pulse profile
    modulateData(simDataRealR, pulseProfile, timeSamples);
    modulateData(simDataImagR, pulseProfile, timeSamples);
    modulateData(simDataRealL, pulseProfile, timeSamples);
    modulateData(simDataImagL, pulseProfile, timeSamples);

    // data vector will contain right polarisation real and imag parts as well as
    // left polarisation real and imag parts, hence length * 4
    size_t dataLength = (simDataRealR.size()) * 4;

    /*
    // combine both polarisation real and imaginary parts to one double vector
    vector<double> simData(dataLength);
    for (size_t i{0}; i < dataLength; i += 4) {
        simData[i + 0] = simDataRealR[i / 4];
        simData[i + 1] = simDataImagR[i / 4];
        simData[i + 2] = simDataRealL[i / 4];
        simData[i + 3] = simDataImagL[i / 4];
    }
    */
    
    // print real and imaginary parts to separate files for testing
    /*
    cout << "Printing right polarisation real parts..." << endl;
    ofstream fileRealR("simRealR.txt");
    for (size_t i{0}; i < simDataRealR.size(); ++i) {
        fileRealR << simDataRealR[i] << endl;
    } 
    fileRealR.close();

    cout << "Printing right polarisation imaginary parts..." << endl;
    ofstream fileImagR("simImagR.txt");
    for (size_t i{0}; i < simDataImagR.size(); ++i) {
        fileImagR << simDataImagR[i] << endl;
    } 
    fileImagR.close();
    */

    /*
    cout << "Printing complex numbers..." << endl;
    ofstream fileComplex("simComplexR.txt");
    for (size_t i{0}; i < simDataImagR.size(); ++i) {
        if (simDataImagR[i] < 0) {
            fileComplex << simDataRealR[i] << simDataImagR[i] << "j" << endl;
        }
        else {
            fileComplex << simDataRealR[i] << "+" << abs(simDataImagR[i]) << "j" << endl;
        }
    } 
    fileComplex.close();
    */
    
    // print magnitudes for testing...
    cout << "Printing magnitudes..." << endl;
    ofstream fileMagsR("simMagsR.txt");
    for (size_t i{0}; i < simDataImagR.size(); ++i)     {
        fileMagsR << sqrt(pow(simDataRealR[i], 2.) + pow(simDataImagR[i], 2.)) << endl;
    } 
    fileMagsR.close();
    
    vector<int8_t> simData(dataLength);
    for (size_t i{0}; i < dataLength; i += 4) {
        simData[i + 0] = static_cast<int8_t>(simDataRealR[i / 4]);
        simData[i + 1] = static_cast<int8_t>(simDataImagR[i / 4]);
        simData[i + 2] = static_cast<int8_t>(simDataRealL[i / 4]); 
        simData[i + 3] = static_cast<int8_t>(simDataImagL[i / 4]);
    }

    cout << "Data elements in one channel: " << simDataRealR.size() << endl;

    // write data to binary file
    binWrite(simData, "random.dada");
  
    return EXIT_SUCCESS;
}      

 
