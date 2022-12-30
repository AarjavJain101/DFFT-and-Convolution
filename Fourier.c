/*
* Name: Aarjav Jain
* Email: Aarjavjain2674@gmail.com
*/

/*
  * (POTENTIAL) TO DO!!
  * Make an FFT for different prime factoring than just 2
*/

// NOTE: IF DATA IS PERSONALLY ADDED, COMMENT OUT FIRST FUNCTION!

#include <stdlib.h>
#include <stdio.h>
#include <math.h>  // Mainly for acos() to get PI.
#include <complex.h> // Standard Library of Complex Numbers 
#include <time.h> // For finding the time taken by my FFT algorithm

// *** LIMITER is changed to be the amount of data points ***
// *** Change FREQRES to the frequency resolution of the FFT'd data ***
// *** SIZEF + SIZEG - 1 MUST = 2^m where m is any positive integer ***
// *** If doing a convolution: CONVOLUTION == 1 else, CONVOLUTION == 0 ***
#define PI acos(-1.0)
#define LIMITER 8
#define DATAFILE "data1.txt"
#define FREQRES 1
#define SIZEF 2
#define SIZEG 1
#define CONVOLUTION 1

void fillArrayFFT(FILE *data, double complex signal[]);
double complex *findFFT(double complex signal[], int N);
void storeDataFFT(double complex signal_fft[], int N);
double complex *inverseFFT(double complex signal[], int N);
void storeDataInverseFFT(double complex *signalInverseFFT_ptr, int N);
void fillArrayInverseFFT(FILE *creal_output, FILE *cimag_output, double complex signalInverse[]);
void fillArrayConvolution(FILE *arrayA_file, FILE *arrayB_file, double complex arrayA[], double complex arrayB[]);
double complex *doConvolution(double complex arrayA[], double complex arrayB[], int sizeA, int sizeB, double complex *convolution, double complex convolution_FFT[]);
void storeDataConvolution(double complex *convolution, int N);
//double complex *findFT(double complex signal[], double complex signalFT[], int N);

int main(void) 
{
  // Do the FFT
  double complex signal[LIMITER];
  fillArrayFFT(fopen(DATAFILE, "r"), signal);
  
  double complex *signalFFT_ptr;
  clock_t FFT_start = clock();    // Start the clock
  signalFFT_ptr = findFFT(signal, LIMITER);
  clock_t FFT_end = clock();  // End the clock
  double timeSpentFFT = ((double)(FFT_end - FFT_start) / CLOCKS_PER_SEC) / 0.001;  // Tells time in milliseconds
  
  storeDataFFT(signalFFT_ptr, LIMITER);
  printf("\ntime spent for FFT = %lf milliseconds\n\n", timeSpentFFT);

  // Do the Inverse FFT
  double complex signalInverse[LIMITER];
  fillArrayInverseFFT(fopen("creal_output.txt", "r"), fopen("cimag_output.txt", "r"), signalInverse);

  double complex *signalInverseFFT_ptr;
  clock_t inverseFFT_start = clock();    // Start the clock
  signalInverseFFT_ptr = inverseFFT(signalInverse, LIMITER);
  clock_t inverseFFT_end = clock();  // End the clock
  double timeSpentInverseFFT = ((double)(inverseFFT_end - inverseFFT_start) / CLOCKS_PER_SEC) / 0.001;  // Tells time in milliseconds

  storeDataInverseFFT(signalInverseFFT_ptr, LIMITER);
  printf("\ntime spent for Inverse FFT = %lf milliseconds\n\n", timeSpentInverseFFT);

  // Do convolution
  double complex arrayA[SIZEF];
  double complex arrayB[SIZEG];
  fillArrayConvolution(fopen("arrayA_data.txt", "r"), fopen("arrayB_data.txt", "r"), arrayA, arrayB);
  double complex *convolution;
  double complex convolution_FFT[SIZEF + SIZEG - 1];
  clock_t  convolution_start = clock();    // Start the clock
  convolution = doConvolution(arrayA, arrayB, SIZEF, SIZEG, convolution, convolution_FFT);
  clock_t  convolution_end = clock();    // End the clock
  double timeSpentConvolution = ((double)(convolution_end - convolution_start) / CLOCKS_PER_SEC) / 0.001;  // Tells time in milliseconds

  storeDataConvolution(convolution, SIZEF + SIZEG - 1); 
  printf("\ntime spent for Convolution = %lf milliseconds\n\n", timeSpentConvolution);

  return 0;
}

void fillArrayFFT(FILE *data, double complex signal[]) // Take data and store into signal
{
  double datas[LIMITER];
  int index = 0;
  while (fscanf(data, "%lf", &datas[index]) == 1) // data points as doubles inside datas[]
  {
    index++;
  }
  fclose(data);

  for (int i = 0; i<LIMITER; i++) // converts double to double complex effectively
  {
    signal[i] = datas[i];
  }
}

double complex *findFFT(double complex signal[], int N) // FFT TIME! I use an in-place algorithm
{
  double complex signalEven[N/2];  // Declare array to store the even indicies
  double complex signalOdd[N/2];  // Declare array to store the odd indicies
  double complex *signalEvenFFT_ptr;  // Declare pointer to take return of the even indicies
  double complex *signalOddFFT_ptr;  // Declare pointer to take return of the odd indicies

  // return the signal if it cannot be broken down further
  if (N == 1)
  {
    return signal;
  }
  
  // Do the splitting (here it is even first then odd)
  for (int i = 0; i<N/2; i++) // Fill in signalEven with the even indexed data
  {
    signalEven[i] = signal[2*i];
  }
  for (int i = 0; i<N/2; i++)  // Fill in signalOdd with the odd indexed data
  {
    signalOdd[i] = signal[2*i +1];
  }

  // Do recursion to keep splitting
  signalEvenFFT_ptr = findFFT(signalEven, N/2);
  signalOddFFT_ptr = findFFT(signalOdd, N/2);
  
  // Bring together the even and odd DFT's and calculate the final.
  for (int i = 0; i<N/2; i++)
  {
    signal[i] = *(signalEvenFFT_ptr+i) + cexp((-2*PI*I)*(i)/N) * *(signalOddFFT_ptr+i);
    signal[N/2 + i] = *(signalEvenFFT_ptr+i) - cexp((-2*PI*I)*(i)/N) * *(signalOddFFT_ptr+i);
  }

  return signal;
}

void fillArrayInverseFFT(FILE *creal_output, FILE *cimag_output, double complex signalInverse[]) // Take FFT and store into signal
{
  double real[LIMITER];
  double imag[LIMITER];
  int index = 0;
  while (fscanf(creal_output, "%lf", &real[index]) == 1) // store into real
  {
    index++;
  }
  index = 0;
  fclose(creal_output);
  while (fscanf(cimag_output, "%lf", &imag[index]) == 1) // store into imag
  {
    index++;
  }
  fclose(cimag_output);

  for (int i = 0; i<LIMITER; i++) // converts double to double complex effectively
  {
    signalInverse[i] = real[i] + imag[i]*I;
  }
}

double complex *inverseFFT(double complex signal[], int N)
{
  double complex signalEven[N/2];  // Declare array to store the even indicies
  double complex signalOdd[N/2];  // Declare array to store the odd indicies
  double complex *signalEvenFFT_ptr;  // Declare pointer to take return of the even indicies
  double complex *signalOddFFT_ptr;  // Declare pointer to take return of the odd indicies

  // return the signal if it cannot be broken down further
  if (N == 1)
  {
    return signal;
  }
  
  // Do the splitting (here it is even first then odd)
  for (int i = 0; i<N/2; i++) // Fill in signalEven with the even indexed data
  {
    signalEven[i] = signal[2*i];
  }
  for (int i = 0; i<N/2; i++)  // Fill in signalOdd with the odd indexed data
  {
    signalOdd[i] = signal[2*i +1];
  }

  // Do recursion to keep splitting
  signalEvenFFT_ptr = inverseFFT(signalEven, N/2);
  signalOddFFT_ptr = inverseFFT(signalOdd, N/2);
  
  // Bring together the even and odd DFT's and calculate the final.
  for (int i = 0; i<N/2; i++)
  {
    signal[i] = *(signalEvenFFT_ptr+i) + cexp((2*PI*I)*(i)/N) * *(signalOddFFT_ptr+i);
    signal[N/2 + i] = *(signalEvenFFT_ptr+i) - cexp((2*PI*I)*(i)/N) * *(signalOddFFT_ptr+i);
  }

  if (N == LIMITER && CONVOLUTION == 1)
  {
    return signal;
  }
  else if (N == LIMITER && CONVOLUTION == 0) // Divide the Sum for each time-point by N
  {
    for (int i = 0; i<N; i++)
    {
      signal[i] = (1.0/N) * signal[i];
    }
  }
  return signal;
}

void storeDataFFT(double complex *signalFFT_ptr, int N)  // Put data in file so I can graph in Excel. The real component goes to creal_output and the imaginary component goes to cimag_output.
{
  FILE *frequency = fopen("frequency.txt", "w");
  FILE *creal_output = fopen("creal_output.txt", "w");
  FILE *cimag_output = fopen("cimag_output.txt", "w");
  for (int i = 0; i<N; i++)
  {
    fprintf(frequency, "%d\n", i);
    fprintf(creal_output, "%lf\n", creal(*(signalFFT_ptr+i)));
    fprintf(cimag_output, "%lf\n", cimag(*(signalFFT_ptr+i)));
  }
  fclose(frequency);
  fclose(creal_output);
  fclose(cimag_output);
}

void storeDataInverseFFT(double complex *signalInverseFFT_ptr, int N)  // Put data in file so I can graph in Excel.
{
  double sampleFrequency = N * FREQRES;
  FILE *inversed = fopen("inversed.txt", "w");
  FILE *radians = fopen("radians.txt", "w");
  for (int i = 0; i<N; i++)
  {
    fprintf(inversed, "%lf\n", creal(*(signalInverseFFT_ptr+i)));
    fprintf(radians, "%lf\n", i*(1/sampleFrequency));
  }
  fclose(inversed);
  fclose(radians);
}

void fillArrayConvolution(FILE *arrayA_file, FILE *arrayB_file, double complex arrayA[], double complex arrayB[]) // Store data for convolution
{
  double realA[SIZEF];
  double realB[SIZEG];
  int indexA = 0;
  while (fscanf(arrayA_file, "%lf", &realA[indexA]) == 1) // store into arrayA
  {
    indexA++;
  }
  fclose(arrayA_file);
  for (int i = 0; i<indexA; i++)
  {
    arrayA[i] = realA[i] + 0.0*I;
  }
  int indexB = 0;
  while (fscanf(arrayB_file, "%lf", &realB[indexB]) == 1) // store into arrayB
  {
    indexB++;
  }
  fclose(arrayB_file);
  for (int i = 0; i<indexB; i++)
  {
    arrayB[i] = realB[i] + 0.0*I;
  }
}

double complex *doConvolution(double complex arrayA[], double complex arrayB[], int sizeA, int sizeB, double complex *convolution, double complex convolution_FFT[]) // Returns the convolution of two arrays of doubles
{
  // Start by padding each array with enough zeroes to match final convolution size
  int N = sizeA + sizeB - 1;
  double complex paddedA[N]; 
  double complex paddedB[N]; 
  for (int i = 0; i<N; i++)
  {
    if (i < sizeA)
    {
      paddedA[i] = arrayA[i];
    }
    else
    {
      paddedA[i] = 0;
    }
  }
  for (int i = 0; i<N; i++)
  {
    if (i < sizeB)
    {
      paddedB[i] = arrayB[i];
    }
    else
    {
      paddedB[i] = 0;
    }
  }

  // Now find the FFT's of both padded arrays to prepare for multipication
  double complex *paddedA_FFT = findFFT(paddedA, N);
  double complex *paddedB_FFT = findFFT(paddedB, N);

  // Do the multipication
  for (int i = 0; i<N; i++)
  {
    convolution_FFT[i] = *(paddedA_FFT+i) * *(paddedB_FFT+i);
  }

  // Now inverse convolution_FFT to get just the convolution and return the pointer
  convolution = inverseFFT(convolution_FFT, N);
  for (int i = 0; i<N; i++)
  {
    *(convolution+i) = 1.0*(*(convolution+i))/N;
  }
  return convolution;
}

void storeDataConvolution(double complex *convolution, int N)  // Put data in file
{
  FILE *convolve = fopen("convolution.txt", "w");
  for (int i = 0; i<N; i++)
  {
    fprintf(convolve, "%lf\n", creal(*(convolution+i)));
  }
  fclose(convolve);
}


/*
double complex *findFT(double complex signal[], double complex signalFT[], int N)
{
  double complex *ptr;
  double complex zetaPowers[N];    // Make an array with what the complex multpliers will be called zetaPowers
  for (int i = 0; i<N; i++)
  {
    zetaPowers[i] = cexp(((-2*PI*I)*i)/N);
  }

  for (int i = 0; i<N; i++)
  {
    signalFT[i] = 0;  // Makes sure that the value you are adding to is 0
    for (int j = 0; j<N; j++) // Computing the sum for the Fourier Transformed data
    {
      signalFT[i] += signal[j] * zetaPowers[(i*j) % N];  // computes vector sum/amount of that frequecy present
    }
  }
  ptr = signalFT;
  return ptr;
}
*/