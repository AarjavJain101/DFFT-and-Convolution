# DFFT-and-Convolution
Welcome to my Fourier Transform Journey!

PURPOSES of program: 

-Given 2^m data points, where m is a positive integer, return the DFFT
  -Store the real, imaginary, and frequency components in respectve files for graphing

-Given 2^m data points, where m is a positive integer, return the Inverse DFFT
  -Store the inversed data and the associated time-points in respective files

-Given two data sets where sizeof(dataA) + sizeof(dataB) - 1 is 2^m data points, and m is a positive integer, return the convolution using the DFFT trick
  -Store the convolution into its respective file

GOALS and TIMELINE of the project:

1. Implement an Algorithm for the Discrete Fourier Transform
  -Use this algorithm to correctly obtain the frequencies of a given wave function (completed  on December 10th, 2022).

2. Implement the Discrete Fourier Transform but using the Fast Fourier Algorithm. Also, make sure to actually understand what is happening mathematically.
  -To my understanding, the most siginifcnt imporvemmnt in the Fast Fourier Algorithm comes   from the fact that the Discrete Fast Fourier Transform requires N^2 operations, however, the Fast Fourier Algorithm requires only N*log_2(N) operations. The reason is because we utilize the conjugaet properties of complex numbers to make use of reqpeat calculations. Essentially, instead of doing (large num)^2 calculations, You just do (large num / 2)^2 + (small num) calculations instead.
  -COMPARISON: 4096 samples and 4096 samples/1 radian, the DFT averaged 67 ms, and the DFFT averaged 1.8 ms. (using O() notation, the DFFT is 341 times faster).

  -COMPARISON: With 65536 samples and 65536 samples/1 radian, the DFT averaged 23 SECONDS, and the DFFT averaged 0.038 SECONDS. (using O() notation, the DFFT is 4096 times faster)

  -COMPARISON: With 131072 samples and 131072 samples/1 radian, the DFT just got core dumped... The DFFT averaged 86.25ms. (using O() notation, the DFFT is 7710 times faster)

  -COMPARISON: With 262144 samples and 262144 samples/1 radian, the DFT probably got core dumped because the DFFT got core dumped this time.

  -COMPLETED on December 24th 2022 (kind of late because of exams and understanding took a while)

3. Implement an algorithm to perform an inverse DFFT
  -This requires inversing the DFT matrix and multiplying by the recipricol of the number of samples
  -COMPLETED on December 26th, 2022

4. Implement an algorithm to perform a convolution of two arrays of numbers
  -COMPLETED on December 27th, 2022
