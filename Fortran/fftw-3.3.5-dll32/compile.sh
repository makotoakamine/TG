gfortran -c -o test.o test.f90
#gfortran -o test_fftw.exe test.o -LC:\fftw-3.3.5-dll32 -lfftw3-3
gfortran -o test_fftw.exe test.o -L/user/local/lib -llibfttw3
./test_fftw

#gfortran -c -o CrFe.o CrFe-spinodal-773K-00.f90
#gfortran -o CrFe.exe CrFe.o -LC:\fftw-3.3.5-dll32 -lfftw3-3
#./CrFe.exe
