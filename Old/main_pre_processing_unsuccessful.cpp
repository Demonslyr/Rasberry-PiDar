///////////////////////////////////////////////////////////////////
// Title: Radar Fourier Transform
//
// Copyright (C) 2014  Daniel Murtha
//
// This file is distributed under GPLv2 Licence
///////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////
// Conditional Compilation Options
///////////////////////////////////////////////////////////////////
#define USE_TIMING_ITER 1
#define BUFSIZE  1048576

///////////////////////////////////////////////////////////////////
// Included Files
///////////////////////////////////////////////////////////////////
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;

///////////////////////////////////////////////////////////////////
// Time difference method
///////////////////////////////////////////////////////////////////

timespec diff(timespec start, timespec end)
{
    timespec temp;
    if ((end.tv_nsec-start.tv_nsec)<0) {
        temp.tv_sec = end.tv_sec-start.tv_sec-1;
        temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    } else {
        temp.tv_sec = end.tv_sec-start.tv_sec;
        temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    return temp;
}

///////////////////////////////////////////////////////////////////
// Main method
///////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{

    ///////////////////////////////////////////////////////////
    // Variables
    ///////////////////////////////////////////////////////////
    char bkt[14]; //for reading in .csv files for testing
      
    int i=0,j=0,h=0,k=0,id=0,y=0;

    double *in;//real input array
    double totalns=0,totals=0,itotalns=0,itotals=0;//timing vars
 
    std::string filepathin = "radar/radar_data/";
    std::string filepathout = "radar/radar_data_processed/";
    std::string filename = "2015-11-18_193605_0001";
    std::string file_type = ".csv";
    //int file_num =0;

    vector<double> valc(2);//vector for real+cplx pairs

    vector< vector<double> > cols;//vector for columns

    vector< vector< vector<double> > > matrix_2d_in_A;//mid calculation 2d array_A
    vector< vector< vector<double> > > matrix_2d_fi_A;//final 2d complex array_A
    vector< vector< vector<double> > > matrix_2d_in_B;//mid calculation 2d array_B
    vector< vector< vector<double> > > matrix_2d_fi_B;//final 2d complex array_B

    struct timespec start, end, tstart, tend,istart, iend, itotal, total;//timing structs

    fftw_plan plan;//left to right and top to bottom
    fftw_complex *out, *final1, *final2;//FFTW in/out arrays

    ///////////////////////////////////////////////////////////
    // Filestreams
    ///////////////////////////////////////////////////////////

    ifstream Matrix_in;
    ofstream mtx_out;

    ///////////////////////////////////////////////////////////
    // Initializing Data Structures
    ///////////////////////////////////////////////////////////

    final1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (506));
    final2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (506));

    for(int i=0;i<128;i++) {
        for(int i=0;i<512;i++) {
            cols.push_back(valc);
        }
        matrix_2d_in_A.push_back(cols);
	matrix_2d_fi_A.push_back(cols);
	matrix_2d_in_B.push_back(cols);
	matrix_2d_fi_B.push_back(cols);
        cols.clear();
    }

    ///////////////////////////////////////////////////////////
    // Opening filestreams
    ///////////////////////////////////////////////////////////
    cout<<filepathin+filename+file_type+"\n";
    Matrix_in.open(filepathin+filename+file_type);

    ///////////////////////////////////////////////////////////
    // FFTW planning
    ///////////////////////////////////////////////////////////

    //ltr_trans = fftw_plan_dft_r2c_1d(1024, in, out, FFTW_EXHAUSTIVE);
    plan = fftw_plan_dft_1d(512, final1, final2, FFTW_FORWARD, FFTW_EXHAUSTIVE);
    //plan = fftw_plan_dft_1d(512, final1, final2, FFTW_FORWARD, FFTW_EXHAUSTIVE);

    // Check to see if all three files opened correctly
    if (Matrix_in.is_open()) {
	std::string line;
        // Read in the 2D matrix file
	for(int i = 0; i < 4; ++i)
	{
	    std::getline(Matrix_in, line);
	}
        for(int col = 0; col < 512; ++col) {
	    for(int row = 0; row < 128; ++row)
	    {
                std::getline(Matrix_in, line);
		if (!Matrix_in.good())
                    break;
                std::stringstream iss(line);
		std::string val;

                std::getline(iss, val, ',');
		if(stoi(val)>=16000)
		matrix_2d_in_A[row][col][0] = 0;
		else
		matrix_2d_in_A[row][col][0] = atof(val.c_str());

                std::getline(iss, val, ',');
		if(stoi(val)>=16000)
		matrix_2d_in_A[row][col][1] = 0;
		else
		matrix_2d_in_A[row][col][1] = atof(val.c_str());
	    }
        }
        for(int col = 0; col < 512; ++col) {
	    for(int row = 0; row < 128; ++row)
	    {
                std::getline(Matrix_in, line);
		if (!Matrix_in.good())
                    break;
                std::stringstream iss(line);
		std::string val;

                std::getline(iss, val, ',');
		if(stoi(val)>=16000)
		matrix_2d_in_B[row][col][0] = 0;
		else
		matrix_2d_in_B[row][col][0] = atof(val.c_str());

                std::getline(iss, val, ',');
		if(stoi(val)>=16000)
		matrix_2d_in_B[row][col][1] = 0;
		else
		matrix_2d_in_B[row][col][1] = atof(val.c_str());
	    }
        }
    } else {
        cout<<"Error opening the file.\n";
        exit(1);
    }

        ///////////////////////////////////////////////////////////
        // Closing unnecessary filestreams
        ///////////////////////////////////////////////////////////

    Matrix_in.close();

        ///////////////////////////////////////////////////////////
        // Input Formatting
        ///////////////////////////////////////////////////////////
    
	std::cout<<"Begin parsing"<<std::endl;
    
	//matrix_2d_in_A.erase(matrix_2d_in_A.begin()+510,matrix_2d_in_A.begin()+511);
	//matrix_2d_in_A.erase(matrix_2d_in_A.begin()+0,matrix_2d_in_A.begin()+3);
	//matrix_2d_in_B.erase(matrix_2d_in_B.begin()+510,matrix_2d_in_B.begin()+511);
	//matrix_2d_in_B.erase(matrix_2d_in_B.begin()+0,matrix_2d_in_B.begin()+3);
    
        
	i=0;
	j=0;
    
    
	std::cout<<"parsing"<<std::endl;
	while(i<506)
	{
    		if ((matrix_2d_in_A[127][i][0] == 0) && (matrix_2d_in_A[127][i][1] == 0))
		{
			std::cout<<"Found a zero column at "<<i<<std::endl;
			while(j<506)
			{
				std::cout<<j<<std::endl;
				if((matrix_2d_in_B[127][j][0] == 0) && (matrix_2d_in_B[127][j][1] == 0))
				{
					std::cout<<"matched a zero column at "<<j<<std::endl;
					matrix_2d_in_A[i].swap(matrix_2d_in_B[j]);
					j++;
					break;
				}
				else
				{
					j++;
				}
			}
		i++;
		}
		else
		{
			i++;
		}
	}
    
    std::cout<<"Done parsing"<<std::endl;
    
#if USE_TIMING_ITER
for(int p=0; p < USE_TIMING_ITER; p++) {
    // Start the performance timing here///

#endif // USE_TIMING_ITER



	///////////////////////////////////////////////////////////
	// FFT Process
    	///////////////////////////////////////////////////////////

	//switches between transforming the A buffer or the B buffer
   

    // Copy data into the in array and perform the second fft
	std::cout<<matrix_2d_in_A.size()<<" "<<matrix_2d_in_A[1].size()<<std::endl;
std::cout<<"Begin FFT matrix A"<<std::endl;
    for(int i=0; i<128; i++){
        for(int j=0; j<506; j++){
            final1[j][0] = matrix_2d_in_A[i][j][0];
            final1[j][1] = matrix_2d_in_A[i][j][1];
        }
	
std::cout<<"Start FFT matrix A"<<std::endl;
        clock_gettime(CLOCK_REALTIME, &tstart);

        fftw_execute(plan);
	
        clock_gettime(CLOCK_REALTIME, &tend);
    	total = diff(tstart, tend);
       	totals+=total.tv_sec;
       	totalns+=total.tv_nsec;
std::cout<<"End FFT matrix A "<<i<<std::endl;
	//////////////////////////////
	///
	///The program seg 
	///faults here on i = 111
	///
	//////////////////////////////
	for(int j=0, k=253; j<253; j++, k++)
	{
	    matrix_2d_fi_A[i][k][0] = final2[j][0];
	    matrix_2d_fi_A[i][k][1] = final2[j][1];
	    matrix_2d_fi_A[i][j][0] = final2[k][0];
	    matrix_2d_fi_A[i][j][1] = final2[k][1];
	}
     }
     
     std::cout<<"Begin FFT matrix B"<<std::endl;
     
     for(int i=0; i<128; i++){
        for(int j=0; j<506; j++){
            final1[j][0] = matrix_2d_in_B[i][j][0];
            final1[j][1] = matrix_2d_in_B[i][j][1];
        }	

        clock_gettime(CLOCK_REALTIME, &tstart);

        fftw_execute(plan);
	
        clock_gettime(CLOCK_REALTIME, &tend);
    	total = diff(tstart, tend);
       	totals+=total.tv_sec;
       	totalns+=total.tv_nsec;

	for(int j=0, k=253; j<253; j++, k++)
	{
	    matrix_2d_fi_B[i][k][0] = final2[j][0];
	    matrix_2d_fi_B[i][k][1] = final2[j][1];
	    matrix_2d_fi_B[i][j][0] = final2[k][0];
	    matrix_2d_fi_B[i][j][1] = final2[k][1];
	}
     }
     
#if USE_TIMING_ITER
    std::cout<<"Begin CSV output"<<std::endl;
    mtx_out.open(filepathout+filename+"(processed)"+file_type);

    for(int i=0;i<128;i++){
        for(int j=0; j<506;j++){
            mtx_out<<matrix_2d_fi_A[i][j][0]<<","<<matrix_2d_fi_A[i][j][1]<<",";
        }
        mtx_out<<"\n";
    }
    for(int i=0;i<128;i++){
        for(int j=0; j<506;j++){
            mtx_out<<matrix_2d_fi_B[i][j][0]<<","<<matrix_2d_fi_B[i][j][1]<<",";
        }
        mtx_out<<"\n";
    }
    mtx_out.close();
    }   // end of timing loop
    
    cout<< "Time: "<< (totals)/(USE_TIMING_ITER);
    cout<<"sec " << (totalns/(USE_TIMING_ITER))/1000000 <<"msec\n";

    #endif // USE_TIMING_ITER

    ///////////////////////////////////////////////////////////
    // Cleaning up
    ///////////////////////////////////////////////////////////

    fftw_destroy_plan(plan);
    fftw_free(final1);
    fftw_free(final2);

    return(0);
}

