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
#include <math.h>
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
    char bkt[16]; //for reading in .csv files for testing
      
    int i=0,j=0,header_lines=0,h=0,y=0;

    double *in;//real input array
    double totalns=0,totals=0,itotalns=0,itotals=0;//timing vars

    std::string input_type = "";
    if(argc > 1){
        input_type = argv[1];
    }
//  std::string filepathin = "radar/radar_data/";
//  std::string filepathout = "radar/radar_data_processed/";
//  std::string filename = "2015-11-18_193605_0001";

    std::string filepathin = "";
    std::string filepathout = "";
    std::string filename = "";
    std::string hann_filename = "Hann506.csv";	
    std::string file_type = ".csv";	
    //int file_num =0;

    ///////////////////////////////////////////////////////////
    // Data Structures
    ///////////////////////////////////////////////////////////
    double hann[506];

    vector<double> valc(2);//vector for real+cplx pairs
    vector<double> mag(128);//vector to hold array megnitudes can be made more efficient by having a 128 vector of 506 vectors but it's left this way for ease of use.

    vector< vector<double> > cols;//vector for columns
    vector< vector<double> >	matrix_2d_A_MAG;//vector to hold columns of magnitudes
    vector< vector<double> > matrix_2d_B_MAG;//vector to hold columns of magnitudes

    vector< vector< vector<double> > > matrix_2d_in_A;//mid calculation 2d array_A
    vector< vector< vector<double> > > matrix_2d_fi_A;//final 2d complex array_A
    vector< vector< vector<double> > > matrix_2d_in_B;//mid calculation 2d array_B
    vector< vector< vector<double> > > matrix_2d_fi_B;//final 2d complex array_B

    ///////////////////////////////////////////////////////////
    // Timing Structures
    ///////////////////////////////////////////////////////////

    struct timespec start, end, tstart, tend,istart, iend, itotal, total;//timing structs

    ///////////////////////////////////////////////////////////
    // Filestreams
    ///////////////////////////////////////////////////////////
    
    fftw_plan plan;//left to right and top to bottom
    fftw_complex *out, *final1, *final2;//FFTW in/out arrays

    ///////////////////////////////////////////////////////////
    // Filestreams
    ///////////////////////////////////////////////////////////

    ifstream Matrix_in;
    ifstream hann_in;
    ofstream mtx_out;
    ofstream mtx_mag_out;

    ///////////////////////////////////////////////////////////
    // Initializing Data Structures
    ///////////////////////////////////////////////////////////

    //final1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (506));
    //final2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (506));
    final1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (512));
    final2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (512));

    for(int i=0;i<512;i++) {
        for(int i=0;i<128;i++) {
            cols.push_back(valc);
        }
        matrix_2d_in_A.push_back(cols);
	matrix_2d_fi_A.push_back(cols);
	matrix_2d_in_B.push_back(cols);
	matrix_2d_fi_B.push_back(cols);
        cols.clear();
    }
    for(int i=0;i<506;i++) {
	matrix_2d_A_MAG.push_back(mag);
	matrix_2d_B_MAG.push_back(mag);
    }
    
    ///////////////////////////////////////////////////////////
    // Configuring Filepaths
    ///////////////////////////////////////////////////////////

	    if(input_type == "PC")
	    {
		    filepathin = "Test_Data/";
		    filename = "fft_result_2015-11-17_103513_716";//PC recording
		    header_lines=2;
		    filepathout = "Test_Data_processed/";
	    }
	    else if(input_type == "Pi2")
	    {		    
		    filepathin = "Test_Data/";
		    filename = "2015-11-17_213141_0004";//Pi2 recording
		    header_lines=4;
		    filepathout = "Test_Data_processed/";
	    }
	    /*
	    else if(input_type == "Cust")
	    {		    	    
		    filepathin = "Test_Data/";
		    filename = argv[2];//Pi2 radar input
		    header_lines=argv[3];
		    filepathout = "Test_Data_processed/";
	    }
	    */
	    else
	    {
		    filepathin = "Test_Data/";
		    filename = "2015-11-17_213141_0004";//Pi2 recording
		    header_lines=4;
		    filepathout = "Test_Data_processed/";
	    }
	    
    ///////////////////////////////////////////////////////////
    // Opening filestreams
    ///////////////////////////////////////////////////////////
    std::cout<<"Input Filepath: "<<filepathin<<filename<<file_type<<"\n";
    Matrix_in.open(filepathin+filename+file_type);
    std::cout<<"Hanning Window Filename: "<<hann_filename<<"\n";
    hann_in.open(hann_filename);

    ///////////////////////////////////////////////////////////
    // FFTW planning
    ///////////////////////////////////////////////////////////

    //plan = fftw_plan_dft_1d(506, final1, final2, FFTW_FORWARD, FFTW_EXHAUSTIVE);
    plan = fftw_plan_dft_1d(512, final1, final2, FFTW_FORWARD, FFTW_EXHAUSTIVE);

    // Check to see if all three files opened correctly
	    
    if(hann_in.is_open())
    {
	for(h=0;h<506;h++)
	{
		assert(hann_in.good());
		hann_in.getline(bkt,16,'\n');
		hann[h]=(double)atof(bkt);
	}
    }
    else
    {
	    std::cout<<"Hanning window file read error! Hanning Window filename: "<<hann_filename<<std::endl;
	    exit(1);	    
    }
	    
    if (Matrix_in.is_open()) {
	std::string line;
	    
	    std::cout<<"Files Opened\n";
        // Read in the 2D matrix file
	for(int i = 0; i < 4; ++i)
	{
	    std::getline(Matrix_in, line);
	}
	//array: [512][128][2]
	//array: [col][row][real+imag]
        for(int col = 0; col < 512; col++) 
	{
	    for(int row = 0; row < 128; row++)
	    {
		//std::cout<<col<<", "<<row<<std::endl;
                std::getline(Matrix_in, line);
		if (!Matrix_in.good())
		{
			break;
		}
                std::stringstream iss(line);
		std::string val;

                std::getline(iss, val, ',');
		if(stoi(val)>=16000)
		{
		matrix_2d_in_A[col][row][0] = 0;
		//matrix_2d_in_A[col][row][1] = 0;
		}
		else
		{
		matrix_2d_in_A[col][row][0] = atof(val.c_str());
		}

                std::getline(iss, val, ',');
		if(stoi(val)>=16000)
		{
		matrix_2d_in_A[col][row][1] = 0;
		//matrix_2d_in_A[col][row][0] = 0;
		}
		else
		{
		matrix_2d_in_A[col][row][1] = atof(val.c_str());
		}
	    }
        }
	
        for(int col = 0; col < 512; col++) 
	{
	    for(int row = 0; row < 128; row++)
	    {
                std::getline(Matrix_in, line);
		if (!Matrix_in.good())
		{
			break;
		}
                std::stringstream iss(line);
		std::string val;
                std::getline(iss, val, ',');
		if(stoi(val)>=16000)
		{
		matrix_2d_in_B[col][row][0] = 0;
		//matrix_2d_in_B[col][row][1] = 0;
		}
		else
		{
		matrix_2d_in_B[col][row][0] = atof(val.c_str());
		}

                std::getline(iss, val, ',');
		if(stoi(val)>=16000)
		{
		matrix_2d_in_B[col][row][1] = 0;
		//matrix_2d_in_B[col][row][0] = 0;
		}
		else
		{
		matrix_2d_in_B[col][row][1] = atof(val.c_str());
		}
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
    
	matrix_2d_in_A.erase(matrix_2d_in_A.begin()+510,matrix_2d_in_A.begin()+512);
	matrix_2d_in_A.erase(matrix_2d_in_A.begin()+0,matrix_2d_in_A.begin()+4);
	matrix_2d_in_B.erase(matrix_2d_in_B.begin()+510,matrix_2d_in_B.begin()+512);
	matrix_2d_in_B.erase(matrix_2d_in_B.begin()+0,matrix_2d_in_B.begin()+4);
    
        
	j=0;
    
	std::cout<<"parsing"<<std::endl;
	
	for(i=0;i<506;i++)
	{
		if((matrix_2d_in_A[i][127] [0]==0) && (matrix_2d_in_A[i][127] [1]==0))
		{
			for(;j<506;j++)
			{
				if((matrix_2d_in_B[j][127][0]==0)&&(matrix_2d_in_B[j][127][1]==0))
				{
					matrix_2d_in_A[i].swap(matrix_2d_in_B[j]);
					j++;
					break;
				}
			}
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

    // Copy data into the in array and perform the fft
	j=0;
	//std::cout<<matrix_2d_in_A.size()<<" "<<matrix_2d_in_A[1].size()<<std::endl;
std::cout<<"Begin FFT matrix A"<<std::endl;
    for(int i=0; i<128; i++){    
        for(int j=0; j<506; j++){

            final1[j][0] = matrix_2d_in_A[j][i][0]*hann[j];
            final1[j][1] = matrix_2d_in_A[j][i][1]*hann[j];
        }
	
        clock_gettime(CLOCK_REALTIME, &tstart);

        fftw_execute(plan);
	
        clock_gettime(CLOCK_REALTIME, &tend);
    	total = diff(tstart, tend);
       	totals+=total.tv_sec;
       	totalns+=total.tv_nsec;


	for(int j=0, k=253; j<253; j++, k++)
	{
	    matrix_2d_fi_A[k][i][0] = final2[j][0];
	    matrix_2d_fi_A[k][i][1] = final2[j][1];
	    matrix_2d_fi_A[j][i][0] = final2[k][0];
	    matrix_2d_fi_A[j][i][1] = final2[k][1];
	}
     }
     
     std::cout<<"Begin FFT matrix B"<<std::endl;
     
     for(int i=0; i<128; i++){
        for(int j=0; j<506; j++){
            final1[j][0] = matrix_2d_in_B[j][i][0]*hann[j];
            final1[j][1] = matrix_2d_in_B[j][i][1]*hann[j];
        }

        clock_gettime(CLOCK_REALTIME, &tstart);

        fftw_execute(plan);
	
        clock_gettime(CLOCK_REALTIME, &tend);
    	total = diff(tstart, tend);
       	totals+=total.tv_sec;
       	totalns+=total.tv_nsec;

	for(int j=0, k=253; j<253; j++, k++)
	{
	    matrix_2d_fi_B[k][i][0] = final2[j][0];
	    matrix_2d_fi_B[k][i][1] = final2[j][1];
	    matrix_2d_fi_B[j][i][0] = final2[k][0];
	    matrix_2d_fi_B[j][i][1] = final2[k][1];
	}
     }
     
#if USE_TIMING_ITER
    std::cout<<"Begin CSV output to file:\n\t"<<filepathout+filename+"(processed"+"_post_"+hann_filename+")"+file_type+"\n\t"+filepathout+filename+"_magnitudes(processed"+"_post_"+hann_filename+")"+file_type<<std::endl;
    mtx_out.open(filepathout+filename+"(processed"+"_post_"+hann_filename+")"+file_type);
    mtx_mag_out.open(filepathout+filename+"_magnitudes(processed"+"_post_"+hann_filename+")"+file_type);

    for(int i=0;i<128;i++){
        for(int j=0; j<506;j++){
            mtx_out<<matrix_2d_fi_A[j][i][0]<<","<<matrix_2d_fi_A[j][i][1]<<",";
	    mtx_mag_out<<sqrt(pow(matrix_2d_fi_A[j][i][0],2.0)+pow(matrix_2d_fi_A[j][i][1],2.0))<<",";
        }
        mtx_out<<"\n";
	mtx_mag_out<<"\n";
    }
    for(int i=0;i<128;i++){
        for(int j=0; j<506;j++){
            mtx_out<<matrix_2d_fi_B[j][i][0]<<","<<matrix_2d_fi_B[j][i][1]<<",";
	    mtx_mag_out<<sqrt(pow(matrix_2d_fi_B[j][i][0],2.0)+(matrix_2d_fi_B[j][i][1],2.0))<<",";
        }
        mtx_out<<"\n";
	mtx_mag_out<<"\n";
    }
    mtx_out.close();
    mtx_mag_out.close();
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

