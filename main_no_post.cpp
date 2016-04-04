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
    char bkt[16]; //for reading in .csv files for testing
      
    int i=0,j=0,header_lines=0,h=0,y=0;

    double *in;//real input array
    double totalns=0,totals=0,FFTns=0,FFTs=0;//timing vars
 
    //Vertical hanning window
    //double hanning_2[512] = {0,3.7797e-05,0.00015118,0.00034013,0.00060463,0.00094463,0.0013601,0.0018509,0.0024171,0.0030584,0.0037749,0.0045665,0.0054329,0.0063741,0.0073899,0.0084803,0.0096449,0.010884,0.012196,0.013583,0.015043,0.016576,0.018182,0.019862,0.021614,0.023438,0.025334,0.027302,0.029341,0.031452,0.033633,0.035885,0.038207,0.040599,0.043061,0.045591,0.04819,0.050858,0.053593,0.056396,0.059266,0.062203,0.065205,0.068274,0.071408,0.074606,0.077869,0.081196,0.084586,0.088038,0.091554,0.09513,0.098769,0.10247,0.10623,0.11004,0.11392,0.11786,0.12185,0.1259,0.13001,0.13417,0.13839,0.14266,0.14699,0.15137,0.1558,0.16029,0.16483,0.16941,0.17405,0.17874,0.18347,0.18826,0.19309,0.19796,0.20288,0.20785,0.21286,0.21792,0.22301,0.22815,0.23333,0.23855,0.24381,0.24911,0.25445,0.25982,0.26523,0.27068,0.27616,0.28167,0.28722,0.2928,0.29841,0.30405,0.30972,0.31542,0.32115,0.32691,0.33269,0.33849,0.34432,0.35018,0.35605,0.36195,0.36787,0.37381,0.37977,0.38574,0.39174,0.39775,0.40377,0.40981,0.41587,0.42193,0.42801,0.4341,0.4402,0.44631,0.45243,0.45855,0.46468,0.47081,0.47695,0.4831,0.48924,0.49539,0.50154,0.50768,0.51383,0.51998,0.52612,0.53225,0.53839,0.54451,0.55063,0.55675,0.56285,0.56894,0.57503,0.5811,0.58716,0.59321,0.59924,0.60526,0.61126,0.61725,0.62321,0.62916,0.63509,0.641,0.64689,0.65275,0.6586,0.66441,0.67021,0.67598,0.68172,0.68743,0.69312,0.69877,0.7044,0.70999,0.71556,0.72109,0.72658,0.73205,0.73748,0.74287,0.74822,0.75354,0.75882,0.76406,0.76926,0.77442,0.77954,0.78462,0.78965,0.79464,0.79958,0.80448,0.80934,0.81414,0.8189,0.82361,0.82827,0.83289,0.83745,0.84196,0.84642,0.85083,0.85518,0.85948,0.86373,0.86792,0.87205,0.87613,0.88015,0.88412,0.88802,0.89187,0.89566,0.89939,0.90306,0.90667,0.91021,0.9137,0.91712,0.92048,0.92377,0.927,0.93017,0.93327,0.9363,0.93927,0.94218,0.94501,0.94778,0.95048,0.95312,0.95568,0.95818,0.96061,0.96296,0.96525,0.96747,0.96961,0.97169,0.97369,0.97562,0.97748,0.97927,0.98099,0.98263,0.9842,0.9857,0.98712,0.98847,0.98975,0.99095,0.99207,0.99313,0.99411,0.99501,0.99584,0.99659,0.99727,0.99788,0.9984,0.99886,0.99923,0.99954,0.99976,0.99991,0.99999,0.99999,0.99991,0.99976,0.99954,0.99923,0.99886,0.9984,0.99788,0.99727,0.99659,0.99584,0.99501,0.99411,0.99313,0.99207,0.99095,0.98975,0.98847,0.98712,0.9857,0.9842,0.98263,0.98099,0.97927,0.97748,0.97562,0.97369,0.97169,0.96961,0.96747,0.96525,0.96296,0.96061,0.95818,0.95568,0.95312,0.95048,0.94778,0.94501,0.94218,0.93927,0.9363,0.93327,0.93017,0.927,0.92377,0.92048,0.91712,0.9137,0.91021,0.90667,0.90306,0.89939,0.89566,0.89187,0.88802,0.88412,0.88015,0.87613,0.87205,0.86792,0.86373,0.85948,0.85518,0.85083,0.84642,0.84196,0.83745,0.83289,0.82827,0.82361,0.8189,0.81414,0.80934,0.80448,0.79958,0.79464,0.78965,0.78462,0.77954,0.77442,0.76926,0.76406,0.75882,0.75354,0.74822,0.74287,0.73748,0.73205,0.72658,0.72109,0.71556,0.70999,0.7044,0.69877,0.69312,0.68743,0.68172,0.67598,0.67021,0.66441,0.6586,0.65275,0.64689,0.641,0.63509,0.62916,0.62321,0.61725,0.61126,0.60526,0.59924,0.59321,0.58716,0.5811,0.57503,0.56894,0.56285,0.55675,0.55063,0.54451,0.53839,0.53225,0.52612,0.51998,0.51383,0.50768,0.50154,0.49539,0.48924,0.4831,0.47695,0.47081,0.46468,0.45855,0.45243,0.44631,0.4402,0.4341,0.42801,0.42193,0.41587,0.40981,0.40377,0.39775,0.39174,0.38574,0.37977,0.37381,0.36787,0.36195,0.35605,0.35018,0.34432,0.33849,0.33269,0.32691,0.32115,0.31542,0.30972,0.30405,0.29841,0.2928,0.28722,0.28167,0.27616,0.27068,0.26523,0.25982,0.25445,0.24911,0.24381,0.23855,0.23333,0.22815,0.22301,0.21792,0.21286,0.20785,0.20288,0.19796,0.19309,0.18826,0.18347,0.17874,0.17405,0.16941,0.16483,0.16029,0.1558,0.15137,0.14699,0.14266,0.13839,0.13417,0.13001,0.1259,0.12185,0.11786,0.11392,0.11004,0.10623,0.10247,0.098769,0.09513,0.091554,0.088038,0.084586,0.081196,0.077869,0.074606,0.071408,0.068274,0.065205,0.062203,0.059266,0.056396,0.053593,0.050858,0.04819,0.045591,0.043061,0.040599,0.038207,0.035885,0.033633,0.031452,0.029341,0.027302,0.025334,0.023438,0.021614,0.019862,0.018182,0.016576,0.015043,0.013583,0.012196,0.010884,0.0096449,0.0084803,0.0073899,0.0063741,0.0054329,0.0045665,0.0037749,0.0030584,0.0024171,0.0018509,0.0013601,0.00094463,0.00060463,0.00034013,0.00015118,3.7797e-05,0,};
    double hann[512];
    std::string input_type = "";
    if(argc > 1){
        input_type = argv[1];
    }
	    
    std::string filepathin = "";
    std::string filepathout = "";
    std::string filename = "";
    std::string hann_filename = "hann2_512.csv";
	    
    //std::string filename = "2015-11-17_213141_0004";header_lines=4;//Pi2 recording
    //std::string filename = "fft_result_2015-11-17_103513_716";header_lines=2;//PC recording
    std::string file_type = ".csv";   
    //int file_num =0;

    vector<double> valc(2);//vector for real+cplx pairs

    vector< vector<double> > cols;//vector for columns

    vector< vector< vector<double> > > matrix_2d_in;//mid calculation 2d array
    vector< vector< vector<double> > > matrix_2d_fi;//final 2d complex array

    struct timespec FFTstart, FFTend, tstart, tend, FFTtotal, total;//timing structs

    fftw_plan plan;//left to right and top to bottom
    fftw_complex *out, *final1, *final2;//FFTW in/out arrays

    ///////////////////////////////////////////////////////////
    // Filestreams
    ///////////////////////////////////////////////////////////

    ifstream Matrix_in;
    ifstream hann_in;
    ofstream mtx_out;

    ///////////////////////////////////////////////////////////
    // Initializing Data Structures
    ///////////////////////////////////////////////////////////

    final1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (512));
    final2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (512));

    for(int i=0;i<256;i++) {
        for(int i=0;i<512;i++) {
            cols.push_back(valc);
        }
        matrix_2d_in.push_back(cols);
	matrix_2d_fi.push_back(cols);
        cols.clear();
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
	    else if(input_type == "Radar")
	    {		    	    
		    filepathin = "radar/radar_data/";
		    filename = argc[2];//Pi2 radar input
		    header_lines=4;
		    filepathout = "radar/radar_data_processed/";
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

    plan = fftw_plan_dft_1d(512, final1, final2, FFTW_FORWARD, FFTW_EXHAUSTIVE);

    // Check to see if files opened correctly
	    
    if(hann_in.is_open())
    {
	for(h=0;h<512;h++)
	{
		assert(hann_in.good());
		hann_in.getline(bkt,16,'\n');
		hann[h]=(double)atof(bkt);
		//std::cout<<hann[h]<<std::endl;
	}
    }
    else
    {
	    std::cout<<"Hanning window file read error! Hanning WIndow filename: "<<hann_filename<<std::endl;
	    exit(1);	    
    }
	    
    if (Matrix_in.is_open()) {
	std::string line;
        // Read in the 2D matrix file
	for(int i = 0; i < header_lines; i++)
	{
	    std::getline(Matrix_in, line);
	}
		//array[256][512][2]
		//array[row][col][real+cplx]
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
				if(stoi(val)>=16000){
					//std::cout<<stoi(val)<<std::endl;
					matrix_2d_in[row][col][0] = 0;
					//matrix_2d_in[row][col][1] = 0;
					break;
				}
				else
				{
					//matrix_2d_in[row][col][0] = atof(val.c_str());
					matrix_2d_in[row][col][0] = hann[col]*(double)atof(val.c_str());
				}

		        std::getline(iss, val, ',');
				if(stoi(val)>=16000)
				{
					//std::cout<<stoi(val)<<std::endl;
					//matrix_2d_in[row][col][0] = 0;
					matrix_2d_in[row][col][1] = 0;
				}
				else
				{
					//matrix_2d_in[row][col][1] = atof(val.c_str());
					matrix_2d_in[row][col][1] = hann[col]*(double)atof(val.c_str());
				}
			
				//std::cout<<matrix_2d_in[row][col][0]<<"+i"<<matrix_2d_in[row][col][1]<<std::endl;
		    }
        }
        for(int col = 0; col < 512; col++) 
        {
		    for(int row = 128; row < 256; row++)
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
					//std::cout<<stoi(val)<<std::endl;
					matrix_2d_in[row][col][0] = 0;
					//matrix_2d_in[row][col][1] = 0;
					break;
				}
				else
				{
					//matrix_2d_in[row][col][0] = atof(val.c_str());
					matrix_2d_in[row][col][0] = hann[col]*(double)atof(val.c_str());
				}

		        std::getline(iss, val, ',');
				if(stoi(val)>=16000)
				{
					//std::cout<<stoi(val)<<std::endl;
					//matrix_2d_in[row][col][0] = 0;
					matrix_2d_in[row][col][1] = 0;
				}
				else
				{
					//matrix_2d_in[row][col][1] = atof(val.c_str());
					matrix_2d_in[row][col][1] = hann[col]*(double)atof(val.c_str());
				}
					//std::cout<<matrix_2d_in[row][col][0]<<"+i"<<matrix_2d_in[row][col][1]<<std::endl;
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

#if USE_TIMING_ITER
for(int p=0; p < USE_TIMING_ITER; p++) {
    // Start the performance timing here///

#endif // USE_TIMING_ITER



	///////////////////////////////////////////////////////////
	// FFT Process
    	///////////////////////////////////////////////////////////
   

    // Copy data into the in array and perform the second fft
    for(int i=0; i<256; i++){
        for(int j=0; j<512; j++){
        	//array[256][512][2]
			//array[row][col][real+cplx]
            final1[j][0] = matrix_2d_in[i][j][0];
            final1[j][1] = matrix_2d_in[i][j][1];
	    //final1[i][0] = 5;
	    //final1[i][1] = 5.5;
		
	    //std::cout<<"Pre FFT: "<< matrix_2d_in[i][j][0]<<"+i"<< matrix_2d_in[i][j][1]<<" recorded as "<<final1[i][0]<<"+i"<<final1[i][1]<<std::endl;
        }	
	
        clock_gettime(CLOCK_REALTIME, &FFTstart);

        fftw_execute(plan);
	
        clock_gettime(CLOCK_REALTIME, &FFTend);
    	FFTtotal = diff(FFTstart, FFTend);
       	FFTs+=FFTtotal.tv_sec;
       	FFTns+=FFTtotal.tv_nsec;
	
	//std::cout<<final2[5][0]<<"+i"<<final2[5][1]<<std::endl;
	//std::cout<<final2[85][0]<<"+i"<<final2[85][1]<<std::endl;
	//std::cout<<final2[175][0]<<"+i"<<final2[175][1]<<std::endl;

	for(int j=0, k=256; j<256; j++, k++)
	{
//	    matrix_2d_fi[i][k][0] = hann[j]*final2[j][0];
//	    matrix_2d_fi[i][k][1] = hann[j]*final2[j][1];
//	    matrix_2d_fi[i][j][0] = hann[k]*final2[k][0];
//	    matrix_2d_fi[i][j][1] = hann[k]*final2[k][1];

	    matrix_2d_fi[i][k][0] = final2[j][0];
	    matrix_2d_fi[i][k][1] = final2[j][1];
	    matrix_2d_fi[i][j][0] = final2[k][0];
	    matrix_2d_fi[i][j][1] = final2[k][1];
	    
	    //std::cout<<"Post FFT: "<<final2[j][0]<<"+i"<<final2[j][1]<<" and "<<final2[k][0]<<"+i"<<final2[k][1]<<std::endl;
	    //std::cout<<final2[j][0]<<"*"<<hann[j]<<"="<<(hann[j]*final2[j][0])<<", This should be the same as "<<matrix_2d_fi[i][k][0]<<std::endl;
	    //std::cout<<"Post FFT: "<<matrix_2d_fi[i][j][0]<<"+i"<<matrix_2d_fi[i][j][1]<<" and "<<matrix_2d_fi[i][k][0]<<"+i"<<matrix_2d_fi[i][k][1]<<std::endl;
	}
     }

#if USE_TIMING_ITER
    std::cout<<"Begin CSV output to file:\n\t"<<filepathout+filename+"(processed"+"_no_post_"+hann_filename+")"+file_type<<std::endl;
    mtx_out.open(filepathout+filename+"(processed"+"_no_post_"+hann_filename+")"+file_type);
    //mtx_out.open(filepathout+filename+"(processed)"+file_type);

    for(int i=0;i<256;i++){
        for(int j=0; j<512;j++){
            mtx_out<<matrix_2d_fi[i][j][0]<<","<<matrix_2d_fi[i][j][1]<<",";
        }
        mtx_out<<"\n";
    }
    mtx_out.close();
    }   // end of timing loop
    
    cout<< "Time: "<< (FFTs)/(USE_TIMING_ITER);
    cout<<"sec " << (FFTns/(USE_TIMING_ITER))/1000000 <<"msec\n";

    #endif // USE_TIMING_ITER

    ///////////////////////////////////////////////////////////
    // Cleaning up
    ///////////////////////////////////////////////////////////

    fftw_destroy_plan(plan);
    fftw_free(final1);
    fftw_free(final2);

    return(0);
}

