#include "pcnn.h"
#include "network.h"
#include <unordered_set>

// #include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <iostream>
#include <stdio.h>

using namespace cv;
using namespace std;

Mat src, edges;
Mat src_gray;

// General TODOs
// - Check if the initial output pulse Y, is set to a zero value matrix in pcnn_simple.cpp
// - calculate the output the each neuron using equation (9), with a threshold larger than 
// the minimum value of U : \Theta = min(U) + 0.01; Y = step(U − \Theta)
// - the edge of image Bin can be obtained by logical operation exclusive disjunction (XOR) 
// on Y0 and Y : Edge = Y0 ⊕ Y

// how to check for inf (std::numeric_limits<double>::max()) ? It will anyway return an error (hardware segfault).

std::vector<double> pcnn_dynamic_get_output(pcnn_dynamic dynamic, int num_osc) {
	// pcnn_dynamic & dynamic = *((pcnn_dynamic *) pointer);

	// pyclustering_package * package = new pyclustering_package((unsigned int) pyclustering_type_data::PYCLUSTERING_TYPE_LIST);
	// package->size = dynamic.size();
	// package->data = new pyclustering_package * [package->size];
	int step = 1;
	pcnn_network_state & current_state = dynamic[step];
	std::vector<double> current_state_result;
	std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"<<std::endl;
	for (int i = 0; i < num_osc; i++) {
		std::cout <<  current_state.m_output[i] << std::endl;
		current_state_result.push_back(current_state.m_output[i]);
		// std::cout <<  current_state.m_output[i] << std::endl;
		// ((pyclustering_package **) package->data)[i] = create_package(&dynamic[i].m_output);
	}

	return current_state_result;
}

double safe_division(double num, double den) 
{
	if (abs(den) > std::numeric_limits<double>::epsilon())
		return num/den;
	else
		return num/std::numeric_limits<double>::epsilon();
}

static void template_dynamic_generation(
		const unsigned int num_osc, 
		const unsigned int steps, 
		const conn_type type_conn, 
		const pcnn_stimulus & stimulus) {

	pcnn_parameters parameters;
	pcnn network(num_osc, type_conn, parameters);

	pcnn_dynamic dynamic;
	network.simulate(steps, stimulus, dynamic);
	
	pcnn_time_signal time_signal;
	dynamic.allocate_time_signal(time_signal);
}

// generates, simulates PCNN. Returns a vector of outputs at a hardcoded timestep 
std::vector<double> pcnn_ensemble_simulate(
		const unsigned int num_osc, 
		const unsigned int steps, 
		const conn_type type_conn, 
		const pcnn_stimulus & stimulus, 
		const pcnn_parameters * const params = nullptr) {

	// use default params if not specified
	pcnn_parameters parameters;
	if (params != nullptr) {
		parameters = *params;
	}

	// class pcnn has a vector of pcnn_oscillators, which are structs which have feeding, linking, output, threshold as members
	// class pcnn also has a member pcnn_params

	pcnn network(num_osc, type_conn, parameters);

	// pcnn_dynamic is a class that represents the dynamic output of pcnn. 
	// It's inherited from templated dynamic_data and uses "pcnn_network_state" (in place of the template)
	// dynamic_data::output_dynamics member stores a std::vector of pcnn_network_state 
	// pcnn_network_state is a struct with members std::vector<double> m_output and double m_time.

	pcnn_dynamic dynamic;

	// Performs static simulation of oscillatory network based on Hodgkin-Huxley neuron model.
	network.simulate(steps, stimulus, dynamic);

	// TODO : return a vector of vectors for all time steps
	std::vector<double> pcnn_state; 
	// pcnn_state.resize(num_osc);
	// dynamic.return_dynamic(3); // time step 3 
	// pcnn_state = dynamic.return_dynamic(3); // time step 3 
	// dynamic.dummy_method_vector();
	dynamic.dummy_method();
	std::cout<<"before call to dynamic method" <<std::endl;
	// dynamic.return_dynamic_test(3); // time step 3
	pcnn_state = pcnn_dynamic_get_output(dynamic, num_osc);
 	std::cout<<pcnn_state.size()<<std::endl;
	std::cout<<"after call to dynamic method"<<std::endl;
	// std::vector<double> temp = {0};
	// return temp;

	return pcnn_state;
	/////* Not needed right now */////

	/* minidoc so that my head doesn't spin after a week reading this code */	

	// pcnn_ensemble data type is basically a typedef of std::vector<unsigned int>;
	// ensemble_data is a class defined in network.h which has 2 member varibles:
	// 1. "value_type" (which is a typedef of a template<typename sync_ensemble_type>)
	// 2. "sync_ensemble" (which is a std::vector<value_type>) 
	// and 3 iterators, ands some functions to facilitate tranversing, pushin, popping, reserving, etc

	// So in the following lines, we are instantiating an ensemble_data object whose "value_type" is this pcnn_ensemble, 
	// whihc btw is basically a std::vector<unsigned int>, dear future self. 
	// Finally pcnn_time_signal is a typedef of typedef std::vector<unsigned int>	

	/* minidoc ends */

	// ensemble_data<pcnn_ensemble> sync_ensembles; /* holder for synchronous oscillators in the current time step */
	// ensemble_data<pcnn_ensemble> spike_ensembles;
	// pcnn_time_signal time_signal;
	
	// // Allocate clusters in line with ensembles of synchronous oscillators where each synchronous ensemble corresponds to only one cluster.
	// dynamic.allocate_sync_ensembles(sync_ensembles);
	
	// // Analyses output dynamic of network and allocates spikes on each iteration as a list of indexes of oscillators.
	// dynamic.allocate_spike_ensembles(spike_ensembles);
	
	// // Analyses output dynamic and calculates time signal (signal vector information) of network output.
	// dynamic.allocate_time_signal(time_signal);

	/////* Not needed right now ends */////
}

int main( int argc, char** argv )
{	
 	//Analyses output dynamic of network and allocates spikes on each iteration as a list of indexes of oscillators.

	const char* filename = argv[1];
	VideoCapture input_video(filename);
	if(!input_video.isOpened())
	  throw "error while reading video";

	VideoWriter output_video;

	output_video.open("output_video.avi", CV_FOURCC('X','V','I','D'), 30,cv::Size(2*input_video.get(CV_CAP_PROP_FRAME_WIDTH), input_video.get(CV_CAP_PROP_FRAME_HEIGHT)), true);
	// arguments are codec, frame per sec, and resolution.
	// Multiply by 2 to accomodate orig and pcnn output

    double red, green, blue, hue, saturation, intensity;
    double num, den;
    int intensity_quantized, intensity_quantized_visualize;	
    int rows, cols; // image size, no. of pixels = no of PCNNs
    
    /* Containers to define the PCNN ensemble */
    pcnn_stimulus pcnn_stimulus_intensity; // note that pcnn_stimulus is std::vector<double>. (auto-conversion)  
    int pcnn_steps = 20; // no of pulses/iterations of PCNN.  

	// for(;;)
	// {
		// Read the current frame from the video
		// input_video >> src;
		src = imread(argv[1], CV_LOAD_IMAGE_COLOR); 
		if(src.empty())
		{ 
			std::cout<<"wtf!";
		}

		// convert current frame from RGB to HSI. The "quantized" I value will be used as a stimulus to the PCNN filter. 
		Mat hsi(src.rows, src.cols, src.type());

		std::cout << src.rows << " "<< src.cols<<endl;
		int test;
		for(int i = 0; i < src.rows; i++)
		{
			for(int j = 0; j < src.cols; j++)
			{
				blue = src.at<Vec3b>(i, j)[0];
				green = src.at<Vec3b>(i, j)[1];
				red = src.at<Vec3b>(i, j)[2];

				intensity = (blue + green + red) / 3; // [0, 255]
				intensity_quantized = floor(intensity/4); // quantized intensity, in 64 levels. [0:1;63]
				intensity_quantized_visualize = 4*intensity_quantized; // throttle up to [0:4:255] to visualize what's going on 
				
				// the intensity_quantized is the input to our PCNN. 
				// std::cout<< i << "   " << j << "    " << intensity_quantized<<endl;
				pcnn_stimulus_intensity.push_back(intensity_quantized);

				int minimum_rgb = 0;
				minimum_rgb = std::min(red, std::min(green, blue));

				saturation = 1 - (minimum_rgb/intensity);

				num = (red - green) + (red - blue);
				den = sqrt( (red - green)*(red - green) + (red - blue)*(green - blue) );
				hue = safe_division(num, den);
				hue = acos(0.5 * hue) * 180/M_PI;

				if(blue > green)
				{
					hue	 = 360 - hue;
				}

				hsi.at<Vec3b>(i, j)[0] = hue; 			 // [0,360]  
				hsi.at<Vec3b>(i, j)[1] = saturation*100; // [0,100]
				hsi.at<Vec3b>(i, j)[2] = intensity_quantized_visualize; // [0:1:63]
			}
		}
		std::cout << pcnn_stimulus_intensity.size()<<endl;

		// for(auto i = pcnn_stimulus_intensity.begin(); i!= pcnn_stimulus_intensity.end(); i++)
		// 	std::cout << * i << std::endl;

		std::vector<double> pcnn_result;
		// pcnn_result.resize(pcnn_stimulus_intensity.size());
		int size_of_result = 153;

		// crop image to square so as to comply with pcnn connection type GRID_EIGHT
		// TODO try to remove this constraint of square images
		pcnn_stimulus_intensity.resize(size_of_result*size_of_result); 


		std::cout<<"main before call to pcnn_ensemble_simulate" <<std::endl;
		pcnn_result = pcnn_ensemble_simulate(pcnn_stimulus_intensity.size(), pcnn_steps, conn_type::GRID_EIGHT, pcnn_stimulus_intensity);
	 	std::cout << pcnn_result.size()<< std::endl;
		std::cout<<"main after call to pcnn_ensemble_simulate" <<std::endl;

		// populate an openCV image to visualize what's going on at hardcoded time step
		Mat pcnn_at_step(size_of_result, size_of_result, CV_64FC1);

		for(int i=0; i < pcnn_at_step.rows; ++i)
			for(int j=0; j < pcnn_at_step.cols; ++j)
			{
				// pcnn_at_step.at<double>(i, j) = pcnn_result.at(i + j);	
				pcnn_at_step.at<Vec3b>(i, j)[0] = pcnn_result.at(i + j);  
				pcnn_at_step.at<Vec3b>(i, j)[1] = pcnn_result.at(i + j);
				pcnn_at_step.at<Vec3b>(i, j)[2] = pcnn_result.at(i + j);
			}

		// pcnn_stimulus stimulus_test { 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
		// template_dynamic_generation(stimulus_test.size(), 20, conn_type::GRID_EIGHT, stimulus_test);
		// pcnn_ensemble_simulate(stimulus_test.size(), 20, conn_type::GRID_EIGHT, stimulus_test);

		// for(auto i = pcnn_result.begin(); i!= pcnn_result.end(); i++)
		// 	std::cout << * i << std::endl;


		// resize(src, src, Size(640, 360), 0, 0, INTER_CUBIC);
		// resize(hsi, hsi, Size(640, 360), 0, 0, INTER_CUBIC);

		namedWindow("RGB image", CV_WINDOW_AUTOSIZE);
		namedWindow("HSI image", CV_WINDOW_AUTOSIZE);
		namedWindow("pcnn result", CV_WINDOW_AUTOSIZE);

		imshow("RGB image", src);
		imshow("HSI image", hsi);
		imshow("pcnn result", pcnn_at_step);

	/*	Mat display = Mat::zeros( src.rows, (2*src.cols) + (20), src.type() );

		//club orginal and pcnn output together in a single window
		src.copyTo(at(display, Rect(0, 0, src.cols, src.rows)));
		pcnn_output.copyTo(Mat(display, Rect(src.cols + 20, 0, src.cols, src.rows)));

		Mat display_resized;
		resize(display, display_resized, Size(), 0.5, 0.5, INTER_CUBIC); // upscale 2x

		imshow(original, display_resized);
		
		// write to output_video
		output_video << display; */

		// waitKey(1);
	// }

	waitKey(0);
	return 0;
}
