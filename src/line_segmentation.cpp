#include "pcnn_simple.h"
#include "network.h"
#include <unordered_set>

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <iostream>
#include <stdio.h>

using namespace cv;
using namespace std; // hashtag bad practice

Mat src, edges;
Mat src_gray;

vector< vector<double> > pcnn_dynamic_get_output(pcnn_dynamic dynamic, int num_osc) 
{
	vector< vector<double> > pcnn_data;
	pcnn_data.resize( dynamic.size() );

	for (int step = 0; step < dynamic.size(); step++)
	{
		pcnn_data[step].resize(num_osc);
		pcnn_network_state & current_state = dynamic[step];

		for (int i = 0; i < num_osc; i++) 
		{
			pcnn_data[step][i] = dynamic.dynamic_oscillator_at(step, i);
			// cout << "step " << step << ", osc_number "<<i <<", output " <<current_state.m_output[i]<<endl;
		}
	}

	return pcnn_data;
}

double safe_division(double num, double den) 
{
	if (abs(den) > std::numeric_limits<double>::epsilon())
		return num/den;
	else
		return num/std::numeric_limits<double>::epsilon();
}

// generates, simulates PCNN. Returns a vector of outputs at a hardcoded timestep 
void pcnn_ensemble_simulate(
		const unsigned int num_osc, 
		const unsigned int steps, 
		const conn_type type_conn, 
		const pcnn_stimulus & stimulus, 
		const size_t height, 
		const size_t width,
		pcnn_dynamic & dynamic,
		const pcnn_parameters * const params = nullptr) 
{
	// use default params if not specified
	pcnn_parameters parameters;
	if (params != nullptr) {
		parameters = *params;
	}

	// class pcnn has a vector of pcnn_oscillators, which are structs which have feeding, linking, output, threshold as members
	// class pcnn also has a member pcnn_params
	pcnn network(num_osc, type_conn, height, width, parameters);

	// pcnn_dynamic is a class that represents the dynamic output of pcnn. 
	// It's inherited from templated dynamic_data and uses "pcnn_network_state" (in place of the template)
	// dynamic_data::output_dynamics member stores a std::vector of pcnn_network_state 
	// pcnn_network_state is a struct with members std::vector<double> m_output and double m_time.
	// pcnn_dynamic dynamic;

	// Performs static simulation of oscillatory network based on Hodgkin-Huxley neuron model.
	network.simulate(steps, stimulus, dynamic);

	// return dynamic;

	/////* Not needed right now ////

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

int animated_display(char* window_name, Mat image, int animation_delay)
{
	imshow(window_name, image);
	int c = waitKey(animation_delay);
	if(c >= 0) {return -1;}
	return 0;
}

int main( int argc, char** argv )
{	
	const char* filename = argv[1];
	VideoCapture input_video(filename);
	if(!input_video.isOpened())
	  throw "error while reading video";

	VideoWriter output_video;

	// arguments are codec, frame per sec, and resolution.
	// Multiply by 2 to accomodate orig and pcnn output
	output_video.open("output_video.avi", CV_FOURCC('X','V','I','D'), 30,cv::Size(2*input_video.get(CV_CAP_PROP_FRAME_WIDTH), input_video.get(CV_CAP_PROP_FRAME_HEIGHT)), true);

    double red, green, blue, hue, saturation, intensity;
    double num, den;
    int intensity_quantized, intensity_quantized_visualize;	
    int rows, cols; // image size, no. of pixels = no of PCNNs
    
    // Containers to define the PCNN ensemble 
    pcnn_stimulus pcnn_stimulus_intensity; // note that pcnn_stimulus is std::vector<double>. (auto-conversion)  
    int pcnn_steps = PCNN_NO_OF_STEPS;     // no of pulses/iterations of PCNN.  


	// for(;;)
	// {
		// Read the current frame from the video
		// input_video >> src;
		src = imread(argv[1], CV_LOAD_IMAGE_COLOR); 
		if(src.empty())
		{ 
			std::cout<<"couldn't read image!"<<std::endl;
		}
		resize(src, src, Size(640, 360), 0, 0, INTER_CUBIC);
		// resize(src, src, Size(640, 360), 0, 0, INTER_CUBIC);
		Mat hsi(src.rows, src.cols, src.type());

		// convert current frame from RGB to HSI. The "quantized" I value will be used as a stimulus to the PCNN filter. 
		std::cout << src.rows << " "<< src.cols<<endl;

		for(int i = 0; i < src.rows; i++)
		{
			for(int j = 0; j < src.cols; j++)
			{
				blue = src.at<Vec3b>(i, j)[0];
				green = src.at<Vec3b>(i, j)[1];
				red = src.at<Vec3b>(i, j)[2];

				intensity = (blue + green + red) / 3; // [0, 255]
				intensity_quantized = floor(intensity/25); // quantized intensity, in 64 levels. [0:1;63]
				intensity_quantized_visualize = 25*intensity_quantized; // throttle up to [0:4:255] to visualize what's going on 
		
				// the intensity_quantized is the input to our PCNN. 
				// std::cout<< i << "   " << j << "    " << intensity<<endl;
				pcnn_stimulus_intensity.push_back(intensity);

				hsi.at<Vec3b>(i, j)[0] = intensity_quantized_visualize; // [0:1:63]
				hsi.at<Vec3b>(i, j)[1] = intensity_quantized_visualize; // [0:1:63]
				hsi.at<Vec3b>(i, j)[2] = intensity_quantized_visualize; // [0:1:63]

			}
		}

		vector< vector<double> > pcnn_result;
		pcnn_dynamic pcnn_dynamic_data_result;
		pcnn_ensemble_simulate(pcnn_stimulus_intensity.size(), pcnn_steps, conn_type::GRID_EIGHT, pcnn_stimulus_intensity, src.rows, src.cols, pcnn_dynamic_data_result);

		pcnn_result = pcnn_dynamic_get_output(pcnn_dynamic_data_result, src.cols*src.rows);

		// populate a vector of fake grayscale openCV image to visualize what's going on.
		// Following block could be optimized. But in the end, we ll need to visualize at a specified time step only, so it doesn't matter
		vector<Mat> pcnn_images;
		vector<double> pcnn_result_current_step;

		cout<<src.rows <<"  " <<src.cols<<endl;
		for(int index = 0; index < pcnn_steps; index++)
		{	
			pcnn_result_current_step = pcnn_result[index];
			Mat pcnn_image_current (src.rows, src.cols, src.type());
			for(int i = 0; i < src.rows; i++)
			{
				for(int j = 0; j < src.cols; j++)
				{
					pcnn_image_current.at<Vec3b>(i,j)[0] = 255 * (pcnn_result_current_step[i*src.cols +j]);
					pcnn_image_current.at<Vec3b>(i,j)[1] = 255 * (pcnn_result_current_step[i*src.cols +j]);
					pcnn_image_current.at<Vec3b>(i,j)[2] = 255 * (pcnn_result_current_step[i*src.cols +j]);
				}
			}

			// medianBlur (pcnn_image_current, pcnn_image_current, 3);
			// cout << index <<endl <<endl << pcnn_image_current<<endl<<endl<<endl<<endl<<endl<<endl<<endl<<endl;
			pcnn_images.push_back(pcnn_image_current);

		    // Find (dark and light) noisy pixels and apply median filter
		}

		namedWindow("RGB image", CV_WINDOW_AUTOSIZE);
		namedWindow("HSI image", CV_WINDOW_AUTOSIZE);

		resize(src, src, Size(640, 360), 0, 0, INTER_CUBIC);
		resize(hsi, hsi, Size(640, 360), 0, 0, INTER_CUBIC);
		
		imshow("RGB image", src);
		imshow("HSI image", hsi);

		namedWindow("pcnn result 1", CV_WINDOW_AUTOSIZE);
		imshow("pcnn result 1", pcnn_images[1]);

		namedWindow("pcnn result 3", CV_WINDOW_AUTOSIZE);
		imshow("pcnn result 3", pcnn_images[3]);

		namedWindow("pcnn result 5", CV_WINDOW_AUTOSIZE);
		imshow("pcnn result 5", pcnn_images[5]);

		namedWindow("pcnn result 7", CV_WINDOW_AUTOSIZE);
		imshow("pcnn result 7", pcnn_images[7]);

		namedWindow("pcnn result 9", CV_WINDOW_AUTOSIZE);
		imshow("pcnn result 9", pcnn_images[9]);

		namedWindow("pcnn result 11", CV_WINDOW_AUTOSIZE);
		imshow("pcnn result 11", pcnn_images[11]);
		
		namedWindow("pcnn result 13", CV_WINDOW_AUTOSIZE);
		imshow("pcnn result 13", pcnn_images[13]);

		namedWindow("pcnn result 15", CV_WINDOW_AUTOSIZE);
		imshow("pcnn result 15", pcnn_images[15]);

		namedWindow("pcnn result 17", CV_WINDOW_AUTOSIZE);
		imshow("pcnn result 17", pcnn_images[17]);

		namedWindow("pcnn result 19", CV_WINDOW_AUTOSIZE);
		imshow("pcnn result 19", pcnn_images[19]);

		// Following for is a basic sketch of what's needed for tiling and scaling the images corresponding
		// to the pulsed output of PCNN in a single window. It's not a priority right now. 

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
