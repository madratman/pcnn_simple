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

// how to check for inf (std::numeric_limits<double>::max()) ? It will anyway return an error (hardware segfault).
double safe_division(double num, double den) 
{
	if (abs(den) > std::numeric_limits<double>::epsilon())
		return num/den;
	else
		return num/std::numeric_limits<double>::epsilon();
}


int main( int argc, char** argv )
{	
	const char* filename = argv[1];
	VideoCapture input_video(filename);
	if(!input_video.isOpened())
	  throw "error while reading video";

	VideoWriter output_video;

	output_video.open("output_video.avi", CV_FOURCC('X','V','I','D'), 30,cv::Size(2*input_video.get(CV_CAP_PROP_FRAME_WIDTH), input_video.get(CV_CAP_PROP_FRAME_HEIGHT)), true);
	// arguments are codec, frame per sec, and resolution.
	// Multiply by 2 to accomodate orig and pcnn output

	for(;;)
	{
		// Read the current frame from the video
		input_video >> src;
		if(src.empty())
		{ 
			std::cout<<"wtf!";
		}

		// convert current frame from RGB to HSI. The "quantized" I value will be used as a stimulus to the PCNN filter. 
		Mat hsi(src.rows, src.cols, src.type());
	    double red, green, blue, hue, saturation, intensity;
	    double num, den;

		for(int i = 0; i < src.rows; i++)
		{
			for(int j = 0; j < src.cols; j++)
			{
				blue = src.at<Vec3b>(i, j)[0];
				green = src.at<Vec3b>(i, j)[1];
				red = src.at<Vec3b>(i, j)[2];

				intensity = (blue + green + red) / 3;

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
				hsi.at<Vec3b>(i, j)[2] = intensity; 	 // [0, 255]
			}
		}

		namedWindow("RGB image", CV_WINDOW_AUTOSIZE);
		namedWindow("HSI image", CV_WINDOW_AUTOSIZE);

		imshow("RGB image", src);
		imshow("HSI image", hsi);

	/*	Mat display = Mat::zeros( src.rows, (2*src.cols) + (20), src.type() );

		//club orginal and pcnn output together in a single window
		src.copyTo(at(display, Rect(0, 0, src.cols, src.rows)));
		pcnn_output.copyTo(Mat(display, Rect(src.cols + 20, 0, src.cols, src.rows)));

		Mat display_resized;
		resize(display, display_resized, Size(), 0.5, 0.5, INTER_CUBIC); // upscale 2x

		imshow(original, display_resized);
		
		// write to output_video
		output_video << display;*/
		waitKey(1);
	}

	waitKey(0);
	return 0;
}
