#ifndef LINE_DETECTOR_H
#define LINE_DETECTOR_H
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <math.h>
#include <sstream>

using namespace cv;

struct Coordinate
{
    double x_;
    double y_;

    Coordinate(){}
    Coordinate(double x, double y){x_=x; y_=y;}
};

struct Line
{
    struct Coordinate start_point_;
    struct Coordinate end_point_;
    double length_;
    double slope_;
    double angle_;
    double intercept_;
    double dist_from_origin_;
    int strength_;

    /* TODO this constructor can prove to be costly.
    Is it better to calculate all this in separate member functions?
    But these are minor calculations anyway, so meh. */
    void evaluate()
    {
        length_ = sqrt(pow(abs(end_point_.x_ - start_point_.x_), 2) + pow(abs(end_point_.y_ - start_point_.x_), 2));
        
        /* TODO a seg fault is bound occur if den = 0  or close to zero. Check what defines close */
        // return(std::min( num/std::max(den, std::numeric_limits<double>::epsilon()) , std::numeric_limits<double>::max() )
        if(!abs(start_point_.x_ - end_point_.x_) < 3) /* these are just pixel values */
        { 
            slope_ = (end_point_.y_ - start_point_.y_) / (end_point_.x_ - start_point_.x_);
            angle_ = std::floor( atan2( (end_point_.x_ - start_point_.x_), (end_point_.y_ - start_point_.y_ ) ) * 180/M_PI);
        }
        else
        { 
            std::cout << "start_point_.x " << start_point_.x_ << "start_point_,y" << start_point_.y_ <<std::endl;
            slope_ = 1000.0; /* not sure about this */
            angle_ = 90;
        }
        intercept_ = (-slope_ * start_point_.x_) + start_point_.y_;    
        dist_from_origin_ = abs(intercept_) / sqrt(1 + pow(slope_,2));
        strength_ = 1; /* one line detected implies strength is 1 */
    }

    Line(){};
    
    Line(Vec4i opencv_line)
    {
        start_point_ = Coordinate(opencv_line[0], opencv_line[1]);
        end_point_ = Coordinate(opencv_line[2], opencv_line[3]);
        evaluate();      
    }

    Line(Coordinate start_point, Coordinate end_point)
    {
        start_point_ = start_point;
        end_point_ = end_point;
        evaluate();      
    }

};

class Line_detector
{
    public:
        Line_detector(std::vector<Vec4i> original_lines); /* ensure a vector is passed if only one line is spitted out by HT */
        cv::Vec4i remove_duplicates(); // returns vector of openCV lines 
        std::vector<cv::Vec4i> from_Lines_struct_to_opencv_lines(std::vector<Line> line_struct_vector); // utility method to convert from variable unique_lines_ to a vector of openCV lines (vector<Vec4i>) 
        cv::Vec4i from_Lines_struct_to_opencv_lines(Line line_struct);

        std::vector<Line> original_lines_;
        std::vector<Line> unique_lines_;

    private:
        double threshold_length_; /* this threshold defines our (fake) "ROI" */  
        // dunno if this is a good idea : Rect ROI_; 
};

#endif /* LINE_DETECTOR_H */