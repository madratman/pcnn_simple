- This repo attempts to implement [Towards automatic power line detection for a UAV surveillance system using pulse coupled neural filter and an improved Hough transforma](http://link.springer.com/article/10.1007%2Fs00138-009-0206-y)

- Here's a summary of the paper I wrote a few weeks back. [Summary](https://docs.google.com/document/d/17NT0FqzACVTZ6PO-D2FGP4kk5K-5-ZAiWekqSQnQ1fk/edit?usp=sharing)

- The first thing to do is to implement a filter based on (a simplified version of) Pulse Coupled Neural Network.
For the same, I am building upon the C++ core of https://github.com/annoviko/pyclustering.(which btw seems like an awesome new library)

- To be precise, the first commit is a direct copy of `pcnn.h, pcnn.cpp, network.h, network.cpp` from https://github.com/annoviko/pyclustering/tree/master/ccore/ccore 

- To compile, do 
```
mkdir build
cd build
cmake ..
make
```

- Usage  
`./line_segmentation <path_to_your_image>`

 - Parameters:
You can adjust the following parameters in pcnn_simple.h
```
  - No of time steps = `PCNN_NO_OF_STEPS`   
   
    // Multiplier for the threshold at the current step.    
    double VT;
    
    // Synaptic weight - neighbours influence on linking compartment
    std::vector<double> W = {sqrt(1/2), 1.0, sqrt(1/2),
							 1.0, 1.0, 1.0, 
						     sqrt(1/2), 1.0, sqrt(1/2)};

    // Linking strength in the network.
    double B; (0.2 from paper)

    // step_value by which each oscillator's threshold is decreased if it 
    double step_value; 
```    
For videos, you can uncomment the for loop in `line_segmentation.cpp`, but right now it is pretty slow.  

*(Adding the following stuff for easier collaboration and tracking)*
### Higher level Roadmap of things to be done (as of 25/10/15)
  * Tune pcnn. 
  * Choose at which time step to perform post processing
  * Locating noisy pixels and applying a median *only* to the noisy pixels, as explained in this reference in the above paper [An Adaptive Method for Image Filtering with
Pulse-coupled Neural Networks ](http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1530009&tag=1)
  * Edge detection (see algorithm 2 in paper)
  * K-means clustering followed by the "smart" Hough transform (the voting procedure) to detect power lines. Read paper for details.
    The voting procedure is better in one of the references : [Real-time line detection through an improved Hough transform voting scheme](http://www.sciencedirect.com/science/article/pii/S0031320307001823)
    For kmeans, I'll again modify the source for pyclustering's ccore.

### Development Log 
  - Week 1
    - Late night 19/10/15 : Found pyclustering. Took out required files ccore and building on top of it.
    - Throughout the week : Getting familiar, adding support 2D weight matrices and rectangular images, methods to read output, executable for image/video line detection, [minor contribution](https://github.com/annoviko/pyclustering/pull/260) to pyclustering, mailing the author of pyclustering - who turned out to be pretty helpful and is already addressing [issue #259](https://github.com/annoviko/pyclustering/issues/259), and of course adding and later removing my own bugs and getting confused here and there. 
    - 25/10 : First part of paper is mostly done. PCNN needs tuning. 
    - *Results* can be seen in `results/preliminary_Oct_25' directory. In the screenshots, "pcnn_result_n.png" is the image at time step n of the network. The "HSI image.png" has a misleading name. It's just a monochrome visualization of quantized intensities(in 64 lines, as written in paper) of the original image. The parameters can be seen from the commit from that day. 

### More ideas for power line detection / segmentation. 
  - *Stacked PCNN, in a ConvNet style architecture*. This is a general idea for segmentation. Dunno about power lines.   
 Here the convolution part is replaced by the PCNN equations. The output of the pixel is fed on to the next layer and so on.    
 You could, in theory, tweak which temporal (or "serial") pulse should be passed on to the next network. Or you could just pass on everything as parallel channels.     
But then, you could also cherry pick wich time step for each layer should be chosen and passed on to the next.   
 This whole thing is similar to (CRFs as RNNs) [http://www.robots.ox.ac.uk/~szheng/CRFasRNN.html]

Daniel's ideas(quoting):
  - Learning-based general purpose edge detection. A recent example based on networks similar to the ones AirLab uses is
[Contour Detection Using Cost-Sensitive Convolutional Neural Networks](http://arxiv.org/pdf/1412.6857v5.pdf)  
  However, if we're assuming wires are lines (which is reasonable most of the time) I think we can do better. In essence, the idea is to do sort of an end-to-end, learning-based Hough detector. Instead of mapping from an input image to some sort of binary label image, we would learn a mapping from an input image to a "hough vote map", i.e. an image where each "pixel" represents belief in a line with those theta, rho parameters.  
Use [Sloth](https://github.com/cvhciKIT/sloth) for labeling  
Also similar to [Real-Time Grasp Detection Using Convolutional Neural Networks](http://pjreddie.com/media/files/papers/grasp_detection_1.pdf)

 - Stacking and autocontext
  http://pages.ucsd.edu/~ztu/publication/cvpr08_autocontext.pdf
  https://www.ri.cmu.edu/pub_files/2010/9/munoz_eccv_10.pdf

#### Resources for CNN
  - [ CS231n: Convolutional Neural Networks for Visual Recognition. ](http://cs231n.github.io/)
  - [Deep Learning for Perception](https://computing.ece.vt.edu/~f15ece6504/).   
  [Homework repos](https://github.com/batra-mlp-lab/)