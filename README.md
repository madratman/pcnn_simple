- This repo attempts to implement [Towards automatic power line detection for a UAV surveillance system using pulse coupled neural filter and an improved Hough transforma](http://link.springer.com/article/10.1007%2Fs00138-009-0206-y)

- Here's a summary of the paper I wrote a few weeks back. [Summary](https://docs.google.com/document/d/17NT0FqzACVTZ6PO-D2FGP4kk5K-5-ZAiWekqSQnQ1fk/edit?usp=sharing)

- The first thing to do is to implement a filter based on (a simplified version of) Pulse Coupled Neural Network.
For the same, I am building upon the C++ core of https://github.com/annoviko/pyclustering.(which btw seems like an awesome new library)

- To be precise, the first commit is a direct copy of `pcnn.h, pcnn.cpp, network.h, network.cpp` from https://github.com/annoviko/pyclustering/tree/master/ccore/ccore 

