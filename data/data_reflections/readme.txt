These files contain the data that was used in the experiments of paper [1]: unary potentials, edge weights (pairwise potentials), and superpixels.
The images are '001_oct' and '020_smallplant' from  

http://melodi.ee.washington.edu/~jegelka/cc/


For each image, there are the following files, where XX = image name:
XX_wts.mat : contains the unary potentials (linear function, or edges to terminal nodes s, t -- this is still in form of a vector and needs to be reshaped into a matrix if you want to view it) and two weight matrices. One is for the horizontal edges, and one for the vertical edges. Using Matlab's imagesc function will help you understand it better.

XX_pYY.mat, where YY is a number. These are superpixel partitions (this is mostly a benchmark for testing algorithms here). YY is the number of superpixels. The file contains a matrix that has the superpixel number for each node / pixel in the image.

The functions that were actually used where of the form
unary(x) + lambda * weights(x) + lambda2 * suppix(x)

The values of lambda and lambda2 we used were:
001_oct  lambda=0,1  lambda2 = 0.002
001_oct  lambda = 0.01  lambda2 = 0.008
020_smallplant  lambda = 0.3  lambda2 = 0.002
020_smallplant  lambda = 0.08  lambda2 = 0.008


If you use this data, cite the following papers:

[1] S. Jegelka, F. Bach and S. Sra. Reflection methods for user-friendly submodular optimization. Neural Information Processing Systems (NIPS), 2013.
[2] S. Jegelka and J. Bilmes. Submodularity beyond Submodular Energies: Coupling Edges in Graph Cuts. IEEE Conference of Computer Vision and Pattern Recognition (CVPR), 2011.

