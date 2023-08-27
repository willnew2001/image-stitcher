# image-stitcher
Matlab program that stitches two images with shared content together to make one larger image.

HOW TO USE:
Please make sure you have 'image1.png' and 'image2.png' in the same directory as each part
to get an output. These images must be taken side by side and have some amount of overlap.
See the example images in this repo. 'image1' must be the left image and 'image2' must be
the right.

You may also need the Image Processing Toolbox to use the 'imgaussfilt' function. Matlab 
should prompt you to install it when you run the program.

RESOLUTION AND RUNTIME
The program takes longer to run on larger images. I've found a good balance of quality to
speed can be found in images of approx size 504x378. You can size your own images down by using the
'scale_factor' variable on line 10. The example scale factor of 8 scales the source images down
by a factor of 8. Regardless, it will take a few moments to run depending on hardware.

RANDOMNESS INVOLVED
A random number generator is used to select optimal point correspondences using RANSAC. 
The program is currently unseeded and therefore the results of the image stitching are
non-deterministic. To set the seed, go to line 17. To set whether or not you even
want it seeded, go to line 18. I found good results with seed 1234 with the default images.
