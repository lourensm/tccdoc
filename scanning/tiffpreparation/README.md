# tccdoc/scanning/tiffpreparation
Analysing and changing scanned tiff images.

# Issues
* I had problems with installing and accessing the libtiff library. I ended up with copying some .h files and using functionality (I do not understand) from libtiff to compile my own code.
* I had some problems trying to use libtiff from perl and python. Made me decide to use C. I might switch to C++ later.



# Problem description

In the first stage, I wish to detect the placing of the text in a scanned
book page. The source is a tai chi book.
The scan is a two page scan. The orientation and offset will vary across
scans.

In a first step I try to determine the rotation and offset:
Detect text lines and their orientation.
Detect vertical and horizontal limits of the text.

The next step is detecting the pagenumber (at the bottom). And
adjust the detected textboudary, assuming that the pagenumbers are placed
consistently.

I am now exploring the possibility of detecting connected blobs in the
image. 

 
# Status Mar 21 2016
e
Last work done Jul 21 2015. I remember that I was heavily disappointed because of some bug
in the libtiff library. Drawing rectangles seemetd to go wrong. I seem to have been able to
calculate a bounding box around text. The next step would have been to rotate pages and do some 
rearranging so that all pages of the "13 treatises" book would align. And then perform OCR.
