
* Need to define different specifications of actions on scans/tiff files.
  For now we have a test tif file "kin.tif" and we work on analysing
  scanned pages of some book.

* We need to pass on to the program how to analyse the tif file.
  For now we simply define the activities in datastructures in the source.
  But we still need an option for this.

* Analysis and cleanup of the source.
  * Let us do all actions starting from a segment.
  * How to deal with the detected "blobs"?
* Details: How to define "blob' actions?
  

First split off rotation and position detection.

We need to have blobs available for the rotation detection.
So the setup of callbacks should be combined.
In steps:

* Encode the blob detection as an LArea based one.
* Implies additional callbacks at end_of_line and end of image

* Improve the list of LActiveAreas:
allocation should be done by analysescan.c functionality.
There should be a free_data callback called at end and the 
LActiveArea elements should be freed in analysescan.c
* separate out the pixelrotation functionality from activeareas.
* make the activeAreas a local datatype?

*Reason for complexity is that we wish to remove the activeareas cleanly.