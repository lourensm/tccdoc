
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
  

TODO?


TODO: document the pixelrotationresults, start encoding the rotation, the
box. (Use imagemagick convert).

We need to have blobs available for the rotation detection later.
So the setup of callbacks should be combined.

Detect the pagenumbers and test whether their position is precise.
Can we recognize, distinguish numbers?

Blobs, blobs, what to do?
Is libspatial useful?
Just try it?

Steps:
1)reliable page positioning, in particular the orientation
2)extract candidate page numbers
3)gather them and recognize them
4)revisit the page positioning, taking the page numebers into consideration/

0)define text rectangles,
1)use imagemagick convert to add the boudary to the original tiff. 

It looks relatively ok:
with img016:
y: 138-8034
x:515 - 5684
rotation: 0.1*97.9/100
But: bit unclear: plot coordinates
But: why is the y range so extended?

Ok, got the rotation ok: see Makefile
I can recognize lines in the text.
First line starts at y =~65  linedistance is about 250



Then, decide on where to place the pages, offset.
First choose thefirst example, rotate around its center of predefined 
window.

Then regenerate page, put border around it.
http://www.codeguru.com/cpp/g-m/bitmap/otherformats/article.php/c4933/Working-with-TIFF-Images.htm
https://workingbarely.wordpress.com
https://workingbarely.wordpress.com/2015/02/27/libtiff-the-night-of-the-missing-tag/

Step by step:
-Can we reopen the tiff?
    First lets try to simply read lines again.
- What datastructure do we need to write a rotated/translated version of
the pages?
    First rescan the whole image, possibly only write a box around the
    text border.
- What origins to use?
- How to deal with the pagenumbers?

A bit awkward:
1)Where to put finalization, combining all results?
Looking at setuppage function:
There combine all data values and define a final call?
But having a global value containing the data is awkward.


http://www.eurasip.org/Proceedings/Eusipco/Eusipco2000/SESSIONS/WEDAM/PO4/CR1860.PDF

tesseract..
http://www.imageprocessingplace.com/downloads_V3/root_downloads/tutorials/contour_tracing_Abeer_George_Ghuneim/square.html
http://www.emanueleferonato.com/2013/03/01/using-marching-squares-algorithm-to-trace-the-contour-of-an-image/

INTERMEZZO
Get acquainted with libspatial.

Seems to be ok.

I would like to add a blob shape..?
Implement minimumDistance..


There is a boost spatial indexing library but also boost polygon?
http://www.boost.org/doc/libs/1_57_0_b1/libs/polygon/doc/index.htm
it seems that libspatialindex only defines rectangular space.

boost:
http://stackoverflow.com/questions/22909171/boostgeometry-nearest-neighbors-using-a-circle


I considered using spatial indexing to further analyse pages.
I would then have a more detailed look at libspatialindex and boost libraries.
Both seem to only distinguish rectangles.

But, our analysis provides horizontal structuring by means of lines.

As a first step I could look for pagenumbers relative to the pageboundaries.
See 13tr.md

Ok start with this approach!

Keep the blobwindows as is. If the expected range is outside the pageblobbox,
warn.

Detect blobs in the expected place.


