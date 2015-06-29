
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
  

TODO FIRST
debug: img012_r.tif
The rectangles do not fit correctly, see o1.tif
Actions:
1)measure xlo,xhi, ylo, yhi
2)measure xtl,ytl  etc

left rectangle?
x1:973, y1:843, x2:6174, y2:8518   14039x9909  334mm x 236
where is middle x: (6174 + 973)/2  x 334/14039 = 85.01 mm
where is middle y: (8518 + 843)/2 x  236/9909 = 111.5 mm
ylow = 175
better: adjust polygon, draw x/2 y/2 and intersection

try:
-draw "polygon 1778.694932,704.989349  6452.062408,1064.345638  5872.322332,8603.766810 1198.954857,8244.410521" -draw "polygon 7730.105749,1001.598210  12926.081824,990.448595  12942.222587,8512.413960 7746.246512,8523.563575"

first draw the intersections with the angles.



Improve on surrounding box, by looking at details of yvalues,xvalues


TODO: document the pixelrotationresults, start encoding the rotation, the
box. (Use imagemagick convert).

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


Drawing line
(x1,y1) - (x2,y2)
width w:

so, in direction 1,2 add w/2
perpendicular, add w/2
define x1,y1 corners:

Draw lines:
a bit more exploring