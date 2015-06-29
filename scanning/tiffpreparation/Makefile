help:
	@echo 'rotatedir           - Rotate tif files 270 DEG'
	@echo 'analyse.fileprefix  - Analyse $(TIFFSOURCEDIR1)/fileprefix.tif'
	@echo 'analyse             - Analyse $(TIFFSOURCEDIR1)/$(TIFFFILE1)'
	@echo 'writecontent        - Testing reading pixels'
	@echo 'setup               - Create data directories'
	@echo 'clean               - Remove .o, ~ files'

clean:
	rm -r *.o *~

SOURCES=analysetiff.c analysescan.c lblobinfo.c pixelrotation.c analysepage.c lprintbox.c ldraw.c
INCLUDES=analysescan.h tiffcommon.h lblobinfo.h pixelrotation.h analysepage.h analysescantests.c lprintbox.h ldraw.h
OBJECTS=$(SOURCES:%.c=%.o)


RESULTDIR=results

spat:
	g++ -std=c++0x libspat.cc -lspatialindex_c -lspatialindex -o libspat
spat.%:
	g++ -std=c++0x libspat$*.cc -lspatialindex_c -lspatialindex -o libspat$*

compiletiff.%:
	gcc -DHAVE_CONFIG_H  -I/usr/local/include   -g -O2 -Wall -W -MT $*.o -MD -MP -MF .deps/$*.Tpo -c -o $*.o $*.c
	mv -f .deps/**.Tpo .deps/$*.Po
	./libtool  --tag=CC   --mode=link gcc  -g -O2 -Wall -W   -o $* $*.o /usr/local/lib/libtiff.la -lz 
	gcc -g -O2 -Wall -W -o $* $*.o  /usr/local/lib/libtiff.dylib  -lz

otoexe.%:
	./libtool  --tag=CC   --mode=link gcc  -g -O2 -Wall -W   -o $* $*.o /usr/local/lib/libtiff.la -lz 
	gcc -g -O2 -Wall -W -o $* $*.o  /usr/local/lib/libtiff.dylib  -lz


setup:
	mkdir -p .deps
	mkdir -p $(TIFFSOURCEDIR1)
	mkdir -p results

TIFF1=$(HOME)/DeskTop/13trorig/img010.tif

TIFFSOURCEDIR1=corrected13tr
ANALYSETIFF=analysetiff
%.o:	%.c $(INCLUDES)
	rm -f $(ANALYSETIFF)
	gcc -DHAVE_CONFIG_H  -I/usr/local/include   -g -O2 -Wall -W -MT $*.o -MD -MP -MF .deps/$*.Tpo -c -o $*.o $*.c


$(ANALYSETIFF):	$(OBJECTS) 
	./libtool  --tag=CC   --mode=link gcc  -g -O2 -Wall -W   -o $(ANALYSETIFF) $(OBJECTS) /usr/local/lib/libtiff.la -lz 
	gcc -g -O2 -Wall -W -o $(ANALYSETIFF) $(OBJECTS)  /usr/local/lib/libtiff.dylib  -lz


l:
	splint +skip-sys-headers +boolint +matchanyintegral analysetiff.c


#TIFFFILE1=img016_r
TIFFFILE1=img012_r

rect:
	convert $(TIFFSOURCEDIR1)/$(TIFFFILE1).tif -fill none -stroke black -strokewidth 10 -draw "polygon 7746, 990 12942,1002 12926,8524 7730,8512" -draw "polygon 1192,1026 5865,667 6445,8206 1772,8566" o1.tif 

a:
	$(MAKE) $(ANALYSETIFF)
	./$(ANALYSETIFF) -t kin1.tif

t:
	$(MAKE) $(ANALYSETIFF)
	./$(ANALYSETIFF) -t test1.tif

s:
	$(MAKE) $(ANALYSETIFF)
	./$(ANALYSETIFF) -t square1.tif


T:
	$(MAKE) $(ANALYSETIFF)
	./$(ANALYSETIFF) -T

analyse:
	$(MAKE) analyse.$(TIFFFILE1)

analyse.%:
	$(MAKE) $(ANALYSETIFF)
	./$(ANALYSETIFF) -t $(TIFFSOURCEDIR1)/$*.tif

writecontent:
	$(MAKE) $(ANALYSETIFF)
	./$(ANALYSETIFF) -u kin.tif


SOURCEDIR=$(HOME)/DeskTop/13trorig

rotatedir:
	mkdir -p $(TIFFSOURCEDIR1)
	#rm -f $(TIFFSOURCEDIR1)/*_r.tif
	ls $(SOURCEDIR)/*.tif
	for f in $(SOURCEDIR)/*.tif;\
	do \
		echo $$f;\
		f1=$(TIFFSOURCEDIR1)/`basename $$f .tif`_r.tif;\
		echo $$f1;\
		if [ -e $$f1 ];\
		then \
			echo FOUND;\
		else \
			echo MISSING;\
			convert $$f -rotate 270 -compress group4 $$f1;\
		fi \
	done



RESULTSDIR=results
DATA=$(RESULTSDIR)/img016_r_l_y_values.dat
DATAX=$(RESULTSDIR)/img016_r_l_x_values.dat
plot0.%:
	rm -f $(RESULTSDIR)/pixelrotation$*.eps
	echo "set output \"$(RESULTSDIR)/pixelrotation$*.eps\";set terminal postscript color enhanced eps;plot [3300:4400] \"$(DATA)\" using 1:$* title \"pixelrotation $*\";quit" | gnuplot
	open $(RESULTSDIR)/pixelrotation$*.eps

RANGEY=0:8500
FY=a
RANGEX=0:10000
FX=a
plotya.%:
	rm -f $(RESULTSDIR)/pixelrotationya$*.eps
	x=`expr $* + 100 + 2`;echo "set output \"$(RESULTSDIR)/pixelrotationya$*.eps\";set terminal postscript color enhanced eps;plot [$(RANGEY)] \"$(DATA)\" using 1:$$x title \"pixelrotation $*\";quit" | gnuplot
	open $(RESULTSDIR)/pixelrotationy$(FY)$*.eps

plotx.%:
	#rm -f $(RESULTSDIR)/pixelrotationx_$*.eps
	x=`expr $* + 100 + 2`;echo "set output \"$(RESULTSDIR)/pixelrotationx_$*.eps\";set terminal postscript color enhanced eps;plot [$(RANGEX)] \"$(DATAX)\" using 1:$$x title \"pixelrotation x $*\";quit" | gnuplot
	open $(RESULTSDIR)/pixelrotationx_$*.eps


t1:
	convert -size 100x100 canvas: -monochrome -compress group4 c.tif
	convert c.tif -fill none -stroke black -strokewidth 3 -draw "polygon 10,8 62, 8 62,85 10,85" c1.tif
	open c1.tif


t2:
	convert -size 10000x5000 canvas: -monochrome -compress group4 c.tif
	convert c.tif -fill none -stroke black -strokewidth 20 -draw "polygon 200,200 4800, 200 4800,4800 200,4800" c1.tif
	open c1.tif

t3:
	convert -size 14000x10000 canvas: -monochrome -compress group4 c.tif
	convert c.tif -fill none -stroke black -strokewidth 10 -draw "polygon 900,800  6000,800  6000,8000 900,8000" c1.tif
