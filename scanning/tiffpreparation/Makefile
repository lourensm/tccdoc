help:
	@echo 'rotatedir           - Rotate tif files 270 DEG'
	@echo 'analyse.fileprefix  - Analyse $(TIFFSOURCEDIR1)/fileprefix.tif'
	@echo 'analyse             - Analyse $(TIFFSOURCEDIR1)/$(TIFFFILE1)'
	@echo 'writecontent        - Testing reading pixels'
	@echo 'setup               - Create data directories'
	@echo 'clean               - Remove .o, ~ files'

clean:
	rm -r *.o *~

RESULTDIR=results

compiletiff.%:
	gcc -DHAVE_CONFIG_H  -I/usr/local/include   -g -O2 -Wall -W -MT $*.o -MD -MP -MF .deps/$*.Tpo -c -o $*.o $*.c
	mv -f .deps/**.Tpo .deps/$*.Po
	./libtool  --tag=CC   --mode=link gcc  -g -O2 -Wall -W   -o $* $*.o /usr/local/lib/libtiff.la -lz 
	gcc -g -O2 -Wall -W -o $* $*.o  /usr/local/lib/libtiff.dylib  -lz


setup:
	mkdir -p .deps
	mkdir -p $(TIFFSOURCEDIR1)
	mkdir -p results

TIFF1=$(HOME)/DeskTop/13trorig/img010.tif

TIFFSOURCEDIR1=corrected13tr
ANALYSETIFF=analysetiff
$(ANALYSETIFF):	$(ANALYSETIFF).c
	$(MAKE) compiletiff.$(ANALYSETIFF)


TIFFFILE1=img016_r
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
	rm -f $(TIFFSOURCEDIR1)/*_r.tif
	ls $(SOURCEDIR)/*.tif
	for f in $(SOURCEDIR)/*.tif;\
	do \
		echo $$f;\
		f1=$(TIFFSOURCEDIR1)/`basename $$f .tif`_r.tif;\
		convert $$f -rotate 270 -compress group4 $$f1;\
	done

#convert $$f -rotate 270 -compress group4 img011rc.tif;\


DATA=$(RESULTSDIR)/pixelrotation.dat
plot0.%:
	rm -f $(RESULTSDIR)/pixelrotation$*.eps
	echo "set output \"$(RESULTSDIR)/pixelrotation$*.eps\";set terminal postscript color enhanced eps;plot [3900:4400] \"$(DATA)\" using 1:$* title \"pixelrotation $*\";quit" | gnuplot
	open pixelrotation$*.eps
plot1.%:
	rm -f pixelrotation1_$*.eps
	echo "set output \"pixelrotation1_$*.eps\";set terminal postscript color enhanced eps;plot [7880:7990] \"$(DATA)\" using 1:$* title \"pixelrotation $*\";quit" | gnuplot
	open pixelrotation1_$*.eps
