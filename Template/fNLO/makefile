FILES = Cards HTML SubProcesses bin lib Source Events FixedOrderAnalysis \
	   index.html README README.lhapdf TemplateVersion.txt MGMEVersion.txt
  
amcatnlo.tar.gz: amcatnlo.tar
	rm -f amcatnlo.tar.gz
	gzip amcatnlo.tar
  
amcatnlo.tar: 
	tar cf amcatnlo.tar $(FILES)
  
clean:
	rm -f *.log
	cd Source; make clean
	rm -f SubProcesses/*.o
	for dir in `ls -d SubProcesses/P*`; do cd $$dir; make clean; cd ../../;\
	    done
