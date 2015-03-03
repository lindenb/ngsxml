# Author:
#    Pierre Lindenbaum PhD @yokofakun
# 
#
XSLTPROC=xsltproc
XMLLINT=xmllint


.PHONY:all clean tests graph

all: tests


tests: test01.mk
	$(MAKE) -f $< -I test/config/01/

# generate Makefile
test01.mk:  test/model01.xml stylesheets/model2make.xsl xsd/schema01.xsd 
	cp test/fastq/uncompressed/*.fastq test/fastq
	gzip -f test/fastq/*.fastq
	# validate schema
	${XMLLINT}  \
			--xinclude \
			--path test/ref \
			--schema $(filter %.xsd,$^) --noout $<
	# generate the makefile
	# the awk command is used to remove the line ending with antislash, because our local version of SGE+qmake is buggy.
	${XSLTPROC} \
		--xinclude \
		--path test/ref \
		$(filter %.xsl,$^) $< |\
	awk '/\\$$/ { L=length($$0); printf("%s",substr($$0,1,L-1)); next;} { printf("%s\n",$$0);}' > $(addsuffix .tmp,$@)
	mv $(addsuffix .tmp,$@) $@

# generate graph with makefile2graph
graph: graph.png 
graph.png : test01.mk
	$(MAKE) -Bdn -f $< -I test/config/01/ |\
		make2graph | dot -Tpng -o $@

clean:
	rm -rf test01.mk OUT test/fastq/*.gz 

	
