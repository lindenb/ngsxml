# Author:
#    Pierre Lindenbaum PhD @yokofakun
# 
#
XSLTPROC=xsltproc
XMLLINT=xmllint


.PHONY:all dryrun clean tests graph

all: tests


tests: test01.mk
	$(MAKE) -f $< -I test/config/01/

dryrun: test01.mk
	$(MAKE) -nd -f $< -I test/config/01/ | grep -E '( is newer than |Must remake)'

# generate Makefile
test01.mk:  test/model01.xml stylesheets/model2make.xsl xsd/schema01.xsd 
	cp test/fastq/uncompressed/*.fastq test/fastq
	gzip -f test/fastq/*.fastq
	#prevent the date of fastqs to be a problem
	touch -t '198001010001'  test/fastq/*.fastq.gz
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
	$(MAKE) --no-builtin-rules -Bdn -f $< -I test/config/01/ |\
		make2graph --basename | dot -Tpng -o $@

clean:
	rm -rf test01.mk OUT test/fastq/*.gz 

	
