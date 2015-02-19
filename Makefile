XSLTPROC=xsltproc
XMLLINT=xmllint
.PHONY:all clean tests graph

all: tests


tests: test01.mk
	$(MAKE) -f $< -I test/config/01/

test01.mk:  test/model01.xml stylesheets/model2make.xsl xsd/schema01.xsd 
	cp test/fastq/uncompressed/*.fastq test/fastq
	gzip -f test/fastq/*.fastq
	${XMLLINT}  --schema $(filter %.xsd,$^) --noout $<
	${XSLTPROC}  $(filter %.xsl,$^) $< |\
	awk '/\\$$/ { L=length($$0); printf("%s",substr($$0,1,L-1)); next;} { printf("%s\n",$$0);}' $(addsuffix .tmp,$@) > $@

graph: graph.png 
graph.png : test01.mk
	$(MAKE) -Bdn -f $< -I test/config/01/ | make2graph | dot -Tpng -o $@

clean:
	rm -rf test01.mk OUT test/fastq/*.gz 

	
