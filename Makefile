XSLTPROC=xsltproc
XMLLINT=xmllint
.PHONY:all clean tests

all: tests


tests: test01.mk
	$(MAKE) -f $< -I test/config/01/

test01.mk:  test/model01.xml stylesheets/model2make.xsl xsd/schema01.xsd 
	cp test/fastq/uncompressed/*.fastq test/fastq
	gzip -f test/fastq/*.fastq
	${XMLLINT}  --schema $(filter %.xsd,$^) --noout $<
	${XSLTPROC} --output $@ $(filter %.xsl,$^) $< 
	 

clean:
	rm -rf test01.mk OUT test/fastq/*.gz 

	
