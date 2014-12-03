XSLTPROC=xsltproc
XMLLINT=xmllint
.PHONY:all clean tests

all: tests


tests: test01.mk
	$(MAKE) -f $< -I test/config/01/ -n

test01.mk:  test/model01.xml stylesheets/model2make.xsl xsd/schema01.xsd \
	test/ref/ref.fa \
	test/fastq/sample_1_01_R2.fastq.gz test/fastq/sample_1_02_R2.fastq.gz \
	test/fastq/sample_2_01_R2.fastq.gz
	-${XMLLINT}  --schema $(filter %.xsd,$^) --noout $<
	${XSLTPROC} --output $@ $(filter %.xsl,$^) $< 
	

test/ref/ref.fa:
	mkdir -p $(dir $@) && \
	rm -f $@  \
	$(foreach  C,chr4_gl000194_random chr1_gl000192_random, && curl "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/$(C).fa.gz" | gunzip -c >> $@ ) 



test/fastq/sample_1_01_R2.fastq.gz: test/fastq/sample_1_01_R1.fastq.gz
test/fastq/sample_1_01_R1.fastq.gz: test/ref/ref.fa
	mkdir -p $(dir $@) && \
	wgsim -N 1000 test/ref/ref.fa test/fastq/sample_1_01_R1.fastq test/fastq/sample_1_01_R2.fastq
	gzip -f test/fastq/sample_1_01_R1.fastq test/fastq/sample_1_01_R2.fastq

test/fastq/sample_1_02_R2.fastq.gz: test/fastq/sample_1_02_R1.fastq.gz
test/fastq/sample_1_02_R1.fastq.gz: test/ref/ref.fa
	mkdir -p $(dir $@) && \
	wgsim -N 1000 test/ref/ref.fa test/fastq/sample_1_02_R1.fastq test/fastq/sample_1_02_R2.fastq
	gzip -f test/fastq/sample_1_02_R1.fastq test/fastq/sample_1_02_R2.fastq

test/fastq/sample_2_01_R2.fastq.gz: test/fastq/sample_2_01_R1.fastq.gz
test/fastq/sample_2_01_R1.fastq.gz: test/ref/ref.fa
	mkdir -p $(dir $@) && \
	wgsim -N 1000 test/ref/ref.fa test/fastq/sample_2_01_R1.fastq test/fastq/sample_2_01_R2.fastq
	gzip -f test/fastq/sample_2_01_R1.fastq test/fastq/sample_2_01_R2.fastq

clean:
	rm -f test01.mk OUT test/fastq test/ref
	
