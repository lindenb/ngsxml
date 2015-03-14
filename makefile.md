Generated Makefile: 
```make

# Name
# 	myProject
# Description:
# 	my project
# 
include config.mk
OUTDIR=OUT
BINDIR=$(abspath ${OUTDIR})/bin


ifeq ($(origin JAVA_HOME), undefined)
JAVA_HOME=/usr
endif
export JAVA_HOME

# if tools are undefined
bwa.version=0.7.12
bwa.exe ?=${BINDIR}/bwa-${bwa.version}/bwa
htslib.version=1.2.1
samtools.version=1.2
samtools.exe ?=${BINDIR}/samtools-${samtools.version}/samtools
bcftools.version=1.2
bcftools.exe ?=${BINDIR}/bcftools-${bcftools.version}/bcftools
tabix.exe ?=$(BINDIR)/htslib-${htslib.version}/tabix
bgzip.exe ?=$(BINDIR)/htslib-${htslib.version}/bgzip
freebayes.version=v0.9.18
freebayes.exe ?=$(BINDIR)/freebayes-${freebayes.version}/bin/freebayes
varscan.version=v2.3.7
varscan.jar=$(BINDIR)/varscan-${varscan.version}/VarScan.${varscan.version}.jar
fastqc.version=0.11.2
fastqc.exe=${BINDIR}/fastqc-${fastqc.version}/fastqc
ant.version=1.8.2
ant.exe=${BINDIR}/apache-ant-${ant-version}/bin/ant
mvn.version=3.2.2
mvn.exe=$(BINDIR)/apache-maven-${mvn.version}/bin/mvn
java.exe ?=$(if ${JAVA_HOME},${JAVA_HOME}/bin/java,java)
picard.version?=1.129
picard.dir=${BINDIR}/picard-tools-${picard.version}
picard.jar=${picard.dir}/picard.jar
gatk.version=3.3
gatk.jar=$(BINDIR)/gatk-protected-${gatk.version}/target/GenomeAnalysisTK.jar
curl.options=-kL
callers= varscan samtools $(sort $(if $(realpath ($(addsuffix /cmake,$(subst :, ,${PATH})))),freebayes,))



define project_dir
$(OUTDIR)/Projects/$(1)
endef

define sample_list
$(call project_dir,$(1))/sample.list
endef


define project_vcf_dir
$(call project_dir,$(1))/VCF
endef

define project_bam_dir
$(call project_dir,$(1))/BAM
endef

# 1 projet
# 2 method
# 3 chrom
# 4 start
# 5 end
define  vcf_segment
$(call project_vcf_dir,$(1))/$(2)/${tmp.prefix}$(1).$(2)_$(3)_$(4)_$(5).vcf.gz
endef

# 1 projet
# 2 method
define vcf_list
$(call project_vcf_dir,$(1))/$(2)/${tmp.prefix}$(1).$(2).vcf.list
endef

# 1 projet
# 2 method
define vcf_final
$(call project_vcf_dir,$(1))/$(2)/$(1).$(2).vcf.gz
endef

# 1 projet
# 2 method
define gather_vcf

$$(call vcf_final,$(1),$(2)) : $$(call vcf_list,$(1),$(2)) ${picard.jar}
	mkdir -p $$(dir $$@) && \
	${java.exe} -jar $$(filter %.jar,$$^) GatherVcfs I=$$< O=$$(addsuffix .tmp.vcf,$$@) && \
	${bgzip.exe} -f $$(addsuffix .tmp.vcf,$$@) && \
	${tabix.exe} -f -p vcf $$(addsuffix .tmp.vcf.gz,$$@) && \
	mv $$(addsuffix .tmp.vcf.gz,$$@) $$@ && \
	mv $$(addsuffix .tmp.vcf.gz.tbi,$$@) $$(addsuffix .tbi,$$@)

endef

# 1 projet
# 2 method
define create_vcf_list


$$(call vcf_list,$(1),$(2)) :  \
	$(call  vcf_segment,$(1),$(2),seg1,0,5101)  \
	$(call  vcf_segment,$(1),$(2),seg2,0,4000) 
	mkdir -p $$(dir $$@)
	rm -f $$(addsuffix .tmp,$$@) 
	echo "$(call  vcf_segment,$(1),$(2),seg1,0,5101)" >> $$(addsuffix .tmp,$$@) 
	echo "$(call  vcf_segment,$(1),$(2),seg2,0,4000)" >> $$(addsuffix .tmp,$$@) 
	mv $$(addsuffix .tmp,$$@) $$@

endef


# 1 target
# 2 fastq files
define run_fastqc

$(1) : $(2) ${fastqc.exe}
	mkdir -p $$(dir $$@) && \
	cat $(2) > $$(addsuffix .tmp.gz,$$@) && \
	${fastqc.exe}   \
		-o $$(dir $$@) \
		-j ${java.exe} \
		 --format fastq  --noextract \
		 $$(addsuffix .tmp.gz,$$@)  && \
	rm $$(addsuffix .tmp.gz,$$@) && \
	mv $$(addsuffix .tmp_fastqc.zip,$$@) $$@
	
endef


.PHONY= all clean all_bams all_vcfs all_quality_controls all_fastqc all_qc_insert_size

all: all_quality_controls all_vcfs 


all_vcfs:  $(foreach method,${callers},$(call vcf_final,Proj1,${method})) 


all_quality_controls : all_fastqc all_qc_insert_size


all_fastqc :  \
	$(call project_dir,Proj1)/Samples/NA12878/FASTQC/NA12878.for_fastq.zip \
	$(call project_dir,Proj1)/Samples/NA12878/FASTQC/NA12878.rev_fastq.zip \
	$(call project_dir,Proj1)/Samples/NA12891/FASTQC/NA12891.for_fastq.zip \
	$(call project_dir,Proj1)/Samples/NA12891/FASTQC/NA12891.rev_fastq.zip \
	$(call project_dir,Proj1)/Samples/NA12892/FASTQC/NA12892.for_fastq.zip \
	$(call project_dir,Proj1)/Samples/NA12892/FASTQC/NA12892.rev_fastq.zip


all_qc_insert_size :  \
	$(call project_dir,Proj1)/Samples/NA12878/SIZE/NA12878.insert_size \
	$(call project_dir,Proj1)/Samples/NA12891/SIZE/NA12891.insert_size \
	$(call project_dir,Proj1)/Samples/NA12892/SIZE/NA12892.insert_size


all_bams:  \
	$(call project_dir,Proj1)/Samples/NA12878/BAM/Proj1_NA12878.bam \
	$(call project_dir,Proj1)/Samples/NA12891/BAM/Proj1_NA12891.bam \
	$(call project_dir,Proj1)/Samples/NA12892/BAM/Proj1_NA12892.bam

$(eval $(foreach C,${callers},$(call gather_vcf,Proj1,$C)))
$(eval $(foreach C,${callers},$(call create_vcf_list,Proj1,$C)))




#
# Create a VCF  for project '' and segment 
#

$(call  vcf_segment,Proj1,samtools,seg1,0,5101) : $(call project_dir,Proj1)/BAM/${tmp.prefix}Proj1.bam.list \
	$(addsuffix .fai,${REFERENCE}) ${samtools.exe}  ${bcftools.exe}
	mkdir -p $(dir $@) && \
	${samtools.exe} mpileup -uf ${REFERENCE} -b $< -r seg1:0-5101 | \
	${bcftools.exe} call  --variants-only --multiallelic-caller --output-type z --output $@


$(call  vcf_segment,Proj1,varscan,seg1,0,5101) : $(call project_dir,Proj1)/BAM/${tmp.prefix}Proj1.bam.list \
	$(addsuffix .fai,${REFERENCE}) ${varscan.jar} ${samtools.exe} $(addsuffix .dict,$(basename ${REFERENCE})) \
	$(call sample_list,Proj1)
	mkdir -p $(dir $@) && \
	${samtools.exe} mpileup -f ${REFERENCE} -b $< -r seg1:0-5101 > $(addsuffix .tmp.mpileup,$@) && \
	${java.exe} -jar $(filter %.jar,$^) mpileup2snp   $(addsuffix .tmp.mpileup,$@) \
		-output-vcf --variants --vcf-sample-list $(call sample_list,Proj1) > $(addsuffix .snp.vcf,$@)  && \
	${java.exe} -jar $(filter %.jar,$^) mpileup2indel $(addsuffix .tmp.mpileup,$@) \
		--output-vcf --variants --vcf-sample-list $(call sample_list,Proj1) > $(addsuffix .indel.vcf,$@)  && \
	head -n 1 $(addsuffix .snp.vcf,$@) > $(addsuffix .tmp.vcf,$@)
	grep -hE '^#' $(addsuffix .snp.vcf,$@) $(addsuffix .indel.vcf,$@) | grep -v "^##fileformat=" | LC_ALL=C sort | uniq >> $(addsuffix .tmp.vcf,$@)
	grep -vhE '^#' $(addsuffix .snp.vcf,$@) $(addsuffix .indel.vcf,$@) | LC_ALL=C sort | uniq >> $(addsuffix .tmp.vcf,$@)
	${java.exe} -jar ${picard.jar} UpdateVcfSequenceDictionary SEQUENCE_DICTIONARY=$(addsuffix .dict,$(basename ${REFERENCE})) I=$(addsuffix .tmp.vcf,$@) O=$(addsuffix .tmp2.vcf,$@)
	${java.exe} -jar ${picard.jar} SortVcf SD=$(addsuffix .dict,$(basename ${REFERENCE})) I=$(addsuffix .tmp2.vcf,$@) O=$(addsuffix .tmp.vcf,$@)
	rm  $(addsuffix .tmp2.vcf,$@) $(addsuffix .snp.vcf,$@) $(addsuffix .indel.vcf,$@) $(addsuffix .tmp.mpileup,$@)
	gzip -f $(addsuffix .tmp.vcf,$@) && mv $(addsuffix .tmp.vcf.gz,$@) $@


$(call  vcf_segment,Proj1,freebayes,seg1,0,5101) : $(call project_dir,Proj1)/BAM/${tmp.prefix}Proj1.bam.list \
	$(addsuffix .fai,${REFERENCE}) ${freebayes.exe}  $(addsuffix .dict,$(basename ${REFERENCE}))
	mkdir -p $(dir $@) && \
	${freebayes.exe} --bam-list $< --vcf $(addsuffix .tmp.vcf,$@) --fasta-reference ${REFERENCE} -r seg1:0-5101
	${java.exe} -jar ${picard.jar} UpdateVcfSequenceDictionary SEQUENCE_DICTIONARY=$(addsuffix .dict,$(basename ${REFERENCE})) I=$(addsuffix .tmp.vcf,$@) O=$(addsuffix .tmp2.vcf,$@)
	rm  $(addsuffix .tmp.vcf,$@)
	gzip -f $(addsuffix .tmp2.vcf,$@) && mv $(addsuffix .tmp2.vcf.gz,$@) $@
	
$(call  vcf_segment,Proj1,unifiedgenotyper,seg1,0,5101) : $(call project_dir,Proj1)/BAM/${tmp.prefix}Proj1.bam.list \
	$(addsuffix .fai,${REFERENCE}) ${gatk.jar}  $(addsuffix .dict,$(basename ${REFERENCE}))
	mkdir -p $(dir $@) && \
	${java.exe} -jar ${gatk.jar} -T UnifiedGenotyper \
	-I $< -R ${REFERENCE}  -glm BOTH -l INFO \
	-L seg1:1-5101 \
	-o $(addsuffix .tmp.vcf,$@)  && \
	gzip -f $(addsuffix .tmp.vcf,$@) && mv $(addsuffix .tmp.vcf.gz,$@) $@



#
# Create a VCF  for project '' and segment 
#

$(call  vcf_segment,Proj1,samtools,seg2,0,4000) : $(call project_dir,Proj1)/BAM/${tmp.prefix}Proj1.bam.list \
	$(addsuffix .fai,${REFERENCE}) ${samtools.exe}  ${bcftools.exe}
	mkdir -p $(dir $@) && \
	${samtools.exe} mpileup -uf ${REFERENCE} -b $< -r seg2:0-4000 | \
	${bcftools.exe} call  --variants-only --multiallelic-caller --output-type z --output $@


$(call  vcf_segment,Proj1,varscan,seg2,0,4000) : $(call project_dir,Proj1)/BAM/${tmp.prefix}Proj1.bam.list \
	$(addsuffix .fai,${REFERENCE}) ${varscan.jar} ${samtools.exe} $(addsuffix .dict,$(basename ${REFERENCE})) \
	$(call sample_list,Proj1)
	mkdir -p $(dir $@) && \
	${samtools.exe} mpileup -f ${REFERENCE} -b $< -r seg2:0-4000 > $(addsuffix .tmp.mpileup,$@) && \
	${java.exe} -jar $(filter %.jar,$^) mpileup2snp   $(addsuffix .tmp.mpileup,$@) \
		-output-vcf --variants --vcf-sample-list $(call sample_list,Proj1) > $(addsuffix .snp.vcf,$@)  && \
	${java.exe} -jar $(filter %.jar,$^) mpileup2indel $(addsuffix .tmp.mpileup,$@) \
		--output-vcf --variants --vcf-sample-list $(call sample_list,Proj1) > $(addsuffix .indel.vcf,$@)  && \
	head -n 1 $(addsuffix .snp.vcf,$@) > $(addsuffix .tmp.vcf,$@)
	grep -hE '^#' $(addsuffix .snp.vcf,$@) $(addsuffix .indel.vcf,$@) | grep -v "^##fileformat=" | LC_ALL=C sort | uniq >> $(addsuffix .tmp.vcf,$@)
	grep -vhE '^#' $(addsuffix .snp.vcf,$@) $(addsuffix .indel.vcf,$@) | LC_ALL=C sort | uniq >> $(addsuffix .tmp.vcf,$@)
	${java.exe} -jar ${picard.jar} UpdateVcfSequenceDictionary SEQUENCE_DICTIONARY=$(addsuffix .dict,$(basename ${REFERENCE})) I=$(addsuffix .tmp.vcf,$@) O=$(addsuffix .tmp2.vcf,$@)
	${java.exe} -jar ${picard.jar} SortVcf SD=$(addsuffix .dict,$(basename ${REFERENCE})) I=$(addsuffix .tmp2.vcf,$@) O=$(addsuffix .tmp.vcf,$@)
	rm  $(addsuffix .tmp2.vcf,$@) $(addsuffix .snp.vcf,$@) $(addsuffix .indel.vcf,$@) $(addsuffix .tmp.mpileup,$@)
	gzip -f $(addsuffix .tmp.vcf,$@) && mv $(addsuffix .tmp.vcf.gz,$@) $@


$(call  vcf_segment,Proj1,freebayes,seg2,0,4000) : $(call project_dir,Proj1)/BAM/${tmp.prefix}Proj1.bam.list \
	$(addsuffix .fai,${REFERENCE}) ${freebayes.exe}  $(addsuffix .dict,$(basename ${REFERENCE}))
	mkdir -p $(dir $@) && \
	${freebayes.exe} --bam-list $< --vcf $(addsuffix .tmp.vcf,$@) --fasta-reference ${REFERENCE} -r seg2:0-4000
	${java.exe} -jar ${picard.jar} UpdateVcfSequenceDictionary SEQUENCE_DICTIONARY=$(addsuffix .dict,$(basename ${REFERENCE})) I=$(addsuffix .tmp.vcf,$@) O=$(addsuffix .tmp2.vcf,$@)
	rm  $(addsuffix .tmp.vcf,$@)
	gzip -f $(addsuffix .tmp2.vcf,$@) && mv $(addsuffix .tmp2.vcf.gz,$@) $@
	
$(call  vcf_segment,Proj1,unifiedgenotyper,seg2,0,4000) : $(call project_dir,Proj1)/BAM/${tmp.prefix}Proj1.bam.list \
	$(addsuffix .fai,${REFERENCE}) ${gatk.jar}  $(addsuffix .dict,$(basename ${REFERENCE}))
	mkdir -p $(dir $@) && \
	${java.exe} -jar ${gatk.jar} -T UnifiedGenotyper \
	-I $< -R ${REFERENCE}  -glm BOTH -l INFO \
	-L seg2:1-4000 \
	-o $(addsuffix .tmp.vcf,$@)  && \
	gzip -f $(addsuffix .tmp.vcf,$@) && mv $(addsuffix .tmp.vcf.gz,$@) $@

	

#
# list of samples, for varscan
#
$(call sample_list,Proj1): 
	mkdir -p $(dir $@)
	rm -f $(addsuffix .tmp,$@) 
	echo "NA12878" >> $(addsuffix .tmp,$@)
	echo "NA12891" >> $(addsuffix .tmp,$@)
	echo "NA12892" >> $(addsuffix .tmp,$@)
	mv $(addsuffix .tmp,$@) $@


#
# create a list of BAMs
#
$(call project_dir,Proj1)/BAM/${tmp.prefix}Proj1.bam.list :  $(addsuffix .bai, $(call project_dir,Proj1)/Samples/NA12878/BAM/Proj1_NA12878.bam $(call project_dir,Proj1)/Samples/NA12891/BAM/Proj1_NA12891.bam $(call project_dir,Proj1)/Samples/NA12892/BAM/Proj1_NA12892.bam)
	mkdir -p $(dir $@)
	rm -f $(addsuffix .tmp,$@) 
	echo "$(call project_dir,Proj1)/Samples/NA12878/BAM/Proj1_NA12878.bam" >> $(addsuffix .tmp,$@) 
	echo "$(call project_dir,Proj1)/Samples/NA12891/BAM/Proj1_NA12891.bam" >> $(addsuffix .tmp,$@) 
	echo "$(call project_dir,Proj1)/Samples/NA12892/BAM/Proj1_NA12892.bam" >> $(addsuffix .tmp,$@) 
	mv $(addsuffix .tmp,$@) $@








#
# index final BAM for Sample 'NA12878'
# 
$(addsuffix .bai, $(call project_dir,Proj1)/Samples/NA12878/BAM/Proj1_NA12878.bam): $(call project_dir,Proj1)/Samples/NA12878/BAM/Proj1_NA12878.bam ${samtools.exe}
	mkdir -p $(dir $@) && \
	${samtools.exe} index  $<
#
# prepare final BAM for Sample 'NA12878'
# 
$(call project_dir,Proj1)/Samples/NA12878/BAM/Proj1_NA12878.bam : $(call project_dir,Proj1)/Samples/NA12878/BAM/${tmp.prefix}rmdup.bam
	mkdir -p $(dir $@) && \
	cp  $< $@


$(call project_dir,Proj1)/Samples/NA12878/BAM/${tmp.prefix}rmdup.bam : $(call project_dir,Proj1)/Samples/NA12878/BAM/${tmp.prefix}merged.bam ${picard.jar}
	mkdir -p $(dir $@) && \
	${java.exe} -jar $(filter %.jar,$^) MarkDuplicates I=$< O=$@ M=$(addsuffix .metrics,$@) AS=true VALIDATION_STRINGENCY=SILENT



#
# merge BAMs 
#
$(call project_dir,Proj1)/Samples/NA12878/BAM/${tmp.prefix}merged.bam :  \
 	$(call project_dir,Proj1)/Samples/NA12878/BAM/${tmp.prefix}1_sorted.bam \
 	$(call project_dir,Proj1)/Samples/NA12878/BAM/${tmp.prefix}2_sorted.bam ${picard.jar}
	mkdir -p $(dir $@) && \
 	${java.exe} -jar $(filter %.jar,$^) MergeSamFiles $(foreach B,$(filter %.bam,$^), I=${B} ) O=$@ AS=true VALIDATION_STRINGENCY=SILENT
  


#
# FastQC
#
$(eval $(call run_fastqc,$(call project_dir,Proj1)/Samples/NA12878/FASTQC/NA12878.for_fastq.zip, test/fastq/NA12878_01_R1.fastq.gz test/fastq/NA12878_02_R1.fastq.gz))
$(eval $(call run_fastqc,$(call project_dir,Proj1)/Samples/NA12878/FASTQC/NA12878.rev_fastq.zip, test/fastq/NA12878_01_R2.fastq.gz test/fastq/NA12878_02_R2.fastq.gz))

#
# Sizes
#
$(call project_dir,Proj1)/Samples/NA12878/SIZE/NA12878.insert_size   : $(call project_dir,Proj1)/Samples/NA12878/BAM/Proj1_NA12878.bam ${picard.jar}
	mkdir -p $(dir $@) && \
 	${java.exe} -jar ${picard.jar} CollectInsertSizeMetrics I=$< O=$@  AS=true H=$(addsuffix .pdf,$@) $(foreach M,ALL_READS SAMPLE LIBRARY READ_GROUP, LEVEL=${M} )




	
#
# Index BAM $(call project_dir,Proj1)/Samples/NA12878/BAM/${tmp.prefix}1_sorted.bam
#
$(addsuffix .bai,$(call project_dir,Proj1)/Samples/NA12878/BAM/${tmp.prefix}1_sorted.bam ): $(call project_dir,Proj1)/Samples/NA12878/BAM/${tmp.prefix}1_sorted.bam ${samtools.exe}
	${samtools} index $<

#
# Align test/fastq/NA12878_01_R1.fastq.gz and test/fastq/NA12878_01_R2.fastq.gz
#
$(call project_dir,Proj1)/Samples/NA12878/BAM/${tmp.prefix}1_sorted.bam : \
	test/fastq/NA12878_01_R1.fastq.gz  \
	test/fastq/NA12878_01_R2.fastq.gz \
	$(addsuffix .bwt,${REFERENCE}) \
	${bwa.exe} ${samtools.exe}
	mkdir -p $(dir $@) && \
	${bwa.exe} mem -R '@RG\tID:idp11858708\tSM:NA12878\tLB:NA12878\tPL:ILLUMINA\tPU:1' \
		${REFERENCE} \
		test/fastq/NA12878_01_R1.fastq.gz \
		test/fastq/NA12878_01_R2.fastq.gz |\
	${samtools.exe} view -uS - |\
	${samtools.exe} sort - $(basename $@) 





	
#
# Index BAM $(call project_dir,Proj1)/Samples/NA12878/BAM/${tmp.prefix}2_sorted.bam
#
$(addsuffix .bai,$(call project_dir,Proj1)/Samples/NA12878/BAM/${tmp.prefix}2_sorted.bam ): $(call project_dir,Proj1)/Samples/NA12878/BAM/${tmp.prefix}2_sorted.bam ${samtools.exe}
	${samtools} index $<

#
# Align test/fastq/NA12878_02_R1.fastq.gz and test/fastq/NA12878_02_R2.fastq.gz
#
$(call project_dir,Proj1)/Samples/NA12878/BAM/${tmp.prefix}2_sorted.bam : \
	test/fastq/NA12878_02_R1.fastq.gz  \
	test/fastq/NA12878_02_R2.fastq.gz \
	$(addsuffix .bwt,${REFERENCE}) \
	${bwa.exe} ${samtools.exe}
	mkdir -p $(dir $@) && \
	${bwa.exe} mem -R '@RG\tID:groupid2\tSM:NA12878\tLB:lib1\tPL:ILMN\tPU:2\tPI:98' \
		${REFERENCE} \
		test/fastq/NA12878_02_R1.fastq.gz \
		test/fastq/NA12878_02_R2.fastq.gz |\
	${samtools.exe} view -uS - |\
	${samtools.exe} sort - $(basename $@) 




#
# index final BAM for Sample 'NA12891'
# 
$(addsuffix .bai, $(call project_dir,Proj1)/Samples/NA12891/BAM/Proj1_NA12891.bam): $(call project_dir,Proj1)/Samples/NA12891/BAM/Proj1_NA12891.bam ${samtools.exe}
	mkdir -p $(dir $@) && \
	${samtools.exe} index  $<
#
# prepare final BAM for Sample 'NA12891'
# 
$(call project_dir,Proj1)/Samples/NA12891/BAM/Proj1_NA12891.bam : $(call project_dir,Proj1)/Samples/NA12891/BAM/${tmp.prefix}rmdup.bam
	mkdir -p $(dir $@) && \
	cp  $< $@


$(call project_dir,Proj1)/Samples/NA12891/BAM/${tmp.prefix}rmdup.bam : $(call project_dir,Proj1)/Samples/NA12891/BAM/${tmp.prefix}merged.bam ${picard.jar}
	mkdir -p $(dir $@) && \
	${java.exe} -jar $(filter %.jar,$^) MarkDuplicates I=$< O=$@ M=$(addsuffix .metrics,$@) AS=true VALIDATION_STRINGENCY=SILENT



#
# merge BAMs 
#
$(call project_dir,Proj1)/Samples/NA12891/BAM/${tmp.prefix}merged.bam :  \
 	$(call project_dir,Proj1)/Samples/NA12891/BAM/${tmp.prefix}1_sorted.bam \
 	$(call project_dir,Proj1)/Samples/NA12891/BAM/${tmp.prefix}2_sorted.bam ${picard.jar}
	mkdir -p $(dir $@) && \
 	${java.exe} -jar $(filter %.jar,$^) MergeSamFiles $(foreach B,$(filter %.bam,$^), I=${B} ) O=$@ AS=true VALIDATION_STRINGENCY=SILENT
  


#
# FastQC
#
$(eval $(call run_fastqc,$(call project_dir,Proj1)/Samples/NA12891/FASTQC/NA12891.for_fastq.zip, test/fastq/NA12891_01_R1.fastq.gz test/fastq/NA12891_02_R1.fastq.gz))
$(eval $(call run_fastqc,$(call project_dir,Proj1)/Samples/NA12891/FASTQC/NA12891.rev_fastq.zip, test/fastq/NA12891_01_R2.fastq.gz test/fastq/NA12891_02_R2.fastq.gz))

#
# Sizes
#
$(call project_dir,Proj1)/Samples/NA12891/SIZE/NA12891.insert_size   : $(call project_dir,Proj1)/Samples/NA12891/BAM/Proj1_NA12891.bam ${picard.jar}
	mkdir -p $(dir $@) && \
 	${java.exe} -jar ${picard.jar} CollectInsertSizeMetrics I=$< O=$@  AS=true H=$(addsuffix .pdf,$@) $(foreach M,ALL_READS SAMPLE LIBRARY READ_GROUP, LEVEL=${M} )




	
#
# Index BAM $(call project_dir,Proj1)/Samples/NA12891/BAM/${tmp.prefix}1_sorted.bam
#
$(addsuffix .bai,$(call project_dir,Proj1)/Samples/NA12891/BAM/${tmp.prefix}1_sorted.bam ): $(call project_dir,Proj1)/Samples/NA12891/BAM/${tmp.prefix}1_sorted.bam ${samtools.exe}
	${samtools} index $<

#
# Align test/fastq/NA12891_01_R1.fastq.gz and test/fastq/NA12891_01_R2.fastq.gz
#
$(call project_dir,Proj1)/Samples/NA12891/BAM/${tmp.prefix}1_sorted.bam : \
	test/fastq/NA12891_01_R1.fastq.gz  \
	test/fastq/NA12891_01_R2.fastq.gz \
	$(addsuffix .bwt,${REFERENCE}) \
	${bwa.exe} ${samtools.exe}
	mkdir -p $(dir $@) && \
	${bwa.exe} mem -R '@RG\tID:idp11860932\tSM:NA12891\tLB:NA12891\tPL:ILLUMINA\tPU:1' \
		${REFERENCE} \
		test/fastq/NA12891_01_R1.fastq.gz \
		test/fastq/NA12891_01_R2.fastq.gz |\
	${samtools.exe} view -uS - |\
	${samtools.exe} sort - $(basename $@) 





	
#
# Index BAM $(call project_dir,Proj1)/Samples/NA12891/BAM/${tmp.prefix}2_sorted.bam
#
$(addsuffix .bai,$(call project_dir,Proj1)/Samples/NA12891/BAM/${tmp.prefix}2_sorted.bam ): $(call project_dir,Proj1)/Samples/NA12891/BAM/${tmp.prefix}2_sorted.bam ${samtools.exe}
	${samtools} index $<

#
# Align test/fastq/NA12891_02_R1.fastq.gz and test/fastq/NA12891_02_R2.fastq.gz
#
$(call project_dir,Proj1)/Samples/NA12891/BAM/${tmp.prefix}2_sorted.bam : \
	test/fastq/NA12891_02_R1.fastq.gz  \
	test/fastq/NA12891_02_R2.fastq.gz \
	$(addsuffix .bwt,${REFERENCE}) \
	${bwa.exe} ${samtools.exe}
	mkdir -p $(dir $@) && \
	${bwa.exe} mem -R '@RG\tID:idp11861588\tSM:NA12891\tLB:NA12891\tPL:ILLUMINA\tPU:1' \
		${REFERENCE} \
		test/fastq/NA12891_02_R1.fastq.gz \
		test/fastq/NA12891_02_R2.fastq.gz |\
	${samtools.exe} view -uS - |\
	${samtools.exe} sort - $(basename $@) 




#
# index final BAM for Sample 'NA12892'
# 
$(addsuffix .bai, $(call project_dir,Proj1)/Samples/NA12892/BAM/Proj1_NA12892.bam): $(call project_dir,Proj1)/Samples/NA12892/BAM/Proj1_NA12892.bam ${samtools.exe}
	mkdir -p $(dir $@) && \
	${samtools.exe} index  $<
#
# prepare final BAM for Sample 'NA12892'
# 
$(call project_dir,Proj1)/Samples/NA12892/BAM/Proj1_NA12892.bam : $(call project_dir,Proj1)/Samples/NA12892/BAM/${tmp.prefix}rmdup.bam
	mkdir -p $(dir $@) && \
	cp  $< $@


$(call project_dir,Proj1)/Samples/NA12892/BAM/${tmp.prefix}rmdup.bam : $(call project_dir,Proj1)/Samples/NA12892/BAM/${tmp.prefix}merged.bam ${picard.jar}
	mkdir -p $(dir $@) && \
	${java.exe} -jar $(filter %.jar,$^) MarkDuplicates I=$< O=$@ M=$(addsuffix .metrics,$@) AS=true VALIDATION_STRINGENCY=SILENT



#
# merge BAMs 
#
$(call project_dir,Proj1)/Samples/NA12892/BAM/${tmp.prefix}merged.bam :  \
 	$(call project_dir,Proj1)/Samples/NA12892/BAM/${tmp.prefix}1_sorted.bam \
 	$(call project_dir,Proj1)/Samples/NA12892/BAM/${tmp.prefix}2_sorted.bam ${picard.jar}
	mkdir -p $(dir $@) && \
 	${java.exe} -jar $(filter %.jar,$^) MergeSamFiles $(foreach B,$(filter %.bam,$^), I=${B} ) O=$@ AS=true VALIDATION_STRINGENCY=SILENT
  


#
# FastQC
#
$(eval $(call run_fastqc,$(call project_dir,Proj1)/Samples/NA12892/FASTQC/NA12892.for_fastq.zip, test/fastq/NA12892_01_R1.fastq.gz test/fastq/NA12892_02_R1.fastq.gz))
$(eval $(call run_fastqc,$(call project_dir,Proj1)/Samples/NA12892/FASTQC/NA12892.rev_fastq.zip, test/fastq/NA12892_01_R2.fastq.gz test/fastq/NA12892_02_R2.fastq.gz))

#
# Sizes
#
$(call project_dir,Proj1)/Samples/NA12892/SIZE/NA12892.insert_size   : $(call project_dir,Proj1)/Samples/NA12892/BAM/Proj1_NA12892.bam ${picard.jar}
	mkdir -p $(dir $@) && \
 	${java.exe} -jar ${picard.jar} CollectInsertSizeMetrics I=$< O=$@  AS=true H=$(addsuffix .pdf,$@) $(foreach M,ALL_READS SAMPLE LIBRARY READ_GROUP, LEVEL=${M} )




	
#
# Index BAM $(call project_dir,Proj1)/Samples/NA12892/BAM/${tmp.prefix}1_sorted.bam
#
$(addsuffix .bai,$(call project_dir,Proj1)/Samples/NA12892/BAM/${tmp.prefix}1_sorted.bam ): $(call project_dir,Proj1)/Samples/NA12892/BAM/${tmp.prefix}1_sorted.bam ${samtools.exe}
	${samtools} index $<

#
# Align test/fastq/NA12892_01_R1.fastq.gz and test/fastq/NA12892_01_R2.fastq.gz
#
$(call project_dir,Proj1)/Samples/NA12892/BAM/${tmp.prefix}1_sorted.bam : \
	test/fastq/NA12892_01_R1.fastq.gz  \
	test/fastq/NA12892_01_R2.fastq.gz \
	$(addsuffix .bwt,${REFERENCE}) \
	${bwa.exe} ${samtools.exe}
	mkdir -p $(dir $@) && \
	${bwa.exe} mem -R '@RG\tID:idp11862556\tSM:NA12892\tLB:NA12892\tPL:ILLUMINA\tPU:1' \
		${REFERENCE} \
		test/fastq/NA12892_01_R1.fastq.gz \
		test/fastq/NA12892_01_R2.fastq.gz |\
	${samtools.exe} view -uS - |\
	${samtools.exe} sort - $(basename $@) 





	
#
# Index BAM $(call project_dir,Proj1)/Samples/NA12892/BAM/${tmp.prefix}2_sorted.bam
#
$(addsuffix .bai,$(call project_dir,Proj1)/Samples/NA12892/BAM/${tmp.prefix}2_sorted.bam ): $(call project_dir,Proj1)/Samples/NA12892/BAM/${tmp.prefix}2_sorted.bam ${samtools.exe}
	${samtools} index $<

#
# Align test/fastq/NA12892_02_R1.fastq.gz and test/fastq/NA12892_02_R2.fastq.gz
#
$(call project_dir,Proj1)/Samples/NA12892/BAM/${tmp.prefix}2_sorted.bam : \
	test/fastq/NA12892_02_R1.fastq.gz  \
	test/fastq/NA12892_02_R2.fastq.gz \
	$(addsuffix .bwt,${REFERENCE}) \
	${bwa.exe} ${samtools.exe}
	mkdir -p $(dir $@) && \
	${bwa.exe} mem -R '@RG\tID:idp11863212\tSM:NA12892\tLB:NA12892\tPL:ILLUMINA\tPU:1' \
		${REFERENCE} \
		test/fastq/NA12892_02_R1.fastq.gz \
		test/fastq/NA12892_02_R2.fastq.gz |\
	${samtools.exe} view -uS - |\
	${samtools.exe} sort - $(basename $@) 




$(addsuffix .fai,${REFERENCE}): ${REFERENCE} ${samtools.exe}
	${samtools.exe} faidx $<

$(addsuffix .bwt,${REFERENCE}): ${REFERENCE} ${bwa.exe}
	${bwa.exe} index $<

$(addsuffix .dict,$(basename ${REFERENCE})): ${REFERENCE} ${picard.jar}
	mkdir -p $(dir $@)  && \
	rm -f $@  && \
	${java.exe} -jar ${picard.jar} CreateSequenceDictionary R=$< O=$@



${bwa.exe}   :	
	rm -rf $(dir $@) && \
	mkdir -p $(BINDIR) && \
	curl ${curl.options} -L -o $(BINDIR)/bwa-${bwa.version}.zip -L "https://github.com/lh3/bwa/archive/${bwa.version}.zip" && \
	unzip $(BINDIR)/bwa-${bwa.version}.zip -d $(BINDIR)  && \
	rm $(BINDIR)/bwa-${bwa.version}.zip && \
	make -C $(dir $@) || make -C $(dir $@) CFLAGS=" -g -Wall -Wno-unused-function -O2 -msse2 "

${bcftools.exe}   : ${BINDIR}/htslib/htslib.mk
	rm -rf $(dir $@) && \
	mkdir -p $(BINDIR) && \
	curl ${curl.options} -L -o $(BINDIR)/bcftools-${bcftools.version}.zip -L "https://github.com/samtools/bcftools/archive/${bcftools.version}.zip" && \
	unzip $(BINDIR)/bcftools-${bcftools.version}.zip -d $(BINDIR)  && \
	rm $(BINDIR)/bcftools-${bcftools.version}.zip && \
	make -C $(dir $@)

${samtools.exe}  :	${BINDIR}/htslib/htslib.mk
	rm -rf $(dir $@) && \
	mkdir -p $(BINDIR) && \
	curl ${curl.options} -L -o $(BINDIR)/samtools-${samtools.version}.zip -L "https://github.com/samtools/samtools/archive/${samtools.version}.zip" && \
	unzip $(BINDIR)/samtools-${samtools.version}.zip -d $(BINDIR)  && \
	rm $(BINDIR)/samtools-${samtools.version}.zip && \
	make -C $(dir $@)

#
# Create symlink because samtools requires a folder named htslib
#
${BINDIR}/htslib/htslib.mk : ${BINDIR}/htslib-${htslib.version}/htslib.mk
	ln -s $(dir $<) ${BINDIR}/htslib

#
# Tabix and bgzip are compiled by htslib
#
${tabix.exe} ${bgzip.exe}: ${BINDIR}/htslib-${htslib.version}/htslib.mk 

#
# Download htslib ${htslib.version}
#
${BINDIR}/htslib-${htslib.version}/htslib.mk :
	rm -rf $(dir $@) && \
	mkdir -p $(BINDIR) && \
	curl ${curl.options} -L -o $(BINDIR)/htslib-${htslib.version}.zip -L "https://github.com/samtools/htslib/archive/${htslib.version}.zip" && \
	unzip $(BINDIR)/htslib-${htslib.version}.zip -d $(BINDIR)  && \
	rm $(BINDIR)/htslib-${htslib.version}.zip && \
	make -C $(dir $@)

#
# Download Fastqc ${fastqc.version}
#
${fastqc.exe} :
	echo "Downloading FASTQC v.${fastqc.version}"  && \
	rm -rf $(dir $@) && \
	mkdir -p $(BINDIR) && \
	curl ${curl.options} -L -o $(BINDIR)/fastqc-${fastqc.version}.zip -L "http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${fastqc.version}.zip" && \
	unzip $(BINDIR)/fastqc-${fastqc.version}.zip -d ${BINDIR}  && \
	mv ${BINDIR}/FastQC $(BINDIR)/fastqc-${fastqc.version} && \
	rm $(BINDIR)/fastqc-${fastqc.version}.zip && \
	chmod +x ${fastqc.exe}


#
# Download picard
#
${picard.jar} :
	echo "DOWNLOADING PICARD version : ${picard.version}"
	rm -rf $(dir $@) && \
	mkdir -p $(BINDIR) && \
	curl ${curl.options} -o ${BINDIR}/picard-tools-${picard.version}.zip -kL "https://github.com/broadinstitute/picard/releases/download/${picard.version}/picard-tools-${picard.version}.zip" && \
	unzip ${BINDIR}/picard-tools-${picard.version}.zip -d ${BINDIR} && \
	rm $(BINDIR)/picard-tools-${picard.version}.zip

#
# Varscan
#
${varscan.jar} :
	echo "DOWNLOADING varscan version : ${varscan.version}"
	rm -rf $(dir $@) && \
	mkdir -p $(dir $@) && \
	curl ${curl.options} -o $@ "http://heanet.dl.sourceforge.net/project/varscan/VarScan.${varscan.version}.jar"

#
# Freebayes
#
${freebayes.exe} :
	echo "DOWNLOADING Freebayes version : ${freebayes.version}"	
	rm -rf $(BINDIR)/freebayes-${freebayes.version} && \
	git clone --recursive "https://github.com/ekg/freebayes.git" $(BINDIR)/freebayes-${freebayes.version}
	(cd $(BINDIR)/freebayes-${freebayes.version} && git checkout ${freebayes.version} && git submodule update --recursive && make )

#
# Ant
#
${ant.exe}:
	echo "DOWNLOADING Ant version : ${ant.version}"	&& \
	mkdir -p ${BINDIR} && \
	rm -rf $(BINDIR)/apache-ant-${ant.version} && \
	curl ${curl.options} -o ant-${ant.version}.zip "http://archive.apache.org/dist/ant/binaries/apache-ant-${ant.version}-bin.zip" && \
	unzip ant-${ant.version}.zip -d ${BINDIR} && \
	rm  ant-${ant.version}.zip && \
	touch --no-create $@

#
# Maven
#
${mvn.exe}:
	echo "DOWNLOADING Maven version : ${mvn.version}" && \
	mkdir -p ${BINDIR} && \
	rm -rf $(BINDIR)/apache-maven-${mvn.version} && \
	curl ${curl.options} -o maven-${mvn.version}.zip "http://archive.apache.org/dist/maven/binaries/apache-maven-${mvn.version}-bin.zip" && \
	unzip maven-${mvn.version}.zip -d ${BINDIR} && \
	rm  maven-${mvn.version}.zip && \
	touch --no-create $@

#
# GATK
# at home I need to run maven twice because sometimes I've got
# a error "could not create JVM"
# 
${gatk.jar} : ${mvn.exe}
	echo "DOWNLOADING GATK version : ${gatk.version}"
	mkdir -p ${BINDIR}/m2 && \
	rm -rf $(BINDIR)/gatk-protected-${gatk.version} gatk-${gatk.version}.zip && \
	curl ${curl.options} -o gatk-${gatk.version}.zip \
		"https://github.com/broadgsa/gatk-protected/archive/${gatk.version}.zip" && \
	unzip gatk-${gatk.version}.zip -d ${BINDIR} && \
	rm gatk-${gatk.version}.zip && \
	(cd $(BINDIR)/gatk-protected-${gatk.version} && ${mvn.exe} verify  -Dmaven.repo.local=../m2)
		
	

clean:
	rm -rf ${BINDIR}

```
