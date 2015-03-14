<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet
	version="1.0"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	>

<xsl:output method='text' />

<xsl:template match='/'>
  <xsl:apply-templates select="model"/>
</xsl:template>

<xsl:template match='model'>
# Name
# 	<xsl:value-of select="@name"/>
# Description:
# 	<xsl:value-of select="@description"/>
# 
include config.mk
OUTDIR=<xsl:value-of select="@directory"/>
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
	mkdir -p $$(dir $$@) &amp;&amp; \
	${java.exe} -jar $$(filter %.jar,$$^) GatherVcfs I=$$&lt; O=$$(addsuffix .tmp.vcf,$$@) &amp;&amp; \
	${bgzip.exe} -f $$(addsuffix .tmp.vcf,$$@) &amp;&amp; \
	${tabix.exe} -f -p vcf $$(addsuffix .tmp.vcf.gz,$$@) &amp;&amp; \
	mv $$(addsuffix .tmp.vcf.gz,$$@) $$@ &amp;&amp; \
	mv $$(addsuffix .tmp.vcf.gz.tbi,$$@) $$(addsuffix .tbi,$$@)

endef

# 1 projet
# 2 method
define create_vcf_list


$$(call vcf_list,$(1),$(2)) : <xsl:for-each select="/model/segments/segment"> \
	$(call  vcf_segment,$(1),$(2),<xsl:value-of select="@chrom"/>,<xsl:value-of select="@start"/>,<xsl:value-of select="@end"/>) </xsl:for-each>
	mkdir -p $$(dir $$@)
	rm -f $$(addsuffix .tmp,$$@) <xsl:for-each select="/model/segments/segment">
	echo "$(call  vcf_segment,$(1),$(2),<xsl:value-of select="@chrom"/>,<xsl:value-of select="@start"/>,<xsl:value-of select="@end"/>)" &gt;&gt; $$(addsuffix .tmp,$$@) </xsl:for-each>
	mv $$(addsuffix .tmp,$$@) $$@

endef


# 1 target
# 2 fastq files
define run_fastqc

$(1) : $(2) ${fastqc.exe}
	mkdir -p $$(dir $$@) &amp;&amp; \
	cat $(2) > $$(addsuffix .tmp.gz,$$@) &amp;&amp; \
	${fastqc.exe}   \
		-o $$(dir $$@) \
		-j ${java.exe} \
		 --format fastq  --noextract \
		 $$(addsuffix .tmp.gz,$$@)  &amp;&amp; \
	rm $$(addsuffix .tmp.gz,$$@) &amp;&amp; \
	mv $$(addsuffix .tmp_fastqc.zip,$$@) $$@
	
endef


.PHONY= all clean all_bams all_vcfs all_quality_controls all_fastqc all_qc_insert_size

all: all_quality_controls all_vcfs 


all_vcfs: <xsl:for-each select="project"> $(foreach method,${callers},$(call vcf_final,<xsl:value-of select="@name"/>,${method})) </xsl:for-each>


all_quality_controls : all_fastqc all_qc_insert_size


all_fastqc : <xsl:for-each select="project/sample"> \
	<xsl:apply-templates select="." mode="fastqc.for"/> \
	<xsl:apply-templates select="." mode="fastqc.rev"/></xsl:for-each>


all_qc_insert_size : <xsl:for-each select="project/sample"> \
	<xsl:apply-templates select="." mode="qc.insert.size"/></xsl:for-each>


all_bams: <xsl:for-each select="project/sample"> \
	<xsl:apply-templates select="." mode="bam.final"/></xsl:for-each>


<xsl:apply-templates select="project"/>

$(addsuffix .fai,${REFERENCE}): ${REFERENCE} ${samtools.exe}
	${samtools.exe} faidx $&lt;

$(addsuffix .bwt,${REFERENCE}): ${REFERENCE} ${bwa.exe}
	${bwa.exe} index $&lt;

$(addsuffix .dict,$(basename ${REFERENCE})): ${REFERENCE} ${picard.jar}
	mkdir -p $(dir $@)  &amp;&amp; \
	rm -f $@  &amp;&amp; \
	${java.exe} -jar ${picard.jar} CreateSequenceDictionary R=$&lt; O=$@



${bwa.exe}   :	
	rm -rf $(dir $@) &amp;&amp; \
	mkdir -p $(BINDIR) &amp;&amp; \
	curl ${curl.options} -L -o $(BINDIR)/bwa-${bwa.version}.zip -L "https://github.com/lh3/bwa/archive/${bwa.version}.zip" &amp;&amp; \
	unzip $(BINDIR)/bwa-${bwa.version}.zip -d $(BINDIR)  &amp;&amp; \
	rm $(BINDIR)/bwa-${bwa.version}.zip &amp;&amp; \
	make -C $(dir $@) || make -C $(dir $@) CFLAGS=" -g -Wall -Wno-unused-function -O2 -msse2 "

${bcftools.exe}   : ${BINDIR}/htslib/htslib.mk
	rm -rf $(dir $@) &amp;&amp; \
	mkdir -p $(BINDIR) &amp;&amp; \
	curl ${curl.options} -L -o $(BINDIR)/bcftools-${bcftools.version}.zip -L "https://github.com/samtools/bcftools/archive/${bcftools.version}.zip" &amp;&amp; \
	unzip $(BINDIR)/bcftools-${bcftools.version}.zip -d $(BINDIR)  &amp;&amp; \
	rm $(BINDIR)/bcftools-${bcftools.version}.zip &amp;&amp; \
	make -C $(dir $@)

${samtools.exe}  :	${BINDIR}/htslib/htslib.mk
	rm -rf $(dir $@) &amp;&amp; \
	mkdir -p $(BINDIR) &amp;&amp; \
	curl ${curl.options} -L -o $(BINDIR)/samtools-${samtools.version}.zip -L "https://github.com/samtools/samtools/archive/${samtools.version}.zip" &amp;&amp; \
	unzip $(BINDIR)/samtools-${samtools.version}.zip -d $(BINDIR)  &amp;&amp; \
	rm $(BINDIR)/samtools-${samtools.version}.zip &amp;&amp; \
	make -C $(dir $@)

#
# Create symlink because samtools requires a folder named htslib
#
${BINDIR}/htslib/htslib.mk : ${BINDIR}/htslib-${htslib.version}/htslib.mk
	ln -s $(dir $&lt;) ${BINDIR}/htslib

#
# Tabix and bgzip are compiled by htslib
#
${tabix.exe} ${bgzip.exe}: ${BINDIR}/htslib-${htslib.version}/htslib.mk 

#
# Download htslib ${htslib.version}
#
${BINDIR}/htslib-${htslib.version}/htslib.mk :
	rm -rf $(dir $@) &amp;&amp; \
	mkdir -p $(BINDIR) &amp;&amp; \
	curl ${curl.options} -L -o $(BINDIR)/htslib-${htslib.version}.zip -L "https://github.com/samtools/htslib/archive/${htslib.version}.zip" &amp;&amp; \
	unzip $(BINDIR)/htslib-${htslib.version}.zip -d $(BINDIR)  &amp;&amp; \
	rm $(BINDIR)/htslib-${htslib.version}.zip &amp;&amp; \
	make -C $(dir $@)

#
# Download Fastqc ${fastqc.version}
#
${fastqc.exe} :
	echo "Downloading FASTQC v.${fastqc.version}"  &amp;&amp; \
	rm -rf $(dir $@) &amp;&amp; \
	mkdir -p $(BINDIR) &amp;&amp; \
	curl ${curl.options} -L -o $(BINDIR)/fastqc-${fastqc.version}.zip -L "http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${fastqc.version}.zip" &amp;&amp; \
	unzip $(BINDIR)/fastqc-${fastqc.version}.zip -d ${BINDIR}  &amp;&amp; \
	mv ${BINDIR}/FastQC $(BINDIR)/fastqc-${fastqc.version} &amp;&amp; \
	rm $(BINDIR)/fastqc-${fastqc.version}.zip &amp;&amp; \
	chmod +x ${fastqc.exe}


#
# Download picard
#
${picard.jar} :
	echo "DOWNLOADING PICARD version : ${picard.version}"
	rm -rf $(dir $@) &amp;&amp; \
	mkdir -p $(BINDIR) &amp;&amp; \
	curl ${curl.options} -o ${BINDIR}/picard-tools-${picard.version}.zip -kL "https://github.com/broadinstitute/picard/releases/download/${picard.version}/picard-tools-${picard.version}.zip" &amp;&amp; \
	unzip ${BINDIR}/picard-tools-${picard.version}.zip -d ${BINDIR} &amp;&amp; \
	rm $(BINDIR)/picard-tools-${picard.version}.zip

#
# Varscan
#
${varscan.jar} :
	echo "DOWNLOADING varscan version : ${varscan.version}"
	rm -rf $(dir $@) &amp;&amp; \
	mkdir -p $(dir $@) &amp;&amp; \
	curl ${curl.options} -o $@ "http://heanet.dl.sourceforge.net/project/varscan/VarScan.${varscan.version}.jar"

#
# Freebayes
#
${freebayes.exe} :
	echo "DOWNLOADING Freebayes version : ${freebayes.version}"	
	rm -rf $(BINDIR)/freebayes-${freebayes.version} &amp;&amp; \
	git clone --recursive "https://github.com/ekg/freebayes.git" $(BINDIR)/freebayes-${freebayes.version}
	(cd $(BINDIR)/freebayes-${freebayes.version} &amp;&amp; git checkout ${freebayes.version} &amp;&amp; git submodule update --recursive &amp;&amp; make )

#
# Ant
#
${ant.exe}:
	echo "DOWNLOADING Ant version : ${ant.version}"	&amp;&amp; \
	mkdir -p ${BINDIR} &amp;&amp; \
	rm -rf $(BINDIR)/apache-ant-${ant.version} &amp;&amp; \
	curl ${curl.options} -o ant-${ant.version}.zip "http://archive.apache.org/dist/ant/binaries/apache-ant-${ant.version}-bin.zip" &amp;&amp; \
	unzip ant-${ant.version}.zip -d ${BINDIR} &amp;&amp; \
	rm  ant-${ant.version}.zip &amp;&amp; \
	touch --no-create $@

#
# Maven
#
${mvn.exe}:
	echo "DOWNLOADING Maven version : ${mvn.version}" &amp;&amp; \
	mkdir -p ${BINDIR} &amp;&amp; \
	rm -rf $(BINDIR)/apache-maven-${mvn.version} &amp;&amp; \
	curl ${curl.options} -o maven-${mvn.version}.zip "http://archive.apache.org/dist/maven/binaries/apache-maven-${mvn.version}-bin.zip" &amp;&amp; \
	unzip maven-${mvn.version}.zip -d ${BINDIR} &amp;&amp; \
	rm  maven-${mvn.version}.zip &amp;&amp; \
	touch --no-create $@

#
# GATK
# at home I need to run maven twice because sometimes I've got
# a error "could not create JVM"
# 
${gatk.jar} : ${mvn.exe}
	echo "DOWNLOADING GATK version : ${gatk.version}"
	mkdir -p ${BINDIR}/m2 &amp;&amp; \
	rm -rf $(BINDIR)/gatk-protected-${gatk.version} gatk-${gatk.version}.zip &amp;&amp; \
	curl ${curl.options} -o gatk-${gatk.version}.zip \
		"https://github.com/broadgsa/gatk-protected/archive/${gatk.version}.zip" &amp;&amp; \
	unzip gatk-${gatk.version}.zip -d ${BINDIR} &amp;&amp; \
	rm gatk-${gatk.version}.zip &amp;&amp; \
	(cd $(BINDIR)/gatk-protected-${gatk.version} &amp;&amp; ${mvn.exe} verify  -Dmaven.repo.local=../m2)
		
	

clean:
	rm -rf ${BINDIR}

</xsl:template>


<xsl:template match='project'>
<xsl:variable name="proj" select="."/>

$(eval $(foreach C,${callers},$(call gather_vcf,<xsl:value-of select="@name"/>,$C)))
$(eval $(foreach C,${callers},$(call create_vcf_list,<xsl:value-of select="@name"/>,$C)))


<xsl:for-each select="/model/segments/segment">

#
# Create a VCF  for project '<xsl:value-of select="@name"/>' and segment <xsl:apply-templates select="." mode="id"/>
#

$(call  vcf_segment,<xsl:value-of select="$proj/@name"/>,samtools,<xsl:value-of select="@chrom"/>,<xsl:value-of select="@start"/>,<xsl:value-of select="@end"/>) : <xsl:apply-templates select="$proj" mode="bam.list"/> \
	$(addsuffix .fai,${REFERENCE}) ${samtools.exe}  ${bcftools.exe}
	mkdir -p $(dir $@) &amp;&amp; \
	${samtools.exe} mpileup -uf ${REFERENCE} -b $&lt; -r <xsl:value-of select="concat(@chrom,':',@start,'-',@end)"/> | \
	${bcftools.exe} call  --variants-only --multiallelic-caller --output-type z --output $@


$(call  vcf_segment,<xsl:value-of select="$proj/@name"/>,varscan,<xsl:value-of select="@chrom"/>,<xsl:value-of select="@start"/>,<xsl:value-of select="@end"/>) : <xsl:apply-templates select="$proj" mode="bam.list"/> \
	$(addsuffix .fai,${REFERENCE}) ${varscan.jar} ${samtools.exe} $(addsuffix .dict,$(basename ${REFERENCE})) \
	$(call sample_list,<xsl:value-of select="$proj/@name"/>)
	mkdir -p $(dir $@) &amp;&amp; \
	${samtools.exe} mpileup -f ${REFERENCE} -b $&lt; -r <xsl:value-of select="concat(@chrom,':',@start,'-',@end)"/> &gt; $(addsuffix .tmp.mpileup,$@) &amp;&amp; \
	${java.exe} -jar $(filter %.jar,$^) mpileup2snp   $(addsuffix .tmp.mpileup,$@) \
		-output-vcf --variants --vcf-sample-list $(call sample_list,<xsl:value-of select="$proj/@name"/>) &gt; $(addsuffix .snp.vcf,$@)  &amp;&amp; \
	${java.exe} -jar $(filter %.jar,$^) mpileup2indel $(addsuffix .tmp.mpileup,$@) \
		--output-vcf --variants --vcf-sample-list $(call sample_list,<xsl:value-of select="$proj/@name"/>) &gt; $(addsuffix .indel.vcf,$@)  &amp;&amp; \
	head -n 1 $(addsuffix .snp.vcf,$@) &gt; $(addsuffix .tmp.vcf,$@)
	grep -hE '^#' $(addsuffix .snp.vcf,$@) $(addsuffix .indel.vcf,$@) | grep -v "^##fileformat=" | LC_ALL=C sort | uniq &gt;&gt; $(addsuffix .tmp.vcf,$@)
	grep -vhE '^#' $(addsuffix .snp.vcf,$@) $(addsuffix .indel.vcf,$@) | LC_ALL=C sort | uniq &gt;&gt; $(addsuffix .tmp.vcf,$@)
	${java.exe} -jar ${picard.jar} UpdateVcfSequenceDictionary SEQUENCE_DICTIONARY=$(addsuffix .dict,$(basename ${REFERENCE})) I=$(addsuffix .tmp.vcf,$@) O=$(addsuffix .tmp2.vcf,$@)
	${java.exe} -jar ${picard.jar} SortVcf SD=$(addsuffix .dict,$(basename ${REFERENCE})) I=$(addsuffix .tmp2.vcf,$@) O=$(addsuffix .tmp.vcf,$@)
	rm  $(addsuffix .tmp2.vcf,$@) $(addsuffix .snp.vcf,$@) $(addsuffix .indel.vcf,$@) $(addsuffix .tmp.mpileup,$@)
	gzip -f $(addsuffix .tmp.vcf,$@) &amp;&amp; mv $(addsuffix .tmp.vcf.gz,$@) $@


$(call  vcf_segment,<xsl:value-of select="$proj/@name"/>,freebayes,<xsl:value-of select="@chrom"/>,<xsl:value-of select="@start"/>,<xsl:value-of select="@end"/>) : <xsl:apply-templates select="$proj" mode="bam.list"/> \
	$(addsuffix .fai,${REFERENCE}) ${freebayes.exe}  $(addsuffix .dict,$(basename ${REFERENCE}))
	mkdir -p $(dir $@) &amp;&amp; \
	${freebayes.exe} --bam-list $&lt; --vcf $(addsuffix .tmp.vcf,$@) --fasta-reference ${REFERENCE} -r <xsl:value-of select="concat(@chrom,':',@start,'-',@end)"/>
	${java.exe} -jar ${picard.jar} UpdateVcfSequenceDictionary SEQUENCE_DICTIONARY=$(addsuffix .dict,$(basename ${REFERENCE})) I=$(addsuffix .tmp.vcf,$@) O=$(addsuffix .tmp2.vcf,$@)
	rm  $(addsuffix .tmp.vcf,$@)
	gzip -f $(addsuffix .tmp2.vcf,$@) &amp;&amp; mv $(addsuffix .tmp2.vcf.gz,$@) $@
	
$(call  vcf_segment,<xsl:value-of select="$proj/@name"/>,unifiedgenotyper,<xsl:value-of select="@chrom"/>,<xsl:value-of select="@start"/>,<xsl:value-of select="@end"/>) : <xsl:apply-templates select="$proj" mode="bam.list"/> \
	$(addsuffix .fai,${REFERENCE}) ${gatk.jar}  $(addsuffix .dict,$(basename ${REFERENCE}))
	mkdir -p $(dir $@) &amp;&amp; \
	${java.exe} -jar ${gatk.jar} -T UnifiedGenotyper \
	-I $&lt; -R ${REFERENCE}  -glm BOTH -l INFO \
	-L <xsl:value-of select="concat(@chrom,':',@start+1,'-',@end)"/> \
	-o $(addsuffix .tmp.vcf,$@)  &amp;&amp; \
	gzip -f $(addsuffix .tmp.vcf,$@) &amp;&amp; mv $(addsuffix .tmp.vcf.gz,$@) $@

</xsl:for-each>	

#
# list of samples, for varscan
#
$(call sample_list,<xsl:value-of select="@name"/>): 
	mkdir -p $(dir $@)
	rm -f $(addsuffix .tmp,$@) <xsl:for-each select="sample">
	echo "<xsl:value-of select="@name"/>" &gt;&gt; $(addsuffix .tmp,$@)</xsl:for-each>
	mv $(addsuffix .tmp,$@) $@


#
# create a list of BAMs
#
<xsl:apply-templates select="." mode="bam.list"/> :  $(addsuffix .bai,<xsl:for-each select="sample"><xsl:text> </xsl:text><xsl:apply-templates select="." mode="bam.final"/> </xsl:for-each>)
	mkdir -p $(dir $@)
	rm -f $(addsuffix .tmp,$@) <xsl:for-each select="sample">
	echo "<xsl:apply-templates select="." mode="bam.final"/>" &gt;&gt; $(addsuffix .tmp,$@) </xsl:for-each>
	mv $(addsuffix .tmp,$@) $@






<xsl:apply-templates select="sample"/>

</xsl:template>


<xsl:template match='sample'>

#
# index final BAM for Sample '<xsl:value-of select="@name"/>'
# 
$(addsuffix .bai, <xsl:apply-templates select="." mode="bam.final"/>): <xsl:apply-templates select="." mode="bam.final"/> ${samtools.exe}
	mkdir -p $(dir $@) &amp;&amp; \
	${samtools.exe} index  $&lt;
#
# prepare final BAM for Sample '<xsl:value-of select="@name"/>'
# 
<xsl:apply-templates select="." mode="bam.final"/><xsl:text> : </xsl:text><xsl:apply-templates select="." mode="bam.rmdup"/>
	mkdir -p $(dir $@) &amp;&amp; \
	cp  $&lt; $@


<xsl:apply-templates select="." mode="bam.rmdup"/><xsl:text> : </xsl:text><xsl:apply-templates select="." mode="bam.merged"/> ${picard.jar}
	mkdir -p $(dir $@) &amp;&amp; \
	${java.exe} -jar $(filter %.jar,$^) MarkDuplicates I=$&lt; O=$@ M=$(addsuffix .metrics,$@) AS=true VALIDATION_STRINGENCY=SILENT


<xsl:if test="count(fastq) &gt; 1">
#
# merge BAMs 
#
<xsl:apply-templates select="." mode="bam.merged"/><xsl:text> : </xsl:text><xsl:for-each select="fastq"> \
 	<xsl:apply-templates select="." mode="bam.sorted"/> </xsl:for-each> ${picard.jar}
	mkdir -p $(dir $@) &amp;&amp; \
 	${java.exe} -jar $(filter %.jar,$^) MergeSamFiles $(foreach B,$(filter %.bam,$^), I=${B} ) O=$@ AS=true VALIDATION_STRINGENCY=SILENT
  
</xsl:if>

#
# FastQC
#
$(eval $(call run_fastqc,<xsl:apply-templates select="." mode="fastqc.for"/>,<xsl:for-each select="fastq/for"><xsl:text> </xsl:text><xsl:apply-templates select="."/></xsl:for-each>))
$(eval $(call run_fastqc,<xsl:apply-templates select="." mode="fastqc.rev"/>,<xsl:for-each select="fastq/rev"><xsl:text> </xsl:text><xsl:apply-templates select="."/></xsl:for-each>))

#
# Sizes
#
<xsl:apply-templates select="." mode="qc.insert.size"/>   : <xsl:apply-templates select="." mode="bam.final"/> ${picard.jar}
	mkdir -p $(dir $@) &amp;&amp; \
 	${java.exe} -jar ${picard.jar} CollectInsertSizeMetrics I=$&lt; O=$@  AS=true H=$(addsuffix .pdf,$@) $(foreach M,ALL_READS SAMPLE LIBRARY READ_GROUP, LEVEL=${M} )

<xsl:apply-templates select="fastq"/>

</xsl:template>


<xsl:template match='fastq'>


	
#
# Index BAM <xsl:apply-templates select="." mode="bam.sorted"/>
#
<xsl:text>$(addsuffix .bai,</xsl:text><xsl:apply-templates select="." mode="bam.sorted"/><xsl:text> ): </xsl:text><xsl:apply-templates select="." mode="bam.sorted"/> ${samtools.exe}
	${samtools} index $&lt;

#
# Align <xsl:value-of select="for"/> and <xsl:value-of select="rev"/>
#
<xsl:apply-templates select="." mode="bam.sorted"/> : \
	<xsl:apply-templates select="for"/>  \
	<xsl:apply-templates select="rev"/> \
	$(addsuffix .bwt,${REFERENCE}) \
	${bwa.exe} ${samtools.exe}
	mkdir -p $(dir $@) &amp;&amp; \
	${bwa.exe} mem -R '@RG\tID:<xsl:apply-templates select="." mode="fastq.id"/>\tSM:<xsl:value-of select="../@name"/>\tLB:<xsl:choose>
		<xsl:when test="@library"><xsl:value-of select="@library"/></xsl:when>
		<xsl:otherwise><xsl:value-of select="../@name"/></xsl:otherwise>
		</xsl:choose>\tPL:<xsl:choose>
		<xsl:when test="@platform"><xsl:value-of select="@platform"/></xsl:when>
		<xsl:otherwise><xsl:text>ILLUMINA</xsl:text></xsl:otherwise>
		</xsl:choose>\tPU:<xsl:choose>
		<xsl:when test="@lane"><xsl:value-of select="@lane"/></xsl:when>
		<xsl:otherwise><xsl:text>1</xsl:text></xsl:otherwise>
		</xsl:choose><xsl:choose>
		<xsl:when test="@median-size">\tPI:<xsl:value-of select="@median-size"/></xsl:when>
		</xsl:choose>' \
		${REFERENCE} \
		<xsl:apply-templates select="for"/> \
		<xsl:apply-templates select="rev"/> |\
	${samtools.exe} view -uS - |\
	${samtools.exe} sort - $(basename $@) 


</xsl:template>

<!-- BAM mapped to reference and sorted -->
<xsl:template match='for|rev'>

<xsl:value-of select="normalize-space(./text())"/>

</xsl:template>



<!-- BAM mapped to reference and sorted -->
<xsl:template match='fastq' mode="bam.sorted">
<xsl:apply-templates select="." mode="bam.dir"/>
<xsl:value-of select="concat('${tmp.prefix}',count(preceding-sibling::fastq) + 1,'_sorted.bam')"/>
</xsl:template>

<!-- directory for the BAMs -->
<xsl:template match='fastq' mode="bam.dir">
<xsl:apply-templates select=".." mode="bam.dir"/>
</xsl:template>

<!-- fastq-id will be used for read-groups in BAM -->
<xsl:template match='fastq' mode="fastq.id">
<xsl:choose>
<xsl:when test="@id"><xsl:value-of select="@id"/></xsl:when>
<xsl:otherwise><xsl:value-of select="generate-id(.)"/></xsl:otherwise>
</xsl:choose>
</xsl:template>


<!-- directory for the BAMs -->
<xsl:template match='sample' mode="bam.dir">
<xsl:apply-templates select="." mode="dir"/>
<xsl:text>/BAM/</xsl:text>
</xsl:template>

<!--Base directory for a sample -->
<xsl:template match='sample' mode="dir">
<xsl:apply-templates select=".." mode="dir"/>
<xsl:text>/Samples/</xsl:text>
<xsl:value-of select="@name"/>
</xsl:template>


<!--Base directory for a sample /FASTQC-->
<xsl:template match='sample' mode="fastqc.dir">
<xsl:apply-templates select="." mode="dir"/>
<xsl:text>/FASTQC/</xsl:text>
</xsl:template>

<!--Base directory for a sample /FASTQC-->
<xsl:template match='sample' mode="fastqc.for">
<xsl:apply-templates select="." mode="fastqc.dir"/>
<xsl:value-of select="@name"/>
<xsl:text>.for_fastq.zip</xsl:text>
</xsl:template>

<!--Base directory for a sample /FASTQC-->
<xsl:template match='sample' mode="fastqc.rev">
<xsl:apply-templates select="." mode="fastqc.dir"/>
<xsl:value-of select="@name"/>
<xsl:text>.rev_fastq.zip</xsl:text>
</xsl:template>


<!-- Print the name BAM with merged sorted BAM -->
<xsl:template match='sample' mode="bam.merged">
<xsl:choose>
  <xsl:when test="count(fastq)=1">
  	<xsl:apply-templates select="fastq" mode="bam.sorted"/>
  </xsl:when>
  <xsl:otherwise>
    <xsl:apply-templates select="." mode="bam.dir"/>
    <xsl:text>${tmp.prefix}merged.bam</xsl:text>
  </xsl:otherwise>
</xsl:choose>
</xsl:template>

<!-- Print the name BAM with marked duplicates for a sample -->
<xsl:template match='sample' mode="bam.rmdup">
  <xsl:apply-templates select="." mode="bam.dir"/>
  <xsl:text>${tmp.prefix}rmdup.bam</xsl:text>
</xsl:template>

<!-- Print the name of final BAM for a sample -->
<xsl:template match='sample' mode="bam.final">
  <xsl:apply-templates select="." mode="bam.dir"/>
  <xsl:value-of select="concat(../@name,'_',@name,'.bam')"/>
</xsl:template>


<xsl:template match='sample' mode="qc.insert.size">
  <xsl:apply-templates select="." mode="dir"/>
  <xsl:value-of select="concat('/SIZE/',@name,'.insert_size')"/>
</xsl:template>





<!-- list of BAM for this project-->
<xsl:template match='project' mode="bam.list">
  <xsl:apply-templates select="." mode="dir"/>
  <xsl:text>/BAM/${tmp.prefix}</xsl:text>
  <xsl:value-of select="concat(@name,'.bam.list')"/>
</xsl:template>



<!-- Print the PROJECT directory -->
<xsl:template match='project' mode="dir">
<xsl:text>$(call project_dir,</xsl:text>
<xsl:value-of select="@name"/>
<xsl:text>)</xsl:text>
</xsl:template>


</xsl:stylesheet>
