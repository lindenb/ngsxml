XSLTPROC=xsltproc
XMLLINT=xmllint
.PHONY:all clean tests

all: tests


tests: test01.mk

test01.mk:  test/model01.xml stylesheets/model2make.xsl xsd/schema01.xsd
	-${XMLLINT}  --schema $(filter %.xsd,$^) --noout $<
	${XSLTPROC} --output $@ $(filter %.xsl,$^) $< 
	$(MAKE) -f $@ -I test/config/01/ -n
	

clean:
	
