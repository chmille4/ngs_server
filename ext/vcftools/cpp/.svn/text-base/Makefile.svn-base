# Make file for vcftools
# Author: Adam Auton
# ($Revision: 230 $)

# Compiler
CPP = g++
# Output executable
EXECUTABLE = vcftools
# Flag used to turn on compilation of PCA routines
ifndef VCFTOOLS_PCA
	VCFTOOLS_PCA = 0
endif
# Compiler flags
CPPFLAGS = -O2 -Wall -Wextra
#CPPFLAGS = -g 
# Included libraries (zlib)
LIB = -lz 
#LIB = -lz -I/opt/local/include/ -L/opt/local/lib/
SOURCES = vcftools.cpp vcf_file.cpp vcf_entry.cpp \
			vcf_entry_getters.cpp vcf_entry_setters.cpp \
			vcf_file_filters.cpp vcf_file_output.cpp \
			vcf_file_format_convert.cpp \
			vcf_file_diff.cpp parameters.cpp \
			vcf_file_index.cpp \
			output_log.cpp

ifeq ($(VCFTOOLS_PCA), 1)
	# Define flag for PCA routine compilation
	CPPFLAGS += -DVCFTOOLS_PCA
	# Add LAPACK library
	LIB += -llapack	
	# Add PCA source code
	SOURCES += dgeev.cpp
endif

all: vcftools

vcftools: $(SOURCES)
	$(CPP) $(CPPFLAGS) $(SOURCES) -o $@ $(LIB)
ifdef BINDIR
	cp $(CURDIR)/$@ $(BINDIR)/$@
endif

clean:
	@rm -f vcftools
	@rm -f $(BINDIR)/vcftools
