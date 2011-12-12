
BIN = vcf-compare fill-aa vcf-annotate vcf-merge vcf-isec vcf-stats vcf-to-tab fill-an-ac \
    vcf-query vcf-convert vcf-subset vcf-validator vcf-concat vcf-sort
MOD = FaSlice.pm Vcf.pm

install:
	    @for i in $(BIN); do cp $(CURDIR)/$$i $(BINDIR)/$$i; done; \
        for i in $(MOD); do cp $(CURDIR)/$$i $(MODDIR)/$$i; done; 

clean:
		@for i in $(BIN); do rm -f $(BINDIR)/$$i; done; \
        for i in $(MOD); do rm -f $(MODDIR)/$$i; done; 
