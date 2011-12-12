
FILES=docs index license options perl_module specs links

install:
	    @for i in $(FILES); do echo "php src/index.php $$i > $$i.html"; php src/index.php $$i > $$i.html; done;

