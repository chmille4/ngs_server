
#dir = Dir.pwd

#Dir.chdir("/Users/chase/Desktop/tmp_workspace/ngs_server/ext/bamtools/build")

file = File.new('MakeFile', 'w')

file. puts "ifndef PREFIX
    export PREFIX = $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
endif
export BINDIR = ${PREFIX}/bin
export MODDIR = ${PREFIX}/lib

DIRS = cpp perl
install:
	    @mkdir -p $(BINDIR); mkdir -p $(MODDIR); \
        for dir in $(DIRS); do cd $$dir && $(MAKE) $(MAKEFLAGS) && cd ..; done

clean:
		@for dir in $(DIRS); do cd $$dir && $(MAKE) clean && cd ..; done"



#Dir.chdir(dir)