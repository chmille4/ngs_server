
#dir = Dir.pwd

#Dir.chdir("/Users/chase/Desktop/tmp_workspace/ngs_server/ext/bamtools/build")
exists = `which cmake`
if (exists == '')
  $stderr.puts "\n\n\n[ERROR]: cmake is required and not installed. Get it here: http://www.cmake.org/\n\n"
end

`cmake .`


#Dir.chdir(dir)