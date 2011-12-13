

def ngsadd(source_path)
  source_path = File.absolute_path(source_path)
  file_name = File.basename(source_path)
  data_dir_path = File.join(File.dirname(__FILE__),'../data', file_name)
  `ln -s #{source_path} #{data_dir_path}`    
end

