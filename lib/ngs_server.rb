require 'rubygems'
require 'sinatra'
require 'sinatra/base'
require 'json'



class MyNgsServer < Sinatra::Base

  use Rack::Deflater
  set :server, 'thin'

  gempath = File.join(File.dirname(__FILE__), "../")
  datapath = File.join(gempath, 'data')
  bamtools_path = "#{gempath}/ext/bamtools/bin/bamtools"
  vcftools_path = "#{gempath}/ext/vcftools/bin/vcf-query"
  
  get '/json/bam/*' do |path|
    
    # turn off json content type so a preflight cors doesn't need to be done
    #content_type :json
    response['Access-Control-Allow-Origin'] = '*'; 
    

    # invoke with eg: base_url/json/bam/subset22-sorted.bam?min=30000000&max=30010000&segment=22   
    json = `#{bamtools_path} convert -in #{datapath}/#{path} -format json -region #{params["segment"]}:#{params["min"]}..#{params["max"]}`
    json       

  end

  get '/json/vcf/*' do |path|
    
    # turn off json content type so a preflight cors doesn't need to be done
    #content_type :json
    response['Access-Control-Allow-Origin'] = '*';

    # invoke with eg: base_url/json/vcf/ALL.2of4intersection.20100804.sites.vcf.gz?min=6992179&max=6992190&segment=1
    json = `#{vcftools_path} #{datapath}/#{path} #{params["segment"]}:#{params["min"]}-#{params["max"]} -f '{"reference":"%CHROM","position":%POS,"ref":"%REF","alt":"%ALT","info":%INFO},'`
    json

  end

  get '/json/sources/*' do |path|

    # turn off json content type so a preflight cors doesn't need to be done
    # content_type :json
    response['Access-Control-Allow-Origin'] = '*';

    extension = File.basename(path)
    dirpath = File.dirname(path)    
    
    private_dir = "private"

    # invoke with eg: base_url/json/vcf/file=ALL.2of4intersection.20100804.sites.vcf.gz?min=6992179&max=6992190&segment=1
    list  = `find -L #{datapath}/#{dirpath} -not \\( -name private -prune \\) -name '*#{extension}' | awk -F#{datapath}/ '{print $2}'`
    list.split("\n").to_json

  end

end  

