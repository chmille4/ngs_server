require'sinatra/base'

class Myserver < Sinatra::Base  
  #set :run, true
  set :public, '.'
 # set :port, opts[:port]
 set :port, 4569
#  set :daemonize, opts[:daemonize]
#  set :daemonize, true

  gempath = File.join(File.dirname(__FILE__), "../")
  datapath = File.join(gempath, 'data')
  bamtools_path = "#{gempath}/ext/bamtools/bin/bamtools"
  vcftools_path = "#{gempath}/ext/vcftools/bin/vcf-query"


  get '/json/bam/*' do |path|

    content_type :json
    response['Access-Control-Allow-Origin'] = '*';

    # invoke with eg: base_url/json/bam/file=subset22-sorted.bam?min=30000000&max=30010000&segment=22   
    json = `#{bamtools_path} convert -in #{datapath}/#{path} -format json -region #{params["segment"]}:#{params["min"]}..#{params["max"]}`
    json       

  end

  get '/json/vcf/*' do |path|

    content_type :json
    response['Access-Control-Allow-Origin'] = '*';

    # invoke with eg: base_url/json/vcf/file=ALL.2of4intersection.20100804.sites.vcf.gz?min=6992179&max=6992190&segment=1
    json = `#{vcftools_path} #{datapath}/#{path} #{params["segment"]}:#{params["min"]}-#{params["max"]} -f '{"reference":"%CHROM","position":%POS,"ref":"%REF","alt":"%ALT","info":%INFO},'`
    json

  end

  get '/json/sources/:extension' do

    content_type :json
    response['Access-Control-Allow-Origin'] = '*';

    # invoke with eg: base_url/json/vcf/file=ALL.2of4intersection.20100804.sites.vcf.gz?min=6992179&max=6992190&segment=1
    list = `find -L #{datapath} -name *#{params['extension']} | awk -F/ '{ for(i=2; i<NF; i++) {printf("%s/", $i)} print $NF }'`
    list

  end

  
end
