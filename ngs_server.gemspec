# -*- encoding: utf-8 -*-
$:.push File.expand_path("../lib", __FILE__)
require "ngs_server/version"

Gem::Specification.new do |s|
  s.name        = "ngs_server"
  s.version     = NgsServer::VERSION
  s.platform    = Gem::Platform::RUBY
  s.authors     = ["Chase Miller"]
  s.email       = ["chmille4@gmail.com"]
  s.homepage    = ""
  s.summary     = %q{"Ultra Lightweight NGS Data Server"}
  s.description = %q{"Converts BAM/VCF files into JSON for consumption by web apps"}
  s.add_dependency('sinatra', '>= 1.2')
  s.add_dependency('thin', '>= 1.2')
  s.add_dependency('rack', '>= 1')

  s.rubyforge_project = "ngs_server"

  s.files         = `git ls-files`.split("\n")
  s.test_files    = `git ls-files -- {test,spec,features}/*`.split("\n")
  s.executables   = `git ls-files -- bin/*`.split("\n").map{ |f| File.basename(f) }
#  s.executables   = "bin/ngs_"
  s.extensions    = ["ext/bamtools/extconf.rb", "ext/vcftools/extconf.rb"]
  s.require_paths = ["lib"]
end
