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
  s.description = %q{"TODO"}

  s.rubyforge_project = "ngs_server"

  s.files         = `git ls-files`.split("\n")
  s.test_files    = `git ls-files -- {test,spec,features}/*`.split("\n")
  s.executables   = `git ls-files -- bin/*`.split("\n").map{ |f| File.basename(f) }
  s.require_paths = ["lib"]
end
