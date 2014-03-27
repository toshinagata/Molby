#
#  startup.rb
#
#  Created by Toshi Nagata.
#  Copyright 2008 Toshi Nagata. All rights reserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation version 2 of the License.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

include Molby
include Math

Deg2Rad = Math::PI / 180.0
Rad2Deg = 180.0 / Math::PI

$startup_dir = Dir.pwd
case RUBY_PLATFORM
  when /mswin|mingw|cygwin|bccwin/
    $platform = "win"
	Encoding.default_external = "shift_jis"
	$home_directory = ENV['USERPROFILE'].gsub(/\\/, "/")
  when /darwin/
    $platform = "mac"
	Encoding.default_external = "utf-8"
	$home_directory = ENV['HOME']
  else
    $platform = "other"
	Encoding.default_external = "locale"
	$home_directory = ENV['HOME']
end

$backtrace = nil

def backtrace
  if $backtrace
    print $backtrace.join("\n")
  end
end

#  Utility methods
module Enumerable
  def sum(&block)
    if block
      self.inject(0) { |sum, v| sum + block.call(v) }
	else
	  self.inject(0) { |sum, v| sum + v }
	end
  end
  def average(&block)
    sum(&block) / Float(self.length)
  end
end

module Kernel
  def filecopy(src, dst)
    fpin = File.open(src, "rb")
    return nil if fpin == nil
    fpout = File.open(dst, "wb")
    if fpout == nil
      fpin.close
      return nil
    end
    a = ""
    while fpin.read(4096, a)
      fpout.write(a)
    end
    fpin.close
    fpout.close
    return true
  end
  def mkdir_recursive(path)
    if FileTest.directory?(path)
	  return 0
    else
	  dir = File.dirname(path)
	  if !FileTest.exist?(dir)
	    mkdir_recursive(dir)
	  end
      return Dir.mkdir(path)
	end
  end
  def create_temp_dir(tag, name = nil)
    #  Create a temporary directory like HomeDirectory/Molby/tag/name.xxxxxx
	name ||= "temp"
	base = $home_directory + "/Molby/" + tag
	mkdir_recursive(base)
	10000.times { |n|
	  p = sprintf("%s/%s.%05d", base, name, n)
	  if !FileTest.exist?(p)
	    Dir.mkdir(p)
		return p
	  end
	}
	raise "Cannot create temporary directory in #{base}"
  end
  def rm_recursive(path)
    if FileTest.directory?(path)
	  Dir.foreach(path) { |en|
	    next if en == "." || en == ".."
		rm_recursive(path + "/" + en)
	  }
	  Dir.rmdir(path)
	else
	  File.unlink(path)
	end
  end
  def cleanup_temp_dir(path, option = nil)
    #  Clean-up temporary directories
	#  If option is nil, then the directory at path is removed
	#  If option is a number, then the directories that are in the same parent directory as path
	#  are removed except for the newest ones
	if path.index($home_directory + "/Molby/") != 0
	  raise "Bad cleanup_temp_dir call: the path does not begin with $HOME/Molby/ (#{path})"
	end
	if !FileTest.directory?(path)
	  raise "Bad cleanup_temp_dir call: the path is not a directory (#{path})"
	end
	option = option.to_i
	if option <= 0
	  rm_recursive(path)
	else
	  base = File.dirname(path)
	  ent = Dir.entries.sort_by { |en| File.mtime("#{base}/#{en}").to_i * (-1) } - [".", ".."]
	  if $platform == "mac"
	    ent -= [".DS_Store"]
	  end
	  ent[0, option] = []  #  Remove newest #{option} entries
	  #  Mark this directory to be ready to remove (see below)
	  open("#{path}/.done", "w") { |fp| }
	  #  Remove the older directories
	  ent.each { |en|
	    #  Check the existence of ".done" file (otherwise, we may accidentarily remove a directory
	    #  that is still in use in other threads)
	    if File.exist?("#{base}/#{en}/.done")
	      rm_recursive("#{base}/#{en}")
		end
	  }
	end
  end
end

class IO
  alias :gets_original :gets
  def gets(rs = $/)
    if rs != $/
	  return gets_original(rs)
	end
    if @end_of_line
	  s = gets_original(@end_of_line)
	  if s && s.chomp!(@end_of_line)
	    s += $/
	  end
	else
	  s = ""
	  while c = getc
	    if c == "\r"
		  #  \r or \r\n
		  if (c = getc) && c != "\n"
		    ungetc(c)
		    @end_of_line = "\r"
		  else
		    @end_of_line = "\r\n"
		  end
		  break
		elsif c == "\n"
		  #  \n
		  @end_of_line = "\n"
		  break
		else
		  s += c
		end
	  end
	  if @end_of_line
	    s += $/
	  end
    end
	return s
  end
end

#  Additional method definitions
load "transform.rb"
load "molecule.rb"
load "loadsave.rb"
load "formula.rb"
load "dialog.rb"

#  Menu commands
load "view.rb"
load "uff.rb"
load "md.rb"
load "mopac6.rb"
load "gamess.rb"
load "crystal.rb"
load "commands.rb"

GC.start
