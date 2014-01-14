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
	$KCODE="SJIS"
	$home_directory = ENV['USERPROFILE'].gsub(/\\/, "/")
  when /darwin/
    $platform = "mac"
	$KCODE="UTF8"
	$home_directory = ENV['HOME']
  else
    $platform = "other"
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
	    if c == 13
		  #  \r or \r\n
		  if (c = getc) && c != 10
		    ungetc(c)
		    @end_of_line = "\r"
		  else
		    @end_of_line = "\r\n"
		  end
		  break
		elsif c == 10
		  #  \n
		  @end_of_line = "\n"
		  break
		else
		  s += c.chr
		end
	  end
	  if @end_of_line
	    s += $/
	  end
    end
	return s
  end
end

load "transform.rb"
load "molecule.rb"
load "loadsave.rb"
load "formula.rb"
load "dialog.rb"
load "commands.rb"
load "md.rb"
load "uff.rb"
load "gamess.rb"
load "crystal.rb"
load "mopac6.rb"

GC.start
