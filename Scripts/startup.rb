# coding: utf-8
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

sdir = get_global_settings("global.scratch_dir")
if sdir == nil || sdir == ""
  sdir = $home_directory
  if $platform == "win" && sdir =~ / /
    #  Try 8.3 name
	tdir = ENV['TEMP']
	if tdir != nil
	  sdir = tdir.gsub(/\\/, "/")
	end
	if sdir =~ /\/AppData\/Local\/Temp$/
	  sdir = Regexp.last_match.pre_match
	end
  end
  if sdir !~ / /
    set_global_settings("global.scratch_dir", sdir)
  end
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

module Math
  def acos_safe(arg)
    if arg <= -1.0
      return PI
    elsif arg >= 1.0
      return 0.0
    else
      return acos(arg)
    end
  end

  def sqrt_safe(arg)
    arg <= 0.0 ? 0.0 : sqrt(arg)
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
    #  Create a temporary directory like %scratch%/Molby/tag/name.xxxxxx
	base = get_global_settings("global.scratch_dir")
	if base == nil || base == ""
	  base = ask_scratch_dir
	  if base == nil
	    raise "Scratch directory is not set."
	  end
	end
	msg = nil
	name ||= "temp"
	base = base + "/Molby/" + tag
	begin
	  mkdir_recursive(base) rescue ((msg = "Cannot create directory #{base}") && raise)
	  10000.times { |n|
	    p = sprintf("%s/%s.%05d", base, name, n)
	    if !FileTest.exist?(p)
	      Dir.mkdir(p) rescue ((msg = "Cannot create directory #{p}") && raise)
		  return p
	    end
	  }
	rescue
	  raise "Cannot create temporary directory in #{base}"
	end
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
#  def cleanup_temp_dir(path, option = nil)
    #  Clean-up temporary directories
	#  If option is nil, then the directory at path is removed
	#  If option is a number, then the directories that are in the same parent directory as path
	#  are removed except for the newest ones
#	base = get_global_settings("global.scratch_dir")
#	if base == nil || base == ""
#	  raise "Scratch directory is not set."
#	end
#	base = base + "/Molby/"
#	if path.index(base) != 0
#	  raise "Bad cleanup_temp_dir call: the path is not a part of the scratch directory (#{path})"
#	end
#	if !FileTest.directory?(path)
#	  raise "Bad cleanup_temp_dir call: the path is not a directory (#{path})"
#	end
#	option = option.to_i
#	if option <= 0
#	  rm_recursive(path)
#	else
#	  base = File.dirname(path)
#	  ent = Dir.entries.sort_by { |en| File.mtime("#{base}/#{en}").to_i * (-1) } - [".", ".."]
#	  if $platform == "mac"
#	    ent -= [".DS_Store"]
#	  end
#	  ent[0, option] = []  #  Remove newest #{option} entries
#	  #  Mark this directory to be ready to remove (see below)
#	  open("#{path}/.done", "w") { |fp| }
#	  #  Remove the older directories
#	  ent.each { |en|
#	    #  Check the existence of ".done" file (otherwise, we may accidentarily remove a directory
#	    #  that is still in use in other threads)
#	    if File.exist?("#{base}/#{en}/.done")
#	      rm_recursive("#{base}/#{en}")
#		end
#	  }
#	end
# end

  def remove_dir(dir)
    entries = Dir.entries(dir)
	entries.each { |en|
	  next if en == "." || en == ".."
	  fname = "#{dir}/#{en}"
	  if File.directory?(fname)
	    remove_dir(fname)
	  else
	    File.unlink(fname)
	  end
	}
	Dir.unlink(dir)
  end
  
  def erase_old_logs(tdir, level = nil, keep_number = 0)
    log_dir = File.dirname(tdir)
	base = get_global_settings("global.scratch_dir")
	if base == nil || base == ""
	  raise "Scratch directory is not set."
	end
	base = base.gsub(/\\/, "/") + "/Molby/"
	tdir = tdir.gsub(/\\/, "/")
	if tdir.index(base) != 0
	  raise "Bad erase_old_logs call: the path is not a part of the scratch directory (#{tdir}, #{base})"
	end
	if level == nil || level == "none"
	  remove_dir(tdir)
	elsif level == "latest"
	  if keep_number == nil
	    keep_number = 5
	  else
	    keep_number = keep_number.to_i
	  end
	  entries = Dir.entries(log_dir).sort_by { |en| File.mtime("#{log_dir}/#{en}").to_i * (-1) } - [".", ".."]
	  if $platform == "mac"
	    entries -= [".DS_Store"]
	  end
	  entries[0, keep_number] = []  #  Remove newest #{keep_number} entries
	  #  Remove the older directories
	  entries.each { |en|
		#  Check the existence of ".in_use" file (otherwise, we may accidentarily remove a directory
		#  that is still in use in other threads)
		if !File.exist?("#{log_dir}/#{en}/.in_use")
		  rm_recursive("#{log_dir}/#{en}")
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
	else
	  return gets_any_eol
	end
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
