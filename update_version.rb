#!/usr/bin/ruby

require 'kconv'

#  Get the version string
#  version = "X.X.X"
version = nil
date = nil
eval IO.read("Version")
ver = version
t = Time.now
year = t.year
month = t.month
day = t.day
d = sprintf("%04d%02d%02d", year, month, day)
# exit 0 if date == d
#File.open("Version", "w") { |fp|
#  fp.print "version = \"#{version}\"\n"
#  fp.print "date = \"#{d}\"\n"
#}
build = "build " + d
# verstr = "v#{ver} #{build}"
verstr = "v#{ver}"
yrange = (year > 2008 ? "2008-#{year}" : "2008")

def modify_file(name, &block)
  ary = IO.readlines(name) rescue return
  modified = false
  ary.each_with_index { |s, i|
    s1 = block.call(s.dup)
    if s1 && s1 != s
      ary[i] = s1
      modified = true
    end
  }
  if modified
    File.rename(name, name + "~")
    open(name, "wb") { |fp|
      ary.each { |s| fp.write(s) }
    }
    File.delete(name + "~")
  end
end

#  Modify Info.plist
nm = "build-xcode/Molby-Info.plist"
version = false
modify_file(nm) { |s|
  if s =~ /Copyright/
    s.sub(/[-0-9]+ Toshi Nagata/, "#{yrange} Toshi Nagata")
  elsif s =~ /Version \d+\.\d+/
    "\t<string>Version #{ver}</string>\n"
  elsif version
    version = false
    "\t<string>#{verstr}</string>\n"
  else
    version = (s =~ /\bCFBundleVersion\b/)
    nil
  end
}

#  Modify_MacLegacy Info.plist
nm = "build-xcode/Molby_MacLegacy-Info.plist"
version = false
modify_file(nm) { |s|
  if s =~ /Copyright/
    s.sub(/[-0-9]+ Toshi Nagata/, "#{yrange} Toshi Nagata")
  elsif s =~ /Version \d+\.\d+/
    "\t<string>Version #{ver}</string>\n"
  elsif version
    version = false
    "\t<string>#{verstr}</string>\n"
  else
    version = (s =~ /\bCFBundleVersion\b/)
    nil
  end
}
  
#  Modify InfoPlist.strings
Dir["xcode-build/*.lproj/InfoPlist.strings"].each { |nm|
  modify_file(nm) { |s|
    s = s.kconv(Kconv::UTF8, Kconv::UTF16)
    if s =~ /Copyright/ && s.sub!(/Toshi Nagata, [-0-9]+/, "Toshi Nagata, #{yrange}")
      s = s.kconv(Kconv::UTF16, Kconv::UTF8)
    else
      nil
    end
  }
}

#  Modify Molby.iss
modify_file("build-win/molby64.iss") { |s|
  if s =~ /AppVerName/ && s.sub!(/\(.*\)*/, "(#{verstr})")
    s
  else
    nil
  end
}

modify_file("build-win32/molby32.iss") { |s|
  if s =~ /AppVerName/ && s.sub!(/\(.*\)*/, "(#{verstr})")
    s
  else
    nil
  end
  }
  
#  Modify MyVersion.c
modify_file("wxSources/MyVersion.c") { |s|
  if s =~ /Version/ && s.sub!(/\".*\"/, "\"#{verstr}\"")
    s
  elsif s =~ /Copyright/ && s =~ /Toshi Nagata/ && s.sub!(/\d\d\d\d(-\d\d\d\d)?/, yrange)
    s
  else
    nil
  end
}

#  Modify doc_source.html
modify_file("Document_Sources/src/doc_source.html") { |s|
  if s =~ /Version/ && s =~ /<!-- version -->/ && s.sub!(/[Vv][-.0-9 A-Za-z_]*/, "Version #{ver}")
    s
  elsif s =~ /<!-- copyright -->/ && s.sub!(/\d\d\d\d(-\d\d\d\d)?/, yrange)
    s
  else 
    nil
  end
}

#  Modify README
modify_file("README") { |s|
  if s =~ /        Version/ && s.sub!(/[Vv][-.0-9 A-Za-z_]*/, "Version #{ver}")
    s
  elsif s =~ /       Copyright/ && s =~ /Toshi Nagata/ && s.sub!(/\d\d\d\d(-\d\d\d\d)?/, yrange)
    s
  else
    nil
  end
}
