#!/usr/bin/ruby

require 'kconv'

#  Get the version string
#  version = "X.X.X"
#  date = "yyyymmdd"
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
File.open("Version", "w") { |fp|
  fp.print "version = \"#{version}\"\n"
  fp.print "date = \"#{d}\"\n"
}
build = "build " + d
# verstr = "v#{ver} #{build}"
verstr = "v#{ver}"
yrange = (year > 2008 ? "2008-#{year}" : "2008")

def modify_file(name, &block)
  ary = IO.readlines(name)
  modified = false
  ary.each_with_index { |s, i|
    s = block.call(s)
    if s
      ary[i] = s
      modified = true
    end
  }
  if modified
    File.rename(name, name + "~")
    open(name, "wb") { |fp|
      ary.each { |s| fp.write(s) }
    }
  end
end

#  Modify Info.plist
nm = "xcode-build/Info.plist"
version = false
modify_file(nm) { |s|
  if version
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
modify_file("msw-build/molby.iss") { |s|
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
  else
    nil
  end
}

#  Modify doc_source.html
modify_file("Documents/src/doc_source.html") { |s|
  if s =~ /Version/ && s.sub!(/[Vv][-.0-9 A-Za-z_]*/, "Version #{ver}")
    s
  else
    nil
  end
}

#  Modify README
modify_file("README") { |s|
  if s =~ /        Version/ && s.sub!(/[Vv][-.0-9 A-Za-z_]*/, "Version #{ver}")
    s
  else
    nil
  end
}
