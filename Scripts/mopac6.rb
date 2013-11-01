#
#  mopac6.rb
#
#  Created by Toshi Nagata on 2013/11/01.
#  Copyright 2013 Toshi Nagata. All rights reserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation version 2 of the License.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

class Molecule

  #  Execute mopac6 (inferior copy of rungms script)
  #  inpname is the input file
  #  mol (optional) is the molecule from which the MOPAC input was built.
  def Molecule.execute_mopac(inpname, mol = nil)

    #  MOPAC executable
    mopexe = "#{ResourcePath}/mopac/mopac606"

    inpbase = File.basename(inpname)
    inpdir = File.dirname(inpname)
    inpbody = inpbase.sub(/\.\w*$/, "")

    #  Prepare the scratch directory in the home directory
    scrdir = document_home.sub(/\/My Documents/, "") + "/mopac"
    n = 0
    while File.exist?(scrdir) && !File.directory?(scrdir)
      if n == 0
        scrdir += ".1"
      else
        scrdir = scrdir.sub(".#{n}", ".#{n + 1}")
      end
      n += 1
    end
    if !File.exist?(scrdir)
      Dir.mkdir(scrdir)
    end
    scrdir = scrdir + "/" + inpbody + "." + $$.to_s + ".0"
    n = 0
    while File.exist?(scrdir)
      scrdir = scrdir.sub(".#{n}", ".#{n + 1}")
      n += 1
    end
    Dir.mkdir(scrdir)

	#  Copy the input file to the scratch directory
	cwd = Dir.pwd
	filecopy(inpname, scrdir + "/FOR005")
	term_callback = proc { |m, n|
	    if n == 0
		    filecopy("#{scrdir}/FOR006", "#{cwd}/#{inpdir}/#{inpbody}.out");
	    end
		Dir.foreach(scrdir) { |fn|
			next if fn == "." || fn == ".."
		    File.delete("#{scrdir}/#{fn}")
		}
		Dir.chdir(cwd)
		Dir.rmdir(scrdir)
	}
	timer_callback = nil
	
	#  Execute mopac
	Dir.chdir(scrdir)
	if mol
	  pid = mol.call_subprocess_async(mopexe, term_callback, timer_callback, "NUL", "NUL")
	else
	  pid = call_subprocess(mopexe, "Running MOPAC", nil, "NUL", "NUL")
	  term_callback.call(nil, pid)
	end
	Dir.chdir(cwd)
	if pid < 0
	  error_message_box("MOPAC failed to run.")
	  return
	end

  end
    
end

Molecule.execute_mopac("test01.inp", nil)
