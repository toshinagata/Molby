#
#  gamess.rb
#
#  Created by Toshi Nagata on 2009/11/22.
#  Copyright 2009 Toshi Nagata. All rights reserved.
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

  def Molecule.read_gamess_basis_sets(fname)
    # $gamess_basis = Hash.new unless $gamess_basis
	$gamess_ecp = Hash.new unless $gamess_ecp
	# $gamess_basis_desc = Hash.new unless $gamess_basis_desc
	# $gamess_basis_keys = [] unless $gamess_basis_keys
	basename = File.basename(fname, ".*")
	keys = []
	descname = nil
    File.open(fname, "r") { |fp|
	  while (s = fp.gets)
	    if s =~ /^\s*!\s*/ && descname == nil
	      #  Get the descriptive name from the first comment line
		  s = Regexp.last_match.post_match
		  if s =~ /  EMSL/
		    descname = Regexp.last_match.pre_match
		  else
		    descname = (s.split)[0]
		  end
		  $gamess_basis_desc[basename] = descname
		  next
		end
		ss, bas = (s.split)[0..1]     #  Tokens delimited by whitespaces
		next if ss == nil || ss == ""
		if ss == "$ECP"
		#  Read the ECP section
		  ecp = $gamess_ecp[basename]
		  if !ecp
		    ecp = []
			$gamess_ecp[basename] = ecp
		  end
		  keys << basename unless keys.include?(basename)
		  while (s = fp.gets)
		    break if s=~ /\$END/
			ary = s.split  #  (PNAME, PTYPE, IZCORE, LMAX+1)
			elem = ary[0].match(/\w{1,2}/).to_a[0]  #  ary[0] should begin with an element name
			ap = Parameter.builtin.elements.find { |p| p.name.casecmp(elem) == 0 }
			raise MolbyError, "the ECP definition does not begin with an element name: #{s}" if !ap
			ecpdef = s
			ln = 1
			(0..Integer(ary[3])).each {
			  s = fp.gets
			  raise MolbyError, "the ECP definition ends unexpectedly at line #{ln} for: #{s}" if !s
			  ln += 1
			  ary = s.split  #  (NGPOT, comment)
			  ecpdef += s
			  (1..Integer(ary[0])).each {
			    s = fp.gets
				raise MolbyError, "the ECP definition ends unexpectedly at line #{ln} for: #{s}" if !s
				ln += 1
				ecpdef += s
			  }
			}
		    ecp[ap.index] = ecpdef
		  end
		elsif ss =~ /\W/
		  #  Comments or other unrecognizable lines
		  next
		elsif (ap = Parameter.builtin.elements.find { |p| p.name.casecmp(ss) == 0 || p.fullname.casecmp(ss) == 0 })
		  #  Valid basis definition
		  if bas == nil || bas =~ /\W/
		    bas = basename
		  end
		  basis = $gamess_basis[bas]
		  if !basis
		    basis = []
			$gamess_basis[bas] = basis
		  end
		  keys << bas unless keys.include?(bas)
		  basdef = ""
		  while (s = fp.gets) && s =~ /\S/
		    basdef += s
		  end
		  basis[ap.index] = basdef
		else
		  raise MolbyError, "ss is not a valid element symbol or name: #{s}"
		end
	  end
    }
	unless $gamess_basis_keys.include?(basename)
	  $gamess_basis_keys.push(basename)
	end
  end

  #  Execute GAMESS (inferior copy of rungms script)
  #  inpname is the input file
  #  mol (optional) is the molecule from which the GAMESS input was built.
  #  If mol is specified and RUNTYP=OPTIMIZE, then the intermediate structures are
  #  displayed real-time.
  def Molecule.execute_gamess(inpname, mol = nil)
    gmsname = get_global_settings("gamess.executable_path")
    gmsdir = nil
    gmsvers = nil
    while 1
      if gmsname == nil || !File.exist?(gmsname)
        gmsname = Dialog.open_panel("Please locate the GAMESS executable")
        exit if gmsname == nil
      end
      gmsbase = File.basename(gmsname)
      gmsdir = File.dirname(gmsname)
      if gmsbase =~ /gamess\.(.*)\.(exe|x)$/i
        gmsvers = $1
        break
      else
        gmsname = nil
        error_message_box(gmsbase + " does not look like a GAMESS executable!")
      end
    end

	if mol == nil
		mol = Molecule.open(inpname)
		if mol == nil
			error_message_box("Cannot open #{inpname} as GAMESS input")
			return
		end
	end
	
    inpbase = File.basename(inpname)
    inpdir = File.dirname(inpname)
    inpbody = inpbase.sub(/\.inp$/, "")
    logbase = inpbody + ".log"

    set_global_settings("gamess.executable_path", gmsname)

    ncpus = get_global_settings("gamess.ncpus").to_i
	if ncpus == 0
	  ncpus = 1
	end
	
    #  Prepare the scratch directory in the home directory
    #  (Not in the document home to avoid space-containing path in Windows)
    scrdir = document_home.sub(/\/My Documents/, "") + "/gamess"
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

    if $platform == "win"
      sep = "\\"
      scrdir.gsub!("/", sep)
      gmsdir.gsub!("/", sep)
    else
      sep = "/"
    end

    #  Get the host name etc.
    hostname = backquote("hostname").chomp
    if $platform == "win"
	  s = backquote("cmd.exe /c dir #{scrdir}")
      freebytes = s.split("\n").pop.match(/([0-9,]+)[^0-9]*$/).to_a[1]
      if freebytes
        freebytes = (freebytes.gsub(",","").to_i / 1024).to_s + " Kbytes"
      else
        freebytes = "(unknown)"
      end
      uname = backquote("cmd.exe /c ver").to_s.gsub("\n", "")
    else
      freebytes = `df -k #{scrdir}`
      uname = `uname`.chomp
    end

    #  Redirect standard output to the log file
    logname = scrdir + sep + logbase
    fpout = File.open(logname, "w")
    fpout.print "----- GAMESS execution script -----\n"
    fpout.print "This job is running on host #{hostname}\n"
    fpout.print "under operating system #{uname} at #{Time.now.to_s}\n"
    fpout.print "Available scratch disk space (Kbyte units) at beginning of the job is\n"
    fpout.print freebytes

    #  Copy the input file
    scrprefix = scrdir + sep + inpbody
    filecopy(inpname, scrprefix + ".F05")

    #  Prepare environmental variables
    auxdir = "#{gmsdir}#{sep}auxdata"
    ENV["ERICFMT"] = "#{auxdir}#{sep}ericfmt.dat"
    ENV["MCPPATH"] = "#{auxdir}#{sep}MCP"
    ENV["BASPATH"] = "#{auxdir}#{sep}BASES"
    ENV["QUANPOL"] = "#{auxdir}#{sep}QUANPOL"
    ENV["EXTBAS"] = "/dev/null"
    ENV["IRCDATA"] = "#{scrprefix}.irc"
    ENV["PUNCH"] = "#{scrprefix}.dat"
    ENV["INPUT"] = "#{scrprefix}.F05"
    ENV["AOINTS"] = "#{scrprefix}.F08"
    ENV["MOINTS"] = "#{scrprefix}.F09"
    ENV["DICTNRY"] = "#{scrprefix}.F10"
    ENV["DRTFILE"] = "#{scrprefix}.F11"
    ENV["CIVECTR"] = "#{scrprefix}.F12"
    ENV["CASINTS"] = "#{scrprefix}.F13"
    ENV["CIINTS"] = "#{scrprefix}.F14"
    ENV["WORK15"] = "#{scrprefix}.F15"
    ENV["WORK16"] = "#{scrprefix}.F16"
    ENV["CSFSAVE"] = "#{scrprefix}.F17"
    ENV["FOCKDER"] = "#{scrprefix}.F18"
    ENV["WORK19"] = "#{scrprefix}.F19"
    ENV["DASORT"] = "#{scrprefix}.F20"
    ENV["DFTINTS"] = "#{scrprefix}.F21"
    ENV["DFTGRID"] = "#{scrprefix}.F22"
    ENV["JKFILE"] = "#{scrprefix}.F23"
    ENV["ORDINT"] = "#{scrprefix}.F24"
    ENV["EFPIND"] = "#{scrprefix}.F25"
    ENV["PCMDATA"] = "#{scrprefix}.F26"
    ENV["PCMINTS"] = "#{scrprefix}.F27"
    ENV["MLTPL"] = "#{scrprefix}.F28"
    ENV["MLTPLT"] = "#{scrprefix}.F29"
    ENV["DAFL30"] = "#{scrprefix}.F30"
    ENV["SOINTX"] = "#{scrprefix}.F31"
    ENV["SOINTY"] = "#{scrprefix}.F32"
    ENV["SOINTZ"] = "#{scrprefix}.F33"
    ENV["SORESC"] = "#{scrprefix}.F34"
    ENV["SIMEN"] = "#{scrprefix}.simen"
    ENV["SIMCOR"] = "#{scrprefix}.simcor"
    ENV["GCILIST"] = "#{scrprefix}.F37"
    ENV["HESSIAN"] = "#{scrprefix}.F38"
    ENV["SOCCDAT"] = "#{scrprefix}.F40"
    ENV["AABB41"] = "#{scrprefix}.F41"
    ENV["BBAA42"] = "#{scrprefix}.F42"
    ENV["BBBB43"] = "#{scrprefix}.F43"
    ENV["MCQD50"] = "#{scrprefix}.F50"
    ENV["MCQD51"] = "#{scrprefix}.F51"
    ENV["MCQD52"] = "#{scrprefix}.F52"
    ENV["MCQD53"] = "#{scrprefix}.F53"
    ENV["MCQD54"] = "#{scrprefix}.F54"
    ENV["MCQD55"] = "#{scrprefix}.F55"
    ENV["MCQD56"] = "#{scrprefix}.F56"
    ENV["MCQD57"] = "#{scrprefix}.F57"
    ENV["MCQD58"] = "#{scrprefix}.F58"
    ENV["MCQD59"] = "#{scrprefix}.F59"
    ENV["MCQD60"] = "#{scrprefix}.F60"
    ENV["MCQD61"] = "#{scrprefix}.F61"
    ENV["MCQD62"] = "#{scrprefix}.F62"
    ENV["MCQD63"] = "#{scrprefix}.F63"
    ENV["MCQD64"] = "#{scrprefix}.F64"
    ENV["NMRINT1"] = "#{scrprefix}.F61"
    ENV["NMRINT2"] = "#{scrprefix}.F62"
    ENV["NMRINT3"] = "#{scrprefix}.F63"
    ENV["NMRINT4"] = "#{scrprefix}.F64"
    ENV["NMRINT5"] = "#{scrprefix}.F65"
    ENV["NMRINT6"] = "#{scrprefix}.F66"
    ENV["DCPHFH2"] = "#{scrprefix}.F67"
    ENV["DCPHF21"] = "#{scrprefix}.F68"
    ENV["GVVPT"] = "#{scrprefix}.F69"

    #    next files are used only during coupled cluster runs, so let's
    #    display the numerous definitions only if they are to be used.
    ENV["CCREST"] = "#{scrprefix}.F70"
    ENV["CCDIIS"] = "#{scrprefix}.F71"
    ENV["CCINTS"] = "#{scrprefix}.F72"
    ENV["CCT1AMP"] = "#{scrprefix}.F73"
    ENV["CCT2AMP"] = "#{scrprefix}.F74"
    ENV["CCT3AMP"] = "#{scrprefix}.F75"
    ENV["CCVM"] = "#{scrprefix}.F76"
    ENV["CCVE"] = "#{scrprefix}.F77"
    ENV["EOMSTAR"] = "#{scrprefix}.F80"
    ENV["EOMVEC1"] = "#{scrprefix}.F81"
    ENV["EOMVEC2"] = "#{scrprefix}.F82"
    ENV["EOMHC1"] = "#{scrprefix}.F83"
    ENV["EOMHC2"] = "#{scrprefix}.F84"
    ENV["EOMHHHH"] = "#{scrprefix}.F85"
    ENV["EOMPPPP"] = "#{scrprefix}.F86"
    ENV["EOMRAMP"] = "#{scrprefix}.F87"
    ENV["EOMRTMP"] = "#{scrprefix}.F88"
    ENV["EOMDG12"] = "#{scrprefix}.F89"
    ENV["MMPP"] = "#{scrprefix}.F90"
    ENV["MMHPP"] = "#{scrprefix}.F91"
    ENV["MMCIVEC"] = "#{scrprefix}.F92"
    ENV["MMCIVC1"] = "#{scrprefix}.F93"
    ENV["MMCIITR"] = "#{scrprefix}.F94"
    ENV["MMNEXM"] = "#{scrprefix}.F95"
    ENV["MMNEXE"] = "#{scrprefix}.F96"
    ENV["MMNREXM"] = "#{scrprefix}.F97"
    ENV["MMNREXE"] = "#{scrprefix}.F98 "
    #
    #     next are for TDHFX code, not used by current GAMESS
    #
    ENV["OLI201"] = "#{scrprefix}.F201"
    ENV["OLI202"] = "#{scrprefix}.F202"
    ENV["OLI203"] = "#{scrprefix}.F203"
    ENV["OLI204"] = "#{scrprefix}.F204"
    ENV["OLI205"] = "#{scrprefix}.F205"
    ENV["OLI206"] = "#{scrprefix}.F206"
    ENV["OLI207"] = "#{scrprefix}.F207"
    ENV["OLI208"] = "#{scrprefix}.F208"
    ENV["OLI209"] = "#{scrprefix}.F209"
    ENV["OLI210"] = "#{scrprefix}.F210"
    ENV["OLI211"] = "#{scrprefix}.F211"
    ENV["OLI212"] = "#{scrprefix}.F212"
    ENV["OLI213"] = "#{scrprefix}.F213"
    ENV["OLI214"] = "#{scrprefix}.F214"
    ENV["OLI215"] = "#{scrprefix}.F215"
    ENV["OLI216"] = "#{scrprefix}.F216"
    ENV["OLI217"] = "#{scrprefix}.F217"
    ENV["OLI218"] = "#{scrprefix}.F218"
    ENV["OLI219"] = "#{scrprefix}.F219"
    ENV["OLI220"] = "#{scrprefix}.F220"
    ENV["OLI221"] = "#{scrprefix}.F221"
    ENV["OLI222"] = "#{scrprefix}.F222"
    ENV["OLI223"] = "#{scrprefix}.F223"
    ENV["OLI224"] = "#{scrprefix}.F224"
    ENV["OLI225"] = "#{scrprefix}.F225"
    ENV["OLI226"] = "#{scrprefix}.F226"
    ENV["OLI227"] = "#{scrprefix}.F227"
    ENV["OLI228"] = "#{scrprefix}.F228"
    ENV["OLI229"] = "#{scrprefix}.F229"
    ENV["OLI230"] = "#{scrprefix}.F230"
    ENV["OLI231"] = "#{scrprefix}.F231"
    ENV["OLI232"] = "#{scrprefix}.F232"
    ENV["OLI233"] = "#{scrprefix}.F233"
    ENV["OLI234"] = "#{scrprefix}.F234"
    ENV["OLI235"] = "#{scrprefix}.F235"
    ENV["OLI236"] = "#{scrprefix}.F236"
    ENV["OLI237"] = "#{scrprefix}.F237"
    ENV["OLI238"] = "#{scrprefix}.F238"
    ENV["OLI239"] = "#{scrprefix}.F239"

    if $platform == "win"
	  if ncpus < 2
	    ncpus = 2
	  end
      fpout.print "Microsoft MPI will be running GAMESS on 1 node.\n"
      fpout.print "The binary kicked off by 'mpiexec' is gamess.#{gmsvers}.exe\n"
      fpout.print "MS-MPI will run #{ncpus} compute process\n"
      #  File containing environmental variables
      envfil = "#{scrprefix}.GMS.ENV"
      fp = File.open(envfil, "w")
      ENV.each { |k, v| fp.print "#{k}=#{v}\n" }
      fp.close
      #  File containing arguments to mpiexec
      procfil = "#{scrprefix}.processes.mpd"
      fp = File.open(procfil, "w")
      fp.print "-env ENVFIL #{envfil} -n #{ncpus} #{gmsdir}#{sep}gamess.#{gmsvers}.exe\n"
      fp.close
    end
    
    fpout.close
    fplog = File.open(logname, "r")
    size = 0
    lines = []
    last_line = ""
    
    #  Callback procs
	term_callback = proc { |m, n|
	  msg = "GAMESS execution on #{inpbase} "
	  hmsg = "GAMESS "
	  if n == 0
	    msg += "succeeded."
		hmsg += "Completed"
		icon = :info
	  else
	    msg += "failed with status #{n}."
		hmsg += "Failed"
		icon = :error
	  end
	  msg += "\n(In directory #{inpdir})"
  	  ext_to_keep = [".dat", ".rst", ".trj", ".efp", ".gamma", ".log"]
	  ext_to_keep.each { |ex|
	    if File.exists?("#{scrprefix}#{ex}")
		  filecopy("#{scrprefix}#{ex}", "#{inpdir}#{sep}#{inpbody}#{ex}")
	    end
	  }
	  Dir.foreach(scrdir) { |file|
	    if file != "." && file != ".." && !ext_to_keep.include?(File.extname(file))
		  File.delete("#{scrdir}#{sep}#{file}")
	    end
	  }
	  erase_old_logs(scrdir, "latest", 5)
	  message_box(msg, hmsg, :ok, icon)
    }
	
    timer_callback = proc { |m, n|
      fplog.seek(0, IO::SEEK_END)
      sizec = fplog.tell
      if sizec > size
        #  Read new lines
        fplog.seek(size, IO::SEEK_SET)
        fplog.each_line { |line|
          if line[-1, 1] == "\n"
            lines.push(last_line + line)
            last_line = ""
          else
            last_line += line
            break
          end
        }
        size = fplog.tell
        last_i = nil
        i = 0
        nserch = -1
        while i < lines.count
          line = lines[i]
          if line =~ /GEOMETRY SEARCH POINT NSERCH= *(\d+)/
            nserch = $1.to_i
            last_i = i
          elsif line =~ /NSERCH:/
          #  print line
			if mol
			  dummy, n, grad = line.match(/NSERCH:[^0-9]*([0-9]+).*GRAD[^0-9]*([-.0-9]+)/).to_a
			  mol.show_text("Search: #{n}\nGradient: #{grad}")
			end
			last_i = i
          elsif nserch > 0 && line =~ /ddikick\.x/
            last_i = -1
            break
          elsif mol && nserch > 0 && line =~ /COORDINATES OF ALL ATOMS/
            #  There should be (natoms + 2) lines
            if i + mol.natoms + 3 <= lines.count
              coords = []
              (i + 3...i + 3 + mol.natoms).each { |j|
                name, charge, x, y, z = lines[j].split
                coords.push(Vector3D[x.to_f, y.to_f, z.to_f])
              }
              mol.create_frame([coords])
              mol.display
              last_i = i + mol.natoms + 2
              i = last_i   #  Skip the processed lines
            end
          end
          i += 1
        end
        if last_i == -1
          lines.clear
          break
        elsif last_i
          lines[0..last_i] = nil
        end
      end
      true
    }

    if $platform == "win"
	  cmdline = "cmd.exe /c \"mpiexec -configfile #{procfil} >>#{logname}\""
	else
	  hosts = "localhost " * ncpus
	  cmdline = "/bin/sh -c '#{gmsdir}/ddikick.x #{gmsdir}/gamess.#{gmsvers}.x #{inpbody} -ddi #{ncpus} #{ncpus} #{hosts} -scr #{scrdir} < /dev/null >>#{logname}'"
	end
	
	if mol
	  pid = mol.call_subprocess_async(cmdline, term_callback, timer_callback)
	  if pid < 0
	    error_message_box("GAMESS failed to start. Please examine GAMESS installation.")
	    return
	  end
	else
	  pid = call_subprocess(cmdline, "Running GAMESS")
	  term_callback(nil, pid)
	end

  end
  
  def export_gamess(fname, hash)

    def reorder_array(ary, ordered_sub_ary)
	  return ordered_sub_ary + (ary - ordered_sub_ary)
	end
	
	now = Time.now.to_s
	basename = File.basename(fname, ".*")
	
	#  Various settings
	icharg = hash["charge"]
	mult = hash["mult"]
	runtyp = ["ENERGY", "PROP", "OPTIMIZE"][Integer(hash["runtype"])]
	scftyp = ["RHF", "ROHF", "UHF"][Integer(hash["scftype"])]
	bssname = hash["basis"]
	bssname2 = hash["secondary_basis"]
	if hash["use_secondary_basis"] != 0 && bssname2 != bssname
	  use_2nd = true
	  element2 = hash["secondary_elements"].split(/[\s,]+/).map { |name| name.capitalize }
	else
	  use_2nd = false
	  bssname2 = bssname
	end
	basis = $gamess_basis[bssname]
	basis2 = $gamess_basis[bssname2]
	if !basis || !basis2
	  raise MolbyError, "Unknown basis set name??? \"#{bssname}\" or \"#{bssname2}\""
	end

	#  Use effective core potentials?
	ecp = $gamess_ecp[bssname]
	ecp2 = $gamess_ecp[bssname2]
	ecp_read = (ecp || ecp2 ? "ECP=READ " : "")

	#  Use only one built-in basis set?
	gbasis = nil
	if !use_2nd
	  case bssname
	  when "PM3"
		gbasis = "PM3 NGAUSS=3"
	  when "STO3G"
		gbasis = "STO NGAUSS=3"
	  when "321G"
		gbasis = "N21 NGAUSS=3"
	  when "631G"
		gbasis = "N31 NGAUSS=6"
	  when "631Gd"
		gbasis = "N31 NGAUSS=6 NDFUNC=1"
	  when "631Gdp"
		gbasis = "N31 NGAUSS=6 NDFUNC=1 NPFUNC=1"
	  when "6311G"
		gbasis = "N311 NGAUSS=6"
	  when "6311Gdp"
		gbasis = "N311 NGAUSS=6 NDFUNC=1 NPFUNC=1"
	  end
	end
    
	#  Count non-dummy atoms
	natoms = 0
	each_atom { |ap|
	  natoms += 1 if ap.atomic_number != 0
	}
	
	#  Fill hash with default values
	h = (hash["CONTRL"] ||= Hash.new)
	h["COORD"] ||= "UNIQUE"
	h["EXETYP"] ||= "RUN"
	h["ICHARG"] ||= icharg.to_s
	h["ICUT"] ||= "20"
	h["INTTYP"] ||= "HONDO"
	h["ITOL"] ||= "30"
	h["MAXIT"] ||= "200"
	h["MOLPLT"] ||= ".T."
	h["MPLEVL"] ||= "0"
	h["MULT"] ||= mult.to_s
	h["QMTTOL"] ||= "1e-08"
	h["RUNTYP"] ||= runtyp
	if hash["dft"] != 0
	  h["DFTTYP"] ||= hash["dfttype"]
	end
	if hash["use_internal"] != 0 && (hash["runtype"] == 2 || h["RUNTYP"] == "OPTIMIZE")
	  nzvar = natoms * 3 - 6  #  TODO: 3N-5 for linear molecules
	  h["NZVAR"] = nzvar.to_s
	else
	  nzvar = 0
	end
	h["SCFTYP"] ||= scftyp
	if ecp_read != ""
	  h["ECP"] ||= "READ"
	end
	h["UNITS"] = "ANGS"
	
	h = (hash["SCF"] ||= Hash.new)
	h["CONV"] ||= "1.0E-06"
	h["DIRSCF"] ||= ".T."
	h["FDIFF"] ||= ".T."
	h["DAMP"] ||= ".T."
	
	h = (hash["STATPT"] ||= Hash.new)
	h["NSTEP"] ||= "400"
	h["OPTTOL"] ||= "1.0E-06"

	h = (hash["SYSTEM"] ||= Hash.new)
	h["MEMDDI"] ||= "0"
	h["MWORDS"] ||= "16"
	h["TIMLIM"] ||= "50000"
	
	h = (hash["GUESS"] ||= Hash.new)
	h["GUESS"] ||= "HUCKEL"
	
	if gbasis
	  h = (hash["BASIS"] ||= Hash.new)
	  h["GBASIS"] ||= gbasis
	end
	
	if nzvar > 0
	  h = (hash["ZMAT"] ||= Hash.new)
	  h["DLC"] ||= ".T."
	  h["AUTO"] ||= ".T."
	end
	
	if hash["esp"] != 0
	  h = (hash["ELPOT"] ||= Hash.new)
	  h["IEPOT"] ||= "1"
	  h["OUTPUT"] ||= "PUNCH"
	  h["WHERE"] ||= "PDC"
	  h = (hash["PDC"] ||= Hash.new)
	  h["CONSTR"] ||= "NONE"
	  h["PTSEL"] ||= "CONNOLLY"
	end
	
    File.open(fname, "wb") { |fp|
	  fp.print "!  GAMESS input\n"
	  fp.print "!  Generated by Molby at #{now}\n"
	  fp.print "!  Basis set: " + $gamess_basis_desc[bssname] + "\n"
	  if use_2nd
	    fp.print "!  [" + element2.join(", ") + "]: " + $gamess_basis_desc[bssname2] + "\n"
	  end
	  controls = reorder_array(hash.keys.select { |k| hash[k].is_a?(Hash) },
		["CONTRL", "SCF", "STATPT", "SYSTEM", "GUESS", "BASIS", "ZMAT", "ELPOT", "PDC"])
	  controls.each { |k|
	    h = hash[k]
		next if h == nil || h.size == 0
		if k == "CONTRL"
		  ordered = ["COORD", "EXETYP", "ICHARG", "ICUT", "INTTYP", "ITOL", "MAXIT", "MOLPLT", "MPLEVL",
		             "MULT", "QMTTOL", "RUNTYP", "SCFTYP", "ECP", "UNITS", "DFTTYP", "NZVAR"]
		elsif k == "SCF"
		  ordered = ["CONV", "DIRSCF", "FDIFF", "DAMP"]
		elsif k == "STATPT"
		  ordered = ["NSTEP", "OPTTOL"]
		elsif k == "SYSTEM"
		  ordered = ["MEMDDI", "MWORDS", "TIMLIM"]
		elsif k == "GUESS"
		  ordered = ["GUESS"]
		elsif k == "BASIS"
		  ordered = ["GBASIS"]
		elsif k == "ZMAT"
		  ordered = ["DLC", "AUTO"]
		elsif k == "ELPOT"
		  ordered = ["IEPOT", "OUTPUT", "WHERE"]
		elsif k == "PDC"
		  ordered = ["CONSTR", "PTSEL"]
		else
		  ordered = []
		end
	    keys = reorder_array(h.keys, ordered)
		n = 0
		keys.each_with_index { |kk, i|
		  v = h[kk]
		  next if v == nil || v == ""
		  if n == 0
		    fp.printf " $%-6s", k
		  elsif n % 3 == 0
		    fp.print "\n        "
		  end
		  fp.printf " %s=%s", kk, h[kk].to_s
		  n += 1
		}
		if n > 0
		  fp.print " $END\n"
		end
	  }
	  fp.print " $DATA\n#{basename}\nC1 0\n"
	  secondary = []
	  each_atom { |ap|
	    next if ap.atomic_number == 0
		fp.printf "%-6s %4d %10.6f %10.6f %10.6f\n", ap.name, ap.atomic_number, ap.r.x, ap.r.y, ap.r.z
		if use_2nd && element2.include?(ap.element)
		  secondary[ap.index] = true
		end
	    if !gbasis
		  #  Basis specification followed by a blank line
		  bas = (secondary[ap.index] ? basis2 : basis)
		  if bas.is_a?(Array)
		    bas = bas[ap.atomic_number]
		  end
		  if !bas
		    raise MolbyError, "Basis set is not defined for atom #{ap.index}, element #{ap.element}"
		  end
		  fp.print bas
		  fp.print "\n"
		end
	  }
	  fp.print " $END\n"
	  ecp_ary = []
	  if ecp || ecp2
	    fp.print " $ECP\n"
		each_atom { |ap|
		  an = ap.atomic_number
		  next if an == 0
		  ecpp = (secondary[ap.index] ? ecp2 : ecp)
		  e = ecp_ary[an] || (ecpp && ecpp[an])
		  if e
		    #  Cache the PNAME of the $ECP entry and re-use it
		    ecp_ary[an] ||= (e.split)[0] + "\n"
		  else
		    e = ap.element.upcase + "-ECP NONE\n"
		  end
		  fp.print e
		}
	    fp.print " $END\n"
	  end
	}
	fname
  end
  
  def cmd_create_gamess_input
    if natoms == 0
      raise MolbyError, "cannot create GAMESS input; the molecule is empty"
    end

	#  Descriptive text and internal string for popup menus
#    bset_desc = ["PM3", "STO-3G", "3-21G", "6-31G", "6-31G(d)", "6-31G(d,p)", "6-311G", "6-311G(d,p)", "LanL2DZ"]
#	bset_internal = ["PM3", "STO3G", "321G", "631G", "631Gd", "631Gdp", "6311G", "6311Gdp", "LanL2DZ"]
    bset_desc = $gamess_basis_keys.map { |key| $gamess_basis_desc[key] }
	dft_desc = ["B3LYP"]
	dft_internal = ["B3LYP"]

	defaults = {"scftype"=>0, "runtype"=>0, "charge"=>"0", "mult"=>"1",
	  "basis"=>4, "use_secondary_basis"=>0, "secondary_elements"=>"",
	  "secondary_basis"=>8, "esp"=>0, "ncpus"=>"1"}

    hash = Dialog.run("GAMESS Export") {
      def load_basis_set_sub(item)
	    fname = Dialog.open_panel("Select a file containing GAMESS basis set:")
		if fname
		  Molecule.read_gamess_basis_sets(fname)
		  bset_desc_new = $gamess_basis_keys.map { |key| $gamess_basis_desc[key] }
		  sel1 = attr("basis", :value)
		  sel2 = attr("secondary_basis", :value)
		  set_attr("basis", :subitems=>bset_desc_new)
		  set_attr("basis", :value=>sel1)
		  set_attr("secondary_basis", :subitems=>bset_desc_new)
		  set_attr("secondary_basis", :value=>sel2)
		end
      end
	  def select_gamess_path(item)
	    while 1
	      fname = Dialog.open_panel("Locate GAMESS executable:")
		  return if fname == nil
		  bname = File.basename(fname)
		  if bname =~ /gamess\.(.*)\.(exe|x)$/i
		    set_attr("executable_path", :value=>fname)
		    return
		  else
		    error_message_box("\"#{bname}\" does not look like a GAMESS executable!  Please try again.")
		  end
		end
	  end
	  layout(4,
		item(:text, :title=>"SCF type"),
		item(:popup, :subitems=>["RHF", "ROHF", "UHF"], :tag=>"scftype"),
	    item(:text, :title=>"Run type"),
		item(:popup, :subitems=>["Energy", "Property", "Optimize"], :tag=>"runtype",
		  :action=>proc { |it| set_attr("use_internal", :enabled=>(it[:value] == 2)) } ),

		item(:checkbox, :title=>"Use internal coordinates for structure optimization", :tag=>"use_internal"),
		-1, -1, -1,

		item(:text, :title=>"Charge"),
		item(:textfield, :width=>80, :tag=>"charge"),
		item(:text, :title=>"Multiplicity"),
		item(:textfield, :width=>80, :tag=>"mult"),

		item(:checkbox, :title=>"Use DFT", :tag=>"dft",
		  :action=>proc { |it| set_attr("dfttype", :enabled=>(it[:value] != 0)) } ),
		-1,
		item(:text, :title=>"DFT type"),
		item(:popup, :subitems=>dft_desc, :tag=>"dfttype"),
		
		item(:line),
		-1, -1, -1,

		item(:text, :title=>"Basis set"),
		item(:popup, :subitems=>bset_desc, :tag=>"basis"),
		-1,
		-1,

		item(:button, :title=>"Load Basis Set...", :action=>:load_basis_set_sub),
		-1, -1, -1,
		
		item(:checkbox, :title=>"Use secondary basis set", :tag=>"use_secondary_basis",
		  :action=>proc { |it|
		    flag = (it[:value] != 0)
			set_attr("secondary_elements", :enabled=>flag)
			set_attr("secondary_basis", :enabled=>flag)
		  }),
		-1, -1, -1,

		item(:text, :title=>"   Elements"),
		item(:textfield, :width=>80, :tag=>"secondary_elements"),
		item(:text, :title=>"Basis set"),
		item(:popup, :subitems=>bset_desc, :tag=>"secondary_basis"),
		
		item(:line),
		-1, -1, -1,
		
		item(:checkbox, :title=>"Calculate electrostatic potential (ESP)", :tag=>"esp"),
		-1, -1, -1,
	
		item(:line),
		-1, -1, -1,
		
		item(:checkbox, :title=>"Execute GAMESS on this machine", :tag=>"execute_local",
		  :action=>proc { |it|
		    flag = (it[:value] != 0)
			set_attr("executable_path", :enabled=>flag)
			set_attr("select_path", :enabled=>flag)
			set_attr("ncpus", :enabled=>flag)
		  }),
		-1, -1, -1,

		item(:text, :title=>"   Path"),
		item(:textfield, :width=>300, :tag=>"executable_path"),
		-1, -1,
		
		-1,
		item(:button, :title=>"Select Path...", :tag=>"select_path", :action=>:select_gamess_path),
		-1, -1,
		
		item(:text, :title=>"   N of CPUs"),
		item(:textfield, :width=>80, :tag=>"ncpus"),
		-1, -1,
		
		item(:line),
		-1, -1, -1
	  )
	  values = Hash.new
	  each_item { |it|
	    tag = it[:tag]
		if tag
		  values[tag] = (get_global_settings("gamess.#{tag}") || defaults[tag])
		  it[:value] = values[tag]
		end
	  }
	  set_attr("secondary_elements", :enabled=>(values["use_secondary_basis"] == 1))
	  set_attr("secondary_basis", :enabled=>(values["use_secondary_basis"] == 1))
	  set_attr("dfttype", :enabled=>(values["dft"] == 1))
	  set_attr("use_internal", :enabled=>(values["runtype"] == 2))
	  set_attr("executable_path", :enabled=>(values["execute_local"] == 1))
	  set_attr("select_path", :enabled=>(values["execute_local"] == 1))
	  set_attr("ncpus", :enabled=>(values["execute_local"] == 1))
	}
	hash.each_pair { |key, value|
	  next if key == :status
	  set_global_settings("gamess.#{key}", value)
	}
	if hash[:status] == 0
	  #  Specify basis by internal keys
	  hash["basis"] = $gamess_basis_keys[hash["basis"]]
	  hash["secondary_basis"] = $gamess_basis_keys[hash["secondary_basis"]]
	  hash["dfttype"] = dft_internal[hash["dfttype"]]
	  basename = (self.path ? File.basename(self.path, ".*") : self.name)
      fname = Dialog.save_panel("Export GAMESS input file:", self.dir, basename + ".inp", "GAMESS input file (*.inp)|*.inp|All files|*.*")
	  return nil if !fname
	  export_gamess(fname, hash)
	  if hash["execute_local"] == 1
	    Molecule.execute_gamess(fname, self)
	  end
	else
	  nil
	end
  end
  
end

$gamess_basis = {
  "PM3"   => " PM3 0\n",
  "STO3G" => " STO 3\n",
  "321G"  => " N21 3\n",
  "631G"  => " N31 6\n" }
$gamess_basis_desc = {
  "PM3"   => "PM3",
  "STO3G" => "STO-3G",
  "321G"  => "3-21G",
  "631G"  => "6-31G" }
$gamess_basis_keys = ["PM3", "STO3G", "321G", "631G"]

["631Gd", "631Gdp", "631+Gd", "631++Gdp", "6311Gdp", "6311+Gd", "6311++Gdp", "6311++G2d2p", "6311++G3df3pd", "LanL2DZ"].each { |n|
  Molecule.read_gamess_basis_sets("basis_sets/#{n}.txt")
}
