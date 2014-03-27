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
    scrdir = $home_directory + "/Molby/gamess"
	begin
	  mkdir_recursive(scrdir)
	rescue
	  error_message_box("Cannot create directory #{scrdir}: " + $!.to_s)
	  return
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
	if $platform == "win"
	  scrbody = inpbody
	else
      scrbody = scrprefix
	end
    filecopy(inpname, scrprefix + ".F05")

    #  Prepare environmental variables
    auxdir = "#{gmsdir}#{sep}auxdata"
    ENV["ERICFMT"] = "#{auxdir}#{sep}ericfmt.dat"
    ENV["MCPPATH"] = "#{auxdir}#{sep}MCP"
    ENV["BASPATH"] = "#{auxdir}#{sep}BASES"
    ENV["QUANPOL"] = "#{auxdir}#{sep}QUANPOL"
    ENV["EXTBAS"] = "/dev/null"
    ENV["IRCDATA"] = "#{scrbody}.irc"
    ENV["PUNCH"] = "#{scrbody}.dat"
    ENV["INPUT"] = "#{scrbody}.F05"
    ENV["AOINTS"] = "#{scrbody}.F08"
    ENV["MOINTS"] = "#{scrbody}.F09"
    ENV["DICTNRY"] = "#{scrbody}.F10"
    ENV["DRTFILE"] = "#{scrbody}.F11"
    ENV["CIVECTR"] = "#{scrbody}.F12"
    ENV["CASINTS"] = "#{scrbody}.F13"
    ENV["CIINTS"] = "#{scrbody}.F14"
    ENV["WORK15"] = "#{scrbody}.F15"
    ENV["WORK16"] = "#{scrbody}.F16"
    ENV["CSFSAVE"] = "#{scrbody}.F17"
    ENV["FOCKDER"] = "#{scrbody}.F18"
    ENV["WORK19"] = "#{scrbody}.F19"
    ENV["DASORT"] = "#{scrbody}.F20"
    ENV["DFTINTS"] = "#{scrbody}.F21"
    ENV["DFTGRID"] = "#{scrbody}.F22"
    ENV["JKFILE"] = "#{scrbody}.F23"
    ENV["ORDINT"] = "#{scrbody}.F24"
    ENV["EFPIND"] = "#{scrbody}.F25"
    ENV["PCMDATA"] = "#{scrbody}.F26"
    ENV["PCMINTS"] = "#{scrbody}.F27"
    ENV["MLTPL"] = "#{scrbody}.F28"
    ENV["MLTPLT"] = "#{scrbody}.F29"
    ENV["DAFL30"] = "#{scrbody}.F30"
    ENV["SOINTX"] = "#{scrbody}.F31"
    ENV["SOINTY"] = "#{scrbody}.F32"
    ENV["SOINTZ"] = "#{scrbody}.F33"
    ENV["SORESC"] = "#{scrbody}.F34"
    ENV["SIMEN"] = "#{scrbody}.simen"
    ENV["SIMCOR"] = "#{scrbody}.simcor"
    ENV["GCILIST"] = "#{scrbody}.F37"
    ENV["HESSIAN"] = "#{scrbody}.F38"
    ENV["SOCCDAT"] = "#{scrbody}.F40"
    ENV["AABB41"] = "#{scrbody}.F41"
    ENV["BBAA42"] = "#{scrbody}.F42"
    ENV["BBBB43"] = "#{scrbody}.F43"
    ENV["MCQD50"] = "#{scrbody}.F50"
    ENV["MCQD51"] = "#{scrbody}.F51"
    ENV["MCQD52"] = "#{scrbody}.F52"
    ENV["MCQD53"] = "#{scrbody}.F53"
    ENV["MCQD54"] = "#{scrbody}.F54"
    ENV["MCQD55"] = "#{scrbody}.F55"
    ENV["MCQD56"] = "#{scrbody}.F56"
    ENV["MCQD57"] = "#{scrbody}.F57"
    ENV["MCQD58"] = "#{scrbody}.F58"
    ENV["MCQD59"] = "#{scrbody}.F59"
    ENV["MCQD60"] = "#{scrbody}.F60"
    ENV["MCQD61"] = "#{scrbody}.F61"
    ENV["MCQD62"] = "#{scrbody}.F62"
    ENV["MCQD63"] = "#{scrbody}.F63"
    ENV["MCQD64"] = "#{scrbody}.F64"
    ENV["NMRINT1"] = "#{scrbody}.F61"
    ENV["NMRINT2"] = "#{scrbody}.F62"
    ENV["NMRINT3"] = "#{scrbody}.F63"
    ENV["NMRINT4"] = "#{scrbody}.F64"
    ENV["NMRINT5"] = "#{scrbody}.F65"
    ENV["NMRINT6"] = "#{scrbody}.F66"
    ENV["DCPHFH2"] = "#{scrbody}.F67"
    ENV["DCPHF21"] = "#{scrbody}.F68"
    ENV["GVVPT"] = "#{scrbody}.F69"

    #    next files are used only during coupled cluster runs, so let's
    #    display the numerous definitions only if they are to be used.
    ENV["CCREST"] = "#{scrbody}.F70"
    ENV["CCDIIS"] = "#{scrbody}.F71"
    ENV["CCINTS"] = "#{scrbody}.F72"
    ENV["CCT1AMP"] = "#{scrbody}.F73"
    ENV["CCT2AMP"] = "#{scrbody}.F74"
    ENV["CCT3AMP"] = "#{scrbody}.F75"
    ENV["CCVM"] = "#{scrbody}.F76"
    ENV["CCVE"] = "#{scrbody}.F77"
    ENV["EOMSTAR"] = "#{scrbody}.F80"
    ENV["EOMVEC1"] = "#{scrbody}.F81"
    ENV["EOMVEC2"] = "#{scrbody}.F82"
    ENV["EOMHC1"] = "#{scrbody}.F83"
    ENV["EOMHC2"] = "#{scrbody}.F84"
    ENV["EOMHHHH"] = "#{scrbody}.F85"
    ENV["EOMPPPP"] = "#{scrbody}.F86"
    ENV["EOMRAMP"] = "#{scrbody}.F87"
    ENV["EOMRTMP"] = "#{scrbody}.F88"
    ENV["EOMDG12"] = "#{scrbody}.F89"
    ENV["MMPP"] = "#{scrbody}.F90"
    ENV["MMHPP"] = "#{scrbody}.F91"
    ENV["MMCIVEC"] = "#{scrbody}.F92"
    ENV["MMCIVC1"] = "#{scrbody}.F93"
    ENV["MMCIITR"] = "#{scrbody}.F94"
    ENV["MMNEXM"] = "#{scrbody}.F95"
    ENV["MMNEXE"] = "#{scrbody}.F96"
    ENV["MMNREXM"] = "#{scrbody}.F97"
    ENV["MMNREXE"] = "#{scrbody}.F98 "
    #
    #     next are for TDHFX code, not used by current GAMESS
    #
    ENV["OLI201"] = "#{scrbody}.F201"
    ENV["OLI202"] = "#{scrbody}.F202"
    ENV["OLI203"] = "#{scrbody}.F203"
    ENV["OLI204"] = "#{scrbody}.F204"
    ENV["OLI205"] = "#{scrbody}.F205"
    ENV["OLI206"] = "#{scrbody}.F206"
    ENV["OLI207"] = "#{scrbody}.F207"
    ENV["OLI208"] = "#{scrbody}.F208"
    ENV["OLI209"] = "#{scrbody}.F209"
    ENV["OLI210"] = "#{scrbody}.F210"
    ENV["OLI211"] = "#{scrbody}.F211"
    ENV["OLI212"] = "#{scrbody}.F212"
    ENV["OLI213"] = "#{scrbody}.F213"
    ENV["OLI214"] = "#{scrbody}.F214"
    ENV["OLI215"] = "#{scrbody}.F215"
    ENV["OLI216"] = "#{scrbody}.F216"
    ENV["OLI217"] = "#{scrbody}.F217"
    ENV["OLI218"] = "#{scrbody}.F218"
    ENV["OLI219"] = "#{scrbody}.F219"
    ENV["OLI220"] = "#{scrbody}.F220"
    ENV["OLI221"] = "#{scrbody}.F221"
    ENV["OLI222"] = "#{scrbody}.F222"
    ENV["OLI223"] = "#{scrbody}.F223"
    ENV["OLI224"] = "#{scrbody}.F224"
    ENV["OLI225"] = "#{scrbody}.F225"
    ENV["OLI226"] = "#{scrbody}.F226"
    ENV["OLI227"] = "#{scrbody}.F227"
    ENV["OLI228"] = "#{scrbody}.F228"
    ENV["OLI229"] = "#{scrbody}.F229"
    ENV["OLI230"] = "#{scrbody}.F230"
    ENV["OLI231"] = "#{scrbody}.F231"
    ENV["OLI232"] = "#{scrbody}.F232"
    ENV["OLI233"] = "#{scrbody}.F233"
    ENV["OLI234"] = "#{scrbody}.F234"
    ENV["OLI235"] = "#{scrbody}.F235"
    ENV["OLI236"] = "#{scrbody}.F236"
    ENV["OLI237"] = "#{scrbody}.F237"
    ENV["OLI238"] = "#{scrbody}.F238"
    ENV["OLI239"] = "#{scrbody}.F239"

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
      fp.print "-env ENVFIL #{envfil} -wdir #{scrdir} -n #{ncpus} #{gmsdir}#{sep}gamess.#{gmsvers}.exe\n"
      fp.close
    end
    
    fpout.close
    fplog = File.open(logname, "r")
    size = 0
    lines = []
    last_line = ""
    
    #  Callback procs
	term_callback = lambda { |m, n|
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
	
    timer_callback = lambda { |m, n|
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
		  :action=>lambda { |it| set_attr("use_internal", :enabled=>(it[:value] == 2)) } ),

		item(:checkbox, :title=>"Use internal coordinates for structure optimization", :tag=>"use_internal"),
		-1, -1, -1,

		item(:text, :title=>"Charge"),
		item(:textfield, :width=>80, :tag=>"charge"),
		item(:text, :title=>"Multiplicity"),
		item(:textfield, :width=>80, :tag=>"mult"),

		item(:checkbox, :title=>"Use DFT", :tag=>"dft",
		  :action=>lambda { |it| set_attr("dfttype", :enabled=>(it[:value] != 0)) } ),
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
		  :action=>lambda { |it|
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
		  :action=>lambda { |it|
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
  
  def cmd_create_cube
    grid = default_MO_grid
	if grid == nil
	  Dialog.run {
	    layout(1,
		  item(:text, :title=>"This molecule does not contain MO information."))
	  }
	  return
	end
    mos = selected_MO
	if mos == nil || mos.length == 0
      Dialog.run {
	    layout(1,
		  item(:text, :title=>"Please select MO(s) in the MO Info table."))
      }
	  return
	end
	hash = Dialog.run {
	  layout(1,
	    item(:text, :title=>"Please specify cube dimensions (in bohr unit):"),
	    layout(4,
		  item(:text, :title=>"Origin"),
		  item(:textfield, :width=>100, :height=>20, :tag=>"originx", :value=>sprintf("%.6f", grid[0].x)),
		  item(:textfield, :width=>100, :height=>20, :tag=>"originy", :value=>sprintf("%.6f", grid[0].y)),
		  item(:textfield, :width=>100, :height=>20, :tag=>"originz", :value=>sprintf("%.6f", grid[0].z)),
		  item(:text, :title=>"Delta"),
		  item(:textfield, :width=>100, :height=>20, :tag=>"deltax", :value=>sprintf("%.6f", grid[1])),
		  item(:textfield, :width=>100, :height=>20, :tag=>"deltay", :value=>sprintf("%.6f", grid[2])),
		  item(:textfield, :width=>100, :height=>20, :tag=>"deltaz", :value=>sprintf("%.6f", grid[3])),
		  item(:text, :title=>"Step"),
		  item(:textfield, :width=>100, :height=>20, :tag=>"stepx", :value=>grid[4].to_s),
		  item(:textfield, :width=>100, :height=>20, :tag=>"stepy", :value=>grid[5].to_s),
		  item(:textfield, :width=>100, :height=>20, :tag=>"stepz", :value=>grid[6].to_s)))
	}
	if hash[:status] == 0
	  path = self.path || self.name
	  dir = self.dir || Dir.pwd
	  origin = Vector3D[hash["originx"], hash["originy"], hash["originz"]]
	  dx = hash["deltax"]
	  dy = hash["deltay"]
	  dz = hash["deltaz"]
	  nx = hash["stepx"]
	  ny = hash["stepy"]
	  nz = hash["stepz"]
	  basename = File.basename(path, ".*")
	  filenames = []
	  mo_type = self.mo_type
	  mos.each { |n|
	    fname1 = fname2 = nil
	    alpha = (mo_type != "UHF" ? "" : "alpha ")
		a = (mo_type != "UHF" ? "" : "a")
	    fname1 = Dialog.save_panel("Cube file name for #{alpha}MO #{n}", dir, basename + "_#{n}#{a}.cube", "Gaussian cube file (*.cube)|*.cube")
		if (mo_type == "UHF")
		  fname2 = Dialog.save_panel("Cube file name for beta MO #{n}", dir, basename + "_#{n}b.cube", "Gaussian cube file (*.cube)|*.cube")
		end
		filenames.push([n, fname1, fname2])
	  }
	  filenames.each { |pair|
	    n = pair[0]
		alpha = (mo_type != "UHF" ? "" : "alpha ")
	    show_progress_panel("Creating cube file for #{alpha}MO #{n}...")
		if pair[1]
		  cubegen(pair[1], n, origin, dx, dy, dz, nx, ny, nz, true)
		end
		if pair[2] && mo_type == "UHF"
		  set_progress_message("Creating cube file for beta MO #{n}...")
		  cubegen(pair[2], n, origin, dx, dy, dz, nx, ny, nz, true, true)
		end
	    hide_progress_panel
 	  }
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

register_menu("QChem\tCreate GAMESS Input...",
  :cmd_create_gamess_input, :non_empty)
register_menu("QChem\tCreate MOPAC6 Input...",
  :cmd_create_mopac_input, :non_empty)    # mopac6.rb
register_menu("QChem\tCreate MO Cube...",
  :cmd_create_cube, lambda { |m| m && m.mo_type } )
