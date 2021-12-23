# coding: utf-8
#
#  loadsave.rb
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

require "scanf.rb"

class Molecule

  def loadcrd(filename)
    if natoms == 0
      raise MolbyError, "cannot load crd; the molecule is empty"
    end
	fp = open(filename, "rb")
    count = 0
	frame = 0
	coords = (0...natoms).collect { Vector3D[0, 0, 0] }
    periodic = (self.box && self.box[4].any? { |n| n != 0 })
	show_progress_panel("Loading AMBER crd file...")
#    puts "sframe = #{sframe}, pos = #{fp.pos}"
    while line = fp.gets
      line.chomp!
	  values = line.scan(/......../)
	  if fp.lineno == 1
	    #  The first line should be skipped. However, if this line seems to contain
		#  coordinates, then try reading them
		if values.size != 10 && (values.size > 10 || values.size != natoms * 3)
		  next  #  Wrong number of coordinates
		end
		if values.each { |v| Float(v) rescue break } == nil
		  #  The line contains non-number 
		  next
	    end
		puts "Loadcrd: The title line seems to be missing"
	  end
      if count + values.size > natoms * 3
        raise MolbyError, sprintf("crd format error - too many values at line %d in file %s; number of atoms = %d, current frame = %d", fp.lineno, fp.path, natoms, frame)
      end
	  values.each { |v|
	    coords[count / 3][count % 3] = Float(v)
	    count += 1
	  }
      if count == natoms * 3
	    #  End of frame
		if frame == 0 && atoms.all? { |ap| ap.r.x == 0.0 && ap.r.y == 0.0 && ap.r.z == 0.0 }
			#  Do not create a new frame
			atoms.each_with_index { |ap, i| ap.r = coords[i] }
		else
			create_frame([coords])
		end
		if periodic
		  #  Should have box information
		  line = fp.gets
		  if line == nil || (values = line.chomp.scan(/......../)).length != 3
		    #  Periodic but no cell information
			puts "Loadcrd: the molecule has a periodic cell but the crd file does not contain cell info"
			periodic = false
			redo
	      end
		  self.cell = [Float(values[0]), Float(values[1]), Float(values[2]), 90, 90, 90]
		end
        count = 0
        frame += 1
		if frame % 5 == 0
		  set_progress_message("Loading AMBER crd file...\n(#{frame} frames completed)")
		end
      end
    end
	fp.close
	if (count > 0)
	  raise MolbyError, sprintf("crd format error - file ended in the middle of frame %d", frame)
	end
	hide_progress_panel
	if frame > 0
	  self.frame = self.nframes - 1
	  return true
	else
	  return false
	end
  end
    
  def savecrd(filename)
    if natoms == 0
      raise MolbyError, "cannot save crd; the molecule is empty"
    end
    fp = open(filename, "wb")
	show_progress_panel("Saving AMBER crd file...")
	fp.printf("TITLE: %d atoms\n", natoms)
	cframe = self.frame
	nframes = self.nframes
	j = 0
	self.update_enabled = false
	begin
		(0...nframes).each { |i|
		  select_frame(i)
		  j = 0
		  while (j < natoms * 3)
			w = atoms[j / 3].r[j % 3]
			fp.printf(" %7.3f", w)
			fp.print("\n") if (j % 10 == 9 || j == natoms * 3 - 1)
			j += 1
		  end
		  if i % 5 == 0
			set_progress_message("Saving AMBER crd file...\n(#{i} frames completed)")
		  end
		}
	ensure
		self.update_enabled = true
	end
	select_frame(cframe)
	fp.close
	hide_progress_panel
	true
  end

  alias :loadmdcrd :loadcrd
  alias :savemdcrd :savecrd

  def sub_load_gamess_log_basis_set(lines, lineno)
    ln = 0
	while (line = lines[ln])
		ln += 1
		break if line =~ /SHELL\s+TYPE\s+PRIMITIVE/
	end
	ln += 1
	i = -1
	nprims = 0
	sym = -10  #  undefined
	ncomps = 0
	clear_basis_set
	while (line = lines[ln])
		ln += 1
		break if line =~ /TOTAL NUMBER OF BASIS SET/
		if line =~ /^\s*$/
		  #  End of one shell
		  add_gaussian_orbital_shell(i, sym, nprims)
		  nprims = 0
		  sym = -10
		  next
		end
		a = line.split
		if a.length == 1
		  i += 1
		  ln += 1  #  Skip the blank line
		  next
		elsif a.length == 5 || a.length == 6
		  if sym == -10
			case a[1]
			when "S"
			  sym = 0; n = 1
			when "P"
			  sym = 1; n = 3
			when "L"
			  sym = -1; n = 4
			when "D"
			  sym = 2; n = 6
			when "F"
			  sym = 3; n = 10
			when "G"
			  sym = 4; n = 15
			else
			  raise MolbyError, "Unknown gaussian shell type at line #{lineno + ln}"
			end
			ncomps += n
		  end
		  if (a.length == 5 && sym == -1) || (a.length == 6 && sym != -1)
			raise MolbyError, "Wrong format in gaussian shell information at line #{lineno + ln}"
		  end
		  exp = Float(a[3])
		  c = Float(a[4])
		  csp = Float(a[5] || 0.0)
		  add_gaussian_primitive_coefficients(exp, c, csp)
		  nprims += 1
		else
		  raise MolbyError, "Error in reading basis set information at line #{lineno + ln}"
		end
	end
	return ncomps
  end

  def sub_load_gamess_log_mo_coefficients(lines, lineno, ncomps)
    ln = 0
	idx = 0
	alpha = true
	while (line = lines[ln]) != nil
		ln += 1
		if line =~ /BETA SET/
			alpha = false
			next
		end
		if line =~ /------------/ || line =~ /EIGENVECTORS/ || line =~ /\*\*\*\* (ALPHA|BETA) SET/
			next
		end
		next unless line =~ /^\s*\d/
		mo_labels = line.split       #  MO numbers (1-based)
		mo_energies = lines[ln].split
		mo_symmetries = lines[ln + 1].split
	#	puts "mo #{mo_labels.inspect}"
		ln += 2
		mo = mo_labels.map { [] }    #  array of *independent* empty arrays
		while (line = lines[ln]) != nil
			ln += 1
			break unless line =~ /^\s*\d/
			5.times { |i|
			  s = line[15 + 11 * i, 11].chomp
			  break if s =~ /^\s*$/
			  mo[i].push(Float(s)) rescue print "line = #{line}, s = #{s}"
			# line[15..-1].split.each_with_index { |s, i|
			#  	mo[i].push(Float(s))
			}
		end
		mo.each_with_index { |m, i|
			idx = Integer(mo_labels[i])
			set_mo_coefficients(idx + (alpha ? 0 : ncomps), Float(mo_energies[i]), m)
		#	if mo_labels[i] % 8 == 1
		#		puts "set_mo_coefficients #{idx}, #{mo_energies[i]}, [#{m[0]}, ..., #{m[-1]}]"
		#	end
		}
#		if line =~ /^\s*$/
#			next
#		else
#			break
#		end
	end
  end
    
  def sub_load_gamess_log(fp)

    if natoms == 0
		new_unit = true
	else
		new_unit = false
    end
	self.update_enabled = false
	mes = "Loading GAMESS log file"
	show_progress_panel(mes)
	energy = nil
	ne_alpha = ne_beta = 0   #  number of electrons
	rflag = nil  #  0, UHF; 1, RHF; 2, ROHF
	mo_count = 0
	search_mode = 0   #  0, no search; 1, optimize; 2, irc
	ncomps = 0   #  Number of AO terms per one MO (sum of the number of components over all shells)
	alpha_beta = nil   #  Flag to read alpha/beta MO appropriately
	nsearch = 0  #  Search number for optimization
	begin
	#	if nframes > 0
	#		create_frame
	#		frame = nframes - 1
	#	end
		n = 0
		while 1
			line = fp.gets
			if line == nil
				break
			end
			line.chomp!
			if line =~ /ATOM\s+ATOMIC\s+COORDINATES/ || line =~ /COORDINATES OF ALL ATOMS ARE/ || line =~ /COORDINATES \(IN ANGSTROM\) FOR \$DATA GROUP ARE/
				set_progress_message(mes + "\nReading atomic coordinates...")
				if line =~ /ATOMIC/
				  first_line = true
				#  if !new_unit
				#    next   #  Skip initial atomic coordinates unless loading into an empty molecule
				#  end
				  line = fp.gets  #  Skip one line
				else
				  first_line = false
				  nsearch += 1
				  if line =~ /COORDINATES OF ALL ATOMS ARE/
				    line = fp.gets  #  Skip one line
			      end
				end
				n = 0
				coords = []
				names = []
				while (line = fp.gets) != nil
					break if line =~ /^\s*$/ || line =~ /END OF ONE/ || line =~ /\$END/
					next if line =~ /-----/
					name, charge, x, y, z = line.split
					v = Vector3D[x, y, z]
					coords.push(v * (first_line ? 0.529177 : 1.0))  #  Bohr to angstrom
					names.push([name, charge])
				end
				if new_unit
					#  Build a new molecule
					names.each_index { |i|
						ap = add_atom(names[i][0])
						ap.atomic_number = names[i][1].to_i
						ap.atom_type = ap.element
						ap.r = coords[i]
					}
					#  Find bonds
					guess_bonds
				#	atoms.each { |ap|
				#		j = ap.index
				#		(j + 1 ... natoms).each { |k|
				#			if calc_bond(j, k) < 1.7
				#				create_bond(j, k)
				#			end
				#		}
				#	}
					new_unit = false
				#	create_frame
				else
				    dont_create = false
					if (search_mode == 1 && nsearch == 1) || first_line
						#  The input coordinate and the first frame for geometry search
						#  can have the same coordinate as the last frame; if this is the case, then
						#  do not create the new frame
						select_frame(nframes - 1)
						dont_create = true
						each_atom { |ap|
						  if (ap.r - coords[ap.index]).length2 > 1e-8
						    dont_create = false
							break
						  end
						}
					end
					if !dont_create
						create_frame([coords])  #  Should not be (coords)
					end
				end
				set_property("energy", energy) if energy
			elsif line =~ /BEGINNING GEOMETRY SEARCH POINT/
				energy = nil   #  New search has begun, so clear the energy
				search_mode = 1
			elsif line =~ /CONSTRAINED OPTIMIZATION POINT/
				energy = nil   #  New search has begun, so clear the energy
				search_mode = 2
			elsif line =~ /FINAL .* ENERGY IS *([-.0-9]+) AFTER/
				if search_mode != 2
					energy = $1.to_f
					set_property("energy", energy)
				end
			elsif line =~ /TOTAL ENERGY += +([-.0-9]+)/
				energy = $1.to_f
			elsif false && line =~ /EQUILIBRIUM GEOMETRY LOCATED/i
				set_progress_message(mes + "\nReading optimized coordinates...")
				fp.gets; fp.gets; fp.gets
				n = 0
				while (line = fp.gets) != nil
					break if line =~ /^\s*$/
					line.chomp
					atom, an, x, y, z = line.split
					ap = atoms[n]
					ap.r = Vector3D[x, y, z]
					n += 1
					break if n >= natoms
				end
			#	if ne_alpha > 0 && ne_beta > 0
			#		#  Allocate basis set record again, to update the atomic coordinates
			#		allocate_basis_set_record(rflag, ne_alpha, ne_beta)
			#	end
			elsif line =~ /ATOMIC BASIS SET/
				lines = []
				lineno = fp.lineno
				while (line = fp.gets)
					break if line =~ /TOTAL NUMBER OF BASIS SET/
					line.chomp!
					lines.push(line)
				end
				ncomps = sub_load_gamess_log_basis_set(lines, lineno)
			elsif line =~ /NUMBER OF OCCUPIED ORBITALS/
				line =~ /=\s*(\d+)/
				n = Integer($1)
				if line =~ /ALPHA/
					ne_alpha = n
				else
					ne_beta = n
				end
			elsif line =~ /SCFTYP=(\w+)/
				scftyp = $1
				if ne_alpha > 0 || ne_beta > 0
					rflag = 0
					case scftyp
					when "RHF"
						rflag = 1
					when "ROHF"
						rflag = 2
					end
				end
			elsif line =~ /(ALPHA|BETA)\s*SET/
				alpha_beta = $1
			elsif line =~ /^\s*(EIGENVECTORS|MOLECULAR ORBITALS)\s*$/
				if mo_count == 0
					clear_mo_coefficients
					set_mo_info(:type=>["UHF", "RHF", "ROHF"][rflag], :alpha=>ne_alpha, :beta=>ne_beta)
				end
				mo_count += 1
				line = fp.gets; line = fp.gets
				lineno = fp.lineno
				lines = []
				set_progress_message(mes + "\nReading MO coefficients...")
				while (line = fp.gets)
					break if line =~ /\.\.\.\.\.\./ || line =~ /----------------/
					line.chomp!
					lines.push(line)
				end
				sub_load_gamess_log_mo_coefficients(lines, lineno, ncomps)
				set_progress_message(mes)
			elsif line =~ /N A T U R A L   B O N D   O R B I T A L   A N A L Y S I S/
				nbo_lines = []
				while (line = fp.gets) != nil
				  break if line =~ /done with NBO analysis/
				  nbo_lines.push(line)
				end
				import_nbo_log(nbo_lines)
			end
		end
	end
	if nframes > 0
	  select_frame(nframes - 1)
	end
	if energy && energy != 0.0
	  set_property("energy", energy)
	end
	hide_progress_panel
	self.update_enabled = true
	(n > 0 ? true : false)
  end
  
  def sub_load_gaussian_log(fp)

    if natoms == 0
		new_unit = true
	else
		new_unit = false
    end
	if nframes > 0
		create_frame
		frame = nframes - 1
	end
	n = 0
	nf = 0
	energy = nil
	use_input_orientation = false
	show_progress_panel("Loading Gaussian out file...")
	while 1
		line = fp.gets
		if line == nil
			break
		end
		line.chomp!
		if line =~ /(Input|Standard) orientation/
			match = $1
			if match == "Input"
				use_input_orientation = true if nf == 0
				next if !use_input_orientation
			else
				next if use_input_orientation
			end
			4.times { line = fp.gets }    #  Skip four lines
			n = 0
			coords = []
			anums = []
			while (line = fp.gets) != nil
				break if line =~ /-----/
				num, charge, type, x, y, z = line.split
				coords.push(Vector3D[x, y, z])
				anums.push(charge)
				n += 1
			end
			if new_unit
				#  Build a new molecule
				anums.each_index { |i|
					ap = add_atom("X")
					ap.atomic_number = anums[i]
					ap.atom_type = ap.element
					ap.name = sprintf("%s%d", ap.element, i)
					ap.r = coords[i]
				}
				#  Find bonds
			#	atoms.each { |ap|
			#		j = ap.index
			#		(j + 1 ... natoms).each { |k|
			#			if calc_bond(j, k) < 1.7
			#				create_bond(j, k)
			#			end
			#		}
			#	}
				guess_bonds
				new_unit = false
			#	create_frame
			else
				create_frame([coords])  #  Should not be (coords)
			end
			if energy
				# TODO: to ensure whether the energy output line comes before
				# or after the atomic coordinates.
				set_property("energy", energy)
			end
			nf += 1
		elsif line =~ /SCF Done: *E\(\w+\) *= *([-.0-9]+)/
			energy = $1.to_f
		end
	end
	if energy
		set_property("energy", energy)
	end
	hide_progress_panel
	(n > 0 ? true : false)
  end

  def loadout(filename)
    retval = false
    fp = open(filename, "rb")
	while s = fp.gets
	  if s =~ /Gaussian/
	    retval = sub_load_gaussian_log(fp)
		break
	  elsif s =~ /GAMESS/
	    retval = sub_load_gamess_log(fp)
		break
	  end
	end
	fp.close
	return retval
  end
  
  alias :loadlog :loadout

  def loadxyz(filename)
	fp = open(filename, "rb")
	n = 0
	coords = []
	names = []
	cell = nil
	while 1
	  line = fp.gets
	  if line == nil
		fp.close
		break
	  end
	  line.chomp
	  toks = line.split
	  if coords.length == 0
	    #  Allow "number of atoms" line or "number of atoms and crystallographic cell parameter"
		#  (Chem3D xyz format)
		next if toks.length == 1
		if toks.length == 7
		  cell = toks[1..6].map { |s| Float(s.sub(/\(\d+\)/, "")) }  #  Remove (xx) and convert to float
		  next
		end
	  end
	  name, x, y, z = line.split
	  next if z == nil
	  x = Float(x.sub(/\(\d+\)/, ""))
	  y = Float(y.sub(/\(\d+\)/, ""))
	  z = Float(z.sub(/\(\d+\)/, ""))
	  r = Vector3D[x, y, z]
	  coords.push(r)
	  names.push(name)
	  n += 1
	end
	celltr = nil
	if cell
	  self.cell = cell
	  celltr = self.cell_transform
	end
	names.each_index { |i|
	  ap = add_atom(names[i])
	  names[i] =~ /^([A-Za-z]{1,2})/
	  element = $1.capitalize
	  ap.element = element
	  ap.atom_type = element
	  if celltr
	    ap.r = celltr * coords[i]
	  else
	    ap.r = coords[i]
	  end
	}
	guess_bonds
	#  Find bonds
#	atoms.each { |ap|
#	  j = ap.index
#	  (j + 1 ... natoms).each { |k|
#		if calc_bond(j, k) < 1.7
#		#  create_bond(j, k)
#		end
#	  }
#	}
#	self.undo_enabled = save_undo_enabled
	(n > 0 ? true : false)
  end
  
  def savexyz(filename)
    open(filename, "wb") { |fp|
	  fp.printf "%d\n", self.natoms
	  each_atom { |ap|
	    fp.printf "%s %.5f %.5f %.5f\n", ap.element, ap.x, ap.y, ap.z
	  }
	}
	return true
  end
  
  def loadzmat(filename)
    self.remove(All)
	open(filename, "rb") { |fp|
	  while (line = fp.gets)
	    line.chomp!
		a = line.split
		an = Molecule.guess_atomic_number_from_name(a[0])
		elm = Parameter.builtin.elements[an].name
		base1 = a[1].to_i - 1
		base2 = a[3].to_i - 1
		base3 = a[5].to_i - 1
		base1 = nil if base1 < 0
		base2 = nil if base2 < 0
		base3 = nil if base3 < 0
		add_atom(a[0], elm, elm, a[2].to_f, base1, a[4].to_f, base2, a[6].to_f, base3)
	  end
	}
	return true
  end
  
  def savezmat(filename)
  end
  
  def loadcom(filename)
#	save_undo_enabled = self.undo_enabled?
#	self.undo_enabled = false
	self.remove(All)
	fp = open(filename, "rb")
	section = 0
	while (line = fp.gets)
	  line.chomp!
	  if section == 0
	    section = 1 if line =~ /^\#/
		next
	  elsif section == 1 || section == 2
	    section += 1 if line =~ /^\s*$/
		next
	  else
	    #  The first line is skipped (charge and multiplicity)
		while (line = fp.gets)
		  line.chomp!
		  break if line =~ /^\s*$/
		  a = line.split(/\s*[ \t,\/]\s*/)
		  r = Vector3D[Float(a[1]), Float(a[2]), Float(a[3])]
		  ap = add_atom(a[0])
		  a[0] =~ /^([A-Za-z]{1,2})/
		  element = $1.capitalize
	      ap.element = element
	      ap.atom_type = element
	      ap.r = r
		end
		break
	  end
	end
	fp.close
	guess_bonds
#	self.undo_enabled = save_undo_enabled
	return true
  end
  
  def loadinp(filename)
#	save_undo_enabled = self.undo_enabled?
#	self.undo_enabled = false
	self.remove(All)
	fp = open(filename, "rb")
	section = 0
	has_basis = false
	while (line = fp.gets)
	  if line =~ /\A \$BASIS/
	    has_basis = true
		next
	  end
	  next if line !~ /\A \$DATA/
	  line = fp.gets   #  Title line
	  line = fp.gets   #  Symmetry line
	  while (line = fp.gets)
#	    puts line
	    line.chomp!
		break if line =~ /\$END/
		a = line.split
		ap = add_atom(a[0])
		ap.atomic_number = Integer(a[1])
		ap.atom_type = ap.element
		r = Vector3D[Float(a[2]), Float(a[3]), Float(a[4])]
		ap.r = r;
		if !has_basis
		  #  Skip until one blank line is detected
		  while (line = fp.gets) && line =~ /\S/
		  end
		end
      end
	  break
	end
	fp.close
	guess_bonds
#	self.undo_enabled = save_undo_enabled
	return true
  end
  
  def saveinp(filename)
    if natoms == 0
      raise MolbyError, "cannot save GAMESS input; the molecule is empty"
    end
    fp = open(filename, "wb")
	now = Time.now.to_s
	fp.print <<end_of_header
!  GAMESS input
!  Generated by Molby at #{now}
 $CONTRL COORD=UNIQUE EXETYP=RUN ICHARG=0
         ICUT=20 INTTYP=HONDO ITOL=30
         MAXIT=200 MOLPLT=.T. MPLEVL=0
         MULT=1 QMTTOL=1e-08 RUNTYP=OPTIMIZE
         SCFTYP=RHF UNITS=ANGS                  $END
 $SCF    CONV=1.0E-06 DIRSCF=.T.                $END
 $STATPT NSTEP=400 OPTTOL=1.0E-06               $END
 $SYSTEM MEMDDI=0 MWORDS=16 TIMLIM=50000        $END
 $BASIS  GBASIS=N31 NDFUNC=1 NGAUSS=6           $END
 $GUESS  GUESS=HUCKEL                           $END
!
 $DATA
 #{name}
 C1
end_of_header
	each_atom { |ap|
		next if ap.atomic_number == 0
		fp.printf " %-6s %4d %10.6f %10.6f %10.6f\n", ap.name, ap.atomic_number, ap.r.x, ap.r.y, ap.r.z
	}
	fp.print " $END\n"
	fp.close
	return true
  end

  def savecom(filename)
    if natoms == 0
      raise MolbyError, "cannot save Gaussian input; the molecule is empty"
    end
    fp = open(filename, "wb")
	base = File.basename(filename, ".*")
	fp.print <<end_of_header
%Chk=#{base}.chk
\# PM3 Opt

 #{name}; created by Molby at #{Time.now.to_s}

 0 1
end_of_header
	each_atom { |ap|
	    next if ap.atomic_number == 0
		fp.printf "%-6s %10.6f %10.6f %10.6f\n", ap.element, ap.r.x, ap.r.y, ap.r.z
	}
	fp.print "\n"
	fp.close
	return true
  end

  alias :loadgjf :loadcom
  alias :savegjf :savecom
  
  def loadcif(filename)
    mol = self
    def getciftoken(fp)
	  while @tokens.length == 0
	    line = fp.gets
	    return nil if !line
		if line[0] == ?;
		  s = line
		  while line = fp.gets
		    break if line[0] == ?;
		    s += line
		  end
		  return s if !line
		  s += ";"
		  @tokens.push(s)
		  line = line[1...line.length]
		end
	    line.strip!
	    line.scan(/'[^\']*'|"[^\"]*"|[^#\s]+|#.*/) do |s|
	      next if s[0] == ?#
		  if s =~ /^data_|loop_|global_|save_|stop_/i
		    s = "#" + s    #  Label for reserved words
		  end
		  @tokens.push(s)
	    end
	  end
	  return @tokens.shift
	end
	def float_strip_rms(str)
	  str =~ /^(-?)(\d*)(\.(\d*))?(\((\d+)\))?$/
	  sgn, i, frac, rms = $1, $2, $4, $6
	  i = i.to_f
	  if frac
		base = 0.1 ** frac.length
		i = i + frac.to_f * base
	  else
	    base = 1.0
	  end
	  if rms
	    rms = rms.to_f * base
	  else
	    rms = 0.0
	  end
	  if sgn == "-"
	    i = -i
	  end
	  return i, rms
	end
	def parse_symmetry_operation(str)
	  if str == "."
	    sym = nil
	  elsif (str =~ /(\d+)_(\d)(\d)(\d)/) || (str =~ /(\d+) +(\d)(\d)(\d)/)
	    sym = [Integer($1) - 1, Integer($2) - 5, Integer($3) - 5, Integer($4) - 5]
	  elsif (str =~ /^(\d+)$/)
	    sym = [Integer($1) - 1, 0, 0, 0]
	  end
	  if sym && (sym[0] == 0 && sym[1] == 0 && sym[2] == 0 && sym[3] == 0)
	    sym = nil
	  end
	  sym
	end
	def find_atom_by_name(mol, name)
	  name = name.delete(" ()")
	  ap = mol.atoms[name] rescue ap = nil
	  return ap
	end
	selfname = self.name
	fp = open(filename, "rb")
	data_identifier = nil
	@tokens = []
	count_up = 1
	while true
	  warn_message = ""
	  verbose = nil
      verbose = true
	  bond_defined = false
	  special_positions = []
	  mol.remove(All)
	  cell = []
	  cell_trans = cell_trans_inv = Transform.identity
	  token = getciftoken(fp)
	  pardigits_re = /\(\d+\)/
	  calculated_atoms = []
	  while token != nil
	    if token =~ /^\#data_/i
		  if data_identifier == nil || mol.natoms == 0
		    #  First block or no atoms yet
            #  Continue processing of this molecule
		    data_identifier = token
			token = getciftoken(fp)
			next
		  else
		    #  Description of another molecule begins here
			data_identifier = token
			break
		  end
	    elsif token =~ /^_cell/
		  val = getciftoken(fp)
		  if token == "_cell_length_a"
		    cell[0], cell[6] = float_strip_rms(val)
		  elsif token == "_cell_length_b"
		    cell[1], cell[7] = float_strip_rms(val)
		  elsif token == "_cell_length_c"
		    cell[2], cell[8] = float_strip_rms(val)
		  elsif token == "_cell_angle_alpha"
		    cell[3], cell[9] = float_strip_rms(val)
		  elsif token == "_cell_angle_beta"
		    cell[4], cell[10] = float_strip_rms(val)
		  elsif token == "_cell_angle_gamma"
		    cell[5], cell[11] = float_strip_rms(val)
		  end
		  if cell.length == 12 && cell.all?
		    mol.cell = cell
            puts "Unit cell is set to #{mol.cell.inspect}." if verbose
		    cell = []
		    cell_trans = mol.cell_transform
		    cell_trans_inv = cell_trans.inverse
		  end
		  token = getciftoken(fp)
		  next
        elsif token.casecmp("#loop_") == 0
	      labels = []
		  while (token = getciftoken(fp)) && token[0] == ?_
		    labels.push(token)
		  end
		  if labels[0] =~ /symmetry_equiv_pos|space_group_symop|atom_site_label|atom_site_aniso_label|geom_bond/
		    hlabel = Hash.new(-10000000)
		    labels.each_with_index { |lb, i|
		  	  hlabel[lb] = i
		    }
		    data = []
		    n = labels.length
		    a = []
		    while true
		   	  break if token == nil || token[0] == ?_ || token[0] == ?#
			  a.push(token)
			  if a.length == n
			    data.push(a)
			    a = []
			  end
			  token = getciftoken(fp)
		    end
		    if labels[0] =~ /^_symmetry_equiv_pos/ || labels[0] =~ /^_space_group_symop/
		      data.each { |d|
			    symstr = d[hlabel["_symmetry_equiv_pos_as_xyz"]] || d[hlabel["_space_group_symop_operation_xyz"]]
			    symstr.delete("\"\'")
			    exps = symstr.split(/,/)
			    sym = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
			    exps.each_with_index { |s, i|
			      terms = s.scan(/([-+]?)(([.0-9]+)(\/([0-9]+))?([xXyYzZ])?|([xXyYzZ]))/)
				  terms.each { |a|
				    #  a[0]: sign, a[2]: numerator, a[4]: denometer
				    if a[4] != nil
				      #  The number part is a[2]/a[4]
				      num = Float(a[2])/Float(a[4])
				    elsif a[2] != nil
				      #  The number part is either integer or a floating point
				      num = Float(a[2])
				    else
				      num = 1.0
				    end
				    num = -num if a[0][0] == ?-
				    xyz = (a[5] || a[6])
				    if xyz == "x" || xyz == "X"
				      sym[i] = num
				    elsif xyz == "y" || xyz == "Y"
				      sym[i + 3] = num
				    elsif xyz == "z" || xyz == "Z"
				      sym[i + 6] = num
				    else
				      sym[9 + i] = num
				    end
				  }
			    }
			    puts "symmetry operation #{sym.inspect}" if verbose
			    mol.add_symmetry(Transform.new(sym))
			  }
			  puts "#{mol.nsymmetries} symmetry operations are added" if verbose
		    elsif labels[0] =~ /^_atom_site_label/
			  #  Create atoms
			  data.each { |d|
			    name = d[hlabel["_atom_site_label"]]
			    elem = d[hlabel["_atom_site_type_symbol"]]
			    fx = d[hlabel["_atom_site_fract_x"]]
			    fy = d[hlabel["_atom_site_fract_y"]]
			    fz = d[hlabel["_atom_site_fract_z"]]
			    uiso = d[hlabel["_atom_site_U_iso_or_equiv"]]
			    biso = d[hlabel["_atom_site_B_iso_or_equiv"]]
			    occ = d[hlabel["_atom_site_occupancy"]]
			    calc = d[hlabel["_atom_site_calc_flag"]]
			    name = name.delete(" ()")
			    if elem == nil || elem == ""
			      if name =~ /[A-Za-z]{1,2}/
				    elem = $&.capitalize
				  else
				    elem = "Du"
				  end
			    end
			    ap = mol.add_atom(name, elem, elem)
			    ap.fract_x, ap.sigma_x = float_strip_rms(fx)
			    ap.fract_y, ap.sigma_y = float_strip_rms(fy)
			    ap.fract_z, ap.sigma_z = float_strip_rms(fz)
			    if biso
			      ap.temp_factor, sig = float_strip_rms(biso)
			    elsif uiso
			      ap.temp_factor, sig = float_strip_rms(uiso)
				  ap.temp_factor *= 78.9568352087149          #  8*pi*pi
			    end
			    ap.occupancy, sig = float_strip_rms(occ)
			    if calc == "c" || calc == "calc"
			      calculated_atoms.push(ap.index)
		        end
			    #  Guess special positions
			    (1...mol.nsymmetries).each { |isym|
			      sr = ap.fract_r
			      sr = (mol.transform_for_symop(isym) * sr) - sr;
				  nx = (sr.x + 0.5).floor
				  ny = (sr.y + 0.5).floor
				  nz = (sr.z + 0.5).floor
				  if (Vector3D[sr.x - nx, sr.y - ny, sr.z - nz].length2 < 1e-6)
				    #  [isym, -nx, -ny, -nz] transforms this atom to itself
				    #  The following line is equivalent to:
				    #    if special_positions[ap.index] == nil; special_positions[ap.index] = []; end;
				    #    special_positions[ap.index].push(...)
				    (special_positions[ap.index] ||= []).push([isym, -nx, -ny, -nz])
				  end
			    }
			    if verbose && special_positions[ap.index]
			      puts "#{name} is on the special position: #{special_positions[ap.index].inspect}"
			    end
			  }
			  puts "#{mol.natoms} atoms are created." if verbose
		    elsif labels[0] =~ /^_atom_site_aniso_label/
		      #  Set anisotropic parameters
			  c = 0
			  data.each { |d|
			    name = d[hlabel["_atom_site_aniso_label"]]
			    ap = find_atom_by_name(mol, name)
			    next if !ap
			    u11 = d[hlabel["_atom_site_aniso_U_11"]]
			    if u11
			      usig = []
			      u11, usig[0] = float_strip_rms(u11)
			      u22, usig[1] = float_strip_rms(d[hlabel["_atom_site_aniso_U_22"]])
			      u33, usig[2] = float_strip_rms(d[hlabel["_atom_site_aniso_U_33"]])
			      u12, usig[3] = float_strip_rms(d[hlabel["_atom_site_aniso_U_12"]])
			      u13, usig[4] = float_strip_rms(d[hlabel["_atom_site_aniso_U_13"]])
			      u23, usig[5] = float_strip_rms(d[hlabel["_atom_site_aniso_U_23"]])
			      ap.aniso = [u11, u22, u33, u12, u13, u23, 8] + usig
				  c += 1
			    end
			  }
			  puts "#{c} anisotropic parameters are set." if verbose
		    elsif labels[0] =~ /^_geom_bond/
		      #  Create bonds
			  exbonds = []
			  data.each { |d|
			    n1 = d[hlabel["_geom_bond_atom_site_label_1"]]
			    n2 = d[hlabel["_geom_bond_atom_site_label_2"]]
			    sym1 = d[hlabel["_geom_bond_site_symmetry_1"]] || "."
			    sym2 = d[hlabel["_geom_bond_site_symmetry_2"]] || "."
			    n1 = find_atom_by_name(mol, n1)
			    n2 = find_atom_by_name(mol, n2)
			    next if n1 == nil || n2 == nil
		        n1 = n1.index
			    n2 = n2.index
			    sym1 = parse_symmetry_operation(sym1)
			    sym2 = parse_symmetry_operation(sym2)
			    if sym1 || sym2
			      exbonds.push([n1, n2, sym1, sym2])
			    else
			      mol.create_bond(n1, n2)
			    end
			    tr1 = (sym1 ? mol.transform_for_symop(sym1) : Transform.identity)
			    tr2 = (sym2 ? mol.transform_for_symop(sym2) : Transform.identity)
			    if special_positions[n1]
				  #  Add extra bonds for equivalent positions of n1
				  special_positions[n1].each { |symop|
				    sym2x = mol.symop_for_transform(tr1 * mol.transform_for_symop(symop) * tr1.inverse * tr2)
				    exbonds.push([n1, n2, sym1, sym2x])
				  }
			    end
			    if special_positions[n2]
				  #  Add extra bonds n2-n1.symop, where symop transforms n2 to self
				  tr = (sym1 ? mol.transform_for_symop(sym1) : Transform.identity)
				  special_positions[n2].each { |symop|
				    sym1x = mol.symop_for_transform(tr2 * mol.transform_for_symop(symop) * tr2.inverse * tr1)
				    exbonds.push([n2, n1, sym2, sym1x])
				  }
			    end				
		      }
              if mol.nbonds > 0
                bond_defined = true
              end
			  puts "#{mol.nbonds} bonds are created." if verbose
			  if calculated_atoms.length > 0
			    #  Guess bonds for calculated hydrogen atoms
			    n1 = 0
			    calculated_atoms.each { |ai|
			      if mol.atoms[ai].connects.length == 0
				    as = mol.find_close_atoms(ai)
				    as.each { |aj|
				      mol.create_bond(ai, aj)
					  n1 += 1
				    }
				  end
			    }
			    puts "#{n1} bonds are guessed." if verbose
			  end
			  if exbonds.length > 0
			    h = Dialog.run("CIF Import: Symmetry Expansion") {
			      layout(1,
				    item(:text, :title=>"There are bonds including symmetry related atoms.\nWhat do you want to do?"),
				    item(:radio, :title=>"Expand only atoms that are included in those extra bonds.", :tag=>"atoms_only"),
				    item(:radio, :title=>"Expand fragments having atoms included in the extra bonds.", :tag=>"fragment", :value=>1),
				    item(:radio, :title=>"Ignore these extra bonds.", :tag=>"ignore")
				  )
				  radio_group("atoms_only", "fragment", "ignore")
			    }
			    if h[:status] == 0 && h["ignore"] == 0
			      atoms_only = (h["atoms_only"] != 0)
				  if !atoms_only
				    fragments = []
				    mol.each_fragment { |f| fragments.push(f) }
				  end
				  debug = nil
				  exbonds.each { |ex|
				    if debug; puts "extra bond #{ex[0]}(#{ex[2].inspect}) - #{ex[1]}(#{ex[3].inspect})"; end
				    ex0 = ex.dup
				    (2..3).each { |i|
				      symop = ex[i]
					  if symop == nil
					    ex[i + 2] = ex[i - 2]
					  else
					    if debug; puts "  symop = #{symop.inspect}"; end
					    #  Expand the atom or the fragment including the atom
					    if atoms_only
						  ig = IntGroup[ex[i - 2]]
						  idx = 0
					    else
						  ig = fragments.find { |f| f.include?(ex[i - 2]) }
						  ig.each_with_index { |n, ii| if n == ex[i - 2]; idx = ii; break; end }
					    end
					    symop[4] = ex[i - 2]  #  Base atom
					    if debug; puts "  expanding #{ig} by #{symop.inspect}"; end
					    a = mol.expand_by_symmetry(ig, symop[0], symop[1], symop[2], symop[3])
					    ex[i + 2] = a[idx]   #  Index of the expanded atom
					  end
				    }
				    if ex[4] && ex[5] && ex[4] != ex[5]
				      if debug; puts "  creating bond #{ex[4]} - #{ex[5]}"; end
				      mol.create_bond(ex[4], ex[5])
				    end
				  }
			    end
			  end
			  puts "#{mol.nbonds} bonds are created." if verbose
		    end
		    next
		  else
		  #  puts "Loop beginning with #{labels[0]} is skipped"
		  end
	    else
	      #  Skip this token
		  token = getciftoken(fp)
	    end
	    #  Skip tokens until next tag or reserved word is detected
	    while token != nil && token[0] != ?_ && token[0] != ?#
		  token = getciftoken(fp)
	    end
	    next
	  end
      if !bond_defined
		mol.guess_bonds
	  end
	  if token != nil && token == data_identifier
	    #  Process next molecule: open a new molecule and start adding atom on that
		mol = Molecule.new
		count_up += 1
		(@aux_mols ||= []).push(mol)
		next
	  end
	  break
	end
	fp.close
#	self.undo_enabled = save_undo_enabled
	return true
  end
  
  def dump(group = nil)
    def quote(str)
	  if str == ""
	    str = "\"\""
	  else
	    str = str.gsub(/%/, "%%")
		str.gsub!(/ /, "%20")
		str.gsub!(/\t/, "%09")
	  end
	  str
	end
	group = atom_group(group ? group : 0...natoms)
	s = ""
	group.each { |i|
	  ap = atoms[i]
	  s += sprintf("%4d %-7s %-4s %-4s %-2s %7.3f %7.3f %7.3f %6.3f [%s]\n",
		ap.index, sprintf("%3s.%d", ap.res_name, ap.res_seq),
		quote(ap.name), quote(ap.atom_type), quote(ap.element),
		ap.r.x, ap.r.y, ap.r.z, ap.charge,
		ap.connects.join(","))
	}
	print s
  end

  def from_dump(arg)
    if natoms > 0
	  raise "Molecule must be empty"
	end
    format = "index residue name atom_type element rx ry rz charge connects"
	keys = []
	resAtoms = Hash.new
	newBonds = []
	#  arg can be either a String or an array of String.
	#  Iterates for each line in the string or each member of the array.
	if arg.is_a?(String)
	  arg = arg.split("\n")
	end
    arg.each { |line|
	  if line =~ /^\#/
	    format = line[1..-1]
		keys = []
	  end
	  if keys.length == 0
	    keys = format.split(" ").collect { |key| key.to_sym }
	  end
	  values = line.chomp.split(" ")
	  next if values == nil || values.length == 0
	  ap = create_atom(sprintf("X%03d", natoms))
	  r = Vector3D[0, 0, 0]
	  keys.each_index { |i|
	    break if (value = values[i]) == nil
		if value == "\"\""
		  value = ""
		else
		  value.gsub(/%09/, "\t")
		  value.gsub(/%20/, " ")
		  value.gsub(/%%/, "%")
		end
		key = keys[i]
		if key == :residue
		  if resAtoms[value] == nil
		    resAtoms[value] = []
		  end
		  resAtoms[value].push(ap.index)
		elsif key == :rx
		  r.x = value.to_f
		elsif key == :ry
		  r.y = value.to_f
		elsif key == :rz
		  r.z = value.to_f
		elsif key == :connects
		  value.scan(/\d+/).each { |i|
		    i = i.to_i
		    if ap.index < i
			  newBonds.push(ap.index)
			  newBonds.push(i)
			end
		  }
		elsif key == :index
		  next
		else
		  ap.set_attr(key, value)
		end
	  }
	  ap.r = r
	}
	resAtoms.each_key { |key|
	  assign_residue(atom_group(resAtoms[key]), key)
	}
	create_bond(*newBonds)
	self
  end

  #  Plug-in for loading mbsf
  def loadmbsf_plugin(s, lineno)
    ""
  end
  
  #  Plug-in for saving mbsf
  def savembsf_plugin
    ""
  end
  
end
