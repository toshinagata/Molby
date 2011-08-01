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
    periodic = (self.box && self.box[4].all? { |n| n != 0 })
	show_progress_panel("Loading AMBER crd file...")
#    puts "sframe = #{sframe}, pos = #{fp.pos}"
    line = fp.gets   #  Skip first line
    while 1
      line = fp.gets
      if line == nil
	    if (count > 0)
		  raise MolbyError, sprintf("crd format error - file ended in the middle of frame %d", frame)
		end
	    fp.close
		break
      end
#      next if line.match(/^TITLE/)
      line.chomp!
      values = line.split(' ')
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
		  if line == nil || (values = line.chomp.split(' ')).length != 3
		    raise "The molecule has a periodic cell but the crd file does not contain cell information"
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

  def loadlog(filename)

    if natoms == 0
		new_unit = true
	else
		new_unit = false
    end
#	save_undo_enabled = self.undo_enabled?
#	self.undo_enabled = false
	self.update_enabled = false
	mes = "Loading GAMESS log file"
	show_progress_panel(mes)
	ne_alpha = ne_beta = 0   #  number of electrons
	rflag = nil  #  0, UHF; 1, RHF; 2, ROHF
	mo_count = 0
	ncomps = 0   #  Number of AO terms per one MO (sum of the number of components over all shells)
	alpha_beta = nil   #  Flag to read alpha/beta MO appropriately
	begin
		fp = open(filename, "rb")
		if nframes > 0
			create_frame
			frame = nframes - 1
		end
		n = 0
		while 1
			line = fp.gets
			if line == nil
				fp.close
				break
			end
			line.chomp!
			if line =~ /ATOM\s+ATOMIC\s+COORDINATES/ || line =~ /COORDINATES OF ALL ATOMS ARE/
				set_progress_message(mes + "\nReading atomic coordinates...")
				first_line = (line =~ /ATOMIC/)
				line = fp.gets    #  Skip one line
				n = 0
				coords = []
				names = []
				while (line = fp.gets) != nil
					break if line =~ /^\s*$/ || line =~ /END OF ONE/
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
					create_frame
				else
					create_frame([coords])  #  Should not be (coords)
				end
			elsif line =~ /EQUILIBRIUM GEOMETRY LOCATED/i
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
				if ne_alpha > 0 && ne_beta > 0
					#  Allocate basis set record again, to update the atomic coordinates
					allocate_basis_set_record(rflag, ne_alpha, ne_beta)
				end
			elsif line =~ /ATOMIC BASIS SET/
				while (line = fp.gets)
					break if line =~ /SHELL\s+TYPE\s+PRIMITIVE/
				end
				line = fp.gets
				i = -1
				nprims = 0
				sym = -10  #  undefined
				ncomps = 0
				while (line = fp.gets)
					break if line =~ /TOTAL NUMBER OF BASIS SET/
					line.chomp!
					if line =~ /^\s*$/
					  #  End of one shell
					  add_gaussian_orbital_shell(sym, nprims, i)
					  # puts "add_gaussian_orbital_shell #{sym}, #{nprims}, #{i}"
					  nprims = 0
					  sym = -10
					  next
					end
					a = line.split
					if a.length == 1
					  i += 1
					  line = fp.gets  #  Skip the blank line
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
						else
						  raise MolbyError, "Unknown gaussian shell type at line #{fp.lineno}"
						end
						ncomps += n
					  end
					  if (a.length == 5 && sym == -1) || (a.length == 6 && sym != -1)
					    raise MolbyError, "Wrong format in gaussian shell information at line #{fp.lineno}"
					  end
					  exp = Float(a[3])
					  c = Float(a[4])
					  csp = Float(a[5] || 0.0)
					  add_gaussian_primitive_coefficients(exp, c, csp)
					  nprims += 1
					  # puts "add_gaussian_primitive_coefficients #{exp}, #{c}, #{csp}"
					else
					  raise MolbyError, "Error in reading basis set information at line #{fp.lineno}"
					end
				end
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
				if ne_alpha > 0 && ne_beta > 0
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
					allocate_basis_set_record(rflag, ne_alpha, ne_beta)
					# puts "allocate_basis_set_record  #{rflag}, #{ne_alpha}, #{ne_beta}"
				end
				mo_count += 1
				idx = 0
				line = fp.gets; line = fp.gets;
				set_progress_message(mes + "\nReading MO coefficients...")
				while (line = fp.gets) != nil
					break unless line =~ /^\s*\d/
					mo_labels = line.split       #  MO numbers (1-based)
					mo_energies = fp.gets.split
					mo_symmetries = fp.gets.split
					mo = mo_labels.map { [] }    #  array of *independent* empty arrays
					while (line = fp.gets) != nil
						break unless line =~ /^\s*\d/
						line[14..-1].split.each_with_index { |s, i|
							mo[i].push(Float(s))
						}
					end
					mo.each_with_index { |m, i|
						idx = Integer(mo_labels[i]) - 1
						set_mo_coefficients(idx + (alpha_beta == "BETA" ? ncomps : 0), Float(mo_energies[i]), m)
					#	if mo_labels[i] % 8 == 1
					#		puts "set_mo_coefficients #{idx}, #{mo_energies[i]}, [#{m[0]}, ..., #{m[-1]}]"
					#	end
					}
					if line =~ /^\s*$/ && idx < ncomps - 1
						next
					else
						break
					end
				end
				set_progress_message(mes)
			end
		end
#	ensure
#		self.undo_enabled = save_undo_enabled
#        hide_progress_panel
#		self.update_enabled = true
	end
	if nframes > 0
	  select_frame(nframes - 1)
	end
	hide_progress_panel
	self.update_enabled = true
	(n > 0 ? true : false)
  end
  
  def loadxyz(filename)
#	save_undo_enabled = self.undo_enabled?
#	self.undo_enabled = false
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
  
  def loadout(filename)

    if natoms == 0
		new_unit = true
 #     raise MolbyError, "cannot load crd; the molecule is empty"
	else
		new_unit = false
    end
#	save_undo_enabled = undo_enabled?
#	self.undo_enabled = false
	fp = open(filename, "rb")
	if nframes > 0
		create_frame
		frame = nframes - 1
	end
	n = 0
	nf = 0
	use_input_orientation = false
	show_progress_panel("Loading Gaussian out file...")
	while 1
		line = fp.gets
		if line == nil
			fp.close
			break
		end
		line.chomp
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
				create_frame
			else
				create_frame([coords])  #  Should not be (coords)
			end
			nf += 1
		end
	end
	hide_progress_panel
#	self.undo_enabled = save_undo_enabled
	(n > 0 ? true : false)
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
		fp.printf "%-6s %10.6f %10.6f %10.6f\n", ap.element, ap.r.x, ap.r.y, ap.r.z
	}
	fp.print "\n"
	fp.close
	return true
  end

  alias :loadgjf :loadcom
  alias :savegjf :savecom
  
  def loadcif(filename)
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
	  return Float(str.sub(/\(\d+\)/, ""))
	end
	@tokens = []
	self.remove(All)
	fp = open(filename, "rb")
	cell = []
	token = getciftoken(fp)
	pardigits_re = /\(\d+\)/
	while token != nil
	  if token =~ /^_cell/
		val = getciftoken(fp)
		if token == "_cell_length_a"
		  cell[0] = float_strip_rms(val)
		elsif token == "_cell_length_b"
		  cell[1] = float_strip_rms(val)
		elsif token == "_cell_length_c"
		  cell[2] = float_strip_rms(val)
		elsif token == "_cell_angle_alpha"
		  cell[3] = float_strip_rms(val)
		elsif token == "_cell_angle_beta"
		  cell[4] = float_strip_rms(val)
		elsif token == "_cell_angle_gamma"
		  cell[5] = float_strip_rms(val)
		end
		if cell.length == 6 && cell.all?
		  self.cell = cell
		  puts "Unit cell is set to #{cell.inspect}."
		  cell = []
		end
		token = getciftoken(fp)
		next
      elsif token.casecmp("#loop_") == 0
	    labels = []
		while (token = getciftoken(fp)) && token[0] == ?_
		  labels.push(token)
		end
		if labels[0] =~ /symmetry_equiv_pos|atom_site_label|atom_site_aniso_label|geom_bond/
		  hlabel = Hash.new(-10000000)
		  labels.each_with_index { |lb, i|
			hlabel[lb] = i
		  }
		  data = []
		  n = labels.length
		  a = []
		  while 1
			break if token == nil || token[0] == ?_ || token[0] == ?#
			a.push(token)
			if a.length == n
			  data.push(a)
			  a = []
			end
			token = getciftoken(fp)
		  end
		  if labels[0] =~ /^_symmetry_equiv_pos/
		    data.each { |d|
			  symstr = d[hlabel["_symmetry_equiv_pos_as_xyz"]]
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
				    sym[i * 3] = num
				  elsif xyz == "y" || xyz == "Y"
				    sym[i * 3 + 1] = num
				  elsif xyz == "z" || xyz == "Z"
				    sym[i * 3 + 2] = num
				  else
				    sym[9 + i] = num
				  end
				}
			  }
			  puts "symmetry operation #{sym.inspect}"
			  add_symmetry(Transform.new(sym))
			}
			puts "#{self.nsymmetries} symmetry operations are added"
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
			  ap = self.add_atom(name, elem, elem)
			  ap.fract_x = float_strip_rms(fx)
			  ap.fract_y = float_strip_rms(fy)
			  ap.fract_z = float_strip_rms(fz)
			  if biso
			    ap.temp_factor = float_strip_rms(biso)
			  elsif uiso
			    ap.temp_factor = float_strip_rms(uiso) * 78.9568352087149 #  8*pi*pi
			  end
			  ap.occupancy = float_strip_rms(occ)
			}
			puts "#{self.natoms} atoms are created."
		  elsif labels[0] =~ /^_atom_site_aniso_label/
		    #  Set anisotropic parameters
			c = 0
			data.each { |d|
			  name = d[hlabel["_atom_site_aniso_label"]]
			  ap = self.atoms[name]
			  next if !ap
			  u11 = d[hlabel["_atom_site_aniso_U_11"]]
			  if u11
			    u11 = float_strip_rms(u11)
			    u22 = float_strip_rms(d[hlabel["_atom_site_aniso_U_22"]])
			    u33 = float_strip_rms(d[hlabel["_atom_site_aniso_U_33"]])
			    u12 = float_strip_rms(d[hlabel["_atom_site_aniso_U_12"]])
			    u13 = float_strip_rms(d[hlabel["_atom_site_aniso_U_13"]])
			    u23 = float_strip_rms(d[hlabel["_atom_site_aniso_U_23"]])
			    ap.aniso = [u11, u22, u33, u12, u13, u23, 8]
				c += 1
			  end
			}
			puts "#{c} anisotropic parameters are set."
		  elsif labels[0] =~ /^_geom_bond/
		    #  Create bonds
			exbonds = []
			data.each { |d|
			  n1 = d[hlabel["_geom_bond_atom_site_label_1"]]
			  n2 = d[hlabel["_geom_bond_atom_site_label_2"]]
			  sym1 = d[hlabel["_geom_bond_site_symmetry_1"]]
			  sym2 = d[hlabel["_geom_bond_site_symmetry_2"]]
			  if sym1 != "." || sym2 != "."
			    exbonds.push([n1, n2, sym1, sym2])
			  else
			    self.create_bond(n1, n2)
			  end
		    }
			if exbonds.length > 0
			  h = Dialog.run {
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
				  self.each_fragment { |f| fragments.push(f) }
				end
				sym_decoder = /(\d+)_(\d)(\d)(\d)/
				debug = nil
				exbonds.each { |ex|
				  #  Convert name to index before expansion
				  if debug
				    ex[4] = ex[0]
					ex[5] = ex[1]
				  end
				  ex[0] = self.atoms[ex[0]].index
				  ex[1] = self.atoms[ex[1]].index
				}
				exbonds.each { |ex|
				  if debug; puts "extra bond #{ex[4]}(#{ex[2]}) - #{ex[5]}(#{ex[3]})"; end
				  (2..3).each { |i|
				    if ex[i] == "."
					  ex[i] = ex[i - 2]    #  No expansion
					elsif ex[i] =~ /(\d+)_(\d)(\d)(\d)/
			          symop = [Integer($1) - 1, Integer($2) - 5, Integer($3) - 5, Integer($4) - 5, ex[i - 2]]
					  if debug; puts "  symop = #{symop.inspect}"; end
					  ap = self.atoms.find { |ap| (s = ap.symop) != nil && s === symop }
					  if ap
					    if debug; puts "  already expanded (atom #{ap.index}, #{ap.name}, #{ap.symop.inspect})"; end
					    ex[i] = ap.index   #  Already expanded
					  else
					    #  Expand the atom or the fragment including the atom
						if atoms_only
						  ig = IntGroup[ex[i - 2]]
						else
						  ig = fragments.find { |f| f.include?(ex[i - 2]) }
						end
					    if debug; puts "  expanding #{ig} by #{symop.inspect}"; end
						self.expand_by_symmetry(ig, symop[0], symop[1], symop[2], symop[3])
						#  Find again the expanded atom
					    ap = self.atoms.find { |ap| (s = ap.symop) != nil && s === symop }
						ex[i] = ap.index
					  end
					else
					  raise "unrecognizable symmetry operation: #{ex[i]}"
					end
				  }
				  if debug; puts "  creating bond #{ex[2]} - #{ex[3]}"; end
				  self.create_bond(ex[2], ex[3])
				}
			  end
			end
			puts "#{self.nbonds} bonds are created."
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
    arg.each { |line|
      #  arg can be either a String or an array of String. If it is a string,
	  #  arg.each iterates for each line in the string. If it is an array,
	  #  arg.each iterates for each member of the array.
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

end
