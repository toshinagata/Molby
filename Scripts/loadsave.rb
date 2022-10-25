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

  def sub_load_psi4_log(fp)
    if natoms == 0
      new_unit = true
    else
      new_unit = false
    end
    n = 0
    nf = 0
    energy = nil
  
    show_progress_panel("Loading Psi4 output file...")

    getline = lambda { @lineno += 1; return fp.gets }

    #  Import coordinates and energies
    vecs = []
    ats = []
    first_frame = nframes
    trans = nil
    hf_type = nil
    nalpha = nil
    nbeta = nil
    while line = getline.call
      if line =~ /==> Geometry <==/
        #  Skip until line containing "------"
        while line = getline.call
          break if line =~ /------/
        end
        vecs.clear
        index = 0
        #  Read atom positions
        while line = getline.call
          line.chomp!
          break if line =~ /^\s*$/
          tokens = line.split(' ')
          if natoms > 0 && first_frame == nframes
            if index >= natoms || tokens[0].upcase != atoms[index].element.upcase
              hide_progress_panel
              raise MolbyError, "The atom list does not match the current structure at line #{@lineno}"
            end
          end
          vecs.push(Vector3D[Float(tokens[1]), Float(tokens[2]), Float(tokens[3])])
          if natoms == 0
            ats.push(tokens[0])
          end
          index += 1
        end
        if natoms == 0
          #  Create molecule from the initial geometry
          ats.each_with_index { |aname, i|
            #  Create atoms
            ap = add_atom(aname)
            ap.element = aname
            ap.atom_type = ap.element
            ap.name = sprintf("%s%d", aname, i)
            ap.r = vecs[i]
          }
          guess_bonds
        else
          if vecs.length != natoms
            break  #  Log file is incomplete
          end
          #  Does this geometry differ from the last one?
          vecs.length.times { |i|
            if (atoms[i].r - vecs[i]).length2 > 1.0e-14
              #  Create a new frame and break
              create_frame
              vecs.length.times { |j|
                atoms[j].r = vecs[j]
              }
              break
            end
          }
        end
        #  end geometry
      elsif line =~ /Final Energy: +([-.0-9]+)/
        #  Energy for this geometry
        energy = Float($1)
        set_property("energy", energy)
        if line =~ /RHF/
          hf_type = "RHF"
        elsif line =~ /UHF/
          hf_type = "UHF"
        elsif line =~ /ROHF/
          hf_type = "ROHF"
        end
      elsif line =~ /^ *Nalpha *= *(\d+)/
        nalpha = Integer($1)
      elsif line =~ /^ *Nbeta *= *(\d+)/
        nbeta = Integer($1)
      end
    end
    hide_progress_panel
    clear_basis_set
    clear_mo_coefficients
    set_mo_info(:type => hf_type, :alpha => nalpha, :beta => nbeta)
    return true
  end

  #  mol.set_mo_info should be set before calling this function
  #  Optional label is for importing JANPA output: "NAO" or "CPLO"
  #  If label is not nil, then returns a hash containing the following key/value pairs:
  #    :atoms => an array of [element_symbol, seq_num, atomic_num, x, y, z] (angstrom)
  #    :gto => an array of an array of [sym, [ex0, c0, ex1, c1, ...]]
  #    :moinfo => an array of [sym, energy, spin (0 or 1), occ]
  #    :mo => an array of [c0, c1, ...]
  def sub_load_molden(fp, label = nil)
    getline = lambda { @lineno += 1; return fp.gets }
    bohr = 0.529177210903
    errmsg = nil
    ncomps = 0  #  Number of components (AOs)
    occ_alpha = 0  #  Number of occupied alpha orbitals
    occ_beta = 0   #  Number of occupied beta orbitals
    if label
      hash = Hash.new
    end
    #  The GTOs (orbital type, contractions and exponents) are stored in gtos[]
    #  and set just before first [MO] is processed.
    #  This is because we do not know whether the orbital type is cartesian or spherical
    #  until we see lines like "[5D]".
    gtos = []
    spherical_d = false
    spherical_f = false
    spherical_g = false
    #  Number of components for each orbital type
    ncomp_hash = { 0=>1, 1=>3, -1=>4, 2=>6, -2=>5, 3=>10, -3=>7, 4=>15, -4=>9 }
    catch :ignore do
      while line = getline.call
        if line =~ /^\[Atoms\]/
          i = 0
          while line = getline.call
            if line =~ /^[A-Z]/
              #  element, index, atomic_number, x, y, z (in AU)
              a = line.split(' ')
              if label
                (hash[:atoms] ||= []).push([a[0], Integer(a[1]), Integer(a[2]), Float(a[3]) * bohr, Float(a[4]) * bohr, Float(a[5]) * bohr])
              else
                if atoms[i].atomic_number != Integer(a[2]) ||
                  (atoms[i].x - Float(a[3]) * bohr).abs > 1e-4 ||
                  (atoms[i].y - Float(a[4]) * bohr).abs > 1e-4 ||
                  (atoms[i].z - Float(a[5]) * bohr).abs > 1e-4
                  errmsg = "The atom list does not match the current molecule."
                  throw :ignore
                end
              end
              i += 1
            else
              break
            end
          end
          redo  #  The next line will be the beginning of the next block
        elsif line =~ /^\[GTO\]/
          shell = 0
          atom_index = 0
          while line = getline.call
            #  index, 0?
            a = line.split(' ')
            break if a.length != 2
            atom_gtos = []  #  [[sym1, [e11, c11, e12, c12, ...], add_exp1], [sym2, [e21, c22, ...], add_exp2], ...]
            #  loop for shells
            while line = getline.call
              #  type, no_of_primitives, 1.00?
              a = line.split(' ')
              break if a.length != 3   #  Terminated by a blank line
              a[0] =~ /^([a-z]+)([0-9]+)?$/
              symcode = $1
              add_exp = ($2 == nil ? 0 : $2.to_i)
              case symcode
              when "s"
                sym = 0
              when "p"
                sym = 1
              when "d"
                sym = 2
              when "f"
                sym = 3
              when "g"
                sym = 4
              else
                raise MolbyError, "Unknown gaussian shell type '#{a[0]}' at line #{@lineno} in MOLDEN file"
              end
              nprimitives = Integer(a[1])
              gtoline = [sym, [], add_exp]
              atom_gtos.push(gtoline)
              nprimitives.times { |i|
                line = getline.call   #  exponent, contraction
                b = line.split(' ')
                gtoline[1].push(Float(b[0]), Float(b[1]))
              }
              #  end of one shell
              shell += 1
            end
            #  end of one atom
            atom_index += 1
            gtos.push(atom_gtos)
          end
          if label
            hash[:gto] = gtos
          end
          redo  #  The next line will be the beginning of the next block
        elsif line =~ /^\[5D\]/ || line =~ /^\[5D7F\]/
          spherical_d = spherical_f = true
        elsif line =~ /^\[5D10F\]/
          spherical_d = true
          spherical_f = false
        elsif line =~ /^\[7F\]/
          spherical_f = true
        elsif line =~ /^\[9G\]/
          spherical_g = true
        elsif line =~ /^\[MO\]/
          #  Add shell info and primitive coefficients to molecule
          gtos.each_with_index { | atom_gtos, atom_index|
            atom_gtos.each { |gtoline|
              sym = gtoline[0]
              #  Change orbital type if we use spherical functions
              sym = -2 if sym == 2 && spherical_d
              sym = -3 if sym == 3 && spherical_f
              sym = -4 if sym == 4 && spherical_g
              gtoline[0] = sym
              coeffs = gtoline[1]
              nprimitives = coeffs.length / 2
              add_exp = gtoline[2]
              ncomps += ncomp_hash[sym]
              if !label
                add_gaussian_orbital_shell(atom_index, sym, nprimitives, add_exp)
                nprimitives.times { |prim|
                  add_gaussian_primitive_coefficients(coeffs[prim * 2], coeffs[prim * 2 + 1], 0.0)
                }
              end
            }
          }
          m = []
          idx_alpha = 1   #  set_mo_coefficients() accepts 1-based index of MO
          idx_beta = 1
          if label
            hash[:mo] = []
            hash[:moinfo] = []
          end
          while true
            #  Loop for each MO
            m.clear
            ene = nil
            spin = nil
            sym = nil   #  Not used in Molby
            occ = nil
            i = 0
            while line = getline.call
              if line =~ /^ *Sym= *(\w+)/
                sym = $1
              elsif line =~ /^ *Ene= *([-+.0-9eE]+)/
                ene = Float($1)
              elsif line =~ /^ *Spin= *(\w+)/
                spin = $1
              elsif line =~ /^ *Occup= *([-+.0-9eE]+)/
                occ = Float($1)
                if occ > 0.0
                  if spin == "Alpha"
                    occ_alpha += 1
                  else
                    occ_beta += 1
                  end
                end
                if label
                  hash[:moinfo].push([sym, ene, (spin == "Alpha" ? 0 : 1), occ])
                end
              elsif line =~ /^ *([0-9]+) +([-+.0-9eE]+)/
                m[i] = Float($2)
                i += 1
                if i >= ncomps
                  if spin == "Alpha"
                    idx = idx_alpha
                    idx_alpha += 1
                  else
                    idx = idx_beta
                    idx_beta += 1
                  end
                  if label
                    hash[:mo].push(m.dup)
                  else
                    set_mo_coefficients(idx, ene, m)
                  end
                  break
                end
              else
                break
              end
            end
            break if i < ncomps  #  no MO info was found
          end
          next
        end #  end if
      end   #  end while
    end     #  end catch
    if errmsg
      message_box("The MOLDEN file was found but not imported. " + errmsg, "Psi4 import info", :ok)
      return (label ? nil : false)
    end
    return (label ? hash : true)
  end

  #  Import the JANPA log and related molden files
  #  Files: inppath.{NAO.molden,CLPO.molden,janpa.log}
  #  If inppath.spherical.molden is available, then clear existing mo info
  #  and load from it (i.e. use the basis set converted by molden2molden)
  def sub_load_janpa_log(inppath)
    begin
      m2name_p = {0=>"x", 1=>"y", -1=>"z"}
      m2name_d = {0=>"zz-rr", 1=>"xz", -1=>"yz", 2=>"xx-yy", -2=>"xy"}
      m2name_f = {0=>"z3-zr2", 1=>"xz2-xr2", -1=>"yz2-yr2", 2=>"x2z-y2z", -2=>"xyz", 3=>"x3-xy2", -3=>"x2y-y3"}
      m2name_g = {0=>"z4-z2r2+r4", 1=>"xz3-xzr2", -1=>"yz3-yzr2", 2=>"x2z2-y2z2", -2=>"xyz2-xyr2", 3=>"x3z-xy2z", -3=>"x2yz-y3z", 4=>"x4-x2y2+y4", -4=>"x3y-xy3"}
      fp = File.open(inppath + ".janpa.log", "rt") rescue fp = nil
      if fp == nil
        hide_progress_panel  #  Close if it is open
        message_box("Cannot open JANPA log file #{inppath + '.janpa.log'}: " + $!.to_s)
        return false
      end
      print("Importing #{inppath}.janpa.log.\n")
      lineno = 0
      getline = lambda { lineno += 1; return fp.gets }
      h = Hash.new
      mfiles = Hash.new
      h["software"] = "JANPA"
      nao_num = nil  #  Set later
      nao_infos = [] #  index=atom_index, value=Hash with key "s", "px", "py" etc.
      #  nao_infos[index][key]: array of [nao_num, occupancy], in the reverse order of appearance
      while line = getline.call
        if line =~ /molden2molden: a conversion tool for MOLDEN/
          while line = getline.call
            break if line =~ /^All done!/
            if line =~ /\.spherical\.molden/
              #  The MOs are converted to spherical basis set
              #  Clear the existing MO and load *.spherical.molden
              sname = inppath + ".spherical.molden"
              fps = File.open(sname, "rt") rescue fps = nil
              if fps != nil
                print("Importing #{sname}.\n")
                @lineno = 0
                type = get_mo_info(:type)
                alpha = get_mo_info(:alpha)
                beta = get_mo_info(:beta)
                clear_basis_set
                set_mo_info(:type=>type, :alpha=>alpha, :beta=>beta)
                #  mol.@hf_type should be set before calling sub_load_molden
                @hf_type = type
                sub_load_molden(fps)
                fps.close
              end
            end
          end
        elsif line =~ /^NAO \#/
          h["NAO"] = []  #  [num, anum, occupied_p, group_num, orb_desc, occ, principle*10+orb, shelltype]
          # num: NAO index (0-based), anum: atom index (0-based),
          # occupied_p: true if marked occupied by JANPA
          # group_num: NAO group number (same for same atom, same principle, same orbital type)
          # orb_desc: like "s", "px|py|pz", "dxy|dyz|dzx|dx2-y2|dz2-r2", etc.
          # occ: occupancy
          # principle: principle number (calculated later from group_num)
          # orb: 0:s, 1:p, 2:d, 3:f, 4:g
          # shelltype: 0: core, 1: valence, 2: rydberg
          nao_infos = []  #  nao_infos[atom_index] = {orb_desc=>[[nao0, occ0], [nao1, occ1] ...]}
          # orb_desc: same as above
          # nao, occ: nao index (0-based) and occupancy
          while line = getline.call
            break if line !~ /^\s*[1-9]/
            num = Integer(line[0, 5]) - 1
            name = line[5, 21]
            occ = Float(line[26, 11])
            #  like A1*: R1*s(0)
            #  atom_number, occupied?, group_number, orb_sym, angular_number
            name =~ /\s*[A-Z]+([0-9]+)(\*?):\s* R([0-9]+)\*([a-z]+)\(([-0-9]+)\)/
            anum = Integer($1) - 1
            occupied = $2
            group_num = Integer($3)
            orb_sym = $4
            ang_num = Integer($5)
            orb_desc = orb_sym
            if orb_desc == "p"
              orb_desc += m2name_p[ang_num]
            elsif orb_desc == "d"
              orb_desc += m2name_d[ang_num]
            elsif orb_desc == "f"
              orb_desc += m2name_f[ang_num]
            elsif orb_desc == "g"
              orb_desc += m2name_g[ang_num]
            end
            h["NAO"].push([num, anum, occupied, group_num, orb_desc, occ, nil, nil])
            nao_num = h["NAO"].length
            ((nao_infos[anum] ||= Hash.new)[orb_desc] ||= []).unshift([num, occ])
          end
          nao_num = h["NAO"].length
          #  Set principle from nao_infos
          val_occ = []  #  val_occ[atom_index][principle*10+orb] = (0: core, 1: valence, 2: ryd)
          nao_infos.each_with_index { |value, atom_index|
            val_occ[atom_index] ||= []
            value.each { |orb_desc, ar|
              ar.each_with_index { |v, idx|
                if v[1] > 1.9
                  val = 0  #  lp
                elsif v[1] > 0.01
                  val = 1  #  valence
                else
                  val = 2  #  ryd
                end
                principle = idx + 1
                orb_sym = orb_desc[0]
                orb = 0
                if orb_sym == "p"
                  orb = 1
                elsif orb_sym == "d"
                  orb = 2
                elsif orb_sym == "f"
                  orb = 3
                elsif orb_sym == "g"
                  orb = 4
                end
                principle += orb
                po = principle * 10 + orb
                val1 = val_occ[atom_index][po]
                if val1 == nil
                  val_occ[atom_index][po] = val
                elsif val1 != val
                  val_occ[atom_index][po] = 1  #  core & val, core & ryd or val & ryd
                else
                  #  If all orbitals in the shell are "lp", then the shell should be core.
                  #  If all orbitals in the shell are "ryd", then the shell should be ryd.
                  #  Otherwise, the shell is valence.
                end
                h["NAO"][v[0]][6] = po
              }
            }
          }
          nao_num.times { |nao_index|
            nao = h["NAO"][nao_index]
            nao[7] = val_occ[nao[1]][nao[6]]  #  shell type
          }
        elsif line =~ /^\s*(C?)LPO\s+D e s c r i p t i o n\s+Occupancy\s+Composition/
          if $1 == "C"
            key = "CLPO"
          else
            key = "LPO"
          end
          h[key] = []  #  [[a1, a2], LP|BD|NB, desc, occ, [h1, h2], lpo_idx, nao_idx]
                       # a1, a2: atom indices (0-based), h1, h2: hybrid indices (0-based)
                       # lpo_idx: orbital index (0-based)
                       # nao_idx: most significant nao (set later)
          if key == "CLPO"
            h["LHO"] = []  # [[a1, a2], lho_idx, nao_idx]
                           # a1, a2: atom indices (0-based)
                           # nao_idx: most significant nao (set later)
          end
          while line = getline.call
            break if line =~ /^\s*$/
            num = Integer(line[0, 5]) - 1
            label1 = line[5, 6].strip
            desc = line[11, 30].strip
            occ = line[41, 11].strip
            comp = line[52, 1000].strip
            desc =~ /\s*([-A-Za-z0-9]+)(,\s*(.*$))?/
            desc1 = $1
            desc2 = ($3 || "")
            if desc2 =~ /^(.*)*\(NB\)\s*$/ && label1 == ""
              label1 = "(NB)"
              desc2 = $1.strip
            end
            label1.gsub!(/[\(\)]/, "")  #  Strip parens
            atoms = desc1.scan(/[A-Za-z]+(\d+)/)   # "C1-H3" -> [["1"], ["3"]]
            atoms = atoms.map { |a| Integer(a[0]) - 1 }  # [0, 2]
            hybrids_a = comp.scan(/h(\d+)@[A-Za-z]+(\d+)/)  #  "h8@C1...h13@H3" -> "[["8", "1"], ["13", "3"]]
            hybrids = []
            hybrids_a.each { |a|
              i = atoms.find_index(Integer(a[1]) - 1)
              if i != nil
                hybrids[i] = Integer(a[0]) - 1
              end
            } # [7, 12]
            #  like ["BD", [0, 2], "Io = 0.2237", occ, [7, 12]]
            #  0, 2 are the atom indices (0-based)
            #  7, 12 are the number of hybrid orbitals (0-based)
            h[key][num] = [atoms, label1, desc2, Float(occ), hybrids, num]
            if key == "CLPO"
              h["LHO"][hybrids[0]] = [atoms, hybrids[0]]
              if hybrids[1]
                h["LHO"][hybrids[1]] = [atoms.reverse, hybrids[1]]
              end
            end
          end
          nao_num.times { |idx|
            if h[key][idx] == nil
              h[key][idx] = [nil, nil, nil, nil, nil, idx]
            end
          }
        elsif line =~ /^ -NAO_Molden_File: (\S*)/
          mfiles["NAO"] = $1
        elsif line =~ /^ -LHO_Molden_File: (\S*)/
          mfiles["LHO"] = $1
        elsif line =~ /^ -CLPO_Molden_File: (\S*)/
          mfiles["CLPO"] = $1
        elsif line =~ /^ -PNAO_Molden_File: (\S*)/
          mfiles["PNAO"] = $1
        elsif line =~ /^ -AHO_Molden_File: (\S*)/
          mfiles["AHO"] = $1
        elsif line =~ /^ -LPO_Molden_File: (\S*)/
          mfiles["LPO"] = $1
        end
      end
      fp.close

      #  Read molden files
      keys = []
      mfiles.each { |key, value|
        fp = Kernel.open(value, "rt") rescue fp = nil
        if fp
          print("Importing #{value}.\n")
          res = sub_load_molden(fp, key)
          if res
            #  Temporarily assign: this array will be reordered later and converted to LAMatrix
            h["AO/#{key}"] = res[:mo]
            #h["AO/#{key}"] = LAMatrix.new(res[:mo])
          end
          fp.close
          keys.push(key)
        end
      }

      #  Reorder orbitals
      natoms = self.natoms
      keys.each { |key|
        if key == "NAO"
          # [num, anum, occupied_p, group_num, orb_desc, occ, principle*10+orb, shelltype]
          # reorder by (1) anum and shelltype (core/val first, ryd later),
          # (2) principle*10+orb, (3) num
          reorder = []
          h["NAO"].each { |nao|
            if nao[7] == 2
              r = (nao[1] + natoms) * 3 + 2
            else
              r = nao[1] * 3 + nao[7]
            end
            r = r * 100 + nao[6]
            r = r * nao_num + nao[0]
            reorder.push(r)
          }
          h["NAOnew2old"] = (0...nao_num).to_a.sort { |a, b| reorder[a] <=> reorder[b] }
          h["NAOold2new"] = []
          h["NAOnew2old"].each_with_index { |i, j|
            h["NAOold2new"][i] = j
          }
        else
          #  LPO, CLPO: [[a1, a2], LP|BD|NB, desc, occ, [h1, h2], lpo_idx, nao_idx, shelltype]
          #  LHO: [[a1, a2], lho_idx, nao_idx, shelltype]
          if h[key]
            #  Set the most significant nao
            mo = LAMatrix.new(h["AO/NAO"]).inverse * LAMatrix.new(h["AO/#{key}"])
            pos = (key == "LHO" ? 2 : 6)  #  position of nao_idx (new index)
            nao_num.times { |o_idx|
              o = (h[key][o_idx] ||= [])
              if o[pos - 1] == nil
                o[pos - 1] = o_idx  #  Original index
              end
              #  Find the NAO that has the largest contribution to this hybrid orbital
              a_max = -1
              i_max = nil
              nao_num.times { |i|
                a = mo[o_idx, i] ** 2
                if a > a_max
                  i_max = i
                  a_max = a
                end
              }
              o[pos] = h["NAOold2new"][i_max]  #  new index
              o[pos + 1] = h["NAO"][i_max][7]  #  0:core, 1:val, 2:ryd
              if o[0] == nil
                #  Set the atom which has the most significant NAO
                o[0] = [h["NAO"][i_max][1]]
              end
            }
            reorder = []
            h[key].each_with_index { |o, o_idx|
              #  Reorder by (1) atom and shelltype (core first, val second, ryd later),
              #  if LP|BD|NB is present, then val/LP and val/BD come earlier than val/NB
              #  (2) second atom if present, and (3) most significant NAO index
              if o[pos + 1] == 2
                r = (o[0][0] + natoms) * 4 + 3
              else
                r = o[0][0] * 4 + o[pos + 1]
                if o[pos + 1] == 1 && (key == "LPO" || key == "CLPO") && o[1] == "NB"
                  r = r + 1
                end
              end
              r = r * (natoms + 1) + (o[0][1] || -1) + 1  #  Second atom index + 1 (0 if none)
              r = r * nao_num + o[pos]
              reorder.push(r)
            }
            h[key + "new2old"] = (0...nao_num).to_a.sort { |a, b| reorder[a] <=> reorder[b] }
            h[key + "old2new"] = []
            h[key + "new2old"].each_with_index { |i, j|
              h[key + "old2new"][i] = j
            }
          end
        end
      }
      #print("h[\"NAO\"] = " + h["NAO"].inspect + "\n")
      #print("h[\"LHO\"] = " + h["LHO"].inspect + "\n")
      #print("h[\"LPO\"] = " + h["LPO"].inspect + "\n")
      #print("h[\"CLPO\"] = " + h["CLPO"].inspect + "\n")

      #  Do reorder and make matrices
      keys.each { |key|
        if key == "NAO" || key == "PNAO"
          old2new = h["NAOold2new"]
        elsif key == "LPO"
          old2new = h["LPOold2new"]
        elsif key == "CLPO"
          old2new = h["CLPOold2new"]
        elsif key == "AHO" || key == "LHO"
          old2new = h["LHOold2new"]
        end
        if old2new
          mo = []
          nao_num.times { |i|
            mo[old2new[i]] = h["AO/" + key][i]
          }
          h["AO/" + key] = mo
        end
        h["AO/" + key] = LAMatrix.new(h["AO/" + key])
      }
      
      #  Create labels
      keys.each { |key|
        old2new = h[key + "old2new"]
        if key == "NAO"
          h["NAO_L"] = []
          nao_num.times { |nao_idx|
            # [num, anum, occupied_p, group_num, orb_desc, occ, principle*10+orb, shelltype]
            nao = h["NAO"][nao_idx]
            aname = self.atoms[nao[1]].name
            orb_desc = nao[4]
            principle = nao[6] / 10
            shell = ["core", "val", "ryd"][nao[7]]
            j = (old2new ? old2new[nao_idx] : nao_idx)
            h["NAO_L"][j] = "#{aname} (#{principle}#{orb_desc}) (#{shell}, \##{nao_idx + 1})"
          }
        elsif key == "LHO"
          h["LHO_L"] = []
          nao_num.times { |lho_idx|
            # [[a1, a2], lho_idx, nao_idx, shelltype]
            lho = h["LHO"][lho_idx]
            aname1 = self.atoms[lho[0][0]].name
            aname2 = ("(" + self.atoms[lho[0][1]].name + ")") rescue ""
            shell = ["core", "val", "ryd"][lho[3]]
            j = (old2new ? old2new[lho_idx] : lho_idx)
            h["LHO_L"][j] = "#{aname1}#{aname2} (#{shell}, \##{lho_idx + 1})"
          }
        elsif key == "LPO" || key == "CLPO"
          h[key + "_L"] = []
          nao_num.times { |key_idx|
            # [[a1, a2], LP|BD|NB, desc, occ, [h1, h2], lpo_idx, nao_idx, shelltype]
            o = h[key][key_idx]
            aname1 = self.atoms[o[0][0]].name
            aname2 = ("-" + self.atoms[o[0][1]].name) rescue ""
            type = o[1]
            shelltype = o[7]
            if shelltype == 2
              type = "RY"  #  Rydberg
            elsif type == "LP" && shelltype == 0
              type = "CR"  #  core
            end
            j = (old2new ? old2new[key_idx] : key_idx)
            h[key + "_L"][j] = "#{aname1}#{aname2} (#{type}, \##{key_idx + 1})"
          }
        end
      }
      
      @nbo = h
      if @nbo["AO/NAO"] && @nbo["AO/LHO"] && @nbo["AO/PNAO"]
        #  Generate PLHO from PNAO, NAO, LHO
        #  This protocol was suggested by the JANPA author in a private commnunication.
        begin
          nao2lho = @nbo["AO/NAO"].inverse * @nbo["AO/LHO"]
          nao2pnao = @nbo["AO/NAO"].inverse * @nbo["AO/PNAO"]
          sign = LAMatrix.diagonal((0...nao2pnao.column_size).map { |i| (nao2pnao[i, i] < 0 ? -1 : 1)})
          @nbo["AO/PLHO"] = @nbo["AO/PNAO"] * sign * nao2lho
        rescue
          @nbo["AO/PLHO"] = nil
        end
      end
      if @nbo["AO/AHO"] && @nbo["LHO_L"]
        #  Copy labels from LHO
        @nbo["AHO_L"] = @nbo["LHO_L"].dup
      end
      return true
    rescue => e
      $stderr.write(e.message + "\n")
      $stderr.write(e.backtrace.inspect + "\n")
    end
  end

  def loadout(filename)
  retval = false
  fp = open(filename, "rb")
  @lineno = 0
  begin
    while s = fp.gets
      @lineno += 1
      if s =~ /Gaussian/
        retval = sub_load_gaussian_log(fp)
        break
      elsif s =~ /GAMESS/
        retval = sub_load_gamess_log(fp)
        break
      elsif s =~ /Psi4/
        retval = sub_load_psi4_log(fp)
        if retval
          #  If .molden file exists, then try to read it
          namepath = filename.gsub(/\.\w*$/, "")
          mname = "#{namepath}.molden"
          if File.exists?(mname)
            fp2 = open(mname, "rb")
            if fp2
              flag = sub_load_molden(fp2)
              fp2.close
              status = (flag ? 0 : -1)
            end
          end
          if File.exists?("#{namepath}.janpa.log")
            flag = sub_load_janpa_log(namepath)
            status = (flag ? 0 : -1)
          end
        end
        break
      end
    end
  rescue
    hide_progress_panel
    raise
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
