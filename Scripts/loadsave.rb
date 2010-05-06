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
	show_progress_panel("Loading a crd file...")
#    puts "sframe = #{sframe}, pos = #{fp.pos}"
    while 1
      line = fp.gets
      if line == nil
	    if (count > 0)
		  raise MolbyError, sprintf("crd format error - file ended in the middle of frame %d", frame)
		end
	    fp.close
		break
      end
      next if line.match(/^TITLE/)
      line.chomp
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
        count = 0
        frame += 1
		if frame % 5 == 0
		  set_progress_message("Loading a crd file...\n(#{frame} frames completed)")
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
	show_progress_panel("Saving a crd file...")
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
			set_progress_message("Saving a crd file...\n(#{i} frames completed)")
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

  def loadlog(filename)

    if natoms == 0
		new_unit = true
	else
		new_unit = false
    end
#	save_undo_enabled = self.undo_enabled?
#	self.undo_enabled = false
	self.update_enabled = false
	show_progress_panel("Loading GAMESS log file...")
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
			line.chomp
			if line =~ /ATOM\s+ATOMIC\s+COORDINATES/ || line =~ /COORDINATES OF ALL ATOMS ARE/
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
			end
		end
	ensure
#		self.undo_enabled = save_undo_enabled
        hide_progress_panel
		self.update_enabled = true
	end
	if nframes > 0
	  select_frame(nframes - 1)
	end
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
  end

  alias :loadgjf :loadcom
  alias :savegjf :savecom
  
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
