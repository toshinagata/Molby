#
#  molecule.rb
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

  def min(group = nil)
    rv = Vector3D[1e20, 1e20, 1e20]
    group = atom_group(group ? group : 0...natoms)
    group.each { |i|
      v = atoms[i].r
      rv.x = v.x if v.x < rv.x
      rv.y = v.y if v.y < rv.y
      rv.z = v.z if v.z < rv.z
    }
    rv
  end
  
  def max(group = nil)
    rv = Vector3D[-1e20, -1e20, -1e20]
    group = atom_group(group ? group : 0...natoms)
    group.each { |i|
      v = atoms[i].r
      rv.x = v.x if v.x > rv.x
      rv.y = v.y if v.y > rv.y
      rv.z = v.z if v.z > rv.z
    }
    rv
  end

  def neutralize(charge = 0.0, group = nil)
    group = atom_group(group ? group : 0...natoms)
    csum = 0.0
    group.each { |i| csum += atoms[i].charge }
    if csum != 0.0
      csum = -csum / group.length
      group.each { |i| atoms[i].charge += csum }
    end
  end

  def all
    return atom_group
  end
  
  def rotate_with_axis(xv, yv, cv = nil, group = nil)
    cv = Vector3D[0, 0, 0] if !cv
    v1 = (xv - cv).normalize
    v3 = v1.cross(yv - cv).normalize
    v2 = v3.cross(v1).normalize
#    tr = Transform[v1, v2, v3, cv - Transform[v1, v2, v3, [0, 0, 0]] * cv]
    tr = Transform[v1, v2, v3, [0,0,0]].inverse * Transform.translation(cv * (-1))
    transform(tr, group)
  end

  #  Get the fragment including the atom n1. If additional atoms are given,
  #  those atoms will not be counted during the search.
#  def fragment(n1, *n2)
#    def fragment_sub(idx, g, gx)
#	  return if gx.include?(idx)
#      g << idx
#      for i in atoms[idx].connects
#        next if g.include?(i)
#        fragment_sub(i, g, gx)
#      end
#    end
#	n1 = atom_index(n1)
#	g = atom_group(n1)
#	gx = atom_group(*n2)
#    for i in atoms[n1].connects
#      fragment_sub(i, g, gx)
#    end
#    g
#  end

  #  Rotate the molecular fragment. The fragment given by "fragment(n2, n1)" is rotated,
  #  with the atom n1 as the center and n2->n1 as the axis.
  def rotate_fragment(n1, n2, angle)
    frag = fragment(n2, n1)
    cen = atoms[n1].r
    axis = (atoms[n2].r - cen).normalize
    rotate(axis, angle, cen, frag)
  end
  
  #  Calculate the bond length defined by two vectors.
  def Molecule.calc_bond(v1, v2)
    (v1 - v2).length
  end
  
  #  Calculate the angle defined by three vectors.
  def Molecule.calc_angle(v1, v2, v3)
    r21 = v1 - v2
    r23 = v3 - v2
    w1 = r21.length
    w2 = r23.length
    if (w1 < 1e-5 || w2 < 1e-5)
      return 0.0
    end
    cs = r21.dot(r23) / (w1 * w2)
    Math.atan2(Math.sqrt(1 - cs*cs), cs) * Rad2Deg
  end
  
  #  Calculate the dihedral angle defined by four vectors.
  def Molecule.calc_dihedral(v1, v2, v3, v4)
    r21 = v1 - v2
    r32 = v2 - v3
    r43 = v3 - v4
    cr1 = r21.cross(r32)
    cr2 = r32.cross(r43)
    cr3 = r32.cross(cr1)
    w1 = cr1.length
    w2 = cr2.length
    w3 = cr3.length
    if (w1 < 1e-5 || w2 < 1e-5 || w3 < 1e-5)
      return 9999.0
    end
    sn = cr3.dot(cr2) / (w3 * w2)
    cs = cr1.dot(cr2) / (w1 * w2)
    Math.atan2(-sn, cs) * Rad2Deg
  end
  
  #  Calculate the bond length between the two atoms.
  def calc_bond(n1, n2)
    n1 = atoms[n1].r if !n1.is_a?(Vector3D)
    n2 = atoms[n2].r if !n2.is_a?(Vector3D)
    self.class.calc_bond(n1, n2)
  end

  #  Calculate the bond angle defined by the three atoms.
  def calc_angle(n1, n2, n3)
    n1 = atoms[n1].r if !n1.is_a?(Vector3D)
    n2 = atoms[n2].r if !n2.is_a?(Vector3D)
    n3 = atoms[n3].r if !n3.is_a?(Vector3D)
    self.class.calc_angle(n1, n2, n3)
  end

  #  Calculate the dihedral angle defined by the four atoms
  def calc_dihedral(n1, n2, n3, n4)
    n1 = atoms[n1].r if !n1.is_a?(Vector3D)
    n2 = atoms[n2].r if !n2.is_a?(Vector3D)
    n3 = atoms[n3].r if !n3.is_a?(Vector3D)
    n4 = atoms[n4].r if !n4.is_a?(Vector3D)
    self.class.calc_dihedral(n1, n2, n3, n4)
  end

  #  Set the dihedral angle n1-n2-n3-n4 to the specified value. The fragment
  #  including n3 and n4 are rotated.
  def set_dihedral(n1, n2, n3, n4, angle)
    dihed = calc_dihedral(n1, n2, n3, n4)
    if (dihed < 9999.0)
      rotate_fragment(n2, n3, dihed - angle)
    end
  end

  #  Intended for use internally. Find the bond vector to be created. If atom to remove (rem) is specified, it is
  #  the vector base->rem. If rem is not specified, look for the "dummy" atom if present.
  #  Otherwise, no atom is removed, and the new bond vector is defined as the most distant
  #  direction from any existing bonds to base. Returns [n, v] where n is the atom to be
  #  removed (or nil), and v is the bond vector.

  def find_dock_vec(base, rem)
    if (!rem)
      atoms.each { |ap|
        if ap.name[0] == ?#
          rem = ap.index
          break
        end
      }
    end
    if (!rem)
      v1 = atoms[base].r
      v2 = Vector3D[0, 0, 0]
      for i in atoms[base].connects
        v2 += (atoms[i].r - v1).normalize
      end
      v2 = (v2.length > 0 ? v2.normalize * (-1.0) : Vector3D[1, 0, 0])
    else
      v2 = (atoms[rem].r - atoms[base].r).normalize
    end
    [rem, v2]
  end

  #  Combine two molecules. A new bond is created between base1 and base2. rem1 and rem2
  #  specify the atom (or fragment) to be removed when creating the new bond. If nil is
  #  specified for either of them, a dummy atom (that has a '#' at the top of its name)
  #  is looked for, and if none is found then no atom is removed and the new bond is
  #  directed so that it is most remote from any of the existing bonds.
  #  rem1 and rem2 can also be Vector3D's that specify the direction of the new bond.
  #  In this case, no atoms are removed.
  #  If len is specified, then the length of the new bond is set to this value. Otherwise,
  #  it is set to 1.5.
  #  If dihed is specified, then the dihedral angle def1-base1-base2-def2 is set to this
  #  value (in degree). If nil is specified for either def1 or def2 (or both), then
  #  the atom which is connected to base1 (base2) and has the smallest index is used.
  #  Returns the group of atom indices that have been added.
  
  def dock(mol, base1, base2, rem1 = nil, rem2 = nil, len = nil, dihed = nil, def1 = nil, def2 = nil)

	if rem1.is_a?(Vector3D)
	  v1 = rem1.normalize
	  rem1 = nil
	else
      rem1, v1 = find_dock_vec(base1, rem1)
	end
	if rem2.is_a?(Vector3D)
	  v2 = rem2.normalize
	  rem2 = nil
	else
      rem2, v2 = mol.find_dock_vec(base2, rem2)
	end
    if (!len)
      #  If new bond length is not specified:
      if (rem1)
        #  If the atom to remove is specified, then the old bond length is kept
        len = (atoms[rem1].r - atoms[base1].r).length
      else
        #  Otherwise, an arbitrary bond length is used
        len = 1.5
      end
    end
    #  Duplicate the second molecule
    mol = mol.dup
    #  Rotate the second molecule so that v2 becomes antiparallel to v1
    cs = -(v1.dot(v2))
    if (cs < 1 - 1e-6)
      if (cs < -1 + 1e-6)
        #  v1 == v2 (180 deg rotation)
        v3 = v2.cross([0, 1, 0])
        if (v3.length < 1e-6)
          v3 = v2.cross([0, 0, 1])
        end
        v3 = v3.normalize
      else
        v3 = v2.cross(v1).normalize
      end
      angle = Math.atan2(Math.sqrt(1.0 - cs*cs), cs) * Rad2Deg
      mol.rotate(v3, angle, mol.atoms[base2].r)
    end
    #  Move the second molecule so that the atom 'base2' is located at atoms[base1].r+v1*len
    mol.translate(atoms[base1].r + v1 * len - mol.atoms[base2].r)
    if (dihed)
      #  Rotate the bond base1-base2 so that the dihedral angle becomes dihed
      #  Find atom of lowest index if def1/2 is not specified
      def1 = (atoms[base1].connects - (rem1 ? [atom_index(rem1)] : [])).min if !def1
      def2 = (mol.atoms[base2].connects - (rem2 ? [mol.atom_index(rem2)] : [])).min if !def2
	  if def1 && def2
        dihed1 = self.calc_dihedral(atoms[def1].r, atoms[base1].r, mol.atoms[base2].r, mol.atoms[def2].r)
        mol.rotate(v1, dihed1 - dihed, mol.atoms[base2].r) if dihed1 < 9999.0
	  end
    end
    #  Calculate the atom indices for the combined molecule
	natoms1 = natoms
	natoms2 = mol.natoms
    base1 = atom_index(base1)
    rem1 = atom_index(rem1) if rem1
    base2 = mol.atom_index(base2) + natoms1
    rem2 = mol.atom_index(rem2) + natoms1 if rem2
	#  Fragment to remove from self
    g = IntGroup[]
	if rem1
	  f1 = fragment(rem1, base1)
	  natoms1 -= f1.length
	  g += f1
	  assign_residue(f1, 0)
	end
	ofs = self.max_residue_number
	self.nresidues = ofs + 1
#	mol.offset_residue(mol.all, ofs)
    #  Merge two molecules, create a bond, and remove fragments if necessary
    add(mol)
    create_bond(base1, base2)
	#  Fragment to remove from added molecule
	if rem2
	  f2 = fragment(rem2, base2)
	  natoms2 -= f2.length
	  g += f2
	end
    remove(g)
#	mol.offset_residue(mol.all, -ofs)
	#  Returns the group of appended atom indices
	IntGroup[natoms1...(natoms1+natoms2)]
  end

  #  Add a new atom. (bond, base1, angle, base2, dihed, base3) defines
  #  the position of the new atom in the Z-matrix style. (Note: angle 
  #  and dihed should be given in *degree*, not in radian). If bond/base1
  #  are specified, a new bond is also created between the new atom and
  #  base1. If bond/base1 are specified but angle/base2 are not, then
  #  the direction of the new bond is assumed as the most distant
  #  direction from the existing bonds.
  #  Returns the reference to the new atom.
  def add_atom(name, atom_type = "c3", element = "C", bond = nil, base1 = nil, angle = nil, base2 = nil, dihed = nil, base3 = nil)
    ap = create_atom(name)
    ap.atom_type = atom_type
    ap.element = element
    if (base1 != nil)
      bp1 = atoms[base1].r
      if ap.res_seq == 0 && (res_seq = atoms[base1].res_seq) > 0
        group = atom_group { |p| p.res_seq == res_seq }
        group << ap.index
        assign_residue(group, atoms[base1].res_name + ".#{res_seq}")
      end
      if base2 == nil
        n, r = find_dock_vec(base1, nil)
        r = bp1 + r.normalize * bond
      else
        v2 = atoms[base2].r - bp1
        if base3 != nil
          v3 = atoms[base3].r - bp1
        else
          v3 = Vector3D[0, 0, 1]
          dihed = 180
        end
        angle *= Deg2Rad
        dihed *= Deg2Rad
        vx = v2.normalize
        vy = v3.cross(v2)
        if vy.length < 1e-8
          n = Math::floor(angle / Math::PI + 0.5)
          if ((angle - n * Math::PI).abs < 1e-8)
            vd = vx * (bond * Math::cos(angle))
          else
            raise "Cannot define dihedral angle for atom #{name}-#{base1}-#{base2}-#{base3}"
          end
        else
          vy = vy.normalize
          vz = vx.cross(vy).normalize
          x = bond * Math::cos(angle)
          y = bond * Math::sin(angle) * Math::sin(dihed)
          z = bond * Math::sin(angle) * Math::cos(dihed)
          vd = Transform[vx, vy, vz, [0, 0, 0]] * Vector3D[x, y, z]
        end
        r = bp1 + vd
      end
      ap.r = r
      create_bond(base1, ap.index)
    end
    return ap
  end

  #  Clean up atom names. All atoms in the group (if specified) are named
  #  with the sequential numbers.
  def guess_names(group = nil)
    count = Hash.new(0)
	atoms.each { |ap|
	  next if group != nil && !group.member?(ap.index)
	  e = ap.element
	  if (e == "Du")
	    ap.element = e = "H"
	  end
	  count[e] += 1
	  ap.name = sprintf("%s%d", e, count[e])
	}
	self
  end
  
  def guess_type_sub(idx)
    ap = atoms[idx]
	if ap.atom_type != ""
	  return ap.atom_type
	end
	e = ap.element
	n = ap.connects.length
	t = nil
	if e == "C"
	  if n == 2
		t = "c1"
	  elsif n == 3
		ap.connects.each { |i|
		  if atoms[i].element == "O" && atoms[i].connects.length == 1
			#  Carbonyl carbon
			t = "c"
			break
		  end
		}
		if t == nil
		  #  Not carbonyl carbon: may be other type, but it is difficult to 
		  #  guess the atom type so assign as "CA"
		  t = "ca"
		end
	  else
		t = "c3"
	  end
	elsif e == "H"
	  tt = guess_type_sub(ap.connects[0])
	  #  Count the number of electronegative atoms connected to the parent atom
	  ne = 0
	  atoms[ap.connects[0]].connects.each { |i|
		ee = atoms[i].element
		if ["N", "O", "F", "P", "S", "Cl", "As", "Se", "Br", "Sb", "Te", "I"].index(ee) != nil
		  ne += 1
		end
	  }
	  te = atoms[ap.connects[0]].element
	  if te == "C"
	    if tt == "ca"
	      if ne == 1
		    t = "h4"
		  elsif ne == 2
		    t = "h5"
		  else
	        t = "ha"
		  end
	    elsif tt == "cz"
	      t = "ha"
	    elsif tt == "c3"
		  if ne == 1
		    t = "h1"
		  elsif ne == 2
		    t = "h2"
		  elsif ne == 3
		    t = "h3"
		  else
		    t = "hc"
		  end
		else
		  t = "hc"
	    end
	  elsif te == "N"
	    t = "hn"
	  elsif te == "O"
		t = "ho"
	  elsif te == "S"
		t = "hs"
	  else
		t = "hc"
	  end
	elsif e == "N"
	  i = ap.connects.find { |j|
	    atoms[j].element == "C" && atoms[j].connects.find { |k|
		  atoms[k].element == "O" && atoms[k].connects.length == 1
		}
	  }
	  if i != nil
	    t = "n"   #  Amide NH
	  elsif n == 4
	    t = "n4"  #  Ammonium group
	  else
	    t = "n3"  #  Amino group
	  end
	elsif e == "O"
	  if n == 1
	    #  Count the number of oxygens connected to the root atom
		no = 0
		atoms[ap.connects[0]].connects.each { |i|
		  if atoms[i].element == "O"
		    no += 1
		  end
		}
		if no == 1
		  t = "o"  #  Carbonyl
		else
		  t = "o"  #  Carboxyl or phosphate
		end
	  else
	    #  Count the number of hydrogens connected to this atom
		nh = 0
		ap.connects.each { |i|
		  if atoms[i].element == "H"
		    nh += 1
		  end
		}
		if nh == 0
		  t = "os"
		elsif nh == 1
		  t = "oh"
		else
		  t = "ow"
		end
	  end
	else
	  t = e.upcase
	end
    ap.atom_type = t
  end

  #  Clean up atom types.
  def guess_types(group = nil)
	atoms.each { |ap|
	  next if group != nil && !group.member?(ap.index)
	  next if ap.atom_type != ""
	  guess_type_sub(ap.index)
	}
  end

  #  Solvate the molecule with the given solvent box.
  #  The first argument (box) must be a Molecule containing a unit cell information.
  #  The second argument defines the size of the solvated system. A positive number represents
  #  an offset to the bounding box of the solute, and a negative number represents an absolute size.
  #  If it is given as an Array (or a Vector), the elements define the sizes in the x/y/z directions.
  #  The third argument represents the limit distance to avoid conflict between the solute and
  #  solvent. The solvent molecule containing atoms within this limit from the solute is removed.
  #  The atom group containing added solvent molecule is returned.
  def solvate(sbox, size = [10.0, 10.0, 10.0], limit = 3.0)

	#  Sanity check
    if sbox.box == nil
	  raise MolbyError, "the solvent box does not have a unit cell information"
    end
	flags = sbox.box[4]
	if flags[0] == 0 || flags[1] == 0 || flags[2] == 0
	  raise MolbyError, "the solvent box does not have three-dimensional periodicity"
	end

	#  Calculate the box size
	b = self.bounds
	bsize = b[1] - b[0]
	if size.kind_of?(Numeric)
	  size = [size, size, size]
	end
	size = size.collect { |s| Float(s) }
	limit = Float(limit)
	(0..2).each do |i|
	  d = (size[i] >= 0 ? bsize[i] + size[i] * 2 : -size[i])
	  if d < bsize[i]
	    raise MolbyError, "the box size is too small"
	  end
	  bsize[i] = d;
	end
#	puts "Box size = #{bsize}"

	#  Translate the solute to the center of the box
	translate((b[0] + b[1]) * -0.5)
	solute_natoms = self.natoms

	#  Add solvents so that the target box is fully covered
	rtr = sbox.cell_transform.inverse
	min = Vector3D[1e30, 1e30, 1e30]
	max = Vector3D[-1e30, -1e30, -1e30]
	[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1]].each do |pt|
	  pt = Vector3D[(pt[0] - 0.5) * bsize[0], (pt[1] - 0.5) * bsize[1], (pt[2] - 0.5) * bsize[2]]
	  rpt = rtr * pt
	  min.x = rpt.x if min.x > rpt.x
	  min.y = rpt.y if min.y > rpt.y
	  min.z = rpt.z if min.z > rpt.z
	  max.x = rpt.x if max.x < rpt.x
	  max.y = rpt.y if max.y < rpt.y
	  max.z = rpt.z if max.z < rpt.z
	end
	xmin = (min.x + 0.5).floor
	xmax = (max.x + 0.5).floor
	ymin = (min.y + 0.5).floor
	ymax = (max.y + 0.5).floor
	zmin = (min.z + 0.5).floor
	zmax = (max.z + 0.5).floor
	sbox_natoms = sbox.natoms
	xv, yv, zv, ov, flags = sbox.box
#	puts "xmin = #{xmin}, ymin = #{ymin}, zmin = #{zmin}"
#	puts "xmax = #{xmax}, ymax = #{ymax}, zmax = #{zmax}"
	(xmin..xmax).each do |x|
	  (ymin..ymax).each do |y|
	    (zmin..zmax).each do |z|
		  add(sbox)
		  translate(xv * x + yv * y + zv * z, IntGroup[self.natoms - sbox_natoms..self.natoms - 1])
		end
	  end
	end
	
	#  Remove out-of-bounds molecules
	g = atom_group do |ap|
	  r = ap.r
	  r.x < -bsize[0] * 0.5 || r.y < -bsize[1] * 0.5 || r.z < -bsize[2] * 0.5 || r.x > bsize[0] * 0.5 || r.y > bsize[1] * 0.5 || r.z > bsize[2] * 0.5
	end
	g = fragment(g)  #  expand by fragment
	remove(g)
#	puts "Removed atoms by bounds: #{g}"
	
	#  Find conflicts
	conf = find_conflicts(limit, IntGroup[0..solute_natoms - 1], IntGroup[solute_natoms..self.natoms - 1])
	g = atom_group(conf.map { |c| c[1] } )    #  atom group containing conflicting atoms
	g = fragment(g)                           #  expand by fragment
	remove(g)
#	puts "Removed atoms by conflicts: #{g}"
	
	#  Renumber residue numbers of solvent molecules
	rseq = max_residue_number(0..solute_natoms - 1)
	rseq = 0 if rseq == nil
	each_fragment do |g|
	  next if g[0] < solute_natoms
	  rseq += 1
	  assign_residue(g, atoms[g[0]].res_name + ".#{rseq}")
	end

    #  Set the unit cell information
	set_box(bsize[0], bsize[1], bsize[2])
	
	return IntGroup[solute_natoms..self.natoms - 1]

  end	
  
  #  Count the (minimum) number of bonds between atoms n1 and n2
  def count_bonds(n1, n2)
    if n1 == n2
	  return 0
    end
    table = Array.new(natoms, 0) #  Results for each atom (offset by +1)
    table[n1] = 1                #  count_bonds(n1, n1) should be 0
    shell = [n1]                 #  The group of atoms with equal count_bonds
    n = 2                        #  Next shell number
    while 1
	  nshell = []
	  shell.each { |i|
	    atoms[i].connects.each { |j|
		  if table[j] == 0
		    if j == n2
			  return n - 1
		    end
		    table[j] = n
		    nshell.push(j)
		  end
	    }
	  }
	  if nshell.length == 0
	    return -1   #  No connection between n1 and n2
	  end
	  shell = nshell
	  n += 1
    end
  end

end
