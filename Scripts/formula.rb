#
#  main.rb
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

class Molecule

  @@known_fragments = Hash.new
  @@known_fragment_names = []
  
  #  The instance variable @dummies is used to retain the indices of the detachable
  #  atoms in the fragment. For example, a methyl group "CH3" is internally 
  #  represented by a 5-atom molecule "Du-CH3" where Du is a dummy atom. When
  #  a molecule is returned from the method Molecule.from_formula, the dummy atom
  #  is replaced with a hydrogen atom but the position of the dummy atom is retained
  #  as @dummies[0].

  attr_accessor :dummies

  #  Find dummy atoms; the atoms whose element is "Du" and name begines with "_" 
  #  are regarded as dummy atoms.
  def find_dummy_atoms(how_many = nil)
    d = []
    atoms.each { |ap|
	  if ap.element == "Du" && ap.name[0] == ?_
	    d.push(ap.index)
		if how_many != nil && d.length >= how_many
		  break
		end
	  end
	}
	d
  end
  
  #  Returns the first atom connected to the specified atom.
  def root_atom(dummy_atom)
    return atoms[dummy_atom].connects[0]
  end
  
  #  Register a fragment in @@known_fragment and rebuild @@known_fragment_names
  def Molecule.register_fragment(name, fragment)
    @@known_fragments[name] = fragment
	idx = @@known_fragment_names.length
	@@known_fragment_names.each_with_index { |n,i|
	  if n.length < name.length
	    idx = i
		break
	  end
	}
	@@known_fragment_names[idx, 0] = name
  end
  
  def Molecule.known_fragment(name)
	@@known_fragments[name].dup
  end
  
  def Molecule.lookup_fragment_from_string(str)
	@@known_fragment_names.each { |name|
	  if str.rindex(name, 0) == 0
		return name
	  end
	}
	nil
  end

  #  Create a fragment containing an atom and appropriate number of dummy atoms.
  #  The argument valence specifies the number of dummy atoms. Length is the
  #  bond length, and degree is the angle Du-X-Du (given in degree, not radian).
  def Molecule.atom_fragment(name, valence, length = nil, degree = nil)
    fragment = Molecule.new
	fragment.add_atom(name, "", name)  #  Atom type will be guessed later
	length = 1.0 if length == nil
	if valence > 0
	  fragment.add_atom("_1", "", "Du", length, 0)
	  if valence > 1
	    if degree == nil
		  degree = 180.0 if valence == 2
		  degree = 120.0 if valence == 3
		  degree = 109.2 if valence == 4
		end
		fragment.add_atom("_2", "", "Du", length, 0, degree, 1)
		if valence > 2
		  if valence == 3
		    sn = Math.sin(degree * Deg2Rad)
			cs = Math.cos(degree * Deg2Rad)
			acs = cs * (1 - cs) / (sn * sn)
			if (acs < -1 || acs > 1)
			  dihed = 180.0
			else
			  dihed = Math.acos(acs) * Rad2Deg
			end
		  elsif valence == 4
		    dihed = 120.0
		  end
		  fragment.add_atom("_3", "", "Du", length, 0, degree, 1, dihed, 2)
		  if valence > 3
		    fragment.add_atom("_4", "", "Du", length, 0, degree, 1, -dihed, 2)
		  end
		end
	  end
	end
	fragment.dummies = fragment.find_dummy_atoms
	fragment
  end

  #  Examine whether the given atom is a carbonyl carbon
  def is_carbonyl(idx)
    if atoms[idx].element == "C"
	  atoms[idx].connects.each { |idx2|
	    if atoms[idx2].element == "O" && atoms[idx2].connects.length == 1
		  return true
		end
	  }
	end
	return false
  end
  
  def guess_bond_length(e1, e2, mult)
    if e1 > e2
	  ee = e1
	  e1 = e2
	  e2 = ee
	end
	if (e1 == "C" && e2 == "C")
	  len = (mult == 3 ? 1.18 : (mult == 2 ? 1.39 : 1.54))
	elsif (e1 == "C" && e2 == "H")
	  len = 1.08
	elsif (e1 == "C" && e2 == "O")
	  len = (mult == 2 ? 1.23 : 1.41)
	elsif (e1 == "C" && e2 == "N")
	  len = (mult == 2 ? 1.29 : 1.48)
	elsif (e1 == "H" && e2 == "O")
	  len = 0.96
	elsif (e1 == "H" && e2 == "N")
	  len = 1.01
	elsif (e1 == "C" && e2 == "Cl")
	  len = 1.79
	elsif (e1 == "C" && e2 == "S")
	  len = 1.82
	elsif (e1 == "O" && e2 == "S")
	  len = (mult == 2 ? 1.60 : 1.80)
	else
	  len = 1.5
	end
	return len
  end
  
  #  Make the 'root' atom trigonal (sp2). The atoms in the dummies array are
  #  considered as one dummy atom.
  def make_sp2(root, dummies, vec)
    nd = atoms[root].connects.select { |n| !dummies.include?(n) }
	return if nd.length != 2
	r0 = atoms[root].r
	r1 = atoms[nd[0]].r - r0
	ax = r1.cross(atoms[nd[1]].r - r0)
	if (ax.length < 1e-4)
	  ax = r1.cross(atoms[dummies[0]].r - r0)
	  if (ax.length < 1e-4)
	    ax = r1.cross(Vector3D[1, 0, 0])
		if (ax.length < 1e-4)
		  ax = r1.cross(Vector3D[0, 1, 0])
		  if (ax.length < 1e-4)
		    raise("Cannot determine axis for rotating fragment")
		  end
		end
	  end
	end
	r2 = atoms[nd[1]].r - r0
	cs = r1.dot(r2) / (r1.length * r2.length)
	ang = Math.acos(cs) * Rad2Deg
	rotate(ax, 120.0 - ang, r0, fragment(nd[0], root))
	v = ((atoms[nd[0]].r - r0).normalize + (atoms[nd[1]].r - r0).normalize) * (-1.0)
	vec.x = v.x
	vec.y = v.y
	vec.z = v.z
	return vec
  end
  
  #  Add a molecular fragment to self. The new connection is created at the position
  #  of the dummy atoms. The 'first' dummy atom in self is kept untouched wherever possible.
  #  Mult (multiplicity) is the number of dummy atoms removed from each fragment. 
  #  If unspecified, it is the same as the number of the dummy atoms in the added fragment.
  def add_formula(mol, mult = nil)
    natoms_s = self.natoms
	natoms_m = mol.natoms
	#  Find list of dummy atoms
    md = mol.find_dummy_atoms
	len = md.length
	if len == 0
	  raise "The given molecule has no dummy atoms"
	end
	if mult == nil
	  mult = len
	elsif mult > len
	  raise "The given molecule has fewer dummy atoms (#{len}) than the bond order (#{mult})"
	else
	  md.slice!(mult..-1)
	end
	sd = self.find_dummy_atoms
	len2 = sd.length
	if len2 < mult
	  raise "Self has too few dummy atoms (#{len2}) to connect the given molecule (#{mult} bond order)"
	elsif len2 > mult
	  sd.shift  #  Keep the first dummy atom untouched
	  sd.slice!(mult..-1)
	end
	#  Find root atoms and the bond direction (which is defined by the average position of all dummy atoms)
	mroot = nil
	sroot = nil
	mvec = Vector3D[0, 0, 0]
	svec = Vector3D[0, 0, 0]
	md.each_index { |i|
	  mr = mol.root_atom(md[i])
	  if mroot == nil
	    mroot = mr
	  elsif mroot != mr
	    raise "The given molecule has multiple dummy atoms bonded to different root atoms"
	  end
	  mvec += mol.atoms[md[i]].r
	  sr = self.root_atom(sd[i])
	  if sroot == nil
	    sroot = sr
	  elsif sroot != sr
	    raise "self has multiple dummy atoms bonded to different root atoms"
	  end
	  svec += self.atoms[sd[i]].r
	}
	mvec *= 1.0 / mult
	svec *= 1.0 / mult
	mvec -= mol.atoms[mroot].r
	svec -= self.atoms[sroot].r
	#  Specify the atoms to define dihedral angles while docking
	#  (Preference: "untouched" dummy atoms > normal atoms > "removed" dummy atoms)
	mdihed = mol.atoms[mroot].connects.sort_by { |n|
	  if md.index(n) != nil
	    n - natoms_m
	  elsif mol.atoms[n].element == "Du"
	    n + natoms_m
	  else
	    n
	  end
	} [-1]
	sdihed = self.atoms[sroot].connects.sort_by { |n|
	  if sd.index(n) != nil
	    n + natoms_s
	  elsif self.atoms[n].element == "Du"
	    n - natoms_s
	  else
	    n
	  end
	} [0]

	#  Determine the bond length
	me = mol.atoms[mroot].element
	se = self.atoms[sroot].element
	len = guess_bond_length(me, se, mult)

	#  Kludge: special treatment for the amide bond
	if (mult == 1 && ((se == "N" && mol.is_carbonyl(mroot)) || (me == "N" && self.is_carbonyl(sroot))))
	  self.make_sp2(sroot, sd, svec)
	  mol.make_sp2(mroot, md, mvec)
	  len = 1.33
	end
	if (mult == 2 && se == "C")
	  self.make_sp2(sroot, sd, svec)
	end
	if (mult == 2 && me == "C")
	  mol.make_sp2(mroot, md, mvec)
	end

	#  Dock the two fragments
	self.dock(mol, sroot, mroot, svec, mvec, len, 180.0, sdihed, mdihed)
	#  Remove the dummy atoms
	md.map! { |item| item + natoms_s }
	self.remove(sd + md)
	#  Cache the list of dummy atoms
	@dummies = self.find_dummy_atoms
	self
  end
  
  #  Internal subroutine for from_formula
  def Molecule.from_formula_sub(str, idx)
    idx0 = idx
    f = nil
	c = str[idx]
	if c == ?(
	  f, idx = from_formula_sub(str, idx + 1)
	  if str[idx] != ?)
	    raise "Missing right parenthesis"
	  end
	  idx += 1
	  if (c = str[idx]) == nil || c < ?0 || c > ?9
	    #  No connection across the right parenthesis
	    return f, idx
	  end
	else
	  s = str[idx..-1]
	  key = lookup_fragment_from_string(s)
	  if key != nil
	    f = known_fragment(key)
		idx += key.length
	  end
	  if f == nil
	    ss = str[0..idx - 1] + "<?>" + str[idx..-1]
	    raise "Cannot find atom/fragment that matches \"#{s}\": #{ss}"
	  end
	end
	c = str[idx]
	valence = f.dummies.length
	if valence == 2 && c != nil && c >= ?0 && c <= ?9
	  #  Concatenate the fragment n times (like (CH2)8)
	  n = c - ?0
	  idx += 1
	  while (c = str[idx]) != nil && c >= ?0 && c <= ?9
		n = n * 10 + c - ?0
		idx += 1
	  end
	  if n == 0
		raise "0 times repeat is not allowed"
	  end
	  if n > 1
		f1 = f.dup
		(n - 1).times { f.add_formula(f1, 1) }
	  end
	end
	#  Add the fragment(s) from the following substring
	while true
	  if (c = str[idx]) == nil || c == ?)
	    #  End of string or parenthesized substring
	    return f, idx
	  end
	  if (valence = f.dummies.length) == 0
		#  This fragment is saturated
		return f, idx
	  elsif valence == 1 && idx0 > 0
	    #  This fragment has only one dummy bond and is not the first one
		return f, idx
	  end
	  ff, idx = from_formula_sub(str, idx)
	  if (c = str[idx]) != nil && (c >= ?0 && c <= ?9)
		#  The same fragment is added multiple times
	    n = c - ?0
	    idx += 1
	    while (c = str[idx]) != nil && c >= ?0 && c <= ?9
		  n = n * 10 + c - ?0
		  idx += 1
	    end
	    if n == 0
		  raise "0 times repeat is not allowed"
	    end
	  else
	    n = 1
	  end
	  n.times { f.add_formula(ff) }
	end # Repeat until the string ends or the fragment gets saturated
  end
  
  #  Convert a string to a molecule.
  def Molecule.from_formula(str, resname = nil)
    f, idx = from_formula_sub(str, 0)
	f.guess_names
#	f.guess_types
#	if f.nresidues != 2 || resname != nil
	  resname = "RES" if resname == nil
	  f.assign_residue(f.all, "#{resname}.1")
#	end
	f
  end
  
  #  Create a fragment from str and replace the group with the generated fragment.
  #  This method is invoked from MainView when a detachable selection is double-clicked
  #  and user enters the formula in the dialog box.
  def dock_formula(str, group = selection)
    if group.length == 1
	  #  Check if the string is an element name
	  Parameter.builtin.elements.each { |par|
		if par.atomic_number > 0 && par.name == str
		  ap = atoms[group[0]]
		  if ap.atomic_number != par.atomic_number
			if ap.name.index(ap.element) == 0
			  ap.name = ap.name.sub(ap.element, par.name)[0..3]
			end
			ap.atomic_number = par.atomic_number
			ap.atom_type = ""
			guess_types(group)
		  end
		  return
		end
	  }
	end
	mol = Molecule.from_formula(str)
	dock_fragment(mol, group)
  end
	
  #  Replace the specified group with the given fragment.
  def dock_fragment(fragment, group = selection)
    n0 = self.natoms
	return if fragment.natoms == 0
    if group.length == 0
	  add(fragment)
	  sel = atom_group(n0...n0 + fragment.natoms)
	else
	  frag = fragment.dup
      bonds_old = bonds_on_border(group)
   	  bonds_new = frag.dummies
	  if bonds_new == nil
	    bonds_new = frag.find_dummy_atoms()
	  end
#	  puts "group = #{group.inspect}, bonds_old = #{bonds_old.inspect}, bonds_new = #{bonds_new.inspect}"
#	  frag.dump
	  if bonds_old.length == 0 || bonds_new.length == 0
	    #  Translate the added fragment to the center of the selection
		v1 = self.center_of_mass(group) - frag.center_of_mass
		frag.translate(v1)
		remove(group)
		n0 = self.natoms
		add(frag)
		sel = atom_group(n0...n0 + frag.natoms)
	  else
	    old1 = bonds_old[0][1]  #  The first "root" atom in the molecule
		old2 = bonds_old[0][0]  #  The first "dummy" atom in the molecule
		new1 = frag.root_atom(bonds_new[0])  #  The first "root" atom in the fragment
		new2 = bonds_new[0]                  #  The first "dummy" atom in the fragment
		if bonds_old.length > 1 && bonds_new.length > 1
		  old3 = bonds_old[1][1]
		  old4 = bonds_old[1][0]
		  new3 = frag.root_atom(bonds_new[1])
		  new4 = bonds_new[1]
		  oldv1 = (atoms[old1].r - atoms[old2].r).normalize * 1.5
		  oldv2 = (atoms[old3].r - atoms[old4].r).normalize * 1.5
		  oldp1 = atoms[old1].r - oldv1
		  oldp2 = atoms[old3].r - oldv2
		  newp1 = frag.atoms[new1].r
		  newp2 = frag.atoms[new3].r
		  newv1 = frag.atoms[new2].r - newp1
		  newv2 = frag.atoms[new4].r - newp2
		  oldo = (oldp1 + oldp2) * 0.5
		  oldx = (oldp1 - oldp2).normalize
		  oldy = oldv1 + oldv2
		  oldz = oldx.cross(oldy).normalize
		  oldy = oldz.cross(oldx).normalize
		  newo = (newp1 + newp2) * 0.5
		  newx = (newp1 - newp2).normalize
		  newy = newv1 + newv2
		  newz = newx.cross(newy).normalize
		  newy = newz.cross(newx).normalize
          tr = Transform[oldx, oldy, oldz, oldo] * (Transform[newx, newy, newz, [0, 0, 0]].inverse) * Transform.translation(newo * (-1))
          frag.transform(tr)
		  add(frag)
		  create_bond(old1, new1 + n0)
		else
		  oldv = atoms[old2].r - atoms[old1].r
		  newv = frag.atoms[new2].r - frag.atoms[new1].r
		  #  Find the most substituted atoms (which will be arranged in trans conformation)
		  if (atoms[old1].connects.length >= 2)
		    olddi1 = atoms[old1].connects.sort_by { |n| (n == old2 ? 0 : atoms[n].connects.length) } [0]
		  else
		    olddi1 = nil
		  end
		  if (frag.atoms[new1].connects.length >= 2)
		    newdi1 = frag.atoms[new1].connects.sort_by { |n| (n == new2 ? 0 : frag.atoms[n].connects.length) } [0]
		  else
		    newdi1 = nil
		  end
		  dock(frag, old1, new1, oldv, newv, 1.5, 180.0, olddi1, newdi1)  #  Add (without removing atoms)
		end
		(1..bonds_old.length - 1).each { |i|
		  break if i >= bonds_new.length
		  #  Create bonds between other "root" atoms
		  create_bond(bonds_old[i][1], frag.root_atom(bonds_new[i]) + n0)
		}
		rem1 = frag.atom_group(bonds_new)
		rem = group + rem1.offset(n0)            #  The atoms to be removed
		sel = atom_group(n0...n0 + frag.natoms) - rem1.offset(n0)  #  The atoms to be selected
		sel = sel.deconvolute(atom_group - rem)  #  Renumber by skipping the atoms in rem
		remove(rem)
	  end
	end
	guess_types(sel)
	set_undoable_selection(sel)
  end
		

  #  Define molecular fragment from dump string
  def self.fragment_from_dump(str)
    f = Molecule.new.from_dump(str)
	f.dummies = f.find_dummy_atoms
	return f
  end

  #  Define atom fragments
  elements = %w(C H Li Be B C N O F Na Mg Al Si P S Cl)
  valences = [4,1,1,2,3,4,3,2,1,1,2,3,4,3,2,1]
  degrees = {"N"=>100, "O"=>110, "P"=>95, "S"=>108}
  elements.each_index { |i|
    name = elements[i]
	angle = degrees[name]
	register_fragment(name, atom_fragment(name, valences[i], nil, angle))
  }
  #  Special fragment to represent "open valence"
  register_fragment("_", atom_fragment("Du", 1, nil, nil))
  
  #  Define common molecular fragments
  formula_name = %w(CO CO2 O2C CH2 CH3 Et Pr iPr Bu tBu)
  formula_desc = ["CO", "C(O)O_", "OC(O)_", "CH2", "CH3", "CH2CH3", "CH2CH2CH3", "CH(CH3)2", "CH2CH2CH2CH3", "C(CH3)3"]
  formula_name.each_index { |i|
    f, idx = from_formula_sub(formula_desc[i], 0)
	#  Convert the "open valence" atoms to the dummy atoms
	f.atoms.select { |ap|
		ap.element == "Du" ? ap : nil
	}.each_with_index { |ap, j|
		ap.name = sprintf("_%d", j)
	}
	f.dummies = f.find_dummy_atoms
	register_fragment(formula_name[i], f)
  }
  register_fragment("Me", known_fragment("CH3"))
  load(MbsfPath + "/fragment_def.rb")
  
  #  Returns an arbitrary unit vector that is orthogonal to v
  def orthogonal_vector(v)
    vx = v.cross(Vector3D[1, 0, 0])
	if vx.length < 1e-3
	  vx = v.cross(Vector3D[0, 1, 0])
	end
	vx.normalize
  end
  
  #  Add missing hydrogen (or other atom) according to the given geometry type.
  #  atype can be one of the following:
  #  "td": tetrahedral, "tr": trigonal, "py": pyramidal (like amine nitrogen), "li": linear
  def add_hydrogen(idx, atype, bond = 1.07, anum = 1)
    nc = atoms[idx].connects.length
	p = atoms[idx].r
	if nc == 0
	  vx = Vector3D[1, 0, 0]
	else
	  v = atoms[idx].connects.map { |i| (atoms[i].r - p).normalize }
	  vx = Vector3D[0, 0, 0]
	  v.each { |v0| vx += v0 }
	  if nc == 3
	    if vx.length > 1e-3
		  vx = vx.normalize
		else
		  #  The three atoms are in trigonal positions
		  vx = v[0].cross(v[1]).normalize
		  if vx.length < 1e-3
		    vx = orthogonal_vector(v[0])
		  end
		end
	  elsif nc == 2
	    if vx.length < 1e-3
		  vx = orthognal_vector(v[0])
		else
		  vx = vx.normalize
		end
	    vy = v[0] - v[1]
		if vy.length < 1e-3
		  vy = orthogonal_vector(vx)
		else
		  vy = vy.normalize
		end
		vz = vx.cross(vy).normalize
	  elsif nc == 1
	    i = atoms[idx].connects[0]
	    #  Find most substituted atom connected to atoms[i]
		vy = nil
		if (atoms[i].connects.length >= 2)
		  j = atoms[i].connects.sort_by { |n| (n == i ? 0 : atoms[n].connects.length) }
		  vy = nil
		  j.each { |j0|
		    vy = vx.cross(atoms[j0].r - atoms[i].r)
		    if vy.length > 1e-3
			  #  The atoms j[0]-i-idx are not linear
		      vz = vx.cross(vy).normalize
			  vy = vy.normalize
			  break
		    end
		  }
		end
	    if !vy
		  vy = orthogonal_vector(vx)
		else
		  vy = vy.normalize
		end
		vz = vx.cross(vy).normalize
	  else
	    vx = Vector3D[1, 0, 0]
		vy = Vector3D[0, 1, 0]
		vz = Vector3D[0, 0, 1]
	  end
	end
	vn = []
	type = nil
    if atype == "td"
	  raise "The atom #{idx} already has #{nc} bonds" if nc >= 4
	  cs = -0.333333  #  cos(109.47)
	  sn =  0.942809  #  sin(109.47) = sqrt(2)*2/3
	  cs2 = -0.577359 #  cos(125.27)
	  sn2 =  0.816490 #  sin(125.27)
	  if nc == 3
		vn << -vx
	  elsif nc == 2
	    vn << vx * cs2 + vz * sn2
		vn << vx * cs2 - vz * sn2
	  else
	    vn << vx if nc == 0
	    vn << vx * cs + vy * sn
		vn << vx * cs - vy * sn * 0.5 + vz * sn * 0.86603
		vn << vx * cs - vy * sn * 0.5 - vz * sn * 0.86603
	  end
	  type = "hc" if anum == 1
	elsif atype == "tr"
	  raise "The atom #{idx} already has #{nc} bonds" if nc >= 3
	  if nc == 2
	    vn << -vx
	  else
	    vn << vx if nc == 0
		vn << vx * (-0.5) + vy * 0.86603
		vn << vx * (-0.5) - vy * 0.86603
	  end
	  type = "ha" if anum == 1
	elsif atype == "py"
	  raise "The atom #{idx} already has #{nc} bonds" if nc >= 3
	  cs = -0.292376  #  cos(107)
	  sn =  0.956303  #  sin(107)
	  cs2 = -0.491527 #  cos(119.4)
	  sn2 =  0.870862 #  sin(119.4)
	  if nc == 2
	    vn << vx * cs2 + vz * sn2
	  else
	    vn << vx if nc == 0
		vn << vx * cs + vy * sn
		vn << vx * cs + vy * (cs * (1 - cs) / sn) + vz * ((1 - cs) / sn * Math.sqrt(1 + 2 * cs))
	  end
	  type = "h2" if anum == 1
	elsif atype == "be"
	  raise "The atom #{idx} already has #{nc} bonds" if nc >= 2
	  cs = -0.258819  #  cos(105)
	  sn =  0.965926  #  sin(105)
	  vn << vx if nc == 0
	  vn << vx * cs + vy * sn
	  type = "ho" if anum == 1
	elsif atype == "li"
	  raise "The atom #{idx} already has #{nc} bonds" if nc >= 2
	  vn << vx if nc == 0
	  vn << -vx
	  type = "hz" if anum == 1
	else
	  raise "Unknown atom type #{atype}"
	end
	aname = Parameter.builtin.element(anum).name
	name = atoms[idx].name
	name = name.gsub(/\A[A-Za-z]*/, aname)[0..2]
	vn.each_with_index { |v, i|
	  ap = create_atom(sprintf("%s%d", name, i + 1), idx + i + 1)
	  ap.atomic_number = anum
	  ap.atom_type = (type || ap.element)
	  assign_residue(atom_group(ap.index), atoms[idx].res_seq)
	  ap.r = atoms[idx].r + v * bond
	  create_bond(idx, ap.index)
	}
	#  Update selection
#	if selection.count > 0
#	  sel = selection.convolute(IntGroup[0..idx, (idx + vn.count + 1)..(natoms - 1)])
#	  self.selection = sel
#	end
  end
  
  def add_hydrogen_on_group(group, atype, bond = 1.07, anum = 1)
    group.reverse_each { |i|
	  add_hydrogen(i, atype, bond, anum)
	}
    self
  end

end
