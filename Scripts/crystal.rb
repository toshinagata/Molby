#
#  crystal.rb
#
#  Created by Toshi Nagata.
#  Copyright 2012 Toshi Nagata. All rights reserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation version 2 of the License.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

#  Definition for use in ORTEP
class AtomRef
  def to_adc
    sym = self.symop
    if sym == nil
      idx = self.index + 1
      symcode = 55501
    else
      idx = sym[4] + 1
      symcode = (sym[1] + 5) * 10000 + (sym[2] + 5) * 1000 + (sym[3] + 5) * 100 + sym[0] + 1
    end
    return idx, symcode
  end
end

class Molecule

def export_ortep(fp)

  #  Create atom list
  hidden = atom_group { |ap| !is_atom_visible(ap.index) }
  hydrogen = self.show_hydrogens
  expanded = self.show_expanded
  atomlist = atom_group { |ap|
    (ap.element != "H" || hydrogen) &&
    (ap.symop == nil || expanded) &&
    (!hidden.include?(ap.index))
  }

  #  Title
  fp.printf "%-78.78s\n", self.name + ": generated by Molby at " + Time.now.to_s

  #  Cell parameters
  cp = self.cell
  fp.printf "%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f\n", cp[0], cp[1], cp[2], cp[3], cp[4], cp[5]
  
  #  Symmetry operations
  syms = self.symmetries
  if syms == nil || syms.length == 0
   fp.print "1             0  1  0  0              0  0  1  0              0  0  0  1\n"
  else
    syms.each_with_index { |s, i|
      a = s.to_a
      fp.printf "%s%14g%3g%3g%3g%15g%3g%3g%3g%15g%3g%3g%3g\n", (i == syms.length - 1 ? "1" : " "), a[9], a[0], a[1], a[2], a[10], a[3], a[4], a[5], a[11], a[6], a[7], a[8]
    }
  end

  #  Atoms (all symmetry unique atoms regardless they are visible or not)
  n = 0
  each_atom { |ap|
    break if ap.symop != nil
    fp.printf " %4.4s%22s%9.4f%9.4f%9.4f%9d\n", ap.name, "", ap.fract_x, ap.fract_y, ap.fract_z, 0
    an = ap.aniso
    if an != nil
      fp.printf " %8.5f%9.6f%9.6f%9.6f%9.6f%9.6f%9d\n", an[0], an[1], an[2], an[3], an[4], an[5], 0
    else
      t = ap.temp_factor
      t = 1.2 if t <= 0
      fp.printf " %8.3f%9g%9g%9g%9g%9g%9d\n", t, 0.0, 0.0, 0.0, 0.0, 0.0, 6
    end
    n += 1
  }
  natoms_tep = n

  #  Special points to specify cartesian axes
  axis, angle = self.get_view_rotation
  tr = Transform.rotation(axis, angle)
  org = self.get_view_center
  x = org + tr.column(0)
  y = org + tr.column(1)
  tr = self.cell_transform.inverse
  org = tr * org
  x = tr * x
  y = tr * y
  fp.printf " CNTR                      %9.4f%9.4f%9.4f        0\n", org.x, org.y, org.z
  fp.printf " %8.3f%9g%9g%9g%9g%9g%9d\n", 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 6
  fp.printf " X                         %9.4f%9.4f%9.4f        0\n", x.x, x.y, x.z
  fp.printf " %8.3f%9g%9g%9g%9g%9g%9d\n", 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 6
  fp.printf " Y                         %9.4f%9.4f%9.4f        0\n", y.x, y.y, y.z
  fp.printf "1%8.3f%9g%9g%9g%9g%9g%9d\n", 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 6

  #  Initialize
  fp.print  "      201\n"

  #  Pen size
  fp.print  "      205        8\n"

  #  Paper size, margin, viewing distance
  fp.print  "      301      6.6      6.6        0      0.0\n"

  #  Sort atoms by symop and index
  acodes = Hash.new
  ig = IntGroup.new   #  Used for collecting contiguous atoms
  each_atom(atomlist) { |ap|
    idx, symcode = ap.to_adc
    acode = symcode * self.natoms + idx - 1
    acodes[acode] = ap.index
    ig.add(acode)
  }
  index2an = []       #  Index in Molecule to Index in ORTEP
  ig.each_with_index { |acode, i|
    index2an[acodes[acode]] = i + 1
  }
  i = 0
  adcs = []
  while (r = ig.range_at(i)) != nil
    s = r.first
    s = (s / self.natoms) + (s % self.natoms + 1) * 100000  #  Rebuild true ADC (atom index is in the upper digit)
    e = r.last
    e = (e / self.natoms) + (e % self.natoms + 1) * 100000
    if s < e
      adcs.push(s)
      adcs.push(-e)
    else
      adcs.push(s)
    end
    i += 1
  end
  k = 0

  #  Atom list
  adcs.each_with_index { |a, i|
    if k == 0
      fp.print "      401"
    end
    fp.printf "%9d", a
    k += 1
    if i == adcs.length - 1 || k == 6 || (k == 5 && i < adcs.length - 2 && adcs[i + 2] < 0)
      fp.print "\n"
      k = 0
    end
  }

  #  Axes
  fp.printf "      501%4d55501%4d55501%4d55501%4d55501%4d55501                 1\n", natoms_tep + 1, natoms_tep + 1, natoms_tep + 2, natoms_tep + 1, natoms_tep + 3
#  fp.print  "      502        1      0.0        2      0.0        3      0.0\n"
  
  #  Autoscale
  fp.print  "      604                               1.538\n"
  
  #  Explicit bonds
  bond_inst = Array.new(6) { [] }   #  Bonds for types 1 to 5
  bonds.each { |b|
    next if !atomlist.include?(b[0]) || !atomlist.include?(b[1])
    #  TODO: determination of bond types should be refined
    an1 = atoms[b[0]].atomic_number
    an2 = atoms[b[1]].atomic_number
    if an1 == 1 || an2 == 1
      btype = 1
    elsif an1 <= 8 && an2 <= 8
      btype = 3
    else
      btype = 5
    end
    bond_inst[btype].push(b[0], b[1])
  }

  #  Output bond specifications
  #  Avoid including too many ADCs in a single 811/821 instruction
  #  (Upper limit is 140. Here we divide at every 36 ADCs)
  output_bonds = proc { |icode|
    bond_inst.each_with_index { |inst, ii|
      next if inst.length == 0
      inst.each_with_index { |b, i|
        if i % 6 == 0
          fp.printf "  %d   %3s", (i >= inst.length - 6 || i % 36 == 30 ? 2 : 1), (i % 36 == 0 ? icode.to_s : "")
        end
        idx, scode = atoms[b].to_adc
        fp.printf "%9d", idx * 100000 + scode
        if i % 6 == 5 || i == inst.length - 1
          fp.print "\n"
          if i == inst.length - 1 || i % 36 == 35
            fp.printf "%21s%3d%12s%6.3f\n", "", ii, "", 0.05
          end
        end
      }
    }
  }

  fp.print "  0  1001     0.02\n"  #  Activate hidden line removal
  output_bonds.call(821)
  
  #  Atom types
  atom_inst = Array.new(5) { IntGroup.new }   #  Atoms for 714 (sphere), 712 (with axes), 711 (shade and axes)
  atomlist.each { |i|
    #  TODO: determination of atom types should be refined
    an1 = atoms[i].atomic_number
    if an1 == 1
      atype = 4
    elsif an1 <= 6
      atype = 2
    else
      atype = 1
    end
    idx, scode = atoms[i].to_adc
    atom_inst[atype].add(idx)
  }
  [4,2,1].each { |ii|
    inst = atom_inst[ii]
    i = 0
    while (r = inst.range_at(i)) != nil
      fp.printf "  1   %3d\n", 710 + ii
      fp.printf "%27s%9d%9d\n", "", r.first, r.last
      i += 1
    end
  }

  output_bonds.call(811)

  #  Close plot
  fp.print "      202\n"
  fp.print "       -1\n"
  
end

def savetep(filename)
  if natoms == 0
	raise MolbyError, "cannot save ORTEP input; the molecule is empty"
  end
  fp = open(filename, "wb")
  export_ortep(fp)
  fp.close
  return true
end

end

#  Best-fit planes
#  Ref. W. C. Hamilton, Acta Cryst. 1961, 14, 185-189
#       T. Ito, Acta Cryst 1981, A37, 621-624

#  An object to described the best-fit plane
#  
#  Contains the plane coefficients and the constant (a, b, c, d for ax + by + cz + d = 0),
#  the error matrix, and the metric tensor.
#
class Molby::Plane
  attr_accessor :molecule, :group, :coeff, :const, :error_matrix, :metric_tensor
  def initialize(mol, group, coeff, const, err, met)
    @molecule = mol
    @group = group
    @coeff = coeff
    @const = const
    @error_matrix = err
    @metric_tensor = met
    self
  end
  def sigma
    [sqrt(@error_matrix[0, 0]), sqrt(@error_matrix[1, 1]), sqrt(@error_matrix[2, 2]), sqrt(@error_matrix[3, 3])]
  end
  def inspect
    s = sprintf("Molby::Plane[\n coeff, const = [[%f, %f, %f], %f],\n", @coeff.x, @coeff.y, @coeff.z, @const)
    s += sprintf(" sigma = [[%10.4e, %10.4e, %10.4e], %10.4e],\n", *self.sigma)
    (0..3).each { |i|
      s += (i == 0 ? " error_matrix = [" : "     ")
      (0..i).each { |j|
        s += sprintf("%12.6e%s", @error_matrix[j, i], (j == i ? (i == 3 ? "],\n" : ",\n") : ","))
      }
    }
    s += sprintf(" molecule = %s\n", @molecule.inspect)
    s += sprintf(" group = %s\n", @group.inspect)
    (0..3).each { |i|
      s += (i == 0 ? " metric_tensor = [" : "     ")
      (0..3).each { |j|
        s += sprintf("%12.6e%s", @metric_tensor[j, i], (j == 3 ? (i == 3 ? "]]\n" : ",\n") : ","))
      }
    }
    s
  end
  def distance(ap)
    if ap.is_a?(AtomRef)
      fr = ap.fract_r
      sig = ap.sigma
    else
      fr = Vector3D[*ap]
      sig = Vector3D[0, 0, 0]
    end
    d = fr.dot(@coeff) + @const
    sig1 = (@coeff.x * sig.x) ** 2 + (@coeff.y * sig.y) ** 2 + (@coeff.z * sig.z) ** 2
    sig2 = LAMatrix.multiply("t", fr, @error_matrix, fr)[0, 0]
    if ap.is_a?(AtomRef) && ap.molecule == @molecule && @group.include?(ap.index)
      #  The atom defines the plane
      sig0 = sig1 - sig2
      sig0 = 0.0 if sig0 < 0.0
    else
      sig0 = sig1 + sig2
    end  
    return d, sqrt(sig0)
  end
  def dihedral(plane)
    e1 = @error_matrix.submatrix(0, 0, 3, 3)
    e2 = plane.error_matrix.submatrix(0, 0, 3, 3)
    m = @metric_tensor.submatrix(0, 0, 3, 3)
    c = plane.coeff
    cos_t = plane.coeff.dot(m * @coeff)
    if cos_t > 1.0
      cos_t = 1.0
    elsif cos_t < -1.0
      cos_t = -1.0
    end
    t = acos(cos_t)
    sig_t = (m * e1).trace + (m * e2).trace
    if sig_t < t * t
      w = 1.0 / sin(t)
      sig_t = w * w * (c.dot(LAMatrix.multiply(m, e1, m, c)) + @coeff.dot(LAMatrix.multiply(m, e2, m, @coeff)))
    end
    t *= 180.0 / PI
    sig_t = sqrt(sig_t) * 180.0 / PI
    return t, sig_t
  end
end

class Molecule

#  Calculate best-fit plane for the given atoms
#  Return value: a Molby::Plane object

def plane(group)

  #  Number of atoms
  dim = group.length

  #  Positional parameters and standard deviations
  x = []
  sig = []
  sig_min = 1e10
  each_atom(group) { |ap|
    x.push(ap.fract_r)
    sig.push(ap.sigma)
    if (s = ap.sigma_x) > 0.0 && s < sig_min
      sig_min = s
    end
    if (s = ap.sigma_y) > 0.0 && s < sig_min
      sig_min = s
    end
    if (s = ap.sigma_z) > 0.0 && s < sig_min
      sig_min = s
    end
  }
  if sig_min == 1e10
    sig_min = 1e-12
  end
  sig.each { |s|
    s.x = sig_min if s.x < sig_min
    s.y = sig_min if s.y < sig_min
    s.z = sig_min if s.z < sig_min
  }

  #  The metric tensor of the reciprocal lattice
  #  g[j, i] = (ai*).dot(aj*), where ai* and aj* are the reciprocal axis vectors
  t = self.cell_transform
  if t.nil?
    t = Transform.identity
  end
  g2inv = LAMatrix[t]
  g2 = g2inv.inverse
  g2[3, 3] = 0.0
  g2inv[3, 3] = 0.0
  g = LAMatrix.multiply("t", g2, g2)

  #  The variance-covariance matrices of the atomic parameters
  #  mm[k][n] is a 3x3 matrix describing the correlation between the atoms k and n,
  #  and its components are defined as: sigma_k[i] * sigma_n[j] * corr[k, i, n, j],
  #  where corr(k, i, n, j) is the correlation coefficients between the atomic parameters
  #  k[i] and n[j].
  mm = Array.new(dim) { Array.new(dim) }
  zero = LAMatrix.zero(3, 3)
  dim.times { |k|
    dim.times { |n|
      mkn = LAMatrix.new(3, 3)
      if k == n
        3.times { |i|
          3.times { |j|
            if i == j
              mkn[j, i] = sig[k][i] * sig[n][j]
            else
              #  Inter-coordinate correlation should be implemented here
            end
          }
        }
      else
        #  Inter-atomic correlation should be implemented here
      end
      mm[k][n] = (mkn == zero ? zero : mkn)
    }
  }

  #  The variance-covariance matrix of the atom-plance distances
  #  m[j, i] = v.transpose * mm[i][j] * v, where v is the plane coefficient vector
  #  The inverse of m is the weight matrix
  m = LAMatrix.new(dim, dim)
  
  #  The matrix representation of the atomic coordinates
  #  y[j, i] = x[i][j] (for j = 0..2), -1 (for j = 3)
  #  y * LAMatrix[a, b, c, d] gives the atom-plane distances for each atom
  y = LAMatrix.new(4, dim)
  dim.times { |i|
    y[0, i] = x[i].x
    y[1, i] = x[i].y
    y[2, i] = x[i].z
    y[3, i] = 1.0
  }

  #  The coefficients to be determined
  n0 = LAMatrix[1, 1, 1, 0]
  v = LAMatrix[1, 1, 1]     #  The coefficient part

  iter = 0
  while iter < 20

    iter += 1

    #  Set zero to the "constant" part, and normalize the "coefficient" part
    n0[0, 3] = 0.0
    n0 = g2 * n0
    n0.multiply!(1.0 / n0.fnorm)
    n0 = g2inv * n0
    3.times { |i| v[0, i] = n0[0, i] }

    #  Build the variance-covariance matrix    
    dim.times { |i|
      dim.times { |j|
        m[j, i] = LAMatrix.multiply("t", v, mm[i][j], v)[0, 0]
      }
    }
    c = LAMatrix.multiply("t", y, "i", m, y)

    #  Invert c: only the inverse is used in the following, so c is inversed destructively
    cinv = c.inverse!
 
    if iter == 1

      #  Determine the tentative solution, which is given by the eigenvector of cinv * g
      #  for the largest eigenvalue
      evals, evecs = (cinv * g).eigenvalues
      4.times { |i| n0[0, i] = evecs[3, i] }

    else

      #  Convert the coefficient vector to the reciprocal space
      h = g * n0
      
      #  Determine multiplier
      #  In this implementation, the sign of delta-n is opposite from that used in
      #  the reference
      lam = 1.0 / (LAMatrix.multiply("t", h, cinv, h)[0, 0])
      
      #  Solve the linearized equation
      #  (Is the equation 21 in the reference really correct? Shouldn't it read
      #   B = 1 - lambda * C.inverse * H* ? )
      b = LAMatrix.multiply(lam, cinv, g)
      b.sub!(LAMatrix.identity(4))

      dn = b * n0
      n0 += dn

      break if dn[0, 0] ** 2 + dn[0, 1] ** 2 + dn[0, 2] ** 2 < 1e-9

    end
  end

  #  Error matrix = b * cinv * b.transpose
  em = LAMatrix.multiply(b, cinv, "t", b)
  coeff = Vector3D[n0[0, 0], n0[0, 1], n0[0, 2]]
  const = n0[0, 3]

  return Molby::Plane.new(self, group, coeff, const, em, g)

end

def cmd_plane
  plane_settings = @plane_settings || Hash.new
  mol = self
  h = Dialog.run("Best-Fit Planes", "Close", nil) {
    refresh_proc = proc { |it|
      n = it[:tag][/\d/].to_i
      g = plane_settings["group#{n}"]
      if g
        str = g.inspect.sub!("IntGroup[", "").sub!("]", "")
        set_value("group#{n}", str)
        if n == 1 || n == 2
          p = mol.plane(g) rescue p = nil
          plane_settings["plane#{n}"] = p
          if p
            coeff = p.coeff
            const = p.const
            sig = p.sigma
            aps = (n == 1 ? "" : "'")
            str = sprintf("a%s = %f(%f)\nb%s = %f(%f)\nc%s = %f(%f)\nd%s = %f(%f)",
                          aps, coeff.x, sig[0],
                          aps, coeff.y, sig[1],
                          aps, coeff.z, sig[2],
                          aps, const, sig[3])
            set_value("result#{n}", str)
          else
            set_value("result#{n}", "")
          end
          p1 = plane_settings["plane1"]
          p2 = plane_settings["plane2"]
          if p1 && p2
            t, sig = p1.dihedral(p2)
            str = sprintf("%f(%f)", t, sig)
            set_value("dihedral", str)
          else
            set_value("dihedral", "")
          end
        else
          p = plane_settings["plane1"]
          if p
            str = ""
            mol.each_atom(g) { |ap|
              d, sig = p.distance(ap)
              str += sprintf("%d %f(%f)\n", ap.index, d, sig)
            }
            str.chomp!
          else
            str = ""
          end
          set_value("result#{n}", str)
        end
      else
        set_value("group#{n}", "")
        set_value("result#{n}", "")
      end
    }
    set_proc = proc { |it|
      n = it[:tag][/\d/].to_i
      sel = mol.selection
      if sel.count > 0
        str = sel.inspect.sub!("IntGroup[", "").sub!("]", "")
        set_value("group#{n}", str)
        plane_settings["group#{n}"] = sel
      else
        plane_settings["group#{n}"] = nil
      end
      refresh_proc.call(it)
    }
    text_proc = proc { |it|
      n = it[:tag][/\d/].to_i
      str = it[:value].gsub(/[^-.,0-9]/, "")  #  Remove unsane characters
      g = eval("IntGroup[#{str}]") rescue g = nil
      plane_settings["group#{n}"] = g
      refresh_proc.call(it)
    }
    layout(3,
      item(:text, :title=>"Plane 1 (ax + by + cz + d = 0)"),
      -1, -1,
      item(:text, :title=>"Atoms"),
      item(:textfield, :width=>240, :height=>32, :tag=>"group1", :action=>text_proc),
      item(:button, :title=>"Set Current Selection", :tag=>"button1", :action=>set_proc),
      item(:text, :title=>"Results"),
      item(:textview, :width=>240, :height=>68, :editable=>false, :tag=>"result1"),     
      item(:button, :title=>"Recalculate", :tag=>"refresh1", :action=>refresh_proc),
      item(:line),
      -1, -1,
      item(:text, :title=>"Plane 2 (a'x + b'y + c'z + d' = 0)"),
      -1, -1,
      item(:text, :title=>"Atoms"),
      item(:textfield, :width=>240, :height=>32, :tag=>"group2", :action=>text_proc),
      item(:button, :title=>"Set Current Selection", :tag=>"button2", :action=>set_proc),
      item(:text, :title=>"Results"),
      item(:textview, :width=>240, :height=>68, :editable=>false, :tag=>"result2"),
      item(:button, :title=>"Recalculate", :tag=>"refresh2", :action=>refresh_proc),
      item(:text, :title=>"Dihedral angle with Plane 1"), -1, -1,
      -1,
      item(:textfield, :width=>240, :height=>16, :tag=>"dihedral"), -1,
      item(:line),
      -1, -1,
      item(:text, :title=>"Distance from Plane 1"), -1, -1,
      item(:text, :title=>"Atoms"),
      item(:textfield, :width=>240, :height=>32, :tag=>"group3", :action=>text_proc),
      item(:button, :title=>"Set Current Selection", :tag=>"button3", :action=>set_proc),
      item(:text, :title=>"Results"),
      item(:textview, :width=>240, :height=>68, :editable=>false, :tag=>"result3"),
      item(:button, :title=>"Recalculate", :tag=>"refresh3", :action=>refresh_proc)
    )
    refresh_proc.call(item_with_tag("refresh1"))
    refresh_proc.call(item_with_tag("refresh2"))
  }
  @plane_settings = plane_settings
end

if lookup_menu("Best-fit Planes...") < 0
  register_menu("", "")
  register_menu("Best-fit Planes...", :cmd_plane)
end

end
