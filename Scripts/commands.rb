# coding: utf-8
#
#  commands.rb
#
#  Created by Toshi Nagata on 2008/06/28.
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

  def cmd_assign_residue
    sel = self.selection
	if sel.length == 0
	  sel = self.atom_group
	end
	atoms = sel.inspect.sub!("IntGroup[", "").sub!("]", "")
    hash = Dialog.run("Assign Residue") {
	  layout(2,
		item(:text, :title=>"New residue name/number\n(like \"RES.1\")\nfor atoms #{atoms}"),
	    item(:textfield, :width=>120, :tag=>"residue"))
    }
    if hash[:status] == 0
	  residue = hash["residue"]
	  assign_residue(sel, residue)
	end
  end

  def cmd_offset_residue
    sel = self.selection
	if sel.length == 0
	  sel = self.atom_group
	end
	atoms = sel.inspect.sub!("IntGroup[", "").sub!("]", "")
    hash = Dialog.run("Offset Residues") {
	  layout(2,
		item(:text, :title=>"Offset residue number:\nfor atoms #{atoms}"),
	    item(:textfield, :width=>120, :tag=>"offset"))
    }
	if hash[:status] == 0
	  offset = hash["offset"].to_i
	  offset_residue(sel, offset)
	end
  end
   
  def cmd_sort_by_residue
    sel = self.selection
	if sel.length == 0
	  sel = self.atom_group
	end
	sorted = sel.sort_by { |i| [self.atoms[i].res_seq, i] }
	ary = []
	j = 0
	(0...natoms).each { |i|
	  if sel.include?(i)
	    ary << sorted[j]
		j += 1
		break if j >= sorted.length
	  else
	    ary << i
	  end
	}
	self.renumber_atoms(ary)
  end

  def cmd_delete_frames
    n = nframes
    return if n == 0
	hash = Dialog.run("Delete Frames") {
	  layout(2,
	    item(:text, :title=>"Start"),
	    item(:textfield, :width=>120, :tag=>"start", :value=>"0"),
		item(:text, :title=>"End"),
		item(:textfield, :width=>120, :tag=>"end", :value=>(n - 1).to_s),
		item(:text, :title=>"Keeping frames every..."),
		-1,
		item(:text, :title=>"Step"),
		item(:textfield, :width=>120, :tag=>"step", :value=>"0"))
	}
	if hash[:status] == 0
	  sframe = Integer(hash["start"])
	  eframe = Integer(hash["end"])
	  step = Integer(hash["step"])
	  return if sframe > eframe
	  eframe = n - 1 if eframe >= n
	  fgroup = IntGroup[sframe..eframe]
	  if step > 0
	    while sframe <= eframe
		  fgroup.delete(sframe)
		  sframe += step
		end
	  end
	  remove_frames(fgroup)
	end
  end

  def cmd_reverse_frames
    n = nframes
    return if n == 0
	hash = Dialog.run("Reverse Frames") {
	  layout(2,
	    item(:text, :title=>"Start"),
	    item(:textfield, :width=>120, :tag=>"start", :value=>"0"),
		item(:text, :title=>"End"),
		item(:textfield, :width=>120, :tag=>"end", :value=>(n - 1).to_s))
	}
	if hash[:status] == 0
	  sframe = Integer(hash["start"])
	  eframe = Integer(hash["end"])
	  return if sframe > eframe
	  eframe = n - 1 if eframe >= n
	  frames = (0...sframe).to_a + (sframe..eframe).to_a.reverse + ((eframe + 1)...n).to_a
	  reorder_frames(frames)
	end
  end

  def cmd_concat_frames
      n = nframes
      return if n == 0
      if Molecule.list.length == 1
          message_box("No other molecule.", "", :ok)
          return
      end
      ls = (0...Molecule.list.length).select { |i|
          m = Molecule[i]
          if m == self || m.natoms != self.natoms
              false
              else
              if nil != (0...m.natoms).find { |n| m.atoms[n].atomic_number != self.atoms[n].atomic_number }
                  false
                  else
                  true
              end
          end
      }
      if ls.length == 0
          message_box("No molecule has the same atomic sequence as the present one.", "", :ok)
          return
      end
      labels = []
      label_hash = Hash.new   #  Check if documents with the same name are present
      ls.each_with_index { |i, idx|  #  i is an index to Molecule, idx is an index to ls
          m = Molecule[i]
          label = m.name
          idx_h = label_hash[label]  #  If non-nil, then duplicate name
          label_hash[label] = idx
          if idx_h
              if labels[idx_h] == label
                  labels[idx_h] = label + " (" + Molecule[ls[idx_h]].dir + ")"
              end
              label = label + " (" + Molecule[i].dir + ")"
          end
          labels.push(label)
      }
      nf = Molecule[ls[0]].nframes
      hash = Dialog.run("Concatenate Frames") {
          layout(2,
                 item(:text, :title=>"From Molecule:"),
                 item(:popup, :subitems=>labels, :tag=>"mol", :width=>240,
                      :action=>lambda { |it|
                      nf = Molecule[ls[it[:value]]].nframes
                      set_attr("start_title", :title=>"Start (0-#{nf})")
                      set_attr("end_title", :title=>"End (0-#{nf})")
                      }),
                 item(:radio, :title=>"All Frames", :tag=>"all_frames", :value=>1,
                      :action=>lambda { |it|
                      flag = (it[:value] == 0)
                      set_attr("start_frame", :enabled=>flag)
                      set_attr("end_frame", :enabled=>flag)
                      }), -1,
                 item(:radio, :title=>"Select Frames", :tag=>"select_frames",
                      :action=>lambda { |it|
                      flag = (it[:value] != 0)
                      set_attr("start_frame", :enabled=>flag)
                      set_attr("end_frame", :enabled=>flag)
                      }), -1,
                 item(:text, :title=>"Start (0-#{nf})", :tag=>"start_title"),
                 item(:textfield, :value=>"0", :tag=>"start_frame", :enabled=>false),
                 item(:text, :title=>"End (0-#{nf})", :tag=>"end_title"),
                 item(:textfield, :value=>"0", :tag=>"end_frame", :enabled=>false))
                 radio_group("all_frames", "select_frames")
      }
      if hash[:status] == 0
          idx_h = Integer(hash["mol"])
          m = Molecule[ls[idx_h]]
          f = m.frame  #  Save the current frame number
          nf = m.nframes
          if hash["all_frame"] != 0
              f1 = 0
              f2 = nf - 1
              else
              f1 = Integer(hash["start_frame"])
              f2 = Integer(hash["end_frame"])
              f1 = 0 if f1 < 0
              f1 = nf - 1 if f1 > nf - 1
              f2 = 0 if f2 < 0
              f2 = nf - 1 if f2 > nf - 1
          end
          if f1 <= f2
              a = (f1..f2).to_a
              else
              a = (f2..f1).to_a.reverse
          end
          prop = m.property_names
          na = m.natoms
          a.each { |n|
              self.create_frame()
              m.select_frame(n)
              na.times { |i|
                  self.atoms[i].r = m.atoms[i].r
              }
              prop.each { |pr|
                  self.set_property(pr, m.get_property(pr))
              }
          }
      end
  end

  def cmd_extra_properties
    mol = self
	get_count = lambda { |it| mol.nframes }
	get_value = lambda { |it, row, col|
	  if col == 0
	    row.to_s
	  else
	    sprintf("%.8f", mol.get_property(col - 1, row)) rescue ""
	  end
    }
    mol.open_auxiliary_window("Extra Props", :resizable=>true, :has_close_box=>true) {
	  names = nil
	  @on_document_modified = lambda { |m|
	    if (n = m.property_names) != names
		  names = n
		  col = [["frame", 40]] + names.map { |nn| [nn, 120] }
		  item_with_tag("table")[:columns] = col
		end
		item_with_tag("table")[:refresh] = true
	  }
	  layout(1,
		item(:table, :width=>320, :height=>300, :flex=>[0,0,0,0,1,1], :tag=>"table",
		  :on_count=>get_count,
		  :on_get_value=>get_value),
	    :flex=>[0,0,0,0,1,1]
	  )
	  set_min_size(320, 200)
	  # listen(mol, "documentModified", on_document_modified)
	  # listen(mol, "documentWillClose", lambda { |m| close } )
	  @on_document_modified.call(mol)
	  show
	}
  end
  
  def cmd_show_energy
    wave = [0.0, 0.0]
    cur = 0
    mol = self
    wave_min = lambda { m = 1e8; wave.each { |x| if x != 0.0 && x < m; m = x; end }; m }
    wave_max = lambda { m = -1e8; wave.each { |x| if x != 0.0 && x > m; m = x; end }; m }
    d = open_auxiliary_window("Energy", :resizable=>true, :has_close_box=>true) {
    graph_item = nil   #  Forward declaration
    target_mol = nil
    draw_graph = lambda { |it|
      clear
      f = graph_item[:frame]
      draw_rectangle(0, 0, f[2], f[3])
      width = f[2] - 25
      height = f[3] - 25
      draw_line(16, 0, 16, height + 12, width + 20, height + 12)
      xx = yy = nil
      min = wave_min.call
      max = wave_max.call
      h = max - min
      h = (h == 0.0 ? 1.0 : height / h)
      c = wave.count
      if c == 0
        w = 1.0
      elsif c == 1
        w = Float(width)
      else
        w = Float(width) / (c - 1)
      end
      lines = []
      a = []
      #  Skip the points that have exactly 0.0 value
      wave.each_with_index { |d, i|
        if d != 0.0
          a.push(i * w + 16)
          a.push(height - (d - min) * h + 12)
        end
        if d == 0.0 || i == wave.length - 1
          #  End of this curve fragment
          if a.length == 0
            #  Do nothing
          else
            if a.length == 2
              if wave.count == 1
                #  If wave has only one point, then draw a horizontal line
                a.push(a[0] + 16)
                a.push(a[1])
              else
                #  Otherwise, draw a zero-length line
                a.push(a[0])
                a.push(a[1])
              end
            end
            lines.push(a)  #  Store this line fragment
            a = []
          end
        end
      }
      lines.each { |a|
        draw_line(a)
      }
      brush(:color=>[0.2, 0.2, 1.0])
      y = wave[cur] || 0.0
      if y < min
        y = min
      elsif y > max
        y = max
      end
      xx = cur * w + 16
      yy = height - (y - min) * h + 12
      draw_ellipse(cur * w + 16, height - (y - min) * h + 12, 6)
    }
	  @on_document_modified = lambda { |m|
		cur = mol.frame
		wave.clear
		if mol.nframes < 2
		  wave = [mol.get_property("energy", 0)] * 2
		else
		  wave = (0...mol.nframes).map { |i| mol.get_property("energy", i) }
		end
		f = graph_item[:frame]
		item_with_tag("energy")[:title] = sprintf("Energy = %.10f hartree, Frame = %d", mol.get_property("energy", cur), cur)
		graph_item.refresh_rect([0, 0, f[2], f[3]])
	  }
	  layout(1,
		item(:text, :title=>"Energy =                 ", :tag=>"energy"),
		item(:view, :frame=>[0, 0, 320, 240], :tag=>"graph", :on_paint=>draw_graph, :flex=>[0,0,0,0,1,1]),
	#	item(:button, :title=>"Update", :action=>doc_modified, :align=>:center, :flex=>[1,1,1,0,0,0]),
	#	item(:button, :title=>"Close", :action=>proc { hide }, :align=>:center, :flex=>[1,1,1,0,0,0]),
		:flex=>[0,0,0,0,1,1]
	  )
	  graph_item = item_with_tag("graph")
	  size = self.size
	  set_min_size(size)
	  set_size(500, 300)
	  @on_document_modified.call(mol)
	  show
	}
	self
  end
  
  def cmd_create_surface
    if @surface_dialog_attr == nil
	  @surface_dialog_attr = Hash.new
	  @surface_dialog_attr["hidden"] = -1
	end
    surface_dialog_attr = @surface_dialog_attr
    mol = self
	mol.open_auxiliary_window("MO Surface", :resizable=>true, :has_close_box=>true) {
	  tags = ["mo_ao", "mo", "color", "opacity", "threshold", "expand", "grid", "reverse", "basis"]
	  motype = mol.get_mo_info(:type)
	  alpha = mol.get_mo_info(:alpha)
	  beta = mol.get_mo_info(:beta)
	  ncomps = mol.get_mo_info(:ncomps)
	  mo_index = 0
	  mo_ao = nil
	  coltable = [[0,0,1], [1,0,0], [0,1,0], [1,1,0], [0,1,1], [1,0,1], [0,0,0], [1,1,1]]
	  mo_ao_items = ["Molecular Orbitals", "Atomic Orbitals"]
	  mo_ao_keys = ["MO", "AO"]
    if (nbo = mol.instance_eval { @nbo }) != nil
	    nbo.keys.each { |key|
        if key[0..2] == "AO/"
          key2 = key[3..-1]
          name2 = key2
          if key2 == "NAO"
            name2 = "Natural Atomic Orbitals"
          elsif key2 == "NBO"
            name2 = "Natural Bond Orbitals"
          elsif key2 == "NHO"
            name2 = "Natural Hybrid Orbitals"
          elsif key2 == "NLMO"
            name2 = "Natural Localized Molecular Orbitals"
          elsif key2 == "PNAO"
            name2 = "Pre-orthogonal Natural Atomic Orbitals"
          elsif key2 == "PNHO"
            name2 = "Pre-orthogonal Natural Hybrid Orbitals"
          elsif key2 == "PNBO"
            name2 = "Pre-orthogonal Natural Bond Orbitals"
          elsif key2 == "AHO"
            name2 = "Atomic Hybrid Orbitals"
          elsif key2 == "LHO"
            name2 = "Lewis Hybrid Orbitals"
          elsif key2 == "PLHO"
            name2 = "Pre-orthogonal Lewis Hybrid Orbitals"
          elsif key2 == "LPO"
            name2 = "Localized Property-optimized Orbitals"
          elsif key2 == "CLPO"
            name2 = "Chemist's Localized Property-optimized Orbitals"
          end
          mo_ao_items.push(name2)
          mo_ao_keys.push(key2)
        end
      }
	  end
	  mult = (motype == "UHF" ? 2 : 1)
	  mo_menu = ["1000 (-00.00000000)"]  #  Dummy entry; create later
	  tabvals = []
    coeffs_matrix = nil  #  Coefficient matrix for showing surface
    coeffs = nil   #  Coefficient array for showing surface
    table_coeffs_matrix = nil  #  Coefficient matrix for table (based on "basis")
	  a_idx_old = -1
	  occ = nil
	  on_update_attr = lambda {
	    tags.each { |tag|
		  surface_dialog_attr[tag] = item_with_tag(tag)[:value]
		}
		it = item_with_tag("hide")
		case surface_dialog_attr["hidden"]
		when -1
		  it[:title] = "Hide"
		  it[:enabled] = false
		when 0
		  it[:title] = "Hide"
		  it[:enabled] = true
		when 1
		  it[:title] = "Show"
		  it[:enabled] = true
		end
	  }
    update_mo_labels = lambda { |mo_ao_index|
      m = []
      if mo_ao_index == 0
        m = (1..(ncomps * mult)).map { |n|
          if motype == "UHF"
            i1 = (n - 1) / 2 + 1
            i2 = n % 2
            c1 = (i2 == 0 ? "B" : "A")
            c2 = (i1 > (i2 == 0 ? beta : alpha) ? "*" : "")
            else
            i1 = n
            i2 = 1
            c1 = ""
            if i1 > beta && i1 > alpha
              c2 = "*"
              elsif i1 > beta || i1 > alpha
              c2 = "S"
              else
              c2 = ""
            end
          end
          en = mol.get_mo_energy(i1 + (i2 == 0 ? ncomps : 0))
          sprintf("%d%s%s (%.8f)", i1, c1, c2, en)
        }
      elsif mo_ao_index == 1
        m = []
        ncomps.times { |i|
          m[i] = sprintf("%d: %s (%s)", i + 1, tabvals[i][2], mol.atoms[tabvals[i][3]].name)
        }
      else
        m = []
        key = mo_ao_keys[mo_ao_index]
        labels = nbo[key + "_L"]
        if labels == nil && key[0] == ?P
          labels = nbo[key[1..-1] + "_L"]
        end
        ncomps.times { |i|
          if labels
            lab = sprintf("%d: %s", i + 1, labels[i])
          else
            lab = sprintf("%s%d", key, i + 1)
          end
          m[i] = lab
        }
      end
      return m
    }
	  on_get_value = lambda { |it, row, col|
    begin
      if !table_coeffs_matrix || !coeffs_matrix
        #  Update tabvals
        tabvals = []
        bs = value("basis")
        if bs == 0  #  AO
          ncomps.times { |i|
            a_idx, s_idx, label = mol.get_gaussian_component_info(i)
            if a_idx_old != a_idx
              a_idx_old = a_idx
              a = a_idx.to_s
              n = mol.atoms[a_idx].name
            else
              a = n = ""
            end
            tabvals.push([a, n, label, a_idx])
          }
        else
          labels = update_mo_labels.call(bs + 1)
          ncomps.times { |i|
            label = labels[i]
            label = label.gsub(/^\d+: /, "")
            tabvals.push([(i + 1).to_s, "", label, i])
          }
        end
        #  Update coeffs_matrix
        if !coeffs_matrix
          #  coeffs_matrix = AO/(orbitals_to_display)
          if mo_ao == 0
            m = []
            ncomps.times { |i|
              mult.times { |j|
                m.push(mol.get_mo_coefficients(i * mult + j + 1))
              }
            }
            coeffs_matrix = LAMatrix.new(m)  #  Matrix AO/MO
          elsif mo_ao == 1
            coeffs_matrix = LAMatrix.identity(ncomps)  #  Matrix AO/AO  (identity)
          else
            coeffs_matrix = nbo["AO/" + mo_ao_keys[mo_ao]]
          end
        end
        #  m2 = AO/(basis)
        if bs == 0   #  AO
          m2 = LAMatrix.identity(ncomps)
        else
          m2 = nbo["AO/" + mo_ao_keys[bs + 1]]
        end
        #  (basis)/(orbitals_to_display)
        table_coeffs_matrix = m2.inverse * coeffs_matrix
      end
      #  Update coefficient array
      if !coeffs
        coeffs = coeffs_matrix.column(mo_index).to_a[0]
      end
      if col < 3
        tabvals[row][col]
      else
        sprintf("%.6f", table_coeffs_matrix[mo_index, row])
      end
    rescue => e
      $stderr.write(e.to_s + "\n")
      $stderr.write(e.backtrace.inspect + "\n")
    end
	  }
	  h = {"mo"=>nil, "color"=>nil, "opacity"=>nil, "threshold"=>nil, "expand"=>nil, "grid"=>nil, "reverse"=>nil, "basis"=>nil }
	  should_update = true
	  on_action = lambda { |it|
	    tag = it[:tag]
		value = it[:value]
		if tag == "color" || tag == "opacity"
		  opac = value("opacity").to_f
		  opac = 0.0 if opac < 0.0
		  opac = 1.0 if opac > 1.0
		  col = value("color")
		  color = coltable[col] + [opac]
		  color0 = [1,1,1,opac]
		  mol.set_surface_attr(:color=>color, :color0=>color0)
		  h[tag] = value
    elsif tag == "threshold" || tag == "reverse"
      thres = item_with_tag("threshold")[:value].to_f
		  thres = 0.001 if thres >= 0.0 && thres < 0.001
		  thres = -0.001 if thres <= 0.0 && thres > -0.001
      if item_with_tag("reverse")[:value] == 1
        thres = -thres
      end
		  mol.set_surface_attr(:thres=>thres)
		  h[tag] = value
		else
	      should_update = false
	      h.each_key { |key|
		    val = value(key)
		    if val && h[key] != val
		      should_update = true
		      break
		    end
		  }
		  item_with_tag("update")[:enabled] = should_update
		end
		on_update_attr.call
	  }
	  on_mo_action = lambda { |it|
	    mo = it[:value]
		if mo_ao == 0
		  if motype == "UHF"
		    mo_index = (mo / 2) + (mo % 2 == 1 ? ncomps : 0)
			if mo_index < alpha || (mo_index >= ncomps && mo_index < ncomps + beta)
			  occ_new = 1
			else
			  occ_new = 0
			end
		  else
		    mo_index = mo
			if mo_index < alpha && mo_index < beta
			  occ_new = 1
			elsif mo_index < alpha || mo_index < beta
			  occ_new = -1
			else
			  occ_new = 0
			end
		  end
		else
		  mo_index = mo
      occ_new = occ
      if mo_menu[mo]
        if mo_menu[mo] =~ /(\(ryd?)|(\(NB)|(\*)/
          occ_new = 0
        else
          occ_new = 1
        end
      end
		end
		coeffs = nil
		if occ_new != occ
		  #  Set default color
		  col = (occ_new == 0 ? 1 : (occ_new == -1 ? 2 : 0))
		  set_value("color", col)
		  h["color"] = col
		  occ = occ_new
		end
		item_with_tag("table")[:refresh] = true
	    on_action.call(it)
		on_update_attr.call
	  }
	  on_set_action = lambda { |it|
	    if mo_ao != it[:value]
        mo_ao = it[:value]
        mo_menu = update_mo_labels.call(mo_ao)
        it0 = item_with_tag("mo")
        val = it0[:value]
        it0[:subitems] = mo_menu
        it0[:value] = val   #  Keep the mo number invariant
        h["mo"] = nil  # "Update" button is forced to be enabled
        coeffs = nil
        coeffs_matrix = nil  #  Matrix needs update
        on_mo_action.call(it0)
      end
      on_update_attr.call
	  }
    on_basis_action = lambda { |it|
      if h["basis"] != value("basis")
        table_coeffs_matrix = nil
        item_with_tag("table")[:refresh] = true
      end
    }
	  on_update = lambda { |it|
	    h.each_key { |key|
		  h[key] = value(key)
		}
		opac = h["opacity"].to_f
		opac = 0.0 if opac < 0.0
		opac = 1.0 if opac > 1.0
		color = coltable[h["color"]] + [opac]
		color0 = [1, 1, 1, opac]
		thres = h["threshold"].to_f
		thres = 0.001 if thres >= 0.0 && thres < 0.001
		thres = -0.001 if thres <= 0.0 && thres > -0.001
    if h["reverse"] == 1
      thres = -thres
    end
		expand = h["expand"].to_f
		expand = 0.01 if expand < 0.01
		expand = 10.0 if expand > 10.0
		grid = h["grid"].to_i
		if grid > 10000000
		  grid = 10000000
		end
		if mo_ao == 0
		  idx = mo_index + 1
		else
		  idx = 0
		  mol.set_mo_coefficients(0, 0.0, coeffs)
		end
		if it[:tag] == "create_cube"
		  basename = File.basename(mol.path || mol.name, ".*")
	      fname1 = Dialog.save_panel("Cube file name", mol.dir || Dir.pwd, basename + ".cube", "Gaussian cube file (*.cube)|*.cube")
		  if fname1
		    mol.cubegen(fname1, idx, grid, true)
		  end
		else
		  mol.create_surface(idx, :npoints=>grid, :color=>color, :thres=>thres, :expand=>expand, :color0=>color0)
		  mol.show_surface
		  on_action.call(it)
		  surface_dialog_attr["hidden"] = 0
		  on_update_attr.call
		end
	  }
	  layout(1,
	    layout(2,
		  item(:text, :title=>"Orbital Set"),
		  item(:popup, :tag=>"mo_ao", :subitems=>mo_ao_items, :action=>on_set_action),
	      item(:text, :title=>"Select"),
	      item(:popup, :tag=>"mo", :subitems=>mo_menu, :action=>on_mo_action)),
		layout(4,
		  item(:text, :title=>"Color"),
		  item(:popup, :tag=>"color", :subitems=>["blue", "red", "green", "yellow", "cyan", "magenta", "black", "white"], :action=>on_action),
		  item(:text, :title=>"Opacity"),
		  item(:textfield, :tag=>"opacity", :width=>80, :value=>"0.8", :action=>on_action),
		  item(:text, :title=>"Threshold"),
		  item(:textfield, :tag=>"threshold", :width=>80, :value=>"0.05", :action=>on_action),
      item(:checkbox, :tag=>"reverse", :title=>"Reverse Phase", :action=>on_action),
      -1),
		layout(4,
		  item(:text, :title=>"Number of Grid Points"),
		  item(:textfield, :tag=>"grid", :width=>80, :value=>"512000", :action=>on_action),
      item(:text, :title=>"Box Limit"),
      item(:textfield, :tag=>"expand", :width=>80, :value=>"1.4", :action=>on_action)),
    #  table coefficients are based on "basis" (AO, NAO, etc.)
    #  MO is not allowed as basis, so the "basis" menu begins with "AO"
    layout(2,
           item(:text, :title=>"Coeffs based on:"),
           item(:popup, :tag=>"basis", :subitems=>mo_ao_keys[1..-1], :value=>1, :action=>on_basis_action)),
    item(:table, :width=>380, :height=>300, :tag=>"table",
		  :columns=>[["Atom", 40], ["Name", 60], ["Label", 140], ["Coeff", 120]],
		  :on_count=> lambda { |it| ncomps },
		  :on_get_value=>on_get_value,
		  :flex=>[0,0,0,0,1,1]),
		layout(3,
		  item(:button, :tag=>"update", :title=>"Update", :action=>on_update),
		  item(:button, :tag=>"clear", :title=>"Clear", :action=>lambda { |it|
		    mol.clear_surface
			item_with_tag("update")[:enabled] = true
			surface_dialog_attr["hidden"] = -1
			on_update_attr.call
		    } ),
		  item(:button, :tag=>"hide", :title=>"Hide", :action=>lambda { |it|
		    case surface_dialog_attr["hidden"]
			when 0
			  surface_dialog_attr["hidden"] = 1
			  mol.hide_surface
			when 1
			  surface_dialog_attr["hidden"] = 0
			  mol.show_surface
			end
			on_update_attr.call
			} ),
		  item(:button, :tag=>"create_cube", :title=>"Create Cube", :action=>on_update),
		  :flex=>[0,1,0,0,0,0]),
		:flex=>[0,0,0,0,1,1]
	  )
	  tags.each { |tag|
	    if (val = surface_dialog_attr[tag]) != nil
		  item_with_tag(tag)[:value] = val
		end
	  }
	  mo_idx = surface_dialog_attr["mo"]
	  on_update_attr.call
	  on_set_action.call(item_with_tag("mo_ao"))
	  item_with_tag("mo")[:value] = surface_dialog_attr["mo"] = mo_idx
	  size = self.size
	  set_min_size(size[0], 250)
	  item_with_tag("table")[:refresh] = true
	  show
    }
  end

  def ask_graphic_export_scale
    scale = get_global_settings("global.export_graphic_scale")
	scale = (scale ? scale.to_i - 1 : 3)
	bg_color = get_global_settings("global.export_background_color")
	bg_color = (bg_color ? bg_color.to_i + 1 : 0)
    hash = Dialog.run("Set Export Properties") {
	  layout(2,
		item(:text, :title=>"Resolution:"),
		item(:popup, :subitems=>["Screen", "Screenx2", "Screenx3", "Screenx4", "Screenx5"], :tag=>"resolution", :value=>scale),
		item(:text, :title=>"Background color:"),
		item(:popup, :subitems=>["Same as Screen", "Transparent", "Black", "White"], :tag=>"bg_color", :value=>bg_color))
	}
	if hash[:status] == 0
	  set_global_settings("global.export_graphic_scale", (hash["resolution"] + 1).to_s)
	  set_global_settings("global.export_background_color", (hash["bg_color"] - 1).to_s)
	  return 0
	else
	  return -1
	end
  end
  
  #  DEBUG
  def cmd_test
    $test_dialog = Dialog.new("Test") { item(:text, :title=>"test"); show }
  end
  
end

module Kernel
  def ask_scratch_dir
    sdir = get_global_settings("global.scratch_dir")
	while 1
      p = Dialog.open_panel("Please select scratch directory", sdir, nil, true)
	  if p
	    if p =~ / /
		  error_message_box("Please avoid path containing a white space.\n" + p.sub(/ /, "<!> <!>"))
		  sdir = p
		  next
		else
		  set_global_settings("global.scratch_dir", p)
		  return p
	    end
	  else
	    return nil
	  end
	end
  end
end

register_menu("Assign residue...", :cmd_assign_residue, :non_empty)
register_menu("Offset residue...", :cmd_offset_residue, :non_empty)
register_menu("Sort by residue", :cmd_sort_by_residue, :non_empty)
register_menu("", "")
register_menu("Delete Frames...", :cmd_delete_frames, lambda { |m| m && m.nframes > 1 } )
register_menu("Reverse Frames...", :cmd_reverse_frames, lambda { |m| m && m.nframes > 1 } )
register_menu("Concatenate Frames...", :cmd_concat_frames, lambda { |m| m && m.nframes > 1 } )
register_menu("", "")
register_menu("Show Energy Window...", :cmd_show_energy, lambda { |m| m && m.property_names.include?("energy") } )
register_menu("Show MO Surface...", :cmd_create_surface, lambda { |m| m && m.get_mo_info(:type) != nil } )
#register_menu("cmd test", :cmd_test)

