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
    mol.open_auxiliary_window("Extra Props", nil, nil, :resizable=>true, :has_close_box=>true) {
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
	  on_document_modified.call(mol)
	  show
	}
  end
  
  def cmd_show_energy
	wave = [0.0, 0.0]
	cur = 0
	mol = self
	d = open_auxiliary_window("Energy", nil, nil, :resizable=>true, :has_close_box=>true) {
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
		min = wave.min
		h = wave.max - min
		h = (h == 0.0 ? 1.0 : height / h)
		w = wave.count
		w = (w == 0 ? 1.0 : Float(width) / w)
		a = []
		wave.each_with_index { |d, i|
		  a.push(i * w + 16)
		  a.push(height - (d - min) * h + 12)
		}
		if wave.count == 1
		  a.push(w + 16)
		  a.push(height - (wave[0] - min) * h + 12)
		end
		draw_line(a)
		brush(:color=>[0.2, 0.2, 1.0])
		y = wave[cur] || 0.0
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
    mol = self
	mol.open_auxiliary_window("MO Surface", nil, nil, :resizable=>true, :has_close_box=>true) {
	  motype = mol.get_mo_info(:type)
	  alpha = mol.get_mo_info(:alpha)
	  beta = mol.get_mo_info(:beta)
	  ncomps = mol.get_mo_info(:ncomps)
	  mo_index = 1
	  mo_ao = nil
	  coltable = [[0,0,1], [1,0,0], [0,1,0], [1,1,0], [0,1,1], [1,0,1], [0,0,0]]
	  if (motype == "RHF")
	    beta = nil
	  end
	  i = (beta ? 2 : 1)
	  mo_menu = []  #  Create later
	  tabvals = []
	  coeffs = nil
	  a_idx_old = -1
	  ncomps.times { |i|
	    a_idx, label, nprims = mol.get_gaussian_shell_info(i)
		if a_idx_old != a_idx
		  a_idx_old = a_idx
		  a = a_idx.to_s
		  n = mol.atoms[a_idx].name
		else
		  a = n = ""
		end
		tabvals.push([a, n, label, a_idx])
	  }
	  on_get_value = lambda { |it, row, col|
	    if col < 3
		  tabvals[row][col]
		else
		  if coeffs == nil
		    if mo_ao == 0
		      coeffs = mol.get_mo_coefficients(mo_index)
		    else
		      coeffs = (0...ncomps).map { |i| (i == mo_index ? 1.0 : 0.0) }
			end
		  end
		  sprintf("%.6f", coeffs[row])
		end
	  }
	  h = {"mo"=>nil, "color"=>nil, "opacity"=>nil, "threshold"=>nil, "expand"=>nil, "grid"=>nil}
	  should_update = true
	  on_action = lambda { |it|
	    should_update = false
	    h.each_key { |key|
		  val = value(key)
		  if val && h[key] != val
		    should_update = true
		    break
		  end
		}
		item_with_tag("update")[:enabled] = should_update
	  }
	  on_mo_action = lambda { |it|
	    mo = it[:value]
		if mo_ao == 0
		  if beta
		    mo_index = (mo / 2) + (mo % 2 == 1 ? ncomps : 0) + 1
		  else
		    mo_index = mo + 1
		  end
		else
		  mo_index = mo
		end
		coeffs = nil
		item_with_tag("table")[:refresh] = true
	    on_action.call(it)
	  }
	  on_set_action = lambda { |it|
	    if mo_ao != it[:value]
		  mo_ao = it[:value]
		  if mo_ao == 0
		    mo_menu = (1..(ncomps * i)).map { |n|
	          if beta
		        i1 = (n - 1) / 2 + 1
		        i2 = n % 2
		        c1 = (i2 == 0 ? "B" : "A")
		        c2 = (i1 > (i2 == 0 ? beta : alpha) ? "*" : "")
		      else
		        i1 = n
		        i2 = 1
		        c1 = ""
		        c2 = (i1 >= alpha ? "*" : "")
		      end
		      en = mol.get_mo_energy(i1 + (i2 == 0 ? ncomps : 0))
		      sprintf("%d%s%s (%.8f)", i1, c1, c2, en)
			}
		  else
		    mo_menu = []
		    ncomps.times { |i|
			  mo_menu[i] = sprintf("AO%d: %s (%s)", i + 1, tabvals[i][2], mol.atoms[tabvals[i][3]].name)
			}
		  end
		  it0 = item_with_tag("mo")
		  it0[:subitems] = mo_menu
		  it0[:value] = 0
		  on_mo_action.call(it0)
		end
	  }
	  on_update = lambda { |it|
	    h.each_key { |key|
		  h[key] = value(key)
		}
		opac = h["opacity"].to_f
		color = coltable[h["color"]] + [opac]
		thres = h["threshold"].to_f
		thres = 0.001 if thres >= 0.0 && thres < 0.001
		thres = -0.001 if thres <= 0.0 && thres > -0.001
		expand = h["expand"].to_f
		expand = 0.01 if expand < 0.01
		expand = 10.0 if expand > 10.0
		grid = h["grid"].to_i
		if grid > 10000000
		  grid = 10000000
		end
		if mo_ao == 0
		  idx = mo_index
		else
		  idx = 0
		  mol.set_mo_coefficients(0, 0.0, coeffs)
		end
		mol.create_surface(idx, :npoints=>grid, :color=>color, :thres=>thres, :expand=>expand)
		on_action.call(it)
	  }
	  layout(1,
	    layout(2,
		  item(:text, :title=>"Orbital Set"),
		  item(:popup, :tag=>"mo_ao", :subitems=>["Molecular Orbitals", "Atomic Orbitals"], :action=>on_set_action),
	      item(:text, :title=>"Select"),
	      item(:popup, :tag=>"mo", :subitems=>mo_menu, :action=>on_mo_action)),
		layout(4,
		  item(:text, :title=>"Color"),
		  item(:popup, :tag=>"color", :subitems=>["blue", "red", "green", "yellow", "cyan", "magenta", "black"], :action=>on_action),
		  item(:text, :title=>"Opacity"),
		  item(:textfield, :tag=>"opacity", :width=>80, :value=>"0.6", :action=>on_action),
		  item(:text, :title=>"Threshold"),
		  item(:textfield, :tag=>"threshold", :width=>80, :value=>"0.05", :action=>on_action),
		  item(:text, :title=>"Box Limit"),
		  item(:textfield, :tag=>"expand", :width=>80, :value=>"1.0", :action=>on_action)),
		layout(2,
		  item(:text, :title=>"Number of Grid Points"),
		  item(:textfield, :tag=>"grid", :width=>120, :value=>"512000", :action=>on_action)),
	    item(:table, :width=>300, :height=>300, :tag=>"table",
		  :columns=>[["Atom", 60], ["Name", 60], ["Label", 60], ["Coeff", 120]],
		  :on_count=> lambda { |it| tabvals.count },
		  :on_get_value=>on_get_value,
		  :flex=>[0,0,0,0,1,1]),
		item(:button, :tag=>"update", :title=>"Update", :action=>on_update, :flex=>[0,1,0,0,0,0]),
		:flex=>[0,0,0,0,1,1]
	  )
	  on_set_action.call(item_with_tag("mo_ao"))
	  size = self.size
	  set_min_size(size[0], 250)
	  item_with_tag("table")[:refresh] = true
	  show
    }
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
register_menu("", "")
register_menu("Show Energy Window...", :cmd_show_energy, lambda { |m| m && m.property_names.include?("energy") } )
register_menu("Show MO Surface...", :cmd_create_surface, lambda { |m| m && m.get_mo_info(:type) != nil } )
#register_menu("cmd test", :cmd_test)

