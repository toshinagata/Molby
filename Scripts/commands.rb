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
    Dialog.new("Extra Props:" + mol.name, nil, nil, :resizable=>true, :has_close_box=>true) {
	  names = nil
	  on_document_modified = lambda { |m|
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
	  listen(mol, "documentModified", on_document_modified)
	  listen(mol, "documentWillClose", lambda { |m| close } )
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
		item(:view, :frame=>[0, 0, 100, 80], :tag=>"graph", :on_paint=>draw_graph, :flex=>[0,0,0,0,1,1]),
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
#register_menu("cmd test", :cmd_test)

