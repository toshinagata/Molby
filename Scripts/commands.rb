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

  def cmd_extra_properties
    mol = self
	on_document_modified = lambda { |m|
	}
	get_count = lambda { |it| mol.nframes }
	get_value = lambda { |it, row, col| sprintf("%.8g", mol.get_property(col, row)) rescue "" }
    Dialog.new("Extra Props:" + mol.name, nil, nil, :resizable=>true, :has_close_box=>true) {
	  columns = mol.property_names.map { |name| [name, 80] }
	  layout(1,
		item(:table, :width=>240, :height=>380, :flex=>[0,0,0,0,1,1], :tag=>"table",
		  :columns=>columns,
		  :on_count=>get_count,
		  :on_get_value=>get_value),
	    :flex=>[0,0,0,0,1,1]
	  )
	  set_min_size(480, 200)
	  listen(mol, "documentModified", on_document_modified)
	  listen(mol, "documentWillClose", lambda { |m| close } )
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
register_menu("", "")
register_menu("Open Extra Properties Window...", :cmd_extra_properties, lambda { |m| m && m.property_names.count > 0 } )
#register_menu("cmd test", :cmd_test)

