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
    hash = Dialog.run {
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
    hash = Dialog.run {
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

  def cmd_create_cube
    grid = default_MO_grid
	if grid == nil
	  Dialog.run {
	    layout(1,
		  item(:text, :title=>"This molecule does not contain MO information."))
	  }
	  return
	end
    mos = selected_MO
	if mos == nil || mos.length == 0
      Dialog.run {
	    layout(1,
		  item(:text, :title=>"Please select MO(s) in the MO Info table."))
      }
	  return
	end
	hash = Dialog.run {
	  layout(1,
	    item(:text, :title=>"Please specify cube dimensions (in bohr unit):"),
	    layout(4,
		  item(:text, :title=>"Origin"),
		  item(:textfield, :width=>100, :height=>20, :tag=>"originx", :value=>sprintf("%.6f", grid[0].x)),
		  item(:textfield, :width=>100, :height=>20, :tag=>"originy", :value=>sprintf("%.6f", grid[0].y)),
		  item(:textfield, :width=>100, :height=>20, :tag=>"originz", :value=>sprintf("%.6f", grid[0].z)),
		  item(:text, :title=>"Delta"),
		  item(:textfield, :width=>100, :height=>20, :tag=>"deltax", :value=>sprintf("%.6f", grid[1])),
		  item(:textfield, :width=>100, :height=>20, :tag=>"deltay", :value=>sprintf("%.6f", grid[2])),
		  item(:textfield, :width=>100, :height=>20, :tag=>"deltaz", :value=>sprintf("%.6f", grid[3])),
		  item(:text, :title=>"Step"),
		  item(:textfield, :width=>100, :height=>20, :tag=>"stepx", :value=>grid[4].to_s),
		  item(:textfield, :width=>100, :height=>20, :tag=>"stepy", :value=>grid[5].to_s),
		  item(:textfield, :width=>100, :height=>20, :tag=>"stepz", :value=>grid[6].to_s)))
	}
	if hash[:status] == 0
	  path = self.path || self.name
	  dir = self.dir || Dir.pwd
	  origin = Vector3D[hash["originx"], hash["originy"], hash["originz"]]
	  dx = hash["deltax"]
	  dy = hash["deltay"]
	  dz = hash["deltaz"]
	  nx = hash["stepx"]
	  ny = hash["stepy"]
	  nz = hash["stepz"]
	  basename = File.basename(path, ".*")
	  filenames = []
	  mo_type = self.mo_type
	  mos.each { |n|
	    fname1 = fname2 = nil
	    alpha = (mo_type != "UHF" ? "" : "alpha ")
		a = (mo_type != "UHF" ? "" : "a")
	    fname1 = Dialog.save_panel("Cube file name for #{alpha}MO #{n}", dir, basename + "_#{n}#{a}.cube", "Gaussian cube file (*.cube)|*.cube")
		if (mo_type == "UHF")
		  fname2 = Dialog.save_panel("Cube file name for beta MO #{n}", dir, basename + "_#{n}b.cube", "Gaussian cube file (*.cube)|*.cube")
		end
		filenames.push([n, fname1, fname2])
	  }
	  filenames.each { |pair|
	    n = pair[0]
		alpha = (mo_type != "UHF" ? "" : "alpha ")
	    show_progress_panel("Creating cube file for #{alpha}MO #{n}...")
		if pair[1]
		  cubegen(pair[1], n, origin, dx, dy, dz, nx, ny, nz, true)
		end
		if pair[2] && mo_type == "UHF"
		  set_progress_message("Creating cube file for beta MO #{n}...")
		  cubegen(pair[2], n, origin, dx, dy, dz, nx, ny, nz, true, true)
		end
	    hide_progress_panel
 	  }
	end
  end

  def cmd_delete_frames
    n = nframes
    return if n == 0
	hash = Dialog.run {
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

  def cmd_show_graphite
    n = self.show_graphite
	flag = self.show_graphite?
	hash = Dialog.run("Show Graphite") {
      layout(1,
	    item(:checkbox, :title=>"Show graphite", :tag=>"show_graphite", :value=>(flag ? 1 : 0),
		  :action=>lambda { |it| set_attr("graphite", :enabled=>(it[:value] == 1)) } ),
	    item(:text, :title=>"Number of graphite rings for each direction:"),
	    item(:textfield, :width=>120, :tag=>"graphite", :value=>n.to_s, :enabled=>flag))
	}
	if hash[:status] == 0
	  self.show_graphite(hash["graphite"])
	  self.show_graphite(hash["show_graphite"] == 1 ? true : false)
	end
  end
  
  #  DEBUG
  def cmd_test
    $test_dialog = Dialog.new("Test") { item(:text, :title=>"test"); show }
  end
  
end

register_menu("Assign residue...", :cmd_assign_residue)
register_menu("Offset residue...", :cmd_offset_residue)
register_menu("Sort by residue", :cmd_sort_by_residue)
register_menu("", "")
register_menu("Delete Frames...", :cmd_delete_frames)
#register_menu("cmd test", :cmd_test)

register_menu("QChem\tCreate GAMESS Input...",
  :cmd_create_gamess_input, :non_empty)   # gamess.rb
register_menu("QChem\tCreate MOPAC6 Input...",
  :cmd_create_mopac_input, :non_empty)    # mopac6.rb
register_menu("QChem\tCreate MO Cube...",
  :cmd_create_cube, lambda { |m| m && m.mo_type } )
