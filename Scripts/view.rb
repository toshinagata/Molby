#
#  view.rb
#
#  Created by Toshi Nagata.
#  Copyright 2014 Toshi Nagata. All rights reserved.
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

def cmd_show_unit_cell
  self.show_unitcell = !self.show_unitcell
end

def cmd_show_hydrogens
  self.show_hydrogens = !self.show_hydrogens
end

def cmd_show_dummy_atoms
  self.show_dummy_atoms = !self.show_dummy_atoms
end

def cmd_show_expanded
  self.show_expanded = !self.show_expanded
end

def cmd_show_ellipsoids
  self.show_ellipsoids = !self.show_ellipsoids
end

def cmd_show_rotation_center
  self.show_rotation_center = !self.show_rotation_center
end

def cmd_line_mode
  self.line_mode = !self.line_mode
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
  
end

register_menu("View\t-", nil)
register_menu("View\t^Show Unit Cell", :cmd_show_unit_cell,
  lambda { |m| [m != nil, m.show_unitcell] } )
register_menu("View\t^Show Hydrogen Atoms", :cmd_show_hydrogens,
  lambda { |m| [m != nil, m.show_hydrogens] } )
register_menu("View\t^Show Dummy Atoms", :cmd_show_dummy_atoms,
  lambda { |m| [m != nil, m.show_dummy_atoms] } )
register_menu("View\t^Show Expanded Atoms", :cmd_show_expanded,
  lambda { |m| [m != nil, m.show_expanded] } )
register_menu("View\t^Show Ellipsoids", :cmd_show_ellipsoids,
  lambda { |m| [m != nil, m.show_ellipsoids] } )
register_menu("View\t^Show Rotation Center", :cmd_show_rotation_center,
  lambda { |m| [m != nil, m.show_rotation_center] } )
register_menu("View\t-", nil)
register_menu("View\t^Show Graphite...", :cmd_show_graphite,
  lambda { |m| [m != nil, m.show_graphite?] } )
register_menu("View\t-", nil)
register_menu("View\t^Line Mode", :cmd_line_mode,
  lambda { |m| [m != nil, m.line_mode] } )
