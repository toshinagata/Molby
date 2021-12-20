# coding: utf-8
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

def cmd_line_mode
  self.line_mode = !self.line_mode
end

def cmd_ball_and_stick_mode
  mol = self
  rad = mol.atom_radius
  res = mol.atom_resolution
  brad = mol.bond_radius
  bres = mol.bond_resolution
  if mol.line_mode == true || mol.atom_radius >= 0.5
    rad = 0.2
	res = 12
	brad = 0.1
	bres = 8
  end
  hash = Dialog.run("Ball and Stick") {
    layout(2,
	  item(:text, :title=>"Atom Radius (0..0.5, default 0.2)"),
	  item(:textfield, :width=>100, :tag=>"atom_radius", :value=>rad.to_s),
	  item(:text, :title=>"Atom Resolution (default 12)"),
	  item(:textfield, :width=>100, :tag=>"atom_resolution", :value=>res.to_s),
	  item(:text, :title=>"Bond Radius (default 0.1)"),
	  item(:textfield, :width=>100, :tag=>"bond_radius", :value=>brad.to_s),
	  item(:text, :title=>"Atom Resolution (default 8)"),
	  item(:textfield, :width=>100, :tag=>"bond_resolution", :value=>bres.to_s))
  }
  if hash[:status] == 0
    f = hash["atom_radius"].to_f
	self.atom_radius = f if f > 0.0 && f < 0.5
    n = hash["atom_resolution"].to_i
	self.atom_resolution = n if n >= 6 && n < 100
    f = hash["bond_radius"].to_f
	self.bond_radius = f if f > 0.0
    n = hash["bond_resolution"].to_i
	self.bond_resolution = n if n >= 4 && n < 100
	self.line_mode = false
  end
end

def cmd_space_filling_mode
  mol = self
  rad = mol.atom_radius
  res = mol.atom_resolution
  if mol.line_mode == true || mol.atom_radius < 0.5
    rad = 1.0
	res = 18
  end
  hash = Dialog.run("Space Filling") {
    layout(2,
	  item(:text, :title=>"Atom Radius (0.5..2.0, default 1.0)"),
	  item(:textfield, :width=>100, :tag=>"atom_radius", :value=>rad.to_s),
	  item(:text, :title=>"Atom Resolution (default 18)"),
	  item(:textfield, :width=>100, :tag=>"atom_resolution", :value=>res.to_s))
  }
  if hash[:status] == 0
    f = hash["atom_radius"].to_f
	self.atom_radius = f if f >= 0.5 && f <= 2.0
    n = hash["atom_resolution"].to_i
	self.atom_resolution = n if n >= 6 && n < 100
	self.line_mode = false
  end
end

end

register_menu("View\t-", nil)
register_menu("View\t^Show Unit Cell", :cmd_show_unit_cell,
  lambda { |m| [m != nil, m && m.show_unitcell] } )
register_menu("View\t^Show Hydrogen Atoms", :cmd_show_hydrogens,
  lambda { |m| [m != nil, m && m.show_hydrogens] } )
register_menu("View\t^Show Dummy Atoms", :cmd_show_dummy_atoms,
  lambda { |m| [m != nil, m && m.show_dummy_atoms] } )
register_menu("View\t^Show Expanded Atoms", :cmd_show_expanded,
  lambda { |m| [m != nil, m && m.show_expanded] } )
register_menu("View\t^Show Ellipsoids", :cmd_show_ellipsoids,
  lambda { |m| [m != nil, m && m.show_ellipsoids] } )
register_menu("View\t^Show Rotation Center", :cmd_show_rotation_center,
  lambda { |m| [m != nil, m && m.show_rotation_center] } )
register_menu("View\t-", nil)
register_menu("View\t^Show Graphite...", :cmd_show_graphite,
  lambda { |m| [m != nil, m && m.show_graphite?] } )
register_menu("View\t-", nil)
register_menu("View\tAppearance\t^Line", :cmd_line_mode,
  lambda { |m| [m != nil, m && m.line_mode] } )
register_menu("View\tAppearance\t^Ball and Stick", :cmd_ball_and_stick_mode,
  lambda { |m| [m != nil, m && !m.line_mode && m.atom_radius < 0.5] } )
register_menu("View\tAppearance\t^Space Filling", :cmd_space_filling_mode,
  lambda { |m| [m != nil, m && !m.line_mode && m.atom_radius >= 0.5] } )
