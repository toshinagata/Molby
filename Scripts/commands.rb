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

  def save_gamess_with_ecp(filename)
    if natoms == 0
      raise MolbyError, "cannot save GAMESS input; the molecule is empty"
    end
    fp = open(filename, "wb")
	now = Time.now.to_s
	fp.print <<end_of_header
!  GAMESS input
!  Generated by Molby at #{now}
 $CONTRL COORD=UNIQUE EXETYP=RUN ICHARG=0
         ICUT=20 INTTYP=HONDO ITOL=30
         MAXIT=200 MOLPLT=.T. MPLEVL=0
         MULT=1 QMTTOL=1e-08 RUNTYP=OPTIMIZE
		 ECP=READ 
         SCFTYP=RHF UNITS=ANGS                  $END
 $SCF    CONV=1.0E-06 DIRSCF=.T. FDIFF=.T. DAMP=.T. $END
 $STATPT NSTEP=400 OPTTOL=1.0E-06               $END
 $SYSTEM MEMDDI=0 MWORDS=16 TIMLIM=50000        $END
 $BASIS  EXTFIL=.T. GBASIS=321LAN               $END
 $GUESS  GUESS=HUCKEL                           $END
!
 $DATA
 #{name}
 C1
end_of_header
	ecp = Hash.new
	each_atom { |ap|
		fp.printf " %-6s %4d %10.6f %10.6f %10.6f\n", ap.name, ap.atomic_number, ap.r.x, ap.r.y, ap.r.z
		if ap.atomic_number >= 11
			ecp[ap.element.upcase] = 1
		end
	}
	fp.print " $END\n"
	#  Read ECP from file
	efp = open($startup_dir + "/gamess_ecp.txt", "rb")
	elem = ""
	ecplines = nil
	while 1
		line = efp.gets
		if line == nil || line =~ /(\w\w)-ECP /
			if ecplines != nil
				print "#{line}"
				ecp[elem] = ecplines
			end
			break if line == nil
			elem = $1.upcase
			if ecp.member?(elem)
				ecplines = ""
			else
				ecplines = nil
			end
			print "#{line}: #{'->ignore' if ecplines == nil}\n"
		end
		if ecplines != nil
			ecplines += line
		end
	end
	efp.close
	#  Create ECP section
	fp.print " $ECP\n"
	each_atom { |ap|
		elem = ap.element.upcase
		if ecp.member?(elem)
			fp.print ecp[elem]
			ecp[elem] = "#{elem}-ECP\n"
		else
			fp.print "#{elem}-ECP NONE\n"
		end
	}
	fp.print " $END\n"	
	fp.close
  end
  
  def cmd_create_gamess
    hash = Dialog.run {
	  layout(1,
	    item(:text, :title=>"Create GAMESS input with LANL2DZ effective core potentials"))
    }
    if hash[:status] == 0
      fname = File.basename(self.path, ".*") + ".inp"
	  fname = Dialog.save_panel(nil, self.dir, fname)
	  if fname
		save_gamess_with_ecp(fname)
	  end
	end
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
	  origin = Vector3D[hash["originx"], hash["originy"], hash["originz"]]
	  dx = hash["deltax"]
	  dy = hash["deltay"]
	  dz = hash["deltaz"]
	  nx = hash["stepx"]
	  ny = hash["stepy"]
	  nz = hash["stepz"]
	  basename = File.basename(self.path, ".*")
	  filenames = []
	  mos.each { |n|
	    fname = Dialog.save_panel("Cube file name for MO #{n}", self.dir, basename + "_#{n}.cube", "Gaussian cube file (*.cube)|*.cube")
		if !fname
		  filenames.clear
		  break
		end
		filenames.push([n, fname])
	  }
	  filenames.each { |pair|
	    n = pair[0]
	    show_progress_panel("Creating cube file for MO #{n}...")
		cubegen(pair[1], n, origin, dx, dy, dz, nx, ny, nz, true)
	    hide_progress_panel
 	  }
	end
  end

  def Dialog.list_remote_files(host, directory)
#    list = `ssh #{host} "env COLUMNS=40 ls -C #{directory}"`
    list = `ssh #{host} "ls #{directory}"`
  end
  
  def Molecule.cmd_load_remote(mol)  #  mol is not used
    hash = Dialog.run {
	  def button_action(it)   #  Action for OK and Cancel buttons
		if it[:index] == 0  #  OK
		  local = File.expand_path(value("local"))
		  if value("local") != ""
		    #  Check whether the local file already exists
			exist = []
			sfile = value("sfile")
			cfile = value("cfile")
			if sfile != "" && FileTest.exist?("#{local}/#{sfile}")
			  exist.push(sfile)
			end
			if cfile != "" && FileTest.exist?("#{local}/#{cfile}")
			  exist.push(cfile)
			end
			if exist.length > 0
			  if exist.length == 1
			    msg = "The file #{exist[0]} already exists"
			  else
			    msg = "The files " + exist.join(", ") + " already exist"
			  end
			  msg += " in directory #{local}. Overwrite?"
			  hash = Dialog.run {
			    layout(1, item(:text, :title=>msg, :width=>240, :height=>60))
			  }
			  return if hash[:status] != 0  #  No call of super -> dialog is not dismissed
			end
		  end
	    end
		end_modal(it)
	  end
	  def text_action(it)
		if value("host") != "" && value("directory") != "" && value("sfile") != ""
		  set_attr(0, :enabled=>true)
		else
		  set_attr(0, :enabled=>false)
		end	  
	  end
	  layout(2,
		item(:text, :title=>"Remote host"),
		item(:textfield, :width=>280, :height=>20, :tag=>"host", :action=>:text_action),
		item(:text, :title=>"Directory"),
		item(:textfield, :width=>280, :height=>20, :tag=>"directory", :action=>:text_action),
		item(:text, :title=>"Structure File"),
		item(:textfield, :width=>280, :height=>20, :tag=>"sfile", :action=>:text_action),
		item(:text, :title=>"Coordinate File"),
		item(:textfield, :width=>280, :height=>20, :tag=>"cfile"),
		item(:text, :title=>"File List"),
		item(:textview, :width=>280, :height=>80, :tag=>"list", :editable=>false),
		nil,
		[ item(:button, :title=>"Update",
			:action=>proc { |it| 
			  list = Dialog.list_remote_files(value("host"), value("directory"))
			  set_value("list", list)
			}
		  ), {:align=>:right} ],
	#	item(:checkbox, :title=>"Copy files to a local directory"),
	#	nil,
		item(:text, :title=>"Local directory"),
		item(:textfield, :width=>280, :height=>20, :tag=>"local"),
		nil,
		[ item(:button, :title=>"Choose...",
		    :action=>proc { |it|
			  dir = Dialog.open_panel(nil, nil, nil, true)
			  if dir
			    set_value("local", dir)
			  end
		    }
		  ), {:align=>:right} ]
	#	layout(1, 2, 1, 0)
	  )
	  set_attr(0, :action=>:button_action)
	  set_attr(1, :action=>:button_action)
	  set_attr(0, :enabled=>false)
	  self.each_item { |it|
		tag = it[:tag]
	    if (type = it[:type]) == :textfield || type == :textview
		  val = get_global_settings("load_remote.#{tag}")
		  if (val != nil)
			set_value(tag, val)
		  end
		end
	  }
    }
	hash.each_pair { |key, value|
	  next if key == :status
	  set_global_settings("load_remote.#{key}", value)
	}
	if hash[:status] == 0
	  sfile = hash["sfile"]
	  cfile = hash["cfile"]
	  host = hash["host"]
	  directory = hash["directory"]
	  local = hash["local"]
	  if local == ""
	    local = `mktemp -d /tmp/MolbyRemoteLoad.XXXXXX`.chomp
		if $? != 0
		  raise "Cannot create a temporary directory"
		end
	  end
	  local = File.expand_path(local)
	  if cfile == "" && sfile =~ /\.psf$/
	    cfile = sfile.sub(/\.psf$/, ".pdb")
	  end
	  if cfile == ""
	    files = sfile
	  else
	    files = "{#{sfile},#{cfile}}"
	  end
	  show_progress_panel("Fetching remote file(s)...")
	  if !system("scp '#{host}:#{directory}/#{files}' '#{local}'")
	    raise "Cannot copy remote files"
	  end
	  hide_progress_panel()
	  mol = Molecule.open("#{local}/#{sfile}")
	  if cfile != "" && FileTest.exist?("#{local}/#{cfile}")
	    mol.undo_enabled = false
	    mol.molload("#{local}/#{cfile}")
		mol.undo_enabled = true
	  end
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

  def cmd_solvate
    #  Find molecule with unit cell defined
    solv = []
	solvnames = []
	Molecule.list.each { |m|
	  if m.cell != nil
	    solv.push m
		solvnames.push m.name
	  end
    }
	if solvnames.length == 0
	  Dialog.run {
	    layout(1,
		  item(:text, :title=>"Please open a molecule file containing a solvent box."))
	  }
	  return
	end
	hash = Dialog.run {
	  layout(1,
	    item(:text, :title=>"Choose solvent box:"),
	    item(:popup, :subitems=>solvnames, :tag=>"solvent"),
		item(:text, :title=>"Box offset\n(Negative numbers for absolute sizes)"),
		layout(2,
		  item(:text, :title=>"x"),
		  item(:textfield, :width=>"120", :tag=>"x", :value=>"10.0"),
		  item(:text, :title=>"y"),
		  item(:textfield, :width=>"120", :tag=>"y", :value=>"10.0"),
		  item(:text, :title=>"z"),
		  item(:textfield, :width=>"120", :tag=>"z", :value=>"10.0")),
		item(:text, :title=>"Exclusion limit distance:"),
		item(:textfield, :width=>"120", :tag=>"limit", :value=>"3.0"))
	}
	if hash[:status] == 0
	  solvate(solv[hash["solvent"]], [hash["x"], hash["y"], hash["z"]], hash["limit"])
	end
  end
     
  def cmd_show_graphite
    n = self.show_graphite
	hash = Dialog.run("Show Graphite") {
      layout(1,
	    item(:text, :title=>"Number of graphite rings for each direction:\n(0 to suppress display)"),
	    item(:textfield, :width=>120, :tag=>"graphite", :value=>n.to_s))
	}
	if hash[:status] == 0
	  self.show_graphite(hash["graphite"])
	end
  end
  
end

register_menu("Assign residue...", :cmd_assign_residue)
register_menu("Offset residue...", :cmd_offset_residue)
register_menu("Sort by residue", :cmd_sort_by_residue)
register_menu("", "")
register_menu("Load remote...", :cmd_load_remote)
#register_menu("Create GAMESS input...", :cmd_create_gamess)
#register_menu("Create Cube file...", :cmd_create_cube)
register_menu("", "")
register_menu("Delete Frames...", :cmd_delete_frames)
register_menu("Solvate...", :cmd_solvate)
