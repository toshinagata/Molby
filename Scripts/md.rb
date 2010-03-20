#
#  md.rb
#
#  Created by Toshi Nagata.
#  Copyright 2009 Toshi Nagata. All rights reserved.
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

  def cmd_md_sub
    arena = self.md_arena
    keys = arena.to_hash.keys
	#  The read-only keys
	read_only = [:step, :coord_frame, :transient_temperature, :average_temperature]
	#  Sort the keys so that they are (a little more) readable
	[:timestep, :temperature, :cutoff, :electro_cutoff, :pairlist_distance,
	 :scale14_vdw, :scale14_elect, :use_xplor_shift, :dielectric,
	 :andersen_freq, :andersen_coupling, :random_seed, :relocate_center,
	 :use_graphite, :surface_probe_radius, :surface_tension, :surface_potential_freq,
	 :gradient_convergence, :coordinate_convergence,
	 :log_file, :coord_file, :vel_file, :force_file, :debug_file, :debug_output_level,
	 :coord_output_freq, :energy_output_freq,
	 :step, :coord_frame, :transient_temperature, :average_temperature].each_with_index { |k, i|
	  keys.delete(k)
	  keys.insert(i, k)
	}
	#  Arrange the items in vertical direction
	n = (keys.count + 1) / 2
	i = 0
	keys = keys.sort_by { |k| i += 1; (i > n ? i - n + 0.5 : i) }
	#  Do dialog
    hash = Dialog.run("Molecular Dynamics Advanced Settings") {
	  items = []
	  keys.each { |k|
	    enabled = !read_only.include?(k)
	    items.push(item(:text, :title=>k.to_s))
		it = item(:textfield, :width=>120, :value=>arena[k].to_s)
		if enabled
		  it[:tag] = k
		else
		  it[:enabled] = false
		end
		items.push(it)
	  }
	  layout(4, *items)
	}
	if hash[:status] == 0
	  hash.keys.each { |k|
	    next if k == :status || read_only.include?(k)
	    arena[k] = hash[k]
	  }
#	  arena.prepare
	  return hash
	else
	  return nil
	end
  end
  
  def prepare_arena
  
	#  Get the MDArena
	arena = self.md_arena
	if !arena.prepare(true)   #  Parameter check only
	  Dialog.run("MD Error") {
		layout(1,
		  item(:text, :title=>"Some parameters are missing. Please open\nthe 'parameters' table and examine."))
		  set_attr(1, :hidden=>true)  #  Hide the cancel button
	  }
	  return nil
    end

    #  Initialize some fields at first invocation
    if !@arena_created
	  if arena.log_file == nil && self.dir != nil
	    arena.log_file = self.dir + "/" + self.name.gsub(/\.\w+$/, ".log")
	  end
	  @md_spf = 500
	  @md_nf = 200
	  @arena_created = true
	end
	
    return arena

  end

  def cmd_md(minimize)

	#  Create arena
	arena = self.prepare_arena
    if !arena
	  return -1
	end
	
    #  Basic parameters
	spf = @md_spf
	nf = @md_nf
	mol = self
	files_save = Hash.new
	[:log_file, :coord_file, :vel_file, :force_file, :debug_file].each { |k|
	  s = arena[k]
	  files_save[k] = s
	  if s != nil
	    arena[k] = File.basename(s)
	  end
	}
	title = (minimize ? "Minimize" : "Molecular Dynamics")
	hash = Dialog.run(title) {
	  items = []
	  if !minimize
	    items.push item(:text, :title=>"Timestep (fs)")
	    items.push item(:textfield, :width=>120, :value=>arena.timestep.to_s, :tag=>"timestep")
	    items.push item(:text, :title=>"Target temperature (K)")
	    items.push item(:textfield, :width=>120, :value=>arena.temperature.to_s, :tag=>"temperature")
	  end
	  items.push item(:text, :title=>"Steps per frame")
	  items.push item(:textfield, :width=>120, :value=>arena.coord_output_freq.to_s, :tag=>"steps_per_frame")
	  items.push item(:text, :title=>"Number of frames")
	  items.push item(:textfield, :width=>120, :value=>nf.to_s, :tag=>"number_of_frames")
	  items.push item(:text, :title=>"Log file")
	  items.push item(:textfield, :width=>120, :value=>(arena.log_file || ""), :tag=>"log_file")
	  items.push item(:button, :title=>"Advanced...",
		     :action=>proc { |it|
			   if mol.cmd_md_sub
			     if !minimize
			       set_value("timestep", arena[:timestep].to_s)
				   set_value("temperature", arena[:temperature].to_s)
				 end
				 set_value("log_file", arena[:log_file])
			   end
			 }
			)
	  layout(2, *items)
	}
	dirstr = (self.dir || document_home)
	status = hash[:status]
	if status == 0
	  arena[:log_file] = hash[:log_file]
	end
	[:log_file, :coord_file, :vel_file, :force_file, :debug_file].each { |k|
	  if status == 0
	    s = arena[k]
	    if s != nil && s != ""
	      arena[k] = dirstr + "/" + s
		else
		  arena[k] = nil
		end
	  else
	    arena[k] = files_save[k]
      end
	}
	if status != 0
	  return -1
	end
	if !minimize
	  arena.temperature = Float(hash["temperature"])
	  arena.timestep = Float(hash["timestep"])
	end
	arena.coord_output_freq = Integer(hash["steps_per_frame"])
	arena.energy_output_freq = arena.coord_output_freq
	return Integer(hash["number_of_frames"])
  end
  
  def cmd_define_unit_cell
    mol = self
    hash = Dialog.run("Define Unit Cell") {
	  @mol = mol
	  def set_box_value(item1)
	    h = Hash.new
		["o0", "o1", "o2", "a0", "a1", "a2", "b0", "b1", "b2", "c0", "c1", "c2"].each { |k|
		  begin
		    s = value(k)
			if s == nil || s == ""
			  h[k] = 0.0
			else
		      h[k] = Float(eval(s))
			  set_value(k, h[k].to_s)
			end
		  rescue
		    mes = "Cannot evaluate #{value(k)}: " + $!.to_s
			Dialog.run("Value Error") {
			  layout(1, item(:text, :title=>mes))
			  set_attr(1, :hidden=>true)
			}
			return nil
		  end
		}
		ax = Vector3D[h["a0"], h["a1"], h["a2"]]
		bx = Vector3D[h["b0"], h["b1"], h["b2"]]
		cx = Vector3D[h["c0"], h["c1"], h["c2"]]
		ox = Vector3D[h["o0"], h["o1"], h["o2"]]
		@mol.set_box(ax, bx, cx, ox)
		return @mol
	  end
#	  def action(item1)
#	    if item1[:index] == 0
#		  if !set_box_value(item1)
#		    return  #  Cannot set box: dialog is not dismissed
#		  end
#		end
#		super
#	  end
	  box = @mol.box
	  layout(4,
	    item(:text, :title=>"Unit cell:"),
		-1, -1, -1,
	    item(:text, :title=>"origin"),
		item(:textfield, :width=>140, :tag=>"o0", :value=>(box ? box[3].x.to_s : "")),
		item(:textfield, :width=>140, :tag=>"o1", :value=>(box ? box[3].y.to_s : "")),
		item(:textfield, :width=>140, :tag=>"o2", :value=>(box ? box[3].z.to_s : "")),
	    item(:text, :title=>"a-axis"),
		item(:textfield, :width=>140, :tag=>"a0", :value=>(box ? box[0].x.to_s : "")),
		item(:textfield, :width=>140, :tag=>"a1", :value=>(box ? box[0].y.to_s : "")),
		item(:textfield, :width=>140, :tag=>"a2", :value=>(box ? box[0].z.to_s : "")),
	    item(:text, :title=>"b-axis"),
		item(:textfield, :width=>140, :tag=>"b0", :value=>(box ? box[1].x.to_s : "")),
		item(:textfield, :width=>140, :tag=>"b1", :value=>(box ? box[1].y.to_s : "")),
		item(:textfield, :width=>140, :tag=>"b2", :value=>(box ? box[1].z.to_s : "")),
	    item(:text, :title=>"c-axis"),
		item(:textfield, :width=>140, :tag=>"c0", :value=>(box ? box[2].x.to_s : "")),
		item(:textfield, :width=>140, :tag=>"c1", :value=>(box ? box[2].y.to_s : "")),
		item(:textfield, :width=>140, :tag=>"c2", :value=>(box ? box[2].z.to_s : "")),
		item(:button, :title=>"Set", :action=>:set_box_value),
		item(:text, :title=>"(Ruby expressions are allowed as the values)"),
		-1, -1)
	  set_attr(0, :action=>proc { |it| set_box_value(it) && end_modal(it) })
	}
  end

  def cmd_show_periodic_image
    mol = self
    hash = Dialog.run("Show Periodic Image") {
	  @mol = mol
	  def set_periodic_image(it)
	    a = []
	    ["amin", "amax", "bmin", "bmax", "cmin", "cmax"].each_with_index { |k, i|
		  s = value(k)
		  if s == nil || s == ""
		    a[i] = 0
		  else
		    a[i] = Integer(s)
		  end
		}
	    @mol.show_periodic_image(a)
	  end
	  pimage = @mol.show_periodic_image
	  layout(4,
	    item(:text, :title=>"Show Periodic Image:"),
		-1, -1, -1,
		item(:text, :title=>"a-axis"),
		item(:textfield, :width=>80, :tag=>"amin", :value=>pimage[0].to_s),
		item(:text, :title=>"to"),
		item(:textfield, :width=>80, :tag=>"amax", :value=>pimage[1].to_s),
		item(:text, :title=>"b-axis"),
		item(:textfield, :width=>80, :tag=>"bmin", :value=>pimage[2].to_s),
		item(:text, :title=>"to"),
		item(:textfield, :width=>80, :tag=>"bmax", :value=>pimage[3].to_s),
		item(:text, :title=>"c-axis"),
		item(:textfield, :width=>80, :tag=>"cmin", :value=>pimage[4].to_s),
		item(:text, :title=>"to"),
		item(:textfield, :width=>80, :tag=>"cmax", :value=>pimage[5].to_s),
		item(:button, :title=>"Set", :action=>:set_periodic_image))
	  set_attr(0, :action=>proc { |it| set_periodic_image(it); end_modal(it) } )
	}
  end

  def cmd_pressure_control
  
    if box == nil
	  Dialog.run("Pressure Control Error") {
		layout(1,
		  item(:text, :title=>"Unit cell is not defined. Please open 'Unit Cell...' \nfrom the MM/MD menu and set up the unit cell."))
		set_attr(1, :hidden=>true)  #  Hide the cancel button
	  }
	  return nil
	end
	
  	#  Create arena
	arena = self.prepare_arena
    if !arena
	  return nil
	end
    mol = self
	
	#  Create pressure arena
	if arena.pressure_freq == nil
	  arena.pressure_freq = 20   #  This assignment will automatically create pressure arena
	end
	
	hash = Dialog.run("Pressure Control") {
	  @mol = mol
	  layout(2,
	    item(:text, :title=>"Frequency"),
		item(:textfield, :width=>100, :tag=>"freq", :value=>arena.pressure_freq.to_s),
		item(:text, :title=>"Coupling"),
		item(:textfield, :width=>100, :tag=>"coupling", :value=>arena.pressure_coupling.to_s),
	    item(:text, :title=>"Pressure a"),
		item(:textfield, :width=>100, :tag=>"pa", :value=>arena.pressure[0].to_s),
	    item(:text, :title=>"Pressure b"),
		item(:textfield, :width=>100, :tag=>"pb", :value=>arena.pressure[1].to_s),
	    item(:text, :title=>"Pressure c"),
		item(:textfield, :width=>100, :tag=>"pc", :value=>arena.pressure[2].to_s),
	    item(:text, :title=>"Cell flexibility a"),
		item(:textfield, :width=>100, :tag=>"fa", :value=>arena.pressure_cell_flexibility[0].to_s),
	    item(:text, :title=>"Cell flexibility b"),
		item(:textfield, :width=>100, :tag=>"fb", :value=>arena.pressure_cell_flexibility[1].to_s),
	    item(:text, :title=>"Cell flexibility c"),
		item(:textfield, :width=>100, :tag=>"fc", :value=>arena.pressure_cell_flexibility[2].to_s),
	    item(:text, :title=>"Cell flexibility bc"),
		item(:textfield, :width=>100, :tag=>"fbc", :value=>arena.pressure_cell_flexibility[3].to_s),
	    item(:text, :title=>"Cell flexibility ca"),
		item(:textfield, :width=>100, :tag=>"fca", :value=>arena.pressure_cell_flexibility[4].to_s),
	    item(:text, :title=>"Cell flexibility ab"),
		item(:textfield, :width=>100, :tag=>"fab", :value=>arena.pressure_cell_flexibility[5].to_s),
	    item(:text, :title=>"Cell flexibility origin"),
		item(:textfield, :width=>100, :tag=>"forig", :value=>arena.pressure_cell_flexibility[6].to_s),
	    item(:text, :title=>"Cell flexibility orientation"),
		item(:textfield, :width=>100, :tag=>"forient", :value=>arena.pressure_cell_flexibility[7].to_s),
	    item(:text, :title=>"Fluctuate cell origin"),
		item(:textfield, :width=>100, :tag=>"orig", :value=>arena.pressure_fluctuate_cell_origin.to_s),
	    item(:text, :title=>"Fluctuate cell orientation"),
		item(:textfield, :width=>100, :tag=>"orient", :value=>arena.pressure_fluctuate_cell_orientation.to_s))
    }
	if hash[:status] == 0
	  arena.pressure_freq = hash["freq"]
	  arena.pressure_coupling = hash["coupling"]
	  arena.pressure = [hash["pa"], hash["pb"], hash["pc"]]
	  arena.pressure_cell_flexibility = [hash["fa"], hash["fb"], hash["fc"], hash["fbc"], hash["fca"], hash["fab"], hash["forig"], hash["forient"]]
	  arena.pressure_fluctuate_cell_origin = hash["orig"]
	  arena.pressure_fluctuate_cell_orientation = hash["orient"]
	end
  end
  
  def ambertools_dialog(tool)
	ante_dir = get_global_settings("antechamber.ante_dir")
	log_dir = get_global_settings("antechamber.log_dir")
	if !ante_dir
	  ante_dir = ($platform == "mac" ? "/Applications/amber10" : ($platform == "win" ? "c:/opt/amber10" : ""))
	end
	if !log_dir
	  log_dir = document_home + "/amber10"
	end
	if $platform == "win"
	  suffix = ".exe"
	else
	  suffix = ""
	end
	log_level = (get_global_settings("antechamber.log_level") || "none")
    hash = Dialog.run("Run " + tool.capitalize) {
	  @toolname = tool + suffix
      def valid_antechamber_dir(s)
	    FileTest.exist?(s + "/" + @toolname)
  	  end
	  layout(2,
		item(:checkbox, :title=>"Optimize structure and calculate charges (may be slow)", :tag=>"calc_charge",
		  :value=>(get_global_settings("antechamber.calc_charge") || 0),
		  :action=>proc { |it| set_attr("nc", :enabled=>(it[:value] != 0)) } ),
		-1,
		item(:text, :title=>"      Net Molecular Charge:"),
		item(:textfield, :width=>"80", :tag=>"nc", :value=>(get_global_settings("antechamber.nc") || "0")),
		item(:checkbox, :title=>"Use the residue information", :tag=>"use_residue",
		  :value=>(get_global_settings("antechamber.use_residue") || 0),
		  :action=>proc { |it| valid_antechamber_dir(it[:value]) && set_attr(0, :enabled=>true) } ),
		-1,
		item(:line),
		-1,
		item(:text, :title=>"Log directory:"),
		[ item(:button, :title=>"Choose...",
			:action=>proc { |it|
			  dir = Dialog.open_panel(nil, nil, nil, true)
			  if dir
				set_value("log_dir", dir)
			  end
			}
		  ), {:align=>:right} ],
		item(:textfield, :width=>360, :height=>40, :tag=>"log_dir", :value=>log_dir),
		-1,
		item(:text, :title=>"Log handling"),
		-1,
		layout(1,
		  item(:radio, :title=>"Do not keep logs", :tag=>"log_none", :value=>(log_level == "none" ? 1 : 0)),
		  item(:radio, :title=>"Keep only when error occurred", :tag=>"log_error_only", :value=>(log_level == "error_only" ? 1 : 0)),
		  layout(3,
		    item(:radio, :title=>"Keep latest ", :tag=>"log_keep_latest", :value=>(log_level == "latest" ? 1 : 0)),
		    item(:textfield, :width=>80, :tag=>"log_keep_number", :value=>(get_global_settings("antechamber.log_keep_number") || "10")),
		    item(:text, :title=>"logs"),
		    {:margin=>0, :padding=>0}),
		  item(:radio, :title=>"Keep All", :tag=>"log_all", :value=>(log_level == "all" ? 1 : 0)),
		  {:margin=>0, :padding=>0}
		),
		-1
	  )
	  set_attr("nc", :enabled=>(attr("calc_charge", :value) == "1"))
    }
	hash.each_pair { |key, value|
	  next if key == :status
	  v = [["log_none", "none"], ["log_error_only", "error_only"], ["log_keep_latest", "latest"], ["log_all", "all"]].assoc(key)
      if v 
	    next if value != 1
	    tag = "log_level"
	    value = v[1]
      end
      set_global_settings("antechamber.#{tag}", value)
    }

	#  The hash values are set in the action() method
#	print "retval = #{hash ? 1 : 0}\n"
	return (hash[:status] == 0 ? 1 : 0)
  end
  
  def create_ante_log_dir(name, key)
    log_dir = get_global_settings("antechamber.log_dir")
	msg = ""
	begin
	  if !FileTest.directory?(log_dir)
	    Dir.mkdir(log_dir) rescue ((msg = "Cannot create directory #{log_dir}") && raise)
	  end
	  n = 1
	  while FileTest.exist?(dname = log_dir + "/#{name}_#{key}.#{n}")
	    n += 1
	  end
	  Dir.mkdir(dname) rescue ((msg = "Cannot create directory #{dname}") && raise)
	  return dname
	rescue
	  error_message_box(msg + ": " + $!.to_s)
	  return
	end
  end
  
  def clean_ante_log_dir(nkeep)
    def rm_recursive(f)
	  if FileTest.directory?(f)
	    Dir.entries(f).each { |file|
		  next if file == "." || file == ".."
		  rm_recursive(f + "/" + file)
		}
		Dir.rmdir(f)
	  else
	    File.delete(f)
	  end
	end
    log_dir = get_global_settings("antechamber.log_dir")
	cwd = Dir.pwd
	count = 0
	begin
	  Dir.chdir(log_dir)
	  #  Get subdirectories and sort by last modified time
	  dirs = Dir.entries(log_dir).reject! { |x|
	    !FileTest.directory?(x) || x == "." || x == ".."
	  }.sort_by { |x| 
	    File.mtime(x)
	  }
	  dirs[0..-(nkeep + 1)].each { |d|
	    rm_recursive(d)
		count += 1
	  }
	rescue
	  error_message_box $!.to_s
	  Dir.chdir(cwd)
	  return
	end
	Dir.chdir(cwd)
	return count
  end
  
  def cmd_antechamber
    return ambertools_dialog("antechamber")
  end
  
  def import_ac(acfile)
    open(acfile, "r") { |fp|
	  while (s = fp.gets)
	    next if s !~ /^ATOM/
		s.chomp!
		idx = Integer(s[4..11]) - 1
		charge = Float(s[54..63])
		type = s[72..-1]
		type.gsub!(/ /, "")
		ap = atoms[idx]
		ap.charge = charge
		ap.atom_type = type
	  end
	}
  end
  
  def import_sqmout(file)
    open(file, "r") { |fp|
	  while (s = fp.gets)
	    next if s !~ /Final Structure/
		s = fp.gets
		s = fp.gets
		s = fp.gets
		idx = 0
		while (s = fp.gets)
		  break if s !~ /QMMM/
		  a = s.split
		  r = Vector3D[Float(a[4]), Float(a[5]), Float(a[6])]
		  atoms[idx].r = r
		  idx += 1
		end
		break
	  end
	}
  end
  
  def import_frcmod(file)
    self.md_arena.prepare(true)  #  Clean up existing parameters
	par = self.parameter
    open(file, "r") { |fp|
	  wtable = Hash.new
	  state = 0
	  while (s = fp.gets)
	    s.chomp!
		case s
		when /^MASS/
		  state = 1
		  next
		when /^BOND/
		  state = 2
		  next
		when /^ANGLE/
		  state = 3
		  next
		when /^DIHE/
		  state = 4
		  next
		when /^IMPR/
		  state = 5
		  next
		when /^NONB/
		  state = 6
		  next
		when ""
		  state = 0
		  next
		else
		  case state
		  when 1
			name, weight = s.split
			wtable[name] = Float(weight)
		  when 2
		    types, k, r0, com = s.split(nil, 4)
		    pp = par.bonds.lookup(types, :local, :missing) || par.bonds.insert
			pp.atom_types = types
			pp.k = k
			pp.r0 = r0
			pp.comment = com
		  when 3
		    types, k, a0, com = s.split(nil, 4)
			pp = par.angles.lookup(types, :local, :missing) || par.angles.insert
			pp.atom_types = types
			pp.k = k
			pp.a0 = a0
			pp.comment = com
		  when 4
		    types, n, k, phi0, period, com = s.split(nil, 6)
			pp = par.dihedrals.lookup(types, :local, :missing) || par.dihedrals.insert
			pp.atom_types = types
			pp.mult = 1
			pp.k = k
			pp.phi0 = phi0
			pp.period = Float(period).round
			pp.comment = com
		  when 5
		    types, k, phi0, period, com = s.split(nil, 5)
			pp = par.impropers.lookup(types, :local, :missing) || par.impropers.insert
			pp.atom_types = types
			pp.mult = 1
			pp.k = k
			pp.phi0 = phi0
			pp.period = Float(period).round
			pp.comment = com
		  when 6
		    name, r_eq, eps, com = s.split(nil, 4)
			pp = par.vdws.lookup(name, :local, :missing) || par.vdws.insert
			pp.atom_type = name
			pp.r_eq = r_eq
			pp.r_eq14 = r_eq
			pp.eps = eps
			pp.eps14 = eps
			if wtable[name]
			  pp.weight = wtable[name]
			end
			pp.comment = com
		  end
		end
	  end
    }
  end
  
  def count_elements
    elements = []
	each_atom { |ap|
	  if (p = elements.assoc(ap.atomic_number))
	    p[1] += 1
	  else
	    elements.push [ap.atomic_number, 1]
	  end
	}
	elements.sort_by { |p| p[0] }
  end
  
  def cmd_run_resp
    if natoms == 0
	  error_message_box "Molecule is empty"
	  return
	elsif nelpots == 0
	  error_message_box "No ESP information is loaded"
	  return
	end
	return unless ambertools_dialog("resp")
	nc = get_global_settings("antechamber.nc")
	ante_dir = get_global_settings("antechamber.ante_dir")

	#  Create the temporary directory
	dname = create_ante_log_dir((self.path ? File.basename(self.path, ".*") : self.name), "rs")
	return unless dname
	cwd = Dir.pwd
	Dir.chdir(dname)

	if true
	  eq = search_equivalent_atoms
	  #  search for methyl (and methylene??) hydrogens
	  methyl = []
	  each_atom { |ap|
	    if ap.element == "C" && ap.connects.length == 4
		  hs = ap.connects.select { |n| atoms[n].element == "H" }
		  if hs.length == 3
		    mc = hs.min
		    hs.each { |n| methyl[n] = mc }  #  methyl hydrogens
			methyl[ap.index] = -1           #  methyl carbon
		  end
		end
	  }
	  for i in [1,2]
	    open("resp.input#{i}", "w") { |fp|
		  fp.print " #{self.name} \n"
		  fp.print " &cntrl\n"
		  fp.print " ioutput=1, IQOPT=#{i}, nmol=1, ihfree=1, irstrnt=1, qwt=#{i*0.0005},\n"
		  fp.print " &end\n"
		  fp.print "  1.0\n"
		  fp.print " #{self.name} \n"
		  fp.printf "%4d%5d\n", nc, self.natoms
		  each_atom { |ap|
		    idx = ap.index
			if i == 1
			  n = (methyl[idx] ? 0 : eq[idx] + 1)
			  n = 0 if n == idx + 1
			else
			  n = (methyl[idx] ? methyl[idx] + 1 : -1)
			  n = 0 if n == idx + 1
			end
			fp.printf "%4d%5d\n", ap.atomic_number, n
		  }
		  fp.print "\n\n\n\n\n\n"
		}
	  end
	else
      #  Export as an antechamber format
  	  c = count_elements
	  formula = c.map { |p| Parameter.builtin.elements[p[0]].name + p[1].to_s }.join(" ")
	  open("respgen_in.ac", "w") { |fp|
	    fp.printf("CHARGE %9.2f ( %d )\n", nc, nc)
	    fp.printf("Formula: #{formula}\n")
	    each_atom { |ap|
	      fp.printf("ATOM %6d  %-4s%-3s%6d%12.3f%8.3f%8.3f%10.6f%8s\n", ap.index + 1, ap.name, ap.res_name, ap.res_seq, ap.r.x, ap.r.y, ap.r.z, ap.charge, ap.atom_type)
	    }
	    bonds.each_with_index { |b, i|
	      fp.printf("BOND%5d%5d%5d%5d  %5s%5s\n", i + 1, b[0] + 1, b[1] + 1, 0, atoms[b[0]].name, atoms[b[1]].name)
	    }
	  }
	
	  #  Create resp input by respgen
	  if !system("#{ante_dir}/respgen -i respgen_in.ac -o resp.input1 -f resp1") \
	  || !system("#{ante_dir}/respgen -i respgen_in.ac -o resp.input2 -f resp2")
	    error_message_box("Cannot run respgen.")
	    Dir.chdir(cwd)
	    return
	  end
	end
	
	#  Create ESP file
	a2b = 1.8897259885   #  angstrom to bohr
	open("resp.esp", "w") { |fp|
	  fp.printf("%5d%5d%5d\n", natoms, nelpots, nc)
	  each_atom { |ap|
	    fp.printf("                %16.7E%16.7E%16.7E\n", ap.r.x * a2b, ap.r.y * a2b, ap.r.z * a2b)
	  }
	  nelpots.times { |i|
	    pos, esp = elpot(i)
	    fp.printf("%16.7E%16.7E%16.7E%16.7E\n", esp, pos.x, pos.y, pos.z)
	  }
	  fp.print "\n\n"
	}

    #  Run resp
	if !system("#{ante_dir}/resp -O -i resp.input1 -o resp.output1 -e resp.esp -t qout_stage1") \
	|| !system("#{ante_dir}/resp -O -i resp.input2 -o resp.output2 -e resp.esp -q qout_stage1 -t qout_stage2")
	  error_message_box("Cannot run resp.")
	  Dir.chdir(cwd)
	  return 
	end
	
	#  Import resp output
	open("punch", "r") { |fp|
	  while (s = fp.gets)
	    next unless s =~ /Point charges/ && s =~ /after optimization/
		s = fp.gets
		i = 0
		while (s = fp.gets)
		  ary = s.split
		  break if ary.count == 0
		  if Integer(ary[1]) != atoms[i].atomic_number
		    error_message_box(sprintf("The atom %d has inconsistent atomic number (%d in mol, %d in resp output)", i, atoms[i].atomic_number, ary[1]))
			Dir.chdir(cwd)
			return
		  end
		  atoms[i].charge = Float(ary[3])
		  i += 1
		end
		break
	  end
    }
	Dir.chdir(cwd)
	return true
  end
  
  def cmd_gamess_resp
	if natoms == 0
	  error_message_box "Molecule is empty"
	  return
	end
	mol = self;
	Dialog.run("#{name}:GAMESS/RESP", nil, nil) {
	  layout(1,
	    item(:text, :title=>"Step 1:\nCreate GAMESS input for ESP calculation"),
		[item(:button, :title=>"Create GAMESS Input...",
		  :action=>proc { |it|
		    esp_save = get_global_settings("gamess.esp")
			set_global_settings("gamess.esp", 1)
			mol.cmd_create_gamess_input
			set_global_settings("gamess.esp", esp_save)
		  }),
		  {:align=>:right}],
		item(:line),
		item(:text, :title=>"Step 2:\nImport GAMESS .dat file"),
		[item(:button, :title=>"Import GAMESS dat...",
		  :action=>proc { |it|
		    fname = Dialog.open_panel("Select GAMESS .dat file", nil, "*.dat")
			if fname
			  errmsg = nil
			  begin
			    mol.loaddat(fname)
				errmsg = "Cannot find ESP results in the dat file." if mol.nelpots == 0
			  rescue
			    errmsg = "Error reading the dat file."
			  end
			  if errmsg
			    message_box(errmsg + "\nYou may want to try another .dat file.", "GAMESS Import Error", :ok, :warning)
			  else
			    set_attr("resp", :enabled=>true)
			  end
			end
		  }),
		  {:align=>:right}],
		item(:line),
		item(:text, :title=>"Step 3:\nRun RESP for charge fitting"),
		[item(:button, :title=>"Run RESP...", :tag=>"resp",
		  :action=>proc { |it|
		    if mol.cmd_run_resp
			  if get_global_settings("antechamber.log_level") == "latest"
			    mol.clean_ante_log_dir(Integer(get_global_settings("antechamber.log_keep_number")))
			  end
			  end_modal(0)
			end
		  }),
		  {:align=>:right}],
		item(:line),
		[item(:button, :title=>"Close", :action=>proc { |it| end_modal(1) }),
		  {:align=>:center}]
	  )
	  if mol.nelpots == 0
	    set_attr("resp", :enabled=>false)
	  end
	}
  end
  
  def cmd_edit_local_parameter_in_mainview(ptype, names, types, value, params)
    #  The parameters are given as space separated strings
	#  e.g. "1.862 200.000"
	case ptype
	when "bond"
	  k = ["r0", "k"]
	  pen = self.parameter.bonds
	when "angle"
	  k = ["a0", "k"]
	  pen = self.parameter.angles
	when "dihedral"
	  k = ["k", "period", "phi0"]
	  pen = self.parameter.dihedrals
	when "improper"
	  k = ["k", "period", "phi0"]
	  pen = self.parameter.impropers
	else
	  return
	end
	p = params.split
	hash = Dialog.run("Edit local parameter") {
      layout(1,
	    item(:text, :title=>"Edit #{ptype} parameter for #{names} (#{types})"),
		item(:text, :title=>"(Current value = #{value})"),
	    layout(4,
	      [item(:text, :title=>"types"), {:align=>:center}],
	      [item(:text, :title=>k[0]), {:align=>:center}],
	      [item(:text, :title=>k[1]), {:align=>:center}],
	      (k[2] ? [item(:text, :title=>k[2]), {:align=>:center}] : -1),
		
		  item(:textfield, :width=>100, :value=>types, :tag=>"types"),
		  item(:textfield, :width=>100, :value=>p[0], :tag=>k[0]),
		  item(:textfield, :width=>100, :value=>p[1], :tag=>k[1]),
		  (k[2] ? item(:textfield, :width=>100, :value=>p[2], :tag=>k[2]) : -1)
		)
	  )
	}
	if hash[:status] == 0
	  pref = pen.lookup(hash["types"], :create, :local, :missing, :nobasetype, :nowildcard)
	  raise "Cannot create new #{ptype} parameter" if !pref
	  k.each { |key| pref.set_attr(key, hash[key]) }
	  self.md_arena.prepare(true)  #  Check parameter only
	  return 1
	else
	  return 0
	end
  end
  
end
