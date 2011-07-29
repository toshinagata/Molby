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
    hash = Dialog.run("Molecular Dynamics Advanced Settings", "Close", nil) {
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
#	if hash[:status] == 0
	  hash.keys.each { |k|
	    next if k == :status || read_only.include?(k)
	    arena[k] = hash[k]
	  }
	  return hash
#	  arena.prepare
#	  return hash
#	else
#	  return nil
#	end
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
	    arena.log_file = self.name.gsub(/\.\w+$/, ".log")
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
#	[:log_file, :coord_file, :vel_file, :force_file, :debug_file].each { |k|
#	  s = arena[k]
#	  files_save[k] = s
#	  if s != nil
#	    arena[k] = File.basename(s)
#	  end
#	}
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
	if !minimize
	  arena.temperature = Float(hash["temperature"])
	  arena.timestep = Float(hash["timestep"])
	end
	arena.coord_output_freq = Integer(hash["steps_per_frame"])
	arena.energy_output_freq = arena.coord_output_freq
	arena.log_file = hash["log_file"]
	if hash[:status] == 0
	  return Integer(hash["number_of_frames"])
	else
	  return -1
	end
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
		@mol.show_periodic_image = (attr("show_flag", :value) == 1 ? true : false)
	  end
	  pimage = @mol.show_periodic_image
	  flag = @mol.show_periodic_image?
	  layout(4,
	    item(:checkbox, :title=>"Show Periodic Image", :tag=>"show_flag", :value=>(flag ? 1 : 0),
		  :action=>proc { |it| @mol.show_periodic_image = (it[:value] == 1 ? true : false) } ),
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
	log_dir = get_global_settings("antechamber.log_dir")
	if !log_dir
	  log_dir = document_home + "/antechamber"
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
	  set_attr("nc", :enabled=>(attr("calc_charge", :value) != 0))
	  set_attr("calc_charge", :enabled=>(tool != "resp"))
    }
	hash.each_pair { |key, value|
	  next if key == :status
	  v = [["log_none", "none"], ["log_error_only", "error_only"], ["log_keep_latest", "latest"], ["log_all", "all"]].assoc(key)
      if v 
	    next if value != 1
	    key = "log_level"
	    value = v[1]
      end
      set_global_settings("antechamber.#{key}", value)
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
		when /^ANGL/
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
		    types, k, r0, com = s.split(" ", 4)
		    pp = par.bonds.lookup(types, :local, :missing) || par.bonds.insert
			pp.atom_types = types
			pp.k = k
			pp.r0 = r0
			pp.comment = com
		  when 3
		    types, k, a0, com = s.split(" ", 4)
			pp = par.angles.lookup(types, :local, :missing) || par.angles.insert
			pp.atom_types = types
			pp.k = k
			pp.a0 = a0
			pp.comment = com
		  when 4
		    types, n, k, phi0, period, com = s.split(" ", 6)
			pp = par.dihedrals.lookup(types, :local, :missing) || par.dihedrals.insert
			pp.atom_types = types
			pp.mult = 1
			pp.k = k
			pp.phi0 = phi0
			pp.period = Float(period).round
			pp.comment = com
		  when 5
		    types, k, phi0, period, com = s.split(" ", 5)
			pp = par.impropers.lookup(types, :local, :missing) || par.impropers.insert
			pp.atom_types = types
			pp.mult = 1
			pp.k = k
			pp.phi0 = phi0
			pp.period = Float(period).round
			pp.comment = com
		  when 6
		    name, r_eq, eps, com = s.split(" ", 4)
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

  def cmd_import_frcmod
    file = Dialog.open_panel("Select AMBER frcmod file", nil, "AMBER frcmod file (*.frcmod)|*.frcmod|All Files (*.*)|*.*")
	import_frcmod(file) if file
  end
  
  def Molecule.import_amberlib(file)
    fp = File.open(file, "r")
	raise MolbyError, "Cannot open file #{file}" if fp == nil
	mols = Hash.new
	while line = fp.gets
	  if line =~ /^!!index/
	    while line = fp.gets
		  break if line =~ /^!/
		  units = line.scan(/"([^\"]*)"/).flatten
		end
	  elsif line =~ /^!entry\.(\w+)\.unit\.(\w+)\s/
	    unit = $1
		cmd = $2
		mol = (mols[unit] ||= Molecule.new)
		if cmd == "atoms"
		  while line = fp.gets
		    break if line =~ /^!/
		    line.chomp!
			en = line.split
			en[0].delete!("\"")
			en[1].delete!("\"")
			ap = mol.add_atom(en[0], en[1])
			elem = Integer(en[6])
		    if elem > 0
			  ap.atomic_number = elem
			else
			  ap.element = en[1][0, 1].upcase
			end
			ap.charge = Float(en[7])
		  end
		elsif cmd == "boundbox"
		  values = []
		  5.times { values.push(Float(fp.gets)) }
		  mol.set_cell([values[2], values[3], values[4], 90, 90, 90])
		  line = fp.gets
		elsif cmd == "connectivity"
		  while line = fp.gets
		    break if line =~ /^!/
			line.chomp!
			con = line.split
			begin
			mol.create_bond(Integer(con[0]) - 1, Integer(con[1]) - 1)
			rescue
			puts "#{$!}: con[0] = #{con[0]}, con[1] = #{con[1]}"
			p "natoms = #{mol.natoms}"
			raise
			end
		  end
		elsif cmd == "positions"
		  index = 0
		  while line = fp.gets
		    break if line =~ /^!/
			line.chomp!
			pos = line.split
			mol.atoms[index].r = [Float(pos[0]), Float(pos[1]), Float(pos[2])]
			index += 1
		  end
		elsif cmd == "residues"
		  resseq = 0
		  name = nil
		  n1 = 0
		  while line = fp.gets
		    if line =~ /^!/
			  values[0] = nil
			  values[3] = mol.natoms + 1
			else
			  values = line.split
			  values[0].delete!("\"")
			  values[3] = Integer(values[3])
			end
			if name
			  mol.assign_residue(IntGroup[n1..values[3] - 2], sprintf("%s.%d", name, resseq))
			end
			name = values[0]
			n1 = values[3] - 1
			resseq += 1
			break if name == nil
		  end
		else
		  while line = fp.gets
			break if line =~ /^!/
		  end
		end
	  end
	  redo if line
	end
	fp.close
	msg = "Imported unit(s):"
	mols.each { |key, val|
	  mol = Molecule.open
	  mol.add(val)
	  if val.cell != nil
		mol.set_cell(val.cell)
	  end
	  if mol.natoms > 1000
	    mol.line_mode(true)
	  end
	  mol.resize_to_fit
	  msg += "\n  #{key} (#{mol.name})"
	  if mol.md_arena.prepare(true) == nil
	    h = Dialog.run("Missing parameters") {
		  layout(1, item(:text, :title=>"Some parameters are missing for unit \"#{key}\". Do you want to import an frcmod file?"))
		}
		if (h[:status] == 0)
		  file = Dialog.open_panel("Select frcmod file", nil, "AMBER frcmod (*.frcmod)|*.frcmod|All Files (*.*)|*.*")
		  if file
		    mol.instance_eval { import_frcmod(file) }
		  end
		end
	  end
	}
	message_box(msg, "AMBER Lib Import Complete", :ok)
	puts msg
  end
  
  def Molecule.cmd_import_amberlib
    file = Dialog.open_panel("Select AMBER lib file", nil, "AMBER lib file (*.lib)|*.lib|All Files (*.*)|*.*")
	import_amberlib(file) if file
  end
  
  def export_ac(acfile)
    open(acfile, "w") { |fp|
	  charge = 0.0
	  each_atom { |ap|
		charge += ap.charge
	  }
	  fp.printf "CHARGE    %6.3f\n", charge
	  form = []
	  each_atom { |ap|
		form[ap.atomic_number] = (form[ap.atomic_number] || 0) + 1
	  }
	  fp.print "Formula: "
	  form.each_with_index { |n, i|
		if n != nil
		  fp.print Parameter.builtin.elements[i].name + n.to_s
		end
	  }
	  fp.print "\n"
	  each_atom { |ap|
		fp.printf "ATOM %6d  %-4s%3.3s     1    %8.3f%8.3f%8.3f%10.5f      %-4s    %s\n", ap.index + 1, ap.name, ap.res_name, ap.x, ap.y, ap.z, ap.charge, ap.atom_type, ap.name
	  }
	  bonds.each_with_index { |bd, i|
		n1, n2 = bd
		fp.printf "BOND%5d%5d%5d%5d   %4s%4s\n", i + 1, n1 + 1, n2 + 1, 0, atoms[n1].name, atoms[n2].name
	  }
	}
  end
  
  def export_frcmod(frcfile)
    p = self.parameter
	open(frcfile, "w") { |fp|
	  fp.print "remark goes here\n"
	  fp.print "MASS\n"
	  p.vdws.each { |np|
		if np.source != "gaff" && np.source != "parm99"
		  fp.printf "%s   %10.5f    %s\n", np.atom_type, np.weight, np.comment
		end
	  }
	  fp.print "\n"
	  fp.print "BOND\n"
	  p.bonds.each { |bp|
		if bp.source != "gaff" && bp.source != "parm99"
		  types = sprintf("%-2s-%-2s", bp.atom_types[0], bp.atom_types[1])
		  fp.printf "%s  %6.3f     %7.3f   %s\n", types, bp.k, bp.r0, bp.comment
		end
	  }
	  fp.print "\n"
	  fp.print "ANGLE\n"
	  p.angles.each { |ap|
		if ap.source != "gaff" && ap.source != "parm99"
		  types = sprintf("%-2s-%-2s-%-2s", ap.atom_types[0], ap.atom_types[1], ap.atom_types[2])
		  fp.printf "%s  %7.3f     %7.3f   %s\n", types, ap.k, ap.a0, ap.comment
		end
	  }
	  fp.print "\n"
	  fp.print "DIHE\n"
	  p.dihedrals.each { |dp|
		if dp.source != "gaff" && dp.source != "parm99"
		  types = sprintf("%-2s-%-2s-%-2s-%-2s", dp.atom_types[0], dp.atom_types[1], dp.atom_types[2], dp.atom_types[3])
		  fp.printf "%s   %d   %6.3f      %7.3f          %6.3f   %s\n", types, dp.mult, dp.k, dp.phi0, dp.period, dp.comment
		end
	  }
	  fp.print "\n"
	  fp.print "IMPROPER\n"
	  p.impropers.each { |dp|
		if dp.source != "gaff" && dp.source != "parm99"
		  types = sprintf("%-2s-%-2s-%-2s-%-2s", dp.atom_types[0], dp.atom_types[1], dp.atom_types[2], dp.atom_types[3])
		  fp.printf "%s      %6.2f      %7.3f          %6.2f   %s\n", types, dp.k, dp.phi0, dp.period, dp.comment
		end
	  }
	  fp.print "\n"
	  fp.print "NONBON\n"
	  p.vdws.each { |np|
		if np.source != "gaff" && np.source != "parm99"
		  types = sprintf("%2s", np.atom_types[0])
		  fp.printf "%s   %10.5f    %10.5f    %s\n", np.atom_type, np.r_eq, np.eps, np.comment
		end
	  }
	  fp.print "\n"
	}
  end
  
  # Create sander input files from current molecule  
  def export_prmtop
  
    def error_dialog(msg)
      Dialog.run("AMBER Export Error") {
        layout(1, item(:text, :title=>msg))
        set_attr(1, :hidden=>true)
      }
    end
    
    def format_print(fp, num, fmt, ary)
      fmts = "%#{fmt}%s"
      ary.each_with_index { |x, i|
        fp.printf(fmts, x, (i % num == num - 1 ? "\n" : ""))
      }
      fp.print "\n" if ary.count % num != 0
    end
  
    par = self.parameter
    
    #  Create residue number table
    restable = []
    resno = -1
    self.each_atom { |ap|
      if ap.res_seq > resno
        restable.push(ap.index)
        resno = ap.res_seq
      elsif ap.res_seq < resno
        error_dialog("The residue number is not in ascending order. Please use 'sort by residue' command before creating AMBER input.")
        return nil
      end
    }
    
    #  Check if periodic box is used
    periodic = (self.box && self.box[4].all? { |n| n != 0 })
    if periodic
      fragments = []
      last_solute = nil
      first_solv_mol = nil
      self.each_fragment { |gr|
        if gr.range_at(1) != nil
          msg = gr.to_s.sub(/IntGroup/, "")
          error_dialog("The atoms #{msg} are one fragment, but the atom indices are not consecutive.")
          return nil
        end
        fragments.push(gr.length)
        if last_solute == nil && self.atoms[gr[0]].seg_name == "SOLV"
          #  Calculate the residue number of the last solute atom (gr[0] - 1)
          n = gr[0] - 1
          if n < 0
            error_dialog("All atoms are labeled as the solvent??")
            return nil
          end
          first_solv_mol = fragments.length
          restable.each_with_index { |m, i|
            n -= m
            if n <= 0
              last_solute = i + 1
              break
            end
          }
        end
      }
      if last_solute == nil
        last_solute = restable.length
        first_solv_mol = fragments.length
      end
    end
    
    #  Count bonds and angles, including and not-including hydrogens
    bonds_h = []
    bonds_a = []
    self.bonds.each_with_index { |e, i|
      if e.any? { |n| self.atoms[n].atomic_number == 1 }
        bonds_h.push(i)
      else
        bonds_a.push(i)
      end
    }
    angles_h = []
    angles_a = []
    self.angles.each_with_index { |e, i|
      if e.any? { |n| self.atoms[n].atomic_number == 1 }
        angles_h.push(i)
      else
        angles_a.push(i)
      end
    }
    
    #  Build exclusion table (before processing dihedral angles)
    exnumbers = []
    extable = []
    exhash = {}
    self.each_atom { |ap|
      ex = ap.exclusion
      exx = []
      ex.each_with_index { |x, i|
        x.each { |n|
          if n > ap.index
            exhash[ap.index * 1000000 + n] = i + 2  #  2 for 1-2, 3 for 1-3, 4 for 1-4
            exx.push(n + 1)
          end
        }
      }
      if exx.length > 0
        exx.sort!
      else
        exx.push(0)
      end
      exnumbers.push(exx.length)
      extable.concat(exx)
    }
    
    #  Count dihedrals and impropers, including and not-including hydrogens
    dihedrals_h = []
    dihedrals_a = []
    self.dihedrals.each_with_index { |e, i|
      flag = 0
      k = (e[0] < e[3] ? e[0] * 1000000 + e[3] : e[3] * 1000000 + e[0])
      if exhash[k] == 4
        exhash[k] = -4   #  Handle 1,4-exclusion for only one dihedrals
      elsif exhash[k] == -4
        flag = 1000000   #  Skip calculation of 1,4-interaction to avoid double counting
      elsif exhash[k] && exhash[k] > 0
        flag = 1000000   #  5-member or smaller ring: e[0]...e[3] is 1-2 or 1-3 exclusion pair
      end
      if e.any? { |n| self.atoms[n].atomic_number == 1 }
        dihedrals_h.push(i + flag)
      else
        dihedrals_a.push(i + flag)
      end
    }
    self.impropers.each_with_index { |e, i|
      if e.any? { |n| self.atoms[n].atomic_number == 1 }
        dihedrals_h.push(i + 2000000)
      else
        dihedrals_a.push(i + 2000000)
      end
    }
    
    #  Create vdw parameter table
    nvdw_pars = par.nvdws
    vdw_pair_table = []
    vdw_access_table = []
    (0..nvdw_pars - 1).each { |i|
      p1 = par.vdws[i]
      (0..i).each { |j|
        p2 = par.vdws[j]
        pp = par.lookup("vdw_pair", [p1.atom_type, p2.atom_type])
        if pp
          #  Vdw pair parameter is defined
          eps = pp.eps
          d = pp.r_eq * 2
        else
          eps = Math.sqrt(p1.eps * p2.eps)
          d = p1.r_eq + p2.r_eq
        end
        vdw_access_table[i * nvdw_pars + j] =  vdw_access_table[j * nvdw_pars + i] = vdw_pair_table.length
        vdw_pair_table.push([(d ** 12) * eps, 2 * (d ** 6) * eps])
      }
    }
        
    basename = (self.path ? File.basename(self.path, ".*") : self.name)
    fname = Dialog.save_panel("AMBER prmtop/inpcrd file name", self.dir, basename + ".prmtop", "All files|*.*")
    return nil if !fname
    
    open(fname, "w") { |fp|
      date = Time.now.localtime.strftime("%m/%d/%y  %H:%M:%S")
      
      fp.print "%VERSION  VERSION_STAMP = V0001.000  DATE = #{date}\n"
    
      fp.print "%FLAG TITLE\n%FORMAT(20a4)\n"
      fp.printf "%-80.80s\n", self.name
    
      fp.print "%FLAG POINTERS\n%FORMAT(10I8)\n"
      fp.printf "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n", self.natoms, par.nvdws, bonds_h.length, bonds_a.length, angles_h.length, angles_a.length, dihedrals_h.length, dihedrals_a.length, 0, 0
      fp.printf "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n", extable.length, restable.length, bonds_a.length, angles_a.length, dihedrals_a.length, par.nbonds, par.nangles, par.ndihedrals + par.nimpropers, par.nvdws, 0
      fp.printf "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n", 0, 0, 0, 0, 0, 0, 0, (periodic ? 1 : 0), restable.max, 0
      fp.printf "%8d\n", 0
      
      fp.print "%FLAG ATOM_NAME\n%FORMAT(20a4)\n"
      format_print(fp, 20, "-4.4s", self.atoms.map { |ap| ap.name })
     
      fp.print "%FLAG CHARGE\n%FORMAT(5E16.8)\n"
      format_print(fp, 5, "16.8E", self.atoms.map { |ap| ap.charge * 18.2223 })
      
      fp.print "%FLAG MASS\n%FORMAT(5E16.8)\n"
      format_print(fp, 5, "16.8E", self.atoms.map { |ap| ap.weight })
      
      fp.print "%FLAG ATOM_TYPE_INDEX\n%FORMAT(10I8)\n"
      format_print(fp, 10, "8d", (0...self.natoms).map { |i| self.vdw_par(i).index + 1 })
      
      fp.print "%FLAG NUMBER_EXCLUDED_ATOMS\n%FORMAT(10I8)\n"
      format_print(fp, 10, "8d", exnumbers)
      
      fp.print "%FLAG NONBONDED_PARM_INDEX\n%FORMAT(10I8)\n"
      format_print(fp, 10, "8d", vdw_access_table.map { |n| n + 1 })
      
      fp.print "%FLAG RESIDUE_LABEL\n%FORMAT(20a4)\n"
      format_print(fp, 20, "-4.4s", restable.map { |n| self.atoms[n].res_name })
      
      fp.print "%FLAG RESIDUE_POINTER\n%FORMAT(10I8)\n"
      format_print(fp, 10, "8d", restable.map { |n| n + 1 })
      
      fp.print "%FLAG BOND_FORCE_CONSTANT\n%FORMAT(5E16.8)\n"
      format_print(fp, 5, "16.8E", par.bonds.map { |p| p.k })
      
      fp.print "%FLAG BOND_EQUIL_VALUE\n%FORMAT(5E16.8)\n"
      format_print(fp, 5, "16.8E", par.bonds.map { |p| p.r0 })
      
      fp.print "%FLAG ANGLE_FORCE_CONSTANT\n%FORMAT(5E16.8)\n"
      format_print(fp, 5, "16.8E", par.angles.map { |p| p.k })
      
      fp.print "%FLAG ANGLE_EQUIL_VALUE\n%FORMAT(5E16.8)\n"
      format_print(fp, 5, "16.8E", par.angles.map { |p| p.a0 * Math::PI / 180.0 })
      
      fp.print "%FLAG DIHEDRAL_FORCE_CONSTANT\n%FORMAT(5E16.8)\n"
      format_print(fp, 5, "16.8E", par.dihedrals.map { |p| p.k } + par.impropers.map { |p| p.k })
      
      fp.print "%FLAG DIHEDRAL_PERIODICITY\n%FORMAT(5E16.8)\n"
      format_print(fp, 5, "16.8E", par.dihedrals.map { |p| p.period } + par.impropers.map { |p| p.period })
      
      fp.print "%FLAG DIHEDRAL_PHASE\n%FORMAT(5E16.8)\n"
      format_print(fp, 5, "16.8E", par.dihedrals.map { |p| p.phi0 * Math::PI / 180.0 } + par.impropers.map { |p| p.phi0 * Math::PI / 180.0 })
      
      fp.print "%FLAG SCEE_SCALE_FACTOR\n%FORMAT(5E16.8)\n"
      format_print(fp, 5, "16.8E", par.dihedrals.map { |p| 1.2 } + par.impropers.map { |p| 0.0 })
      
      fp.print "%FLAG SCNB_SCALE_FACTOR\n%FORMAT(5E16.8)\n"
      format_print(fp, 5, "16.8E", par.dihedrals.map { |p| 2.0 } + par.impropers.map { |p| 0.0 })
      
      fp.print "%FLAG SOLTY\n%FORMAT(5E16.8)\n"
      format_print(fp, 5, "16.8E", (0...par.nvdws).map { 0.0 } )
    
      fp.print "%FLAG LENNARD_JONES_ACOEF\n%FORMAT(5E16.8)\n"
      format_print(fp, 5, "16.8E", vdw_pair_table.map { |x| x[0] })
    
      fp.print "%FLAG LENNARD_JONES_BCOEF\n%FORMAT(5E16.8)\n"
      format_print(fp, 5, "16.8E", vdw_pair_table.map { |x| x[1] })
    
      fp.print "%FLAG BONDS_INC_HYDROGEN\n%FORMAT(10I8)\n"
      format_print(fp, 10, "8d", bonds_h.map { |n|
        x = self.bonds[n]
        [x[0] * 3, x[1] * 3, self.bond_par(n).index + 1] }.flatten)
      
      fp.print "%FLAG BONDS_WITHOUT_HYDROGEN\n%FORMAT(10I8)\n"
      format_print(fp, 10, "8d", bonds_a.map { |n|
        x = self.bonds[n]
        [x[0] * 3, x[1] * 3, self.bond_par(n).index + 1] }.flatten)
      
      fp.print "%FLAG ANGLES_INC_HYDROGEN\n%FORMAT(10I8)\n"
      format_print(fp, 10, "8d", angles_h.map { |n|
        x = self.angles[n]
        [x[0] * 3, x[1] * 3, x[2] * 3, self.angle_par(n).index + 1] }.flatten)
      
      fp.print "%FLAG ANGLES_WITHOUT_HYDROGEN\n%FORMAT(10I8)\n"
      format_print(fp, 10, "8d", angles_a.map { |n|
        x = self.angles[n]
        [x[0] * 3, x[1] * 3, x[2] * 3, self.angle_par(n).index + 1] }.flatten)
      
      [dihedrals_h, dihedrals_a].each { |dihed|
        if dihed == dihedrals_h
          fp.print "%FLAG DIHEDRALS_INC_HYDROGEN\n"
        else
          fp.print "%FLAG DIHEDRALS_WITHOUT_HYDROGEN\n"
        end
        fp.print "%FORMAT(10I8)\n"
        format_print(fp, 10, "8d", dihed.map { |n|
          if n < 2000000
            x = self.dihedrals[n % 1000000]
            #  Note: if n >= 1000000, then the 1-4 interaction should be ignored to avoid double calculation. This is specified by the negative index for the third atom, but if the third atom is 0 we have a problem. To avoid this situation, the atom array is reversed.
            if n >= 1000000 && x[2] == 0
              x = x.reverse
            end
            k = self.dihedral_par(n % 1000000).index + 1
            [x[0] * 3, x[1] * 3, (n >= 1000000 ? -x[2] : x[2]) * 3, x[3] * 3, k]
          else
            x = self.impropers[n % 1000000]
            #  The improper torsion is specified by the negative index for the fourth atom, but if the fourth atom is 0 then we have a problem. To avoid this situation, the first and fourth atoms are exchanged.
            if x[3] == 0
              x = [x[3], x[1], x[2], x[0]]
            end
            k = self.improper_par(n % 1000000).index + 1 + par.ndihedrals
            [x[0] * 3, x[1] * 3, -x[2] * 3, -x[3] * 3, k]
          end
        }.flatten)
      }  #  end each [dihedral_h, dihedral_a]
      
      fp.print "%FLAG EXCLUDED_ATOMS_LIST\n%FORMAT(10I8)\n"
      format_print(fp, 10, "8d", extable)
    
      fp.print "%FLAG HBOND_ACOEF\n%FORMAT(5E16.8)\n"
      fp.print "\n"
    
      fp.print "%FLAG HBOND_BCOEF\n%FORMAT(5E16.8)\n"
      fp.print "\n"
    
      fp.print "%FLAG HBCUT\n%FORMAT(5E16.8)\n"
      fp.print "\n"
    
      fp.print "%FLAG AMBER_ATOM_TYPE\n%FORMAT(20a4)\n"
      format_print(fp, 20, "-4.4s", self.atoms.map { |ap| ap.atom_type })
    
      fp.print "%FLAG TREE_CHAIN_CLASSIFICATION\n%FORMAT(20a4)\n"
      format_print(fp, 20, "-4.4s", self.atoms.map { |ap| "M" })
    
      fp.print "%FLAG JOIN_ARRAY\n%FORMAT(10I8)\n"
      format_print(fp, 10, "8d", self.atoms.map { |ap| 0 })
    
      fp.print "%FLAG IROTAT\n%FORMAT(10I8)\n"
      format_print(fp, 10, "8d", self.atoms.map { |ap| 0 })
    
      fp.print "%FLAG RADIUS_SET\n%FORMAT(1a80)\n"
      fp.print "modified Bondi radii (mbondi)\n"
      
      fp.print "%FLAG RADII\n%FORMAT(5E16.8)\n"
      format_print(fp, 5, "16.8E", self.atoms.map { |ap|
        (ap.atomic_number == 1 ? 1.30 : 1.70) })
      
      fp.print "%FLAG SCREEN\n%FORMAT(5E16.8)\n"
      format_print(fp, 5, "16.8E", self.atoms.map { |ap|
        (ap.atomic_number == 1 ? 0.85 : 0.72) })
      
      if periodic
        fp.print "%FLAG SOLVENT_POINTERS\n%FORMAT(3I8)\n"
        fp.printf "%8d%8d%8d\n", last_solute, fragments.length, first_solv_mol
        
        fp.print "%FLAG ATOMS_PER_MOLECULE\n%FORMAT(10I8)\n"
        format_print(fp, 10, "8d", fragments)
        
        fp.print "%FLAG BOX_DIMENSIONS\n%FORMAT(5E16.8)\n"
        box = self.box
        fp.printf "%16.8E%16.8E%16.8E%16.8E\n", 0, box[0].length, box[1].length, box[2].length
      end
    }
    
    fname2 = fname.sub(/\.prmtop$/, '.inpcrd')
    fname2 += '.inpcrd' if fname == fname2
    
    open(fname2, "w") { |fp|
      fp.printf "%-80.80s\n", self.name
      fp.printf "%6d\n", self.natoms
      format_print(fp, 6, "12.7f", self.atoms.map { |ap| [ap.x, ap.y, ap.z] }.flatten )
	  if periodic
	    fp.printf "%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f\n", box[0].length, box[1].length, box[2].length, 90, 90, 90
	  end
    }
  
    return true

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
	ante_dir = MolbyResourcePath + "/amber11/bin"

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
	    open("resp.input#{i}", "wb") { |fp|
		  fp.print " #{self.name} \n"
		  fp.print " &cntrl\n"
		  fp.print " ioutopt=1, IQOPT=#{i}, nmol=1, ihfree=1, irstrnt=1, qwt=#{i*0.0005},\n"
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
	  open("respgen_in.ac", "wb") { |fp|
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
	  if !call_subprocess("\"#{ante_dir}/respgen\" -i respgen_in.ac -o resp.input1 -f resp1", "respgen (stage 1)") \
	  || !call_subprocess("\"#{ante_dir}/respgen\" -i respgen_in.ac -o resp.input2 -f resp2", "respgen (stage 2)")
	    error_message_box("Cannot run respgen.")
	    Dir.chdir(cwd)
	    return
	  end
	  hide_progress_panel
	end
	
	#  Create ESP file
	a2b = 1.8897259885   #  angstrom to bohr
	open("resp.esp", "wb") { |fp|
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
	if !call_subprocess("\"#{ante_dir}/resp\" -O -i resp.input1 -o resp.output1 -e resp.esp -t qout_stage1", "resp (stage 1)") \
	|| !call_subprocess("\"#{ante_dir}/resp\" -O -i resp.input2 -o resp.output2 -e resp.esp -q qout_stage1 -t qout_stage2", "resp (stage 2)")
	  error_message_box("Cannot run resp.")
	  Dir.chdir(cwd)
	  return 
	end
	hide_progress_panel
	
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
