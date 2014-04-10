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
	read_only = [:step, :coord_frame, :transient_temperature, :average_temperature]
	keys = [:timestep, :temperature, :cutoff, :electro_cutoff, :pairlist_distance,
	 :switch_distance, :scale14_vdw, :scale14_elect, :use_xplor_shift, :dielectric,
	 :andersen_freq, :andersen_coupling, :random_seed, :relocate_center,
	 :use_graphite, :minimize_cell,
	 :gradient_convergence, :coordinate_convergence,
	 :log_file, :coord_file, :vel_file, :force_file, :debug_file, :debug_output_level,
	 :coord_output_freq, :energy_output_freq,
	 :step, :coord_frame, :transient_temperature, :average_temperature]
	#  Arrange the items in vertical direction
	n = (keys.count + 1) / 2
	i = 0
	keys = keys.sort_by { |k| i += 1; (i > n ? i - n + 0.5 : i) }
	#  Do dialog
    hash = Dialog.run("MM/MD Advanced Settings", "Close", nil) {
	  items = []
	  keys.each { |k|
	    enabled = !read_only.include?(k)
	    items.push(item(:text, :title=>k.to_s + (enabled ? "" : "*")))
		it = item(:textfield, :width=>120, :value=>arena[k].to_s)
		if enabled
		  it[:tag] = k
		else
		  it[:enabled] = false
		end
		items.push(it)
	  }
	  items.push(item(:text, :title=>"(*: read-only parameters)"), -1, -1, -1)
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
	  mol = self
	  h = Dialog.run("MM/MD Parameter Missing", nil, nil) {
		layout(1,
		  item(:text, :title=>"Some MM/MD parameters are missing.\nPlease select one of the following options."),
		  item(:button, :title=>"Auto Guess", :align=>:center,
		    :action=>lambda { |it| end_modal(2) } ),
		  item(:button, :title=>"Run Antechamber Manually...", :align=>:center,
		    :action=>lambda { |it| end_modal(3) } ),
		  item(:button, :title=>"Cancel", :align=>:center,
		    :action=>lambda { |it| end_modal(1) } ))
	  }
	  status = h[:status]
	  if status == 2
	    n = invoke_antechamber(false)
		return nil if n != 0
	  elsif status == 3
	    n = invoke_antechamber(true)
		return nil if n != 0
	  else
	    return nil
	  end
    end

    #  Initialize some fields at first invocation
    if !@arena_created
	  # if arena.log_file == nil && self.dir != nil
	  #  arena.log_file = self.name.gsub(/\.\w+$/, ".log")
	  # end
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
		     :action=>lambda { |it|
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
	f = hash["log_file"].strip
	arena.log_file = (f == "" ? nil : f)
	if hash[:status] == 0
	  return Integer(hash["number_of_frames"])
	else
	  return -1
	end
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
  
  def ambertools_dialog(tool, msg = nil, block = nil)
	log_dir = get_global_settings("antechamber.log_dir")
	if !log_dir
	  log_dir = $home_directory + "/Molby/antechamber"
	end
	if $platform == "win"
	  suffix = ".exe"
	else
	  suffix = ""
	end
	log_level = (get_global_settings("antechamber.log_level") || "none")
	if msg == nil
	  if tool == "antechamber"
	    msg = "Auto Guess MM/MD Parameters (Antechamber)"
	  else
	    msg = "Run " + tool.capitalize
	  end
    end
    hash = Dialog.run(msg) {
	  @toolname = tool + suffix
      def valid_antechamber_dir(s)
	    FileTest.exist?(s + "/" + @toolname)
  	  end
	  layout(2,
	    (block ?
		  layout(3,
		    item(:text, :title=>"Atoms to process: "),
			item(:text, :title=>block, :width=>120),
			item(:button, :title=>"Skip this block", :action=>lambda { |it| end_modal(2) } )) :
		  -1),
		-1,
		item(:text, :title=>"Net Molecular Charge:"),
		item(:textfield, :width=>"80", :tag=>"nc", :value=>(get_global_settings("antechamber.nc") || "0")),
		(tool == "resp" ?
		  -1 :
		  item(:checkbox, :title=>"Calculate partial charges", :tag=>"calc_charge",
		    :action=>lambda { |it|
			  set_attr("optimize_structure", :enabled=>(it[:value] != 0)) } )),
		-1,
		(tool == "resp" ?
		  -1 :
		  layout(2,
		    item(:view, :width=>"12"),
			item(:checkbox, :title=>"Optimize structure before calculating charges (may be slow)",
			  :tag=>"optimize_structure",
			  :value=>(get_global_settings("antechamber.optimize_structure") || 0)))),
		-1,
		(tool == "resp" ?
		  -1 :
		  item(:checkbox, :title=>"Guess atom types", :tag=>"guess_atom_types",
		    :value=>(get_global_settings("antechamber.guess_atom_types") || 1))),
		-1,
		item(:checkbox, :title=>"Use the residue information for connection analysis", :tag=>"use_residue",
		  :value=>(get_global_settings("antechamber.use_residue") || 0),
		  :action=>lambda { |it| valid_antechamber_dir(it[:value]) && set_attr(0, :enabled=>true) } ),
		-1,
		item(:line),
		-1,
		item(:text, :title=>"Log directory:"),
		[ item(:button, :title=>"Choose...",
			:action=>lambda { |it|
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
	  it = item_with_tag("calc_charge")
	  if it
	    it[:value] = (get_global_settings("antechamber.calc_charge") || 0)
		it[:action].call(it)
	  end
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
	return 1 - hash[:status]  #  1: OK, 0: Cancel, -1: Skip
  end
  
  def create_ante_log_dir(name, key)
    log_dir = get_global_settings("antechamber.log_dir")
	msg = ""
	begin
	  if !FileTest.directory?(log_dir)
	    mkdir_recursive(log_dir) rescue ((msg = "Cannot create directory #{log_dir}") && raise)
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
  
  def invoke_antechamber(ask_options = true, msg = nil)
    #  Find the ambertool directory
	ante_dir = "#{ResourcePath}/amber11/bin"
	#  Ask for antechamber options and log directory
	if ask_options
  	  n = ambertools_dialog("antechamber", msg)
	  return 1 if n == 0
	end
	nc = get_global_settings("antechamber.nc").to_i
	calc_charge = get_global_settings("antechamber.calc_charge").to_i
	guess_atom_types = get_global_settings("antechamber.guess_atom_types").to_i
	optimize_structure = get_global_settings("antechamber.optimize_structure").to_i
	use_residue = get_global_settings("antechamber.use_residue").to_i
	#  Create log directory
	name = (self.name || "unknown").sub(/\.\w*$/, "").sub(/\*/, "")  #  Remove the extension and "*"
	log_dir = get_global_settings("antechamber.log_dir")
	if log_dir == nil
	  log_dir = document_home + "/Molby/antechamber"
	end
	if !File.directory?(log_dir)
	  mkdir_recursive(log_dir)
	end
	tdir = nil
	1000.times { |i|
	  tdir = sprintf("%s/%s.%04d", log_dir, name, i)
	  if !File.exists?(tdir) && (Dir.mkdir(tdir) == 0)
	    break
	  end
	  tdir = nil
	}
	if tdir == nil
	  error_message_box("Cannot create log directory in #{log_dir}.")
	  return -1
	end
	cwd = Dir.pwd
	Dir.chdir(tdir)
	mol2 = self.dup
	if use_residue == 0
	  mol2.assign_residue(mol2.all, 1)
	end
	#  Rename the molecule (antechamber assumes the atom names begin with element symbol)
	mol2.each_atom { |ap|
	  ap.name = ap.element
	  if ap.name.length == 1
	    ap.name += sprintf("%03d", ap.index % 1000)
	  else
	    ap.name += sprintf("%02d", ap.index % 100)
	  end
	}
	mol2.savepdb("./mol.pdb")
	#  Set environmental variable
	p = ENV["AMBERHOME"]
	if p == nil
	  amberhome = "#{ResourcePath}/amber11"
	  if $platform == "win"
	    amberhome.gsub!("/", "\\")
	  end
	  ENV["AMBERHOME"] = amberhome
	end
	if calc_charge != 0
	  opt = "-nc #{nc} -c bcc"
	  if optimize_structure == 0
	    opt += " -ek 'maxcyc=0'"
	  end
	else
	  opt = ""
	end
	if guess_atom_types == 0
	  opt += " -j 0"
    end
	n = call_subprocess("#{ante_dir}/antechamber -i mol.pdb -fi pdb -o mol.ac -fo ac #{opt}", "Antechamber")
	if n != 0
	  error_message_box("Antechamber failed: status = #{n}.")
	  Dir.chdir(cwd)
	  return n
	else
	  if guess_atom_types != 0
	    n = call_subprocess("#{ante_dir}/parmchk -i mol.ac -f ac -o frcmod", "Parmchk")
	    if n != 0
	      error_message_box("Parmchk failed: status = #{n}.")
		  Dir.chdir(cwd)
		  return n
	    end
	  end
	end
	Dir.chdir(cwd)
	n = import_ac("#{tdir}/mol.ac", calc_charge, guess_atom_types)
	if n != 0
	  error_message_box("Cannot import antechamber output.")
	  return n
	end
	if calc_charge != 0
	  n = import_sqmout("#{tdir}/sqm.out")
	  if n != 0
	    error_message_box("Cannot import sqm output.")
		return n
	  end
	end
	if guess_atom_types != 0
	  n = import_frcmod("#{tdir}/frcmod")
	  if n != 0
	    error_message_box("Cannot import parmchk output.")
	    return n
	  end
	  if self.nimpropers > 0
	    remove_improper(IntGroup[0...self.nimpropers])
	  end
	end
	log_level = get_global_settings("antechamber.log_level")
	log_keep_number = get_global_settings("antechamber.log_keep_number")
	erase_old_logs(tdir, log_level, log_keep_number)
	return 0
	
  end
  
  def import_ac(acfile, read_charge, read_type)
    open(acfile, "r") { |fp|
	  while (s = fp.gets)
	    next if s !~ /^ATOM/
		s.chomp!
		idx = Integer(s[4..11]) - 1
		ap = atoms[idx]
		if read_charge != 0
		  charge = Float(s[54..63])
  		  ap.charge = charge
		end
		if read_type != 0
		  type = s[72..-1]
		  type.gsub!(/ /, "")
		  ap.atom_type = type
		end
	  end
	}
	return 0
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
	return 0
  end
  
  def import_frcmod(file)
    self.md_arena.prepare(true)  #  Clean up existing parameters
	par = self.parameter
	if !FileTest.exist?(file)
	  error_message_box("#{file} does not exist?? really??")
	end
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
		    types = s[0..4]
			k, r0, com = s[5..-1].split(" ", 3)
		    pp = par.bonds.lookup(types, :local, :missing) || par.bonds.insert
			pp.atom_types = types
			pp.k = k
			pp.r0 = r0
			pp.comment = com
		  when 3
		    types = s[0..7]
			k, a0, com = s[8..-1].split(" ", 3)
			pp = par.angles.lookup(types, :local, :missing) || par.angles.insert
			pp.atom_types = types
			pp.k = k
			pp.a0 = a0
			pp.comment = com
		  when 4
		    types = s[0..10]
			n, k, phi0, period, com = s[11..-1].split(" ", 5)
			pp = par.dihedrals.lookup(types, :local, :missing) || par.dihedrals.insert
			pp.atom_types = types
			pp.mult = 1
			pp.k = k
			pp.phi0 = phi0
			pp.period = Float(period).round
			pp.comment = com
		  when 5
		    types = s[0..10]
		    k, phi0, period, com = s[11..-1].split(" ", 4)
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
	return 0
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
  
  def Molecule.cmd_import_amberlib(mol)
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
  def export_prmtop(filename = nil)
  
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
          eps = Math.sqrt_safe(p1.eps * p2.eps)
          d = p1.r_eq + p2.r_eq
        end
        vdw_access_table[i * nvdw_pars + j] =  vdw_access_table[j * nvdw_pars + i] = vdw_pair_table.length
        vdw_pair_table.push([(d ** 12) * eps, 2 * (d ** 6) * eps])
      }
    }
        
	if filename
	  fname = filename
    else
      basename = (self.path ? File.basename(self.path, ".*") : self.name)
      fname = Dialog.save_panel("AMBER prmtop/inpcrd file name", self.dir, basename + ".prmtop", "AMBER prmtop (*.prmtop)|*.prmtop|All files|*.*")
      return nil if !fname
	  hash = Dialog.run("Select prmtop format") {
	    layout(2,
		  item(:text, :title=>"Select prmtop format:", :tag=>"prmtopname"),
		  item(:popup, :subitems=>["AMBER 8/NAMD", "AMBER 11"], :tag=>"ambertype"))
	    set_attr("ambertype", :value=>1)
	  }
	  if hash[:status] == 0
	    case hash["ambertype"]
		when 0
		  ambertype = 8
		when 1
		  ambertype = 11
	    end
	  else
	    return nil
      end
    end

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
      format_print(fp, 10, "8d", (0...self.natoms).map { |i| par.vdws.lookup(i).index + 1 })
      
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
      
	  if ambertype > 8
        fp.print "%FLAG SCEE_SCALE_FACTOR\n%FORMAT(5E16.8)\n"
        format_print(fp, 5, "16.8E", par.dihedrals.map { |p| 1.2 } + par.impropers.map { |p| 0.0 })
      
        fp.print "%FLAG SCNB_SCALE_FACTOR\n%FORMAT(5E16.8)\n"
        format_print(fp, 5, "16.8E", par.dihedrals.map { |p| 2.0 } + par.impropers.map { |p| 0.0 })
      end
	  
      fp.print "%FLAG SOLTY\n%FORMAT(5E16.8)\n"
      format_print(fp, 5, "16.8E", (0...par.nvdws).map { 0.0 } )
    
      fp.print "%FLAG LENNARD_JONES_ACOEF\n%FORMAT(5E16.8)\n"
      format_print(fp, 5, "16.8E", vdw_pair_table.map { |x| x[0] })
    
      fp.print "%FLAG LENNARD_JONES_BCOEF\n%FORMAT(5E16.8)\n"
      format_print(fp, 5, "16.8E", vdw_pair_table.map { |x| x[1] })
    
      fp.print "%FLAG BONDS_INC_HYDROGEN\n%FORMAT(10I8)\n"
      format_print(fp, 10, "8d", bonds_h.map { |n|
        x = self.bonds[n]
        [x[0] * 3, x[1] * 3, par.bonds.lookup(x).index + 1] }.flatten)
      
      fp.print "%FLAG BONDS_WITHOUT_HYDROGEN\n%FORMAT(10I8)\n"
      format_print(fp, 10, "8d", bonds_a.map { |n|
        x = self.bonds[n]
        [x[0] * 3, x[1] * 3, par.bonds.lookup(x).index + 1] }.flatten)
      
      fp.print "%FLAG ANGLES_INC_HYDROGEN\n%FORMAT(10I8)\n"
      format_print(fp, 10, "8d", angles_h.map { |n|
        x = self.angles[n]
        [x[0] * 3, x[1] * 3, x[2] * 3, par.angles.lookup(x).index + 1] }.flatten)
      
      fp.print "%FLAG ANGLES_WITHOUT_HYDROGEN\n%FORMAT(10I8)\n"
      format_print(fp, 10, "8d", angles_a.map { |n|
        x = self.angles[n]
        [x[0] * 3, x[1] * 3, x[2] * 3, par.angles.lookup(x).index + 1] }.flatten)
      
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
            k = par.dihedrals.lookup(x).index + 1
            [x[0] * 3, x[1] * 3, (n >= 1000000 ? -x[2] : x[2]) * 3, x[3] * 3, k]
          else
            x = self.impropers[n % 1000000]
            #  The improper torsion is specified by the negative index for the fourth atom, but if the fourth atom is 0 then we have a problem. To avoid this situation, the first and fourth atoms are exchanged.
            if x[3] == 0
              x = [x[3], x[1], x[2], x[0]]
            end
            k = par.impropers.lookup(x).index + 1 + par.ndihedrals
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
      
	  if ambertype > 8
        fp.print "%FLAG RADIUS_SET\n%FORMAT(1a80)\n"
        fp.print "modified Bondi radii (mbondi)\n"
      
        fp.print "%FLAG RADII\n%FORMAT(5E16.8)\n"
        format_print(fp, 5, "16.8E", self.atoms.map { |ap|
          (ap.atomic_number == 1 ? 1.30 : 1.70) })

        fp.print "%FLAG SCREEN\n%FORMAT(5E16.8)\n"
        format_print(fp, 5, "16.8E", self.atoms.map { |ap|
          (ap.atomic_number == 1 ? 0.85 : 0.72) })
      end
	  
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
  
  def cmd_antechamber
    if natoms == 0
	  error_message_box "Molecule is empty"
	  return
	end
	return if invoke_antechamber(true) != 0
	message_box("Antechamber succeeded.", "Antechamber Success", :ok)
  end
  
  def cmd_run_resp
    if natoms == 0
	  error_message_box "Molecule is empty"
	  return
	elsif nelpots == 0
	  error_message_box "No ESP information is loaded"
	  return
	end
	return if ambertools_dialog("resp") == 0
	nc = get_global_settings("antechamber.nc")
	ante_dir = Molby::ResourcePath + "/amber11/bin"

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
	  status = call_subprocess("\"#{ante_dir}/respgen\" -i respgen_in.ac -o resp.input1 -f resp1", "respgen (stage 1)")
	  if status == 0
		status = call_subprocess("\"#{ante_dir}/respgen\" -i respgen_in.ac -o resp.input2 -f resp2", "respgen (stage 2)")
	  end
	  if status != 0
	    error_message_box("Cannot run respgen: status = #{status}")
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
	status = call_subprocess("\"#{ante_dir}/resp\" -O -i resp.input1 -o resp.output1 -e resp.esp -t qout_stage1", "resp (stage 1)")
	if status == 0
	  status = call_subprocess("\"#{ante_dir}/resp\" -O -i resp.input2 -o resp.output2 -e resp.esp -q qout_stage1 -t qout_stage2", "resp (stage 2)")
	  if status == 255 && File.exist?("punch")
		status = 0   #  Ignore error at the second stage
	  end
	end
	if status != 0
	  error_message_box("Cannot run resp: status = #{status}")
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
		  :action=>lambda { |it|
		    esp_save = get_global_settings("gamess.esp")
			set_global_settings("gamess.esp", 1)
			mol.cmd_create_gamess_input
			set_global_settings("gamess.esp", esp_save)
		  }),
		  {:align=>:right}],
		item(:line),
		item(:text, :title=>"Step 2:\nImport GAMESS .dat file"),
		[item(:button, :title=>"Import GAMESS dat...",
		  :action=>lambda { |it|
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
		  :action=>lambda { |it|
		    if mol.cmd_run_resp
			  if get_global_settings("antechamber.log_level") == "latest"
			    mol.clean_ante_log_dir(Integer(get_global_settings("antechamber.log_keep_number")))
			  end
			  end_modal(0)
			end
		  }),
		  {:align=>:right}],
		item(:line),
		[item(:button, :title=>"Close", :action=>lambda { |it| end_modal(1) }),
		  {:align=>:center}]
	  )
	  if mol.nelpots == 0
	    set_attr("resp", :enabled=>false)
	  end
	}
  end
  
  def cmd_edit_local_parameter_in_mainview(ptype, indices, names, types, value, partypes, params)
    #  The parameters are given as space separated strings
	#  e.g. "1.862 200.000"
	mol = self
	case ptype
	when "bond"
	  k = ["k", "r0"]
	  pen = self.parameter.bonds
	when "angle"
	  k = ["k", "a0"]
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
	    item(:text, :title=>"Edit #{ptype} parameter for #{indices} (#{names}, #{types})"),
		item(:text, :title=>"(Current value = #{value})"),
	    layout(4,
	      [item(:text, :title=>"types"), {:align=>:center}],
	      [item(:text, :title=>k[0]), {:align=>:center}],
	      [item(:text, :title=>k[1]), {:align=>:center}],
	      (k[2] ? [item(:text, :title=>k[2]), {:align=>:center}] : -1),
		
		  item(:textfield, :width=>100, :value=>partypes, :tag=>"types"),
		  item(:textfield, :width=>100, :value=>p[0], :tag=>k[0]),
		  item(:textfield, :width=>100, :value=>p[1], :tag=>k[1]),
		  (k[2] ? item(:textfield, :width=>100, :value=>p[2], :tag=>k[2]) : -1)
		),
		(ptype == "bond" || ptype == "angle" ?
		  item(:button, :title=>"Guess k by UFF...", :align=>:right,
		    :action=>lambda { |it|
			  guess, cval = mol.guess_uff_parameter_dialog(value(k[1]), indices)
			  if guess
			    set_value("k", guess)
				if (ptype == "angle")
				  set_value("a0", cval)
				end
			  end
			}) : 
	      nil)
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
  
  #  Solvate the molecule with the given solvent box.
  #  The first argument (box) must be a Molecule containing a unit cell information.
  #  The second argument defines the size of the solvated system. A positive number represents
  #  an offset to the bounding box of the solute, and a negative number represents an absolute size.
  #  If it is given as an Array (or a Vector), the elements define the sizes in the x/y/z directions.
  #  The third argument represents the limit distance to avoid conflict between the solute and
  #  solvent. The solvent molecule containing atoms within this limit from the solute is removed.
  #  The atom group containing added solvent molecule is returned.
  def solvate(sbox, size = [10.0, 10.0, 10.0], limit = 3.0)

	#  Sanity check
    if sbox.box == nil
	  raise MolbyError, "the solvent box does not have a unit cell information"
    end
	flags = sbox.box[4]
	if flags[0] == 0 || flags[1] == 0 || flags[2] == 0
	  raise MolbyError, "the solvent box does not have three-dimensional periodicity"
	end

	show_progress_panel("Adding a solvent box...")

	#  Calculate the box size
	b = self.bounds
	bsize = b[1] - b[0]
	if size.kind_of?(Numeric)
	  size = [size, size, size]
	end
	size = size.collect { |s| Float(s) }
	limit = Float(limit)
	(0..2).each do |i|
	  d = (size[i] >= 0 ? bsize[i] + size[i] * 2 : -size[i])
	  if d < bsize[i]
	    raise MolbyError, "the box size is too small"
	  end
	  bsize[i] = d;
	end
#	puts "Box size = #{bsize}"

	#  Translate the solute to the center of the box
	translate((b[0] + b[1]) * -0.5)
	solute_natoms = self.natoms

	#  Add solvents so that the target box is fully covered
	set_progress_message("Duplicating the solvent box...")
	rtr = sbox.cell_transform.inverse
	min = Vector3D[1e30, 1e30, 1e30]
	max = Vector3D[-1e30, -1e30, -1e30]
	[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1]].each do |pt|
	  pt = Vector3D[(pt[0] - 0.5) * bsize[0], (pt[1] - 0.5) * bsize[1], (pt[2] - 0.5) * bsize[2]]
	  rpt = rtr * pt
	  min.x = rpt.x if min.x > rpt.x
	  min.y = rpt.y if min.y > rpt.y
	  min.z = rpt.z if min.z > rpt.z
	  max.x = rpt.x if max.x < rpt.x
	  max.y = rpt.y if max.y < rpt.y
	  max.z = rpt.z if max.z < rpt.z
	end
	xmin = (min.x + 0.5).floor
	xmax = (max.x + 0.5).floor
	ymin = (min.y + 0.5).floor
	ymax = (max.y + 0.5).floor
	zmin = (min.z + 0.5).floor
	zmax = (max.z + 0.5).floor
	sbox_natoms = sbox.natoms
	xv, yv, zv, ov, flags = sbox.box
#	puts "xmin = #{xmin}, ymin = #{ymin}, zmin = #{zmin}"
#	puts "xmax = #{xmax}, ymax = #{ymax}, zmax = #{zmax}"
    newbox = Molecule.new
	(xmin..xmax).each do |x|
	  (ymin..ymax).each do |y|
	    (zmin..zmax).each do |z|
		  newbox.add(sbox)
		  newbox.translate(xv * x + yv * y + zv * z, IntGroup[newbox.natoms - sbox_natoms..newbox.natoms - 1])
		end
	  end
	end
	
	#  Remove out-of-bounds molecules
	set_progress_message("Removing out-of-bounds molecules...")
	g = newbox.atom_group do |ap|
	  r = ap.r
	  r.x < -(bsize[0] - limit) * 0.5 || r.y < -(bsize[1] - limit) * 0.5 || r.z < -(bsize[2] - limit) * 0.5 ||
	  r.x > (bsize[0] - limit) * 0.5 || r.y > (bsize[1] - limit) * 0.5 || r.z > (bsize[2] - limit) * 0.5
	end
	g = newbox.fragment(g)  #  expand by fragment
	newbox.remove(g)
	
	#  Add solvent molecules
	self.line_mode(true)
	self.add(newbox)
#	puts "Removed atoms by bounds: #{g}"
	
	#  Find conflicts
	set_progress_message("Removing conflicting molecules...")
	conf = find_conflicts(limit, IntGroup[0..solute_natoms - 1], IntGroup[solute_natoms..self.natoms - 1])
	g = atom_group(conf.map { |c| c[1] } )    #  atom group containing conflicting atoms
	g = fragment(g)                           #  expand by fragment
	remove(g)
#	puts "Removed atoms by conflicts: #{g}"
	
	#  Renumber residue numbers of solvent molecules
	rseq = max_residue_number(0..solute_natoms - 1)
	rseq = 0 if rseq == nil
	each_fragment do |g|
	  next if g[0] < solute_natoms
	  rseq += 1
	  assign_residue(g, atoms[g[0]].res_name + ".#{rseq}")
	end

	#  Set the segment number and name (SOLV)
	seg = (0..solute_natoms - 1).map { |n| self.atoms[n].seg_seq }.max + 1
	each_atom(solute_natoms..self.natoms - 1) { |ap|
	  ap.seg_seq = seg if ap.seg_seq < seg
	  ap.seg_name = "SOLV"
	}

    #  Set the unit cell information
	set_box(bsize[0], bsize[1], bsize[2])
	
	hide_progress_panel
	resize_to_fit
	
	return IntGroup[solute_natoms..self.natoms - 1]

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
	  Dialog.run(" ") {
	    layout(1,
		  item(:text, :title=>"Please open a molecule file containing a solvent box."))
	  }
	  return
	end
	hash = Dialog.run("Solvate") {
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
     
end

register_menu("MM/MD\tImport AMBER Lib...", :cmd_import_amberlib)
register_menu("MM/MD\t-", nil)
register_menu("MM/MD\tGuess MM/MD Parameters...", :cmd_antechamber, :non_empty)
register_menu("MM/MD\tGuess UFF Parameters...", :guess_uff_parameters, :non_empty) # uff.rb
register_menu("MM/MD\tGAMESS and RESP...", :cmd_gamess_resp, :non_empty)
register_menu("MM/MD\tCreate SANDER Input...", :export_prmtop, :non_empty)
register_menu("MM/MD\tImport AMBER Frcmod...", :cmd_import_frcmod, :non_empty)
register_menu("MM/MD\t-", nil)
register_menu("MM/MD\tSolvate...", :cmd_solvate, :non_empty)
