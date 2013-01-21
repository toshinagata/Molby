#
#  gamess.rb
#
#  Created by Toshi Nagata on 2009/11/22.
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

  def Molecule.read_gamess_basis_sets(fname)
    # $gamess_basis = Hash.new unless $gamess_basis
	$gamess_ecp = Hash.new unless $gamess_ecp
	# $gamess_basis_desc = Hash.new unless $gamess_basis_desc
	# $gamess_basis_keys = [] unless $gamess_basis_keys
	basename = File.basename(fname, ".*")
	keys = []
	descname = nil
    File.open(fname, "r") { |fp|
	  while (s = fp.gets)
	    if s =~ /^\s*!\s*/ && descname == nil
	      #  Get the descriptive name from the first comment line
		  s = Regexp.last_match.post_match
		  if s =~ /  EMSL/
		    descname = Regexp.last_match.pre_match
		  else
		    descname = (s.split)[0]
		  end
		  $gamess_basis_desc[basename] = descname
		  next
		end
		ss, bas = (s.split)[0..1]     #  Tokens delimited by whitespaces
		next if ss == nil || ss == ""
		if ss == "$ECP"
		#  Read the ECP section
		  ecp = $gamess_ecp[basename]
		  if !ecp
		    ecp = []
			$gamess_ecp[basename] = ecp
		  end
		  keys << basename unless keys.include?(basename)
		  while (s = fp.gets)
		    break if s=~ /\$END/
			ary = s.split  #  (PNAME, PTYPE, IZCORE, LMAX+1)
			elem = ary[0].match(/\w{1,2}/).to_a[0]  #  ary[0] should begin with an element name
			ap = Parameter.builtin.elements.find { |p| p.name.casecmp(elem) == 0 }
			raise MolbyError, "the ECP definition does not begin with an element name: #{s}" if !ap
			ecpdef = s
			ln = 1
			(0..Integer(ary[3])).each {
			  s = fp.gets
			  raise MolbyError, "the ECP definition ends unexpectedly at line #{ln} for: #{s}" if !s
			  ln += 1
			  ary = s.split  #  (NGPOT, comment)
			  ecpdef += s
			  (1..Integer(ary[0])).each {
			    s = fp.gets
				raise MolbyError, "the ECP definition ends unexpectedly at line #{ln} for: #{s}" if !s
				ln += 1
				ecpdef += s
			  }
			}
		    ecp[ap.index] = ecpdef
		  end
		elsif ss =~ /\W/
		  #  Comments or other unrecognizable lines
		  next
		elsif (ap = Parameter.builtin.elements.find { |p| p.name.casecmp(ss) == 0 || p.fullname.casecmp(ss) == 0 })
		  #  Valid basis definition
		  if bas == nil || bas =~ /\W/
		    bas = basename
		  end
		  basis = $gamess_basis[bas]
		  if !basis
		    basis = []
			$gamess_basis[bas] = basis
		  end
		  keys << bas unless keys.include?(bas)
		  basdef = ""
		  while (s = fp.gets) && s =~ /\S/
		    basdef += s
		  end
		  basis[ap.index] = basdef
		else
		  raise MolbyError, "ss is not a valid element symbol or name: #{s}"
		end
	  end
    }
	unless $gamess_basis_keys.include?(basename)
	  $gamess_basis_keys.push(basename)
	end
  end
  
  def export_gamess(fname, hash)

    def reorder_array(ary, ordered_sub_ary)
	  return ordered_sub_ary + (ary - ordered_sub_ary)
	end
	
	now = Time.now.to_s
	basename = File.basename(fname, ".*")
	
	#  Various settings
	icharg = hash["charge"]
	mult = hash["mult"]
	runtyp = ["ENERGY", "PROP", "OPTIMIZE"][Integer(hash["runtype"])]
	scftyp = ["RHF", "ROHF", "UHF"][Integer(hash["scftype"])]
	bssname = hash["basis"]
	bssname2 = hash["secondary_basis"]
	if hash["use_secondary_basis"] != 0 && bssname2 != bssname
	  use_2nd = true
	  element2 = hash["secondary_elements"].split(/[\s,]+/).map { |name| name.capitalize }
	else
	  use_2nd = false
	  bssname2 = bssname
	end
	basis = $gamess_basis[bssname]
	basis2 = $gamess_basis[bssname2]
	if !basis || !basis2
	  raise MolbyError, "Unknown basis set name??? \"#{bssname}\" or \"#{bssname2}\""
	end

	#  Use effective core potentials?
	ecp = $gamess_ecp[bssname]
	ecp2 = $gamess_ecp[bssname2]
	ecp_read = (ecp || ecp2 ? "ECP=READ " : "")

	#  Use only one built-in basis set?
	gbasis = nil
	if !use_2nd
	  case bssname
	  when "PM3"
		gbasis = "PM3"
	  when "STO3G"
		gbasis = "STO NGAUSS=3"
	  when "321G"
		gbasis = "N21 NGAUSS=3"
	  when "631G"
		gbasis = "N31 NGAUSS=6"
	  when "631Gd"
		gbasis = "N31 NGAUSS=6 NDFUNC=1"
	  when "631Gdp"
		gbasis = "N31 NGAUSS=6 NDFUNC=1 NPFUNC=1"
	  when "6311G"
		gbasis = "N311 NGAUSS=6"
	  when "6311Gdp"
		gbasis = "N311 NGAUSS=6 NDFUNC=1 NPFUNC=1"
	  end
	end
    
	#  Count non-dummy atoms
	natoms = 0
	each_atom { |ap|
	  natoms += 1 if ap.atomic_number != 0
	}
	
	#  Fill hash with default values
	h = (hash["CONTRL"] ||= Hash.new)
	h["COORD"] ||= "UNIQUE"
	h["EXETYP"] ||= "RUN"
	h["ICHARG"] ||= icharg.to_s
	h["ICUT"] ||= "20"
	h["INTTYP"] ||= "HONDO"
	h["ITOL"] ||= "30"
	h["MAXIT"] ||= "200"
	h["MOLPLT"] ||= ".T."
	h["MPLEVL"] ||= "0"
	h["MULT"] ||= mult.to_s
	h["QMTTOL"] ||= "1e-08"
	h["RUNTYP"] ||= runtyp
	if hash["dft"] != 0
	  h["DFTTYP"] ||= hash["dfttype"]
	end
	if hash["use_internal"] != 0 && (hash["runtype"] == 2 || h["RUNTYP"] == "OPTIMIZE")
	  nzvar = natoms * 3 - 6  #  TODO: 3N-5 for linear molecules
	  h["NZVAR"] = nzvar.to_s
	else
	  nzvar = 0
	end
	h["SCFTYP"] ||= scftyp
	if ecp_read != ""
	  h["ECP"] ||= "READ"
	end
	h["UNITS"] = "ANGS"
	
	h = (hash["SCF"] ||= Hash.new)
	h["CONV"] ||= "1.0E-06"
	h["DIRSCF"] ||= ".T."
	h["FDIFF"] ||= ".T."
	h["DAMP"] ||= ".T."
	
	h = (hash["STATPT"] ||= Hash.new)
	h["NSTEP"] ||= "400"
	h["OPTTOL"] ||= "1.0E-06"

	h = (hash["SYSTEM"] ||= Hash.new)
	h["MEMDDI"] ||= "0"
	h["MWORDS"] ||= "16"
	h["TIMLIM"] ||= "50000"
	
	h = (hash["GUESS"] ||= Hash.new)
	h["GUESS"] ||= "HUCKEL"
	
	if gbasis
	  h = (hash["BASIS"] ||= Hash.new)
	  h["GBASIS"] ||= gbasis
	end
	
	if nzvar > 0
	  h = (hash["ZMAT"] ||= Hash.new)
	  h["DLC"] ||= ".T."
	  h["AUTO"] ||= ".T."
	end
	
	if hash["esp"] != 0
	  h = (hash["ELPOT"] ||= Hash.new)
	  h["IEPOT"] ||= "1"
	  h["OUTPUT"] ||= "PUNCH"
	  h["WHERE"] ||= "PDC"
	  h = (hash["PDC"] ||= Hash.new)
	  h["CONSTR"] ||= "NONE"
	  h["PTSEL"] ||= "CONNOLLY"
	end
	
    File.open(fname, "wb") { |fp|
	  fp.print "!  GAMESS input\n"
	  fp.print "!  Generated by Molby at #{now}\n"
	  fp.print "!  Basis set: " + $gamess_basis_desc[bssname] + "\n"
	  if use_2nd
	    fp.print "!  [" + element2.join(", ") + "]: " + $gamess_basis_desc[bssname2] + "\n"
	  end
	  controls = reorder_array(hash.keys.select { |k| hash[k].is_a?(Hash) },
		["CONTRL", "SCF", "STATPT", "SYSTEM", "GUESS", "BASIS", "ZMAT", "ELPOT", "PDC"])
	  controls.each { |k|
	    h = hash[k]
		next if h == nil || h.size == 0
		if k == "CONTRL"
		  ordered = ["COORD", "EXETYP", "ICHARG", "ICUT", "INTTYP", "ITOL", "MAXIT", "MOLPLT", "MPLEVL",
		             "MULT", "QMTTOL", "RUNTYP", "SCFTYP", "ECP", "UNITS", "DFTTYP", "NZVAR"]
		elsif k == "SCF"
		  ordered = ["CONV", "DIRSCF", "FDIFF", "DAMP"]
		elsif k == "STATPT"
		  ordered = ["NSTEP", "OPTTOL"]
		elsif k == "SYSTEM"
		  ordered = ["MEMDDI", "MWORDS", "TIMLIM"]
		elsif k == "GUESS"
		  ordered = ["GUESS"]
		elsif k == "BASIS"
		  ordered = ["GBASIS"]
		elsif k == "ZMAT"
		  ordered = ["DLC", "AUTO"]
		elsif k == "ELPOT"
		  ordered = ["IEPOT", "OUTPUT", "WHERE"]
		elsif k == "PDC"
		  ordered = ["CONSTR", "PTSEL"]
		else
		  ordered = []
		end
	    keys = reorder_array(h.keys, ordered)
		n = 0
		keys.each_with_index { |kk, i|
		  v = h[kk]
		  next if v == nil || v == ""
		  if n == 0
		    fp.printf " $%-6s", k
		  elsif n % 3 == 0
		    fp.print "\n        "
		  end
		  fp.printf " %s=%s", kk, h[kk].to_s
		  n += 1
		}
		if n > 0
		  fp.print " $END\n"
		end
	  }
	  fp.print " $DATA\n#{basename}\nC1 0\n"
	  secondary = []
	  each_atom { |ap|
	    next if ap.atomic_number == 0
		fp.printf "%-6s %4d %10.6f %10.6f %10.6f\n", ap.name, ap.atomic_number, ap.r.x, ap.r.y, ap.r.z
		if use_2nd && element2.include?(ap.element)
		  secondary[ap.index] = true
		end
	    if !gbasis
		  #  Basis specification followed by a blank line
		  bas = (secondary[ap.index] ? basis2 : basis)
		  if bas.is_a?(Array)
		    bas = bas[ap.atomic_number]
		  end
		  if !bas
		    raise MolbyError, "Basis set is not defined for atom #{ap.index}, element #{ap.element}"
		  end
		  fp.print bas
		  fp.print "\n"
		end
	  }
	  fp.print " $END\n"
	  ecp_ary = []
	  if ecp || ecp2
	    fp.print " $ECP\n"
		each_atom { |ap|
		  an = ap.atomic_number
		  next if an == 0
		  ecpp = (secondary[ap.index] ? ecp2 : ecp)
		  e = ecp_ary[an] || (ecpp && ecpp[an])
		  if e
		    #  Cache the PNAME of the $ECP entry and re-use it
		    ecp_ary[an] ||= (e.split)[0] + "\n"
		  else
		    e = ap.element.upcase + "-ECP NONE\n"
		  end
		  fp.print e
		}
	    fp.print " $END\n"
	  end
	}
	fname
  end
  
  def cmd_create_gamess_input
    if natoms == 0
      raise MolbyError, "cannot create GAMESS input; the molecule is empty"
    end

	#  Descriptive text and internal string for popup menus
#    bset_desc = ["PM3", "STO-3G", "3-21G", "6-31G", "6-31G(d)", "6-31G(d,p)", "6-311G", "6-311G(d,p)", "LanL2DZ"]
#	bset_internal = ["PM3", "STO3G", "321G", "631G", "631Gd", "631Gdp", "6311G", "6311Gdp", "LanL2DZ"]
    bset_desc = $gamess_basis_keys.map { |key| $gamess_basis_desc[key] }
	dft_desc = ["B3LYP"]
	dft_internal = ["B3LYP"]

	defaults = {"scftype"=>0, "runtype"=>0, "charge"=>"0", "mult"=>"1",
	  "basis"=>4, "use_secondary_basis"=>0, "secondary_elements"=>"",
	  "secondary_basis"=>8, "esp"=>0}

    hash = Dialog.run("GAMESS Export") {
      def load_basis_set_sub(item)
	    fname = Dialog.open_panel("Select a file containing GAMESS basis set:")
		if fname
		  Molecule.read_gamess_basis_sets(fname)
		  bset_desc_new = $gamess_basis_keys.map { |key| $gamess_basis_desc[key] }
		  sel1 = attr("basis", :value)
		  sel2 = attr("secondary_basis", :value)
		  set_attr("basis", :subitems=>bset_desc_new)
		  set_attr("basis", :value=>sel1)
		  set_attr("secondary_basis", :subitems=>bset_desc_new)
		  set_attr("secondary_basis", :value=>sel2)
		end
      end
	  layout(4,
		item(:text, :title=>"SCF type"),
		item(:popup, :subitems=>["RHF", "ROHF", "UHF"], :tag=>"scftype"),
	    item(:text, :title=>"Run type"),
		item(:popup, :subitems=>["Energy", "Property", "Optimize"], :tag=>"runtype",
		  :action=>proc { |it| set_attr("use_internal", :enabled=>(it[:value] == 2)) } ),

		item(:checkbox, :title=>"Use internal coordinates for structure optimization", :tag=>"use_internal"),
		-1, -1, -1,

		item(:text, :title=>"Charge"),
		item(:textfield, :width=>80, :tag=>"charge"),
		item(:text, :title=>"Multiplicity"),
		item(:textfield, :width=>80, :tag=>"mult"),

		item(:checkbox, :title=>"Use DFT", :tag=>"dft",
		  :action=>proc { |it| set_attr("dfttype", :enabled=>(it[:value] != 0)) } ),
		-1,
		item(:text, :title=>"DFT type"),
		item(:popup, :subitems=>dft_desc, :tag=>"dfttype"),
		
		item(:line),
		-1, -1, -1,

		item(:text, :title=>"Basis set"),
		item(:popup, :subitems=>bset_desc, :tag=>"basis"),
		-1,
		-1,

		item(:button, :title=>"Load Basis Set...", :action=>:load_basis_set_sub),
		-1, -1, -1,
		
		item(:checkbox, :title=>"Use secondary basis set", :tag=>"use_secondary_basis",
		  :action=>proc { |it|
		    flag = (it[:value] != 0)
			set_attr("secondary_elements", :enabled=>flag)
			set_attr("secondary_basis", :enabled=>flag)
		  }),
		-1, -1, -1,

		item(:text, :title=>"   Elements"),
		item(:textfield, :width=>80, :tag=>"secondary_elements"),
		item(:text, :title=>"Basis set"),
		item(:popup, :subitems=>bset_desc, :tag=>"secondary_basis"),
		
		item(:line),
		-1, -1, -1,
		
		item(:checkbox, :title=>"Calculate electrostatic potential (ESP)", :tag=>"esp"),
		-1, -1, -1,
	
		item(:line),
		-1, -1, -1
	  )
	  values = Hash.new
	  each_item { |it|
	    tag = it[:tag]
		if tag
		  values[tag] = (get_global_settings("gamess.#{tag}") || defaults[tag])
		  it[:value] = values[tag]
		end
	  }
	  set_attr("secondary_elements", :enabled=>(values["use_secondary_basis"] == 1))
	  set_attr("secondary_basis", :enabled=>(values["use_secondary_basis"] == 1))
	  set_attr("dfttype", :enabled=>(values["dft"] == 1))
	  set_attr("use_internal", :enabled=>(values["runtype"] == 2))
	}
	hash.each_pair { |key, value|
	  next if key == :status
	  set_global_settings("gamess.#{key}", value)
	}
	if hash[:status] == 0
	  #  Specify basis by internal keys
	  hash["basis"] = $gamess_basis_keys[hash["basis"]]
	  hash["secondary_basis"] = $gamess_basis_keys[hash["secondary_basis"]]
	  hash["dfttype"] = dft_internal[hash["dfttype"]]
	  basename = (self.path ? File.basename(self.path, ".*") : self.name)
      fname = Dialog.save_panel("GAMESS input file name", self.dir, basename + ".inp", "GAMESS input file (*.inp)|*.inp|All files|*.*")
	  return nil if !fname
	  export_gamess(fname, hash)
	else
	  nil
	end
  end
  
end

$gamess_basis = {
  "PM3"   => " PM3 0\n",
  "STO3G" => " STO 3\n",
  "321G"  => " N21 3\n",
  "631G"  => " N31 6\n" }
$gamess_basis_desc = {
  "PM3"   => "PM3",
  "STO3G" => "STO-3G",
  "321G"  => "3-21G",
  "631G"  => "6-31G" }
$gamess_basis_keys = ["PM3", "STO3G", "321G", "631G"]

["631Gd", "631Gdp", "631+Gd", "631++Gdp", "6311Gdp", "6311+Gd", "6311++Gdp", "6311++G2d2p", "6311++G3df3pd", "LanL2DZ"].each { |n|
  Molecule.read_gamess_basis_sets("basis_sets/#{n}.txt")
}
