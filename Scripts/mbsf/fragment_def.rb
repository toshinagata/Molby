#
#  This script is loaded from formula.rb
#

class Molecule

  #  "fragment name", "mbsf name" | "other fragment name", "dummy atoms"
  table = [
	"Alicyclic",
	["cyclopropane", "alicyclic/cyclopropane.mbsf"],
	["C3H6", "cyclopropane"],
	["cyclopropyl", "cyclopropane", "H11"],
	["C3H5", "cyclopropyl"],

	["cyclobutane", "alicyclic/cyclobutane.mbsf"],
	["C4H8", "cyclobutane"],
	["cyclobutyl", "cyclobutane", "H11"],
	["C4H7", "cyclobutyl"],
	
	["cyclopentane", "alicyclic/cyclopentane.mbsf"],
	["C5H10", "cyclopentane"],
	["cyclopentyl", "cyclopentane", "H11"],
	["C5H9", "cyclopentyl"],

	["cyclohexane", "alicyclic/cyclohexane.mbsf"],
	["C6H12", "cyclohexane"],
	["cyclohexyl", "cyclohexane", "H11"],
	["C6H11", "cyclohexyl"],
	["Cy", "cyclohexyl"],
	
	["cycloheptane", "alicyclic/cycloheptane.mbsf"],
	["C7H14", "cycloheptane"],
	["cycloheptyl", "cycloheptane", "H11"],
	["C7H13", "cycloheptyl"],
	
	["cyclooctane", "alicyclic/cyclooctane.mbsf"],
	["C8H16", "cyclooctane"],
	["cyclooctyl", "cyclooctane", "H11"],
	["C8H15", "cyclooctyl"],

	["cyclohexane (twist boat)", "alicyclic/cyclohexane-twist-boat.mbsf"],

    "Aromatic",
    ["benzene", "aromatic/benzene.mbsf"],
    ["C6H6", "benzene"],
    ["C6H5", "benzene", "H1"],
    ["Ph", "C6H5"],
	["phenyl", "C6H5"],
    ["C6H4", "benzene", "H1", "H4"],
    ["C6H3", "benzene", "H1", "H3", "H5"],
	
	["cyclopentadienyl", "aromatic/cyclopentadienyl.mbsf"],
	["C5H5", "cyclopentadienyl"],
	["Cp", "cyclopentadienyl"],
	["C5H4", "cyclopentadienyl", "H1"],
	
	["cycloheptatrienyl", "aromatic/cycloheptatrienyl.mbsf"],
	["C7H7", "cycloheptatrienyl"],
	["C7H6", "cycloheptatrienyl", "H1"],
	
    ["C60 fullerene", "aromatic/c60.mbsf"],
    
	"Heterocyclic",
	["furan", "heterocyclic/furan.mbsf"],
	["imidazole", "heterocyclic/imidazole.mbsf"],
	["oxazole", "heterocyclic/oxazole.mbsf"],
	["phthalocyanine", "heterocyclic/phthalocyanine.mbsf"],
	["porphine", "heterocyclic/porphine.mbsf"],
	["pyrazine", "heterocyclic/pyrazine.mbsf"],
	["pyrazole", "heterocyclic/pyrazole.mbsf"],
	["pyridazine", "heterocyclic/pyridazine.mbsf"],
	["pyridine", "heterocyclic/pyridine.mbsf"],
	["pyrimidine", "heterocyclic/pyrimidine.mbsf"],
	["pyrrole", "heterocyclic/pyrrole.mbsf"],
	["thiazole", "heterocyclic/thiazole.mbsf"],
	["thiophene", "heterocyclic/thiophene.mbsf"],
	
	"Coordination",
	["MX2 linear", "coordination/MX2.mbsf"],
	["MX3 trigonal", "coordination/MX3y.mbsf"],
	["MX3 T-shape", "coordination/MX3t.mbsf"],
	["MX4 tetrahedral", "coordination/MX4t.mbsf"],
	["MX4 square-planar", "coordination/MX4p.mbsf"],
	["MX5 trigonal bipyramidal", "coordination/MX5.mbsf"],
	["MX6 octahedral", "coordination/MX6.mbsf"]
  ]

  $named_fragments = []
  
  table.each { |t|
	if t.is_a?(String)
	  $named_fragments.push([t, "-"])  #  Subtitle
	  next
	end
    fname = t[1]
    if fname =~ /\.mbsf$/
      f = new(MbsfPath + "/" + fname)  #  Molecule.new
	  $named_fragments.push([t[0], t[1]]) #  Registered fragment
    else
      f = known_fragment(fname)
    end
    (2...t.length).each { |i|
      p = f.atoms[t[i]]
      if p
        p.name = "_#{i - 1}"
        p.atom_type = ""
        p.element = "Du"
      end
    }
	f.dummies = f.find_dummy_atoms
    register_fragment(t[0], f)
  }

  $named_fragments.concat([            #  Non-registered fragments
    ["Solvent boxes", "-"],
	["CHCl3", "solvents/chcl3box.mbsf"],
	["DMSO", "solvents/dmsobox.mbsf"],
	["MeOH", "solvents/meohbox.mbsf"],
	["N-methylacetamide", "solvents/nmabox.mbsf"],
	["H2O (tip3p)", "solvents/tip3pbox.mbsf"]])

end
