#!/usr/bin/perl
#  amberparm2molby: read amber .dat file and convert to a Molby parameter file
#  input (stdin): a .dat file
#  output (stdout): a Molby parm file
#  Written as amberparm2namd.pl by Toshi Nagata, Dec 2 2004
#  Updated for Molby, Nov 13 2009

$title = <>;
print "! ", $title;

$blanks = ' ' x 100;

sub format_com {
	my ($com) = @_;
	$com =~ s/(^\s*)|(\s*$)//g;
	$com =~ s/^\!\s*//;
	if ($com eq "!") {
		$com = "";
	} elsif ($com ne "") {
		$com = " ! " . $com;
	}
	return $com;
}

while (<>) {
	chomp;
	last if /^\s*$/;
	($aname, $molw, $pol, $com) = unpack("A2 x1 A10 A10 A77", $_ . $blanks);
	push @atoms, $aname;
	$molws{$aname} = $molw;
	$comments{$aname} = format_com($com);
}

#  Hydrophilic atoms
$line = <>;

#  Bonds
while (<>) {
	chomp;
	last if /^\s*$/;
	($a1, $a2, $fc, $len, $com) = unpack("A2 x1 A2 A10 A10 A75", $_ . $blanks);
	printf "bond %-2s   %-2s   %7.1f %5.2f%s\n", $a1, $a2, $fc, $len, format_com($com);
}
print "\n";

#  Angles
while (<>) {
	chomp;
	last if /^\s*$/;
	($a1, $a2, $a3, $fc, $ang, $com) = unpack("A2 x1 A2 x1 A2 A10 A10 A72", $_ . $blanks);
	printf "angle %-2s   %-2s   %-2s   %7.1f %5.2f%s\n", $a1, $a2, $a3, $fc, $ang, format_com($com);
}
print "\n";

#  Dihedrals
while (<>) {
	chomp;
	last if /^\s*$/;
	($a1, $a2, $a3, $a4, $div, $fc, $ang, $per, $com) = unpack("A2 x1 A2 x1 A2 x1 A2 A4 A15 A15 A15 A40", $_ . $blanks);
	if ($div > 0) {
		$fc /= $div;
	}
	next if $per < 0;
	printf "dihe %-2s   %-2s   %-2s   %-2s %7.2f %7d %7.2f%s\n", $a1, $a2, $a3, $a4, $fc, $per, $ang, format_com($com);
}
print "\n";

#  Impropers
while (<>) {
	chomp;
	last if /^\s*$/;
	($a1, $a2, $a3, $a4, $div, $fc, $ang, $per, $com) = unpack("A2 x1 A2 x1 A2 x1 A2 A4 A15 A15 A15 A40", $_ . $blanks);
	printf "impr %-2s   %-2s   %-2s   %-2s %7.2f %7d %7.2f%s\n", $a1, $a2, $a3, $a4, $fc, $per, $ang, format_com($com);
}
print "\n";

$line = <>;  #  H-bond param
$line = <>;  #  blank

#  Atom equivalences
while (<>) {
	chomp;
	last if /^\s*$/;
	$label = "";
	foreach $i (unpack("A2 x2" x 20, $_ . $blanks)) {
		if ($label eq "") {
			$label = $i;
		} elsif ($i !~ /^\s*$/) {
			push @{$equil{$label}}, $i;
		}
	}
}

while (<>) {
	chomp;
	last if /^\s*$/;
	($label, $kind) = unpack("A4 x6 A2", $_. $blanks);
	if ($label eq "END") {
		last;
	}
	if ($kind ne "RE") {
		print STDERR "The vdW parameter other than 'RE' format ($kind) is not supported yet.\n";
		print STDERR "\$_='$_'\n\$label='$label'\n\$kind='$kind'\n";
		exit(1);
	}
	while (<>) {
		chomp;
		last if /^\s*$/;
		($sym, $r, $edep) = unpack("x2 A2 x6 A10 A10", $_ . $blanks);
		$sigma = 2 * $r * 0.890899;  #  sigma = 2R * (1/2)**(1/6)
		printf "vdw %-2s %7.4f %7.4f %7.4f %7.4f %d %7.3f%s\n", $sym, $edep, $r, $edep, $r, 0, $molws{$sym}, $comments{$sym};
		foreach $sym2 (@{$equil{$sym}}) {
			printf "vdw %-2s %7.4f %7.4f %7.4f %7.4f %d %7.3f%s\n", $sym2, $edep, $r, $edep, $r, 0, $molws{$sym2}, $comments{$sym2};
		}
	}
}
print "\n";
