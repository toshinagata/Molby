#!/usr/bin/perl

foreach $f (@ARGV) {
  open(FIN, $f) || die "Cannot open $f: $!";
  ($ff = $f) =~ s/\.in$//;
  open(FOUT, ">$ff") || die "Cannot create $ff: $!";
  while ($s = <FIN>) {
    $s =~ s/\{([^}]*)\}/unpack("V", $1."\0\0\0\0\0\0\0\0")/eg;
    print FOUT $s;
  }
}
