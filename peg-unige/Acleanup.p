#!/usr/local/bin/perl -w
# Acleanup.p  reads an alist and writes an identical alist matrix
#              in which the mlist is correct
# usage
# Acleanup.p  inA > outA
#
# 99 12 29
#

$degreesprovided = 1 ; # a normal Alist does include the degrees of rows and columns, though we ignore them here.
# Since they are ignored, some files we receive might not have them!
$rowsnotcols = 0 ; # my convention is to supply the columns, and clean up based on them.
# but the IBM files supply the rows.

# usage in /home/mackay/code/IBMPEG/codes:
# Acleanup.p rowsnotcols=1 degreesprovided=0 Reg504x1008.dat > PEGReg504x1008
# Acleanup.p rowsnotcols=1 degreesprovided=0 Reg252x504.dat > PEGReg252x504
# Acleanup.p rowsnotcols=1 degreesprovided=0 irReg252x504.dat > PEGirReg252x504
# Acleanup.p rowsnotcols=1 degreesprovided=0 irReg504x1008.dat > PEGirReg504x1008
# Acleanup.p rowsnotcols=1 degreesprovided=0 irUppTriang1030x2048.dat > PEGirUppTriang1030x2048
# Acleanup.p rowsnotcols=1 degreesprovided=0 irUppTriang518x1024.dat  > PEGirUppTriang518x1024

eval "\$$1=\$2"  while @ARGV && $ARGV[0]=~ /^([\w\[\]]+)=(.*)/ && shift ;

print STDERR "--------Acleanup.p----------- DJCM 99/12/29\n" ; 
# first, read the old alist in

if($degreesprovided) {
    $_ = <>;
    ($L,$M) = split ;
    $_ = <>; # max numbers
    $_ = <>;
    @oldnumnlist = split ; 
    $_ = <>;
# @oldnummlist = split ; 
} else { # IBM files
    $_ = <>;
    ($L) = split ;
    $_ = <>;
    ($M) = split ;
    $_ = <>; # max number
}    
print STDERR "------------\nRead in $L  $M\n" ;
$N = $L ; 
for ( $n = 1 ; $n <= $N ; $n ++ ) {
    $nlist[$n] = "" ;	# 
    $numnlist[$n] = 0 ;
}
for ( $m = 1 ; $m <= $M ; $m ++ ) {
    $mlist[$m] = "" ; 
    $nummlist[$m] = 0 ;
}

if ($rowsnotcols) {
    for ( $l = 1 ; $l <= $M ; $l ++ ) {
	$_ = <>;
	@y = split ;
	for ( $u = 0 ; $u <= $#y ; $u ++ ){
	    $ty = $y[$u] ;
	    if ( $ty > 0 ) {
		$n = $ty ; $m = $l ; 
		&mnadd($m,$n);
	    }
	}
    }                               #
} else { # normal mode
    for ( $l = 1 ; $l <= $L ; $l ++ ) {
	$_ = <>;
	@y = split ;
	for ( $u = 0 ; $u <= $#y ; $u ++ ){
	    $ty = $y[$u] ;
	    if ( $ty > 0 ) {
		$m = $ty ; $n = $l ; 
		&mnadd($m,$n);
	    }
	}
    }                               #
}

print STDERR "done N = $N ; M = $M\n" ; 

if ( 1 ) {
    $head = "$N $M\n" ; 
    print STDERR "Writing code.\n" ;

# find $t (max) and $tr (max)

    $t = 0 ; $tr = 0 ; 
    for ( $n = 1 ; $n <= $N ; $n ++ ) {
	if ( $numnlist[$n] > $t ) {
	    $t = $numnlist[$n] ;
	}
    }
    for ( $m = 1 ; $m <= $M ; $m ++ ) {
	if ( $nummlist[$m] > $tr ) {
	    $tr = $nummlist[$m] ;
	}
    }
    $head .= "$t $tr\n" ;
# append zeros to short lines
    for ( $n = 1 ; $n <= $N ; $n ++ ) {
	$thisn = $numnlist[$n] ;
	if ( $thisn ) {
	    $head .= "$thisn " ; 
	    while ( $thisn < $t ) {
		$nlist[$n] .= "0\t" ;
		$thisn ++ ;
	    }
	} # otherwise it is a dead column
    }
    $head .= "\n" ;
    
    for ( $m = 1 ; $m <= $M ; $m ++ ) {
	$thism = $nummlist[$m] ; 
	$head .= "$thism " ; 
	while ( $thism < $tr ) {
	    $mlist[$m] .= "0\t" ;
	    $thism ++ ;
	}
    }
    $head .= "\n" ;
    
    print $head ; 
#    print STDERR $head ; 
# finish lists 
    for ( $n = 1 ; $n <= $N ; $n ++ ) {
	$thisn = $numnlist[$n] ;
	if ( $thisn ) {
	    $nlist[$n] =~ s/\t$/\n/ ; 
	    print $nlist[$n] ; 
	}
    }
    for ( $m = 1 ; $m <= $M ; $m ++ ) {
	$mlist[$m] =~ s/\t$/\n/ ; 
	print $mlist[$m] ; 
    }
}


sub mnadd {
    local($m,$n)=@_;

    if ($m>0) {				# 
	# put m and n on each others lists 
	$nlist[$n] .= $m."\t" ;
	$mlist[$m] .= $n."\t" ;
	$numnlist[$n] ++ ;
	$nummlist[$m] ++ ;
    }
}

# remove from the mlist only
sub mnsubtract {
    local($m,$n)=@_;

    $count = ( $mlist[$m] =~ s/\b$n\t// ) ;
    if ( $count ) {
	$nummlist[$m] -- ;
    } else {
	print STDERR "error, unable to remove $n from mlist $m\n $mlist[$m]\n" ; 
    }
}






