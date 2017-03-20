#!/usr/bin/awk -f

# Usage: smoothed.real.dmr.metrics.awk -v n=<number_fields> <region.file + methylation.file>
# Output: refion.file +, sum_of_differences, max_difference, n_CpG, average_difference

# methylation.file: chr_1, start_1, end_1, strand_1, u_1, m_1, s_1, ls_1, us_1, chr_2, start_2, end_2, strand_2, u_2, m_2, s_2, ls_2, us_2 

# define absolute value
function abs(x)
{
	return x<0 ? -x : x
}

# define max
function max(a, b)
{
	if(a > b){
		return a
	}else{
		return b
	}
}


# set initial values
BEGIN{
	FS = "\t"
	OFS = "\t"
	dsum = 0
	dmax = 0
	dn = 0
}
{
	# get key (first n fields)
	key=$1
	for(i=2; i<= n; i++) key = key "\t" $i

	# if the key has changed, print and clear
	if( key != oldkey && NR > 1){
		if( dsum == 0){
			dmax = 0
			dmean = 0
		}else{
			dmax = dmax*dsum/abs(dsum)
			dmean = dsum/dn
		}
		# print output
		print oldkey, dsum, dmax, dn, dmean
		# clear variables
		dsum = 0
		dmax = 0
		dn = 0
	}
		
	# get smoothed methylation
	l1 = $(n + 7)
	l2 = $(n + 16)

	# compute metrics
	if( l1 != "." && l2 != "."  && l1 != "" && l2 != ""){
		d = l1 - l2
		dsum += d
		dmax = max(abs(d), abs(dmax))
		dn += 1
	}
	
	# store key
	oldkey = key
}
END{
	if( dsum == 0){
		dmax = 0
		dmean = 0
	}else{
		dmax = dmax*dsum/abs(dsum)
		dmean = dsum/dn
	}
	# print output
	print oldkey, dsum, dmax, dn, dmean
}
