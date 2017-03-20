#!/usr/bin/awk -f

# Usage: hmt.summary.awk -v n=<number_fields> <region.file + methylation.file>
# Output: region.file + nCpG, u, m, sumSmo, sumL

# methylation.file: chr_1, start_1, end_1, strand, u, m , smo, ll, ul

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
	u = 0
	m = 0
	smo = 0
	lsum = 0
	ncpg = 0
}
{
	# get key (first n fields)
	key=$1
	for(i=2; i<= n; i++) key = key "\t" $i

	# if the key has changed, print and clear
	if( key != oldkey && NR > 1){
		# print output
		print oldkey, ncpg, u, m, smo, lsum
		# clear variables
		u = 0
		m = 0
		smo = 0
		lsum = 0
		ncpg = 0
	}
		
	# get methylation counts
	l = 0
	t = $(n + 5) + $(n + 6)
	if( t > 0){ l = $(n + 6)/t
	u += $(n + 5)
	m += $(n + 6)
	smo += $(n + 7)
	lsum += l
	ncpg += 1
	}
	# store key
	oldkey = key
}
END{
	print oldkey, ncpg, u, m, smo, lsum
}
