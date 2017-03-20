#!/usr/bin/awk -f

# -v K =  number of fields in the key
# -v TH = min. number og register per bed region
# -v FIRST =  first field with aggreable info

BEGIN{
	OFS="\t"
	out=""
}
{
	# set key
	key=$1
	for(i=2; i <= K; i ++){
		key = key "\t" $i
	}

	# if new bed region, print and clear
	if(key != oldkey && NR > 1) {
		# set up output string
		if(n[FIRST] > TH && m[FIRST] != "NA"){
			out = m[FIRST]/n[FIRST]
			p += 1
		}else{
			out = "NA"
		}
		for(j = FIRST + 1; j <= (NF - 1); j++) {
			if(n[j] > TH && m[j] != "NA") {
				out = out "\t" m[j]/n[j]
				p += 1
			}else{
				out = out "\t" "NA"
			}
		}
		# print result
		if(p > 0){
			print oldkey, out
		}
		# clear arrays
		split("", n)
		split("", m)
		p = 0
	}

	# fill arrays
	for(j = FIRST; j <= (NF - 1); j++) {
		if($j != "NA" && $j != "." && $j != -1){
			m[j] += $j; n[j] += 1
		}
	}
	
	# store key
	oldkey=key
}
