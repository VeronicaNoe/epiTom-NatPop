metabolite=( $(cat metaTarget) )
marker=( $(cat markerTarget) )
kinship=( $(cat kinshipTarget) )

for meta in "${metabolite[@]}"; do
	for mm in "${marker[@]}"; do
		for mmKin in "${marker[@]}"; do
			for kin in "${kinship[@]}"; do
    				echo "${meta}_${mm}_${mmKin}_${kin}" >> toTarget
	  		done
		done
	done
done
