# modify .Rev files 
for i in `seq 1 50` 
do 
	cp templates/unconstrained_tree_inference.Rev tree_inference_$i.Rev
	sed -i -e s/bears.nex/rep_$i.nex/g tree_inference_$i.Rev
	sed -i -e s/mk/rep_$i/g tree_inference_$i.Rev
done

