# modify .Rev files

mkdir exp_alpha_est
cd exp_alpha_est

rev=unconstrained_tree_inference.Rev 
nchar=300 

for i in `seq 1 50`
do
	file=tree_inference_n${nchar}_${i}.Rev
	cp ../templates/$rev $file
	sed -i '' "s/data\\/bears.nex/..\\/data\\/n$nchar\\/rep_$i.nex/g" $file
	sed -i '' "s/mk/rep_n${nchar}_${i}/g" $file
done

nchar=30

for i in `seq 1 50`
do
	file=tree_inference_n${nchar}_${i}.Rev
	cp ../templates/$rev $file
	sed -i '' "s/data\\/bears.nex/..\\/data\\/n$nchar\\/rep_$i.nex/g" $file
	sed -i '' "s/mk/rep_n${nchar}_${i}/g" $file
done

cd ../


## fixed exponential

mkdir exp_alpha_fixed
cd exp_alpha_fixed

rev=unconstrained_tree_inference_fixed_exp.Rev

nchar=300

for i in `seq 1 50`
do
	file=tree_inference_n${nchar}_${i}.Rev
	cp ../templates/$rev $file
	sed -i '' "s/data\\/bears.nex/..\\/data\\/n$nchar\\/rep_$i.nex/g" $file
	sed -i '' "s/mk/rep_n${nchar}_${i}/g" $file
done

nchar=30

for i in `seq 1 50`
do
	file=tree_inference_n${nchar}_${i}.Rev
	cp ../templates/$rev $file
	sed -i '' "s/data\\/bears.nex/..\\/data\\/n$nchar\\/rep_$i.nex/g" $file
	sed -i '' "s/mk/rep_n${nchar}_${i}/g" $file
done

cd ../


## dirichlet

mkdir dirichlet
cd dirichlet

rev=unconstrained_tree_inference_dir.Rev

nchar=300

for i in `seq 1 50`
do
	file=tree_inference_n${nchar}_${i}.Rev
	cp ../templates/$rev $file
	sed -i '' "s/data\\/bears.nex/..\\/data\\/n$nchar\\/rep_$i.nex/g" $file
	sed -i '' "s/mk/rep_n${nchar}_${i}/g" $file
done

nchar=30

for i in `seq 1 50`
do
	file=tree_inference_n${nchar}_${i}.Rev
	cp ../templates/$rev $file
	sed -i '' "s/data\\/bears.nex/..\\/data\\/n$nchar\\/rep_$i.nex/g" $file
	sed -i '' "s/mk/rep_n${nchar}_${i}/g" $file
done

cd ../


