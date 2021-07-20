
#!/bin/bash
N=1000

stride=100
#50

nsteps=10000000

#t1s="10.00 20.00 40.00 100.00 200.00 400.00 1000.00"
#"400.00 1000.00 2000.00 4000.00 10000.00"
t1s="148.03"

t2s="100.49"
#t2s="5.00 10.00 15.00 20.00 25.00 30.00 35.00 40.00 45.00 50.00"

as="-0.00294"
#as="0.000 0.001 0.005 0.010"

#bs="0.00 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10"
bs="-0.09886"

md=$(pwd)

for a in $as

do

for b in $bs

do

for t1 in $t1s

do

for t2 in $t2s

do

folder=$md"/results_mem_ohs/results_a=${a}_b=${b}_t1=${t1}_t2=${t2}/"

mkdir ${folder}

cd ${folder}
eval "cat ../../spectra_mem_template.sh | sed -e s/XXXNXXX/"$N"/g | sed -e s/XXXstrideXXX/"$stride"/g | sed -e s/XXXnstepsXXX/"$nsteps"/g | sed -e s/XXXt1XXX/"$t1"/g | sed -e s/XXXt2XXX/"$t2"/g | sed -e s/XXXaXXX/"$a"/g | sed -e s/XXXbXXX/"$b"/g > submit.sh"
sbatch submit.sh
#sh submit.s
sleep 0.1

done

done

done

done


