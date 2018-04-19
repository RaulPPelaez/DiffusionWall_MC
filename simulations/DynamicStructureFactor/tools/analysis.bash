#!/bin/bash
#This script will analyse results output by mc. It will average and fit the dynamic structure factors to exponentials and store the relaxation times for each wave number. Also it wil lcompute rdf (if possible) and height probability distribution.

#It expects to find a data.main, hydroGridOptions.nml, uammd.log and the results of the simulation in the .. folder.
#It needs xmgrace to run

dataFolder=..
#Check requirements
if ! ls $dataFolder/data.main >/dev/null 2>&1
then
    echo "Could not find data.main! Be sure to put the simulation in the $dataFolder folder"
exit 1
fi
if ! ls $dataFolder/hydro*.nml >/dev/null 2>&1
then
    echo "Could not find hydroGridOptions.nml! Be sure to put the simulation in the $dataFolder folder"
exit 1
fi
if ! ls $dataFolder/uammd.log >/dev/null 2>&1
then
    echo "Could not find uammd.log! Be sure to put the simulation in the $dataFolder folder"
exit 1
fi

#Check datamash version
useDatamash=false
if which datamash >/dev/null 2>&1
then
    useDatamash=true
    version=$(datamash --version | head -1 | awk '{print $4}' | sed 's+\.++g')
    if [ $version -lt 110 ]
    then
	echo "I need datamash version >= 1.1.0!"
	echo "I wont use datamash"
	useDatamash=false
    fi

fi

if ! which xmgrace >/dev/null 2>&1
then
    echo "I need xmgrace to fit curves!, please install it"
    exit 1
fi

#Averages a stream of lines, grouping them by the value of the first column.
#Takes two arguments, the first and last column to average
#This function is equivalent to datamash -W --sort groupby 1 mean $1-$2 sstdev $1-$2
function average_files(){
    if [ $# -eq 1 ]
    then
	startCol=$1
	endCol=$startCol
    else
	startCol=$1
	endCol=$2
    fi
    if $useDatamash
    then
	if [ $# -eq 1 ]
	then
	    cols=$(echo "$startCol")
	else
	    cols=$(echo "$startCol-$endCol")
	fi
	cat /dev/stdin | datamash -W --sort groupby 1 mean $cols sstdev $cols
    else
	
	command=$(echo 'BEGIN{
        $ini='$startCol'-1;
        $end='$endCol'-1;
        }
        {
        $counter{$F[0]}++;
        $n=$counter{$F[0]};
        for $i ($ini..$end){
           $delta= $F[$i] - $mean{$F[0]}[$i];
           $mean{$F[0]}[$i] += $delta/$n;
           $delta2= $F[$i] - $mean{$F[0]}[$i];
           $M2{$F[0]}[$i] += $delta*$delta2
        }
        }
        END{
        foreach $val (keys %mean){
           printf "$val ";
           for $j ($ini..$end){
             printf "%e ",$mean{$val}[$j];
           }
           for $j ($ini..$end){
             if($counter{$val}<=1){
               printf "nan ";
             }
             else{
               printf "%e ",sqrt($M2{$val}[$j]/($counter{$val}-1));
             }
           }

           printf"\n"}
        }')

	command=$(echo $command | tr '\n' ' ');

	cat /dev/stdin | perl -lane  "$command"
    fi
}


if ! which xmgrace >/dev/null 2>&1
then
   echo "ERROR! I need xmgrace to fit functions!"
   exit  1
fi


outputname=$(grep -E '^\s*outputname' $dataFolder/data.main | awk '{print $3}')
Lx=$(cat $dataFolder/data.main | grep "lx"  | awk '{print $3}')
Ly=$(cat $dataFolder/data.main | grep "ly"  | awk '{print $3}')

i2kx=$(echo $Lx | awk '{print 2*3.141592/$1}')
i2ky=$(echo $Ly | awk '{print 2*3.141592/$1}')

#Compute wavenumber modulus for each column in the files
cat $dataFolder/uammd.log |
    grep "Tracking" |
    awk '{print $4,$6,$7}' |
    awk 'BEGIN{a=1}{if($1<a)exit; a=$1; print}' |
    awk '{print $1+1, sqrt(($2*'$i2kx')^2+($3*'$i2ky')^2)}' |
    sort -g -k2 > column2Kmod.dat


ncolumns=$(cat column2Kmod.dat | wc -l | awk '{print $1+1}')

meanFile=$outputname.S_k_t.pair=1.mean.dat
if(! ls $meanFile >/dev/null 2>&1)
then
    echo "Averaging S(k,w) time windows"
    #Average all S(k,w) files (they should be identical)
    ls $dataFolder/$outputname.*.S_k_t.pair=1.dat --hide "*.00000000.*" |
	xargs cat |
	#datamash -W --sort groupby 1 mean 2-$ncolumns sstdev 2-$ncolumns |
	average_files 2 $ncolumns |
	sort -g -k1 > $meanFile
fi    

tmpDir=$(mktemp -d)

#mkdir -p tmp

echo "Splitting S(k,w) file"
awk '{for(i=2;i<='$ncolumns';i++) print $1, $i >> "'$tmpDir'/" sprintf("%d", i)}' $meanFile

echo "Fitting S(k,w)"
tools=$(pwd)
cd $tmpDir
#Fit each Sk
for i in $(ls)
do

    kmod=$(awk '$1=='$i'{print $2}' $tools/column2Kmod.dat)           
    #cat $i >> $kmod

    initialA=$(head -1 $i | awk '{print $2}')
    Gamma=0.001

    correlation=0
    minCorrelation=0.9
    maxIter=20
    iter=0
    endT=$(tail -1 $i | awk '{print $1}')
    
    maxTime=$(awk 'BEGIN{t=0} $2<=('$initialA'/100){t=$1*4; print t; exit} END{ if(t==0){ t='$endT'; print t}}' $i)

    
    tail -n+2 $i | awk '$1<='$maxTime'' > fitfile
    size=$(cat fitfile | wc -l)

    while (( $(echo "$correlation < $minCorrelation" | bc -l) ))
    do
	let iter++

	if [ $iter -gt 15 ]
	then
	    minCorrelation=0.89
	    tail -n+2 $i | awk '$1<='$maxTime'' | awk 'NR<='$size'/2' > fitfile
	fi

	if [ $iter -gt $maxIter ]
	then
	    echo "I could not fit k=$kmod!, current fit correlation: $correlation"
	    break;
	fi

	cat $tools/fitExp.bat |
	    sed -r 's/(a0\s+=).*/\1 '$initialA'/g' |
	    sed -r 's/(a1\s+=).*/\1 '$Gamma'/g' |
	    gracebat fitfile -bat -  | tee fit.log |
	    grep -E '^\s+a.*' |	    
	    tail -2 > fit.tmp

	correlation=$(cat fit.log | grep Correlation | tail -1 | cut -d":" -f2)

	initialA=$(cat fit.tmp | grep "a0" | cut -d"=" -f2 | awk '{print $1}')
	Gamma=$(cat fit.tmp | grep "a1" | cut -d"=" -f2 | awk '{print $1}')
	
#	echo $correlation $kmod $initialA $Gamma
#	exit

	rm -f fit.log fit.tmp *ps
    done
    rm -f fitfile

    echo $kmod $(echo $Gamma | awk '{print 1.0/$1}') >> tauvsk.dat

    rm $i
done
cd -

#Average tau results for each kmod
cat $tmpDir/tauvsk.dat |
    #    datamash -W groupby 1 mean 2 sstdev 2 > tauvsk.dat
    average_files 2 | sort -g -k1 > tauvsk.dat

rm -rf $tmpDir


#Compute rdf
if which rdf >/dev/null 2>&1 && ls $dataFolder/$outputname.particle.pos >/dev/null 2>&1
then
    echo "Compute rdf"
    N=$(tail -n+2 $dataFolder/$outputname.particle.pos | awk 'substr($1,0,1)=="#"{print NR-1; exit}')
    Nsnapshots=$(grep -c '#' $dataFolder/$outputname.particle.pos)
    lx=$(grep '^lx' $dataFolder/data.main | awk '{print $3}')
    ly=$(grep '^ly' $dataFolder/data.main | awk '{print $3}')
    lz=$(grep '^lz' $dataFolder/data.main | awk '{print $3}')
    rcut=$(echo $lx | awk '{print $1/10}')
    cat $dataFolder/$outputname.particle.pos |
	rdf -N $N -Nsnapshots $Nsnapshots -Lx $lx -Ly $ly -Lz $lz  -nbins 100 -rcut $rcut -dim q2D > $outputname.rdf.dat
fi

echo "Compute height probability distribution"
#Compute height probability distribution
if ls $dataFolder/$outputname.particle.pos >/dev/null 2>&1
then
    lx=$(grep '^lz' $dataFolder/data.main | awk '{print $3}')
    ini=$(echo $lx | awk '{print -$1/2}')
    end=$(echo $lx | awk '{print $1/2}')

    nbins=100
    grep -v "#"  $dataFolder/$outputname.particle.pos |	
	awk '{print $3}' |
	awk '{print int(($1-('$ini'))/('$end'-('$ini'))*'$nbins')}' |
	sort -g -k1 |
	uniq -c |
	awk '{print ($2/'$nbins'.0)*(('$end')-('$ini'))+('$ini'), $1}' |
	sort -g -k1 > kk
    
    sum=$(cat kk | datamash -W sum 2)

    awk '{print $1, $2/'$sum'}' kk > $outputname.heightProbabilityDistribution.dat;
    rm kk
fi

echo "Storing results"

resultFolder=$outputname.results
mkdir -p $resultFolder

mv $outputname.heightProbabilityDistribution.dat $resultFolder
mv $outputname.rdf.dat $resultFolder
mv tauvsk.dat $resultFolder/$outputname.tauvsk.dat


