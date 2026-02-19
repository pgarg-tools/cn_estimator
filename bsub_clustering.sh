file=$1
prefix=`basename $file | sed -e 's/.rds//'`
prefix="./test/logs/"${prefix}
bsub="bsub\
  -P acc_sharpa01a\
  -q premium\
  -J "$prefix"\
  -n 48 -R span[hosts=1] \
  -R "rusage[mem=21000]" -R himem \
  -W 7:59 \
  -oo $prefix.stdout \
  -eo $prefix.stderr"
echo $bsub "\"ml R/4.2.0; Rscript scripts/cnv_clustering.r $file ./test/\""
