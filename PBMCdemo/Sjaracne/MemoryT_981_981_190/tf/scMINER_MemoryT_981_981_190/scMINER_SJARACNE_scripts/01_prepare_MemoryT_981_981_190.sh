perl -pe 's/\r\n|\n|\r/\n/g' /home/cqian/PBMCdemo/Sjaracne//MemoryT_981_981_190/tf/MemoryT_110_110_190_tf.txt > /home/cqian/PBMCdemo/Sjaracne//MemoryT_981_981_190/tf/MemoryT_110_110_190_tf.txt.tmp
rm /home/cqian/PBMCdemo/Sjaracne//MemoryT_981_981_190/tf/MemoryT_110_110_190_tf.txt
mv /home/cqian/PBMCdemo/Sjaracne//MemoryT_981_981_190/tf/MemoryT_110_110_190_tf.txt.tmp /home/cqian/PBMCdemo/Sjaracne//MemoryT_981_981_190/tf/MemoryT_110_110_190_tf.txt
perl -pe 's/\r\n|\n|\r/\n/g' /home/cqian/PBMCdemo/Sjaracne//MemoryT_981_981_190/MemoryT_981_981_190.exp > /home/cqian/PBMCdemo/Sjaracne//MemoryT_981_981_190/MemoryT_981_981_190.exp.tmp
rm /home/cqian/PBMCdemo/Sjaracne//MemoryT_981_981_190/MemoryT_981_981_190.exp
mv /home/cqian/PBMCdemo/Sjaracne//MemoryT_981_981_190/MemoryT_981_981_190.exp.tmp /home/cqian/PBMCdemo/Sjaracne//MemoryT_981_981_190/MemoryT_981_981_190.exp
