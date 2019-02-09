perl -pe 's/\r\n|\n|\r/\n/g' /home/cqian/PBMCdemo/Sjaracne//NaiveT_884_884_156/tf/NaiveT_92_92_156_tf.txt > /home/cqian/PBMCdemo/Sjaracne//NaiveT_884_884_156/tf/NaiveT_92_92_156_tf.txt.tmp
rm /home/cqian/PBMCdemo/Sjaracne//NaiveT_884_884_156/tf/NaiveT_92_92_156_tf.txt
mv /home/cqian/PBMCdemo/Sjaracne//NaiveT_884_884_156/tf/NaiveT_92_92_156_tf.txt.tmp /home/cqian/PBMCdemo/Sjaracne//NaiveT_884_884_156/tf/NaiveT_92_92_156_tf.txt
perl -pe 's/\r\n|\n|\r/\n/g' /home/cqian/PBMCdemo/Sjaracne//NaiveT_884_884_156/NaiveT_884_884_156.exp > /home/cqian/PBMCdemo/Sjaracne//NaiveT_884_884_156/NaiveT_884_884_156.exp.tmp
rm /home/cqian/PBMCdemo/Sjaracne//NaiveT_884_884_156/NaiveT_884_884_156.exp
mv /home/cqian/PBMCdemo/Sjaracne//NaiveT_884_884_156/NaiveT_884_884_156.exp.tmp /home/cqian/PBMCdemo/Sjaracne//NaiveT_884_884_156/NaiveT_884_884_156.exp
