perl -pe 's/\r\n|\n|\r/\n/g' /home/cqian/PBMCdemo/Sjaracne//CD8em_821_821_89/tf/CD8em_87_87_89_tf.txt > /home/cqian/PBMCdemo/Sjaracne//CD8em_821_821_89/tf/CD8em_87_87_89_tf.txt.tmp
rm /home/cqian/PBMCdemo/Sjaracne//CD8em_821_821_89/tf/CD8em_87_87_89_tf.txt
mv /home/cqian/PBMCdemo/Sjaracne//CD8em_821_821_89/tf/CD8em_87_87_89_tf.txt.tmp /home/cqian/PBMCdemo/Sjaracne//CD8em_821_821_89/tf/CD8em_87_87_89_tf.txt
perl -pe 's/\r\n|\n|\r/\n/g' /home/cqian/PBMCdemo/Sjaracne//CD8em_821_821_89/CD8em_821_821_89.exp > /home/cqian/PBMCdemo/Sjaracne//CD8em_821_821_89/CD8em_821_821_89.exp.tmp
rm /home/cqian/PBMCdemo/Sjaracne//CD8em_821_821_89/CD8em_821_821_89.exp
mv /home/cqian/PBMCdemo/Sjaracne//CD8em_821_821_89/CD8em_821_821_89.exp.tmp /home/cqian/PBMCdemo/Sjaracne//CD8em_821_821_89/CD8em_821_821_89.exp
