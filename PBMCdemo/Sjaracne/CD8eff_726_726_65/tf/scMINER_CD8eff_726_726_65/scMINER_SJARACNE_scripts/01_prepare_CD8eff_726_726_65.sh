perl -pe 's/\r\n|\n|\r/\n/g' /home/cqian/PBMCdemo/Sjaracne//CD8eff_726_726_65/tf/CD8eff_71_71_65_tf.txt > /home/cqian/PBMCdemo/Sjaracne//CD8eff_726_726_65/tf/CD8eff_71_71_65_tf.txt.tmp
rm /home/cqian/PBMCdemo/Sjaracne//CD8eff_726_726_65/tf/CD8eff_71_71_65_tf.txt
mv /home/cqian/PBMCdemo/Sjaracne//CD8eff_726_726_65/tf/CD8eff_71_71_65_tf.txt.tmp /home/cqian/PBMCdemo/Sjaracne//CD8eff_726_726_65/tf/CD8eff_71_71_65_tf.txt
perl -pe 's/\r\n|\n|\r/\n/g' /home/cqian/PBMCdemo/Sjaracne//CD8eff_726_726_65/CD8eff_726_726_65.exp > /home/cqian/PBMCdemo/Sjaracne//CD8eff_726_726_65/CD8eff_726_726_65.exp.tmp
rm /home/cqian/PBMCdemo/Sjaracne//CD8eff_726_726_65/CD8eff_726_726_65.exp
mv /home/cqian/PBMCdemo/Sjaracne//CD8eff_726_726_65/CD8eff_726_726_65.exp.tmp /home/cqian/PBMCdemo/Sjaracne//CD8eff_726_726_65/CD8eff_726_726_65.exp
