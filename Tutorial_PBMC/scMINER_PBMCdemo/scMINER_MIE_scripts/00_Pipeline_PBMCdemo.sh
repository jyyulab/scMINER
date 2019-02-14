psub -K -P PBMCdemo -J PBMCdemo_MIE_Calc -q compbio -M 2000 -i /home/cqian/PBMCdemo/scMINER_PBMCdemo/scMINER_MIE_scripts/03_Calc_PBMCdemo.sh -oo /home/cqian/PBMCdemo/scMINER_PBMCdemo/scMINER_MIE_log/PBMCdemo_MIE_Calc.%J.%I.out -eo /home/cqian/PBMCdemo/scMINER_PBMCdemo/scMINER_MIE_log/PBMCdemo_MIE_Calc.%J.%I.err 
sleep 30
psub -K -P PBMCdemo -J PBMCdemo_MIE_Merge -q compbio -M 2000 -i /home/cqian/PBMCdemo/scMINER_PBMCdemo/scMINER_MIE_scripts/04_Merge_PBMCdemo.sh -oo /home/cqian/PBMCdemo/scMINER_PBMCdemo/scMINER_MIE_log/PBMCdemo_MIE_Merge.%J.%I.out -eo /home/cqian/PBMCdemo/scMINER_PBMCdemo/scMINER_MIE_log/PBMCdemo_MIE_Merge.%J.%I.err 
sleep 30
psub -K -P PBMCdemo -J PBMCdemo_MIE_Norm -q compbio -M 2000 -i /home/cqian/PBMCdemo/scMINER_PBMCdemo/scMINER_MIE_scripts/05_Norm_PBMCdemo.sh -oo /home/cqian/PBMCdemo/scMINER_PBMCdemo/scMINER_MIE_log/PBMCdemo_MIE_Norm.%J.%I.out -eo /home/cqian/PBMCdemo/scMINER_PBMCdemo/scMINER_MIE_log/PBMCdemo_MIE_Norm.%J.%I.err 
