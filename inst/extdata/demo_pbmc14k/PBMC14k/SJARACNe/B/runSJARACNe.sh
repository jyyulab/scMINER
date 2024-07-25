#BSUB -P runSJARACNe
#BSUB -n 1
#BSUB -M 8000
#BSUB -oo runSJARACNe.out -eo runSJARACNe.err
#BSUB -J runSJARACNe
#BSUB -q rhel8_standard

#sjaracne lsf -e /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/B/B.8572_1902.exp.txt -g /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/B/TF/B.835_1902.tf.txt -o /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/B/TF/bt100_pc001 -n 100 -pc 1e-2 -j /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/B/config_cwlexec.json
sjaracne lsf -e /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/B/B.8572_1902.exp.txt -g /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/B/SIG/B.4148_1902.sig.txt -o /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/B/SIG/bt100_pc001 -n 100 -pc 1e-2 -j /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/B/config_cwlexec.json
