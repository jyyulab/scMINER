#BSUB -P runSJARACNe
#BSUB -n 1
#BSUB -M 8000
#BSUB -oo runSJARACNe.out -eo runSJARACNe.err
#BSUB -J runSJARACNe
#BSUB -q rhel8_standard

#sjaracne lsf -e /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD8TN/CD8TN.8568_1984.exp.txt -g /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD8TN/TF/CD8TN.825_1984.tf.txt -o /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD8TN/TF/bt100_pc001 -n 100 -pc 1e-2 -j /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD8TN/config_cwlexec.json
sjaracne lsf -e /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD8TN/CD8TN.8568_1984.exp.txt -g /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD8TN/SIG/CD8TN.4145_1984.sig.txt -o /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD8TN/SIG/bt100_pc001 -n 100 -pc 1e-2 -j /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD8TN/config_cwlexec.json
