#BSUB -P runSJARACNe
#BSUB -n 1
#BSUB -M 8000
#BSUB -oo runSJARACNe.out -eo runSJARACNe.err
#BSUB -J runSJARACNe
#BSUB -q rhel8_standard

#sjaracne lsf -e /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/Monocyte/Monocyte.8573_1837.exp.txt -g /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/Monocyte/TF/Monocyte.832_1837.tf.txt -o /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/Monocyte/TF/bt100_pc001 -n 100 -pc 1e-2 -j /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/Monocyte/config_cwlexec.json
sjaracne lsf -e /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/Monocyte/Monocyte.8573_1837.exp.txt -g /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/Monocyte/SIG/Monocyte.4191_1837.sig.txt -o /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/Monocyte/SIG/bt100_pc001 -n 100 -pc 1e-2 -j /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/Monocyte/config_cwlexec.json
