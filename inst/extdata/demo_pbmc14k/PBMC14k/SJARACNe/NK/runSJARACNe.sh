#BSUB -P runSJARACNe
#BSUB -n 1
#BSUB -M 8000
#BSUB -oo runSJARACNe.out -eo runSJARACNe.err
#BSUB -J Test
#BSUB -q rhel8_standard

#sjaracne lsf -e /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/NK/NK.8654_1936.exp.txt -g /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/NK/TF/NK.841_1936.tf.txt -o /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/NK/TF/bt100_pc001 -n 100 -pc 1e-2 -j /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/NK/config_cwlexec.json
sjaracne lsf -e /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/NK/NK.8654_1936.exp.txt -g /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/NK/SIG/NK.4202_1936.sig.txt -o /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/NK/SIG/bt100_pc001 -n 100 -pc 1e-2 -j /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/NK/config_cwlexec.json
