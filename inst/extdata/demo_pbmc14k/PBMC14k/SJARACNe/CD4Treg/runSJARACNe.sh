#BSUB -P runSJARACNe
#BSUB -n 1
#BSUB -M 8000
#BSUB -oo runSJARACNe.out -eo runSJARACNe.err
#BSUB -J runSJARACNe
#BSUB -q rhel8_standard

#sjaracne lsf -e /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD4Treg/CD4Treg.8659_1985.exp.txt -g /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD4Treg/TF/CD4Treg.838_1985.tf.txt -o /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD4Treg/TF/bt100_pc001 -n 100 -pc 1e-2 -j /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD4Treg/config_cwlexec.json
sjaracne lsf -e /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD4Treg/CD4Treg.8659_1985.exp.txt -g /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD4Treg/SIG/CD4Treg.4214_1985.sig.txt -o /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD4Treg/SIG/bt100_pc001 -n 100 -pc 1e-2 -j /research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD4Treg/config_cwlexec.json
