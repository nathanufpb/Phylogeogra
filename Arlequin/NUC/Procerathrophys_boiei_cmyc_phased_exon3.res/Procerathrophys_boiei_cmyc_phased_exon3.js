
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:///Z:/home/nathan/Documents/programs/arlequin31/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (Procerathrophys_boiei_cmyc_phased_exon3.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 27/11/24 at 15:30:56", "Procerathrophys_boiei_cmyc_phased_exon3.htm#27_11_24at15_30_56"))
	insDoc(aux1, gLnk("R", "Settings", "Procerathrophys_boiei_cmyc_phased_exon3.htm#27_11_24at15_30_56_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure", "Procerathrophys_boiei_cmyc_phased_exon3.htm#27_11_24at15_30_56_gen_struct"))
		insDoc(aux2, gLnk("R", "Euclidean distances", "Procerathrophys_boiei_cmyc_phased_exon3.htm#27_11_24at15_30_56_amova_dm"))
		insDoc(aux2, gLnk("R", "AMOVA", "Procerathrophys_boiei_cmyc_phased_exon3.htm#27_11_24at15_30_56_amova"))
	aux1 = insFld(foldersTree, gFld("Run of 27/11/24 at 15:38:23", "Procerathrophys_boiei_cmyc_phased_exon3.htm#27_11_24at15_38_23"))
	insDoc(aux1, gLnk("R", "Settings", "Procerathrophys_boiei_cmyc_phased_exon3.htm#27_11_24at15_38_23_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure", "Procerathrophys_boiei_cmyc_phased_exon3.htm#27_11_24at15_38_23_gen_struct"))
		insDoc(aux2, gLnk("R", "Euclidean distances", "Procerathrophys_boiei_cmyc_phased_exon3.htm#27_11_24at15_38_23_amova_dm"))
		insDoc(aux2, gLnk("R", "AMOVA", "Procerathrophys_boiei_cmyc_phased_exon3.htm#27_11_24at15_38_23_amova"))
	aux1 = insFld(foldersTree, gFld("Run of 27/11/24 at 21:40:01", "Procerathrophys_boiei_cmyc_phased_exon3.htm#27_11_24at21_40_01"))
	insDoc(aux1, gLnk("R", "Settings", "Procerathrophys_boiei_cmyc_phased_exon3.htm#27_11_24at21_40_01_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure", "Procerathrophys_boiei_cmyc_phased_exon3.htm#27_11_24at21_40_01_gen_struct"))
		insDoc(aux2, gLnk("R", "Euclidean distances", "Procerathrophys_boiei_cmyc_phased_exon3.htm#27_11_24at21_40_01_amova_dm"))
		insDoc(aux2, gLnk("R", "AMOVA", "Procerathrophys_boiei_cmyc_phased_exon3.htm#27_11_24at21_40_01_amova"))
