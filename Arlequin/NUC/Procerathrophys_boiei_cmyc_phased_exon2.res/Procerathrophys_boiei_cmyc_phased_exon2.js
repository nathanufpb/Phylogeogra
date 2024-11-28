
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:///Z:/home/nathan/Documents/programs/arlequin31/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (Procerathrophys_boiei_cmyc_phased_exon2.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 27/11/24 at 15:19:47", "Procerathrophys_boiei_cmyc_phased_exon2.htm#27_11_24at15_19_47"))
	insDoc(aux1, gLnk("R", "Settings", "Procerathrophys_boiei_cmyc_phased_exon2.htm#27_11_24at15_19_47_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure", "Procerathrophys_boiei_cmyc_phased_exon2.htm#27_11_24at15_19_47_gen_struct"))
		insDoc(aux2, gLnk("R", "Euclidean distances", "Procerathrophys_boiei_cmyc_phased_exon2.htm#27_11_24at15_19_47_amova_dm"))
		insDoc(aux2, gLnk("R", "AMOVA", "Procerathrophys_boiei_cmyc_phased_exon2.htm#27_11_24at15_19_47_amova"))
	aux1 = insFld(foldersTree, gFld("Run of 27/11/24 at 21:41:44", "Procerathrophys_boiei_cmyc_phased_exon2.htm#27_11_24at21_41_44"))
	insDoc(aux1, gLnk("R", "Settings", "Procerathrophys_boiei_cmyc_phased_exon2.htm#27_11_24at21_41_44_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure", "Procerathrophys_boiei_cmyc_phased_exon2.htm#27_11_24at21_41_44_gen_struct"))
		insDoc(aux2, gLnk("R", "Euclidean distances", "Procerathrophys_boiei_cmyc_phased_exon2.htm#27_11_24at21_41_44_amova_dm"))
		insDoc(aux2, gLnk("R", "AMOVA", "Procerathrophys_boiei_cmyc_phased_exon2.htm#27_11_24at21_41_44_amova"))
