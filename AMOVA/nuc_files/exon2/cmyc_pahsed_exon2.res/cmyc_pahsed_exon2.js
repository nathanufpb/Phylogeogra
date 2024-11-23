
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:///Z:/home/nathan/Documents/programs/arlecore_win/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (cmyc_pahsed_exon2.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 22/11/24 at 23:50:00", "cmyc_pahsed_exon2.htm#22_11_24at23_50_00"))
	insDoc(aux1, gLnk("R", "Settings", "cmyc_pahsed_exon2.htm#22_11_24at23_50_00_run_information"))
		aux2 = insFld(aux1, gFld("Shared haplotypes", "cmyc_pahsed_exon2.htm#22_11_24at23_50_00_shared%20haplotypes"))
		insDoc(aux2, gLnk("R", "PE_N_S", "cmyc_pahsed_exon2.htm#22_11_24at23_50_00_gr_shared0"))
		insDoc(aux2, gLnk("R", "PE_S_N", "cmyc_pahsed_exon2.htm#22_11_24at23_50_00_gr_shared1"))
		aux2 = insFld(aux1, gFld("Samples", ""))
		insDoc(aux2, gLnk("R", "PE_N_S", "cmyc_pahsed_exon2.htm#22_11_24at23_50_00_group0"))
		insDoc(aux2, gLnk("R", "PE_S_N", "cmyc_pahsed_exon2.htm#22_11_24at23_50_00_group1"))
		aux2 = insFld(aux1, gFld("Within-samples summary", ""))
		insDoc(aux2, gLnk("R", "Basic indices", "cmyc_pahsed_exon2.htm#22_11_24at23_50_00_comp_sum_Basic"))
		insDoc(aux2, gLnk("R", "Heterozygosity", "cmyc_pahsed_exon2.htm#22_11_24at23_50_00_comp_sum_het"))
		insDoc(aux2, gLnk("R", "No. of alleles", "cmyc_pahsed_exon2.htm#22_11_24at23_50_00_comp_sum_numAll"))
		insDoc(aux2, gLnk("R", "Molecular diversity", "cmyc_pahsed_exon2.htm#22_11_24at23_50_00_comp_sum_moldiv"))
