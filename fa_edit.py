#!/usr/bin/python
########################################################################################
#A general purpose batch FASTA header editing script. Can relabel and trim FASTA headers.
#Originally created as part of my tree making workflow.
#
#Dependencies: core
#
#Gregg Thomas, Summer 2015
########################################################################################

import sys, os, random, argparse
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/corelib/"))
import core

#BABOON
#"ENSBTA:cow,ENSCJA:marmoset,ENSCAF:dog,ENSECA:horse,ENSP00:human,ENSMMU:macaque,ENSMUS:mouse,ENSNLE:gibbon,ENSPTR:chimp,ENSPAN:baboon,ENSPPY:orang,ENSRNO:rat,ENSMOD:possum"

#ARANEAE
#"LHESP:Latrodectus_hesperus WesternBlackWidow,LRECL:Loxosceles_reclusa BrownRecluseSpider,PTEPI:Parasteatoda_tepidariorum CommonHouseSpider,407821:Stegodyphus_mimosarum AfricanSocialVelvetSpider"

#COLEOPTERA
#"217634:Anoplophora_glabripennis AsianLonghornedBeetle,APLAN:Agrilus_planipennis EmeraldAshBorer,77166:Dendroctonus_ponderosae MountainPineBeetle,7539:Leptinotarsa_decemlineata ColoradoPotatoBeetle,OTAUR:Onthophagus_taurus BullHeadedDungBeetle,7070:Tribolium_castaneum RedFlourBeetle"

#DIPTERA
#"7159:Aedes_egypti YellowfeverMosquito,7167:Anopheles_albimanus NewWorldMalariaMosquito,62324:Anopheles_funestus AfricanMalariaMosquito,7165:Anopheles_gambiae AfricanMalariaMosquito,7213:Ceratitis_capitata MediterraneanFruitFly,7176:Culex_quinquefasciatus SouthernHouseMosquito,7222:Drosophila_grimshawi FruitFly,7227:Drosophila_melanogaster FruitFly,7237:Drosophila_pseudoobscura FruitFly,7394:Glossina_morsitans SavannahTsetseFly,LCUP2:Lucilia_cuprina SheepBlowfly,7200:Lutzomyia_longipalpis SandFly,39758:Mayetiola_destructor HessianFly,7370:Musca_domestica HouseFly"

#HEMIPTERA
#"HVITR:Homalodisca_vitripennis GlassyWingedSharpshooter,HHALY:Halyomorpha_halys BrownMarmoratedStinkBug,OFAS2:Oncopeltus_fasciatus MilkweedBug,PVENU:Pachypsylla_venusta HackberryPetioleGallPsyllid,GBUEN:Gerris_buenoi WaterStrider,7029:Acyrthosiphon_pisu PeaAphid,79782:Cimex_lectularius BedBug"

#HYMENOPTERA
#"12957:Atta_cephalotes LeafcutterAnt,103372:Acromyrmex_echinatior PanamanianLeafcutterAnt,7463:Apis_florea LittleHoneyBee,7460:Apis_mellifera HoneyBee,37344:Athalia_rosae TurnipSawfly,132113:Bombus_impatiens CommonEasternBumbleBee,30195:Bombus_terrestris Buff-tailedBumbleBee,211228:Cephus_cinctus WheatStemSawfly,104421:Camponotus_floridanus FloridaCarpenterAnt,286306:Cardiocondyla_obscurior TrampAnt,COPFL:Copidosoma_floridanum Wasp,178035:Dufourea_novaeangliae SolitaryUnivoltineBee,516756:Eufriesea_mexicana OrchidBee,597456:Habropoda_laboriosa SoutheasternBlueberryBee,610380:Harpegnathos_saltator JerdonsJumpingAnt,88501:Lasioglossum_albipes HalictidBee,83485:Linepithema_humile ArgentineAnt,166423:Melipona_quadrifasciata NeotropicalStinglessBee,143995:Megachile_rotundata AlfalfaLeafcuttingBee,7425:Nasonia_vitripennis ParasiticJewelWasp,222816:Orussus_abietinus ParasiticWoodWasp,144034:Pogonomyrmex_barbatus RedHarvesterAnt,13686:Solenopsis_invicta RedFireAnt,7493:Trichogramma_pretiosum ParasiticStinglessWasp"

#LEPIDOPTERA
#"7091:Bombyx_mori DomesticSilkwormMoth,13037:Danaus_plexippus MonarchButterfly,34740:Heliconius_melpomene PostmanButterfly,7130:Manduca_sexta TobaccoHornwormMoth,51655:Plutella_xylostella DiamondbackMoth"

#ARTHROPODA
#"LHESP:Latrodectus_hesperus WesternBlackWidow,LRECL:Loxosceles_reclusa BrownRecluseSpider,PTEPI:Parasteatoda_tepidariorum CommonHouseSpider,407821:Stegodyphus_mimosarum AfricanSocialVelvetSpider,217634:Anoplophora_glabripennis AsianLonghornedBeetle,APLAN:Agrilus_planipennis EmeraldAshBorer,77166:Dendroctonus_ponderosae MountainPineBeetle,7539:Leptinotarsa_decemlineata ColoradoPotatoBeetle,OTAUR:Onthophagus_taurus BullHeadedDungBeetle,7070:Tribolium_castaneum RedFlourBeetle,7159:Aedes_egypti YellowfeverMosquito,7167:Anopheles_albimanus NewWorldMalariaMosquito,62324:Anopheles_funestus AfricanMalariaMosquito,7165:Anopheles_gambiae AfricanMalariaMosquito,7213:Ceratitis_capitata MediterraneanFruitFly,7176:Culex_quinquefasciatus SouthernHouseMosquito,7222:Drosophila_grimshawi FruitFly,7227:Drosophila_melanogaster FruitFly,7237:Drosophila_pseudoobscura FruitFly,7394:Glossina_morsitans SavannahTsetseFly,LCUP2:Lucilia_cuprina SheepBlowfly,7200:Lutzomyia_longipalpis SandFly,39758:Mayetiola_destructor HessianFly,7370:Musca_domestica HouseFly,HVITR:Homalodisca_vitripennis GlassyWingedSharpshooter,HHALY:Halyomorpha_halys BrownMarmoratedStinkBug,OFAS2:Oncopeltus_fasciatus MilkweedBug,PVENU:Pachypsylla_venusta HackberryPetioleGallPsyllid,GBUEN:Gerris_buenoi WaterStrider,7029:Acyrthosiphon_pisu PeaAphid,79782:Cimex_lectularius BedBug,12957:Atta_cephalotes LeafcutterAnt,103372:Acromyrmex_echinatior PanamanianLeafcutterAnt,7463:Apis_florea LittleHoneyBee,7460:Apis_mellifera HoneyBee,37344:Athalia_rosae TurnipSawfly,132113:Bombus_impatiens CommonEasternBumbleBee,30195:Bombus_terrestris Buff-tailedBumbleBee,211228:Cephus_cinctus WheatStemSawfly,104421:Camponotus_floridanus FloridaCarpenterAnt,286306:Cardiocondyla_obscurior TrampAnt,COPFL:Copidosoma_floridanum Wasp,178035:Dufourea_novaeangliae SolitaryUnivoltineBee,516756:Eufriesea_mexicana OrchidBee,597456:Habropoda_laboriosa SoutheasternBlueberryBee,610380:Harpegnathos_saltator JerdonsJumpingAnt,88501:Lasioglossum_albipes HalictidBee,83485:Linepithema_humile ArgentineAnt,166423:Melipona_quadrifasciata NeotropicalStinglessBee,143995:Megachile_rotundata AlfalfaLeafcuttingBee,7425:Nasonia_vitripennis ParasiticJewelWasp,222816:Orussus_abietinus ParasiticWoodWasp,144034:Pogonomyrmex_barbatus RedHarvesterAnt,13686:Solenopsis_invicta RedFireAnt,7493:Trichogramma_pretiosum ParasiticStinglessWasp,7091:Bombyx_mori DomesticSilkwormMoth,13037:Danaus_plexippus MonarchButterfly,34740:Heliconius_melpomene PostmanButterfly,7130:Manduca_sexta TobaccoHornwormMoth,51655:Plutella_xylostella DiamondbackMoth,BGERM:Blatella_germanica GermanCockroachm,CAQUI:Catajapyx_aquilonaris SilvestrisNorthernForcepstail,CSCUL:Centruroides_sculpturatus BarkScorpion,6669:Daphnia_pulex WaterFlea,EAFFI:Eurytemora_affinis CalanoidCopepod,1049336:Ephemera_danica GreenDrakeMayfly,133901:Frankliniella_occidentalis WesternFlowerThrip,6945:Ixodes_scapularis BlackLeggedDeerTick,123851:Ladona_fulva ScarceChaserDragonfly,LLUNA:Limnephilus_lunatus Caddisfly,34638:Metaseiulus_occidentalis WesternPredatoryMite,121225:Pediculus_humanus BodyLouse,126957:Strigamia_maritima GeophilomorphCentipede,32264:Tetranychus_urticae TwoSpottedSpiderMite,136037:Zootermopsis_nevadensis DampwoodTermite,HAZTE:Hyalella_azteca FreshwaterEpibenthicAmphipod"

aas = ["G","A","V","L","I","P","F","Y","W","S","T","C","M","N","Q","K","R","H","D","E", "X"];
nts = ["A","T","C","G"];

############################################
#Function Definitions
############################################
def IO_fileParse():
#This function handles the command line options.

	parser = argparse.ArgumentParser(description="A general purpose FASTA editing script.");

	parser.add_argument("-i", dest="input", help="A directory containing FASTA formatted files or a single FASTA file.");
	parser.add_argument("-r", dest="relabel_opt", help="Boolean to tell the script whether to relabel the FASTA headers (1) or not (0). Default: 1", type=int, default=1);
	parser.add_argument("-s", dest="spec_dict", help="A string formatted as a Python dictionary with the current species ID as the key and the label to add to the beginning of the FASTA header as the value. Must be provided if -r set to 1.");
	parser.add_argument("-t", dest="trim_opt", help="Boolean to tell the script whether to trim the FASTA headers (1) or not (0). Default: 0", type=int, default=0);
	parser.add_argument("-d", dest="trim_delim", help="The character string at which to trim the FASTA headers if -t is set to 1. Default: \" \"", default=" ");
	parser.add_argument("-p", dest="ss_opt", help="Boolean to tell the script whether to remove start and stops from the alignment (1) or not (0). Default: 0", type=int, default=0);
	parser.add_argument("-m", dest="replacement", help="This option will replace all characters in each sequence with another character. For example, AB will replace all As with Bs. If the input is an alignment, if A: is entered, all As will be replaced with another AA that is not present in the column.", default="");
	parser.add_argument("-o", dest="output", help="The directory or file to which the relabeled and/or trimmed FASTA sequences are written.");

	args = parser.parse_args();

	if args.input == None or args.output == None:
		parser.print_help();
		sys.exit();

	if args.relabel_opt not in [0,1]:
		print " ------------------------------------------------";
		print "|**Error 1: -r must take values of either 0 or 1 |";
		print " ------------------------------------------------";
		parser.print_help();
		sys.exit();
	elif args.relabel_opt == 1:
		if args.spec_dict == None:
			print " -------------------------------------------------------";
			print "|**Error 2: With -r set to 1, -s must also be specified |";
			print " -------------------------------------------------------";
			parser.print_help();
			sys.exit();
		else:
			specs = args.spec_dict.split(",");
			sd = {};
			for each in specs:
				spec = each.split(":");
				sd[spec[0]] = spec[1];
	else:
		sd = "";

	if args.trim_opt not in [0,1]:
		print " ------------------------------------------------";
		print "|**Error 3: -t must take values of either 0 or 1 |";
		print " ------------------------------------------------";
		parser.print_help();
		sys.exit();

	if args.ss_opt not in [0,1]:
		print " ------------------------------------------------";
		print "|**Error 4: -p must take values of either 0 or 1 |";
		print " ------------------------------------------------";
		parser.print_help();
		sys.exit();

	replacement = args.replacement;
	if args.replacement != "":
		replacement = replacement.upper();

		if len(replacement) > 2 or (replacement[1] not in aas and replacement[1] != ":"):
			print " -----------------------------------------------------------------------------------------";
			print "|**Error 5: For -m, the second character entered must be a valid amino acid symbol or ':' |";
			print " -----------------------------------------------------------------------------------------";
			parser.print_help();
			sys.exit();


	return args.input, args.relabel_opt, sd, args.trim_opt, args.trim_delim, args.ss_opt, replacement, args.output;

############################################
#Main Block
############################################

ins, r, specdict, t, td, ss, repl, outs = IO_fileParse();
suffix = "";

print "==============================================================================================";
print "\t\t\tFASTA editing";
print "\t\t\t" + core.getDateTime();
if os.path.isfile(ins):
	print "INPUT    | Editing single file:\t\t" + ins;
	filelist = [ins];
else:
	print "INPUT    | Editing files in directory:\t\t\t\t" + ins;
	filelist = os.listdir(ins);
if r == 0 and t == 0 and ss == 0 and repl == "":
	print "INFO     | NOT trimming OR relabeling FASTA headers OR removing start/stop positions.\nSimply re-writing the sequences.";
	suffix = ".cp";
else:
	if r == 1:
		print "---";
		print "INFO     | Relabeling FASTA header(s).";
		print "EXISTING LABEL\t\tADDED LABEL";
		for sid in specdict:
			print sid + "\t\t\t" + specdict[sid];
		print "---";
		suffix = suffix + ".l";
	else:
		print "INFO     | NOT relabeling FASTA header(s).";
	if t == 1:
		if td == " ":
			print "INFO     | Trimming FASTA header(s) at first occurrence of:\t\" \"";
		else:
			print "INFO     | Trimming FASTA header(s) at first occurrence of:\t" + td;
		suffix = suffix + ".t";
	else:
		print "INFO     | NOT trimming FASTA header(s).";
	if ss == 1:
		print "INFO     | Removing start and stop positions.";
		suffix = suffix + ".ss";
	else:
		print "INFO     | NOT removing start and stop positions.";
	if repl != "" and repl[1] != ":":
		print "INFO     | Replacing all " + repl[0] + " symbols in input sequences with " + repl[1] + ".";
		suffix = suffix + ".r";
	elif repl != "" and repl[1] == ":":
		print "INFO     | Replacing all " + repl[0] + " symbols in input sequences with a random AA symbol not present in that column.";
		suffix = suffix + ".r";
if os.path.isfile(ins):
	print "OUTPUT   | Writing output to the following file:\t\t" + outs;
else:
	if outs[len(outs)-1] != "/":
		outs = outs + "/";
	print "OUTPUT   | Writing output to the following directory:\t\t" + outs;
	if not os.path.exists(outs):
		print "+Creating output directory.";
		os.system("mkdir " + outs);
print "---------------------------------------------";
#sys.exit();

print core.getTime() + " Starting...";

numfiles = len(filelist);
numbars = 0;
donepercent = [];
i = 0;

for each in filelist:
	if each.find(".fa") == -1:
		continue;
	
	if os.path.isfile(ins):
		print ins;
		infilename = ins;
		outfilename = outs;

	else:
		numbars, donepercent = core.loadingBar(i, numfiles, donepercent, numbars);
		i = i + 1;
		infilename = ins + each;
		outfilename = outs + each[:each.index(".fa")] + suffix + ".fa";

	outfile = open(outfilename, "w");
	outfile.write("");
	outfile.close();

	inseqs = core.fastaGetDict(infilename);
	#Reads the sequences from the input file.

	if r == 1:
	#Re-labels the FASTA headers.
		rseqs = {};
		for title in inseqs:
			curid = title[1:title.index(":")];
			for specid in specdict:
				if curid.find(specid) != -1:
					newtitle = ">" + specdict[specid] + " " + title[1:];
					rseqs[newtitle] = inseqs[title];
					break;
		inseqs = rseqs;

	if t == 1:
	#Trims the FASTA headers.
		tseqs = {};
		for title in inseqs:
			newtitle = title[:title.index(td)];
			tseqs[newtitle] = inseqs[title];
		inseqs = tseqs;

	if ss == 1:
	#Removes start and stop codons.
		for title in inseqs:
			if inseqs[title][0] == "M":
				inseqs[title] = inseqs[title][1:];
			#if inseqs[title][len(inseqs[title])-1] == "*":
			#	inseqs[title] = inseqs[title][:len(inseqs[title])-1];
			if inseqs[title].find("*") != -1:
				inseqs[title] = inseqs[title].replace("*","");

	if repl != "":
	#Replaces symbols in the sequences.
		r1 = repl[0];
		r2 = repl[1];
		seqlen = len(inseqs[inseqs.keys()[0]]);
		if r2 == ":":
			pseqs = {};
			for x in xrange(0,seqlen):
				present = [];
				for title in inseqs:
					if x == 0:
						pseqs[title] = "";
					if inseqs[title][x] not in present:
						present.append(inseqs[title][x]);
				if r1 in present:
					newaa = present[0];
					while newaa in present:
						newaa = random.choice(aas);
					for title in inseqs:
						if inseqs[title][x] == r1:
							pseqs[title] = pseqs[title] + newaa;
						else:
							pseqs[title] = pseqs[title] + inseqs[title][x];
				else:
					for title in inseqs:
						pseqs[title] = pseqs[title] + inseqs[title][x];
			inseqs = pseqs;
		else:
			for title in inseqs:
				inseqs[title] = inseqs[title].replace(r1,r2);		

	for title in inseqs:
	#Writes the sequences to the output file.
		core.writeSeqOL(outfilename, inseqs[title], title);

if not os.path.isfile(ins):
	pstring = "100.0% complete.";
	sys.stderr.write('\b' * len(pstring) + pstring);
print "\n" + core.getTime() + " Done!";
print "==============================================================================================";
