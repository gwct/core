#!/usr/bin/python
########################################################################################
#A general purpose batch FASTA header editing script. Can relabel and trim FASTA headers.
#Originally created as part of my tree making workflow.
#
#Dependencies: core
#
#Gregg Thomas, Summer 2015
########################################################################################

import sys, os, argparse
import core

#BABOON
#"ENSBTA:cow,ENSCJA:marmoset,ENSCAF:dog,ENSECA:horse,ENSP00:human,ENSMMU:macaque,ENSMUS:mouse,ENSNLE:gibbon,ENSPTR:chimp,ENSPAN:baboon,ENSPPY:orang,ENSRNO:rat"

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

############################################
#Function Definitions
############################################
def IO_fileParse():
#This function handles the command line options.

	parser = argparse.ArgumentParser(description="A general purpose FASTA editing script.");

	parser.add_argument("-i", dest="input", help="A directory containing FASTA formatted alignments.");
	parser.add_argument("-r", dest="relabel_opt", help="Boolean to tell the script whether to relabel the FASTA headers (1) or not (0). Default: 1", type=int, default=1);
	parser.add_argument("-s", dest="spec_dict", help="A string formatted as a Python dictionary with the current species ID as the key and the label to add to the beginning of the FASTA header as the value. Must be provided if -r set to 1.");
	parser.add_argument("-t", dest="trim_opt", help="Boolean to tell the script whether to trim the FASTA headers (1) or not (0). Default: 0", type=int, default=0);
	parser.add_argument("-d", dest="trim_delim", help="The character string at which to trim the FASTA headers if -t is set to 1. Default: \" \"", default=" ");
	parser.add_argument("-p", dest="ss_opt", help="Boolean to tell the script whether to remove start and stops from the alignment (1) or not (0). Default: 0", type=int, default=0);
	parser.add_argument("-o", dest="output", help="The directory to which the relabeled and/or trimmed FASTA sequences are written.");

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


	return args.input, args.relabel_opt, sd, args.trim_opt, args.trim_delim, args.ss_opt, args.output;

############################################
#Main Block
############################################

ins, r, specdict, t, td, ss, outs = IO_fileParse();
suffix = "";

print "==============================================================================================";
print "\t\t\tFASTA editing";
print "\t\t\t" + core.getDateTime();
if os.path.isfile(ins):
	print "Editing single file:\t\t" + ins;
	filelist = [ins];
else:
	print "Editing files in directory:\t\t\t\t" + ins;
	filelist = os.listdir(ins);
if r == 0 and t == 0 and ss == 0:
	print "NOT trimming OR relabeling FASTA headers OR removing start/stop positions.\nSimply re-writing the sequences.";
	suffix = ".cp";
else:
	if r == 1:
		print "---";
		print "Relabeling FASTA header(s).";
		print "EXISTING LABEL\t\tADDED LABEL";
		for sid in specdict:
			print sid + "\t\t\t" + specdict[sid];
		print "---";
		suffix = suffix + ".r";
	else:
		print "NOT relabeling FASTA header(s).";
	if t == 1:
		if td == " ":
			print "Trimming FASTA header(s) at first occurrence of:\t\" \"";
		else:
			print "Trimming FASTA header(s) at first occurrence of:\t" + td;
		suffix = suffix + ".t";
	else:
		print "NOT trimming FASTA header(s).";
	if ss == 1:
		print "Removing start and stop positions.";
		suffix = suffix + ".ss";
	else:
		print "NOT removing start and stop positions.";
if os.path.isfile(ins):
	print "Writing output to the following file:\t\t" + outs;
else:
	print "Writing output to the following directory:\t\t" + outs;
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
			for specid in specdict:
				if title.find(specid) != -1:
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
			if inseqs[title][len(inseqs[title])-1] == "*":
				inseqs[title] = inseqs[title][:len(inseqs[title])-1];

	for title in inseqs:
	#Writes the sequences to the output file.
		core.writeSeqOL(outfilename, inseqs[title], title);

if not os.path.isfile(ins):
	pstring = "100.0% complete.";
	sys.stderr.write('\b' * len(pstring) + pstring);
print "\n" + core.getTime() + " Done!";
