#!/usr/bin/python
########################################################################################
#Part of my Ensembl tree making pipeline. This takes a list of orthologous sequence IDs
#(ie the output from orth_combine.py) and combines the actual sequences into individual
#fasta files for future alignment.
#
#Gregg Thomas, Spring 2015
########################################################################################

#"ENSTTR:Tursiops_truncatus.turTru1.74.cds.all.fa_isofiltered,ENSMLU:Myotis_lucifugus.Myoluc2.0.74.cds.all.fa_isofiltered,ENSMUS:Mus_musculus.GRCm38.74.cds.all.fa_isofiltered,ENSPVA:Pteropus_vampyrus.pteVam1.74.cds.all.fa_isofiltered,ENSBTA:Bos_taurus.UMD3.1.74.cds.all.fa_isofiltered,ENST00:Homo_sapiens.GRCh37.74.cds.all.fa_isofiltered,ENSLAF:Loxodonta_africana.loxAfr3.74.cds.all.fa_isofiltered,ENSCJA:Callithrix_jacchus.C_jacchus3.2.1.74.cds.all.fa_isofiltered,ENSVPA:Vicugna_pacos.vicPac1.74.cds.all.fa_isofiltered,ENSMOD:Monodelphis_domestica.BROADO5.75.cds.all.fa_isofiltered"
#"ENSPFO:molly_pfor_filtered.fa,ENSAMX:cavefish_amex_filtered.fa,ENSDAR:zebrafish_drer_filtered.fa,ENSTNI:puffer_tnir_filtered.fa,ENSGAC:stickleback_gacu_filtered.fa,ENSONI:tilapia_onil_filtered.fa,ENSORL:medaka_olat_filtered.fa,ENSXMA:platyfish_xmac_filtered.fa"
#"ENSPFO:Poecilia_formosa.PoeFor_5.1.2.pep.all.fa,ENSAMX:Astyanax_mexicanus.AstMex102.pep.all.fa,ENSDAR:Danio_rerio.Zv9.pep.all.fa,ENSTNI:Tetraodon_nigroviridis.TETRAODON8.pep.all.fa,ENSGAC:Gasterosteus_aculeatus.BROADS1.pep.all.fa,ENSONI:Oreochromis_niloticus.Orenil1.0.pep.all.fa,ENSORL:Oryzias_latipes.MEDAKA1.pep.all.fa,ENSXMA:Xiphophorus_maculatus.Xipmac4.4.2.pep.all.fa"
#"ENSBTA:cow_peptides_filtered.fa,ENSCJA:marmoset_peptides_filtered.fa,ENSCAF:dog_peptides_filtered.fa,ENSECA:horse_peptides_filtered.fa,ENSP00:human_peptides_filtered.fa,ENSMMU:macaque_peptides_filtered.fa,ENSMUS:mouse_peptides_filtered.fa,ENSNLE:gibbon_peptides_filtered.fa,ENSPTR:chimp_peptides_filtered.fa,ENSPAN:baboon_peptides_filtered.fa,ENSPPY:orang_peptides_filtered.fa,ENSRNO:rat_peptides_filtered.fa,ENSMOD:opossum_peptides_filtered.fa"
#"ENSBTA:cow_peptides_filtered.fa,ENSCJA:marmoset_peptides_filtered.fa,ENSP00:human_peptides_filtered.fa,ENSMMU:macaque_peptides_filtered.fa,ENSMUS:mouse_peptides_filtered.fa,ENSPTR:chimp_peptides_filtered.fa,ENSPAN:baboon_peptides_filtered.fa,ENSPPY:orang_peptides_filtered.fa"
#"chimp:chimp-cds.fa,cow:cow-cds.fa,gibbon:gibbon-cds.fa,gorilla:gorilla-cds.fa,human:human-cds.fa,macaque:macaque-cds.fa,mlemur:mlemur-cds.fa,mouse:mouse-cds.fa,orang:orang-cds.fa,owlmonkey:owlmonkey-cds.fa,vervet:vervet-cds.fa"

import sys, argparse, os
import core

############################################
#Function Definitions
############################################

def optParse():
#This function handles the command line options.

	parser = argparse.ArgumentParser();

	parser.add_argument("-i", dest="input_file", help="A file with tab delimited lists of orthologs (The output from orth_combine.py).");
	parser.add_argument("-s", dest="seq_dir", help="A directory with full fasta sequences from all the species found in the ortholog list.");
	parser.add_argument("-d", dest="spec_dict", help="A necessary option to associate filename with species identifier... format exactly as follows for ALL species: \"spec1ID:spec1.fa,spec2ID:spec2.fa\"");
	parser.add_argument("-delim", dest="title_delim", help="A character to trim the fasta titles at. To trim at the first space, enter 'space'. Default: take whole title. ", default=False)
	parser.add_argument("-m", dest="rem_start", help="A boolean option to either remove the start Methionine (1) or not (0). Default: 0", type=int, default=0);
	parser.add_argument("-o", dest="output_dir", help="Output directory where all combined sequences will be written.");

	args = parser.parse_args();

	if args.input_file == None or args.output_dir == None or args.seq_dir == None or args.spec_dict == None:
		sys.exit(core.errorOut(1, "-i, -o, -s, and -d must all be defined"));

	if args.rem_start not in [0,1]:
		sys.exit(core.errorOut(2, "-m must take values of either 0 or 1"));

	if args.title_delim == 'space':
		args.title_delim = ' ';

	return args.input_file, args.seq_dir, args.spec_dict.split(","), args.title_delim, args.rem_start, args.output_dir;

############################################
#Main Block
############################################

infilename, seqdir, speclist, delim, remstart, outdir = optParse();

#idmap = { 'chimp':'ENSPTR', 'cow':'ENSBTA', 'gibbon':'ENSNLE', 'gorilla':'ENSGGO', 'human':'ENSG00', 'macaque':'ENSMMU', 'mlemur':'ENSMIC', 'mouse':'ENSMUS', 'orang':'ENSPPY', 'owlmonkey':'ENSANA', 'vervet':'ENSCSA' }
#idmap = { 'ENSPTR':'chimp', 'ENSBTA':'cow', 'ENSNLE':'gibbon', 'ENSGGO':'gorilla', 'ENSG00':'human', 'ENSMMU':'macaque', 'ENSMIC':'mlemur', 'ENSMUS':'mouse', 'ENSPPY':'orang', 'ENSANA':'owlmonkey', 'ENSCSA':'vervet' }

print("# =======================================================================");
print("# \t\t\tRetrieving sequences in FASTA format");
print("# \t\t\t" + core.getDateTime());
print("# Retrieving ortholog sequences from:\t", infilename);
print("# Sequence directory\t\t\t" + seqdir);
if delim:
	print("# Trimming FASTA headers at:\t", delim);
else:
	print("# Retaining full FASTA headers.")
print("# Writing combined sequence files to:\t", outdir);
if remstart == 1:
	print("# Removing start Methionines (-m 1)");
else:
	print("# NOT removing start Methionines (-m 0)");
print("# Reminder: Please ensure your species dictionary was entered correctly.");
print("# Note: The script will skip any lines that do not have all species.");
print("# -------------------------------------");
print("# " + core.getTime() + " Preparing species dictionary...");

specdict = {};
for each in speclist:
	current = each.split(":");
	specdict[current[0]] = current[1];

#print specdict;
print("# -------------------------------------");
print("# " + core.getTime() + " Reading source files and extracting sequence IDs...");
#tmp_seq_dict = {};
main_seq_dict = {};
for spec in specdict:
	main_seq_dict[spec] = core.fastaGetDict(os.path.join(seqdir,specdict[spec]));

#main_seq_dict = {};
#for spec in tmp_seq_dict:
#	main_seq_dict[spec] = {};
#	for title in tmp_seq_dict[spec]:
#		#new_title = title[1:title.index(" ")];
#		main_seq_dict[spec][new_title] = tmp_seq_dict[spec][title];

#del tmp_seq_dict;

print("# -------------------------------------");
count = core.getFileLen(infilename);
print("# " + core.getTime() + " Combining", count, "orthologs...");

i = 0;
numbars = 0;
donepercent = [];

nonorth = 0;

for line in open(infilename):
	numbars, donepercent = core.loadingBar(i, count, donepercent, numbars);

	tmpline = line.replace("\n","").split("\t");
	if i == 0:
		numspec = len(tmpline);
		i = i + 1;
		continue;

	if len(tmpline) != numspec:
		print(line);
		nonorth = nonorth + 1;
		continue;

	sn = 0;
	finalseqs = {};
	for gid in tmpline:
		specid = gid[:6];
		#specid = idmap[specid]

		for title in main_seq_dict[specid]:

			if gid in title:
				curseq = main_seq_dict[specid][title];
				if remstart == 1 and curseq[0] == "M":
					curseq = curseq[1:];

				new_title = title.replace("_", " ", 1);
				finalseqs[new_title] = curseq;

				if sn == 0:
					tid = title[title.index("_")+1:title.index(" ")];
					outfilename = os.path.join(outdir, tid + ".fa");
					#outfilename = outdir + tid + ".fa";
					#print outfilename;

		sn = sn + 1;

	if len(finalseqs) == numspec:
		core.filePrep(outfilename, "");
		for title in finalseqs:
			if delim:
				outtitle = title[:title.index(delim)];
			else:
				outtitle = title;
			core.writeSeqOL(outfilename,finalseqs[title],outtitle);

	i = i + 1;

pstring = "100.0% complete.";
sys.stderr.write('\b' * len(pstring) + pstring);
print("\n# " + core.getTime() + " Done!");
print(nonorth, "lines skipped.");
print("# =======================================================================");
