#!/usr/bin/env python3
#----------------------------------------------------------------------------
# Author  :  Aurelie Kapusta
# email   :  4urelie.k@gmail.com
#----------------------------------------------------------------------------
import sys
import argparse
import os
import random
from Bio import SeqIO
import re

chlog = """
#	v1.0 = 2019 09 19
#	v2.0 = 2019 09 26
#		The printing was wrong -> wasn't making the columns
#		Fixing that required loading things differently
#		Also add -f (to force rerunning things)
#	v3.0 = 2019 10 18
#		Weird behavior with -max_target_seq from blast can change the hits 
#          => remove it and just load the best hit only
#		add -b option for blast path
#		bug fix for -f option
#       print all groups, not just 3 max, and add -n option
#       add -a option, to add the Roary group descriptions
#	v3.1 = 2019 10 20
#		Bug fix if no best hit
#	v3.2 = 2019 10 20
#		Add option -m

#	TO DO - put back the ".csv"

"""
version = "3.1"
scriptname = "Roary_check_groups.py"

#Deal with formatting the arguments & stuff:
parser = argparse.ArgumentParser(
usage=scriptname+" -i -r [-f] [-d] [-h] [-v]",
description="Will check groups from different roary runs to find which ones go together.",
#formatter_class=argparse.RawDescriptionHelpFormatter #Use this if formatting in the description
)
parser._optionals.title = "Arguments" #I prefer avoiding positional arguments
parser.add_argument('-i', metavar='REQ', type=str, action='store', help="""The *.extract folders (extracted by group) in -i directory (outputs of Prokka_Roary_extract.pl). Should be named with the Roary run ID (e.g. "gene_presence_absence_90.csv"), because that folder name will be the "ID" in the output of this script.""", required=True)
parser.add_argument('-r', metavar='REQ', type=str, action='store', help="""The csv file name to start from (what is before the ".extract" in the -i directory for that roary run). For each group of that csv, the corresponding groups in the other XXX.extract folders will be listed.""", required=True)
parser.add_argument('-m', help="Comma separated percentage and minimum number of sequences of a group to extract (Default: 10,10)",action='store_true', default="10,10") 
parser.add_argument('-f', help="To force re running the extraction & blasts even if the files are there",action='store_true', default=False) 
parser.add_argument('-a', help="To add the descriptions in the final output",action='store_true', default=False) 
parser.add_argument('-b', metavar='', type=str, action='store', help="Path to blast binaries if not in PATH", default=False)
parser.add_argument('-n', metavar='', type=int, action='store', help="To set a maximum number of corresponding groups to print", default=False)
parser.add_argument('-d', help="Debugger mode on (will print a bunch of stuff)",action='store_true', default=False) 
parser.add_argument('-v', help="Verbose mode on",action='store_true', default=False) 
parser.add_argument('--version', help="get version",action='version', version=scriptname+' v'+version, default=False) 
args = parser.parse_args()

#Clean up the / at the end if it's there
if args.i.endswith("/"):
	fp = args.i.strip("/")
else:
	fp = args.i
if args.b.endswith("/"):
	bp = args.b.strip("/")
else:
	bp = args.b

if args.v:
	print(scriptname+" started, v"+version)
	print("With")
	print("   -i: "+fp)
	print("   -r: "+args.r)
	print()

#-------------------------------------------------------------------------------
#--- FUNCTIONS
#-------------------------------------------------------------------------------
def load_dirs_and_files():
	fadirs = []
	fafiles = []
	groups = []
	for d in os.listdir(fp):
		if os.path.isdir(os.path.join(fp, d)):
			csv = d.replace(".extract", "", 1)
			if csv == args.r:
				for f in os.listdir(fp+"/"+d):
					if f.endswith(".fa"):
						groups.append(fp+"/"+d+"/"+f)
			else:
				faname=outdir+"/"+csv+".fasta"
				fafiles.append(faname)
				if not os.path.isfile(faname) or args.f:
					cmd = "cat "+fp+"/"+d+"/*.fa > "+faname
					if args.d:
						print("   with command: "+cmd)
					os.system(cmd)
				elif args.d and os.path.isfile(faname) and not args.f:
					print("    fasta output exists and -f not set: skipping")

#Even for just -d, that's just huge lists to print, remove				
# 	if args.d:
# 		print("groups:")
# 		print(groups)
# 		print("fasta files:")
# 		print(fafiles)			
	return(groups,fafiles)

#-------------------------------------------------------------------------------
def load_groups_description():
	for d in os.listdir(fp):
		efp = os.path.join(fp, d)
		if not os.path.isdir(efp):
			continue
		csv=d.replace(".csv.extract","",1)
		if not csv in desc:
			desc[csv] = {}
		for fa in os.listdir(efp):
			if not fa.endswith(".fa"):
				continue
			faefp = os.path.join(efp, fa)
			records = list(SeqIO.parse(faefp, "fasta"))
			for rec in records:
				fields=rec.description.split("\t") #rec.descriptio is the full header
				seqid = fields[0]
				seqid=re.sub("##.*","",seqid)
				try:
					thisdesc = fields[1]
				except IndexError:
					thisdesc = "na"
				thisdesc=re.sub(".*ROARY_DESC=","",thisdesc)
				if seqid in desc[csv]:
					continue
				desc[csv][seqid]=thisdesc
	return(desc)

#-------------------------------------------------------------------------------
def make_blast_databases(fafiles):
	for fa in fafiles:
		if not os.path.isfile(fa+".nhr") or args.f:
			if args.b:
				cmd =bp+"/makeblastdb"
			else:
				cmd ="makeblastdb"
			cmd = cmd+" -in "+fa+" -dbtype nucl -out "+fa+" &> "+fa+".makeblastdb.log"
			if args.d: 
				print("   for file: "+fa)
				print("   with command: "+cmd)
			os.system(cmd)
		elif args.d and os.path.isfile(fa+".nhr") and not args.f:
			print("   blastdb already exists for file: "+fa)
	return()

#-------------------------------------------------------------------------------
def set_random_ids(start, end, num): 
	res = [] 
	for j in range(num): 
		index = random.randint(start, end)
		check=0
		while index in res:
			index = random.randint(start, end)
			check+=1
			if check == 100:
				break
		res.append(index) 
	return(res) 

#-------------------------------------------------------------------------------
def extract_X_percent_from_fasta(fa,perc,min):
	(path,name)=os.path.split(fa)
	faname=name.replace(".fa",".10.fa",1)
	fa_extract=os.path.join(outdirfa, faname) 

	#load fasta file
	records = list(SeqIO.parse(fa, "fasta"))
	#check how many to extract
	recnb=len(records)
	toextract=int(recnb/100*perc)
	if recnb < min:
		toextract = recnb
	elif toextract < min:
		toextract = min		
	if not os.path.isfile(fa_extract) or args.f:
		if args.d:
			print("   has %i records =>  will blast %i (if > 10, randomly extracted)" % (recnb,toextract))		
		#Get random indexes 
		rand_list=set_random_ids(0, recnb-1, toextract)	#do this only if toextract < recnb		
		with open(fa_extract,"w") as fho:
			for i, record in enumerate(SeqIO.parse(fa, "fasta")):
				if i in rand_list:
					fho.write(record.format("fasta"))
		if args.d: 
			print("   %i sequences extracted in: %s" % (toextract,fa_extract))
	elif args.d and os.path.isfile(fa_extract) and not args.f: 
			print("    %i sequences already extracted and -f not set: skipping" % toextract)		
	return(fa_extract)
	
#-------------------------------------------------------------------------------
def run_blastn(groups,fafiles):
	blastout=[]
	for gall in groups:
		if args.d: 
			print()
			print("   For file: "+gall)
		g = extract_X_percent_from_fasta(gall,args.m)
		for fa in fafiles:
			(gpath,gname)=os.path.split(g)
			(fapath,faname)=os.path.split(fa)
			fanameshort=faname.replace(".csv.fasta", "", 1)
			outdirblast=outdir+"/"+fanameshort+".blastout"
			if not os.path.isdir(outdirblast):
				os.makedirs(outdirblast)
			out=outdirblast+"/"+gname+".out"
			if not os.path.isfile(out) or args.f:
				eval="10e-20"
				blastfmt="6 qseqid qstart qend qlen qcovs qcovhsp sseqid sstart send sstrand slen evalue bitscore score length pident nident mismatch ppos positive gapopen gaps"
				if args.b:
					cmd =bp+"/blastn"
				else:
					cmd ="blastn"
				cmd = cmd+" -query "+g+" -out "+out+" -db "+fa+" -evalue "+eval+" -outfmt \""+blastfmt+"\" &> "+out+".log"
				os.system(cmd)
				if args.d: 
					print("   with command: "+cmd)
			elif args.d and os.path.isfile(out) and not args.f: 
				print("    blast output exists and -f not set: skipping")		
			blastout.append(out)
	return(blastout)
		
#-------------------------------------------------------------------------------
def parse_blast_outputs(blastout):
	choices = {}
	maxiter = {}
	for f in blastout:
		if args.d:
			print("   For file: "+f)
		(fpath,fname)=os.path.split(f)
		grpi=fname.replace(".10.fa.out","",1)
		#if group seen for the first time, define the dicts:
		if grpi not in choices:
			choices[grpi]={}
			maxiter[grpi]=1
		#get blastoutput (bo) name
		(path,bo)=os.path.split(fpath)
		prevquery=""
		try:
			with open(f,"r") as fhi:
				counts = {}
				tot=0
				for line in fhi:
					line = line.rstrip()
					dat=line.split("\t")
					(grpo,splo)=dat[6].split("##")
					#sorted by bitscores so I can just keep only the top one
					if dat[0] == prevquery:
						continue
					else:
						prevquery = dat[0]
						if grpo in counts:
							counts[grpo]+=1
						else:
							counts[grpo]=1
						tot+=1
				
				#Now sort counts to get the corresponding group
				#Also keep the % that it represents
				iter = 1
				for g, c in sorted(counts.items(), key=lambda item: item[1], reverse=True):
					p = c / tot *100
					if args.d:
						print("      group %s: %i hits against %s => %.2f %%" % (grpi, c, g, p))
					if iter not in choices[grpi]:
						choices[grpi][iter]={}
					choices[grpi][iter][bo]="%s\t%i\t%.2f" % (g, c, p)
					if maxiter[grpi] < iter:
						maxiter[grpi]=iter
					if args.n and iter == args.n:
						break
					iter+=1
		except IOError:    
			print("Can't open file to read:" , file,"- skip")
			continue

	#initialize values to get same amount of 'iter' for all bo
	for g in choices:
		for i in range(1, maxiter[g]+1):
			for f in blastout:
				(fpath,fname)=os.path.split(f)
				(path,bo)=os.path.split(fpath)
				if i not in choices[g]:
					choices[g][i]={}
				if bo not in choices[g][i]:
					choices[g][i][bo]=".\t.\t."
	return(choices,maxiter)

#-------------------------------------------------------------------------------
def print_all_group_choices(dir,out,choices,maxiter):
	#choices, 3 columns for each: choices[grpi][iter][bo]
	#Add descriptions if args.a
	#files=csvs_gene_presence_absence_90.csv.check_groups/gene_presence_absence_95.blastout/aknK.fa.outchoice.tab		
	with open(out, 'w') as fho:
		#first print the header = the ids
		fho.write(args.r)
		if args.a:
			fho.write("\t"+args.r)
		for bo in os.listdir(dir):
			if os.path.isdir(os.path.join(dir, bo)) and bo.endswith(".blastout"):
				boid = bo.replace(".blastout","",1)
				fho.write("\t%s\t%s\t%s" % (boid,boid,boid))		
				if args.a:
						fho.write("\t%s" % boid)	
		if args.a:
			fho.write("\ninitial_group\tinitial_group_desc\tcorresponding_group\tcorresponding_group_desc\tnumber_with_that_hit\tpercentage_with_that_hit\n")
		else:
			fho.write("\ninitial_group\tcorresponding_group\tnumber_with_that_hit\tpercentage_with_that_hit\n")
	
		for gr in choices:
			for i in range(1, maxiter[gr]+1):
				#loop on the bo files same order:
				line=""
				for bo in os.listdir(dir): 
					if os.path.isdir(os.path.join(dir, bo)) and bo.endswith(".blastout"):
						if args.a:
							#split the columns to insert the description
							boid = bo.replace(".blastout","",1)
							dat=choices[gr][i][bo].split("\t")
							if dat[0] == ".":
								thisdesc = "."
							else:
								thisdesc = desc[boid][dat[0]]
							line=line+"\t"+dat[0]+"\t"+thisdesc+"\t"+dat[1]+"\t"+dat[2]
						else:
							line=line+"\t"+choices[gr][i][bo]
				if args.a:
					id=args.r.replace(".csv","",1) 
					fho.write(gr+"\t"+desc[id][gr]+line+"\n")
				else:
					fho.write(gr+line+"\n")
	return()

#-------------------------------------------------------------------------------
#--- MAIN
#-------------------------------------------------------------------------------
#create output directory if needed
outdir=fp+"_"+args.r+".check_groups"
if not os.path.isdir(outdir):
	os.makedirs(outdir)

#Load files from the subdirs
if args.v:
	print("Loading group IDs from "+args.i+args.r+".extract")
	print("And creating fasta files for all other Roary runs")
(groups,fafiles) = load_dirs_and_files()
if args.v:
	print()

#Load descriptions -f -a is set
desc = {}
if args.a:
	if args.v:
		print("Loading group descriptions from all .extract folders")
	desc = load_groups_description() #keys=csv name, and then seqid
	if args.v:
		print()

#Now make the bastable databases
if args.v:
	print("Creating blastable databases")
make_blast_databases(fafiles)
if args.v:
	print()

if args.v:
	print("Blasting groups from -r file against other runs")
outdirfa=outdir+"/extracted_fa"
if not os.path.isdir(outdirfa):
	os.makedirs(outdirfa)
blastout = run_blastn(groups,fafiles)
if args.v:
	print()

if args.v:
	print("Parsing the blast outputs")
(choices,maxiter) = parse_blast_outputs(blastout)
if args.v:
	print()

if args.v:
	print("Printing the results (including pre-existing ones)")
	if args.a:
		print("(adding roary descriptions)")
finalout = outdir+".choices.tab"
print_all_group_choices(outdir,finalout,choices,maxiter)
if args.v:
	print("   -> "+finalout)
if args.v:
	print()




