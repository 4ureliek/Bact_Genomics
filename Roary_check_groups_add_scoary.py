#!/usr/bin/env python3
#----------------------------------------------------------------------------
# Author  :  Aurelie Kapusta
# email   :  4urelie.k@gmail.com
#----------------------------------------------------------------------------
import sys
import argparse
import os
import re
import csv

chlog = """
#	v1.0 = 2019 10 XX - TO DO
#		
"""
version = "1.0"
scriptname = "Roary_check_groups_add_scoary.py"

#Deal with formatting the arguments & stuff:
parser = argparse.ArgumentParser(
usage=scriptname+" -i -r [-f] [-d] [-h] [-v]",
description="From the output of Roary_check_groups.py, now load some columns from Scoary.",
#formatter_class=argparse.RawDescriptionHelpFormatter #Use this if formatting in the description
)
parser._optionals.title = "Arguments" #I prefer avoiding positional arguments
parser.add_argument('-g', metavar='REQ', type=str, action='store', help="Output file of Roary_check_groups.py (extract.*.csv.check_groups.choices.tab).", required=True)
parser.add_argument('-i', metavar='REQ', type=str, action='store', help="Config file, 2 columns tab delimited: \"ID_of_the_csv_in_-g_file \\t full_path_to_corresponding_scoary_folder\". If you have headers, use # at the beginning of the line.", required=True)
parser.add_argument('-c', metavar='', type=str, action='store', help="Wanted columns from Scoary output: comma separated indexes (0 based).Default: 3,4,5,6", default="3,4,5,6")
parser.add_argument('-v', help="Verbose mode on",action='store_true', default=False) 
parser.add_argument('-d', help="Debugger mode on (will print a bunch more stuff)",action='store_true', default=False) 
parser.add_argument('--version', help="get version",action='version', version=scriptname+' v'+version, default=False) 
args = parser.parse_args()
		
if args.v:
	print(scriptname+" started, v"+version)
	print("With")
	print("   -g: "+args.g)
	print("   -i: "+args.i)
	print("   -c: "+args.c)
	print()

#-------------------------------------------------------------------------------
#--- FUNCTIONS
#-------------------------------------------------------------------------------
def load_config_file():
	try:
		with open(args.i,"r") as fhi:
			for line in fhi:
				line = line.rstrip()
				if re.search(r'^#',line):
					continue		
				dat=line.split("\t")
				if dat[0] not in FILES:
					FILES[dat[0]]=dat[1]
	except EnvironmentError:
		print("ERROR - %s does not exist? Exiting..." % args.i)
		exit()
	return()

#-------------------------------------------------------------------------------
def load_scoary():
	for f in FILES:
		if not os.path.isdir(FILES[f]):
			print("   WARN: %s is not a directory? Skipping" % FILES[f])
		else:			
			#gather the scoary outputs:
			for sc in os.listdir(FILES[f]):
				if not sc.endswith(".results.csv"):
					continue
				if f not in SCOARY:
					SCOARY[f]={}
				phen = re.sub("_.*","",sc)
				if phen not in PHEN:
					PHEN[phen]=1
				if phen not in SCOARY[f]:
					SCOARY[f][phen] = {}		
							
				with open(os.path.join(FILES[f], sc)) as fhi:
					csv_reader = csv.reader(fhi, delimiter=',')
					line_count = 0
					for dat in csv_reader:
						g = dat[0]
						if g == "Gene":
							continue
						if g not in SCOARY[f][phen]:
							SCOARY[f][phen][g] = []	
						for c in COLS:
							SCOARY[f][phen][g].append(dat[int(c)])

	return()



#-------------------------------------------------------------------------------
def edit_group_choice_file(out,phen):
	"""
	Insert the values from Scoary:
		SCOARY[csv_id][phenotype][group_id]
	One file per medical phenotype found
	"""
	h=[] #headers
	i = 0	
	with open(args.g,"r") as fhi, open(out,"w") as fho:
		for line in fhi:
			line = line.rstrip()
			v=line.split("\t")
			fn=int((len(v)+2)/4)
			p=2
			if i == 0:		
				#header: insert len(COLS) times the header value
				for f in range(0,fn):
					toinsert = []
					#set what to insert
					for t in range(0,len(COLS)):
						toinsert.append(v[p-1])
					#now insert it
					for toins in toinsert:
						v.insert(p,toins)
					h.append(v[p]) #save files in same order
					p+=len(COLS)+4
			elif i == 1:	
				toinsert = []
				for c in COLS:
					scoaryid = "scoary_colum_%i" % int(c)
					v.insert(p,scoaryid)
					v.append(scoaryid)
					p+=1
			else:
				#same thing as the header for the positions
				for f in h:
					if f not in SCOARY: #no data for this phenotype
						for t in range(0,len(COLS)):
							v.insert(p,".")
							p+=1
					else:
						gi=p-4
						if gi < 0:
							gi = 0
						g=v[gi]
						if phen not in SCOARY[f]:
							SCOARY[f][phen]={}
						if g not in SCOARY[f][phen]:
							SCOARY[f][phen][g]=[]
							for t in range(0,len(COLS)):
								SCOARY[f][phen][g].append(".")				
						for val in SCOARY[f][phen][g]:
							v.insert(p,val)
							p+=1
					p+=4				
			i+=1
			data="\t"
			data = data.join(v)
			fho.write(data+"\n")

	return()



#-------------------------------------------------------------------------------
#--- MAIN
#-------------------------------------------------------------------------------
#Load fa files of interest
if args.v:
	print("Loading config file: "+args.i)
FILES={}
load_config_file()
if args.d:
	print(FILES)
print()

#Load the scoary outputs to insert them in the choice file
PHEN = {} #save a record of phenotypes in case several scoary outputs
COLS = args.c.split(",") #columns to add
SCOARY = {}
if args.v:
	print("Loading Scoary data")
load_scoary()
if args.d:
	for f in SCOARY:
		print("   Patient phenotypes for file: %s" % f)
		for phen in SCOARY[f]:
			 print("      "+phen)
print()

#Now insert the data
if args.v:
	print("Inserting SCOARY data in: "+args.g)
for phen in PHEN:
	out=re.sub("tab",phen+".tab",args.g)
	edit_group_choice_file(out,phen)
if args.v:
	print("   -> "+out)
print()




