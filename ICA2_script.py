#/usr/bin/python3
import os,sys,subprocess,string,re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

print("Welcome!")
# get user's input
prot_fam = input("Please enter your protein family: \n")
tax_group = input("Please enter your taxonomy group: \n")
print("\nYou have provided the following information:\n\tProtein Family: ",prot_fam,"\n\tTaxonomy Group: ",tax_group)
def check():
    answer = input("Is this right?\n")
    if answer.upper()[0] != "Y":
        print("Alright, exiting... You can restart and have another try!")
        exit()
    else:
        print("Ok, thank you!")
check()

# obtain the relevant protein sequence data
fa_name = prot_fam+"_"+tax_group+".fasta"
cmd = "esearch -db protein -query \"("+prot_fam+"[prot])"+" AND ("+tax_group+"[organism])\"|efetch -format fasta > "+fa_name
print("Now searching for protein sequences...\n")
subprocess.call(cmd,shell=True)
print("Complete!")
fasta = open(fa_name)
fasta_contents = fasta.read()
fasta.close()
if fasta_contents.count(">")>0: 
    print("Your fasta file "+fa_name+" is ready, "+str(fasta_contents.count(">"))+" sequences found.\n")
else:    
    print("No sequence found, please check your protein family and/or taxonomy group again.\nExiting...")
    os.remove(fa_name)
    exit()
def check1():
    print("Recommend: The sequences set shouldn't have more than 1000 sequences.")
    answer = input("Do you want to continue?\n")
    if answer.upper()[0] != "Y":
        print("Alright, exiting... You can restart and have another try!")
        exit()
    else:
        print("Ok, thank you!")
check1()

# determine and plot the level of conservation between the protein sequences
cmd = "emma -sequence "+fa_name+" -outseq aligned_"+fa_name+" -dendoutfile clustalW.tree > align.log"
print("Processing multiple sequences alignment...\n")
subprocess.call(cmd,shell=True)
subprocess.call("infoalign -sequence aligned_"+fa_name+" -nousa -outfile info_alignment.txt", shell=True)
print("Complete! Your output file is aligned_"+fa_name+", and you can check for more information in info_alignment.txt and align.log\n")



def choose():
    everything = 42
    # just an empty funtion to avoid error in choose1()
def choose1():
    # I tried dealing with the info_alignment.txt in python, unfortunately it use multiple space instead of some tab, even using replace(" ","\t") leads to "\t\t\t" in the header line... so I turn to bash
    cmd = "sed 's/ /\t/' info_alignment.txt| sed 's/\t$//' > tmp.txt"
    subprocess.call(cmd,shell=True)
    df = pd.read_csv('tmp.txt', sep='\t')
    # acturally there is still a little bug: the column of name have the colname "#", and the column of sequence length have the colname "Name      SeqLen", anyway, our manipulation of the dataframe is just for filtering, if the user don't check the code and it will be okay...
    os.remove('tmp.txt')
    print("\nWe have identical/different/similar sequences (similar means score>0 in the comparison matrix but not identical)\nWhich indicator would you like to filter by[1/2/3/4]:\n\t1 identical sequences percentage\n\t2 different sequences percentage\n\t3 don't care the indicators(take subset from the beginning)\n\t4 back to the previous choose\n")
    while True:
        answer = input("Make a choice!\n")
        if answer == '1':
            df['IdentPer'] = 100-df['% Change']
            print("The identical percentage of your sequences ranges from "+str(df['IdentPer'].min())+" to "+str(df['IdentPer'].max()))
            per = input("\nThe minimum identical percentage you want:\n")
            try:
                if float(per)>=df['IdentPer'].min() and float(per)<=df['IdentPer'].max():
                    pass
            except:
                f'Not valid n!'
            ids = df[ df['IdentPer']>=float(per) ]["#"]
            id_list = ids.tolist()
            # the item in id_list don't have the ">" and "\n"
            id_tmp = ['>' + item for item in id_list]
            target = '\n,'.join(id_tmp).split(sep=",")
            # as the sequences are aligned, each sequence have same number of lines in fasta file
            a_fasta = open('aligned_'+fa_name)
            a_fasta_contents = a_fasta.read()
            a_fasta.close()
            nline = a_fasta_contents.count("\n")/a_fasta_contents.count(">")
            input_f = open('aligned_'+fa_name, 'r') 
            input_f_con = input_f.readlines()
            output_f = open('subset_aligned_'+fa_name, 'w')
            output = []
            for line in input_f_con:
                if line in target:
                    output.append(line)
                    nr =  nline - 1
                while nr>0 :
                    next_line = input_f_con[int(input_f_con.index(line)+nline-nr)]
                    output.append(next_line)
                    nr -= 1
            output_f.write(''.join(output))
            input_f.close()
            output_f.close()
            cmd1 = "plotcon -sequences subset_aligned_"+fa_name+" -winsize 4 -graph x11"
            cmd2 = "plotcon -sequences subset_aligned_"+fa_name+" -winsize 4 -graph png -goutfile \"ident>"+per+"%_conservation\""
            subprocess.call(cmd1,shell=True)
            subprocess.call(cmd2,shell=True)
            print("Your output file is ident>"+per+"%_conservation.1.png\n")
            break
        elif answer == '2':
            df['DiffPer'] = df['Differ']/df['Name        SeqLen']*100
            print("The different percentage of your sequences ranges from "+str(df['DiffPer'].min())+" to "+str(df['DiffPer'].max()))
            per = input("\nThe maximum different percentage you want:\n")
            try:
                if float(per)>=df['DiffPer'].min() and float(per)<=df['DiffPer'].max():
                    pass
            except:
                f'Not valid n!'
            ids = df[ df['DiffPer']<=float(per) ]["#"]
            id_list = ids.tolist()
            id_tmp = ['>' + item for item in id_list]
            target = '\n,'.join(id_tmp).split(sep=",")
            a_fasta = open('aligned_'+fa_name)
            a_fasta_contents = a_fasta.read()
            a_fasta.close()
            nline = a_fasta_contents.count("\n")/a_fasta_contents.count(">")
            input_f = open('aligned_'+fa_name, 'r')
            input_f_con = input_f.readlines()
            output_f = open('subset_aligned_'+fa_name, 'w')
            output = []
            for line in input_f_con:
                if line in target:
                    output.append(line)
                    nr =  nline - 1
                while nr>0 :
                    next_line = input_f_con[int(input_f_con.index(line)+nline-nr)]
                    output.append(next_line)
                    nr -= 1
            output_f.write(''.join(output))
            input_f.close()
            output_f.close()
            cmd1 = "plotcon -sequences subset_aligned_"+fa_name+" -winsize 4 -graph x11"
            cmd2 = "plotcon -sequences subset_aligned_"+fa_name+" -winsize 4 -graph png -goutfile \"differ<"+per+"%_conservation\""
            subprocess.call(cmd1,shell=True)
            subprocess.call(cmd2,shell=True)
            print("Your output file is differ<"+per+"%_conservation.1.png\n")
            break
        elif answer == '3':
            a_fasta = open('aligned_'+fa_name)
            a_fasta_contents = a_fasta.read()
            a_fasta.close()
            n = input("Please give me the number of sequences in the subset:\n")
            try:
                if x>0 and x<fasta_contents.count(">"):
                    pass
            except:
                f'Not valid n!'
            lines = int(n)*a_fasta_contents.count("\n")/a_fasta_contents.count(">")
            with open('subset_aligned_'+fa_name, 'a') as sub:
                sub.write('\n'.join(a_fasta_contents.split()[:int(lines)]))
            cmd1 = "plotcon -sequences subset_aligned_"+fa_name+" -winsize 4 -graph x11"
            cmd2 = "plotcon -sequences subset_aligned_"+fa_name+" -winsize 4 -graph png -goutfile \""+n+"_conservation\""
            subprocess.call(cmd1,shell=True)
            subprocess.call(cmd2,shell=True)
            print("Your output file is "+n+"_conservation.1.png")
            break
        elif answer == '4':
            print("Alright, hope this time you'll choose what you really want to do...")
            choose()
            break
def choose():
    print("Now you can choose from [1/2/3/4]:\n\t1 plot conservation level directly and save\n\t2 choose a subset to plot and save\n\t3 skip and go scanning motifs of interest\n\t4 quit (files will be saved)\n")
    while True:
        answer = input("What do you want next?\n")
        if answer == '1':
            cmd1 = "plotcon -sequences aligned_"+fa_name+" -winsize 4 -graph x11"
            cmd2 = "plotcon -sequences aligned_"+fa_name+" -winsize 4 -graph png -goutfile \"conservation\""
            subprocess.call(cmd1,shell=True)
            subprocess.call(cmd2,shell=True)
            print("Your output file is conservation.1.png")
            break
        elif answer == '2':
            choose1()
            break
        elif answer == '3':
            print("Ok, let's go!")
            break
        elif answer == '4':
            print("Ok, thank you, bye!")
            exit()
        else:
            print("Please choose one from [1/2/3/4]...")
choose()

# scan protein sequence(s) of interest with motifs from the PROSITE database, to determine whether any known motifs (domains) are associated with this subset of sequences
print("\nGenerating the consensus sequence from aligned sequences...\n")
cmd = "cons -sequence aligned_"+fa_name+" -outseq consensus.fasta -name consensus"
subprocess.call(cmd,shell=True)
print("Complete! Your consensus sequence is consensus.fasta\n")
print("Scanning consensus sequence for motifs...\n")
cmd = "patmatmotifs -sequence consensus.fasta -rformat2 motif -outfile "+prot_fam+"_"+tax_group+".motif"
subprocess.call(cmd,shell=True)
print("Complete!")
f = open(prot_fam+"_"+tax_group+".motif")
f_con = f.read()
f.close()
motif = re.findall("Motif: (\w+)",f_con)
sequence = re.findall("Sequence: (?!consensus)(\w+)",f_con)
motif_dict = dict(zip(motif,sequence))
for key,value in motif_dict.items():
    print("Found motif: "+key+" (sequence: "+value+") in consensus sequence.\n")
df = pd.DataFrame()
df['motif'] = motif
df['sequences'] = sequence
start = re.findall("residues (\w+)-",f_con)
end = re.findall("->(\w+)",f_con)
df['start'] = start
df['end'] = end
END = re.findall("to: (\w+)",f_con)
numend = int(''.join(END))
print("Now plotting motif positions...\n")
fig,ax=plt.subplots()
ax.axhline(y=0, xmin=0, xmax=numend, color='black')
rect = patches.Rectangle((int(df['start']), -1), int(df['end'])-int(df['start']), 2, linewidth=1, edgecolor='r', facecolor='none')
ax.add_patch(rect)
ax.set_xlim(0,numend)
ax.set_ylim(-1,1.3)
plt.title('motif position')
plt.text(int(df['start']), 1.2, df['motif'].to_string(index=False, header=False),ha='center',size='x-small')
plt.text(int(df['start']), 1.1, df['sequences'].to_string(index=False, header=False),ha='center',size='x-small')
plt.savefig("motif_position.png",transparent=True)
plt.show()
print("Plot saved as motif_position.png\n")
print("More detail information of motifs is in "+prot_fam+"_"+tax_group+".motif\n")
print("Thank you for using! Bye!")
