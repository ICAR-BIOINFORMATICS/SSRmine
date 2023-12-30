import pandas as pd
import warnings
import sys
warnings.filterwarnings("ignore")

path = sys.argv[1]

print("\n\n*************************************************************************\n\n\
 	\t\tWelcome to \n\n\
   \t\tSSRmine: SSR Mining Tool\n\n\
*************************************************************************\n\n\
    Documentation Link: https://github.com/ICAR-BIOINFORMATICS/SSRmine \n\n\
*************************************************************************\n\n\
    Institute: ICAR-NIPB & ICAR-IASRI, New Delhi, India\n\n\
    Developed by: \n\
    \tDr. Shbana Begam, ICAR-NIPB, New Delhi, India \n\
    \tDr. Samarth Godara, ICAR-IASRI, New Delhi, India\n\n\
*************************************************************************\n\n")

def extract_ssr_info(t, thresh, seq, id, fp, rp, idx):
    ssr_start=0
    i=0
    while i < len(seq)-(t*2):
        if seq[i:i+t] == seq[i+t:i+(t*2)]:
            ssr_type = seq[i:i+t]
            ssr_start = i
            j = i
            while j < len(seq)-(t*2):
                if seq[j:j+t] == seq[j+t:j+(t*2)]:
                    j = j+t
                else:
                    break
            ssr_end = j+t
            i=j
            if ssr_end - ssr_start >= (thresh*t):
                id_list.append(id)
                ssr_type_list.append(ssr_type)
                ssr_start_list.append(ssr_start+idx)
                ssr_end_list.append(ssr_end+idx)
                motif_repeats_list.append((ssr_end - ssr_start)/t)
                try:
                    ssr_for_list.append(seq[ssr_start-fp:ssr_start])
                    ssr_rev_list.append(seq[ssr_end:ssr_end+rp])
                except:
                    ssr_for_list.append(seq[ssr_start:ssr_start+1])
                    ssr_rev_list.append(seq[ssr_end-1:ssr_end])
                i=j-t
            #i=j+t
        else:
            i+=1

def extract_info(id, seq, fp, rp, idx):
    global err
    for ssr_t in range(1,7):
        try:
            extract_ssr_info(ssr_t,thresh_dict[ssr_t], seq, id, fp, rp, idx)
        except:
            print("Error processing ID = ", id, "MOTIF Length = ", ssr_t)
            err+=1

err=0
fp=150
rp=150
seq= ''
id=''
idx=0
id_list = []
ssr_type_list = []
ssr_start_list = []
ssr_end_list = []
motif_repeats_list = []
ssr_for_list = []
ssr_rev_list = []
total_chars = len(open(path, 'r').read())
read_chars = 0
thresh_dict={1:10, 2:6, 3:5, 4:5, 5:5, 6:5}

print("Opening file :", path)

def preprocess(p):
	new_path = p[:p.rfind('.')]+"_pre."+p[p.rfind('.')+1:]
	print("Preprocessing file : ", p)
	print("Saving file : ", new_path)
	
	with open(p, "r") as input_file:
		with open(new_path, "w") as output_file:
			for input_line in input_file:
				if input_line[0]!='>':
					output_file.write(input_line.upper())
				else:
					output_file.write(input_line)
	return new_path
	
user_input = input("Do you want to perform preprocessing? (yes/no): ")
if user_input.lower() == "yes":
	path = preprocess(path)

print("\n\nProcessing completed (%)")

with open(path, "r") as f:
  for line in f:
    read_chars += len(line)
    if line[0]=='>' or len(seq)>10000:
        if id=='':
            id = line[1:-1]
        else:
            print("\b\b\b\b\b", end='')
            print(round((read_chars/total_chars)*100, 2), end='')
            extract_info(id, seq, fp, rp, idx)
            if len(seq)<=10000:
              id = line[1:-1]
              seq=''
              idx=0
       	    else:
       	      seq=seq[-50:]
       	      idx+=9950
    else:
        seq=seq+line[:-1]

extract_info(id, seq, fp, rp, idx)

print("\nCalculating SSR occurances...")

print("Error counts:", err)

ssr_df = pd.DataFrame()

ssr_df['ID'] = id_list
ssr_df['SSR_type'] = ssr_type_list
ssr_df['SSR_start'] = ssr_start_list
ssr_df['SSR_end'] = ssr_end_list
ssr_df['MOTIF_repeats'] = motif_repeats_list
ssr_df['Forward_sequence'] = ssr_for_list
ssr_df['Reverse_sequence'] = ssr_rev_list

count_list =  ssr_df.SSR_type.value_counts()

c_vals = count_list.values
c_keys = count_list.index.tolist()

def fail_motif(m):

    if len(m)%2==0:
        fh = m[:int(len(m)/2)]
        sh = m[int(len(m)/2):]

        if fh==sh:
            return True

        if len(m)%3==0 :
            ft = m[:int(len(m)/3)]
            st = m[int(len(m)/3):int(len(m)*(2/3))]
            tt = m[int(len(m)*(2/3)):]

            if ft==st and st==tt:
                return True

    elif len(m)>1:
        for j in range(1,len(m)):
            if m[0]!=m[j]:
                return False

        return True

    return False

for i in range(len(c_keys)):
    if fail_motif(c_keys[i]):
        c_vals[i]=0

ssr_dict = {c_keys[i]: c_vals[i] for i in range(len(c_keys))}

def get_ssr_rep(rec):
    return ssr_dict[rec['SSR_type']]

ssr_df['Occurances'] = ssr_df.apply(get_ssr_rep, axis=1)

ssr_df = ssr_df[ssr_df.Occurances > 0]

new_path = path[:path.rfind('.')]+"_ssr.csv"

print("\nSaving output file :", new_path)

ssr_df.to_csv(new_path, index=False)

