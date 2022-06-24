import sys
import pandas as pd
import gzip

argv = sys.argv[1:]
hapfile_path = argv[0]
pgsnp_path = argv[1]
vcf_path = argv[2]
chromosome = argv[3]
viewpointpos = int(argv[4])

print(chromosome)

def recalculate_all(hap1, hap2):
    allpositions = pd.concat([hap1['pos'], hap2['pos']], axis=0).drop_duplicates()
    for position in allpositions:
        searchID1 = hap1[hap1['pos']==position]['ID']
        searchID2 = hap2[hap2['pos']==position]['ID']
        
        #check for the position with the searchID how strong it's linked to the current haplotype 1
        if len(hap1[hap1['pos']==position])==1:
            #extract all lines with where ID1 matches the searchID
            links = hapfile[hapfile['ID1'].isin(searchID1)]
            #extract all of these where the other ID is in either haplotype (split)
            linkstohap1 = links[links['ID2'].isin(hap1['ID'])]
            linkstohap2 = links[links['ID2'].isin(hap2['ID'])]

            #get according weight/plexity values
            weight_1 = linkstohap1['weight'].sum()
            ambig_1 = linkstohap2['weight'].sum()
            plexity_1 = len(linkstohap1)
            amplexity_1 = len(linkstohap2)
        else:
            weight_1 = 0
            ambig_1 = 0
            plexity_1 = 0
            amplexity_1 = 0
        
        #check for the position with the searchID how strong it's linked to the current haplotype 2
        if len(hap2[hap2['pos']==position]==1):
            #extract all lines with where ID1 matches the searchID
            links = hapfile[hapfile['ID1'].isin(searchID2)]
            #extract all of these where the other ID is in either haplotype (split)
            linkstohap2 = links[links['ID2'].isin(hap2['ID'])]
            linkstohap1 = links[links['ID2'].isin(hap1['ID'])]
            
            #get according weight/plexity values
            weight_2 = linkstohap2['weight'].sum()
            ambig_2 = linkstohap1['weight'].sum()
            plexity_2 = len(linkstohap2)
            amplexity_2 = len(linkstohap1)
        else:
            weight_2 = 0
            ambig_2 = 0
            plexity_2 = 0
            amplexity_2 = 0
        
        #get final recalculation values for the current snp
        weight = weight_1 + weight_2
        ambig = ambig_1 + ambig_2
        plexity = plexity_1 + plexity_2
        amplexity = amplexity_1 + amplexity_2
        ratio = round(ambig/weight, 3)

        #assign recalculated values depending on hap1 or hap2
        hap1.loc[hap1.pos == position, 'weight'] = weight
        hap1.loc[hap1.pos == position, 'ambig'] = ambig
        hap1.loc[hap1.pos == position, 'plexity'] = plexity
        hap1.loc[hap1.pos == position, 'amplexity'] = amplexity
        hap1.loc[hap1.pos == position, 'ratio'] = ratio
        if (weight_1 > 0 and weight_2 > 0 and(hap1.loc[hap1.pos == position, 'class'].values[0] != "VP")):
            hap1.loc[hap1.pos == position, 'class'] = "Double_link"

        hap2.loc[hap2.pos == position, 'weight'] = weight
        hap2.loc[hap2.pos == position, 'ambig'] = ambig
        hap2.loc[hap2.pos == position, 'plexity'] = plexity
        hap2.loc[hap2.pos == position, 'amplexity'] = amplexity
        hap2.loc[hap2.pos == position, 'ratio'] = ratio
        if (weight_1 > 0 and weight_2 > 0 and(hap2.loc[hap2.pos == position, 'class'].values[0] != "VP")):
            hap2.loc[hap2.pos == position, 'class'] = "Double_link"

    return(hap1, hap2)
            

#Create names
basename = hapfile_path.replace(".hap", "")
hap1_out_file =  "".join(basename.rsplit("temp/", 1)) + "_hap1.txt"
hap2_out_file =  "".join(basename.rsplit("temp/", 1)) + "_hap2.txt"
haps_out_file =  "".join(basename.rsplit("temp/", 1)) + "_haps.txt"
netwtable_file = "".join(basename.rsplit("temp/", 1)) + "_network.txt"
vcf_out_file =   "".join(basename.rsplit("temp/", 1)) + "_variant.vcf"

hapfile = pd.read_csv(hapfile_path, sep="\t", lineterminator="\n", header = None)
hapfile.columns = ['var1', 'pos1', 'var2', 'pos2', 'weight']

pgsnp = pd.read_csv(pgsnp_path, sep="\t", lineterminator="\n", header = None, skiprows=1)
pgsnp.columns = ['chr', 'pos', 'pos+1', 'vars', 'var#', 'reads', '0,0']

#If there's no viewpoint override, autodefine best VP based on how many reads with links were found for each position
weight_frame = hapfile[['pos1','weight']]
weight_frame['total']=weight_frame.groupby(['pos1'])['weight'].transform('sum')
weight_frame = weight_frame[['pos1','total']].drop_duplicates()
max_weight = weight_frame[weight_frame['total'] == weight_frame['total'].max()]

if viewpointpos == 0:
    viewpointpos = int(max_weight['pos1'].head(1))
    
#find the two possible NTs at the viewpoint pos
vars_opts = pgsnp[pgsnp['pos']==viewpointpos]['vars'].values
var1 = vars_opts[0].split('/')[0]
var2 = vars_opts[0].split('/')[1]
viewpointvars=var1,var2

#create hap1 and assign the VP values
hap1 = pd.DataFrame(columns=['pos', 'var', 'ID', 'weight', 'linkedweight', 'ambig', 'class', 'ratio', 'round', 'plexity', 'amplexity'])

hap1 = hap1.append({
    'pos' : viewpointpos,
    'var' : viewpointvars[0],
    'ID' : (str(viewpointpos) + ":" + viewpointvars[0]),
    'weight' : 10000,
    'linkedweight' : 10000,
    'ambig' : 0,
    'class' : "VP",
    'ratio' : 0,
    'round' : 0,
    'plexity' : 0,
    'amplexity' : 0
    }, ignore_index=True)

#create hap2 and assign the VP values
hap2 = pd.DataFrame(columns=['pos', 'var', 'ID', 'weight', 'linkedweight', 'ambig', 'class', 'ratio', 'round', 'plexity', 'amplexity'])
hap2 = hap2.append({
    'pos' : viewpointpos,
    'var' : viewpointvars[1],
    'ID' : (str(viewpointpos) + ":" + viewpointvars[1]),
    'weight' : 10000,
    'linkedweight' : 10000,
    'ambig' : 0,
    'class' : "VP",
    'ratio' : 0,
    'round' : 0,
    'plexity' : 0,
    'amplexity' : 0
    }, ignore_index=True)

#initiate seedhap for first round
seedhap1 = hap1
seedhap2 = hap2

#create IDlist with all positional values of the SNVs
IDlist = hapfile['pos1'].drop_duplicates()
IDlist = IDlist[IDlist != viewpointpos]
IDlistnew = IDlist

hapfile['ID1']=hapfile["pos1"].astype(str)+":"+hapfile['var1']
hapfile['ID2']=hapfile["pos2"].astype(str)+":"+hapfile['var2']

iterations = 25
iteration = iter(range(iterations))
average_ratio_1 = 0

#Loop to build on haplotype, each time using the haplotypes as seeds
for i in iteration:
    loopnr = i+1


    if IDlistnew.size > 0:
        badsnps = {"isolated": [], "single_allele": [], "ambigious": []}
        for snp in IDlistnew:
            #all lines from hapfile which contain the current SNP position were itterating over
            links = hapfile[hapfile['pos1']==snp]
            #all lines from that which have the other position in the current haplotype 1
            hap_links = links[links['pos2'].isin(seedhap1['pos'])]
            
            if (len(links['ID1'].unique())==2):
                #get all links that belong to one or the other ID for the linked variant
                var1 = links['ID1'].drop_duplicates().reset_index(drop=True)[0]
                var2 = links['ID1'].drop_duplicates().reset_index(drop=True)[1]
                var1_links = hap_links[hap_links['ID1'] == var1]
                var2_links = hap_links[hap_links['ID1'] == var2]

                #determine weight for possibilities of how vars can be linked to haplotypes
                #posibility 1: allele/hap 1/1 - 2/2
                var1_hap1 = var1_links[var1_links["ID2"].isin(seedhap1['ID'])]
                var2_hap2 = var2_links[var2_links["ID2"].isin(seedhap2['ID'])]
                p1_weight = var1_hap1['weight'].sum() + var2_hap2['weight'].sum()
                p1_plexity = len(var1_hap1) + len(var2_hap2)

                #posibility 1: allele/hap 1/2 - 2/1
                var1_hap2 = var1_links[var1_links["ID2"].isin(seedhap2['ID'])]
                var2_hap1 = var2_links[var2_links["ID2"].isin(seedhap1['ID'])]
                p2_weight = var1_hap2['weight'].sum() + var2_hap1['weight'].sum()
                p2_plexity = len(var1_hap2) + len(var2_hap1)

                if(p1_weight == 0 and p2_weight == 0):
                    badsnps["isolated"].append(snp)

                if(p1_weight == p2_weight):
                    badsnps["ambigious"].append(snp)
                
                #if var1 is in hap1 en vice versa
                if(p1_weight>p2_weight):
                    weight = p1_weight
                    ambig = p2_weight
                    plexity = p1_plexity
                    amplexity = p2_plexity
                
                    #check if it's a double of single link
                    if (var1_hap1['weight'].sum() > 0 and var2_hap2['weight'].sum() > 0):
                        snv_class_1 = "Double_link"
                        snv_class_2 = "Double_link"
                    elif (var1_hap1['weight'].sum() > 0 and var2_hap2['weight'].sum() == 0): 
                        snv_class_1 = "Single_link"
                        snv_class_2 = "Indirect"
                    elif (var2_hap2['weight'].sum() > 0 and var1_hap1['weight'].sum() == 0): 
                        snv_class_2 = "Single_link"
                        snv_class_1 = "Indirect"

                    #write line for haplotype 1
                    hap1 = hap1.append({
                                'pos' : snp,
                                'var' : var1[-1:],
                                'ID' : var1,
                                'weight' : weight,
                                'linkedweight' : weight,
                                'ambig' : ambig,
                                'class' : snv_class_1,
                                'ratio' : round(ambig/weight, 3),
                                'round' : i+1,
                                'plexity' : plexity,
                                'amplexity' : amplexity
                                }, ignore_index=True)

                    #write line for haplotype 2
                    hap2 = hap2.append({
                                'pos' : snp,
                                'var' : var2[-1:],
                                'ID' : var2,
                                'weight' : weight,
                                'linkedweight' : weight,
                                'ambig' : ambig,
                                'class' : snv_class_2,
                                'ratio' : round(ambig/weight, 3),
                                'round' : i+1,
                                'plexity' : plexity,
                                'amplexity' : amplexity
                                }, ignore_index=True)

                #if var1 is in hap2 en vice versa
                if(p2_weight>p1_weight):
                    weight = p2_weight
                    ambig = p1_weight
                    plexity = p2_plexity
                    amplexity = p1_plexity
                
                    #check if it's a double of single link
                    if (var1_hap2['weight'].sum() > 0 and var2_hap1['weight'].sum() > 0):
                        snv_class_1 = "Double_link"
                        snv_class_2 = "Double_link"
                    elif (var2_hap1['weight'].sum() > 0 and var1_hap2['weight'].sum() == 0): 
                        snv_class_1 = "Single_link"
                        snv_class_2 = "Indirect"
                    elif (var1_hap2['weight'].sum() > 0 and var2_hap1['weight'].sum() == 0): 
                        snv_class_2 = "Single_link"
                        snv_class_1 = "Indirect"

                    #write line for haplotype 1
                    hap1 = hap1.append({
                                'pos' : snp,
                                'var' : var2[-1:],
                                'ID' : var2,
                                'weight' : weight,
                                'linkedweight' : weight,
                                'ambig' : ambig,
                                'class' : snv_class_1,
                                'ratio' : round(ambig/weight,3),
                                'round' : loopnr,
                                'plexity' : plexity,
                                'amplexity' : amplexity
                                }, ignore_index=True)

                    #write line for haplotype 2
                    hap2 = hap2.append({
                                'pos' : snp,
                                'var' : var1[-1:],
                                'ID' : var1,
                                'weight' : weight,
                                'linkedweight' : weight,
                                'ambig' : ambig,
                                'class' : snv_class_2,
                                'ratio' : round(ambig/weight,3),
                                'round' : loopnr,
                                'plexity' : plexity,
                                'amplexity' : amplexity
                                }, ignore_index=True)

            else:
                badsnps["single_allele"].append(snp)
    
    #Subselect and recalculate with decreasing stringency per iteration
    if loopnr < 5:
        quantileweight = 30-(5*loopnr)
        hap1 = hap1[(hap1['ratio']<0.1)|(hap1['class']=="VP")]
        hap2 = hap2[(hap2['ratio']<0.1)|(hap2['class']=="VP")]

        hap1, hap2 = recalculate_all(hap1, hap2)

        hap1 = hap1[(hap1['linkedweight']>quantileweight)&(hap1['class']!="Single_link")&(hap1['class']!="Indirect")|(hap1['class']=="VP")]
        hap2 = hap2[(hap2['linkedweight']>quantileweight)&(hap2['class']!="Single_link")&(hap2['class']!="Indirect")|(hap2['class']=="VP")]

    elif loopnr < (iterations-4): 
        hap1 = hap1[hap1['ratio']<0.3]
        hap2 = hap2[hap2['ratio']<0.3]
        
        hap1, hap2 = recalculate_all(hap1, hap2)

        hap1 = hap1[(hap1['class']!="Single_link")&(hap1['class']!="Indirect")]
        hap2 = hap2[(hap2['class']!="Single_link")&(hap2['class']!="Indirect")]
        
    else:
        hap1, hap2 = recalculate_all(hap1, hap2)

    #Update the seedhap and IDlist for the next round
    seedhap1 = hap1
    seedhap2 = hap2

    IDlistnew = IDlist[~IDlist.isin(seedhap1['pos'])]
    IDlistnew = IDlistnew[~IDlistnew.isin(seedhap2['pos'])]

#Re-evaluate bad snps
final_badsnps = {"Single_linked": [], "Single_ambiguous": [], "Single_isolated": [], "Isolated": []}

hap1_b = pd.DataFrame(columns=['pos', 'var', 'ID', 'weight', 'linkedweight', 'ambig', 'class', 'ratio', 'round', 'plexity', 'amplexity'])
hap2_b = pd.DataFrame(columns=['pos', 'var', 'ID', 'weight', 'linkedweight', 'ambig', 'class', 'ratio', 'round', 'plexity', 'amplexity'])

for snp in badsnps["single_allele"]:
    #get all links that with the snp that link to hap1
    links = hapfile[hapfile['pos1'] == snp]
    hap_links = links[links['pos2'].isin(hap1['pos'])]
    
    #Split all these links into whether they would fit in hap1 or hap2
    hap1_links = hap_links[hap_links['ID2'].isin(hap1['ID'])]
    hap2_links = hap_links[hap_links['ID2'].isin(hap2['ID'])]

    #get the NT for each to use when adding to the dataframe
    hap1_var = hap1_links['var1'].drop_duplicates().reset_index(drop=True)
    hap2_var = hap2_links['var1'].drop_duplicates().reset_index(drop=True)

    #get weight and plexity for each    
    weight_1 = hap1_links['weight'].sum()
    weight_2 = hap2_links['weight'].sum()        
    plexity_1 = len(hap1_links)
    plexity_2 = len(hap2_links)

    #determine whether the badsnp fits into one, both or none of the haplotypes and process accordingly
    if weight_1 > 0 and weight_2 == 0:
        final_badsnps["Single_linked"].append(snp)
        hap1 = hap1.append({
                    'pos' : snp,
                    'var' : hap1_var[0],
                    'ID' : str(snp) + ":" + hap1_var[0],
                    'weight' : weight_1,
                    'linkedweight' : weight_1,
                    'ambig' : weight_2,
                    'class' : "Single_link",
                    'ratio' : round(weight_2/weight_1,3),
                    'round' : 26,
                    'plexity' : plexity_1,
                    'amplexity' : plexity_2
                    }, ignore_index=True)

    if weight_2 > 0 and weight_1 == 0:
        final_badsnps["Single_linked"].append(snp)
        hap2 = hap2.append({
                    'pos' : snp,
                    'var' : hap2_var[0],
                    'ID' : str(snp) + ":" + hap2_var[0],
                    'weight' : weight_2,
                    'linkedweight' : weight_2,
                    'ambig' : weight_1,
                    'class' : "Single_link",
                    'ratio' : round(weight_1/weight_2,3),
                    'round' : 26,
                    'plexity' : plexity_2,
                    'amplexity' : plexity_1
                    }, ignore_index=True)
                    
    if weight_1 > 0 and weight_2 > 0:
        final_badsnps["Single_ambiguous"].append(snp)

    if weight_1 == 0 and weight_2 == 0:
        final_badsnps["Single_isolated"].append(snp)

#Filter bad snps for acceptable ones after recalculation
hap1_c = pd.concat([hap1, hap1_b], axis=0)    
hap2_c = pd.concat([hap2, hap2_b], axis=0)    

hap1_c, hap2_c = recalculate_all(hap1=hap1_c, hap2=hap2_c)

hap1_c = hap1_c[hap1_c['ratio']<0.5]
hap2_c = hap2_c[hap2_c['ratio']<0.5]

hap1_b = hap1_b[hap1_b['ID'].isin(hap1_c['ID'])]
hap2_b = hap2_b[hap2_b['ID'].isin(hap2_c['ID'])]

hap1 = pd.concat([hap1, hap1_b], axis=0)
hap2 = pd.concat([hap2, hap2_b], axis=0)

#Final recalculation and filtering for acceptable ratio
hap1, hap2 = recalculate_all(hap1, hap2)

#get all single_linked_bad positions
hap1_sl_pos = hap1[(hap1['class']=="Single_link")&(~hap1['pos'].isin(hap2['pos']))]['pos']
hap2_sl_pos = hap2[(hap2['class']=="Single_link")&(~hap2['pos'].isin(hap1['pos']))]['pos']

for pos in hap1_sl_pos:
    #get both possible vars at position from pgsnp file
    vars_opts = pgsnp[pgsnp['pos']==pos]['vars'].values
    var1 = vars_opts[0].split('/')[0]
    var2 = vars_opts[0].split('/')[1]
    var_hap1 = hap1[hap1['pos']==pos]['var'].values[0]
    
    #find which is var for each
    if var_hap1 == var1:
        var_hap2 = var2
    if var_hap1 == var2:
        var_hap2 = var1

    #copy hap2 line and adjust to make hap1 indirect linked snp
    hap2_line = hap1[hap1['pos']==pos]
    hap2_line['var'] = var_hap2
    hap2_line['ID'] = str(pos) + ":" + var_hap2
    hap2_line['class'] = "Indirect"

    #add to hap1
    hap2 = pd.concat([hap2, hap2_line])

for pos in hap2_sl_pos:
    #get both possible vars at position from pgsnp file
    vars_opts = pgsnp[pgsnp['pos']==pos]['vars'].values
    var1 = vars_opts[0].split('/')[0]
    var2 = vars_opts[0].split('/')[1]
    var_hap2 = hap2[hap2['pos']==pos]['var'].values[0]
    
    #find which is var for each
    if var_hap2 == var1:
        var_hap1 = var2
    if var_hap2 == var2:
        var_hap1 = var1

    #copy hap2 line and adjust to make hap1 indirect linked snp
    hap1_line = hap2[hap2['pos']==pos]
    hap1_line['var'] = var_hap1
    hap1_line['ID'] = str(pos) + ":" + var_hap1
    hap1_line['class'] = "Indirect"

    #add to hap1
    hap1 = pd.concat([hap1, hap1_line])

hap1 = hap1[hap1['ratio']<0.5]
hap2 = hap2[hap2['ratio']<0.5]

hap1.drop_duplicates(subset=['pos'], inplace=True, keep='last')
hap2.drop_duplicates(subset=['pos'], inplace=True, keep='last')

hap1['hap'] = 1
hap2['hap'] = 2

#Write output files
hap1_out = hap1[['pos', 'var']]
hap2_out = hap2[['pos', 'var']]
haps_out = pd.concat([hap1, hap2], axis=0)

hap1_out.to_csv(hap1_out_file, sep=":", index=False, header=False)
hap2_out.to_csv(hap2_out_file, sep=":", index=False, header=False)
haps_out.to_csv(haps_out_file, sep="\t", index=False)

#make network files

hapfile['node1added'] = "X"
hapfile['node2added'] = "X"

hap1_links_nw = ""
hap2_links_nw = ""

for index, row in hap1.iterrows():
    connectables = hap1[hap1['round']<row['round']]['ID']
    links_nw = hapfile[hapfile['ID2']==row['ID']]
    links_nw = links_nw[links_nw['ID1'].isin(connectables)]

    if len(links_nw)>0:
        links_nw['used'] = row['round']
        links_nw['node2added'] = row['round']
        for index, row in links_nw.iterrows():
            
            linkID = links_nw.loc[[index],['ID1']].values[0][0]
            if len(hap1[hap1['ID']==linkID]['round']) > 0:
                node1added = hap1[hap1['ID']==linkID]['round'].values[0]
                links_nw.loc[[index],['node1added']] = node1added

        if len(hap1_links_nw) > 0:
            hap1_links_nw = pd.concat([hap1_links_nw, links_nw], axis=0)
        else:
            hap1_links_nw = links_nw  
    else:
        hap1.loc[[index],['class']] = "indirect"

for index, row in hap2.iterrows():
    connectables = hap2[hap2['round']<row['round']]['ID']
    links_nw = hapfile[hapfile['ID1']==row['ID']]
    links_nw = links_nw[links_nw['ID2'].isin(connectables)]

    if len(links_nw)>0:
        links_nw['used'] = row['round']
        links_nw['node2added'] = row['round']
        for index, row in links_nw.iterrows():
            linkID = links_nw.loc[[index],['ID1']].values[0][0]
        
            if len(hap2[hap2['ID']==linkID]['round']) > 0:
                node1added = hap2[hap2['ID']==linkID]['round'].values[0]
                links_nw.loc[[index],['node1added']] = node1added

        if len(hap2_links_nw) > 0:
            hap2_links_nw = pd.concat([hap2_links_nw, links_nw], axis=0)
        else:
            hap2_links_nw = links_nw  
    else:
        hap2.loc[[index],['class']] = "indirect"

hap1_links_nw['hap'] = 1
hap1_links_nw['hap.1'] = 1
hap1_links_nw['interaction'] = "pp"

hap2_links_nw['hap'] = 2
hap2_links_nw['hap.1'] = 2
hap2_links_nw['interaction'] = "pp"

all_links = pd.concat([hap1_links_nw, hap2_links_nw], axis = 0)

both = pd.concat([hap1, hap2])


netwtable = all_links.iloc[:, [5,6,12,9,7,8,10,11,1,3,4]]

netwtable['weight1'] = "X"
netwtable['weight2'] = "X"
netwtable['ambig1'] = "X"
netwtable['ambig2'] = "X"
netwtable['plexity1'] = "X"
netwtable['plexity2'] = "X"
netwtable['amplexity1'] = "X"
netwtable['amplexity2'] = "X"

for index, row in netwtable.iterrows():
    hapdata1 = both[both['ID'] == row['ID1']].reset_index(drop=True)
    if len(hapdata1) == 1:
        netwtable.loc[[index],['weight1']] = hapdata1['weight'].values[0]
        netwtable.loc[[index],['ambig1']] = hapdata1['ambig'].values[0]
        netwtable.loc[[index],['plexity1']] = hapdata1['plexity'].values[0]
        netwtable.loc[[index],['amplexity1']] = hapdata1['amplexity'].values[0]

    
    hapdata2 = both[both['ID'] == row['ID2']].reset_index(drop=True)
    if len(hapdata2) == 1:
        netwtable.loc[[index],['weight2']] = hapdata2['weight'].values[0]
        netwtable.loc[[index],['ambig2']] = hapdata2['ambig'].values[0]
        netwtable.loc[[index],['plexity2']] = hapdata2['plexity'].values[0]
        netwtable.loc[[index],['amplexity2']] = hapdata2['amplexity'].values[0]


netwtable['edgeweight'] = netwtable['weight']
netwtable.drop(columns = ['weight'], inplace = True)

netwtable.columns = ["ID1", "ID2", "interaction","interactionused", "added", "added", "hap", "hap",
                    "pos", "pos","weight", "weight", "ambig", "ambig", "plexity", "plexity", "amplexity", "amplexity", 
                    "interactionweight"]

netwtable.to_csv(netwtable_file, sep="\t", index=False)

class VCF_file:
    
    def load(self,vcf,chromosome):
        if vcf.endswith('.gz'):
            file = gzip.open(vcf,'rt')
        else:
            file = open(vcf)
        header_count = 0
        self.header = []
        for line in file.readlines():
            if line.startswith("##"):
                header_count += 1
                self.header.append(line)
            else:
                break
        file.close()
        variants = pd.read_csv(vcf,sep='\t',header=header_count)
        variants.index = variants['POS']
        if chromosome != None:
            self.data = variants.loc[variants['#CHROM'] == chromosome]
        else:
            self.data = variants
        
    def __init__(self,vcf,chromosome=None):
        self.load(vcf,chromosome)
        self.haplotype = pd.DataFrame(columns=['hap1','hap2'])
        
    def phase(self,haps_out):
        haps_sorted = haps_out.reset_index()
        for row in  haps_sorted.index:
            if haps_sorted.loc[row]['pos'] in self.data.index:
                if self.data.loc[haps_sorted.loc[row]['pos'],'REF'] == haps_sorted.loc[row]['var']:
                    self.haplotype.loc[haps_sorted.loc[row]['pos'],f"hap{haps_sorted.loc[row]['hap']}"] = 0
                elif self.data.loc[haps_sorted.loc[row]['pos'],'ALT'] == haps_sorted.loc[row]['var']:
                    self.haplotype.loc[haps_sorted.loc[row]['pos'],f"hap{haps_sorted.loc[row]['hap']}"] = 1
        for pos in self.haplotype.index:
            sample = self.data.columns[9]
            sample_split = self.data.loc[pos,sample].split(':')
            if len([x for x in self.haplotype.loc[pos] if type(x) != int]) == 0:
                sample_split[0] = '|'.join([str(x) for x in self.haplotype.loc[pos]])
                self.data.loc[pos,sample] = ':'.join(sample_split)
            
    def write(self,out):
        with open(out,'w') as fh:
            for line in self.header:
                fh.write(line)
            fh.write('\t'.join([str(x) for x in self.data.columns])+'\n')
            for row in self.data.index:
                fh.write('\t'.join([str(x) for x in self.data.loc[row]])+'\n')

variants_vcf = VCF_file(vcf_path,chromosome.split(':')[0])
variants_vcf.phase(haps_out)
variants_vcf.write(vcf_out_file)
