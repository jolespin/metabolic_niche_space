import pandas as pd
import numpy as np
# import scipy.stats
from scipy.stats import hypergeom
# from re import search
import csv
import os.path

import sys
# sys.path.append('../GOCAM_Project/dev')
import os

import rpy2
from rpy2.robjects.packages import importr
BiasedUrn = importr('BiasedUrn')


def csv2dict(file, sep = ','):
    d={}
    with open(file, "r") as f:
        csv_reader = csv.reader(f, delimiter = sep)
        for row in csv_reader:
            key = row[0]            
            val = list(row[1:])
            d[key] = val
    return d

def map_dict_vals(mapping_dict,dict2bemapped):
    dictmapped = {}
    for k, vals in dict2bemapped.items():
        mapped_vals = set(pd.Series(list(vals)).map(mapping_dict))
        # dictmapped[k] = inputs
        dictmapped[k] = mapped_vals

    return dictmapped

def dict2csv(d,file, sep = ','):
    array = []
    for key, val in d.items():
        temp = [key]
        temp.extend(list(val))
        #trun = [i.split('/')[-1] for i in temp]
        #trun[0] = trun[0].split('-')[-1]
        array.append(temp) #change to trun

    with open(file, "w") as f:
        writer = csv.writer(f, delimiter = sep)
        writer.writerows(array)
        
def update_dict(d,df,key,val):
    for index, row in df.iterrows():
        if (row[key] in d) == False:
            d[row[key]]={row[val]}
        else:
            prev = d.get(row[key])
            prev.add(row[val])
            d[row[key]] = prev
    return d


def reverse_dict(original_dict):
    reversed_dict = {}
    
    for key, values in original_dict.items():
        for value in values:
            if value not in reversed_dict:
                reversed_dict[value] = []
            reversed_dict[value].append(key)
    
    return reversed_dict

def kegg_make_EC_reaction():
    # Read the input file
    with open('../data/kegg/wol2/reaction_enzyme.txt', 'r') as file:
        lines = file.readlines()

    # Process the lines to create the desired output format
    output_lines = []
    for line in lines:
        identifier, ec_numbers = line.strip().split('\t')
        ec_numbers = ec_numbers.strip("[]").replace("'", "").split(', ')
        output_line = identifier + '\t' + '\t'.join(ec_numbers)
        output_lines.append(output_line)

    # Write the output to a TSV file
    with open('../data/kegg/wol2/reaction-to-enzyme.map', 'w') as file:
        for output_line in output_lines:
            file.write(output_line + '\n')

    print("reaction-to-enzyme.map created from reaction_enzyme.txt")
    
    re_ec = csv2dict('../data/kegg/wol2/reaction-to-enzyme.map', sep = '\t')
    ec_re = reverse_dict(re_ec)
    dict2csv(ec_re, '../data/kegg/wol2/enzyme-to-reaction.map', sep = '\t')
    print("enzyme-to-reaction.map created")
    
def kegg_make_pathway_reaction(backend_type, enrich_against):
    
    be_pa = csv2dict(f'../data/kegg/wol2/{backend_type}-to-{enrich_against}.map', sep = '\t')
    pa_be = reverse_dict(be_pa)
    dict2csv(pa_be, f'../data/kegg/wol2/{enrich_against}-to-{backend_type}.map', sep = '\t')
    print(f"{enrich_against}-to-{backend_type}.map created")
    
def get_sizes (data): #data= dataframe with gocam IDs and gene identifiers as columns
    """get number of entities in each gocam"""
    return data['gocam'].value_counts()
    
def get_sets (gene_list):
    """map list of genes to all sets that contain members of that list"""
    sets = []
    not_in_a_set = []
    members2setID = csv2dict('../data/members2setID.csv')
    setID2members_input = {}
    for g in gene_list:
        s = members2setID.get(g)
        if s != None:
            sets = sets +s
            for i in s:
                if (i in setID2members_input) == False:
                    setID2members_input[i]={g}
                else:
                    prev = setID2members_input.get(i)
                    prev.add(g)
                    setID2members_input[i] = prev
        else:
            not_in_a_set.append(g)
    return not_in_a_set, list(set(sets)),setID2members_input #remove duplicates

def filter_gene_list(gene_list, Dict):
    """remove members of gene_list that are not in Dict.
    use function to filter a user's input list of genes based on those that appear at least 
    once in the gocam model database"""
    filtered_gene_list = []
    filtered_out = []
    for gene in gene_list:
        if gene in Dict:
            filtered_gene_list.append(gene)
        else:
            filtered_out.append(gene)
    return filtered_out, filtered_gene_list

def count_genes(gene_list, Dict):
    """ count number of genes in user's gene_list that are in each gocam"""
    gocam_counts = {} #key=gocam, value=list of genes in gocam that are also in the user's list
    for g in gene_list:
        gocams = Dict.get(g)
        for gocam in gocams:
            if (gocam in gocam_counts) == False:
                gocam_counts[gocam]=[g]
            else:
                prev = gocam_counts.get(gocam)
                prev.append(g)
                gocam_counts[gocam] = prev
    return gocam_counts

def remove_non_pathways_kegg(counts, enrich_against):
    """kegg has a number of pathways that aren't meaningful pathways, such as ko01100 'Metabolic Pathways'. These are also very large and slow down the code"""
    if enrich_against == 'pathways':
        for p in ['ko01100','ko01110','ko01120', 'ko01200']:
            counts.pop(p, None)
    elif enrich_against == 'modules':
        pass
    return counts

def construct_dicts(custom):
    pass

def convert_IDs(genes, input_type, backend_type = 'reaction', enrich_against = None):
    """converts input type to uniprot IDs for gocams"""
    # if input_type == 'ko' or input_type == 'enzyme': #kegg
    file = f'../data/kegg/wol2/{input_type}-to-{backend_type}.map'
    if input_type == 'enzyme' and os.path.isfile(file) == False:
        kegg_make_EC_reaction()
    d = csv2dict(file, sep = '\t')
    genes[backend_type]=genes['g'].apply(lambda x: d.get(x))
    t = genes[genes[backend_type].isna()].copy()        
    not_converted = t[t[backend_type].isna()]
    genes = genes.dropna()
    genes = pd.concat([genes,t[t[backend_type].notna()]])
    temp = genes.explode(backend_type)
    temp.drop_duplicates(subset=backend_type,inplace = True)
    
    #remove reactions that are not annotated to a pathway
    d2 = csv2dict(f'../data/kegg/wol2/{backend_type}-to-{enrich_against}.map', sep = '\t')
    mask = temp[backend_type].apply(lambda x: True if d2.get(x) != None else False)
    temp = temp[mask]
    
    reactions2input = pd.Series(temp.g.values, index=temp[backend_type]).to_dict()
    reaction_list = list(temp[backend_type].values)
    return reaction_list, reactions2input, not_converted

    

#Benjamini Hochberg correction
def correct_pval_and_format(enriched_gocams, background_num_gocams,FDR,kegg, enrich_against = None):
    """performs Benjamini Hochberg correction to control the false discovery rate and formats output for display"""
    df = pd.DataFrame(enriched_gocams, columns =['url', 'pval (uncorrected)', '# entities in list','#entities in model','shared entities in gocam'])
    df.sort_values('pval (uncorrected)',inplace=True)
    df.reset_index(drop=True, inplace=True)
    df['FDR_val'] = (df.index+1)*FDR/background_num_gocams
    df['Less_than'] = (df['pval (uncorrected)'] < df['FDR_val'])
    index = df.Less_than.where(df.Less_than==True).last_valid_index()
    df_significant = df
    
    df_significant = df.loc[0:index].copy()
    if index == None:
        df_significant = pd.DataFrame(columns =['url', 'pval (uncorrected)', '# entities in list','#entities in model','shared entities in gocam'])
    df_display = df_significant[['url','pval (uncorrected)', '# entities in list', '#entities in model','shared entities in gocam']].copy()
    #modelID2title = pd.read_csv('../data/modelID2title_mouse.csv')
    
    temp = pd.read_csv('../data/modelID2title_mouse.csv',header = 0,names=['pathway','title'])
    if kegg:
        temp = pd.read_csv(f'../data/kegg/wol2/{enrich_against}_name.txt',header = None, sep = '\t', names=['pathway','title'])
    modelID2title = pd.Series(temp.title.values,index=temp.pathway).to_dict()
    
    df_display['title'] = df_display['url'].map(modelID2title)
    cols = df_display.columns.to_list()
    cols[0]='title'
    cols[-1]='url'
    df_display = df_display[cols]
    return df_display




def get_M_wM(setID2members, ID2pathway):
    """ returns M, the number of entities in the background, and w_M, the mean size of entities in the background"""
    
    l = []
    for s,m in setID2members.items():
        l.append(len(m))
    l = np.array(l)
    l = np.sort(l)
    num_empty_sets = np.sum(l==0)
    
    l = l[l!=0]
    mean = np.mean(l)#l[4:-4]) 1% trimmed mean?
    num_sets = len(l)
    bg = len(ID2pathway)
    M = bg-num_empty_sets
    
    w_M = np.round(((M-num_sets)+num_sets*mean)/M,decimals=2)
    return M, w_M

def make_initial_vectors(gocam2ID,setID2members, gc, M,w_M, kegg = False):
    """initializes counts vector (m) and weights vector (w), where each entity gets its own element in the arrays
- values in m only take on 0 (if there is no solo proteins) or 1
- values in w correspond to the weight of each element in m (weighted by the # genes in a set or 1 for solo proteins)"""
    w_gc = [1] #initialize with 1 as the weight of single proteins (irrespective of whether there are any)
    m_gc = [0] #initialize with 0 single proteins
    num_protein = 0
    for i in gocam2ID.get(gc):
        if "sset:" in i or kegg: #all reactions treated as sets in kegg
            try:
                w_i = len(setID2members[i])
            except KeyError:
                print(f'{i} not found in dictionary')
                continue
            if kegg and w_i == 1:
                num_protein +=1
                continue
                
            w_gc.append(w_i)
            m_gc.append(1)
        else:
            num_protein+=1
    m_gc[0] = num_protein
    m_gc.append(M-np.sum(m_gc)) #entities not in the gocam (roughly)
    w_gc.append(w_M) #weight for entities not in the gocam (all weighted as w_M (the mean))
    return w_gc, m_gc


def make_new_vectors(w_gc,m_gc,M,w_M):
    """compress the m and w vectors by grouping elements according to their weights
- w is the ordered set of unique weights for entities of the gocam + the background bin
- m[i] is the number of entities in the pathway with the weight specified in w[i] + the background bin"""
    w_temp = w_gc[:-1]
    if w_temp[0] != 1:
        print('Possible bug: w_temp[0] != 1',w_temp)
        
    w_new, m_temp = np.unique(w_temp, return_counts=True)
    m_temp[0]=m_gc[0] #w_gc and m_gc have weight 1 as w_gc[0] and the number of single proteins as m_gc[0]
    m_new = np.append(m_temp,np.array([M-np.sum(m_temp)]))
    w_new = np.append(np.unique(w_temp),np.array([w_M]))
    return w_new, m_new




def ncHGT_sf(XT,m,N,w):
    """survival function, sums PMF for all possibilities where K >= k by calling BiasedUrn"""
    #l = len(XT)/len(m)
    if len(XT) == 0:
        print('len(XT) = 0')
        return -1
    pval = 0
    #np.seterr(under='warn')
    ### This could be optimized by setting a threshold and stopping the for loop when the sum exceeds some threshold###
    for i in range(len(XT)):
        x = rpy2.robjects.IntVector(XT[i])
        pval = pval + BiasedUrn.dMFNCHypergeo(x,m,N,w, precision = 1e-10)[0]
    return pval


def enumerate_possibilities(m_new,i,prev_array):
    """enumerate all possible counts vectors"""
    first = True
    for j in range(m_new[i]+1):
        xt = prev_array.copy()
        xt[0][i] = j
        
        #recursion
        if (i < len(m_new)-1):
            xt = enumerate_possibilities(m_new, i+1, xt) #will return matrix (array of arrays)
            
        #combining results into matrix
        if not first:
            XT = np.concatenate([XT,xt], axis = 0)
        else:
            XT = xt
            first = False
    return XT


def do_ncHGT(k,gc,M,N, input_type = '', kegg = False, backend_type = '', enrich_against = None):
    setID2members = '../data/setID2members.csv'
    gocam2ID = '../data/gocam2ID_mouse.csv'
    ID2gocam = '../data/ID2gocam_mouse.csv'
    sep = ','
    if kegg:
        setID2members = f'../data/kegg/wol2/reaction-to-{input_type}.map'
        gocam2ID = f'../data/kegg/wol2/{enrich_against}-to-reaction.map'
        ID2gocam = f'../data/kegg/wol2/reaction-to-{enrich_against}.map'
        sep = '\t'
        if os.path.isfile(f'../data/kegg/wol2/{enrich_against}-to-{backend_type}.map') == False:
            kegg_make_pathway_reaction(backend_type, enrich_against)
             
    setID2members = csv2dict(setID2members, sep = sep)
    gocam2ID = csv2dict(gocam2ID, sep = sep)
    ID2gocam = csv2dict(ID2gocam, sep = sep)
    
    M, w_M = get_M_wM(setID2members, ID2gocam)
    
    #make weight (w) and bin size (m) vectors where each entity in the gocam gets its own entry
    w_in, m_in = make_initial_vectors(gocam2ID, setID2members, gc, M,w_M, kegg = kegg)

    #update m and w vectors by grouping sets of the same size
    w_new , m_new= make_new_vectors(w_in,m_in,M,w_M)

    
    #make XT matrix, an enumeration of all possible arangements of balls in bins based on m_new and w_new
    m_gc = m_new[:-1] #don't pass the background bin to XT
    XT = enumerate_possibilities(m_gc,0,np.zeros(shape=(1,len(m_gc))))
    
    complement = False
    if kegg and k < 5: #initial testing was extremely slow with kegg
        complement = True
        #filter XT to only include the region of the sample space < k (complement of what we want to sum probabilities over)
        mask1 = (np.sum(XT, axis=1) < k)
        XT = XT[mask1]
    else:
        #filter XT to only include the region of the sample space >= k (which is what we want to sum probabilities over)
        mask1 = (np.sum(XT, axis=1) >= k)
        XT = XT[mask1]

    #filter XT to ensure that more than N entities are not picked
    mask2 = (np.sum(XT, axis=1) <= N)
    XT = XT[mask2]

    #add the remaining entities to the m+1th bin (non gocam bin)
    x_mp1_vec = N- np.sum(XT, axis = 1) #number of balls to be drawn from the last bin (the non-gocam background)
    XT = np.concatenate((XT,x_mp1_vec.reshape(len(x_mp1_vec),1)), axis = 1)
    
    m = rpy2.robjects.IntVector(m_new)
    w = rpy2.robjects.FloatVector(w_new)
    pval = ncHGT_sf(XT,m,N,w)
    if complement:
        pval = 1 - pval
    return pval


#ncHGT is either False (indicating that regular HGT should be done) or a positive integer denoting N for ncHGT
def hgt(counts, gocam_sizes, FDR, gene_list_size, background_gene_list_size, ncHGT = False, 
        kegg = False, input_type = None, backend_type = None, enrich_against = None):
    """ performs either the hypergeometric test or our introduced test using Fisher's noncentral hypergeometric dist.
    Whether our unweighted set enrichment or the standard HGT is performed is determined upstream based on what
    Dict of gocams->entities and filtered gene_list are passed into count_genes().
    ncHGT is either False (for set or standard methods) or corresponds to N """
    results = []
    iterator = tqdm.tqdm(counts.items())
    for gocam, gene_list in iterator:
        count = len(gene_list) 
        gocam_size = gocam_sizes[gocam]
        pvalue = None
        if ncHGT:
            if count <=1: #avoid unnecessary calls to BiasedUrn due to computation time
                pvalue = 1
            else:
                #changed from count -1 to count on 7/15/24. This error would not invalidate any results in the paper. It would make them more significant, as 
                #p-values are the probability of obtaining something "as extreme as" k. Survival functions give P(K > k) which is why we use hypergeom.sf(count-1) 
                #in the next code block. However, do_ncHGT computes P(K >= k). This error led to reporting values as slightly less significant by also adding 
                #P(K = k - 1) to the sum
                pvalue = noncentralHGT.do_ncHGT(count,gocam,background_gene_list_size,ncHGT, 
                                                kegg = kegg, input_type = input_type, backend_type = backend_type, enrich_against = enrich_against)
        else: #set or standard methods
            pvalue = hypergeom.sf(count-1, background_gene_list_size,  gocam_size, gene_list_size) 
        if pvalue < 1: #FDR:
            r = (gocam, pvalue, count, gocam_size, gene_list )
            results.append(r)
    return results


#Dict can only contain 1 instance of each gene per gocam (no duplicates)
def enrich(gene_list, uni_list,uniprot2input,pathway_sizes, Dict, ncHGT=False,FDR=.05, 
           kegg = False, input_type = None, backend_type = None, enrich_against = None):
    """uni_list is the list of uniprot IDs, because the backend dictionary, Dict, is gocam_id-> list(uniprot id's).
    uniprot2input is a dictionary keeping track of which of the user's inputs mapped to which uniprot id's so results can be 
    displayed in the user's inputted format, as the mapping is not always 1:1."""
    background_gene_list_size = len(Dict)
    if kegg == False:
        if ncHGT: 
        #we consider the background size to be equal to the total # of genes 
        #(the sum of the weights of all entities would double count genes that occur in multiple sets
        #... is this the right thing to do though?
            background_gene_list_size = len(csv2dict('../data/ID2gocam_mouse_ff.csv'))

        not_in_a_set, sets, setID2members_input_uni = get_sets(uni_list)

        setID2members_input = map_dict_vals(uniprot2input, setID2members_input_uni)

        filtered_out1, set_list_filtered = filter_gene_list(sets,Dict)
        filtered_out2, gene_list_filtered = filter_gene_list(uni_list, Dict) #need to clean gene_list to only include genes in the gocam


        filtered_list = gene_list_filtered + set_list_filtered
        gene_list_size = len(filtered_list)

        flist2input = {**uniprot2input, **setID2members_input}
        filtered_list_as_genes = set(pd.Series(list(filtered_list)).map(flist2input).explode())
        filtered_out_genes = set(gene_list) - filtered_list_as_genes

        counts = count_genes(filtered_list, Dict)

        N_ncHGT = False
        if ncHGT == True:
            N_ncHGT = len(gene_list)-len(filtered_out_genes)
            if N_ncHGT <= 0:
                return "error no genes found in gocams"

        enriched_gocams = hgt(counts, pathway_sizes, FDR, gene_list_size, background_gene_list_size, ncHGT=N_ncHGT)
        background_num_gocams = len(pathway_sizes)
        df_display = correct_pval_and_format(enriched_gocams, background_num_gocams,FDR, kegg)
        return filtered_out_genes, filtered_list, setID2members_input_uni, setID2members_input, df_display
    
    else:
        if ncHGT: 
            background_gene_list_size = len(csv2dict(f'../data/kegg/wol2/{input_type}-to-reaction.map', sep = '\t')) #should adjust this when filtering is applied
        gene_list_size = len(uni_list) #should be named backend_list. number of reactions
        counts = count_genes(uni_list, Dict) #number of reactions per pathway
        counts = remove_non_pathways_kegg(counts, enrich_against)
        N_ncHGT = False
        if ncHGT == True:
            N_ncHGT = len(gene_list) #number of inputted genes #CHECK
        enriched_pathways = hgt(counts, pathway_sizes, FDR, gene_list_size, background_gene_list_size, ncHGT=N_ncHGT, 
                                kegg = kegg, input_type = input_type, backend_type = backend_type, enrich_against = enrich_against)
        background_num_pathways = len(pathway_sizes)
        df_display = correct_pval_and_format(enriched_pathways, background_num_pathways,FDR, kegg, enrich_against = enrich_against)
        
        return ['no filtering with kegg'],uni_list,{'not done with kegg':'not done with kegg'}, uniprot2input, df_display
        

def enrich_wrapper(filename, input_type, method = 'set', return_all = False, FDR=.05,fpath= '../test_data', display_gene_symbol = True, display_input = False, 
                   kegg = False, backend_type = 'reaction', enrich_against = 'module', custom = {}):
    """ wrapper to perform enrichment given a filename, gene ID type, enrichment method, and false discovery rate.
    other parameters:
    
    return_all: 
        if false, only returns the dataframe displaying results. 
        if true: returns (gene_list, filtered_out_genes, filtered_list, setID2members_input_uni, setID2members_input, df_display)
        return_all = True is not just for debugging. User may want to know which of their input genes were filtered out as well as how
        the IDs were mapped, as uniprot IDs can sometimes map to more than one HGNC gene symbol
    display_gene_symbol: if true, display HGNC symbols on output regardless of input ID type"""
    if (backend_type != 'reaction' or enrich_against != 'module') and kegg == False:
        raise ValueError("backend_type, enrich_against argument(s) should only be used with kegg")
    #set method files
    #gcs = '../data/gocam_sizes_mouse.csv'
    id2g = '../data/ID2gocam_mouse.csv'
    sep = ','
    #standard method files
    if method == 'standard':
        gcs = gcs[:-4] + '_ff.csv'   #'../data/gocam_sizes_mouse_ff.csv'
        id2g = id2g[:-4] + '_ff.csv' #'../data/ID2gocam_mouse_ff.csv'
        
    if kegg:
        #gcs = '../data/kegg/pathway_sizes_kegg.csv'
        print('kegg functionality: map kegg orthologs or ECs to reactions and enrich reactions against pathways or modules. Reactions are treated as "sets." Custom addition or removal of entities not implemented yet. Weighted enrichment against pathways is slow due to large sizes; enrich_against for this is set to "module" as default. Certain non-meaningful pathways such as k01100 "Metabolic Pathways" are removed from enrichment.')
        sep = '\t'
        id2g = f'../data/kegg/wol2/{backend_type}-to-{enrich_against}.map'
        if method == 'standard':
            raise ValueError('standard mapping not implemented yet for kegg. use set or ncHGT')
        
        
    if custom != {}:
        id2g = construct_dicts(custom)
        
    gene_list = pd.read_csv(os.path.join(fpath,filename),header=None,names = ['g'])
    
    #normally not needed, but I found a bug where HSPA1A and HSPA1B are listed as synonyms, both in Simplemine and official sources like the Alliance
    gene_list.drop_duplicates(inplace = True) 
    
    gene_list_converted = []
    backend2input = {}
    not_converted = []
    
    #conversion to uniprot IDs not needed for a list of uniprot IDs
    #kegg does not use conversion. input must be ko or enzyme, and they are converted to reactions as an analog for sets
    if input_type in ['uniprot']:
        gene_list_converted = gene_list.g
        backend2input = pd.Series(gene_list_converted.values,index=gene_list_converted).to_dict()
    else:
        gene_list_converted, backend2input, not_converted = convert_IDs(gene_list,input_type, backend_type = backend_type, enrich_against = enrich_against)
    #read in dictionary and the gocam sizes
    #x = pd.read_csv(gcs)
    #gocam_sizes = pd.Series(x.sizes.values,index=x.gocam)
    Dict = csv2dict(id2g,sep = sep)
    pathways2members = reverse_dict(Dict)
    sizes = []
    for k,v in pathways2members.items():
        sizes.append(len(v))
    
    gocam_sizes = pd.Series(sizes, index = pathways2members.keys())    #call enrich()
    ncHGT = False
    if method == 'ncHGT':
        ncHGT = True
    #results: (filtered_out_genes, filtered_list, setID2members_input_uni, setID2members_input, df_display)
    results = enrich(list(gene_list.g), gene_list_converted, backend2input, gocam_sizes, Dict, ncHGT = ncHGT, FDR=FDR, 
                     kegg = kegg, input_type = input_type, backend_type = backend_type,enrich_against= enrich_against)
    
    if display_gene_symbol == True:
        results[4]['shared entities in gocam'] = backend2gene(results[4]['shared entities in gocam'])
        results[4]['shared entities in gocam'] = results[4]['shared entities in gocam'].apply(lambda x: [x_.replace('sset:','set:') for x_ in x])
    if display_input == True:
        #overrides default of display_gene_symbol
        results[4]['shared entities in gocam'] = results[4]['shared entities in gocam'].apply(lambda x: [backend2input.get(x_) for x_ in x])

    if method == 'set' or method == 'ncHGT':
        print(f"Analysis run on {len(results[1])} entities from {len(gene_list)-len(results[0])} out of {len(gene_list)} input genes")
    

    if return_all:
        return (gene_list, *results)
    else:
        return results[4]
    