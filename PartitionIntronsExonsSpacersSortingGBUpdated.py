import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation
import copy 

##### Input

annotation_file_input = "/Users/SJAnnaldasula/Documents/BGBM/Emma/Melicytus ramiflorus MW238813.gb"
input_gb = SeqIO.read(annotation_file_input, 'genbank')

sequence_file_input = "/Users/SJAnnaldasula/Documents/BGBM/Emma/Melicytus ramiflorus MW238813.fasta"
input_seq = SeqIO.read(sequence_file_input, 'fasta')
#input_seq = input_gb # Uncomment this line and comment the upper two lines if you want to use the sequence provided by the GenBank file

workdir = ""
sample = "Test"

trna = False # Set False if want to remove tRNA from introns and exons  
rrna = False # Set False if want to remove rRNA from introns and exons

            
##### Exons and Introns

def GenBankProperSortingFormat(feature):
    # Define the sorting order: gene first, then CDS
    if feature.type == 'source':
        return (feature.location.start,0)
    elif feature.type == 'misc_feature':
        return (feature.location.start,1)
    elif feature.type == 'gene':
        return (feature.location.start,2)
    elif feature.type == 'CDS':
        return (feature.location.start,3)
    else:
        return (feature.location.start,4)
    
def FastaProperSortingFormat(feature):
    if ("<" in feature.name or ">" in feature.name):
        feature.name = feature.name.replace("<","")
        feature.name = feature.name.replace(">","")
    return int(feature.name)
    
features_to_select = ["CDS"]
if (trna): features_to_select += ["tRNA"]
if (rrna): feature_to_select += ["rRNA"]
    
count = 0
exons = []
introns = []

input_gb.features = sorted(input_gb.features, key=GenBankProperSortingFormat)

for feature_num in range(len(input_gb.features)):
    feature_curr = input_gb.features[feature_num]
    if (feature_curr.type == "source"):
        start,end = feature_curr.location.start,feature_curr.location.end
    
    elif ((feature_curr.type == "gene") and ("pseudo" not in feature_curr.qualifiers)):
        feature = input_gb.features[feature_num+1]
        
        if (feature.type in features_to_select):
            try:
                gene_curr = feature.qualifiers["gene"][0]
            except:
                try:
                    gene_curr = feature.qualifiers["locus_tag"][0]
                except:
                    gene_curr = feature.qualifiers["standard_name"][0]

            exon_prev = None
            for exon_curr in feature.location.parts:
                # Exception for if a gene spans the starting location
                if (exon_prev != None and exon_curr.start == start and exon_prev.end == end):
                    feature_loc_1 = FeatureLocation(exon_prev.start,exon_prev.end,strand=exon_prev.strand) 
                    feature_loc_2 = FeatureLocation(exon_curr.start,exon_curr.end,strand=exon_curr.strand)
                    record = SeqRecord(
                        feature_loc_1.extract(input_seq).seq + feature_loc_2.extract(input_seq).seq,
                        id="%s" %gene_curr,
                        name=str(exon_curr.start),
                        description="Exon %d-%d %s" %(exon_prev.start,exon_curr.end,sample),
                    )
                    exons.append(record)
                else:
                    # Exon
                    if (exon_prev == None and exon_curr.end == end):
                        pass
                    else:
                        feature_loc = FeatureLocation(exon_curr.start,exon_curr.end,strand=exon_curr.strand)
                        record = SeqRecord(
                            feature_loc.extract(input_seq).seq,
                            id="%s" %gene_curr,
                            name=str(exon_curr.start),
                            description="Exon %d-%d %s" %(exon_curr.start,exon_curr.end,sample),
                        )
                        exons.append(record) 

                    #Intron
                    # If RPS12,19, dont make intron
                    if (exon_prev != None and ("rps12" not in feature.qualifiers["gene"][0])):
                        if (exon_curr.strand == 1):
                            feature_loc = FeatureLocation(exon_prev.end,exon_curr.start,strand=exon_curr.strand)
                            record = SeqRecord(
                                feature_loc.extract(input_seq).seq,
                                id="%s" %gene_curr,
                                name=str(exon_curr.start),
                                description="Intron %d-%d %s" %(exon_prev.end,exon_curr.start,sample),
                            )
                        else:
                            feature_loc = FeatureLocation(exon_curr.end,exon_prev.start,strand=exon_curr.strand)
                            record = SeqRecord(
                                feature_loc.extract(input_seq).seq,
                                id="%s" %gene_curr,
                                name=str(exon_curr.end),
                                description="Intron %d-%d %s" %(exon_curr.end,exon_prev.start,sample),
                            )
                        introns.append(record) 

                    exon_prev = exon_curr

exons = sorted(exons, key=FastaProperSortingFormat)
exons_seq = os.path.join(workdir,'%s_Exons.fasta' %sample) 
with open(exons_seq, "w+") as result_file:
    for record in exons:
        SeqIO.write(record, result_file, "fasta")

introns = sorted(introns, key=FastaProperSortingFormat)
introns_seq = os.path.join(workdir,'%s_Introns.fasta' %sample)
with open(introns_seq, "w+") as result_file:
    for record in introns:
        SeqIO.write(record, result_file, "fasta")
        
##### Seperate RPS12 stuff     

rps12_gene_segments_positions = []
rps12_gene_segments = []

remove_features = []
for feature in input_gb.features:
    if((feature.type == "gene") and ("gene" in feature.qualifiers) and ("rps12" in feature.qualifiers["gene"][0])):
        if (len(feature.location.parts) > 1):
            for rps12_gene_segment_postions in feature.location.parts:
                if (rps12_gene_segment_postions not in rps12_gene_segments_positions):
                    rps12_gene_segment = copy.deepcopy(feature)
                    rps12_gene_segment.location = rps12_gene_segment_postions
                    rps12_gene_segments_positions.append(rps12_gene_segment_postions)
                    rps12_gene_segments.append(rps12_gene_segment)
            remove_features.append(feature)

for feature in rps12_gene_segments:
    input_gb.features.append(feature)

for feature in remove_features:
    input_gb.features.remove(feature)

input_gb.features = sorted(input_gb.features, key=GenBankProperSortingFormat)

##### Spacers
 
gene_prev = "None"
gene_prev_end = 0

count = 0
spacers = []
for feature in input_gb.features:
    if (feature.type == "gene"):
        gene_curr_start = feature.location.parts[0].start
        gene_curr_end = feature.location.parts[-1].end
        try:
            gene_curr = feature.qualifiers["gene"][0]
        except:
            try:
                gene_curr = feature.qualifiers["locus_tag"][0]
            except:
                gene_curr = feature.qualifiers["standard_name"][0]
                
        if (gene_prev_end < gene_curr_start and gene_curr_start < gene_curr_end):
            feature_loc = FeatureLocation(gene_prev_end,gene_curr_start,strand=1)
            record = SeqRecord(
                feature_loc.extract(input_seq).seq,
                id="%s+++%s" %(gene_prev.replace(" ",""),gene_curr.replace(" ","")),
                name=str(gene_prev_end),
                description="Spacer %d-%d %s" %(gene_prev_end,gene_curr_start,sample)
            )
            spacers.append(record)
            gene_prev = gene_curr
            gene_prev_end = gene_curr_end
        elif (gene_curr_end < gene_prev_end):
            gene_prev_end = gene_prev_end
            gene_prev = gene_prev
        else:
            gene_prev_end = gene_curr_end
            gene_prev = gene_curr

spacers = sorted(spacers, key=FastaProperSortingFormat)
spacers_seq = os.path.join(workdir,'%s_Spacers.fasta' %sample) 
with open(spacers_seq, "w+") as result_file:
    for record in spacers:
        SeqIO.write(record, result_file, "fasta")