import torch

def load_mapping(mapping_file, mapping_type):
    mapping = {}
    file_handle = open(mapping_file)
    for line in file_handle:
        line = line.rstrip().split()
        mapping[line[0]] = int(line[1])

    file_handle.close()
    print('Total number of {} = {}'.format(mapping_type, len(mapping)))
    return mapping

# build mask: matrix (nrows = number of relevant gene set, ncols = number all genes)
# elements of matrix are 1 if the corresponding gene is one of the relevant genes
def create_term_mask(term_direct_gene_map, gene_dim, cuda_id):
    term_mask_map = {}
    for term, gene_set in term_direct_gene_map.items():
        mask = torch.zeros(len(gene_set), gene_dim, requires_grad=False).to(cuda_id)
        for i, gene_id in enumerate(gene_set):
            mask[i, gene_id] = 1
        term_mask_map[term] = mask
    return term_mask_map

def build_input_vector(inputdata, gene_id_mapping):
    inputdata = inputdata.to(torch.int64) # for scatter to work index has to have dtype int64
    genedim = len(gene_id_mapping)
    sample_num = inputdata.size()[0]
    features = torch.zeros((sample_num, genedim), requires_grad=False)
    features.scatter_(1, inputdata, 1)
    return features

def pearson_corr(x, y):
    xx = x - torch.mean(x)
    yy = y - torch.mean(y)

    return torch.sum(xx*yy) / (torch.norm(xx, 2)*torch.norm(yy,2))
        