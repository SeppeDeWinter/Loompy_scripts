import argparse

def average_expression(args):
    import numpy as np
    import pandas as pd
    import loompy
    import itertools
    from tqdm import tqdm

    row_attr = None
    genes_file_path = None
    
    loom_file_path = args.input
    column_attrs = args.column_attributes
    row_attr = args.row_attribute
    out_file_path = args.output
    genes_file_path = args.genes
    batch_s = args.batch_size
    lyrs = args.layers
    ra_use = args.ra_use


    if(row_attr != None and genes_file_path != None):
        raise Exception('Either select row attribute OR a list of genes names (OR none) but not both!')

    with loompy.connect(loom_file_path, 'r', validate = False) as ds:
        cas = [np.unique(ds.ca[column_attrs[x]]) for x in range(len(column_attrs))]
        mean_table = pd.DataFrame(index = [("%s "*len(x)%x).strip() for x in list(itertools.product(*cas))], columns = ds.ra[ra_use])

        if(row_attr != None):
            print("using: ", row_attr)
            Gene_use = ds.ra[ra_use][ds.ra.Selected.astype(bool)]
            print(len(Gene_use))
            for (ix, selection, view) in tqdm(ds.scan(axis=0,items=ds.ra[row_attr].astype(bool), batch_size=batch_s, layers=lyrs), total=int((ds.shape[0]/batch_s)+0.5)):
                for t in itertools.product(*cas):
                    selected_columns = np.prod([ds.ca[column_attrs[x]] == t[x] for x in range(len(column_attrs))], axis = 0).astype(bool)
                    mean_table.loc[("%s "*len(t)%t).strip()][selection] = np.mean(view[:, selected_columns], axis = 1)
        elif(genes_file_path != None):
            with open(genes_file_path) as f:
                selected_genes = f.read().splitlines()
            row_attr = np.in1d(ds.ra[ra_use], selected_genes)
            print("Using: ", sum(row_attr), "of: ",len(selected_genes), "selected genes")
            Gene_use = ds.ra[ra_use][row_attr]

            print(len(Gene_use))
            for (ix, selection, view) in tqdm(ds.scan(axis=0, items=row_attr, batch_size=batch_s, layers=lyrs), total=int((ds.shape[0]/batch_s)+0.5)):
                for t in itertools.product(*cas):
                    selected_columns = np.prod([ds.ca[column_attrs[x]] == t[x] for x in range(len(column_attrs))], axis = 0).astype(bool)
                    mean_table.loc[("%s "*len(t)%t).strip()][selection] = np.mean(view[:, selected_columns], axis = 1)
        else:
            print("Using all genes")
            Gene_use = ds.ra[ra_use]
            print(len(Gene_use))
            for (ix, selection, view) in tqdm(ds.scan(axis=0, batch_size=batch_s, layers=lyrs), total=int((ds.shape[0]/batch_s)+0.5)):
                for t in itertools.product(*cas):
                    selected_columns = np.prod([ds.ca[column_attrs[x]] == t[x] for x in range(len(column_attrs))], axis = 0).astype(bool)
                    mean_table.loc[("%s "*len(t)%t).strip()][selection] = np.mean(view[:, selected_columns], axis = 1)
        
        mean_table = mean_table[Gene_use]
        mean_table.to_csv(out_file_path)

def subsample_cells(args):
    import loompy
    import numpy as np
    from tqdm import tqdm
    n_cells = args.num_cells
    new_file = args.output
    loom_file  = args.input

    with loompy.new(new_file) as dsout:
        with loompy.connect(loom_file, 'r') as ds:
            if(n_cells > ds.shape[1]):
                raise Exception('Number of cells <--num_cells> needs to be smaller to the number of cells in the original loom file!')

            print('subsampling to: ', int((n_cells/ds.shape[1]) * 100), ' percent of original dataset')
            choosen_cells = np.random.choice(np.arange(ds.shape[1]), size = n_cells, replace = False)
            for (ix, selection, view) in tqdm(ds.scan(items=choosen_cells, axis=1), total=int((ds.shape[1]/512)+0.5)):
                dsout.add_columns(view.layers, col_attrs=view.ca, row_attrs=view.ra)

def map_ra(args):
    import loompy
    import pandas as pd
    import numpy as np
    from tqdm import tqdm
    loom_file = args.input
    mapping = args.mapping
    ra = args.row_attr
    new_ra_name = args.new_ra_name

    mapping = pd.read_csv(mapping, header = None, index_col = 0)

    with loompy.connect(loom_file, 'r+', validate=False) as ds:
        source = ds.ra[ra]
        source_in_mapping_bool = np.in1d(source, mapping.iloc[:, 0])
        print('Of the: ', len(source), ' row attributes in the loom file, ', sum(source_in_mapping_bool), ' are in the provided mapping dataframe ')
        source_in_mapping = source[source_in_mapping_bool]
        source_to_target_in_mapping = {x: mapping[mapping.iloc[:, 0] == x].iloc[:, 1].values[0] for x in tqdm(source_in_mapping, total = len(source_in_mapping))}
        target = [source_to_target_in_mapping[x] if x in source_to_target_in_mapping.keys() else 'hs_'+x for x in source]
        ds.ra[new_ra_name] = target

def create_argument_parser():
    parser = argparse.ArgumentParser(description='CLI for doing operations on loom files')

    subparsers = parser.add_subparsers(help='sub-command help')

    # --------------------------------------------------------
    # create the parser for the "average_expression" command
    # --------------------------------------------------------

    parser_avg = subparsers.add_parser('average_expression',
                                            help='Calculates average gene expression for given combination of column attributes')
    parser_avg.add_argument('-i', '--input',
                                 help='<Required> Path to loom file',
                                 required=True)
    parser_avg.add_argument('-o', '--output',
                                 type=argparse.FileType('w'),
                                 help='<Required> Path to write output csv file',
                                 required=True)
    parser_avg.add_argument('-ca', '--column_attributes',
                                 nargs = '+',
                                 help='<Required> Column attribute or combination of column attributes to calculate average gene expression on. For example -ca ClusterID to calculate average gene expression across all clusters. For example -ca ClusterID TimePoint to calculate average gene expression across all combinations of ClusterIDs and timepoints (e.g. cluster_1 timepoint_1, cluster_1 timepoint_2, cluster_2 timepoint_1 and cluster_2 timepoint_2 when the loom files has two ClusterIDs and two TimePoints).',
                                 required=True)
    parser_avg.add_argument('-ra', '--row_attribute',
                                 help='<Optional> row attribute specifying on which genes average expression should be calculated, leave empty for all genes',
                                 required=False)
    parser_avg.add_argument('-g', '--genes',
                                 help='<Optional> file containing gene names from which the average expression has to be calculated, leave empty for all genes',
                                 required=False)
    parser_avg.add_argument('-b', '--batch_size',
                                 help='<Optional> batch size for scanning through loom file, default is 512',
                                 required=False,
                                 default=512,
                                 type = int)
    parser_avg.add_argument('-ly', '--layers',
                                 nargs = '+',
                                 help='<Optional> Layers to calculate average expression on',
                                 required=False,
                                 default=[""])
    parser_avg.add_argument('-r', '--ra_use',
                                help='<Required> which row attribute is used (e.g. gene)',
                                required=True)                                                                                                                               
    parser_avg.set_defaults(func=average_expression)

    # --------------------------------------------------------
    # create the parser for the "subsample_cells" command
    # --------------------------------------------------------
    parser_sample = subparsers.add_parser('subsample_cells',
                                            help='Subsamples loom file to a certain number of cells')
    parser_sample.add_argument('-i', '--input',
                                 help='<Required> Path to loom file',
                                 required=True)
    parser_sample.add_argument('-o', '--output',
                                 help='<Required> Path to write output loom file',
                                 required=True)
    parser_sample.add_argument('-n', '--num_cells',
                                 type = int,
                                 required = True,
                                 help = '<Required> number of cells to sample')                                
    parser_sample.set_defaults(func=subsample_cells)

    # --------------------------------------------------------
    # Create the parser for the "map_ra" command
    # --------------------------------------------------------

    parser_map_ra = subparsers.add_parser('map_ra',
                                            help='Map row attributes of the loom file to some other terms, e.g. ortholog genes')
    parser_map_ra.add_argument('-i', '--input',
                                 help='<Required> Path to loom file (Loom file shoud be writable)',
                                 required=True)  
    parser_map_ra.add_argument('-ra', '--row_attr',
                                 help='<Required> Row attribute to map from (e.g. Gene)',
                                 required=True) 
    parser_map_ra.add_argument('-m', '--mapping',
                                 help='<Required> Path to csv file containing mapping (A 2-column data frame defining variable name mapping. First column is source variable name and second column is target variable name',
                                 required=True)
    parser_map_ra.add_argument('-n', '--new_ra_name',
                                 help='<Required> New of row attribute for the new mapping',
                                 required=True)
    parser_map_ra.set_defaults(func=map_ra)                                                                                               

    return(parser)

def main(argv=None):
    parser = create_argument_parser()
    args = parser.parse_args(args=argv)
    if not hasattr(args, 'func'):
        parser.print_help()
    else:
        args.func(args)


if __name__ == '__main__':
    main()
