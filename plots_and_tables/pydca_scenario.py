from pydca.meanfield_dca import meanfield_dca

test_file = '/users/bag/hlq763/hbv_covar3/analysis/sim_seq/scenarios/single_mutate_10tip_10.fasta'

mfdca_inst = meanfield_dca.MeanFieldDCA(
    test_file,
    'rna',
    pseudocount=0.5,
    seqid=0.8,
)

mfdca_FN_APC = mfdca_inst.compute_sorted_FN_APC()

for site_pair, score in mfdca_FN_APC[:10]:
    print(site_pair)
    print(score)

for site_pair, score in mfdca_FN_APC[:100]:
    output = ",".join(list(site_pair) + [score])
    print(output)



