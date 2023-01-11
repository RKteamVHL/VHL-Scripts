def create_cluster_figures(directory, dfs):
    for df_type, df_out in dfs.items():

        clustered = dataframe_snf(df_out).fillna(0)

        create_cluster_summaries(directory, clustered, df_type)
        create_cluster_phenotype_summaries(directory, clustered, df_type)

def create_cluster_summaries(directory, df, analysis_type):
    base_path = os.path.join(directory, analysis_type, "cluster")
    if not os.path.isdir(base_path):
        os.makedirs(base_path)

    for clust_type in ["cluster_labels_best", "cluster_labels_second"]:
        fig_path = os.path.join(base_path, clust_type)
        if not os.path.isdir(fig_path):
            os.makedirs(fig_path)

        properties = ["generalized_phenotype", "grouped_mutation_type", "domain", "aa_change"]
        for prop in properties:
            prop_df = df.set_index(clust_type)[COMPUTED_COLUMNS[prop]].fillna(0)
            plot_cluster_property(prop_df, cluster_column=clust_type, figure_path=fig_path, property_name=prop,
                                  save_csv=True)

        prop = "age"
        prop_df = df.set_index(clust_type)[COMPUTED_COLUMNS[prop]]
        plot_cluster_property(prop_df, cluster_column=clust_type, figure_path=fig_path, property_name=prop,
                              use_mean=True)
        df["codon_start"] = df["codon_start"].astype(int)
        codon = df.set_index("codon_start")


        codon = codon[codon['generalized_mutant_type.missense_variant'] != 0]
        # codon = codon.dropna(subset=['generalized_mutant_type.missense_variant'])
        codon = pd.get_dummies(codon[codon.index.notnull()][clust_type]).sort_index()
        codon = codon.groupby(codon.index).sum()
        codons = pd.DataFrame(index=set(range(0, 214)), columns=codon.columns)
        codons[:] = 0
        codons.loc[set(codon.index)] = codon
        axs = codons.plot(kind="bar", xticks=[], subplots=True, figsize=(12, 10), title=["" for v in codon.columns])
        # for ax in axs:
        #     ax.set_xticks(DOMAIN_TICKS)
        # plt.savefig(os.path.join(fig_path, f'clustered_codon_start.pdf'))
        plt.savefig(os.path.join(fig_path, f'clustered_codon_start.eps'), format='eps')

        df.to_csv(os.path.join(fig_path, f"clustered_out.tsv"), sep='\t')
        plt.close('all')


def create_cluster_phenotype_summaries(directory, df, analysis_type):
    base_path = os.path.join(directory, analysis_type, "cluster")
    for clust_type in ["cluster_labels_best", "cluster_labels_second"]:
        fig_path = os.path.join(base_path, clust_type)
        if not os.path.isdir(fig_path):
            os.makedirs(fig_path)

        df = df.set_index(clust_type)
        for phenotype in COMPUTED_COLUMNS["generalized_phenotype"]:
            fig_path_pheno = os.path.join(fig_path, phenotype)

            df_phen = df[df[phenotype] >= 1]
            properties = ["grouped_mutation_type", "domain", "aa_change"]
            for prop in properties:
                prop_df = df_phen[COMPUTED_COLUMNS[prop]].fillna(0)
                plot_cluster_property(prop_df, cluster_column=clust_type, figure_path=fig_path_pheno,
                                      property_name=prop, ratio_type="across")

            prop_df = df_phen[phenotype].fillna(0)
            df_counts = prop_df.groupby(clust_type).sum()
            df_ratio = df_counts.divide(df_counts.sum())
            ax = df_ratio.plot(kind="bar", figsize=(12, 8))
            autolabel(ax, ax.patches)
            plt.savefig(os.path.join(fig_path, f'clustered_phenotype_ratios.eps'), format='eps')
            # plt.savefig(os.path.join(fig_path, f'clustered_phenotype_ratios.pdf'))
            plt.close('all')